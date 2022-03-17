#include <sstream>
#include <omp.h>
#include <ctime>
#include "LeesEdwardsBoundaryWorld.h"

#pragma warning(disable : 4996)

using std::stringstream;

void LeesEdwardsBoundaryWorld::prepare(const string& path)
{
	stringstream md;
	md.str("");
	md << "md " << path.c_str();
	system(md.str().c_str());
	checkTimeIncrement();
	writeParticleInformation2Files(path);
	createParticleTimeHistoryFiles(path);
}

void LeesEdwardsBoundaryWorld::takeTimeStep()
{
	modifyParticlePosition();
	applyBackgroudVelocity();
	findPossibleContacts();
	updateContacts();
	removeBackgroudVelocity();
	collectForceAndTorque();
	takeTimeIntegral();
	updateLeesEdwardBoundary();
}


void LeesEdwardsBoundaryWorld::findPossibleContacts()
{
	// Update grid size and numbers
	_xlen = _xmax - _xmin;
	_ylen = _ymax - _ymin;
	int numGridX = static_cast<int>(floor(_xlen / _maxParticleSize));
	int numGridY = static_cast<int>(floor(_ylen / _maxParticleSize));
	int numGrid = numGridX * numGridY;
	double lenGridX = _xlen / numGridX;
	double lenGridY = _ylen / numGridY;
	// Locate particles into grids
	int numGridPerThread = numGrid / NUM_THREADS + 1;
	int numGridPerThreadLE = numGridX / NUM_THREADS + 1;
	/*
	2021.04.23 Jiang comment:
	Besides the normal box to which particles belong to,
	another two boxes are need to apply Lees-Edward boundary condition.
	Contact detection will also deal with particle interactions near L-E boundary seperately.
	*/
	vector<vector<int>> box(numGrid, vector<int>(0));
	vector<vector<int>> boxLE_top(numGridX, vector<int>(0));
	vector<vector<int>> boxLE_bottom(numGridX, vector<int>(0));
	vector<array<int, 2>> particleInGrid(_numParticle, array<int, 2>{-1, -1});
	vector<int> particleInGridLE_top(_numParticle, -1);
	vector<int> particleInGridLE_bottom(_numParticle, -1);
#pragma omp parallel num_threads(NUM_THREADS)
	{
		// Determine particle belongs to which grid
#pragma omp for schedule(dynamic, SIZE_CHUNK) 
		for (int i = 0; i < _numParticle; ++i) {
			int ix = static_cast<int>(floor((_particleHandlers[i]->getPosX() - _xmin) / lenGridX));
			int iy = static_cast<int>(floor((_particleHandlers[i]->getPosY() - _ymin) / lenGridY));
			particleInGrid[i] = { ix, iy };
			// 2021.04.22 LEES EDWARDS BOUNDARY CONDITION ONLY...
			/*
			2021.04.23 Jiang comment:
			If a particle is detected at the very bottom grid of simulation domain, say iy==0,
			this particle (i) should be labelled as "particleInGridLE_top[i]" = ix,
			[!!] note that this "ix" is totally different from the above one, it has a shift due to L-E boundary.
			The same process for particle on the very top grid of simulation domain will be done.
			*/
			if (iy == 0) {
				double tempPosX = _particleHandlers[i]->getPosX() + _totalShift;
				if (tempPosX < _xmin) tempPosX += _xlen;
				else if (tempPosX > _xmax) tempPosX -= _xlen;
				ix = static_cast<int>(floor((tempPosX - _xmin) / lenGridX));
				particleInGridLE_top[i] = ix;
			}
			else if (iy == numGridY - 1) {
				double tempPosX = _particleHandlers[i]->getPosX() - _totalShift;
				if (tempPosX < _xmin) tempPosX += _xlen;
				else if (tempPosX > _xmax) tempPosX -= _xlen;
				ix = static_cast<int>(floor((tempPosX - _xmin) / lenGridX));
				particleInGridLE_bottom[i] = ix;
			} 
			// 2021.04.22 LEES EDWARDS BOUNDARY CONDITION ONLY...
		}
		// Put particles into grid box
		int threadID = omp_get_thread_num();
		int boxStart = threadID * numGridPerThread;
		int boxEnd = (threadID + 1) * numGridPerThread - 1;
		for (int i = 0; i < _numParticle; ++i) {
			int ib = particleInGrid[i][1] * numGridX + particleInGrid[i][0];
			if (ib >= boxStart && ib <= boxEnd) {
				box[ib].push_back(i);
			}
		}
		// LEES 2021.04.22 EDWARDS BOUNDARY CONDITION ONLY...
		/*
		2021.04.23 Jiang comment:
		Considering the convenience for parallel computation,
		particles are first labelled, and after that, located into the box.
		[!!] particleInGridLE_top    -> boxLE_top
			 particleInGridLE_bottom -> boxLE_bottom 
		*/
		boxStart = threadID * numGridPerThreadLE;
		boxEnd = (threadID + 1) * numGridPerThreadLE - 1;
		for (int i = 0; i < _numParticle; ++i) {
			int ix_top = particleInGridLE_top[i];
			if (ix_top > -1) {
				if (ix_top >= boxStart && ix_top <= boxEnd) {
					boxLE_top[ix_top].push_back(i);
				}
			}
			int ix_bottom = particleInGridLE_bottom[i];
			if (ix_bottom > -1) {
				if (ix_bottom >= boxStart && ix_bottom <= boxEnd) {
					boxLE_bottom[ix_bottom].push_back(i);
				}
			}
		}
		// 2021.04.22 LEES EDWARDS BOUNDARY CONDITION ONLY...
	}
	// Find possible contacts and construct possible contact list
#pragma omp parallel num_threads(NUM_THREADS)
	{
#pragma omp for schedule(dynamic, SIZE_CHUNK) 
		for (int i = 0; i < _particleHandlers.size(); ++i) {
			// Todo: check the performance with <vector>
			// Done: <array> performs definitely better than <vector>
			auto pList = _possibleContacts[i].begin();
			auto pList_top = _possibleContactsLE_top[i].begin();
			auto pList_bottom = _possibleContactsLE_bottom[i].begin();
			BaseParticle* ip = _particleHandlers[i];
			int ix = particleInGrid[i][0];
			int iy = particleInGrid[i][1];
			int ib = iy * numGridX + ix;
			for (int jx0 = -1; jx0 <= 1; ++jx0) {
				int jx = ix + jx0;
				if (jx == -1) { jx = numGridX - 1; }
				else if (jx == numGridX) { jx = 0; }
				for (int jy0 = -1; jy0 <= 1; ++jy0) {
					int jy = iy + jy0;
					// 2021.04.22 LEES EDWARDS BOUNDARY CONDITION ONLY...
					/*
					2021.04.23 Jiang Comment:
					[!!] For normal contact detection (also normal periodic contact detection)
					Main particle i, in box ib [ix,iy]
					Sub  particle j, in box jb [jx,jy]
					[!!] But for L-E boundary condition
					Main particle i, in box ib [ix,iy]
					Sub  particle j, in "boxLE_bottom[jx]" and in "boxLE_top[jx]"
					When the main particle i is detected in the bottom box, say iy == 0,
					AND its neighbour box label jy == -1,
					then, instead of normal "box", "boxLE_bottom[jx]" should be used.
					*/
					if (jy == -1 || jy == numGridY) {
						if (jy == -1) {
							if (boxLE_bottom[jx].empty()) continue;
							for (int j = 0; j < boxLE_bottom[jx].size(); ++j) {
								BaseParticle* jp = _particleHandlers[boxLE_bottom[jx][j]];
								if (ip->getID() >= jp->getID()) continue;
								if (ip->getUpLevel() == jp->getUpLevel()) continue;
								*pList_bottom = jp;
								pList_bottom++;
							}
						}
						else if (jy == numGridY) {
							if (boxLE_top[jx].empty()) continue;
							for (int j = 0; j < boxLE_top[jx].size(); ++j) {
								BaseParticle* jp = _particleHandlers[boxLE_top[jx][j]];
								if (ip->getID() >= jp->getID()) continue;
								if (ip->getUpLevel() == jp->getUpLevel()) continue;
								*pList_top = jp;
								pList_top++;
							}
						}
					}
					// 2021.04.22 LEES EDWARDS BOUNDARY CONDITION ONLY...
					else {
						int jb = jy * numGridX + jx;
						if (box[jb].empty()) continue;
						for (int j = 0; j < box[jb].size(); ++j) {
							BaseParticle* jp = _particleHandlers[box[jb][j]];
							if (ip->getID() >= jp->getID()) continue;
							if (ip->getUpLevel() == jp->getUpLevel()) continue;
							*pList = jp;
							pList++;
							if (pList == _possibleContacts[i].end()) {
								printf("##__FUNCTION__exceed array range...\n");
								exit(-2);
							}
						}
					}
				}
			}
			*pList = nullptr;
			*pList_top = nullptr;
			*pList_bottom = nullptr;
		}
	}
}

void LeesEdwardsBoundaryWorld::updateContacts()
{
	//printf("Function called ... LeesEdwardsBoundaryWorld::updateContactsWithLeesEdwardBoundary\n");
#pragma omp parallel num_threads(NUM_THREADS)
	{
#pragma omp for schedule(dynamic, SIZE_CHUNK)
		for (int i = 0; i < _numParticle; ++i) {
			BaseParticle* ip = _particleHandlers[i];
			auto pList = _possibleContacts[i].begin();
			while (*pList != nullptr) {
				BaseParticle* jp = *pList;
				vector2d offset(0.0, 0.0);
				if (jp->getPosX() - ip->getPosX() > _xlen / 2.0) offset.setX(-_xlen);
				else if (jp->getPosX() - ip->getPosX() < -_xlen / 2.0) offset.setX(_xlen);
				if (jp->getPosY() - ip->getPosY() > _ylen / 2.0) offset.setY(-_ylen);
				else if (jp->getPosY() - ip->getPosY() < -_ylen / 2.0) offset.setY(_ylen);
				ip->calculateContactForceAndTorque(jp, _dt, offset);
				pList++;
			}
			// 2021.04.22 LEES EDWARDS BOUNDARY CONDITION ONLY...
			/*
			2021.04.23 Jiang comment:
			Different from normal possible contact detection process,
			possible contact detection process with L-E boundary condition need position and velocity shifts.
			Next code block finds the amount of shifts.
			*/
			auto pList_top = _possibleContactsLE_top[i].begin();
			while (*pList_top != nullptr) {
				BaseParticle* jp = *pList_top;
				vector2d offset(_totalShift, _ylen);
				double tempPosX = offset.x() + jp->getPosX();
				if (tempPosX < _xmin) {
					offset.moveX(_xlen);
					tempPosX += _xlen;
				}
				else if (tempPosX > _xmax) {
					offset.moveX(-_xlen);
					tempPosX -= _xlen;
				}
				if (tempPosX - ip->getPosX() > _xlen / 2.0) {
					offset.moveX(-_xlen);
				}
				else if (tempPosX - ip->getPosX() < -_xlen / 2.0) {
					offset.moveX(_xlen);
				}
				//vector2d offsetVel(-2.0 * jp->getVelX(), 0.0);
				vector2d offsetVel(_shiftVel, 0.0);
				ip->calculateContactForceAndTorque(jp, _dt, offset, offsetVel);
				pList_top++;
			}
			auto pList_bottom = _possibleContactsLE_bottom[i].begin();
			while (*pList_bottom != nullptr) {
				BaseParticle* jp = *pList_bottom;
				vector2d offset(-_totalShift, -_ylen);
				double tempPosX = offset.x() + jp->getPosX();
				if (tempPosX < _xmin) {
					offset.moveX(_xlen);
					tempPosX += _xlen;
				}
				else if (tempPosX > _xmax) {
					offset.moveX(-_xlen);
					tempPosX -= _xlen;
				}
				if (tempPosX - ip->getPosX() > _xlen / 2.0) {
					offset.moveX(-_xlen);
				}
				else if (tempPosX - ip->getPosX() < -_xlen / 2.0) {
					offset.moveX(_xlen);
				}
				//vector2d offsetVel(-2.0 * jp->getVelX(), 0.0);
				vector2d offsetVel(-_shiftVel, 0.0);
				ip->calculateContactForceAndTorque(jp, _dt, offset, offsetVel);
				pList_bottom++;
			}// 2021.04.22 LEES EDWARDS BOUNDARY CONDITION ONLY...
		}
	}
}

void LeesEdwardsBoundaryWorld::updateLeesEdwardBoundary()
{
	_shiftVel = _shiftStrainRate * _ylen;
	_totalStrain += _shiftStrainRate * _dt;
	_totalShift += _shiftVel * _dt;
	if (_totalShift > _xmax) {
		_totalShift -= _xlen;
	}
	else if (_totalShift < _xmin) {
		_totalShift += _xlen;
	}
}


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// It is not necessary to modify particle velocity here, 
// when a overall background velocity modification is applied.
// If both particle velocity modification and background velocity are applied,
// particles near the top-bottom boundaries would be modified twice!!
void LeesEdwardsBoundaryWorld::modifyParticlePosition()
{
#pragma omp parallel num_threads(NUM_THREADS)
	{
#pragma omp for schedule(static)
		for (int i = 0; i < _numParticle; ++i) {
			BaseParticle* ip = _particleHandlers[i];
			// 2021.04.22 LEES EDWARDS BOUNDARY CONDITION ONLY...
			if (ip->getPosY() < _ymin) {
				ip->shiftPosY(_ylen);
				ip->shiftPosX(_totalShift);
				//ip->reverseVelX();
				//ip->setVelX(ip->getVelX() + _shiftVel);
				// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				// It is not necessary to modify particle velocity here, 
				// when a overall background velocity modification is applied.
				// If both particle velocity modification and background velocity are applied,
				// particles near the top-bottom boundaries would be modified twice!!
			}
			else if (ip->getPosY() > _ymax) {
				ip->shiftPosY(-_ylen);
				ip->shiftPosX(-_totalShift);
				//ip->reverseVelX();
				//ip->setVelX(ip->getVelX() - _shiftVel);
			}
			if (ip->getPosX() < _xmin) ip->shiftPosX(_xlen);
			else if (ip->getPosX() > _xmax) ip->shiftPosX(-_xlen);
		}
	}
}

void LeesEdwardsBoundaryWorld::createParticleTimeHistoryFiles(const string& path)
{
	stringstream fileName;

	fileName << path << "\\output_velocity.dat";
	FILE* file_velocity = fopen(fileName.str().c_str(), "w");
	files.push_back(file_velocity);
	fileName.str("");

	fileName << path << "\\output_position.dat";
	FILE* file_position = fopen(fileName.str().c_str(), "w");
	files.push_back(file_position);
	fileName.str("");

	fileName << path << "\\output_omega.dat";
	FILE* file_orientationVelocity = fopen(fileName.str().c_str(), "w");
	files.push_back(file_orientationVelocity);
	fileName.str("");

	fileName << path << "\\output_rotation.dat";
	FILE* file_orientation = fopen(fileName.str().c_str(), "w");
	files.push_back(file_orientation);
	fileName.str("");

	fileName << path << "\\output_contact.dat";
	FILE* file_contact = fopen(fileName.str().c_str(), "w");
	files.push_back(file_contact);
	fileName.str("");

	fileName << path << "\\output_boundary.dat";
	FILE* file_boundary = fopen(fileName.str().c_str(), "w");
	files.push_back(file_boundary);
	fileName.str("");

	fileName << path << "\\output_kineticEnergy.dat";
	FILE* file_kineticEnergy = fopen(fileName.str().c_str(), "w");
	files.push_back(file_kineticEnergy);
	fileName.str("");

	fileName << path << "\\output_boundaryShift.dat";
	FILE* file_boundaryShift = fopen(fileName.str().c_str(), "w");
	files.push_back(file_boundaryShift);
	fileName.str("");

	fileName << path << "\\output_stress.dat";
	FILE* file_stress = fopen(fileName.str().c_str(), "w");
	files.push_back(file_stress);
	fileName.str("");
}

void LeesEdwardsBoundaryWorld::writeParticleTimeHistory2Files()
{
	for (auto it : _particleHandlers) it->print2File_velocity(files[0]);
	fprintf(files[0], "\n");

	for (auto it : _particleHandlers) it->print2File_position(files[1]);
	fprintf(files[1], "\n");

	for (auto it : _particleHandlers) it->print2File_rotationVelocity(files[2]);
	fprintf(files[2], "\n");

	for (auto it : _particleHandlers) it->print2File_rotation(files[3]);
	fprintf(files[3], "\n");

	fprintf(files[4], "%d\n", countActucalContacts());
	for (auto it : _particleHandlers) it->print2File_contact(files[4]);

	print2File_boundary(files[5]);

	print2File_kineticEnergy(files[6]);

	print2File_boundaryShift(files[7]);

	print2File_stress(files[8]);
}

void LeesEdwardsBoundaryWorld::writeParticleTimeHistory2SeparateFiles(const string& path)
{
	stringstream md;
	md << "md " << path.c_str();
	system(md.str().c_str());

	stringstream fileName;

	fileName << path << "\\output_velocity.dat";
	FILE* file_velocity = fopen(fileName.str().c_str(), "w");
	fileName.str("");
	for (auto it : _particleHandlers) it->print2File_velocity(file_velocity);
	fprintf(file_velocity, "\n");
	fclose(file_velocity);

	fileName << path << "\\output_position.dat";
	FILE* file_position = fopen(fileName.str().c_str(), "w");
	fileName.str("");
	for (auto it : _particleHandlers) it->print2File_position(file_position);
	fprintf(file_position, "\n");
	fclose(file_position);

	fileName << path << "\\output_omega.dat";
	FILE* file_orientationVelocity = fopen(fileName.str().c_str(), "w");
	fileName.str("");
	for (auto it : _particleHandlers) it->print2File_rotationVelocity(file_orientationVelocity);
	fprintf(file_orientationVelocity, "\n");
	fclose(file_orientationVelocity);

	fileName << path << "\\output_rotation.dat";
	FILE* file_orientation = fopen(fileName.str().c_str(), "w");
	fileName.str("");
	for (auto it : _particleHandlers) it->print2File_rotation(file_orientation);
	fprintf(file_orientation, "\n");
	fclose(file_orientation);

	fileName << path << "\\output_contact.dat";
	FILE* file_contact = fopen(fileName.str().c_str(), "w");
	fileName.str("");
	fprintf(file_contact, "%d\n", countActucalContacts());
	for (auto it : _particleHandlers) it->print2File_contact(file_contact);
	fclose(file_contact);

	fileName << path << "\\output_boundary.dat";
	FILE* file_boundary = fopen(fileName.str().c_str(), "w");
	fileName.str("");
	print2File_boundary(file_boundary);
	fclose(file_boundary);

	fileName << path << "\\output_kineticEnergy.dat";
	FILE* file_kineticEnergy = fopen(fileName.str().c_str(), "w");
	fileName.str("");
	print2File_kineticEnergy(file_kineticEnergy);
	fclose(file_kineticEnergy);

	fileName << path << "\\output_stress.dat";
	FILE* file_stress = fopen(fileName.str().c_str(), "w");
	fileName.str("");
	print2File_stress(file_stress);
	fclose(file_stress);

	fileName << path << "\\output_boundaryShift.dat";
	FILE* file_boundaryShift = fopen(fileName.str().c_str(), "w");
	fileName.str("");
	print2File_boundaryShift(file_boundaryShift);
	fclose(file_boundaryShift);
}

void LeesEdwardsBoundaryWorld::print2File_boundaryShift(FILE* file)
{
	fprintf(file, "%le\n", _totalShift);
}


void LeesEdwardsBoundaryWorld::setShiftStrainRate(const double& shiftStrainRate)
{
	_shiftStrainRate = shiftStrainRate;
	_shiftVel = _shiftStrainRate * _ylen;
}

void LeesEdwardsBoundaryWorld::setParticleVelocityAccordingToPosition()
{
	for (auto& p : _particleHandlers) {
		double scale = p->getPosY() / _ymax;
		p->setVelX(scale * _shiftVel);
	}
}

void LeesEdwardsBoundaryWorld::moveParticlePositionAccordingToPosition()
{
	for (auto& p : _particleHandlers) {
		//double scale = p->getPosY() / _ymax;
		//p->shiftPosX(scale * _shiftVel * _dt);
		p->shiftPosX(p->getPosY() * _shiftStrainRate * _dt);
	}
}

void LeesEdwardsBoundaryWorld::applyBackgroudVelocity()
{
#pragma omp parallel num_threads(NUM_THREADS)
	{
#pragma omp for schedule(static) 
		for (int i = 0; i < _numParticle; ++i) {
			BaseParticle* p = _particleHandlers[i];
			p->setVelX(p->getVelX() + _shiftStrainRate * p->getPosY());
		}
	}
}

void LeesEdwardsBoundaryWorld::removeBackgroudVelocity()
{
#pragma omp parallel num_threads(NUM_THREADS)
	{
#pragma omp for schedule(static) 
		for (int i = 0; i < _numParticle; ++i) {
			BaseParticle* p = _particleHandlers[i];
			p->setVelX(p->getVelX() - _shiftStrainRate * p->getPosY());
		}
	}
}



void LeesEdwardsBoundaryWorld::print2Screen_worldState(const int& iStep)
{
	static clock_t timerTotal = clock();
	static clock_t timer = clock();
	timer = clock() - timer;
	printf("Current step = %d, Time used per 1000 steps = %lf, Total time used = %lf\n", iStep,
		(static_cast<double>(timer) / CLOCKS_PER_SEC),
		((static_cast<double>(clock()) - timerTotal) / CLOCKS_PER_SEC));
	timer = clock();
	printf("\tCurrent solid fraction = %lf\n", getSolidFraction());
	printf("\tPossible contact number = %d\n", countPossibleContacts());
	printf("\tActual contact number = %d\n", countActucalContacts());
	printf("\tSystem kinetic energy = %le\n", getKineticEnergyPerParticle());
	printf("\tSystem hydro stress = %le\n", getStress().trace() / 2.0);
	printf("\tCurrent shear strain = %le\n", _totalShift / _ylen);
}
