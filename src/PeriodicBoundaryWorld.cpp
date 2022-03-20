#include <omp.h>
#include <vector>
#include <sstream>
#include <ctime>
#include "PeriodicBoundaryWorld.h"

#pragma warning(disable : 4996)

using std::stringstream;
using std::vector;

void PeriodicBoundaryWorld::findPossibleContacts()
{
	// Update grid size and numbers
	_xlen = _xmax - _xmin;
	_ylen = _ymax - _ymin;
	int numGridX = static_cast<int>(floor(_xlen / _maxParticleSize));
	int numGridY = static_cast<int>(floor(_ylen / _maxParticleSize));
	int numGrid = numGridX * numGridY;
	//printf("Logger msg: Total number of Grid %d\n", numGrid);
	double lenGrixX = _xlen / numGridX;
	double lenGrixY = _ylen / numGridY;
	// Locate particles into grids
	int numGridPerThread = numGrid / NUM_THREADS + 1;
	vector<vector<int>> box(numGrid, vector<int>(0));
	vector<array<int, 2>> particleInGrid(_numParticle, array<int, 2>{-1, -1});
#pragma omp parallel num_threads(NUM_THREADS)
	{
		// Determine particle belongs to which grid
#pragma omp for schedule(dynamic, SIZE_CHUNK) 
		for (int i = 0; i < _numParticle; ++i) {
			int ix = static_cast<int>(floor((_particleHandlers[i]->getPosX() - _xmin) / lenGrixX));
			int iy = static_cast<int>(floor((_particleHandlers[i]->getPosY() - _ymin) / lenGrixY));
			particleInGrid[i] = { ix, iy };
		}
		// Put particles into grid box
		int threadID = omp_get_thread_num();
		int boxStart = threadID * numGridPerThread;
		int boxEnd = (threadID + 1) * numGridPerThread - 1;
		for (int i = 0; i < _numParticle; ++i) {
			int ib = particleInGrid[i][1] * numGridX + particleInGrid[i][0];
			if (ib < boxStart || ib > boxEnd) continue;
			box[ib].push_back(i);
		}
	}
	// Find possible contacts and construct possible contact list
#pragma omp parallel num_threads(NUM_THREADS)
	{
#pragma omp for schedule(dynamic, SIZE_CHUNK) 
		for (int i = 0; i < _particleHandlers.size(); ++i) {
			// Todo: check the performance with <vector>
			auto pList = _possibleContacts[i].begin();
			BaseParticle* ip = _particleHandlers[i];
			if ((ip->getSettings()) & SETTING::SKIP_CONTACT_DETECTION) continue;
			int ix = particleInGrid[i][0];
			int iy = particleInGrid[i][1];
			int ib = iy * numGridX + ix;
			for (int jx0 = -1; jx0 <= 1; ++jx0) {
				int jx = ix + jx0;
				if (jx == -1) { jx = numGridX - 1; }
				else if (jx == numGridX) { jx = 0; }
				for (int jy0 = -1; jy0 <= 1; ++jy0) {
					int jy = iy + jy0;
					if (jy == -1) { jy = numGridY - 1; }
					else if (jy == numGridY) { jy = 0; }
					int jb = jy * numGridX + jx;
					if (box[jb].empty()) continue;
					for (int j = 0; j < box[jb].size(); ++j) {
						BaseParticle* jp = _particleHandlers[box[jb][j]];
						if ((jp->getSettings()) & SETTING::SKIP_CONTACT_DETECTION) continue;
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
			*pList = nullptr;
		}
	}
}

void PeriodicBoundaryWorld::updateContacts()
{
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
		}
	}
}

void PeriodicBoundaryWorld::prepare(const string& path)
{
	stringstream md;
	md << "md " << path.c_str();
	system(md.str().c_str());
	checkTimeIncrement();
	createParticleTimeHistoryFiles(path);
	writeParticleInformation2Files(path);
}

void PeriodicBoundaryWorld::takeTimeStep()
{
	modifyParticlePosition();
	findPossibleContacts();
	updateContacts();
	collectForceAndTorque();
	takeTimeIntegral();
	//resetVelocityOfRattlers();
}

void PeriodicBoundaryWorld::modifyParticlePosition()
{
#pragma omp parallel num_threads(NUM_THREADS)
	{
#pragma omp for schedule(static)
		for (int i = 0; i < _numParticle; ++i) {
			BaseParticle* ip = _particleHandlers[i];
			if (ip->getPosX() < _xmin) ip->shiftPosX(_xlen);
			else if (ip->getPosX() > _xmax) ip->shiftPosX(-_xlen);
			if (ip->getPosY() < _ymin) ip->shiftPosY(_ylen);
			else if (ip->getPosY() > _ymax) ip->shiftPosY(-_ylen);
		}
	}
}

void PeriodicBoundaryWorld::createParticleTimeHistoryFiles(const string& path)
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

	fileName << path << "\\output_stress.dat";
	FILE* file_stress = fopen(fileName.str().c_str(), "w");
	files.push_back(file_stress);
	fileName.str("");
}

void PeriodicBoundaryWorld::writeParticleTimeHistory2Files()
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

	print2File_stress(files[7]);
}

void PeriodicBoundaryWorld::writeParticleTimeHistory2SeparateFiles(const string& path)
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

}

void PeriodicBoundaryWorld::scaleBoundarySizeToGivenSolidFraction(const double& initSolidFrac)
{
	double currentSolidFrac = getSolidFraction();
	printf("Solid fraction before scaling = %lf\n", getSolidFraction());
	double scale = sqrt(currentSolidFrac / initSolidFrac);
	scaleAllParticlePosition(scale, scale);
	scaleBoundarySize(scale);
	printf("Solid fraction after scaling  = %lf\n", getSolidFraction());
}

void PeriodicBoundaryWorld::scaleBoundarySizeAndParticlePosition(const double& scale)
{
	scaleAllParticlePosition(scale, scale);
	scaleBoundarySize(scale);
	printf("Solid fraction after scaling  = %lf\n", getSolidFraction());
}

void PeriodicBoundaryWorld::print2Screen_worldState(const int& iStep)
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
	printf("\tSystem stress = %le\n", updateTotalStress().trace()/2.0);
}

void PeriodicBoundaryWorld::print2File_stress(FILE* file)
{
	fprintf(file, "%le\t%le\t%le\t%le\n", _totalStress.xx(), _totalStress.xy(), _totalStress.yx(), _totalStress.yy());
}

