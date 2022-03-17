#include "QuadronAnalysis.h"

#include <iostream>
using std::cout;
using std::endl;
#include <sstream>
using std::stringstream;
#include <numeric> 
using std::iota;
#include <algorithm> 
using std::find;
using std::sort;

//#define DEBUG_MESSAGE
#ifdef DEBUG_MESSAGE
#define LOGGER(x) x
#else
#define LOGGER(x)
#endif // DEBUG_MESSAGE

#define STRING_LENGTH 6

double QuadronSystem::_xmin = 0.0;
double QuadronSystem::_xmax = 0.0;
double QuadronSystem::_ymin = 0.0;
double QuadronSystem::_ymax = 0.0;
double QuadronSystem::_xlen = 0.0;
double QuadronSystem::_ylen = 0.0;
double QuadronSystem::_xlenHalf = 0.0;
double QuadronSystem::_ylenHalf = 0.0;
double QuadronSystem::_offset = 0.0;

unordered_map<string, pair<int, int>> QuadronSystem::oldContactCellList;
unordered_map<string, pair<int, int>> QuadronSystem::newContactCellList;
unordered_map<string, QuadronSystem::Polygon> QuadronSystem::oldCellList;
unordered_map<string, QuadronSystem::Polygon> QuadronSystem::newCellList;
unordered_map<string, pair<int, int>> QuadronSystem::cellSurviveTime;

int QuadronSystem::NUM_THREADS = 4;
int QuadronSystem::SIZE_CHUNK = 100;

void QuadronSystem::setPeriodicBoundarySize(const double& xmin, const double& xmax, const double& ymin, const double& ymax)
{
	_xmin = xmin;
	_xmax = xmax;
	_ymin = ymin;
	_ymax = ymax;
	_xlen = _xmax - _xmin;
	_ylen = _ymax - _ymin;
	_xlenHalf = _xlen / 2.0;
	_ylenHalf = _ylen / 2.0;
}

void QuadronSystem::setLeesEdwardsBoundaryOffset(const double& offset)
{
	_offset = offset;
}

// As an example: reference position (particle), target position (contact)
void QuadronSystem::modifyPosition(double& tempX, double& tempY, const double& x, const double& y)
{
	if (y - tempY > _ylenHalf) { // Particle: top, Contact: bottom.
		tempY += _ylen;
		tempX += _offset;
	}
	else if (y - tempY < -_ylenHalf) { // particle: bottom, contact: top.
		tempY -= _ylen;
		tempX -= _offset;
	}
	if (x - tempX > _xlenHalf) { // particle: right, contact: left.
		tempX += _xlen;
	}
	else if (x - tempX < -_xlenHalf) { // particle: left, contact: right.
		tempX -= _xlen;
	}

	if (x - tempX > _xlenHalf) { // particle: right, contact: left.
		tempX += _xlen;
	}
	else if (x - tempX < -_xlenHalf) { // particle: left, contact: right.
		tempX -= _xlen;
	}
}

void QuadronSystem::print()
{
	printf("xmin = %le\nxmax = %le\nymin = %le\nymax = %le\n", _xmin, _xmax, _ymin, _ymax);
}

void QuadronSystem::quadronAnalysisInThisStep(const string path, vector<BaseParticle*>& particleHandlers)
{
	stringstream md;
	md << "md " << path.c_str();
	system(md.str().c_str());

	std::stringstream fileName;
	fileName << path << "quadronAnalysis_cell.dat";
	FILE* fileCell = fopen(fileName.str().c_str(), "w");
	fileName.str("");
	fileName << path << "quadronAnalysis_quadron.dat";
	FILE* fileQuadron = fopen(fileName.str().c_str(), "w");
	fileName.str("");
	fileName << path << "quadronAnalysis_neighborCell.dat";
	FILE* fileNeighborCell = fopen(fileName.str().c_str(), "w");
	fileName.str("");
	fileName << path << "quadronAnalysis_coordinationNum.dat";
	FILE* fileCoordinationNum = fopen(fileName.str().c_str(), "w");
	QuadronSystem::QuadronAnalysis quadronAnalysis(particleHandlers);
	quadronAnalysis.runAnalysis();
	quadronAnalysis.printNeighborCell2File(fileNeighborCell);
	quadronAnalysis.printCells2File(fileCell);
	quadronAnalysis.printCoordinationNum2File(fileCoordinationNum);
	quadronAnalysis.printQuadrons2File(fileQuadron);
	fclose(fileCell);
	fclose(fileQuadron);
	fclose(fileNeighborCell);
	fclose(fileCoordinationNum);
}

QuadronSystem::Particle::Particle()
{
}

QuadronSystem::Particle::Particle(const int& id, const double& x, const double& y, const double& vol) :
	_id(id), _x(x), _y(y), _volume(vol)
{}

void QuadronSystem::Particle::pushBackSurroundingParticle(QuadronSystem::Particle* p)
{
	_pParticles.push_back(p);
}

void QuadronSystem::Particle::pushBackSurroundingContact(QuadronSystem::Contact* c)
{
	_pContacts.push_back(c);
}

void QuadronSystem::Particle::setAsInactiveParticle()
{
	LOGGER(printf("Particle %d is labelled as rattler particle...\n", _id));
	_isActive = false;
}

void QuadronSystem::Particle::setAllSurroundingContactsAsInactive()
{
	for (QuadronSystem::Contact* contact : _pContacts) {
		// check all surrounding particles if it has only one "ACTIVE" contact,
		// if it does, inactivate it.
		if (contact->_pParticle1 == this) {
			int count = 0;
			for (auto& contact0 : contact->_pParticle2->_pContacts) {
				if (contact0->_isActive) count++;
			}
			if (count == 1) contact->_pParticle2->setAsInactiveParticle();
		}
		else if (contact->_pParticle2 == this) {
			int count = 0;
			for (auto& contact0 : contact->_pParticle1->_pContacts) {
				if (contact0->_isActive) count++;
			}
			if (count == 1) contact->_pParticle1->setAsInactiveParticle();
		}
		contact->setAsInactiveContact();
	}
}

//void QuadronSystem::Particle::setAllSurroundingContactsAsInactive()
//{
//	for (QuadronSystem::Contact* contact : _pContacts) {
//		contact->setAsInactiveContact();
//		// check all surrounding particles if it has only one contact,
//		// if it does, inactivate it.
//		if (contact->_pParticle1 == this) {
//			if (contact->_pParticle2->_pContacts.size() == 1) {
//				contact->_pParticle2->setAsInactiveParticle();
//			}
//		}
//		else if (contact->_pParticle2 == this) {
//			if (contact->_pParticle1->_pContacts.size() == 1) {
//				contact->_pParticle1->setAsInactiveParticle();
//			}
//		}
//	}
//}

void QuadronSystem::Particle::rearrangeContactsClockwisely()
{
	//if (_pContacts.size() < 3) return;
	vector<double> localTheta(_pContacts.size());
	for (int i = 0; i < _pContacts.size(); ++i) {
		double tempX = _pContacts[i]->_x;
		double tempY = _pContacts[i]->_y;
		// Modify contact position (tempX, tempY) with respect to particle position (_x, _y)
		QuadronSystem::modifyPosition(tempX, tempY, _x, _y);
		// Calculate angle of contact point
		localTheta[i] = std::atan2(tempY - _y, tempX - _x);
	}
	vector<int> indices(localTheta.size());
	iota(indices.begin(), indices.end(), 0);
	sort(indices.begin(), indices.end(), 
		[&](int A, int B) -> bool 
		{ return localTheta[A] < localTheta[B]; });
	for (int i = 0; i < indices.size(); ++i) {
		if (_pContacts[indices[i]]->_isActive) {
			_pContactsSorted.push_back(_pContacts[indices[i]]);
			_pParticlesSorted.push_back(_pParticles[indices[i]]);
			_link.push_back(1);
		}
	}
	constructPolygon();
}

void QuadronSystem::Particle::constructPolygon()
{
	LOGGER(printf("Function called...\n"));
	for (QuadronSystem::Contact* contact : _pContactsSorted) {
		double tempX = contact->_x;
		double tempY = contact->_y;
		modifyPosition(tempX, tempY, _x, _y);
		_polygon.pushBackVertex(tempX, tempY);
		LOGGER(printf("Vertec (%le, %le) is pushed back\n", tempX, tempY));
	}
	if (_polygon._x.size() > 0) {
		_polygon.calculateProperties();
		LOGGER(printf("Polygon area = %le\n", _polygon._area));
	}
}

void QuadronSystem::Particle::calculateStress()
{
	_stress = matrix2d(0.0, 0.0, 0.0, 0.0);
	for (QuadronSystem::Contact* contact : _pContactsSorted) {
		double tempX = contact->_x;
		double tempY = contact->_y;
		QuadronSystem::modifyPosition(tempX, tempY, _x, _y);
		vector2d branchVector(tempX - _x, tempY - _y);
		if (contact->_pParticle1 == this) {
			_stress += matrix2d::dyadic(branchVector, contact->_force);
		}
		else {
			_stress += matrix2d::dyadic(branchVector, -contact->_force);
		}
	}
	for (QuadronSystem::Quadron* quadron : _pQuadrons) {
		_totalArea += quadron->_polygon._area;
	}
	_stress /= _totalArea;
}

QuadronSystem::Contact::Contact(const int& id, const double& x, const double& y, 
	QuadronSystem::Particle* p1, QuadronSystem::Particle* p2,
	const vector2d& force) :
	_id(id), _x(x), _y(y), _pParticle1(p1), _pParticle2(p2), _force(force)
{}

void QuadronSystem::Contact::setAsInactiveContact()
{
	_isActive = false;
}

QuadronSystem::Quadron::Quadron(const int id, 
	QuadronSystem::Particle* p, QuadronSystem::Cell* c, 
	QuadronSystem::Contact* c1, QuadronSystem::Contact* c2):
	_id(id), _pParticle(p), _pCell(c), _pContact1(c1), _pContact2(c2)
{

}

void QuadronSystem::Quadron::constructPolygon()
{
	const double refX = _pCell->_polygon._cx;
	const double refY = _pCell->_polygon._cy;
	double tempX, tempY;
	
	tempX = _pContact1->_x;
	tempY = _pContact1->_y;
	modifyPosition(tempX, tempY, refX, refY);
	_polygon.pushBackVertex(tempX, tempY);

	tempX = _pParticle->_polygon._cx;
	tempY = _pParticle->_polygon._cy;
	modifyPosition(tempX, tempY, refX, refY);
	_polygon.pushBackVertex(tempX, tempY);

	tempX = _pContact2->_x;
	tempY = _pContact2->_y;
	modifyPosition(tempX, tempY, refX, refY);
	_polygon.pushBackVertex(tempX, tempY);

	_polygon.pushBackVertex(_pCell->_polygon._cx, _pCell->_polygon._cy);
	_polygon.calculateProperties();
}

void QuadronSystem::Quadron::setPointerToParticlesAndCells()
{
	_pParticle->_pQuadrons.push_back(this);
	_pCell->_pQuadrons.push_back(this);
}



QuadronSystem::Cell::Cell(QuadronSystem::Particle* particle, QuadronSystem::Contact* contact)
{
	_pParticles.push_back(particle);
	_pContacts.push_back(contact);
}

void QuadronSystem::Cell::calculateStress()
{
	_stress = matrix2d(0.0, 0.0, 0.0, 0.0);
	for (QuadronSystem::Particle* particle : _pParticles) {
		_stress += (particle->_stress * particle->_totalArea);
		_totalArea += particle->_totalArea;
	}
	_stress /= _totalArea;
}

void QuadronSystem::Cell::pushBackNewParticle(QuadronSystem::Particle* particle)
{
	_pParticles.push_back(particle);
}

void QuadronSystem::Cell::pushBackNewContact(QuadronSystem::Contact* contact)
{
	_pContacts.push_back(contact);
}

void QuadronSystem::Cell::pushBackNewEdge(QuadronSystem::Particle* particle, QuadronSystem::Contact* contact)
{
	_pParticles.push_back(particle);
	_pContacts.push_back(contact);
}

void QuadronSystem::Cell::setPointerToContacts()
{
	for (int i = 0; i < _pContacts.size() - 1; ++i) {
		_pContacts[i]->_pCell.push_back(this);
	}
	/*for (QuadronSystem::Contact* contact : _pContacts) {
		contact->_pCell.push_back(this);
	}*/
}

void QuadronSystem::Cell::constructPolygon()
{
	double refX = _pParticles[0]->_x;
	double refY = _pParticles[0]->_y;
	for (int i = 0; i < _pContacts.size() - 1; ++i) {
		double tempX = _pContacts[i]->_x;
		double tempY = _pContacts[i]->_y;
		modifyPosition(tempX, tempY, refX, refY);
		_polygon.pushBackVertex(tempX, tempY);
	}
	if (_polygon._x.size() > 0) {
		_polygon.calculateProperties();
		LOGGER(printf("Polygon area = %le\n", _polygon._area));
	}
}

QuadronSystem::QuadronAnalysis::QuadronAnalysis(const vector<BaseParticle*>& P)
{
	// Generate particles from DEM particles...
	for (const auto& p : P) particles.push_back(Particle(p->getID(), p->getPosX(), p->getPosY(),p->getVolume()));

	// Generate contacts from DEM contacts...
	int count_c = 0;
	for (const auto& p : P) {
		int p1 = p->getID();
		for (const auto& c0 : p->getContactInfo()) {
			int p2 = c0.first;
			const auto& c = c0.second;
			contacts.push_back(Contact(
				count_c, c._contactPoint.x(), c._contactPoint.y(), 
				&particles[p1], &particles[p2], c._force));
			count_c++;
		}
	}

	// Set surrounding particles and contacts to each particle...
	for (QuadronSystem::Contact& contact : contacts) {
		contact._pParticle1->pushBackSurroundingContact(&contact);
		contact._pParticle1->pushBackSurroundingParticle(contact._pParticle2);
		contact._pParticle2->pushBackSurroundingContact(&contact);
		contact._pParticle2->pushBackSurroundingParticle(contact._pParticle1);
	}
}

QuadronSystem::QuadronAnalysis::QuadronAnalysis(const string path)
{
	// Generate particles from output files...
	ifstream filePos;
	filePos.open(path + "output_position.dat");
	ifstream fileRadius;
	fileRadius.open(path + "output_radius.dat");
	ifstream fileAspectRatio;
	fileAspectRatio.open(path + "output_aspectRatio.dat");
	int count_p = 0;
	while (1) {
		double x, y, r, k;
		filePos >> x >> y;
		fileRadius >> r;
		fileAspectRatio >> k;
		if (filePos.eof()) break;
		particles.push_back(Particle(count_p, x, y, M_PI * r * r / k));
		count_p++;
		//cout << count_p << " " << x << " " << y << endl;
	}
		//std::cin.get();
	filePos.close();
	fileRadius.close();
	fileAspectRatio.close();
	cout << particles.size() << " particles are generated...\n";

	// Generate contacts from output contacts...
	ifstream fileContact;
	fileContact.open(path + "output_contact.dat");
	int numContact;
	fileContact >> numContact;
	for (int i = 0; i < numContact; ++i) {
		int p1, p2;
		fileContact >> p1 >> p2;
		double x, y;
		fileContact >> x >> y;
		double fnx, fny, ftx, fty;
		fileContact >> fnx >> fny >> ftx >> fty;
		vector2d f(fnx + ftx, fny + fty);
		contacts.push_back(Contact(i, x, y, &particles[p1], &particles[p2], f));
		//cout << i << " " << x << " " << y << endl;
	}
		//std::cin.get();
	cout << contacts.size() << " contacts are generated...\n";

	// Set surrounding particles and contacts to each particle...
	for (QuadronSystem::Contact& contact : contacts) {
		contact._pParticle1->pushBackSurroundingContact(&contact);
		contact._pParticle1->pushBackSurroundingParticle(contact._pParticle2);
		contact._pParticle2->pushBackSurroundingContact(&contact);
		contact._pParticle2->pushBackSurroundingParticle(contact._pParticle1);
	}
}



QuadronSystem::QuadronAnalysis::QuadronAnalysis(vector<double>& vol, ifstream& filePos, ifstream& fileContact)
{
	for (int i = 0; i < vol.size(); ++i)
	{
		double x, y;
		filePos >> x >> y;
		particles.push_back(Particle(i, x, y, vol[i]));
	}

	// Generate contacts from DEM contacts...
	int numContact;
	fileContact >> numContact;
	for (int i = 0; i < numContact; ++i) {
		int p1, p2;
		fileContact >> p1 >> p2;
		double x, y;
		fileContact >> x >> y;
		double fnx, fny, ftx, fty;
		fileContact >> fnx >> fny >> ftx >> fty;
		vector2d f(fnx + ftx, fny + fty);
		contacts.push_back(Contact(i, x, y, &particles[p1], &particles[p2], f));
		//cout << i << " " << x << " " << y << endl;
	}
	//std::cin.get();
	//cout << contacts.size() << " contacts are generated...\n";

	// Set surrounding particles and contacts to each particle...
	for (QuadronSystem::Contact& contact : contacts) {
		contact._pParticle1->pushBackSurroundingContact(&contact);
		contact._pParticle1->pushBackSurroundingParticle(contact._pParticle2);
		contact._pParticle2->pushBackSurroundingContact(&contact);
		contact._pParticle2->pushBackSurroundingParticle(contact._pParticle1);
	}
}


void QuadronSystem::QuadronAnalysis::runAnalysis()
{
	inactivateRattlerParticles();
	removeOpenChainParticles();
	rearrangeContactsClockwisely();
	constructCellsFromParticlesAndContacts();
	giveCellID();
	constructQuadronsFromCells();

	calculateParticleStress();
	calculateCellStress();
}

// This function construct the cell list using HASH MAP...	
// Each cell has a hash key: a string...
// The string starts from the smallest particle ID in this cell and rotate clockwisely...
// The hash key is unique so it can be used to find the difference of cells between 2 time steps...
void QuadronSystem::QuadronAnalysis::makeCellList()
{
	oldCellList = newCellList;
	newCellList.clear();
	newCellList.reserve(cells.size());
	for (Cell& cell : cells) {
		int minID = cell._pParticles[0]->_id;
		int index = 0;
		for (int i = 0; i < cell._cellOrder; ++i) {
			if (cell._pParticles[i]->_id < minID) {
				minID = cell._pParticles[i]->_id;
				index = i;
			}
		}
		string key;
		for (int i = index; i < cell._cellOrder; ++i) {
			string key0 = to_string(cell._pParticles[i]->_id);
			key += string(STRING_LENGTH - key0.length(), '-') + key0;
		}
		for (int i = 0; i < index; ++i) {
			string key0 = to_string(cell._pParticles[i]->_id);
			key += string(STRING_LENGTH - key0.length(), '-') + key0;
		}
		newCellList.emplace(key, cell._polygon);
	}
}


void QuadronSystem::QuadronAnalysis::makeContactCellList()
{
	oldContactCellList = newContactCellList;
	newContactCellList.clear();
	newContactCellList.reserve(contacts.size());
	for (Contact& contact : contacts) {
		if (!contact._isActive) continue;
		string key1 = to_string(contact._pParticle1->_id);
		string key2 = to_string(contact._pParticle2->_id);
		string key = string(STRING_LENGTH - key1.length(), '0') + key1 + " " + string(STRING_LENGTH - key2.length(), '0') + key2;
		newContactCellList.emplace(key, pair<int, int>(contact._pCell.front()->_cellOrder, contact._pCell.back()->_cellOrder));
	}
}

void QuadronSystem::QuadronAnalysis::compareContactCellLists(FILE* fileB, FILE* fileC, const bool& isPrint2File)
{
	vector<string> forwardBrokenContacts, backwardBrokenContacts;
	vector<pair<int, int>> forwardCombinedCells, backwardCombinedCells;
	for (auto& contactCell : oldContactCellList) {
		if (!newContactCellList.count(contactCell.first)) {
			forwardBrokenContacts.push_back(contactCell.first);
			forwardCombinedCells.push_back(contactCell.second);
		}
	}
	for (auto& contactCell : newContactCellList) {
		if (!oldContactCellList.count(contactCell.first)) {
			backwardBrokenContacts.push_back(contactCell.first);
			backwardCombinedCells.push_back(contactCell.second);
		}
	}
	/*printf("%d contacts broke from previous step to current step...\n", forwardBrokenContacts.size());
	printf("%d contacts broke from current step to previous step...\n", backwardBrokenContacts.size());
	printf("current number of cells = %d\n", cells.size());
	int previousNumOfCell = cells.size() - backwardBrokenContacts.size() + forwardBrokenContacts.size();
	printf("previous number of cells = %d\n", previousNumOfCell);*/
	
	if (isPrint2File) {
		for (auto& forword : forwardCombinedCells) {
			fprintf(fileB, "%d\t%d\n", forword.first, forword.second);
		}
		for (auto& backword : backwardCombinedCells) {
			fprintf(fileC, "%d\t%d\n", backword.first, backword.second);
		}
	}
	fflush(fileB);
	fflush(fileC);
}


// This function compares the previous cells with current cells...
// Use 'key' to check if current(previous) cells exist in previous(current) time step...
// Check all pairs (previous, current) if the previous cell center is inside current cell...
// If yes, then link these two cells...
// !! IMPORTANT: If the number of linked cells is larger than 1, then this cell is the bigger one...
// The cells that linked to this cell are the smaller one...
void QuadronSystem::QuadronAnalysis::compareCellLists(FILE* fileB, FILE* fileC)
{
	vector<Polygon> oldCells, newCells;
	vector<string> oldCellLabels, newCellLabels;
	for (auto& oldCell : oldCellList) {
		if (!newCellList.count(oldCell.first)) {
			oldCells.push_back(oldCell.second);
			oldCellLabels.push_back(oldCell.first);
		}
	}
	for (auto& newCell : newCellList) {
		if (!oldCellList.count(newCell.first)) {
			newCells.push_back(newCell.second);
			newCellLabels.push_back(newCell.first);
		}
	}

	vector<vector<Polygon*>> oldCellLinker(oldCells.size(), vector<Polygon*>());
	vector<vector<Polygon*>> newCellLinker(newCells.size(), vector<Polygon*>());
	for (int i = 0; i < oldCells.size(); ++i) {
		for (int j = 0; j < newCells.size(); ++j) {
			double tempX, tempY;
			tempX = oldCells[i]._cx;
			tempY = oldCells[i]._cy;
			modifyPosition(tempX, tempY, newCells[j]._cx, newCells[j]._cy);
			if (newCells[j].pointInsidePolygon(tempX, tempY)) {
				newCellLinker[j].push_back(&oldCells[i]);
			}
			tempX = newCells[j]._cx;
			tempY = newCells[j]._cy;
			modifyPosition(tempX, tempY, oldCells[i]._cx, oldCells[i]._cy);
			if (oldCells[i].pointInsidePolygon(tempX, tempY)) {
				oldCellLinker[i].push_back(&newCells[j]);
			}
		}
	}
	for (int i = 0; i < oldCells.size(); ++i) {
		if (oldCellLinker[i].size() == 2) {
			fprintf(fileB, "%d\t", static_cast<int>(oldCells[i]._x.size()));
			for (auto& newCell : oldCellLinker[i]) {
				fprintf(fileB, "%d\t", static_cast<int>(newCell->_x.size()));
			}
			fprintf(fileB, "\n");
		}
	}
	for (int i = 0; i < newCells.size(); ++i) {
		if (newCellLinker[i].size() == 2) {
			fprintf(fileC, "%d\t", static_cast<int>(newCells[i]._x.size()));
			for (auto& oldCell : newCellLinker[i]) {
				fprintf(fileC, "%d\t", static_cast<int>(oldCell->_x.size()));
			}
			fprintf(fileC, "\n");
		}
	}

	fflush(fileB);
	fflush(fileC);
}

void QuadronSystem::QuadronAnalysis::checkCellSurviveTime(FILE* file)
{
	// add new cells and accumulate old cells' surviving time
	for (auto& cell : newCellList) {
		string key = cell.first;
		auto existingCell = cellSurviveTime.find(key);
		if (existingCell != cellSurviveTime.end()) {
			existingCell->second.first++;
		}
		else {
			cellSurviveTime.emplace(key, pair<int, int>(1, cell.second._x.size()));
		}
	}
	// find dying cells and output
	for (auto& cell : cellSurviveTime) {
		string key = cell.first;
		auto dyingCell = newCellList.find(key);
		if (dyingCell == newCellList.end()) {
			fprintf(file, "%d\t%d\n", cell.second.first, cell.second.second);
			cellSurviveTime.erase(key);
		}
	}
}

void QuadronSystem::QuadronAnalysis::inactivateRattlerParticles()
{
	for (QuadronSystem::Particle& particle : particles) {
		if (particle._pContacts.empty()) particle.setAsInactiveParticle();
	}
}

QuadronSystem::Particle* QuadronSystem::QuadronAnalysis::findParticleWithOnlyOneContact()
{
	for (QuadronSystem::Particle& particle : particles) {
		if (!particle._isActive) continue;
		int count = 0;
		for (QuadronSystem::Contact* contact : particle._pContacts) {
			if (contact->_isActive) count++;
		}
		if (count == 1) return &particle;
	}
	return nullptr;
}

QuadronSystem::Particle* QuadronSystem::QuadronAnalysis::findParticleWithLessThanThreeContacts()
{
	for (QuadronSystem::Particle& particle : particles) {
		if (!particle._isActive) continue;
		int count = 0;
		for (QuadronSystem::Contact* contact : particle._pContacts) {
			if (contact->_isActive) count++;
		}
		if (count < 3) return &particle;
	}
	return nullptr;
}

void QuadronSystem::QuadronAnalysis::calculateCellStress()
{
	for (QuadronSystem::Cell& cell : cells) {
		cell.calculateStress();
	}
}

double QuadronSystem::QuadronAnalysis::calculateMeanCoordinationNumber()
{
	int numParticle = 0;
	int numContact = 0;
	for (const auto& particle : particles) if (particle._isActive) numParticle++;
	for (const auto& contact : contacts) if (contact._isActive) numContact++;
	if (numParticle == 0)
		return 0.0;
	else
		return static_cast<double>(numContact) / numParticle * 2.0;
}

double QuadronSystem::QuadronAnalysis::calculateMeanOriginalCoordinationNumber()
{
	int numParticle = 0;
	int numContact = 0;
	for (const auto& particle : particles) numParticle++;
	for (const auto& contact : contacts) numContact++;
	return static_cast<double>(numContact) / numParticle * 2.0;
}

double QuadronSystem::QuadronAnalysis::calculateRattlerFreeSolidFraction()
{
	double solidFrac = 0.0;
	for (const auto& p : particles) {
		if (p._isActive) solidFrac += p._volume;
	}
	solidFrac /= (QuadronSystem::_xlen * QuadronSystem::_ylen);
	return solidFrac;
}

double QuadronSystem::QuadronAnalysis::calculateSolidFraction()
{
	double solidFrac = 0.0;
	for (const auto& p : particles) {
		solidFrac += p._volume;
	}
	solidFrac /= (QuadronSystem::_xlen * QuadronSystem::_ylen);
	return solidFrac;
}

vector<double> QuadronSystem::QuadronAnalysis::calculateDetailCoordinationNumber()
{
	vector<double> cn(10, 0);
	for (const auto& particle : particles) {
		int numContact = 0;
		for (const auto& contact : particle._pContacts) {
			if (contact->_isActive) numContact++;
		}
		if (numContact >= 10) continue;
		cn[numContact]++;
	}
	return cn;
}

vector<double> QuadronSystem::QuadronAnalysis::calculateDetailOriginalCoordinationNumber()
{
	vector<double> cn(10, 0);
	for (const auto& particle : particles) {
		int numContact = 0;
		for (const auto& contact : particle._pContacts) {
			numContact++;
		}
		if (numContact >= 10) continue;
		cn[numContact]++;
	}
	return cn;
}

double QuadronSystem::QuadronAnalysis::checkQudronArea()
{
	double sumArea = 0.0;
	for (QuadronSystem::Quadron& quadron : quadrons) {
		sumArea += quadron._polygon._area;
	}
	sumArea /= (QuadronSystem::_xlen * QuadronSystem::_ylen);
	//printf("Total area of quadrons / Total space area = %le\n", sumArea);
	return sumArea;
}

double QuadronSystem::QuadronAnalysis::checkParticleArea()
{
	double sumArea = 0.0;
	for (QuadronSystem::Particle& particle : particles) {
		sumArea += particle._polygon._area;
	}
	sumArea /= (QuadronSystem::_xlen * QuadronSystem::_ylen);
	return sumArea;
}

double QuadronSystem::QuadronAnalysis::checkCellArea()
{
	double sumArea = 0.0;
	for (const QuadronSystem::Cell& cell : cells) {
		sumArea += cell._polygon._area;
	}
	sumArea /= (QuadronSystem::_xlen * QuadronSystem::_ylen);
	return sumArea;
}

int QuadronSystem::QuadronAnalysis::getCellNum()
{
	return cells.size();
}

void QuadronSystem::QuadronAnalysis::printQuadrons2File(const std::string& path)
{
	std::stringstream fileName;
	fileName << path << "\\quadronAnalysis_quadron.dat";
	FILE* file = fopen(fileName.str().c_str(), "w");
	printQuadrons2File(file);
	fclose(file);
}

void QuadronSystem::QuadronAnalysis::printQuadrons2File(FILE* file)
{
	fprintf(file, "%d\n", static_cast<int>(quadrons.size()));
	for (QuadronSystem::Quadron& quadron : quadrons) {
		fprintf(file, "%d\t%le\n", quadron._cellOrder, quadron._polygon._area);
		fprintf(file, "%le\t", quadron._polygon._cx);
		for (double& x : quadron._polygon._x) fprintf(file, "%le\t", x);
		fprintf(file, "\n");
		fprintf(file, "%le\t", quadron._polygon._cy);
		for (double& y : quadron._polygon._y) fprintf(file, "%le\t", y);
		fprintf(file, "\n");
	}
	fflush(file);
}

void QuadronSystem::QuadronAnalysis::printCells2File(const std::string& path)
{
	std::stringstream fileName;
	fileName << path << "\\quadronAnalysis_cell.dat";
	FILE* file = fopen(fileName.str().c_str(), "w");
	printCells2File(file);
	fclose(file);
}

void QuadronSystem::QuadronAnalysis::printCells2File(FILE* file)
{
	double dummyX = 0;
	double dummyY = 0;
	fprintf(file, "%d\n", static_cast<int>(cells.size()));
	for (QuadronSystem::Cell& cell : cells) {
		fprintf(file, "%d\t%le\n", cell._cellOrder, cell._polygon._area);
		if (isnormal(cell._stress.xx())) {
			fprintf(file, "%le\t%le\t%le\t%le\n",
				cell._stress.xx(), cell._stress.xy(),
				cell._stress.yx(), cell._stress.yy());
		}
		else {
			fprintf(file, "%le\t%le\t%le\t%le\n",
				0.0, 0.0, 0.0, 0.0);
		}
		
		fprintf(file, "%le\t", cell._polygon._cx);
		for (double& x : cell._polygon._x) fprintf(file, "%le\t", x);
		for (int i = 0; i < cell._cellOrder; ++i) {
			double tempX = cell._pParticles[i]->_x;
			modifyPosition(tempX, dummyY, cell._polygon._cx, 0);
			fprintf(file, "%le\t", tempX);
		}
		fprintf(file, "\n");
		fprintf(file, "%le\t", cell._polygon._cy);
		for (double& y : cell._polygon._y) fprintf(file, "%le\t", y);
		for (int i = 0; i < cell._cellOrder; ++i) {
			double tempY = cell._pParticles[i]->_y;
			modifyPosition(dummyX, tempY, 0, cell._polygon._cy);
			fprintf(file, "%le\t", tempY);
		}
		fprintf(file, "\n");
	}
	fflush(file);
}

void QuadronSystem::QuadronAnalysis::printParticle2File(FILE* file)
{
	double dummyX = 0;
	double dummyY = 0;
	fprintf(file, "%d\n", static_cast<int>(particles.size()));
	for (QuadronSystem::Particle& p : particles) {
		if (p._pContactsSorted.size() == 0) {
			fprintf(file, "0\t0\t0\n");
			continue;
		}
		fprintf(file, "%d\t%le\t%le\n", static_cast<int>(p._pContactsSorted.size()), p._polygon._area, p._totalArea);
		fprintf(file, "%le\t%le\t%le\t%le\n",
			p._stress.xx(), p._stress.xy(),
			p._stress.yx(), p._stress.yy());
		fprintf(file, "%le\t", p._x);
		for (double& x : p._polygon._x) fprintf(file, "%le\t", x);
		for (auto& q : p._pQuadrons) {
			double tempX = q->_polygon._x[3];
			modifyPosition(tempX, dummyY, p._polygon._cx, 0);
			fprintf(file, "%le\t", tempX);
		}
		fprintf(file, "\n");
		fprintf(file, "%le\t", p._y);
		for (double& y : p._polygon._y) fprintf(file, "%le\t", y);
		for (auto& q : p._pQuadrons) {
			double tempY = q->_polygon._y[3];
			modifyPosition(dummyX, tempY, 0, p._polygon._cy);
			fprintf(file, "%le\t", tempY);
		}
		fprintf(file, "\n");
		for (auto& c : p._pContactsSorted) {
			if (c->_pParticle1 == &p) {
				fprintf(file, "%le\t%le\n", c->_force.x(), c->_force.y());
			}
			else {
				fprintf(file, "%le\t%le\n", -c->_force.x(), -c->_force.y());
			}
		}
	}
}

void QuadronSystem::QuadronAnalysis::printNeighborCell2File(FILE* file)
{
	fprintf(file, "%d\t%d\n", 0, 0);
	for (auto& contact : contacts) {
		if (contact._pCell.size() == 2) {
			fprintf(file, "%d\t%d\t%d\t%d\n", 
				contact._pCell[0]->_cellOrder, contact._pCell[1]->_cellOrder, 
				contact._pCell[0]->_id, contact._pCell[1]->_id);
		}
		//else {
			//printf("contact with the number of neighboring cell: %d\n", contact._pCell.size());
		//}
	}
	fflush(file);
}

void QuadronSystem::QuadronAnalysis::printCoordinationNum2File(FILE* file)
{
	fprintf(file, "%le\t%le\t", calculateRattlerFreeSolidFraction(), calculateMeanCoordinationNumber());
	const vector<double> coordinationNum = calculateDetailCoordinationNumber();
	for (const auto& cn : coordinationNum) {
		fprintf(file, "%le\t", cn);
	}
	fprintf(file, "\n");
	fflush(file);
	return;
}

void QuadronSystem::QuadronAnalysis::printOriginalCoordinationNum2File(FILE* file)
{
	//cout << "calculate original coordination number" << endl;
	fprintf(file, "%le\t%le\t", calculateSolidFraction(), calculateMeanOriginalCoordinationNumber());
	vector<double> coordinationNum = calculateDetailOriginalCoordinationNumber();
	for (const auto& cn : coordinationNum) {
		fprintf(file, "%le\t", cn);
	}
	fprintf(file, "\n");
	fflush(file);
	//cout << "Finish: calculate original coordination number" << endl;
}

void QuadronSystem::QuadronAnalysis::printRattlerParticle2File(FILE* file)
{
	for (const auto& p : particles) {
		if (p._isActive) fprintf(file, "1\n");
		else if (!p._isActive) fprintf(file, "0\n");
	}
}

void QuadronSystem::QuadronAnalysis::printCellMemberParticle2File(FILE* file)
{
	fprintf(file, "%d\n", static_cast<int>(cells.size()));
	for (const auto& cell : cells) {
		for (const auto& particle : cell._pParticles) {
			fprintf(file, "%d\t", particle->_id);
		}
		fprintf(file, "\n");
	}
}


void QuadronSystem::QuadronAnalysis::removeOpenChainParticles()
{
	while (Particle* particle = findParticleWithOnlyOneContact()) {
		particle->setAsInactiveParticle();
		particle->setAllSurroundingContactsAsInactive();
	}
}

void QuadronSystem::QuadronAnalysis::rearrangeContactsClockwisely()
{
	for (QuadronSystem::Particle& particle : particles) {
		particle.rearrangeContactsClockwisely();
	}
}

 //===============================================
 //MAIN ALGORITHM OF CELL(QUADRON) CONSTRUCTION...
 //===============================================
void QuadronSystem::QuadronAnalysis::constructCellsFromParticlesAndContacts_anotherVersion()
{
	for (int ip = 0; ip < particles.size(); ++ip) {
		Particle* particle = &particles[ip];
		for (size_t ic = 0; ic < particle->_pContactsSorted.size(); ++ic) {
			LOGGER(printf("=========== new cell constructed ===========\n"));
			LOGGER(printf("particle id = %d\n", particle->_id));
			// Initialize a new cell...
			// push back one particle and one contact
			Contact* contact = particle->_pContactsSorted[ic];
			Cell cell(particle, contact);
			// Check if the first path is accessible...
			// if is, then push back the second particle
			// Note that the second contact is not pushed back here
			if (particle->_link[ic]) {
				// Break the link from 1st particle to 2nd particle...
				// This should be done before set the 'particle' pointer to the next one
				particle->_link[ic] = 0;
				size_t index = ic - 1;
				if (index < 0) index = particle->_pContactsSorted.size() - 1;
				contact = particle->_pContactsSorted[index];
				cell.pushBackNewContact(contact);
			}
			else {
				continue;
			}
			// Find next path until the cell head is the same as the cell tail...
			// This loop begins with 2nd particle and 1st contact
			while (cell._pContacts.back() != cell._pContacts.front()) {

				if (contact->_pParticle1 == particle) particle = contact->_pParticle1;
				else if (contact->_pParticle2 == particle) particle = contact->_pParticle2;
				cell.pushBackNewParticle(particle);
				
				vector<QuadronSystem::Contact*>::iterator it = find(
					particle->_pContactsSorted.begin(), particle->_pContactsSorted.end(), contact);
				
				size_t index = it - particle->_pContactsSorted.begin();
				if (particle->_link[index]) {
					// Break the link from 1st particle to 2nd particle...
					// This should be done before set the 'particle' pointer to the next one
					particle->_link[index] = 0;
					index--;
					if (index < 0) index = particle->_pContactsSorted.size() - 1;
					contact = particle->_pContactsSorted[index];
					cell.pushBackNewContact(contact);
				}
			}
			// Here, the new cell is constructed...
			// But the length of particle list and contact list are not the same
			// For convienence, push one more contact, which is the same as 1st one, back to the contact list
			cell._cellOrder = static_cast<int>(cell._pParticles.size());
			cell.pushBackNewParticle(cell._pParticles.front());
			// Before push back the new "cell" to "cells" vector, check if it is closed...
			// And also check if its cell order is greater than 2
			if (cell._pParticles.back() == cell._pParticles.front() && cell._cellOrder >= 2) {
				cell.constructPolygon();
				cells.push_back(cell);
			}
		}
	}
	LOGGER(printf("[QuadronAnalysis::Message]__%d cells are constructed...\n", static_cast<int>(cells.size())));
}
void QuadronSystem::QuadronAnalysis::giveCellID()
{
	int count = 0;
	for (auto& cell : cells) {
		cell._id = count;
		count++;
	}
}
void QuadronSystem::QuadronAnalysis::constructCellsFromParticlesAndContacts()
{
	for (int ip = 0; ip < particles.size(); ++ip) {
		Particle* particle = &particles[ip];
		for (int ic = 0; ic < particle->_pContactsSorted.size(); ++ic) {
			LOGGER(printf("=========== new cell constructed ===========\n"));
			LOGGER(printf("particle id = %d\n", particle->_id));
			// Initialize a new cell...
			// push back one particle and one contact
			Contact* contact = particle->_pContactsSorted[ic];
			Cell cell(particle, contact);
			// Check if the first path is accessible...
			// if is, then push back the second particle
			// Note that the second contact is not pushed back here
			if (particle->_link[ic]) {
				// Break the link from 1st particle to 2nd particle...
				// This should be done before set the 'particle' pointer to the next one
				particle->_link[ic] = 0;
				particle = particle->_pParticlesSorted[ic];
				cell.pushBackNewParticle(particle);
				LOGGER(printf("particle id = %d\n", particle->_id));
			}
			else {
				continue;
			}
			// Find next path until the cell head is the same as the cell tail...
			// This loop begins with 2nd particle and 1st contact
			while (cell._pParticles.back() != cell._pParticles.front()) {
				// Find the index of 1st contact on 2nd particle...
				// std::find returns the iterator of 1st contact on 2nd particle
				// So it = it + 1 means, find the 2nd contact on 2nd particle
				// Note that:
				// From the point of view of 2nd particle center, find the contact clockwisely
				// From the point of view of this cell, fin the contact anti-clockwisely
				vector<QuadronSystem::Contact*>::iterator it = find(
					particle->_pContactsSorted.begin(), particle->_pContactsSorted.end(), contact) + 1;
				// !! IMPORTANT !! modify the iterator if it points to the last contact in the list...
				if (it == particle->_pContactsSorted.end()) {
					it = particle->_pContactsSorted.begin();
				}
				// Find the actual index of the list...
				size_t index = it - particle->_pContactsSorted.begin();
				LOGGER(printf("index = %d, length of contactSortedList = %d\n", index, particle->_pContactsSorted.size()));
				// Change the contact pointer from 1st contact to 2nd contact and push it back...
				contact = particle->_pContactsSorted[index];
				cell.pushBackNewContact(contact);
				// Check if the new path is accessible...
				if (particle->_link[index]) {
					// Break the link from 2nd particle to 3rd particle...
					// This should be done before set the 'particle' pointer to the next one
					particle->_link[index] = 0;
					particle = particle->_pParticlesSorted[index];
					// Push back 3rd particle...
					cell.pushBackNewParticle(particle);
					LOGGER(printf("particle id = %d\n", particle->_id));
				}
				else {
					// If the path is inaccessible, break the while loop...
					break;
				}
			}
			// Here, the new cell is constructed...
			// But the length of particle list and contact list are not the same
			// For convienence, push one more contact, which is the same as 1st one, back to the contact list
			cell._cellOrder = static_cast<int>(cell._pContacts.size());
			cell.pushBackNewContact(cell._pContacts.front());
			// Before push back the new "cell" to "cells" vector, check if it is closed...
			// And also check if its cell order is greater than 2
			if (cell._pParticles.back() == cell._pParticles.front() && cell._cellOrder >= 2) {
				cell.constructPolygon();
				cells.push_back(cell);
			}
		}
	}
	for (QuadronSystem::Cell& cell : cells) {
		cell.setPointerToContacts();
	}
	LOGGER(printf("[QuadronAnalysis::Message]__%d cells are constructed...\n", static_cast<int>(cells.size())));
}

void QuadronSystem::QuadronAnalysis::constructQuadronsFromCells()
{
	int countID = 0;
	for (QuadronSystem::Cell& cell : cells) {
		LOGGER(printf("Start construct quadron...\n"));
		LOGGER(printf("%d quadrons will be constructed\n", cell._cellOrder));
		for (size_t i = 0; i < cell._cellOrder; ++i) {
			// Construct new quadron...
			QuadronSystem::Quadron quadron(countID++, cell._pParticles[i + 1], &cell, cell._pContacts[i], cell._pContacts[i + 1]);
			quadron.constructPolygon();
			quadron._cellOrder = cell._cellOrder;
			// Push new quadron back to "_quadron" list...
			quadrons.push_back(quadron);
		}
	}
	for (QuadronSystem::Quadron& quadron : quadrons) {
		quadron.setPointerToParticlesAndCells();
	}
}

void QuadronSystem::QuadronAnalysis::calculateParticleStress()
{
	for (QuadronSystem::Particle& particle : particles) {
		particle.calculateStress();
	}
}

QuadronSystem::Polygon::Polygon()
{
}

void QuadronSystem::Polygon::pushBackVertex(const double& x, const double& y)
{
	_x.push_back(x);
	_y.push_back(y);
}

void QuadronSystem::Polygon::calculateProperties()
{
	if (_x.size() == 2) {
		_cx = (_x.front() + _x.back()) / 2.0;
		_cy = (_y.front() + _y.back()) / 2.0;
		_area = 0.0;
	}
	else {
		for (size_t i = 0; i < _x.size() - 1; ++i) {
			double temp = _x[i] * _y[i + 1] - _x[i + 1] * _y[i];
			_area += temp;
			_cx += (_x[i] + _x[i + 1]) * temp;
			_cy += (_y[i] + _y[i + 1]) * temp;
		}
		double temp = (_x.back() * _y.front() - _x.front() * _y.back());
		_area += temp;
		_cx += (_x.front() + _x.back()) * temp;
		_cy += (_y.front() + _y.back()) * temp;
		_area *= 0.5;
		_cx /= (6.0 * _area);
		_cy /= (6.0 * _area);
		_area = abs(_area);
	}
	LOGGER(printf("Polygon area = %le\n", _area));
	LOGGER(printf("Polygon center = (%le, %le)\n", _cx, _cy));
}

int QuadronSystem::Polygon::pointInsidePolygon(const double& x, const double& y)
{
	int nvert = static_cast<int>(_x.size());
	int i, j, c = 0;
	for (i = 0, j = nvert - 1; i < nvert; j = i++) {
		if (((_y[i] > y) != (_y[j] > y)) &&
			(x < (_x[j] - _x[i]) * (y - _y[i]) / (_y[j] - _y[i]) + _x[i]))
			c = !c;
	}
	return c;
}

void QuadronSystem::Polygon::print2File(FILE* file)
{
	for (double& x : _x) {
		fprintf(file, "%le\t", x);
	}
	fprintf(file, "\n");
	for (double& y : _y) {
		fprintf(file, "%le\t", y);
	}
	fprintf(file, "\n");
}
