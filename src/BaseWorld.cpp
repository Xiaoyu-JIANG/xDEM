#include <omp.h>
#include <sstream>
#include "BaseWorld.h"

using std::stringstream;

void BaseWorld::setWorldBoundary(const double& xmin, const double& xmax, const double& ymin, const double& ymax)
{
	_xmin = xmin;
	_xmax = xmax;
	_ymin = ymin;
	_ymax = ymax;
	_xlen = _xmax - _xmin;
	_ylen = _ymax - _ymin;
}

void BaseWorld::collectForceAndTorque()
{
	_totalForce.clear();
	_totalForce.resize(_numParticle, vector2d(0.0, 0.0));
	_totalTorque.clear();
	_totalTorque.resize(_numParticle, 0.0);
	_numContactPerParticle.clear();
	_numContactPerParticle.resize(_numParticle, 0);
	setAllParticlesAsRattler();
#pragma omp parallel num_threads(NUM_THREADS)
	{
		vector<vector2d> threadForce(_numParticle, vector2d(0.0, 0.0));
		vector<double> threadTorque(_numParticle, 0.0);
		vector<int> threadNumContactPerParticle(_numParticle, 0);
		// Todo: compare the performence with coarse grain algorithm
#pragma omp for schedule(dynamic, SIZE_CHUNK)
		for (int p1 = 0; p1 < _numParticle; ++p1) {
			const unordered_map<int, ContactInformation>& contactInfomation = _particleHandlers[p1]->getContactInfo();
			for (auto& iter : contactInfomation) {
				const int& p2 = iter.first;
				const ContactInformation& contactInfo = iter.second;
				threadForce[p1] += contactInfo._forceOnP1;
				threadForce[p2] += contactInfo._forceOnP2;
				threadTorque[p1] += contactInfo._torqueOnP1;
				threadTorque[p2] += contactInfo._torqueOnP2;
				threadNumContactPerParticle[p1]++;
				threadNumContactPerParticle[p2]++;
				_particleHandlers[p1]->setAsNonRattler();
				_particleHandlers[p2]->setAsNonRattler();
			}
		}
#pragma omp critical (BaseWorld_collectForceAndTorque)
		{
			for (int i = 0; i < _numParticle; ++i) {
				_totalForce[i] += threadForce[i];
				_totalTorque[i] += threadTorque[i];
				_numContactPerParticle[i] += threadNumContactPerParticle[i];
			}
		}
	}
}

void BaseWorld::takeTimeIntegral()
{
#pragma omp parallel num_threads(NUM_THREADS)
	{
#pragma omp for schedule(static)
		for (int i = 0; i < _numParticle; ++i) {
			BaseParticle* p = _particleHandlers[i];
			if (p->getSettings() & SETTING::SKIP_TIME_INTEGRAL) continue;
			/*
			2021.04.26 Jiang comment
			the stl iterator begin() function points to the first element,
			this coincide well with the particle ID, which also begins with zero.
			*/
			auto iterForce = _totalForce.begin() + p->getID();
			auto iterTorque = _totalTorque.begin() + p->getID();
			p->takeTimeIntegral(_dt, iterForce, iterTorque);
		}
	}
}

void BaseWorld::resetVelocityOfRattlers()
{
#pragma omp parallel num_threads(NUM_THREADS)
	{
#pragma omp for schedule(static)
		for (int i = 0; i < _numParticle; ++i) {
			if (_particleHandlers[i]->getSettings() & SETTING::SKIP_MACRO_SUM) continue;
			if (_particleHandlers[i]->isRattlerParticle()) {
				_particleHandlers[i]->resetVel();
				_particleHandlers[i]->resetRotVel();
			}
		}
	}
}

void BaseWorld::scaleVelocityOfRattlers(const double scale)
{
#pragma omp parallel num_threads(NUM_THREADS)
	{
#pragma omp for schedule(static)
		for (int i = 0; i < _numParticle; ++i) {
			if (_particleHandlers[i]->getSettings() & SETTING::SKIP_MACRO_SUM) continue;
			if (_particleHandlers[i]->isRattlerParticle()) {
				_particleHandlers[i]->scaleVelocity(scale);
			}
		}
	}
}

void BaseWorld::setAllParticlesAsRattler()
{
#pragma omp parallel num_threads(NUM_THREADS)
	{
#pragma omp for schedule(static)
		for (int i = 0; i < _numParticle; ++i) {
			_particleHandlers[i]->setAsRattler();
		}
	}
}

void BaseWorld::scaleAllParticlePosition(const double& scaleX, const double& scaleY)
{
	for (auto& iter : _particleHandlers) {
		iter->scalePosition(scaleX, scaleY);
	}
}

void BaseWorld::scaleAllParticleVelocity(const double& scale)
{
	for (auto& iter : _particleHandlers) {
		iter->scaleVelocity(scale);
	}
}

void BaseWorld::scaleBoundarySize(const double& scale)
{
	_xmin *= scale;
	_xmax *= scale;
	_ymin *= scale;
	_ymax *= scale;
	_xlen = _xmax - _xmin;
	_ylen = _ymax - _ymin;
}

void BaseWorld::updateAllParticleProperties(const double& density)
{
	for (auto& iter : _particleHandlers) {
		iter->calculateParticleProperties(density);
	}
}

const matrix2d BaseWorld::updateTotalStress()
{
	_totalStress = matrix2d(0.0, 0.0, 0.0, 0.0);
#pragma omp parallel num_threads(NUM_THREADS)
	{
		matrix2d threadStress(0.0, 0.0, 0.0, 0.0);
#pragma omp for schedule(dynamic, SIZE_CHUNK)
		for (int p1 = 0; p1 < _numParticle; ++p1) {
			if (_particleHandlers[p1]->getSettings() & SETTING::SKIP_MACRO_SUM) continue;
			threadStress += _particleHandlers[p1]->getStressFromContactList();
		}
#pragma omp critical (BaseWorld_updateTotalStress)
		{
			_totalStress += threadStress;
		}
	}
	_totalStress /= (_xlen * _ylen);
	return _totalStress;
}

void BaseWorld::setTimeIncrement(const double& dt)
{
	_dt = dt;
}

void BaseWorld::set_OMP_NumThreadsAndSizeChunk(const int& num_threads, const int& size_chunk)
{
	NUM_THREADS = num_threads;
	SIZE_CHUNK = size_chunk;
}

void BaseWorld::updateTotalParticleVolume()
{
	double totalParticleVolume = 0.0;
#pragma omp parallel for num_threads(NUM_THREADS) reduction(+:totalParticleVolume) schedule(static)
	for (int i = 0; i < _numParticle; ++i) {
		if (_particleHandlers[i]->getSettings() & SETTING::SKIP_MACRO_SUM) continue;
		totalParticleVolume += _particleHandlers[i]->getVolume();
	}
	_totalParticleVolume = totalParticleVolume;
	printf("Logger msg: Total particle volume = %lf\n", _totalParticleVolume);
}

void BaseWorld::updateMaximumParticleSize()
{
#pragma omp parallel num_threads(NUM_THREADS) 
	{
		double maxParticleSize = 0.0;
#pragma omp for schedule(static)
		for (int i = 0; i < _numParticle; ++i) {
			if (_particleHandlers[i]->getOuterRadius() > maxParticleSize)
				maxParticleSize = _particleHandlers[i]->getOuterRadius();
		}
#pragma omp critical(BaseWorld_updateMaximumParticleSize)
		{
			if (_maxParticleSize < maxParticleSize)
				_maxParticleSize = maxParticleSize;
		}
	}
	_maxParticleSize *= 2.01;
	printf("Logger msg: Maximum particle size = %lf\n", _maxParticleSize);
}

double BaseWorld::checkTimeIncrement()
{
	double dtMax = 0.0;
	double dtMin = 1e100;
#pragma omp parallel num_threads(NUM_THREADS)
	{
		double dtMax0 = 0.0;
		double dtMin0 = 1e100;
#pragma omp for schedule(static)
		for (int i = 0; i < _numParticle; ++i) {
			for (int j = 0; j < _numParticle; ++j) {
				auto& p1 = _particleHandlers[i];
				auto& p2 = _particleHandlers[j];
				double dtc = _contactModelHandler->checkTimeIncrement(p1->getMass(), p2->getMass());
				if (dtMax0 < dtc) dtMax0 = dtc;
				if (dtMin0 > dtc) dtMin0 = dtc;
			}
		}
#pragma omp critical(BaseWorld_checkTimeIncrement)
		{
			if (dtMax < dtMax0) dtMax = dtMax0;
			if (dtMin > dtMin0) dtMin = dtMin0;
		}
	}
	printf("\n[MESSAGE]\nNumber of collision step checked...\nmin = %d, max = %d\n", 
		static_cast<int>(floor(dtMin / _dt)), static_cast<int>(floor(dtMax / _dt)));

	return dtMin;
}

const int BaseWorld::countPossibleContacts() const
{
	int sum = 0;
	for (int i = 0; i < _possibleContacts.size(); ++i) {
		for (auto pList = _possibleContacts[i].begin(); pList != _possibleContacts[i].end() && *pList != nullptr; ++pList) {
			sum++;
		}
	}
	return sum;
}

const int BaseWorld::countActucalContacts() const
{
	size_t sum = 0;
	for (int i = 0; i < _particleHandlers.size(); ++i) {
		sum += _particleHandlers[i]->getContactInfo().size();
	}
	return static_cast<int>(sum);
}

const double BaseWorld::getKineticEnergy() const
{
	double kineticEnergy = 0.0;
#pragma omp parallel for num_threads(NUM_THREADS) reduction(+:kineticEnergy) schedule(static)
	for (int i = 0; i < _numParticle; ++i) {
		kineticEnergy += _particleHandlers[i]->getKineticEnergy();
	}
	return kineticEnergy;
}

const double BaseWorld::getKineticEnergyPerParticle() const
{
	return getKineticEnergy() / _numParticle;
}

const double BaseWorld::getKineticEnergyPerNonRattlerParticle() const
{
	double kineticEnergy = 0.0;
	int count = 0;
#pragma omp parallel for num_threads(NUM_THREADS) reduction(+:kineticEnergy,count) schedule(static)
	for (int i = 0; i < _numParticle; ++i) {
		if (_particleHandlers[i]->isRattlerParticle()) continue;
		kineticEnergy += _particleHandlers[i]->getKineticEnergy();
		count++;
	}
	return kineticEnergy / count;
}

const double BaseWorld::getSolidFraction() const
{
	return _totalParticleVolume / (_xlen * _ylen);
}

void BaseWorld::writeParticleInformation2Files(const string& path)
{
	stringstream md;
	md << "md " << path.c_str();
	system(md.str().c_str());

	stringstream fileName;
	vector<FILE*> files;
	fileName << path << "\\output_radius.dat";
	FILE* file_radius = fopen(fileName.str().c_str(), "w"); fileName.str("");
	if (!file_radius) {
		printf("can not create radius output file...\n");
		std::cin.get();
	}
	files.push_back(file_radius);
	fileName << path << "\\output_aspectRatio.dat";
	FILE* file_aspectRatio = fopen(fileName.str().c_str(), "w"); fileName.str("");
	files.push_back(file_aspectRatio);
	for (auto it : _particleHandlers) it->print2File_shape(files);
	fclose(file_radius);
	fclose(file_aspectRatio);
}

void BaseWorld::print2File_boundary(FILE* file)
{
	fprintf(file, "%le\t%le\t%le\t%le\n", _xmin, _xmax, _ymin, _ymax);
}

void BaseWorld::print2File_kineticEnergy(FILE* file)
{
	fprintf(file, "%le\n", getKineticEnergyPerParticle());
}

void BaseWorld::print2File_allParticles(FILE* file)
{
	fprintf(file, "%d\n", static_cast<int>(_particleHandlers.size()));
	for (auto& iter : _particleHandlers) {
		iter->BaseParticle::print2File_particle(file);
	}
}

void BaseWorld::closeParticleTimeHistoryFiles()
{
	for (FILE* file : files) {
		fclose(file);
	}
}

void BaseWorld::flushAllFiles()
{
	for (auto& file : files) {
		fflush(file);
	}
}

void BaseWorld::applyBodyForce(const vector2d& bodyForce)
{
#pragma omp parallel num_threads(NUM_THREADS)
	{
#pragma omp for schedule(dynamic, SIZE_CHUNK)
		for (int i = 0; i < _totalForce.size(); ++i) {
			_totalForce[i] += bodyForce;
		}
	}
}

