#include <sstream>
#include <ctime>
#include "LongRangeForceWorld.h"

using std::stringstream;

void LongRangeForceWorld::prepare(const string& path)
{
	stringstream md;
	md << "md " << path.c_str();
	system(md.str().c_str());
	createParticleTimeHistoryFiles(path);
	writeParticleInformation2Files(path);
}

void LongRangeForceWorld::takeTimeStep()
{
	updateContacts();
	collectForceAndTorque();
	takeTimeIntegral();
}

void LongRangeForceWorld::createParticleTimeHistoryFiles(const string& path)
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

	fileName << path << "\\output_kineticEnergy.dat";
	FILE* file_kineticEnergy = fopen(fileName.str().c_str(), "w");
	files.push_back(file_kineticEnergy);
	fileName.str("");
}

void LongRangeForceWorld::writeParticleTimeHistory2Files()
{
	for (auto it : _particleHandlers) it->print2File_velocity(files[0]);
	fprintf(files[0], "\n");

	for (auto it : _particleHandlers) it->print2File_position(files[1]);
	fprintf(files[1], "\n");

	print2File_kineticEnergy(files[2]);
}

void LongRangeForceWorld::writeParticleTimeHistory2SeparateFiles(const string& path)
{
	
}

void LongRangeForceWorld::print2Screen_worldState(const int& iStep)
{
	static clock_t timerTotal = clock();
	static clock_t timer = clock();
	timer = clock() - timer;
	printf("Current step = %d, Time used per 1000 steps = %lf, Total time used = %lf\n", iStep,
		(static_cast<double>(timer) / CLOCKS_PER_SEC),
		((static_cast<double>(clock()) - timerTotal) / CLOCKS_PER_SEC));
	timer = clock();
}

void LongRangeForceWorld::updateContacts()
{
#pragma omp parallel num_threads(NUM_THREADS)
	{
#pragma omp for schedule(dynamic, SIZE_CHUNK)
		for (int i = 0; i < _numParticle; ++i) {
			BaseParticle* ip = _particleHandlers[i];
			for (int j = 0; j < _numParticle; ++j) {
				if (i <= j) continue;
				BaseParticle* jp = _particleHandlers[j];
				ip->calculateContactForceAndTorque(jp, _dt);
			}
		}
	}
}
