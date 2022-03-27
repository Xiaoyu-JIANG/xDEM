#pragma once

#include <vector>
#include <array>
#include <string>

#include "MathVector.h"
#include "MathMatrix2d.h"
#include "BaseParticle.h"

#define MAXIMUM_POSSIBLE_CONTACT_LIST_SIZE 200

using std::vector;
using std::array;
using std::string;

class BaseWorld
{
protected:
	int NUM_THREADS = 1;
	int SIZE_CHUNK = 100;
	int _numParticle = 0;
	double _dt = 0.0;
	double _xmin = 0.0;
	double _xmax = 0.0;
	double _ymin = 0.0;
	double _ymax = 0.0;
	double _xlen = 0.0;
	double _ylen = 0.0;
	vector<BaseParticle*> _particleHandlers;
	ContactModel* _contactModelHandler = nullptr;
	vector<array<BaseParticle*, MAXIMUM_POSSIBLE_CONTACT_LIST_SIZE>> _possibleContacts;
	vector<vector2d> _totalForce;
	vector<double> _totalTorque;
	vector<int> _numContactPerParticle;
	matrix2d _totalStress;
	double _totalParticleVolume = 0.0;
	double _maxParticleSize = 0.0;
	vector<FILE*> files;

public:
	BaseWorld() {}

	BaseWorld(const double& dt,
		const vector<BaseParticle*> particleHandlers,
		ContactModel* contactModelHandler) :
		_dt(dt),
		_particleHandlers(particleHandlers),
		_contactModelHandler(contactModelHandler)
	{
		_numParticle = static_cast<int>(_particleHandlers.size());
		array<BaseParticle*, MAXIMUM_POSSIBLE_CONTACT_LIST_SIZE> initialArray; initialArray.fill(nullptr);
		_possibleContacts.resize(_particleHandlers.size(), initialArray);
		updateMaximumParticleSize();
		updateTotalParticleVolume();
	}

	// Main interfaces in the time loop
	virtual void prepare(const string& path) = 0;
	virtual void takeTimeStep() = 0;
	double checkTimeIncrement();
	void collectForceAndTorque();
	void takeTimeIntegral();
	void closeParticleTimeHistoryFiles();
	void flushAllFiles();
	void applyBodyForce(const vector2d& bodyForce);

	// Kernel methods: contact detection
	virtual void findPossibleContacts() = 0;
	virtual void updateContacts() = 0;
	virtual void modifyParticlePosition() = 0;
	

	// Get methods
	const double getXmin() const { return _xmin; }
	const double getXmax() const { return _xmax; }
	const double getYmin() const { return _ymin; }
	const double getYmax() const { return _ymax; }
	const double getXlen() const { return _xlen; }
	const double getYlen() const { return _ylen; }
	const int countActucalContacts() const;
	const double getElasticEnergy() const;
	const double getElasticEnergyPerContact() const;
	const int countPossibleContacts() const;
	const double getKineticEnergy() const;
	const double getKineticEnergyPerParticle() const;
	const double getKineticEnergyPerNonRattlerParticle() const;
	const double getSolidFraction() const;
	const matrix2d getStress() const { return _totalStress; };

	// Set methods
	void setTimeIncrement(const double& dt);
	void set_OMP_NumThreadsAndSizeChunk(const int& num_threads, const int& size_chunk);
	void setWorldBoundary(const double& xmin, const double& xmax, const double& ymin, const double& ymax);
	void resetVelocityOfRattlers();
	void scaleVelocityOfRattlers(const double scale);
	void setAllParticlesAsRattler();
	void scaleAllParticlePosition(const double& scaleX, const double& scaleY);
	void scaleAllParticleVelocity(const double& scale);
	void resetAllParticleVelocity();
	void scaleBoundarySize(const double& scale);

	// Update methods
	void updateAllParticleProperties(const double& density);
	const matrix2d updateTotalStress();
	void updateTotalParticleVolume();
	void updateMaximumParticleSize();

	// Print methods
	virtual void print2Screen_worldState(const int& iStep) = 0;
	virtual void createParticleTimeHistoryFiles(const string& path) = 0;
	virtual void writeParticleTimeHistory2Files() = 0;
	virtual void writeParticleTimeHistory2SeparateFiles(const string& path) = 0;
	void writeParticleInformation2Files(const string& path);
	void print2File_boundary(FILE* file);
	void print2File_kineticEnergy(FILE* file);
	void print2File_allParticles(FILE* file);


	
	
};