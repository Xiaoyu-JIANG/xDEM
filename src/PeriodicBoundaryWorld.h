#pragma once

#include <string>
#include "BaseWorld.h"

using std::string;

class PeriodicBoundaryWorld : public BaseWorld
{
public:
	PeriodicBoundaryWorld() : BaseWorld() {}

	PeriodicBoundaryWorld(const double& dt,
		const vector<BaseParticle*> particleHandlers,
		ContactModel* contactModelHandler) :
		BaseWorld(dt, particleHandlers, contactModelHandler)
	{
		
	}

	// Main interfaces in the time loop
	virtual void prepare(const string& path) override;
	virtual void takeTimeStep() override;
	// Kernel methods: contact detection
	virtual void findPossibleContacts() override;
	virtual void updateContacts() override;
	virtual void modifyParticlePosition() override;

	// Other useful interfaces
	void scaleBoundarySizeToGivenSolidFraction(const double& initSolidFrac);
	void scaleBoundarySizeAndParticlePosition(const double& scale);

	// Get methods

	// Set methods

	// Update methods

	// Print methods
	virtual void print2Screen_worldState(const int& iStep) override;
	virtual void createParticleTimeHistoryFiles(const string& path) override;
	virtual void writeParticleTimeHistory2Files() override;
	virtual void writeParticleTimeHistory2SeparateFiles(const string& path) override;
	void print2File_stress(FILE* file);

	

};