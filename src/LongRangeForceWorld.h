#pragma once

#include "BaseWorld.h"

class LongRangeForceWorld : public BaseWorld
{
public:
	LongRangeForceWorld() : BaseWorld() {}

	LongRangeForceWorld(const double& dt,
		const vector<BaseParticle*> particleHandlers,
		ContactModel* contactModelHandler) :
		BaseWorld(dt, particleHandlers, contactModelHandler)
	{}

	// Main interfaces in the time loop
	virtual void prepare(const string& path) override;
	virtual void takeTimeStep() override;
	virtual void updateContacts() override;
	virtual void findPossibleContacts() override {}
	virtual void modifyParticlePosition() override {}

	// Print methods
	virtual void print2Screen_worldState(const int& iStep) override;
	virtual void createParticleTimeHistoryFiles(const string& path) override;
	virtual void writeParticleTimeHistory2Files() override;
	virtual void writeParticleTimeHistory2SeparateFiles(const string& path) override;


	
};