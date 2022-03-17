#pragma once

#include "PeriodicBoundaryWorld.h"

class PeriodicControlWorld : public PeriodicBoundaryWorld
{
protected:
	double _mass = 0.0;
	double _velX = 0.0;
	double _velY = 0.0;
	double _maxStrainRate = 0.0;

public:
	PeriodicControlWorld() {}

	PeriodicControlWorld(const double& dt,
		const vector<BaseParticle*> particleHandlers,
		ContactModel* contactModelHandler) :
		PeriodicBoundaryWorld(dt, particleHandlers, contactModelHandler)
	{}

	// Main interfaces in the time loop
	virtual void prepare(const string& path) override;
	virtual void prepare2(const string& path);
	virtual void takeTimeStep() override;

	// Get methods

	// Set methods
	void setVirtualBoundaryMass();
	void setVirtualBoundaryMass(const double& mass);
	void setMaxStrainRate(const double maxStrainRate);

	// Update methods
	void updatePeriodicBoundary_strainControl(const double& strainRateX, const double& strainRateY);
	void updatePeriodicBoundary_stressControlUniform(const double& stress);
	void updatePeriodicBoundary_stressControlUniform_strainLimit(const double& stress);
	void updatePeriodicBoundary_stressControlX(const double& stress);
	void updatePeriodicBoundary_stressControlY(const double& stress);

	// Print methods
	virtual void print2Screen_worldState(const int& iStep) override;

};