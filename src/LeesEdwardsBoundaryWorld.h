#pragma once

#include "PeriodicControlWorld.h"

class LeesEdwardsBoundaryWorld : public PeriodicControlWorld
{
protected:
	double _shiftStrainRate = 0.0;
	double _shiftVel = 0.0;
	double _totalShift = 0.0;
	double _totalStrain = 0.0;
	vector<array<BaseParticle*, MAXIMUM_POSSIBLE_CONTACT_LIST_SIZE>> _possibleContactsLE_top;
	vector<array<BaseParticle*, MAXIMUM_POSSIBLE_CONTACT_LIST_SIZE>> _possibleContactsLE_bottom;

public:
	LeesEdwardsBoundaryWorld() {}

	LeesEdwardsBoundaryWorld(const double& dt,
		const vector<BaseParticle*> particleHandlers,
		ContactModel* contactModelHandler) :
		PeriodicControlWorld(dt, particleHandlers, contactModelHandler)
	{
		array<BaseParticle*, MAXIMUM_POSSIBLE_CONTACT_LIST_SIZE> initialArray; initialArray.fill(nullptr);
		_possibleContactsLE_top.resize(_particleHandlers.size(), initialArray);
		_possibleContactsLE_bottom.resize(_particleHandlers.size(), initialArray);
	}

	// Main interfaces in the time loop
	virtual void prepare(const string& path) override;
	virtual void takeTimeStep() override;
	// Kernel methods: contact detection
	virtual void findPossibleContacts() override;
	virtual void updateContacts() override;
	virtual void modifyParticlePosition() override;
	// Special kernel methods
	void updateLeesEdwardBoundary();

	// Get methods
	const double getTotalShift() const { return _totalShift; }
	const double getTotalStrain() const { return _totalStrain; }

	// Set methods
	void setShiftStrainRate(const double& shiftStrainRate);
	void setParticleVelocityAccordingToPosition();
	void moveParticlePositionAccordingToPosition();
	void applyBackgroudVelocity();
	void removeBackgroudVelocity();

	// Update methods

	// Print methods
	virtual void print2Screen_worldState(const int& iStep) override;
	virtual void createParticleTimeHistoryFiles(const string& path) override;
	virtual void writeParticleTimeHistory2Files() override;
	virtual void writeParticleTimeHistory2SeparateFiles(const string& path) override;
	void print2File_boundaryShift(FILE* file);


	
};