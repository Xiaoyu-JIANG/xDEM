#pragma once

#include "BaseParticle.h"

class PointParticle : public BaseParticle
{
public:
	PointParticle() : BaseParticle() {}

	PointParticle(const int& id,
		const vector2d& position,
		const vector2d& velocity,
		ContactModel* contactModel) :
		BaseParticle(id, 0.0, position, velocity, 0.0, 0.0, contactModel)
	{}

	virtual void calculateParticleProperties(
		const double& density) override {}

	virtual bool calculateContactGeometry(
		const BaseParticle* other,
		ContactGeometry& contactGeometry,
		const vector2d& offset = vector2d(0.0, 0.0),
		const vector2d& offsetVel = vector2d(0.0, 0.0)) override;

	virtual bool calculateContactForceAndTorque(
		const BaseParticle* const other,
		const double& dt,
		const vector2d& offset = vector2d(0.0, 0.0),
		const vector2d& offsetVel = vector2d(0.0, 0.0)) override;


	virtual const bool pointInsideParticle(const vector2d& point) const override { return false; }


	// Print and read methods
	virtual void print2File_shape(vector<FILE*>& files) override {}
private:
	bool contactWithPointParticle(
		const PointParticle* other,
		ContactGeometry& contactGeometry,
		const vector2d& offset = vector2d(0.0, 0.0),
		const vector2d& offsetVel = vector2d(0.0, 0.0));
};