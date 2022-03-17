#pragma once

#include "BaseParticle.h"

class GluedParticle : public BaseParticle
{
protected:
	vector<BaseParticle*> _subParticles;
	vector<vector2d> _branchVectors;

public:
	GluedParticle() {}

	GluedParticle(const int& id,
		const vector2d& velocity,
		const double& rotationVel)
	{
		_id = id;
		_upLevel = _id;
		_velocity = velocity;
		_rotationVel = rotationVel;
	}

	// Over-riden functions
	virtual void calculateParticleProperties(const double& density) override;

	virtual bool calculateContactGeometry(
		const BaseParticle* other,
		ContactGeometry& contactGeometry,
		const vector2d& offset = vector2d(0.0, 0.0),
		const vector2d& offsetVel = vector2d(0.0, 0.0)) override;

	virtual void takeTimeIntegral(
		const double& dt, 
		vector<vector2d>::iterator iterForce, 
		vector<double>::iterator iterTorque);

	// Set methods
	void pushBackSubParticle(BaseParticle* subParticle) { _subParticles.push_back(subParticle); }

	// Other methods
	virtual const bool pointInsideParticle(const vector2d& point) const override { return false; }
	virtual const matrix2d getStressFromContactList() const override;
	virtual const bool isRattlerParticle() const override;
	void generateMatlabCodeForCheck() const;

	// Print and read methods
	virtual void print2File_shape(vector<FILE*>& files) override;
	void print2File_radius(FILE* file) { fprintf(file, "%lf\n", 0.0); }
	void print2File_aspectRatio(FILE* file) { fprintf(file, "%lf\n", 1.0); }
};
