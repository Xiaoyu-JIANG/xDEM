#pragma once

#include <fstream>
#include "BaseParticle.h"
#include "ContactGeometry.h"
#include "SegmentObject.h"

class SegmentObject;

class CircularParticle : public BaseParticle
{
protected:
	double _radius = -1.0;

public:
	CircularParticle() : BaseParticle() {}

	CircularParticle(const int& id,
		const double& outerRadius,
		const vector2d& position,
		const vector2d& velocity,
		const double& rotation,
		const double& rotationVel,
		ContactModel* contactModel) :
		BaseParticle(id, outerRadius, position, velocity, rotation, rotationVel, contactModel),
		_radius(outerRadius)
	{}

	CircularParticle (std::ifstream& is, const bool& ifReadContacts);



	// Over-riden functions
	virtual void calculateParticleProperties(const double& density) override;
	virtual bool calculateContactGeometry(
		const BaseParticle* other,
		ContactGeometry& contactGeometry,
		const vector2d& offset = vector2d(0.0, 0.0),
		const vector2d& offsetVel = vector2d(0.0, 0.0)) override;

	// Get methods
	inline const double getRadius() const { return _radius; }

	// Set methods
	inline void setRadius(const double& radius) { _radius = radius; }
	inline void scaleRadius(const double& scale) { _radius *= scale; }

	// Other methods
	virtual const bool pointInsideParticle(const vector2d& point) const override;

	// Print methods
	virtual void print2File_shape(vector<FILE*>& files) override;
	virtual void print2File_particle(FILE* file) override;
	void print2File_radius(FILE* file) { fprintf(file, "%lf\n", _radius); }
	void print2File_aspectRatio(FILE* file) { fprintf(file, "%lf\n", 1.0); }

private:
	bool contactWithCircularParticle(
		const CircularParticle* other,
		ContactGeometry& contactGeometry,
		const vector2d& offset = vector2d(0.0, 0.0),
		const vector2d& offsetVel = vector2d(0.0, 0.0));

	bool contactWithSegmentObject(
		const SegmentObject* other,
		ContactGeometry& contactGeometry,
		const vector2d& offset = vector2d(0.0, 0.0),
		const vector2d& offsetVel = vector2d(0.0, 0.0));

};
