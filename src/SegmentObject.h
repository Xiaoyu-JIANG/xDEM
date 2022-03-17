#pragma once

#include <fstream>
#include "BaseParticle.h"
#include "ContactGeometry.h"
#include "CircularParticle.h"


class CircularParticle;

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// This segment object is only for fixed object
class SegmentObject : public BaseParticle
{
private:
	double _radius = -1.0;
	vector2d _vertex1;
	vector2d _vertex2;
	vector2d _dirVector;

public:
	SegmentObject() {}

	SegmentObject(const int& id,
		const vector2d& vertex1,
		const vector2d& vertex2,
		const vector2d& velocity,
		const double& rotationVel,
		ContactModel* contactModel) :
		BaseParticle(id, 0.0, vector2d(0.0, 0.0), velocity, 0.0, rotationVel, contactModel),
		_vertex1(vertex1),
		_vertex2(vertex2)
	{
		_outerRadius = _vertex1.getDistance(_vertex2) / 2.0;
		_radius = _outerRadius;
		_position = (_vertex1 + _vertex2) / 2.0;
		vector2d branchVector = _vertex2 - _vertex1;
		_rotation = atan2(branchVector.y(), branchVector.x());
		_dirVector = branchVector.getUnitVector();
		appendSetting(SETTING::SKIP_TIME_INTEGRAL);
	}

	// Over-riden functions
	virtual void calculateParticleProperties(const double& density) override;
	virtual bool calculateContactGeometry(
		const BaseParticle* other,
		ContactGeometry& contactGeometry,
		const vector2d& offset = vector2d(0.0, 0.0),
		const vector2d& offsetVel = vector2d(0.0, 0.0)) override;

	// Get methods
	inline const double getRadius() const { return _radius; }
	inline const vector2d getVertex1() const { return _vertex1; }
	inline const vector2d getVertex2() const { return _vertex2; }
	inline const vector2d getDirVector() const { return _dirVector; }

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
};