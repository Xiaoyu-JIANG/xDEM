#pragma once

#include <fstream>
#include "BaseParticle.h"
#include "CircularParticle.h"
#include "ContactGeometry.h"

class CircularParticle;

class EllipticParticle : public BaseParticle
{
protected:
	double _radius = -1.0;
	double _aspectRatio = 1.0;

public:
	EllipticParticle() : BaseParticle() {}

	EllipticParticle(const int& id,
		const double& outerRadius,
		const double& aspectRatio,
		const vector2d& position,
		const vector2d& velocity,
		const double& rotation,
		const double& rotationVel,
		ContactModel* contactModel) :
		BaseParticle(id, outerRadius, position, velocity, rotation, rotationVel, contactModel),
		_radius(outerRadius),
		_aspectRatio(aspectRatio)
	{}

	EllipticParticle(const CircularParticle& p, 
		const double& aspectRatio, 
		ContactModel* contactModel) :
		BaseParticle(p.getID(), p.getOuterRadius(), p.getPos(), p.getVel(), p.getRot(), p.getRotVel(), contactModel),
		_radius(_outerRadius),
		_aspectRatio(aspectRatio)
	{}

	EllipticParticle(std::ifstream& is, const double& aspectRatio, ContactModel* contactModel, const bool& ifReadContact)
	{
		read4File_particle(is, ifReadContact);
		_upLevel = _id;
		_radius = _outerRadius;
		_aspectRatio = aspectRatio;
		_contactModel = contactModel;
		printf("%d-th particle read finish...\n", _id);
	}

	// Over-riden functions
	virtual void calculateParticleProperties(const double& density) override;
	virtual bool calculateContactGeometry(const BaseParticle* other, 
		ContactGeometry& contactGeometry, 
		const vector2d& offset = vector2d(0.0, 0.0),
		const vector2d& offsetVel = vector2d(0.0, 0.0)) override;

	// Get methods
	inline const double getRadius() const { return _radius; }
	inline const double getAspectRatio() const { return _aspectRatio; }
	
	// Set methods
	void scaleRadius(const double& scale) { _radius = _outerRadius * scale; }

	// Other methods
	virtual const bool pointInsideParticle(const vector2d& point) const override;

	// Print methods
	virtual void print2File_shape(vector<FILE*>& files) override;
	virtual void print2File_particle(FILE* file) override;
	void print2File_radius(FILE* file) { fprintf(file, "%lf\n", _radius); }
	void print2File_aspectRatio(FILE* file) { fprintf(file, "%lf\n", _aspectRatio); }


private:
	bool contactWithEllipticParticle(
		const EllipticParticle* other, 
		ContactGeometry& contactGeometry, 
		const vector2d& offset = vector2d(0.0, 0.0),
		const vector2d& offsetVel = vector2d(0.0, 0.0));

	bool contactWithEllipticParticle0(
		const EllipticParticle* other, 
		ContactGeometry& contactGeometry, 
		const vector2d& offset = vector2d(0.0, 0.0),
		const vector2d& offsetVel = vector2d(0.0, 0.0));

	bool solveQuarticGeneral(double& a, double& b, double& c, double& d, double& e, // input parameters
		double& alphaInit, double& betaInit, double& gammaInit, double& deltaInit);
	double fastQuarticSolver(double& a, double& b, double& c, double& d,
		double& alphaInit, double& betaInit, double& gammaInit, double& deltaInit,
		double& alpha, double& beta, double& gamma, double& delta);
	double NewtonRaphsonForQuartic(const double A, const double B, const double C, const double D, const double E, const double init);
};