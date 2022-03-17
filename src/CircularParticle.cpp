#include <ctime>
#include <random>
#define _USE_MATH_DEFINES
#include <math.h>
#include "CircularParticle.h"

#define FORMAT_CircularParticle \
"%7d\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%4d\n"

CircularParticle::CircularParticle(std::ifstream& is, const bool& ifReadContact)
{
	read4File_particle(is, ifReadContact);
	_radius = _outerRadius;
}

void CircularParticle::calculateParticleProperties(const double& density)
{
	_volume = M_PI * _radius * _radius;
	_mass = _volume * density;
	_invMass = 1.0 / _mass;
	_inertia = _mass * _radius * _radius / 2.0;
	_invInertia = 1.0 / _inertia;
}

bool CircularParticle::calculateContactGeometry(
	const BaseParticle* otherParticle, 
	ContactGeometry& contactGeometry, 
	const vector2d& offset, 
	const vector2d& offsetVel)
{
	if (const CircularParticle* other = dynamic_cast<const CircularParticle*>(otherParticle)) {
		return contactWithCircularParticle(other, contactGeometry, offset, offsetVel);
	}
	else if (const SegmentObject* other = dynamic_cast<const SegmentObject*>(otherParticle)) {
		return contactWithSegmentObject(other, contactGeometry, offset, offsetVel);
	}
	else {
		printf("Todo: more particle types...\n"); exit(-1);
	}
}

const bool CircularParticle::pointInsideParticle(const vector2d& point) const
{
	return _position.getDistance(point) < _radius;
}

void CircularParticle::print2File_shape(vector<FILE*>& files)
{
	print2File_radius(files[0]);
	print2File_aspectRatio(files[1]);
}

void CircularParticle::print2File_particle(FILE* file)
{
	fprintf(file, FORMAT_CircularParticle,
		_id, _outerRadius, _radius,
		_position.x(), _position.y(),
		_velocity.x(), _velocity.y(),
		_rotation, _rotationVel,
		_volume, _mass, _invMass, _inertia, _invInertia,
		static_cast<int>(_contactInfo.size()));
	for (auto& iter : _contactInfo) {
		iter.second.print2File_contactInfo(file);
	}
}

bool CircularParticle::contactWithCircularParticle(
	const CircularParticle* other, 
	ContactGeometry& contactGeometry, 
	const vector2d& offset, 
	const vector2d& offsetVel)
{
	vector2d branchVector = (other->getPos() + offset) - _position;
	double distance = branchVector.getLength();
	double sumOfRadius = _radius + other->getRadius();
	contactGeometry.normalOverlap = sumOfRadius - distance;
	if (contactGeometry.normalOverlap < 0.0) {
		contactGeometry.isContacted = false;
		return false;
	}
	else {
		contactGeometry.isContacted = true;
		contactGeometry.contactNormal = branchVector.getUnitVector();
		vector2d p1 = _position + contactGeometry.contactNormal * _radius;
		vector2d p2 = (other->getPos() + offset) - contactGeometry.contactNormal * other->getRadius();
		contactGeometry.contactPoint = (p1 + p2) / 2.0;
		contactGeometry.branchVectorFromP1 = contactGeometry.contactPoint - _position;
		contactGeometry.branchVectorFromP2 = contactGeometry.contactPoint - (other->getPos() + offset);
		calculateRelativeVelocityAtContactPoint(other, contactGeometry, offsetVel);
		return true;
	}
}

bool CircularParticle::contactWithSegmentObject(
	const SegmentObject* other, 
	ContactGeometry& contactGeometry, 
	const vector2d& offset, 
	const vector2d& offsetVel)
{
	double projectLength = ((_position + offset) - other->getPos()).dot(other->getDirVector());
	vector2d refCenter;
	if (projectLength > other->getRadius()) {
		refCenter = other->getVertex2();
	}
	else if (projectLength < -other->getRadius()) {
		refCenter = other->getVertex1();
	}
	else {
		refCenter = other->getPos() + other->getDirVector() * projectLength;
	}

	vector2d branchVector = (refCenter + offset) - _position;
	double distance = branchVector.getLength();
	double sumOfRadius = _radius;
	contactGeometry.normalOverlap = sumOfRadius - distance;
	if (contactGeometry.normalOverlap < 0.0) {
		contactGeometry.isContacted = false;
		return false;
	}
	else {
		contactGeometry.isContacted = true;
		contactGeometry.contactNormal = branchVector.getUnitVector();
		contactGeometry.contactPoint = refCenter;
		contactGeometry.branchVectorFromP1 = contactGeometry.contactPoint - _position;
		contactGeometry.branchVectorFromP2 = contactGeometry.contactPoint - (refCenter + offset);
		calculateRelativeVelocityAtContactPoint(other, contactGeometry, offsetVel);
		return true;
	}
}

