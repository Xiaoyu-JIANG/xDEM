#include "SegmentObject.h"
#include "CircularParticle.h"

void SegmentObject::calculateParticleProperties(const double& density)
{
	_volume = _radius * 2.0; // linear volume
	_mass = _volume * density;
	_invMass = 1.0 / _mass;
	_inertia = _mass * _radius * _radius / 3.0;
	_invInertia = 1.0 / _inertia;
}

bool SegmentObject::calculateContactGeometry(
	const BaseParticle* otherParticle,
	ContactGeometry& contactGeometry, 
	const vector2d& offset, 
	const vector2d& offsetVel)
{
	if (const CircularParticle* other = dynamic_cast<const CircularParticle*>(otherParticle)) {
		return contactWithCircularParticle(other, contactGeometry, offset, offsetVel);
	}
	else if (const SegmentObject* other = dynamic_cast<const SegmentObject*>(otherParticle)) {
		return 0;
	}
	else {
		printf("Todo: more particle types...\n"); exit(-1);
	}
}

const bool SegmentObject::pointInsideParticle(const vector2d& point) const
{
	return false;
}

void SegmentObject::print2File_shape(vector<FILE*>& files)
{
	print2File_radius(files[0]);
	print2File_aspectRatio(files[1]);
}

void SegmentObject::print2File_particle(FILE* file)
{
}

bool SegmentObject::contactWithCircularParticle(
	const CircularParticle* other, 
	ContactGeometry& contactGeometry, 
	const vector2d& offset, 
	const vector2d& offsetVel)
{
	double projectLength = ((other->getPos() + offset) - _position).dot(_dirVector);
	vector2d refCenter;
	if (projectLength > _radius) {
		refCenter = _vertex2;
	}
	else if (projectLength < -_radius) {
		refCenter = _vertex1;
	}
	else {
		refCenter = _position + _dirVector * projectLength;
	}

	vector2d branchVector = (other->getPos() + offset) - refCenter;
	double distance = branchVector.getLength();
	double sumOfRadius = other->getRadius();
	contactGeometry.normalOverlap = sumOfRadius - distance;
	if (contactGeometry.normalOverlap < 0.0) {
		contactGeometry.isContacted = false;
		return false;
	}
	else {
		contactGeometry.isContacted = true;
		contactGeometry.contactNormal = branchVector.getUnitVector();
		contactGeometry.contactPoint = refCenter;
		contactGeometry.branchVectorFromP1 = contactGeometry.contactPoint - refCenter;
		contactGeometry.branchVectorFromP2 = contactGeometry.contactPoint - (other->getPos() + offset);
		calculateRelativeVelocityAtContactPoint(other, contactGeometry, offsetVel);
		return true;
	}
}
