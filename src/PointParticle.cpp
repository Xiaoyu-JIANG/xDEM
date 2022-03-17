#include "PointParticle.h"

bool PointParticle::calculateContactGeometry(const BaseParticle* otherParticle, ContactGeometry& contactGeometry, const vector2d& offset, const vector2d& offsetVel)
{
	if (const PointParticle* other = dynamic_cast<const PointParticle*>(otherParticle)) {
		return contactWithPointParticle(other, contactGeometry, offset, offsetVel);
	}
	else {
		printf("Todo: more particle types...\n"); exit(-1);
	}
}

bool PointParticle::calculateContactForceAndTorque(
	const BaseParticle* const other, 
	const double& dt, 
	const vector2d& offset, 
	const vector2d& offsetVel)
{
	ContactGeometry contactGeometry;
	// If there is a contact...
	calculateContactGeometry(other, contactGeometry, offset, offsetVel);
	int key = other->getID();
	if (!_contactInfo.count(key)) {
		_contactInfo.emplace(key,
			ContactInformation(other->getID(), contactGeometry.contactNormal, 
				vector2d(_mass, other->getMass()), // The class member _contactPoint is used to save mass information...
				vector2d(0.0, 0.0), vector2d(0.0, 0.0)));
	}
	else {
		auto& contactInfo = _contactInfo.at(key);
		contactInfo._contactNormal = contactGeometry.contactNormal;
	}
	auto& contactInfo = _contactInfo.at(key);
	_contactModel->calculateNormalForce(contactGeometry, contactInfo);
	contactInfo.calculateTotalForceAndTorque();
	return true;
}

bool PointParticle::contactWithPointParticle(const PointParticle* other, ContactGeometry& contactGeometry, const vector2d& offset, const vector2d& offsetVel)
{
	vector2d branchVector = (other->getPos() + offset) - _position;
	contactGeometry.normalOverlap = branchVector.getLength();
	contactGeometry.contactNormal = branchVector.getUnitVector();
	//printf("[DEBUG] Function called...PointParticle::contactWithPointParticle\n");
	//std::cout << branchVector << std::endl;
	//std::cout << contactGeometry.normalOverlap << std::endl;
	return true;
}
