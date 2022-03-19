#include "BaseParticle.h"
#include "MathMatrix2d.h"

#pragma warning(disable : 4996)

#define FORMAT_BaseParticle \
"%7d\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%4d\n"

const double BaseParticle::_GlobalDamping_ = 0.1;

void BaseParticle::calculateRelativeVelocityAtContactPoint(
	const BaseParticle* const other, 
	ContactGeometry& contactGeometry,
	const vector2d& offsetVel)
{
	contactGeometry.relVel =
		(_velocity + contactGeometry.branchVectorFromP1.cross(_rotationVel)) -
		(other->getVel() + offsetVel + contactGeometry.branchVectorFromP2.cross(other->getRotVel()));

	contactGeometry.normalRelVel =
		contactGeometry.contactNormal * contactGeometry.relVel.dot(contactGeometry.contactNormal);

	contactGeometry.shearRelVel =
		contactGeometry.relVel - contactGeometry.normalRelVel;
}

bool BaseParticle::calculateContactForceAndTorque(
	const BaseParticle* const other, 
	const double& dt, 
	const vector2d& offset, 
	const vector2d& offsetVel)
{
	ContactGeometry contactGeometry;
	// If there is a contact...
	if (calculateContactGeometry(other, contactGeometry, offset, offsetVel)) {
		int key = other->getID();
		// If the contact is a new one...
		if (!_contactInfo.count(key)) {
			// Add new contact into undered map.
			_contactInfo.emplace(key,
				ContactInformation(other->getID(), contactGeometry.contactNormal, contactGeometry.contactPoint,
					contactGeometry.branchVectorFromP1, contactGeometry.branchVectorFromP2));
			auto& contactInfo = _contactInfo.at(key);
			// Calculate initial overlap in shear dirction.
			if (contactGeometry.relVel.dot(contactGeometry.contactNormal) == 0.0) {
				contactInfo._shearOverlap = contactGeometry.shearRelVel * dt;
			}
			else {
				contactInfo._shearOverlap = contactGeometry.shearRelVel *
					std::min(dt, contactGeometry.normalOverlap / contactGeometry.relVel.dot(contactGeometry.contactNormal));
			}
			// Calculate damping constant for this contact.
			if (!(_setting & SKIP_TIME_INTEGRAL) && !(other->getSettings() & SKIP_TIME_INTEGRAL)) {
				contactInfo._dampingConstant = _contactModel->calculateDampingConstact(_mass, other->getMass());
			}
			else if (!(_setting & SKIP_TIME_INTEGRAL) && (other->getSettings() & SKIP_TIME_INTEGRAL)) {
				contactInfo._dampingConstant = _contactModel->calculateDampingConstact(_mass);
			}
			else if ((_setting & SKIP_TIME_INTEGRAL) && !(other->getSettings() & SKIP_TIME_INTEGRAL)) {
				contactInfo._dampingConstant = _contactModel->calculateDampingConstact(other->getMass());
			}
			else {
				contactInfo._dampingConstant = _contactModel->calculateDampingConstact(_mass, other->getMass());
			}
		}
		// If the contact has previous data...
		else {
			auto& contactInfo = _contactInfo.at(key);
			// Modify direction of shear overlap.
			double dot = contactInfo._contactNormal.dot(contactGeometry.contactNormal);
			double det = contactInfo._contactNormal.cross(contactGeometry.contactNormal);
			double angle = atan2(det, dot);
			matrix2d R(std::cos(angle), -std::sin(angle),
					   std::sin(angle),  std::cos(angle));
			contactInfo._shearOverlap = R.dot(contactInfo._shearOverlap);
			contactInfo._shearOverlap += contactGeometry.shearRelVel * dt;
			contactInfo._contactNormal = contactGeometry.contactNormal;
			contactInfo._contactPoint = contactGeometry.contactPoint;
			contactInfo._branchVectorFromP1 = contactGeometry.branchVectorFromP1;
			contactInfo._branchVectorFromP2 = contactGeometry.branchVectorFromP2;
			//std::cout << contactInfo._shearOverlap << "  " << _position << std::endl;
		}
		// Calculate normal and shear forces.
		auto& contactInfo = _contactInfo.at(key);
		contactInfo._branchVector = _position - (other->getPos() + offset);
		_contactModel->calculateNormalForce(contactGeometry, contactInfo);
		_contactModel->calculateShearForce(contactGeometry, contactInfo);
		contactInfo.calculateTotalForceAndTorque();
		return true;
		//contactInfo.print();
	}
	// If there is no contact...
	else {
		// Delete old contact information.
		_contactInfo.erase(other->getID());
		return false;
	}
}

/*
2021.04.26 Jiang comment
Instead of passing references of force and torque into this function,
a std::vector iterator is passed, this may cause a little bit performance lose.
But for implement of glued particle model, 
this change will well keep the uniformity of the time integral.
*/
void BaseParticle::takeTimeIntegral(
	const double& dt, 
	vector<vector2d>::iterator iterForce, 
	vector<double>::iterator iterTorque)
{
	_velocity = (_velocity * (1.0 - _GlobalDamping_ * dt) + (*iterForce * _invMass) * dt)
		/ (1.0 + _GlobalDamping_ * dt);
	_position += _velocity * dt;
	_rotationVel = (_rotationVel * (1.0 - _GlobalDamping_ * dt) + *iterTorque * _invInertia * dt) 
		/ (1.0 + _GlobalDamping_ * dt);
	_rotation += _rotationVel * dt;
}

void BaseParticle::takeTimeIntegral(
	const double& dt, 
	const double& velXcontrol, 
	vector<vector2d>::iterator iterForce, 
	vector<double>::iterator iterTorque)
{
	_velocity = (_velocity * (1.0 - _GlobalDamping_ * dt) + *iterForce * _invMass * dt)
		/ (1.0 + _GlobalDamping_ * dt);
	_position.moveX(velXcontrol * dt);
	_position.moveY(_velocity.y() * dt);
	_rotationVel = (_rotationVel * (1.0 - _GlobalDamping_ * dt) + *iterTorque * _invInertia * dt)
		/ (1.0 + _GlobalDamping_ * dt);
	_rotation += _rotationVel * dt;
}

//void BaseParticle::takeTimeIntegral(const double& dt, const vector2d& force, const double& torque)
//{
//	_velocity += force * _invMass * dt;
//	_position += _velocity * dt;
//	_rotationVel += torque * _invInertia * dt;
//	_rotation += _rotationVel * dt;
//}

const matrix2d BaseParticle::getStressFromContactList() const
{
	matrix2d stress(0.0, 0.0, 0.0, 0.0);
	for (auto& iter : _contactInfo) {
		const ContactInformation& contactInfo = iter.second;
		stress += matrix2d::dyadic(contactInfo._force, contactInfo._branchVector);
	}
	return stress;
}

const bool BaseParticle::isRattlerParticle() const
{
	return _isRattler;
}

void BaseParticle::print2File_particle(FILE* file)
{
	fprintf(file, FORMAT_BaseParticle,
		_id, _outerRadius,
		_position.x(), _position.y(),
		_velocity.x(), _velocity.y(),
		_rotation, _rotationVel,
		_volume, _mass, _invMass, _inertia, _invInertia,
		static_cast<int>(_contactInfo.size()));
	for (auto& iter : _contactInfo) {
		iter.second.print2File_contactInfo(file);
	}
}

void BaseParticle::read4File_particle(std::ifstream& is, const bool& ifReadContact)
{
	int numContact;
	is >> _id >> _outerRadius 
		>> _position >> _velocity
		>> _rotation >> _rotationVel
		>> _volume >> _mass >> _invMass >> _inertia >> _invInertia
		>> numContact;
	for (int i = 0; i < numContact; ++i) {
		ContactInformation contactInfo(is);
		if (ifReadContact) {
			_contactInfo.emplace(contactInfo._otherID, contactInfo);
		}
	}
}

void BaseParticle::print2File_contact(FILE* file)
{
	for (auto& iter : _contactInfo) {
		auto& contact = iter.second;
		fprintf(file, "%d\t%d\t%le\t%le\t%le\t%le\t%le\t%le\n",
			_id, contact._otherID,
			contact._contactPoint.x(), contact._contactPoint.y(),
			contact._normalForce.x(), contact._normalForce.y(),
			contact._shearForce.x(), contact._shearForce.y());
	}
}
