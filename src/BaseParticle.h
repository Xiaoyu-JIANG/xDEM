#pragma once

// Todo list:
// 1. Add displacement control and force control interfaces...
//		remove control method from world to BaseParticle class member function...

#include <vector>
#include <unordered_map>
#include <iostream>
#include "MathVector.h"
#include "MathMatrix2d.h"
#include "ContactGeometry.h"
#include "ContactInformation.h"
#include "ContactModelList.h"

using std::unordered_map;
using std::vector;

enum SETTING {
	NONE = 0,
	SKIP_TIME_INTEGRAL = 1 << 0,
	SKIP_CONTACT_DETECTION = 1 << 1,
	SKIP_MACRO_SUM = 1 << 2
};
constexpr SETTING operator|(SETTING lhs, SETTING rhs) { return static_cast<SETTING>(static_cast<int>(lhs) | static_cast<int>(rhs)); }
constexpr SETTING operator&(SETTING lhs, SETTING rhs) { return static_cast<SETTING>(static_cast<int>(lhs) & static_cast<int>(rhs)); }
constexpr SETTING& operator|=(SETTING& lhs, SETTING rhs) { lhs = lhs | rhs; return lhs; }


class BaseParticle
{
protected:
	const static double _GlobalDamping_;
	int			_upLevel = -1;
	int			_id = -1;
	bool		_isRattler = false;
	double		_outerRadius = -1.0;
	vector2d	_position;
	vector2d	_velocity;
	double		_rotation = 0.0;
	double		_rotationVel = 0.0;
	double		_volume = -1.0;
	double		_mass = -1.0;
	double		_invMass = -1.0;
	double		_inertia = -1.0;
	double		_invInertia = -1.0;
	unordered_map<int, ContactInformation> _contactInfo;
	ContactModel* _contactModel = nullptr;
	SETTING		_setting = SETTING::NONE;

public:
	BaseParticle() {}

	BaseParticle(const int& id,
		const double& outerRadius,
		const vector2d& position,
		const vector2d& velocity,
		const double& rotation,
		const double& rotationVel,
		ContactModel* contactModel) :
		_id(id),
		_outerRadius(outerRadius),
		_position(position),
		_velocity(velocity),
		_rotation(rotation),
		_rotationVel(rotationVel),
		_contactModel(contactModel)
	{
		_upLevel = _id;
		//_contactInfo.reserve(30);
	}

	~BaseParticle() {}

	// Virtual functions
	virtual void calculateParticleProperties(
		const double& density) = 0;

	virtual bool calculateContactGeometry(
		const BaseParticle* other,
		ContactGeometry& contactGeometry,
		const vector2d& offset = vector2d(0.0, 0.0),
		const vector2d& offsetVel = vector2d(0.0, 0.0)) = 0;

	void calculateRelativeVelocityAtContactPoint(
		const BaseParticle* const other,
		ContactGeometry& contactGeometry,
		const vector2d& offsetVel = vector2d(0.0, 0.0));

	virtual bool calculateContactForceAndTorque(
		const BaseParticle* const other,
		const double& dt,
		const vector2d& offset = vector2d(0.0, 0.0),
		const vector2d& offsetVel = vector2d(0.0, 0.0));

	//virtual void takeTimeIntegral(const double& dt, const vector2d& force, const double& torque);
	virtual void takeTimeIntegral(
		const double& dt,
		vector<vector2d>::iterator iterForce,
		vector<double>::iterator iterTorque);
	virtual void takeTimeIntegral(
		const double& dt,
		const double& velXcontrol,
		vector<vector2d>::iterator iterForce,
		vector<double>::iterator iterTorque);

	// get methods
	inline const int& getUpLevel() const { return _upLevel; }
	inline const int& getID() const { return _id; }
	inline const double& getOuterRadius() const { return _outerRadius; }
	inline const vector2d& getPos() const { return _position; }
	inline const double getPosX() const { return _position.x(); }
	inline const double getPosY() const { return _position.y(); }
	inline const vector2d& getVel() const { return _velocity; }
	inline const double getVelX() const { return _velocity.x(); }
	inline const double getVelY() const { return _velocity.y(); }
	inline const double& getRot() const { return _rotation; }
	inline const double& getRotVel() const { return _rotationVel; }
	inline const double& getVolume() const { return _volume; }
	inline const double& getMass() const { return _mass; }
	inline const double& getInvMass() const { return _invMass; }
	inline const double& getInertia() const { return _inertia; }
	inline const double& getInvInertia() const { return _invInertia; }
	inline const unordered_map<int, ContactInformation> getContactInfo() const { return _contactInfo; }
	inline const SETTING& getSettings() const { return _setting; }
	inline const double getKineticEnergy() const { return 0.5 * (_mass * _velocity.getLengthSquared() + _inertia * _rotationVel * _rotationVel); }

	inline const bool takeCoarseDetection(const BaseParticle* const other) const {
		return (_outerRadius + other->getOuterRadius()) > _position.getDistance(other->getPos()); }
	inline const bool takeCoarseDetection(const BaseParticle* const other, const vector2d& offset) const {
		return (_outerRadius + other->getOuterRadius()) > _position.getDistance(other->getPos() + offset); }

	// set methods
	inline void setUpLevel(const int& level) { _upLevel = level; }
	inline void setMass(const double& mass) { _mass = mass; _invMass = 1.0 / _mass; }
	inline void setInertia(const double& inertia) { _inertia = inertia; _invInertia = 1.0 / _inertia; }
	inline void setOuterRadius(const double& outerRadius) { _outerRadius = outerRadius; }
	inline void scaleOuterRadius(const double& scale) { _outerRadius *= scale; }
	inline void setPos(const vector2d& pos) { _position = pos; }
	inline void shiftPosX(const double& offset) { _position.moveX(offset); }
	inline void shiftPosY(const double& offset) { _position.moveY(offset); }
	inline void setRot(const double& rot) { _rotation = rot; }
	inline void accumulateRot(const double& value) { _rotation += value; }
	inline void resetVel() { _velocity = vector2d(0.0, 0.0); }
	inline void resetRotVel() { _rotationVel = 0.0; }
	inline void scalePosition(const double& scale) { _position *= scale; }
	inline void scalePosition(const double& scaleX, const double& scaleY) { _position.scaleX(scaleX); _position.scaleY(scaleY); }
	inline void scaleVelocity(const double& scale) { _velocity *= scale; _rotationVel *= scale; }
	inline void reverseVelX() { _velocity.setX(-_velocity.x()); }
	inline void setVel(const vector2d& velocity) { _velocity = velocity; }
	inline void setRotVel(const double& rotVel) { _rotationVel = rotVel; }
	inline void setVelX(const double& value) { _velocity.setX(value); }
	inline void setVelY(const double& value) { _velocity.setY(value); }
	inline void appendSetting(const SETTING& setting) { _setting |= setting; }
	inline void clearSetting() { _setting = SETTING::NONE; }
	inline void setAsRattler() { _isRattler = true; }
	inline void setAsNonRattler() { _isRattler = false; }

	// Other methods
	virtual const bool pointInsideParticle(const vector2d& point) const = 0;
	virtual const matrix2d getStressFromContactList() const;
	virtual const bool isRattlerParticle() const;

	// Print and read methods
	virtual void print2File_shape(vector<FILE*>& files) = 0;
	virtual void print2File_particle(FILE* file);
	void read4File_particle(std::ifstream& is, const bool& ifReadContact);
	void print2File_position(FILE* file) { fprintf(file, "%le\t%le\t", _position.x(), _position.y()); }
	void print2File_velocity(FILE* file) { fprintf(file, "%le\t%le\t", _velocity.x(), _velocity.y()); }
	void print2File_rotation(FILE* file) { fprintf(file, "%le\t", _rotation); }
	void print2File_rotationVelocity(FILE* file) { fprintf(file, "%le\t", _rotationVel); }
	void print2File_contact(FILE* file);
};
