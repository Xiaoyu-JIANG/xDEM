#pragma once

#include <fstream>
#include "MathVector.h"

#pragma warning(disable : 4996)
#define FORMAT_ContactInformation_print2File_contactInfo \
"\t%7d\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t\
%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t\
%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\n"

struct ContactInformation
{
public:
	int _otherID = -1;
	vector2d _contactNormal;
	vector2d _contactPoint;
	vector2d _branchVector;
	vector2d _branchVectorFromP1;
	vector2d _branchVectorFromP2;
	vector2d _shearOverlap;
	vector2d _normalForce;
	vector2d _shearForce;
	vector2d _force;
	vector2d _forceOnP1;
	vector2d _forceOnP2;
	double _torqueOnP1 = 0.0;
	double _torqueOnP2 = 0.0;
	double _dampingConstant = -1.0;

	ContactInformation() {}

	ContactInformation(const int& otherID,
		const vector2d contactNormal,
		const vector2d contactPoint,
		const vector2d branchVectorFromP1,
		const vector2d branchVectorFromP2) :
		_otherID(otherID),
		_contactNormal(contactNormal),
		_contactPoint(contactPoint),
		_branchVectorFromP1(branchVectorFromP1),
		_branchVectorFromP2(branchVectorFromP2),
		_shearOverlap(vector2d(0.0, 0.0)),
		_normalForce(vector2d(0.0, 0.0)),
		_shearForce(vector2d(0.0, 0.0)),
		_force(vector2d(0.0, 0.0)),
		_forceOnP1(vector2d(0.0, 0.0)),
		_forceOnP2(vector2d(0.0, 0.0))
	{}

	ContactInformation(std::istream& is) {
		is >> _otherID
			>> _contactNormal >> _contactPoint
			>> _branchVector >> _branchVectorFromP1 >> _branchVectorFromP2
			>> _shearOverlap
			>> _normalForce >> _shearForce >> _force
			>> _forceOnP1 >> _forceOnP2
			>> _torqueOnP1 >> _torqueOnP2
			>> _dampingConstant;
	}

	void calculateTotalForceAndTorque() {
		_force = _normalForce + _shearForce;
		_forceOnP1 = _force;
		_forceOnP2 = -_force;
		_torqueOnP1 = _branchVectorFromP1.cross(_forceOnP1);
		_torqueOnP2 = _branchVectorFromP2.cross(_forceOnP2);
	}

	void print() {
		printf("Contact Information...\n");
		printf("%+25s : %d\n", "Contact with", _otherID);
		printf("%+25s : (%lf, %lf)\n", "Contact normal", _contactNormal.x(), _contactNormal.y());
		printf("%+25s : (%lf, %lf)\n", "Contact point", _contactPoint.x(), _contactPoint.y());
		printf("%+25s : (%lf, %lf)\n", "Shear overlap", _shearOverlap.x(), _shearOverlap.y());
		printf("%+25s : (%lf, %lf)\n", "Normal force", _normalForce.x(), _normalForce.y());
		printf("%+25s : (%lf, %lf)\n", "Shear force", _shearForce.x(), _shearForce.y());
		printf("%+25s : %lf\n", "Damping constant", _dampingConstant);
	}

	void print2File_contactInfo(FILE* file) {
		fprintf(file, FORMAT_ContactInformation_print2File_contactInfo,
			_otherID,
			_contactNormal.x(), _contactNormal.y(),
			_contactPoint.x(), _contactPoint.y(),
			_branchVector.x(), _branchVector.y(),
			_branchVectorFromP1.x(), _branchVectorFromP1.y(),
			_branchVectorFromP2.x(), _branchVectorFromP2.y(),
			_shearOverlap.x(), _shearOverlap.y(),
			_normalForce.x(), _normalForce.y(),
			_shearForce.x(), _shearForce.y(),
			_force.x(), _force.y(),
			_forceOnP1.x(), _forceOnP1.y(),
			_forceOnP2.x(), _forceOnP2.y(),
			_torqueOnP1,
			_torqueOnP2,
			_dampingConstant);
	}

	friend std::ostream& operator<<(std::ostream& os, const ContactInformation& c) {
		os 
			<< c._otherID << "\t"
			<< c._contactNormal << "\t"
			<< c._contactPoint << "\t"
			<< c._branchVector << "\t"
			<< c._branchVectorFromP1 << "\t"
			<< c._branchVectorFromP2 << "\t"
			<< c._shearOverlap << "\t"
			<< c._normalForce << "\t"
			<< c._shearForce << "\t"
			<< c._force << "\t"
			<< c._forceOnP1 << "\t"
			<< c._forceOnP2 << "\t"
			<< c._torqueOnP1 << "\t"
			<< c._torqueOnP2 << "\t"
			<< c._dampingConstant << "\n";
		return os;
	}

	//friend std::istream& operator>>(std::istream& is, ContactInformation& c) {
	//	is
	//		>> c._otherID
	//		>> c._contactNormal
	//		>> c._contactPoint
	//		>> c._branchVector
	//		>> c._branchVectorFromP1
	//		>> c._branchVectorFromP2
	//		>> c._shearOverlap
	//		>> c._normalForce
	//		>> c._shearForce
	//		>> c._force
	//		>> c._forceOnP1
	//		>> c._forceOnP2
	//		>> c._torqueOnP1
	//		>> c._torqueOnP2
	//		>> c._dampingConstant;
	//	return is;
	//}
};