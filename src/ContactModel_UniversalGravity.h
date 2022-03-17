#pragma once
#define _USE_MATH_DEFINES
#include <math.h>
#include "BaseContactModel.h"

class ContactModel_UniversalGravity : public BaseContactModel
{
public:
	ContactModel_UniversalGravity() : BaseContactModel() {}

	virtual const void calculateNormalForce(const ContactGeometry& contactGeometry, ContactInformation& contactInformation) const override;
	virtual const void calculateShearForce(const ContactGeometry& contactGeometry, ContactInformation& contactInformation) const override;

	const double calculateDampingConstact(const double& mass1, const double& mass2) const { return 0.0; }
	double checkTimeIncrement(const double& mass1, const double& mass2) { return 0.0; }
private:
	const double _gravitationalConstant = 6.674e-11;

};