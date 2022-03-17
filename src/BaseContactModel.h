#pragma once

#include "MathVector.h"
#include "ContactGeometry.h"
#include "ContactInformation.h"

class BaseContactModel
{
public:
	BaseContactModel() {}
	
	virtual const void calculateNormalForce(const ContactGeometry& contactGeometry, ContactInformation& contactInformation) const = 0;
	virtual const void calculateShearForce(const ContactGeometry& contactGeometry, ContactInformation& contactInformation) const = 0;
};

