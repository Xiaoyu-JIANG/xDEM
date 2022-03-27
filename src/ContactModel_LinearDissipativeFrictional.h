#pragma once
#define _USE_MATH_DEFINES
#include <math.h>
#include "BaseContactModel.h"

class ContactModel_LinearDissipativeFrictioanl : public BaseContactModel
{
public:
	ContactModel_LinearDissipativeFrictioanl() : BaseContactModel() {}

	ContactModel_LinearDissipativeFrictioanl(const double& normalStiffness,
		const double& shearStiffness,
		const double& restitutionCoeff,
		const double& frictionalCoeff) : 
		BaseContactModel(),
		_normalStiffness(normalStiffness),
		_shearStiffness(shearStiffness),
		_restitutionCoeff(restitutionCoeff),
		_frictionalCoeff(frictionalCoeff)
	{
		double logOfRessitutionCoeff = log(_restitutionCoeff);
		_factor = 2.0 * sqrt(_normalStiffness) * abs(logOfRessitutionCoeff) 
			/ sqrt(M_PI * M_PI + logOfRessitutionCoeff * logOfRessitutionCoeff);
	}
	
	virtual const void calculateNormalForce(const ContactGeometry& contactGeometry, ContactInformation& contactInformation) const override;
	virtual const void calculateShearForce(const ContactGeometry& contactGeometry, ContactInformation& contactInformation) const override;
	const double calculateDampingConstact(const double& mass1, const double& mass2) const;
	const double calculateDampingConstact(const double& mass) const;

	double checkTimeIncrement(const double& mass1, const double& mass2);

	void setFriction(const double& fric);
	void setRestitutionCoeff(const double& restitutionCoeff);
	void recalculateFactor();
	const double getNormalStiffness() const { return _normalStiffness; }

private:
	inline const vector2d getElasticNormalForce(const ContactGeometry& contactGeometry) const;
	inline const vector2d getDissipativeNormalForce(const ContactGeometry& contactGeometry, ContactInformation& contactInformation) const;
	inline const vector2d getElasticShearForce(ContactInformation& contactInformation) const;
	inline const vector2d getDissipativeShearForce(const ContactGeometry& contactGeometry, ContactInformation& contactInformation) const;

	double _normalStiffness = 0.0;
	double _shearStiffness = 0.0;
	double _restitutionCoeff = 0.0;
	double _frictionalCoeff = 0.0;
	double _factor = 0.0;
};