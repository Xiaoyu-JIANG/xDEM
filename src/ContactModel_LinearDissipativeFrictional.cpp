#include "ContactModel_LinearDissipativeFrictional.h"

const void ContactModel_LinearDissipativeFrictioanl::calculateNormalForce(const ContactGeometry& contactGeometry, ContactInformation& contactInformation) const
{
    contactInformation._normalForce =
        getElasticNormalForce(contactGeometry) + getDissipativeNormalForce(contactGeometry, contactInformation);
}

const void ContactModel_LinearDissipativeFrictioanl::calculateShearForce(const ContactGeometry& contactGeometry, ContactInformation& contactInformation) const
{
    //printf("Function called ... ContactModel_LinearDissipativeFrictioanl::calculateShearForc\n");
    double maximumShearForce = contactInformation._normalForce.getLength() * _frictionalCoeff;
    vector2d elasticShearForce = getElasticShearForce(contactInformation);
    vector2d dissipativeShearForce = getDissipativeShearForce(contactGeometry, contactInformation);
    vector2d shearForce = elasticShearForce + dissipativeShearForce;
    if (shearForce.getLength() > maximumShearForce) {
        shearForce = shearForce.getUnitVector() * maximumShearForce;
        if (elasticShearForce.getLength() > maximumShearForce) {
            contactInformation._shearOverlap = -shearForce / _shearStiffness;
        }
        else {
            contactInformation._shearOverlap = -elasticShearForce / _shearStiffness;
        }
    }
    contactInformation._shearForce = shearForce;
}

inline const vector2d ContactModel_LinearDissipativeFrictioanl::getElasticNormalForce(const ContactGeometry& contactGeometry) const
{
    return - contactGeometry.contactNormal * _normalStiffness * contactGeometry.normalOverlap;
}

inline const vector2d ContactModel_LinearDissipativeFrictioanl::getDissipativeNormalForce(const ContactGeometry& contactGeometry, ContactInformation& contactInformation) const
{
    return - contactGeometry.normalRelVel * contactInformation._dampingConstant;
}

inline const vector2d ContactModel_LinearDissipativeFrictioanl::getElasticShearForce(ContactInformation& contactInformation) const
{
    return - contactInformation._shearOverlap * _shearStiffness;
}

inline const vector2d ContactModel_LinearDissipativeFrictioanl::getDissipativeShearForce(const ContactGeometry& contactGeometry, ContactInformation& contactInformation) const
{
    return - contactGeometry.shearRelVel * contactInformation._dampingConstant;
}

double ContactModel_LinearDissipativeFrictioanl::checkTimeIncrement(const double& mass1, const double& mass2)
{
    double mass0 = (mass1 * mass2) / (mass1 + mass2);
    double dampingConstant = sqrt(mass0) * _factor;
    double dtCritical = M_PI / sqrt(_normalStiffness / mass0 - pow(dampingConstant / (2.0 * mass0), 2));
    return dtCritical;
}

void ContactModel_LinearDissipativeFrictioanl::setFriction(const double& fric)
{
    _frictionalCoeff = fric;
}

void ContactModel_LinearDissipativeFrictioanl::setRestitutionCoeff(const double& restitutionCoeff)
{
    _restitutionCoeff = restitutionCoeff;
}

void ContactModel_LinearDissipativeFrictioanl::recalculateFactor()
{
    double logOfRestitutionCoeff = log(_restitutionCoeff);
    _factor = 2.0 * sqrt(_normalStiffness) * abs(logOfRestitutionCoeff)
        / sqrt(M_PI * M_PI + logOfRestitutionCoeff * logOfRestitutionCoeff);
}

const double ContactModel_LinearDissipativeFrictioanl::calculateDampingConstact(const double& mass1, const double& mass2) const
{
    double mass0 = (mass1 * mass2) / (mass1 + mass2);
    return sqrt(mass0) * _factor;
}

const double ContactModel_LinearDissipativeFrictioanl::calculateDampingConstact(const double& mass) const
{
    return sqrt(mass) * _factor;
}
