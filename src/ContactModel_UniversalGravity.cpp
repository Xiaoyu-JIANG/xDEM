#include "ContactModel_UniversalGravity.h"

const void ContactModel_UniversalGravity::calculateNormalForce(const ContactGeometry& contactGeometry, ContactInformation& contactInformation) const
{
    contactInformation._normalForce =
        contactGeometry.contactNormal *
        _gravitationalConstant *
        contactInformation._contactPoint.x() * // The class member _contactPoint is used to save mass information...
        contactInformation._contactPoint.y() / // The class member _contactPoint is used to save mass information...
        (contactGeometry.normalOverlap * contactGeometry.normalOverlap);
    //printf("[DEBUG] Function called...ContactModel_UniversalGravity::calculateNormalForce\n");
    //printf("        Calculated normal force = (%le, %le)\n", 
    //    contactInformation._normalForce.x(), contactInformation._normalForce.y());
    //printf("        Contact normal = (%le, %le)\n",
    //    contactGeometry.contactNormal.x(), contactGeometry.contactNormal.y());
    //printf("        Object mass = (%le, %le)\n",
    //    contactInformation._contactPoint.x(), contactInformation._contactPoint.y());
    //printf("        Distance = %le\n", contactGeometry.normalOverlap);
}

const void ContactModel_UniversalGravity::calculateShearForce(const ContactGeometry& contactGeometry, ContactInformation& contactInformation) const
{
    return void();
}
