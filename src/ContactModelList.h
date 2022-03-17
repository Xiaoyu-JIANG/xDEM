#pragma once

#include "BaseContactModel.h"

// Add new contact model here...

#include "ContactModel_LinearDissipativeFrictional.h"
#include "ContactModel_UniversalGravity.h"
class ContactModel_LinearDissipativeFrictioanl;
class ContactModel_UniversalGravity;



// Tell which contact model would you like to use in the simulation, defined as "ContactModel"

typedef ContactModel_LinearDissipativeFrictioanl ContactModel;
//typedef ContactModel_UniversalGravity ContactModel;