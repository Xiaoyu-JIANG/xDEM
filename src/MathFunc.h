#pragma once

#define _USE_MATH_DEFINES
#include <math.h>

#define pow2(x) ((x) * (x))
#define pow3(x) ((x) * (x) * (x))
#define pow4(x) ((x) * (x) * (x) * (x))



namespace MathFunc
{
    // power of complex number
    inline void cpow(const double Re0, const double Im0, double& Re, double& Im, const double power) 
    {
        double magnitude = std::pow(sqrt(pow2(Re0) + pow2(Im0)), power);
        double theta = atan2(Im0, Re0) * power;
        Re = magnitude * cos(theta);
        Im = magnitude * sin(theta);
    }
}