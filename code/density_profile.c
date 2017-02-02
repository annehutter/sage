#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"

double getConcentration_Bullock(double Mvir, double z)
{
    double result = 5. / (1.+z) * pow(Mvir / 1.3e14, -0.13);

    return result;
}

double getBvalue_Makino(double conc, double gamma)
{
    double result = 2.*conc / (9. * gamma * (Log(1. + conc) - conc / (1. + conc)));

    return result;
}

double getIntegral_Makino(double conc, double bvalue)
{
    int i;

    for(i=0; i<num; i++)
    {
        x = i * conc / num;
        dx = conc / num;
        integral += x * x * pow(1. + x, 13.5 * bvalue / x) * dx;
    }

    return integral;
}

double getRho0_Makino(double Mhot, double Rvir, double conc, double bvalue)
{
    double rs = Rvir / conc;
    double factor = Mhot / (4. * PI * rs * rs * rs) * exp(13.5 * bvalue);
    double integral = getIntegral_Makino(conc, bvalue);

    double result = factor / integral;

    return result;
}

double getRho_Makino(double r, double Mvir, double Rvir, double z)
{
    double conc = getConcentration_Bullock(Mvir, z);
    double bvalue = getBvalue_Makino(conc, 1.5);
    double rho0 = getRho0_Makino(conc, bvalue);

    double A = -0.178 * bvalue + 0.982;
    double rceff = 0.22 * Rvir / conc;
    double betaeff = 0.9 * bvalue;

    double tmp = r / rceff;
    double result = rho0 * A / pow(1. + tmp * tmp, 1.5 * betaeff);

    return result;
}
