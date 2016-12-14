#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"

// Eddington limited case
void accreteOnBH_EddingtonLimited(double BHmass, double BHaccrete, double rad_efficiency, double dt, double *BHaccretionRate, double *BHaccretionMass, double *Luminosity)
{
    double dt_peak;
    double fEdd = 1.;
    double BHaccretionMassEdd = getBHaccretionMass_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt);
    if(BHaccretionMassEdd - BHmass > BHaccrete * (1. - rad_efficiency))
    {
      dt_peak = getPeakTime_EddingtonLimited(BHmass, fEdd, rad_efficiency, BHaccrete * (1. - rad_efficiency));
    }else{
      dt_peak = dt;
    }
    *BHaccretionRate = getBHaccretionRate_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt_peak);
    *BHaccretionMass = getBHaccretionMass_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt_peak);
    *Luminosity = getLuminosity_radEfficient(*BHaccretionRate, rad_efficiency);
}

void accreteOnBH_EddingtonLimited_redshift(double BHmass, double BHaccrete, double redshift, double rad_efficiency, double dt, double *BHaccretionRate, double *BHaccretionMass, double *Luminosity)
{
    double fEdd;
    double dt_peak;

    if(redshift < 3.0)
    {
      fEdd = 0.3 * pow(0.25*(1 + redshift), 0.25);
    }else{
      fEdd = 0.3;
    }

    double BHaccretionMassEdd = getBHaccretionMass_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt);
    if(BHaccretionMassEdd - BHmass > BHaccrete * (1. - rad_efficiency))
    {
      dt_peak = getPeakTime_EddingtonLimited(BHmass, fEdd, rad_efficiency, BHaccrete * (1. - rad_efficiency));
    }else{
      dt_peak = dt;
    }
    *BHaccretionRate = getBHaccretionRate_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt_peak);
    *BHaccretionMass = getBHaccretionMass_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt_peak);
    *Luminosity = getLuminosity_radEfficient(*BHaccretionRate, rad_efficiency);
}


void accreteOnBH_Hopkins(double BHmass, double BHaccrete, double rad_efficiency, double dt, double *BHaccretionRate, double *BHaccretionMass, double *Luminosity)
{
    double Mpeak, dt_peak;

    double fEdd = 1.;
    double F = 0.7;
    double BHaccretionMassEdd = getBHaccretionMass_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt);

    if(BHmass <= 0.){
      Mpeak = BHaccretionMassEdd;
    }else{
      Mpeak = BHmass + F * BHaccrete * (1. - rad_efficiency);
    }

    if(BHaccretionMassEdd <= Mpeak)
    {
      *BHaccretionRate = getBHaccretionRate_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt);
      *BHaccretionMass = getBHaccretionMass_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt);
      *Luminosity = getLuminosity_radEfficient(*BHaccretionRate, rad_efficiency);
    }else{
      dt_peak = getPeakTime_EddingtonLimited(BHmass, fEdd, rad_efficiency, F * BHaccrete * (1. - rad_efficiency));
      *BHaccretionRate = getBHaccretionRate_Hopkins(BHmass, rad_efficiency, BHaccrete, F, dt-dt_peak);
      *BHaccretionMass = getBHaccretionMass_Hopkins(BHmass, rad_efficiency, BHaccrete, F, dt-dt_peak);
      *Luminosity = getLuminosity_Hopkins(BHmass, rad_efficiency, BHaccrete, F, dt-dt_peak);
      // printf("Mpeak = %e\t dt = %e\t dt_peak = %e\n", Mpeak, dt, dt_peak);
      // printf("BHmass = %e\t BHaccrete = %e BHaccretionMass = %e\n\n", BHmass, BHaccrete, BHaccretionMass-BHmass);
    }
}

// Hirschmann 2014 model
void accreteOnBH_(double mass, double BHmass, double BHaccrete, double temperature, double gamma, double mu, double rad_efficiency, double dt, double *BHaccretionRate, double *BHaccretionMass, double *Luminosity)
{
    double dt_peak;
    if(BHmass <= 0.){
        BHmass = 1.e-8 * Hubble_h;
    }
    double Tvir = 1.9e4 * pow(mass, 0.6666) * ((1.+20.) / 21.);
    printf("Tvir = %e\n", Tvir);
    double rho = mu * PROTONMASS * 9.5e5 * (temperature/1.e4)*(temperature/1.e4) / (BHmass*BHmass*1.e12) * Tvir*1.e-4;
    double fEdd = compute_fEdd(BHmass, rho, temperature, gamma, mu, 1., rad_efficiency, dt);
    double BHaccretionMass_tmp = getBHaccretionMass(BHmass, rho, temperature, gamma, mu, 1., rad_efficiency, dt);
    if(fEdd > 5000.){
        dt_peak = getPeakTime_Bondi(BHmass, rho, temperature, gamma, mu, BHaccrete);
    }else{
        dt_peak = getPeakTime_EddingtonLimited(BHmass, fEdd, rad_efficiency, BHaccrete);
    }
    if(dt_peak < dt)  dt = dt_peak;

    *BHaccretionRate = getBHaccretionRate(BHmass, rho, temperature, gamma, mu, 1., rad_efficiency, dt);
    *BHaccretionMass = getBHaccretionMass(BHmass, rho, temperature, gamma, mu, 1., rad_efficiency, dt);
    *Luminosity = getLuminosity(*BHaccretionRate, rad_efficiency, BHmass, rho, temperature, gamma, mu, fEdd);
}

// void accreteOnBH_Pezulli()
// {
//     LuminosityEdd = ;
//     BHaccretionRateEdd = 16. * LuminosityEdd / (C * C);
//
//     tmp = BHaccretionRateEdd / ;
//     Luminosity = LuminosityEdd * A * (0.985 / (tmp + B) + 0.015 / (tmp + C));
// }
