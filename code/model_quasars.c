#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"


//--------------------------------------------------
// EDDINGTON LIMITED & BONDI MODEL (MARULLI 2009 & HIRSCHMANN 2014)
//--------------------------------------------------

// ACCRETION MODEL

double getBHaccretionRate(double BHmass, double rho, double temperature, double gamma, double mu, double fEdd, double rad_efficiency, double dt)
{
    double BondiRate = getBHaccretionRate_Bondi(BHmass, rho, temperature, gamma, mu);
    double EddRate = getBHaccretionRate_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt);

    if(BondiRate < EddRate){
        return BondiRate;
    }else{
        return EddRate;
    }
}

double getBHaccretionMass(double BHmass, double rho, double temperature, double gamma, double mu, double fEdd, double rad_efficiency, double dt)
{
    double BondiRate = getBHaccretionRate_Bondi(BHmass, rho, temperature, gamma, mu);
    double EddRate = getBHaccretionRate_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt);

    if(BondiRate < EddRate){
        return getBHaccretionMass_Bondi(BHmass, rho, temperature, gamma, mu, dt);
    }else{
        return getBHaccretionMass_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt);
    }
}

double getBondiRate(double BHmass, double rho, double temperature, double gamma, double mu)
{
    double cs = sqrt(gamma * GAS_CONST * temperature / mu);
    return 4. * PI * GRAVITY * GRAVITY * BHmass * BHmass * rho / (cs * cs * cs);
}

double compute_fEdd(double BHmass, double rho, double temperature, double gamma, double mu, double fEdd, double rad_efficiency, double dt)
{
    double BondiRate = getBondiRate(BHmass * 1.e10 /Hubble_h, rho, temperature, gamma, mu);
    double EddRate = getBHaccretionRate_EddingtonLimited(BHmass * 1.e10 /Hubble_h, fEdd, rad_efficiency, dt);

    return BondiRate / EddRate;
}

// Eddington limited Accretion

double getBHaccretionRate_EddingtonLimited(double BHmass, double fEdd, double rad_efficiency, double dt)
{
    double tmp = fEdd / (TEDD_YR * rad_efficiency) * (1. - rad_efficiency);
    double dt_yr = dt / (SEC_PER_YEAR * Hubble_h * CM_PER_KM) * CM_PER_MPC;
    if(BHmass == 0.) BHmass = 1.e-8 * Hubble_h;
    return BHmass * 1.e10 / Hubble_h * tmp * exp(tmp * dt_yr);
}

double getBHaccretionMass_EddingtonLimited(double BHmass, double fEdd, double rad_efficiency, double dt)
{
    double tmp = fEdd / (TEDD_YR * rad_efficiency) * (1. - rad_efficiency);
    double dt_yr = dt / (SEC_PER_YEAR * Hubble_h * CM_PER_KM) * CM_PER_MPC;
    if(BHmass == 0.) BHmass = 1.e-8 * Hubble_h;
    // printf("Eddington: BHmass = %e\t fEdd = %e\t rad_efficiency = %e\t t = %e\n", BHmass*1.e10/Hubble_h, fEdd, rad_efficiency, dt_yr);
    // printf("Eddington: BHmass = %e\n", BHmass * exp(tmp * dt_yr)*1.e10/Hubble_h);
    printf("tmp = %e\t dt = %e\t fEdd = %e\t BHmass = %e\n", tmp, dt_yr, fEdd, BHmass);
    return BHmass * exp(tmp * dt_yr);
}

// Bondi Accretion

double getBHaccretionRate_Bondi(double BHmass, double rho, double temperature, double gamma, double mu)
{
    return getBondiRate(BHmass * 1.e10 / Hubble_h, rho, temperature, gamma, mu);
}

double getBHaccretionMass_Bondi(double BHmass, double rho, double temperature, double gamma, double mu, double dt)
{
    double cs = sqrt(gamma * GAS_CONST * temperature / mu);
    double tmp = 4. * PI * GRAVITY * GRAVITY * rho / (cs * cs * cs);

    if(dt * tmp * BHmass * 1.e10 / Hubble_h < 1.)
    {
        return BHmass / ( 1. - tmp * dt * BHmass * 1.e10 / Hubble_h );
    }else{
        return 1.e20;
    }
}

// LUMINOSITY MODEL

double getLuminosity(double BHaccretionRate, double rad_efficiency, double BHmass, double rho, double temperature, double gamma, double mu, double fEdd)
{
    if(fEdd > 0.1)
    {
        return getLuminosity_radEfficient(BHaccretionRate, rad_efficiency);
    }else{
        return getLuminosity_radInefficient(BHmass, fEdd);
    }
}

double getEddingtonLuminosity(double BHmass)
{
    return 4. * PI * GRAVITY * BHmass * PROTONMASS * C / THOMSON_CS;
}

double getLuminosity_radEfficient(double BHaccretionRate, double rad_efficiency)
{
    return rad_efficiency / (1. - rad_efficiency) * BHaccretionRate * C * C;// * SOLAR_MASS / SEC_PER_YEAR;
}

double getLuminosity_radInefficient(double BHmass, double fEdd)
{
    double EddRate = getEddingtonLuminosity(BHmass * 1.e10 / Hubble_h);
    return 0.1 * EddRate * fEdd * fEdd * 100.;
}


//--------------------------------------------------
// MARULLI 2009 - HOPKINS MODEL
//--------------------------------------------------

double getPeakTime(double BHmass, double fEdd, double rad_efficiency, double BHaccretionMassEdd)
{
  double tmp = fEdd / (TEDD_YR * rad_efficiency) * (1. - rad_efficiency);
  double peakTime_yr = log((BHmass + BHaccretionMassEdd) / BHmass) / tmp;
  double peakTime = peakTime_yr * (SEC_PER_YEAR * Hubble_h * CM_PER_KM) / CM_PER_MPC;
  return peakTime;
}

double getBHaccretionRate_Hopkins(double BHmass, double rad_efficiency, double BHaccrete, double F, double dt)
{
  double dt_yr = dt / (SEC_PER_YEAR * Hubble_h * CM_PER_KM) * CM_PER_MPC;
  double BHpeak = BHmass + F * BHaccrete * (1.-rad_efficiency);
  double Lpeak = getEddingtonLuminosity(BHpeak * 1.e10 / Hubble_h);
  double alpha = -0.95 + 0.32 * log10(Lpeak * SOLAR_MASS/ (1.e12 * SOLAR_LUM));
  if(alpha > -0.2) alpha = -0.2;

  double tmp = pow(Lpeak * SOLAR_MASS/(1.e9 * SOLAR_LUM), alpha);
  double tmp2 = 1./alpha;

  double Avar = (1. - rad_efficiency) / (rad_efficiency * TEDD_YR) * BHpeak * 1.e10 / Hubble_h;
  double Bvar = tmp2 + 1.;
  double Cvar = tmp / 1.e9;

  double BHaccretionRate = Avar * pow(1. + Cvar * dt_yr, Bvar - 1.);
  return BHaccretionRate;
}

double getBHaccretionMass_Hopkins(double BHmass, double rad_efficiency, double BHaccrete, double F, double dt)
{
  double dt_yr = dt / (SEC_PER_YEAR * Hubble_h * CM_PER_KM) * CM_PER_MPC;
  double BHpeak = BHmass + F * BHaccrete * (1.-rad_efficiency);
  double Lpeak = getEddingtonLuminosity(BHpeak * 1.e10 / Hubble_h);
  double alpha = -0.95 + 0.32 * log10(Lpeak * SOLAR_MASS / (1.e12 * SOLAR_LUM));
  if(alpha > -0.2 || Lpeak * SOLAR_MASS > 1.e14 * SOLAR_LUM) alpha = -0.2;

  double tmp = pow(Lpeak * SOLAR_MASS /(1.e9 * SOLAR_LUM), alpha);
  double tmp2 = 1./alpha;

  double Avar = (1. - rad_efficiency) / (rad_efficiency * TEDD_YR) * BHpeak;
  double Bvar = tmp2 + 1.;
  double Cvar = tmp / 1.e9;

  double BHaccretionMass = BHpeak + Avar / (Bvar * Cvar) * (pow(1. + Cvar * dt_yr, Bvar) - 1.);
  return BHaccretionMass;
}

double getLuminosity_Hopkins(double BHmass, double rad_efficiency, double BHaccrete, double F, double dt)
{
  double dt_yr = dt / (SEC_PER_YEAR * Hubble_h * CM_PER_KM) * CM_PER_MPC;
  double BHpeak = BHmass + F * BHaccrete * (1.-rad_efficiency);
  double Lpeak = getEddingtonLuminosity(BHpeak * 1.e10 / Hubble_h);
  double alpha = -0.95 + 0.32 * log10(Lpeak * SOLAR_MASS / (1.e12 * SOLAR_LUM));
  if(alpha > -0.2) alpha = -0.2;

  double tmp = pow(Lpeak * SOLAR_MASS /(1.e9 * SOLAR_LUM), alpha);
  double tmp2 = 1./alpha;

  double Bvar = tmp2 + 1.;
  double Cvar = tmp / 1.e9;

  double luminosity = Lpeak * pow(1. + Cvar * dt_yr, Bvar - 1.) * SEC_PER_YEAR;
  return luminosity;
}

// //--------------------------------------------------
// // PEZULLI 2016 MODEL
// //--------------------------------------------------
//
// double Sadowski_fit(double MEdd_Maccrete, double BHspin)
// {
//     double A;
//     double B;
//     double C;
//
//     A = pow(0.9663 - 0.9292*BHspin, -0.5639);
//     B = pow(4.627 - 4.445*BHspin, -0.5524);
//     C = pow(827.3 - 718.1*BHspin, -0.7060);
//
//     return A * ( 0.985 / (MEdd_Maccrete + B) + 0.015 / (MEdd_Maccrete * C) );
// }
//
// double Pezulli_BHaccretion()
// {
//     // Pezulli 2006 Model
//     Reff = 2.884031503126606e-06*pow(Gal[merger_centralgal].BulgeMass, 0.56);
//     rb = Reff/1.8153;
//     BulgeGasMass = Gal[merger_centralgal].HotGas;
//
//     rd = Gal[merger_centralgal].DiskScaleRadius;
//     x = Reff/rd;
//     GasSurfaceDensity = (Gal[merger_centralgal].ColdGas + (Gal[merger_centralgal].StellarMass - Gal[merger_centralgal].BulgeMass))/(2.*PI * rd * rd);
//     DiskGasMass = Gal[merger_centralgal].ColdGas + (Gal[merger_centralgal].StellarMass - Gal[merger_centralgal].BulgeMass);
//
//     cvir =
//     rs = Rvir/cvir;
//     A = log(1. + cvir) - cvir / (1. + cvir);
//
//     DiskGasMass_t = DiskGasMass * (1.-exp(-Reff/rd)*(Reff+rd)/rd);
//     BulgeGasMass_t = BulgeGasMass * Reff*Reff / ((Reff + rb) * (Reff + rb));
//     DMmass_t = Gal[merger_centralgal].Mvir * (log((Reff + rs) / rs) - Reff/(Reff + rs));
//
//     HaloMass = DiskGasMass + BulgeGasMass + Gal[merger_centralgal].Mvir;
//     fgal = (DiskGassMass + BulgeGasMass) / HaloMass;
//
//     vd_sq = PI * GRAVITY * GasSurfaceDensity * x * ();
//     vb_sq = GRAVITY * (Gal[merger_centralgal].BulgeMass + Gal[merger_centralgal].HotGas) * Reff /((Reff+rb)*(Reff+rb));
//     vDM_sq = GRAVITY * (1.-fgal) * Gal[merger_centralgal].Mvir / Reff;
//     vc_Reff = sqrt(vd_sq*vd_sq + vb_sq*vb_sq + vDM_sq*vDM_sq);
//     taub = Reff/vc_Reff;
//
//     faccrete = beta;
//     Maccrete = faccrete * Gal[merger_centralgal].ColdGas/taub;
//
//     LEdd = ;
//     Lbol = LEdd * Sadowski_fit(16.*LEdd/(C*C)/Maccrete, 0.98);
//
//     BHaccrete = Maccrete - Lbol / (C * C));
//
//     return BHaccrete;
// }
//
//
// //--------------------------------------------------
// // RYU 2016 MODEL
// //--------------------------------------------------
//
// double Ryu_BHaccretion()
// {
//     Maccrete_Bondi = ;
//     Maccrete_Edd = LEdd / (C*C);
//     Maccrete_in = 0.5;  // in Msun/yr
//     rad_efficiency = 0.1;
//
//     if(Maccrete_Bondi < Maccrete_in) min = Maccrete_Bondi;
//     else min = Maccrete_in;
//
//     if(min<3000.)
//     {
//         if(Maccrete_Edd/rad_efficiency < min) min = Maccrete_Edd/rad_efficiency;
//     }
//
//     BHaccrete = min;
//
//     return BHaccrete;
// }
//
// //--------------------------------------------------
// // PARK 2016 MODEL
// //--------------------------------------------------
//
// //--------------------------------------------------
// // SIMPLE MODEL
// //--------------------------------------------------
//
// void grow_black_hole(int merger_centralgal, double mass_ratio)
// {
//   double BHaccrete, metallicity;
//
//   if(Gal[merger_centralgal].ColdGas > 0.0)
//   {
//
//     // Simple Model
//     BHaccrete = BlackHoleGrowthRate * mass_ratio /
//       (1.0 + pow(280.0 / Gal[merger_centralgal].Vvir, 2.0)) * Gal[merger_centralgal].ColdGas;
//
//
//     // cannot accrete more gas than is available!
//     if(BHaccrete > Gal[merger_centralgal].ColdGas)
//       BHaccrete = Gal[merger_centralgal].ColdGas;
//
//     metallicity = get_metallicity(Gal[merger_centralgal].ColdGas, Gal[merger_centralgal].MetalsColdGas);
//     Gal[merger_centralgal].BlackHoleMass += BHaccrete;
//     Gal[merger_centralgal].ColdGas -= BHaccrete;
//     Gal[merger_centralgal].MetalsColdGas -= metallicity * BHaccrete;
//
//     quasar_mode_wind(merger_centralgal, BHaccrete);
//   }
// }
