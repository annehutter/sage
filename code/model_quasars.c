#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"

//--------------------------------------------------
//  SUB EDDINGTON DENSITY (where Bondi rate == Eddington rate)
//--------------------------------------------------

double getSubEddDensity(double BHmass, double cs, double gamma, double mu, double rad_efficiency)
{
    double massFact = 1.e10 / Hubble_h;

    double lambdaB = 1.12;
    if(gamma > 1.1) lambdaB = 0.25;
    lambdaB = lambdaB / 1.12;

    double temperature = 1.e-4 * (cs * cs * mu / (gamma * GAS_CONST));
    BHmass = BHmass * massFact;
    rad_efficiency = rad_efficiency * 10.;

    // result in cm^-3
    double result = 4.e6 / (BHmass * temperature * sqrt(temperature) * rad_efficiency);
    return result;
}

//--------------------------------------------------
// EDDINGTON DENSITY (where Bondi rate (with radiative feedback) == Eddington rate)
//--------------------------------------------------

double getEddingtonDensity(double BHmass, double cs, double gamma, double mu, double meanPhotEnergy, double rad_efficiency)
{
    double massFact = 1.e10 / Hubble_h;

    double lambdaB = 1.12;
    if(gamma > 1.1) lambdaB = 0.25;
    lambdaB = lambdaB / 1.12;

    double temperature = 1.e-4 * (cs * cs * mu / (gamma * GAS_CONST));
    BHmass = BHmass * massFact;
    rad_efficiency = rad_efficiency * 10.;

    // result in cm^-3
    double result = 4.e8 / (BHmass * temperature * rad_efficiency);
    return result;
}

//--------------------------------------------------
// CRITICAL DENSITY (where  Stroemgren radius approaches the effective accretion radius)
// (period of accretion burst correlates with Stroemgren radius; below gas depletion is dominated by accretion onto BH; above outflows dominate)
//--------------------------------------------------

double getCriticalDensity(double BHmass, double cs, double gamma, double mu, double meanPhotEnergy)
{
    double massFact = 1.e10 / Hubble_h;
    double lambdaB = 1.12;
    if(gamma > 1.1) lambdaB = 0.25;
    lambdaB = lambdaB / 1.12;

    meanPhotEnergy = meanPhotEnergy / 41.;
    double temperature = 1.e-4 * (cs * cs * mu / (gamma * GAS_CONST));
    BHmass = BHmass * massFact;

    // result in cm^-3
    double result = 5.e8 / (BHmass * pow(meanPhotEnergy, 0.5625));// * sqrt(lambdaB)) * temperature;
    return result;
}

//--------------------------------------------------
//  SUPER EDDINGTON DENSITY (where Bondi rate == Eddington rate)
//--------------------------------------------------

double getSuperEddDensity(double BHmass, double cs, double gamma, double mu, double meanPhotEnergy, double rad_efficiency)
{
    double massFact = 1.e10 / Hubble_h;

    double lambdaB = 1.12;
    if(gamma > 1.1) lambdaB = 0.25;
    lambdaB = lambdaB / 1.12;

    meanPhotEnergy = meanPhotEnergy / 41.;
    double temperature = 1.e-4 * (cs * cs * mu / (gamma * GAS_CONST));
    BHmass = BHmass * massFact;
    rad_efficiency = rad_efficiency * 10.;

    // result in cm^-3
    double result = 1.e10 / (BHmass * rad_efficiency * temperature * sqrt(temperature)) * pow(meanPhotEnergy, 0.42);    // only part not from Park but from Inayoshi 2016
    return result;
}

//--------------------------------------------------
// DENSITY AT BONDI RADIUS
//--------------------------------------------------

double getRhoAtBondiRadius(double BHmass, double GasMass, double Mvir, double cs, double gamma, double parameter)
{
    double massFact = 1.e10 / Hubble_h;
    double rvir = GRAVITY * Mvir * massFact / (cs * cs);
    double rBondi_div_rc = BHmass / (parameter * Mvir);

    double factor = 1. / (pow(parameter, 3) * (1./parameter - atan(1./parameter)));
    
    double lambdaB = 1.12;
    if(gamma > 1.1) lambdaB = 0.25;

    // result in g cm^-3
    double result = lambdaB * GasMass * massFact / (4. * PI * pow(rvir, 3) * SOLAR_MASS * SOLAR_MASS * (1. + rBondi_div_rc * rBondi_div_rc)) * factor;
    return result;
}

double getDensityAtBondiRadius(double BHmass, double GasMass, double Mvir, double cs, double gamma, double mu, double parameter)
{
    // result in g cm^-3
    double result = getRhoAtBondiRadius(BHmass, GasMass, Mvir, cs, gamma, parameter);
    // result in cm^-3
    result = result / (mu * PROTONMASS);
    return result;
}

//--------------------------------------------------
// BONDI RATE
//--------------------------------------------------

double getBondiRate(double BHmass, double GasMass, double Mvir, double cs, double gamma, double parameter)
{
    double rho = getRhoAtBondiRadius(BHmass, GasMass, Mvir, cs, gamma, parameter);
    double massFact = SOLAR_MASS * 1.e10 / Hubble_h;
    double Grav_times_BHmass = GRAVITY * BHmass * massFact;

    // result in g s^-1 * Msun^2
    double result = 4. * PI * Grav_times_BHmass * Grav_times_BHmass * rho / (cs * cs * cs);
    // result in Msun yr^-1
    return result / SOLAR_MASS * SEC_PER_YEAR;
}

//--------------------------------------------------
// EDDINGTON RATE
//--------------------------------------------------

double getEddRate(double BHmass)
{
    double massFact = 1.e10 / Hubble_h;

    // result in g s^-1 * Msun
    double result = 4. * PI * GRAVITY * BHmass * massFact * PROTONMASS / (THOMSON_CS * C);
    // result in Msun yr^-1
    return result * SEC_PER_YEAR;
}


//--------------------------------------------------
// INFLOW RATE
//--------------------------------------------------

double getInflowRate(double cs)
{
    // result in g s^-1
    double result = (cs * cs * cs) / GRAVITY;
    // result in Msun yr^-1
    return result / SOLAR_MASS * SEC_PER_YEAR;
}

//#####################################################################################################

//--------------------------------------------------
// MASS WITHIN THE BONDI RADIUS
//--------------------------------------------------

double getBondiMass(double BHmass, double GasMass, double Mvir, double parameter)
{
    double BHmass_div_Mvir = BHmass / (Mvir * parameter);
    double tmp = 1. / parameter;
    double factor = BHmass_div_Mvir - atan(BHmass_div_Mvir);
    double factor2 = tmp - atan(tmp);

    double result = GasMass * factor / factor2;
    // result in 1.e10/h Msun
    return result;
}

//#####################################################################################################

//--------------------------------------------------
// EDDINGTON LIMITED ACCRETION MODEL
//--------------------------------------------------

double getBHaccretionRate_EddingtonLimited(double BHmass, double fEdd, double rad_efficiency, double dt)
{
    double tmp = fEdd / (TEDD_YR * rad_efficiency) * (1. - rad_efficiency);
    double dt_yr = dt / (SEC_PER_YEAR * Hubble_h * CM_PER_KM) * CM_PER_MPC;
    if(BHmass == 0.) BHmass = 1.e-8 * Hubble_h;
    double result = BHmass * 1.e10 / Hubble_h * tmp * exp(tmp * dt_yr);
    // result in 1.e10/h Msun
    return result;
}

double getBHaccretionMass_EddingtonLimited(double BHmass, double fEdd, double rad_efficiency, double dt)
{
    double tmp = fEdd / (TEDD_YR * rad_efficiency) * (1. - rad_efficiency);
    double dt_yr = dt / (SEC_PER_YEAR * Hubble_h * CM_PER_KM) * CM_PER_MPC;
    if(BHmass == 0.) BHmass = 1.e-8 * Hubble_h;
    // printf("Eddington: BHmass = %e\t fEdd = %e\t rad_efficiency = %e\t t = %e\n", BHmass*1.e10/Hubble_h, fEdd, rad_efficiency, dt_yr);
    // printf("Eddington: BHmass = %e\n", BHmass * exp(tmp * dt_yr)*1.e10/Hubble_h);
    printf("tmp = %e\t dt_yr = %e\n", tmp, dt_yr);
    double result = BHmass * (exp(tmp * dt_yr) - 1.);   
    // result in 1.e10/h Msun
    return result;
}

double getPeakTime_EddingtonLimited(double BHmass, double fEdd, double rad_efficiency, double BHaccrete)
{
  double tmp = fEdd / (TEDD_YR * rad_efficiency) * (1. - rad_efficiency);
  double peakTime_yr = log((BHmass + BHaccrete) / BHmass) / tmp;
  double peakTime = peakTime_yr * (SEC_PER_YEAR * Hubble_h * CM_PER_KM) / CM_PER_MPC;
  return peakTime;
}


//#####################################################################################################

//--------------------------------------------------
// BONDI ACCRETION MODEL
//--------------------------------------------------

double getBHaccretionRate_Bondi(double BHmass, double GasMass, double Mvir, double cs, double gamma, double parameter, double rad_efficiency, double lambdaRad)
{
    // in Msun yr^-1
    double result = (1. - rad_efficiency) * lambdaRad * getBondiRate(BHmass, GasMass, Mvir, cs, gamma, parameter);
    return result;
}

double getBHaccretionMass_Bondi(double BHmass, double GasMass, double Mvir, double cs, double gamma, double parameter, double rad_efficiency, double lambdaRad, double dt)
{
    double massFact = 1.e10 / Hubble_h;
    double BondiRate = getBondiRate(BHmass, GasMass, Mvir, cs, gamma, parameter);
    double dt_yr = dt / (SEC_PER_YEAR * Hubble_h * CM_PER_KM) * CM_PER_MPC;

    // in Msun
    double result = (1. - rad_efficiency) * lambdaRad * BondiRate * dt_yr;
    // in 1.e10 h^-1 Msun
    result = result / massFact;
    return result;
}

double getPeakTime_Bondi(double BHmass, double GasMass, double Mvir, double cs, double gamma, double parameter, double rad_efficiency, double lambdaRad, double BHaccrete)
{
    double BondiRate = getBondiRate(BHmass, GasMass, Mvir, cs, gamma, parameter);
    double massFact = 1.e10 / Hubble_h;

    double peakTime_yr = BHaccrete * massFact / ((1. - rad_efficiency) * lambdaRad * BondiRate);
    double peakTime = peakTime_yr * (SEC_PER_YEAR * Hubble_h * CM_PER_KM) / CM_PER_MPC;
    printf("dt_peak = %e\n", peakTime);
    return peakTime;
}

//#####################################################################################################

//--------------------------------------------------
// INFLOW ACCRETION MODEL
//--------------------------------------------------

double getBHaccretionRate_Inflow(double cs, double rad_efficiency)
{
    // in Msun yr^-1
    double result = (1. - rad_efficiency) * getInflowRate(cs);
    return result;
}

double getBHaccretionMass_Inflow(double BHmass, double cs, double rad_efficiency, double dt)
{
    double massFact = 1.e10 / Hubble_h;
    double InflowRate = getInflowRate(cs);
    double dt_yr = dt / (SEC_PER_YEAR * Hubble_h * CM_PER_KM) * CM_PER_MPC;

    // in Msun
    double result = (1. - rad_efficiency) * InflowRate * dt_yr;
    // in 1.e10 h^-1 Msun
    result = result / massFact;
    return result;
}

double getPeakTime_Inflow(double cs, double rad_efficiency, double BHaccrete)
{
    double InflowRate = getInflowRate(cs);
    double massFact = 1.e10 / Hubble_h;

    double peakTime_yr = BHaccrete * massFact / ((1. - rad_efficiency) * InflowRate);
    double peakTime = peakTime_yr * (SEC_PER_YEAR * Hubble_h * CM_PER_KM) / CM_PER_MPC;
    printf("dt_peak = %e\n", peakTime);
    return peakTime;
}

//#####################################################################################################

//--------------------------------------------------
// PARK ACCRETION MODELS
//--------------------------------------------------

double getMeanPhotEnergy(double spectralIndex)
{
    double maxEnergy = 0.2e3;

    double result;
    if(spectralIndex > 1.){
          result = 13.6 * spectralIndex / (spectralIndex - 1.);
      }else if(spectralIndex == 1.){
          result = 13.6 * log(maxEnergy / 13.6);
      }else{
          result = 13.6 * spectralIndex / (1. - spectralIndex) * pow(maxEnergy / 13.6, spectralIndex);
      }
    return result;
}

double getLambda_rad(double BHmass, double density, double cs, double gamma, double mu, double meanPhotEnergy)
{
    double massFact = 1.e10 / Hubble_h;
    BHmass = BHmass * massFact;                                              // in Msun
    double temperature = 1.e-4 * (cs * cs * mu / (gamma * GAS_CONST));       // in K
    meanPhotEnergy = meanPhotEnergy / 41.;                            // in eV

    double result = 0.01 * pow(temperature, 2.5) / meanPhotEnergy;
    if(BHmass * density <= 1.e7){
        result = result * sqrt(density * 1.e-5);
    }
    return result;
}

double getFduty(double density, double density_crit, double cs, double gamma, double mu)
{
    double temperature = 1.e-4 * (cs * cs * mu / (gamma * GAS_CONST));       // in K

    double result;
    if(density > density_crit){
        result = 0.5;
    }else{
        result = 0.06 * sqrt(temperature);
    }
    return result;
}

double getTcycle(double BHmass, double density, double density_crit, double meanPhotEnergy, double rad_efficiency)
{
    double massFact = 1.e10 / Hubble_h;
    BHmass = BHmass * massFact;                                              // in Msun
    rad_efficiency = rad_efficiency * 10.;
    meanPhotEnergy = meanPhotEnergy / 41.;                                   // in eV

    double result;
    if(density > density_crit){
        result = 1.e9 * rad_efficiency / density * pow(meanPhotEnergy, -0.875);                               // in yr
    }else{
        result = 1.e5 * pow(BHmass * BHmass * rad_efficiency / density, 1./3.) * pow(meanPhotEnergy, -0.75);  // in yr
    }
    return result;
}

double getBHmassBoundary(double GasMass, double Mvir, double cs, double mu, double parameter, double BoundaryValue)
{
    double massFact = 1.e10 / Hubble_h;
    double tmp = (cs * cs) / GRAVITY;
    double factor = (1./parameter - atan(1./parameter));
    
    // in Msun
    double result = GasMass * tmp * tmp * tmp * factor / (4. * PI * Mvir * mu * PROTONMASS  * BoundaryValue * SOLAR_MASS * SOLAR_MASS);
//     printf("result = %e\t BoundrayValue = %e\t tmp = %e\t factor = %e\n", result, BoundaryValue, tmp, factor);
    // in 1.e10 h^-1 Msun
    result = result / massFact;
    return result;
}

double getLmaxDivLEdd(double BHmass, double density, double density_crit, double cs, double gamma, double mu, double meanPhotEnergy, double rad_efficiency)
{
    double massFact = 1.e10 / Hubble_h;
    BHmass = BHmass * massFact;                                              // in Msun
    double temperature = 1.e-4 * (cs * cs * mu / (gamma * GAS_CONST));       // in K
    rad_efficiency = rad_efficiency * 10.;
    meanPhotEnergy = meanPhotEnergy / 41.;                                   // in eV

    double result;
    if(density > density_crit){
        result = 5.e-6 * density * rad_efficiency * sqrt(temperature) / meanPhotEnergy;  
    }else{
        result = 6.e-7 * density * rad_efficiency * temperature / meanPhotEnergy; 
    }

    if(result > 1.) result = 1.;
    
    return result;
}

//#####################################################################################################

//--------------------------------------------------
// ACCRETION PARAMETER
//--------------------------------------------------

double getMaccr(double BHaccretionRate, double BHmass)
{
    double EddRate = getEddRate(BHmass);
    double result = BHaccretionRate / EddRate;
    return result;
}

//#####################################################################################################

//--------------------------------------------------
// LUMINOSITY MODELS
//--------------------------------------------------

double getEddLuminosity(double BHmass)
{
   double massFact = 1.e10 / Hubble_h;

    // result in erg s^-1 Msun^-1
    double result = 4. * PI * GRAVITY * BHmass * massFact * PROTONMASS * C / THOMSON_CS;
    return result * SEC_PER_YEAR;
}

double getLuminosity_radEfficient(double BHaccretionRate, double rad_efficiency)
{
    double result = rad_efficiency / (1. - rad_efficiency) * BHaccretionRate * C * C;// * SOLAR_MASS / SEC_PER_YEAR;
    return result;
}

double getLuminosity_radIneff_log(double BHaccretionRate, double maccr, double boundary)
{
    double rad_efficiency;
    if(maccr > boundary) rad_efficiency = 2./maccr * (1. + log(maccr / boundary));
    else rad_efficiency = 1./maccr;
    
    double result = rad_efficiency / (1. - rad_efficiency) * BHaccretionRate * C * C;
    return result;
}

double getLuminosity_radIneff_quot(double BHaccretionRate, double maccr, double boundary)
{
    double rad_efficiency = boundary / (10. + boundary * maccr);
    
    double result = rad_efficiency / (1. - rad_efficiency) * BHaccretionRate * C * C;
    return result;
}

double getLuminosity_radIneff_lowMaccr(double BHmass, double maccr)
{
    double EddRate = getEddLuminosity(BHmass * 1.e10 / Hubble_h);
    return 0.1 * EddRate * maccr * maccr;
}

double getLuminosity_oscillations(double BHmass, double density, double density_crit, double cs, double gamma, double mu, double meanPhotEnergy, double rad_efficiency)
{
    double LEdd = getEddLuminosity(BHmass);
    double Lmax = getLmaxDivLEdd(BHmass, density, density_crit, cs, gamma, mu, meanPhotEnergy, rad_efficiency) * LEdd;
    double fduty = getFduty(density, density_crit, cs, gamma, mu);
    double tcycle = getTcycle(BHmass, density, density_crit, meanPhotEnergy, rad_efficiency);
    
    srand((unsigned)time(NULL));
    double t = tcycle * ((double)rand()/(double)RAND_MAX);
    
    double result = Lmax * exp(t / (fduty * tcycle));
    printf("L = %e\n", result);
    
    return result;
}

double getLuminosity_Sadowski(double BHmass, double BHspin, double BHaccretionRate)
{
    double Avar = pow(0.9663 - 0.9292 * BHspin, -0.5639);
    double Bvar = pow(4.627 - 4.445 * BHspin, -0.5524);
    double Cvar = pow(827.3 - 718.1 * BHspin, -0.7060);

    double tmp = (16. * getEddRate(BHmass)) / BHaccretionRate;

    return Avar * ( 0.985 / (tmp + Bvar) + 0.015 / (tmp + Cvar)) * getEddLuminosity(BHmass);
}

//#####################################################################################################

//--------------------------------------------------
// MARULLI 2009 - HOPKINS MODEL
//--------------------------------------------------
double getBHaccretionRate_Hopkins(double BHmass, double rad_efficiency, double BHaccrete, double F, double dt)
{
  double dt_yr = dt / (SEC_PER_YEAR * Hubble_h * CM_PER_KM) * CM_PER_MPC;
  double BHpeak = BHmass + F * BHaccrete * (1.-rad_efficiency);
  double Lpeak = getEddLuminosity(BHpeak * 1.e10 / Hubble_h);
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
  double Lpeak = getEddLuminosity(BHpeak * 1.e10 / Hubble_h);
  double alpha = -0.95 + 0.32 * log10(Lpeak * SOLAR_MASS / (1.e12 * SOLAR_LUM));
  if(alpha > -0.2 || Lpeak * SOLAR_MASS > 1.e14 * SOLAR_LUM) alpha = -0.2;

  double tmp = pow(Lpeak * SOLAR_MASS /(1.e9 * SOLAR_LUM), alpha);
  double tmp2 = 1./alpha;

  double Avar = (1. - rad_efficiency) / (rad_efficiency * TEDD_YR) * BHpeak;
  double Bvar = tmp2 + 1.;
  double Cvar = tmp / 1.e9;

  double BHaccretionMass = BHpeak + Avar / (Bvar * Cvar) * (pow(1. + Cvar * dt_yr, Bvar) - 1.) - BHmass;
  return BHaccretionMass;
}

double getLuminosity_Hopkins(double BHmass, double rad_efficiency, double BHaccrete, double F, double dt)
{
  double dt_yr = dt / (SEC_PER_YEAR * Hubble_h * CM_PER_KM) * CM_PER_MPC;
  double BHpeak = BHmass + F * BHaccrete * (1.-rad_efficiency);
  double Lpeak = getEddLuminosity(BHpeak * 1.e10 / Hubble_h);
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
