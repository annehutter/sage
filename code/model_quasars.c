#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"

// double getBondiRate()
// {
//     return ;
// }

double getBHaccretionRate_EddingtonLimited(double BHmass, double fEdd, double rad_efficiency, double dt)
{
    double tmp = fEdd / (TEDD_YR * rad_efficiency) * (1. - rad_efficiency);
    double dt_yr = dt / (SEC_PER_YEAR * Hubble_h * CM_PER_KM) * CM_PER_MPC;
    return (BHmass + 1.e-8) * 1.e10 * tmp * exp(tmp * dt_yr);
}

double getBHaccretionMass_EddingtonLimited(double BHmass, double fEdd, double rad_efficiency, double dt)
{
    double tmp = fEdd / (TEDD_YR * rad_efficiency) * (1. - rad_efficiency);
    double dt_yr = dt / (SEC_PER_YEAR * Hubble_h * CM_PER_KM) * CM_PER_MPC;
    printf("tmp = %e\t dt_yr = %e\t BHmass = %e\n", tmp, dt_yr, BHmass);
    return (BHmass + 1.e-8) * exp(tmp * dt_yr);
}

double getLuminosity_radEfficient(double BHaccretionRate, double rad_efficiency)
{
    return rad_efficiency / (1. - rad_efficiency) * BHaccretionRate * C * C;// * SOLAR_MASS / SEC_PER_YEAR;
}

// double getLuminosity_radInefficient(double BHaccretionRate, double rad_efficiency)
// {
//     double tmp = BHacrretionRate_Bondi / BHaccretionRate_Edd;
//     return 0.1 * LEdd * tmp * tmp * 100.;
// }



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
