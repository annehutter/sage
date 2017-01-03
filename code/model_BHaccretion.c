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
    if(BHaccretionMassEdd > BHaccrete * (1. - rad_efficiency))
    {
      dt_peak = getPeakTime_EddingtonLimited(BHmass, fEdd, rad_efficiency, BHaccrete * (1. - rad_efficiency));
    }else{
      dt_peak = dt;
    }
//     printf("dt = %e\t dt_peak = %e\t BHmass = %e\n", dt, dt_peak, BHmass);

    *BHaccretionRate = getBHaccretionRate_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt_peak);
    *BHaccretionMass = getBHaccretionMass_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt_peak);
//     printf("BHaccretonRate = %e\n", *BHaccretionRate);
//     printf("BHaccretonMass = %e\n", *BHaccretionMass);
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
    if(BHaccretionMassEdd > BHaccrete * (1. - rad_efficiency))
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

    if(BHaccretionMassEdd + BHmass <= Mpeak)
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

// Ryu / Inayoshi 2016
void accreteOnBH_Ryu(double BHmass, double GasMass, double Mvir, double cs, double cs_inflow, double gamma, double rad_efficiency, double parameter, double BHaccrete, double dt, double *BHaccretionRate, double *BHaccretionMass, double *Luminosity)
{
    if(BHmass <= 0.){
        BHmass = 1.e-8 * Hubble_h;
    }

    double fEdd = 1.;
    int type_acc = 0;
    double BondiRate = 0., InflowRate = 0., EddRate = 0.;
    double m_acc = 1., minRate = 0.;
    double dt_peak;

    double BondiMass = getBondiMass(BHmass, GasMass, Mvir, parameter);

    if(BondiMass > BHmass){
        dt_peak = getPeakTime_Inflow(cs_inflow, rad_efficiency, BHaccrete);
        if(dt_peak < dt) dt = dt_peak;
        type_acc = 4;

        *BHaccretionRate = getBHaccretionRate_Inflow(cs_inflow, rad_efficiency);
        *BHaccretionMass = getBHaccretionMass_Inflow(BHmass, cs_inflow, rad_efficiency, dt);
    }else{
        BondiRate = getBondiRate(BHmass, GasMass, Mvir, cs, gamma, parameter);
        InflowRate = getInflowRate(cs_inflow);
        EddRate = getEddRate(BHmass);

        if(InflowRate > BondiRate){
            minRate = BondiRate;
            type_acc = 1;
        }else{
            minRate = InflowRate;
            type_acc = 2;
        }
        
//         printf("BondiRate = %e \t InflowRate = %e \t EddRate = %e\n", BondiRate, InflowRate, EddRate);

        m_acc = minRate / EddRate * rad_efficiency;

        if(minRate < 3.e3 * EddRate){
            if(minRate > EddRate / rad_efficiency) 
            {
                minRate = EddRate / rad_efficiency;
                type_acc = 3;
            }
        }

        *BHaccretionRate = minRate;

        switch(type_acc)
        {
          case 1:
              dt_peak = getPeakTime_Bondi(BHmass, GasMass, Mvir, cs, gamma, parameter, rad_efficiency, 1., BHaccrete);
              if(dt_peak < dt) dt = dt_peak;
              *BHaccretionMass = getBHaccretionMass_Bondi(BHmass, GasMass, Mvir, cs, gamma, parameter, rad_efficiency, 1., dt);
              break;
          case 2:
              dt_peak = getPeakTime_Inflow(cs_inflow, rad_efficiency, BHaccrete);
              if(dt_peak < dt) dt = dt_peak;
              *BHaccretionMass = getBHaccretionMass_Inflow(BHmass, cs_inflow, rad_efficiency, dt);
              break;
          case 3:
              dt_peak = getPeakTime_EddingtonLimited(BHmass, fEdd, rad_efficiency, (1.-rad_efficiency)*BHaccrete);
              if(dt_peak < dt) dt = dt_peak;
//               printf("dt = %e\t dt_peak = %e\t BHmass = %e\n", dt, dt_peak, BHmass);
              *BHaccretionMass = getBHaccretionMass_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt);
              break;
          default:
            break;
        }

    }
    
    m_acc = getMaccr(*BHaccretionRate, BHmass);
    *Luminosity = getLuminosity_radIneff_log(*BHaccretionRate, m_acc, 20.);
    printf("%d: maccr = %e\t BHmass = %e\n", type_acc, m_acc, BHmass);
//     printf("BHaccretonRate = %e\n", *BHaccretionRate);
//     printf("BHaccretonMass = %e\n", *BHaccretionMass);
//     *Luminosity = getLuminosity_radIneff_quot(*BHaccretionRate, m_acc, 1.);
//     *Luminosity = getLuminosity_radEfficient(*BHaccretionRate, rad_efficiency);
}


void accreteOnBH_Park(double BHmass, double GasMass, double Mvir, double cs, double cs_inflow, double gamma, double mu, double rad_efficiency, double parameter, double BHaccrete, double dt, double *BHaccretionRate, double *BHaccretionMass, double *Luminosity)
{
    double tmpBHaccretionRate, tmpBHaccretionMass;
    double BHmaxRegime, BHmaxRegime2;
    double dt_tmp, dt_peak;
    double density;
    double maccr, accrRegime, densSubEddRegime, densCritRegime, densEddRegime, densSuperEddRegime;
    double fEdd = 1.;
    double lambdaRad = 1.;
    int type_acc = 0;
    
    double spectralIndex = QuasarSpectralIndex;
    
    if(BHmass <= 0.){
        BHmass = 1.e-8 * Hubble_h;
    }
    double massFact = 1.e10 / Hubble_h;

    // get mean photon energy to estimate temperature within Stroemgren sphere
    double meanPhotEnergy = getMeanPhotEnergy(spectralIndex);
    
    // compute n_H^crit and n_H^Edd
    densSubEddRegime = BHmass * massFact * getSubEddDensity(BHmass, cs, gamma, mu, rad_efficiency);
    densCritRegime = BHmass * massFact * getCriticalDensity(BHmass, cs, gamma, mu, meanPhotEnergy);
    densEddRegime = BHmass * massFact * getEddingtonDensity(BHmass, cs, gamma, mu, meanPhotEnergy, rad_efficiency);
    densSuperEddRegime = BHmass * massFact * getSuperEddDensity(BHmass, cs, gamma, mu, meanPhotEnergy, rad_efficiency);

    tmpBHaccretionMass = 0.;
    tmpBHaccretionRate = 0.;
    dt_tmp = 0.;
    
//     if(GasMass <= 0.) printf("\n\n\n GasMass = %e\n\n\n", GasMass);
    
    while(dt_tmp < dt && tmpBHaccretionMass < BHaccrete && GasMass > 0.)
    {
//         printf("dt_tmp = %e\t dt = %e\t tmpBHaccretionMass = %e\t BHaccrete = %e\n", dt_tmp, dt, tmpBHaccretionMass, BHaccrete);
        
        // determine accretion rate
        density = getDensityAtBondiRadius(BHmass, GasMass, Mvir, cs, gamma, mu, parameter);
//         printf("check density: BHmass = %e\t GasMass = %e\t Mvir = %e\t cs = %e\t gamma = %e\t mu = %e\t parameter = %e\n", BHmass, GasMass, Mvir, cs, gamma, mu, parameter);
        accrRegime = density * (BHmass + tmpBHaccretionMass) * massFact;
        
//         printf("regimes:\t subEdd: %e\t densCrit: %e\t densEdd: %e \t superEdd: %e\n", densSubEddRegime, densCritRegime, densEddRegime, densSuperEddRegime);
//         printf("%e\t%e\t%e\n", getSubEddDensity(BHmass, cs, gamma, mu, rad_efficiency), getCriticalDensity(BHmass, cs, gamma, mu, meanPhotEnergy), getEddingtonDensity(BHmass, cs, gamma, mu, meanPhotEnergy, rad_efficiency));
//         printf("density = %e\t accrRegime = %e\t BHmass = %e\t %e\n", density, accrRegime, (BHmass) * massFact, tmpBHaccretionMass);
        
        if(accrRegime > densSuperEddRegime)
        {
            BHmaxRegime = getBHmassBoundary(GasMass, Mvir, cs, mu, parameter, 1.e32);
            BHmaxRegime2 = getBHmassBoundary(GasMass, Mvir, cs, mu, parameter, densSuperEddRegime);
//             printf("BHmass = %e\t BHmaxRegime = %e\t BHmaxRegime2 = %e\n", BHmass, BHmaxRegime, BHmaxRegime2);
            if(BHmaxRegime2 > BHmaxRegime) BHmaxRegime = BHmaxRegime2;
//             printf("BHmax = %e\n", BHmaxRegime);
            
            dt_peak = getPeakTime_Bondi(BHmass, GasMass, Mvir, cs, gamma, parameter, rad_efficiency, 1., BHaccrete-tmpBHaccretionMass);
            
            if(dt_peak >= dt - dt_tmp) dt_peak = dt - dt_tmp;
            
            tmpBHaccretionRate = getBHaccretionRate_Bondi(BHmass, GasMass, Mvir, cs, gamma, parameter, rad_efficiency, 1.);
            tmpBHaccretionMass += getBHaccretionMass_Bondi(BHmass, GasMass, Mvir, cs, gamma, parameter, rad_efficiency, 1., dt_peak);
            
            type_acc = 1;
        }
        
        if(accrRegime > densEddRegime && accrRegime <= densSuperEddRegime)
        {
            BHmaxRegime = getBHmassBoundary(GasMass, Mvir, cs, mu, parameter, densSuperEddRegime);
            BHmaxRegime2 = getBHmassBoundary(GasMass, Mvir, cs, mu, parameter, densEddRegime);
//             printf("BHmass = %e\t BHmaxRegime = %e\t BHmaxRegime2 = %e\n", BHmass, BHmaxRegime, BHmaxRegime2);
            if(BHmaxRegime2 > BHmaxRegime) BHmaxRegime = BHmaxRegime2;
//             printf("BHmax = %e\n", BHmaxRegime);

            dt_peak = getPeakTime_EddingtonLimited(BHmass, fEdd, rad_efficiency, BHmaxRegime-(BHmass+tmpBHaccretionMass));
            
//             printf("dt = %e\t dt_peak = %e\n", dt, dt_peak);

            if(dt_peak >= dt - dt_tmp) dt_peak = dt - dt_tmp;

            tmpBHaccretionRate = getBHaccretionRate_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt_peak);
            tmpBHaccretionMass += getBHaccretionMass_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt_peak);
            
            type_acc = 2;
        }
        
        if(accrRegime > densSubEddRegime && accrRegime <= densEddRegime)
        {
            BHmaxRegime = getBHmassBoundary(GasMass, Mvir, cs, mu, parameter, densEddRegime);
            BHmaxRegime2 = getBHmassBoundary(GasMass, Mvir, cs, mu, parameter, densSubEddRegime);
//             printf("BHmass = %e\t BHmaxRegime = %e\t BHmaxRegime2 = %e\n", BHmass, BHmaxRegime, BHmaxRegime2);
            if(BHmaxRegime2 > BHmaxRegime) BHmaxRegime = BHmaxRegime2;
//             printf("BHmax = %e\n", BHmaxRegime);

            lambdaRad = getLambda_rad(BHmass, density, cs, gamma, mu, meanPhotEnergy);
            
            dt_peak = getPeakTime_Bondi(BHmass, GasMass, Mvir, cs, gamma, parameter, rad_efficiency, lambdaRad, BHmaxRegime-(BHmass+tmpBHaccretionMass));
            
            if(dt_peak >= dt - dt_tmp) dt_peak = dt - dt_tmp;
            
            tmpBHaccretionRate = getBHaccretionRate_Bondi(BHmass, GasMass, Mvir, cs, gamma, parameter, rad_efficiency, lambdaRad);
            tmpBHaccretionMass += getBHaccretionMass_Bondi(BHmass, GasMass, Mvir, cs, gamma, parameter, rad_efficiency, lambdaRad, dt_peak);
            
            type_acc = 3;
        }
        
        if(accrRegime <= densSubEddRegime)
        {
            BHmaxRegime = getBHmassBoundary(GasMass, Mvir, cs, mu, parameter, densSubEddRegime);
//             printf("BHmax = %e\n", BHmaxRegime);

            dt_peak = getPeakTime_Bondi(BHmass, GasMass, Mvir, cs, gamma, parameter, rad_efficiency, 1., BHmaxRegime-(BHmass+tmpBHaccretionMass));
            
            if(dt_peak >= dt - dt_tmp) dt_peak = dt - dt_tmp;
            
            tmpBHaccretionRate = getBHaccretionRate_Bondi(BHmass, GasMass, Mvir, cs, gamma, parameter, rad_efficiency, 1.);
            tmpBHaccretionMass += getBHaccretionMass_Bondi(BHmass, GasMass, Mvir, cs, gamma, parameter, rad_efficiency, 1., dt_peak);
            
            type_acc = 4;
        }

        dt_tmp += dt_peak;
        
        printf("type_acc = %d\n", type_acc);
    }
    
    *BHaccretionRate = tmpBHaccretionRate;
    *BHaccretionMass = tmpBHaccretionMass;
    
    maccr = getMaccr(*BHaccretionRate, BHmass);
    
    if(type_acc == 4)
    {
        *Luminosity = getLuminosity_radIneff_lowMaccr(BHmass, maccr);
    }
    else if(type_acc == 3 || type_acc == 2)
    {
        *Luminosity = getLuminosity_oscillations(BHmass, density, densCritRegime/(BHmass * massFact), cs, gamma, mu, meanPhotEnergy, rad_efficiency);
    }
    else if(type_acc == 1)
    {
        *Luminosity = getLuminosity_radIneff_quot(*BHaccretionRate, maccr, 1.);
    }
    else
    {
        *Luminosity = 0.;
    }
    
    printf("\n");
    
}

// void accreteOnBH_Pezulli()
// {
//     LuminosityEdd = ;
//     BHaccretionRateEdd = 16. * LuminosityEdd / (C * C);
//
//     tmp = BHaccretionRateEdd / ;
//     Luminosity = LuminosityEdd * A * (0.985 / (tmp + B) + 0.015 / (tmp + C));
// }
