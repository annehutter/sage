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
    if(BHaccrete > EPSILON * BHmass)
    {
        if(BHaccretionMassEdd > BHaccrete * (1. - rad_efficiency))
        {
          dt_peak = getPeakTime_EddingtonLimited(BHmass, fEdd, rad_efficiency, BHaccrete * (1. - rad_efficiency));
          *BHaccretionRate = 0.;//getBHaccretionRate_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt_peak);
          *BHaccretionMass = getBHaccretionMass_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt_peak);
          *Luminosity = 0.;//getLuminosity_radEfficient(*BHaccretionRate, rad_efficiency);
        }else{
          dt_peak = dt;
          *BHaccretionRate = getBHaccretionRate_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt_peak);
          *BHaccretionMass = getBHaccretionMass_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt_peak);
          *Luminosity = getLuminosity_radEfficient(*BHaccretionRate, rad_efficiency);
        }
    }
}

void accreteOnBH_EddingtonLimited_redshift(double BHmass, double BHaccrete, double redshift, double rad_efficiency, double dt, double *BHaccretionRate, double *BHaccretionMass, double *Luminosity)
{
    double fEdd;
    double dt_peak;

    if(redshift < 3.0)
    {
      fEdd = 0.8 * pow((1. + redshift)/4., 1.4);
    }else{
      fEdd = 0.8;
    }

    double BHaccretionMassEdd = getBHaccretionMass_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt);
    if(BHaccretionMassEdd > BHaccrete * (1. - rad_efficiency))
    {
      dt_peak = getPeakTime_EddingtonLimited(BHmass, fEdd, rad_efficiency, BHaccrete * (1. - rad_efficiency));
      *BHaccretionRate = 0.;//getBHaccretionRate_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt_peak);
      *BHaccretionMass = getBHaccretionMass_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt_peak);
      *Luminosity = 0.;//getLuminosity_radEfficient(*BHaccretionRate, rad_efficiency);
    }else{
      dt_peak = dt;
      *BHaccretionRate = getBHaccretionRate_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt_peak);
      *BHaccretionMass = getBHaccretionMass_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt_peak);
      *Luminosity = getLuminosity_radEfficient(*BHaccretionRate, rad_efficiency);
    }
}


void accreteOnBH_Hopkins(double BHmass, double BHaccrete, double rad_efficiency, double dt, double *BHaccretionRate, double *BHaccretionMass, double *Luminosity)
{
    double Mpeak, dt_peak;

    double fEdd = 1.;
    double F = 0.7;
    double BHaccretionMassEdd = getBHaccretionMass_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt);

    Mpeak = BHmass + F * BHaccrete * (1. - rad_efficiency);

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
    }
}

// Ryu / Inayoshi 2016
void accreteOnBH_Ryu(double BHmass, double GasMass, double Mvir, double Vvir, double cs, double cs_inflow, double gamma, double rad_efficiency, double parameter, double BHaccrete, double dt, double *BHaccretionRate, double *BHaccretionMass, double *Luminosity, int *type_acc_global)
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

    if(BHaccrete > EPSILON * BHmass)
    {
        if(BondiMass > BHmass){
            dt_peak = getPeakTime_Inflow(cs_inflow, rad_efficiency, BHaccrete);
            if(dt_peak < dt) dt = dt_peak;
            type_acc = 4;

            *BHaccretionRate = getBHaccretionRate_Inflow(cs_inflow, rad_efficiency);
            *BHaccretionMass = getBHaccretionMass_Inflow(BHmass, cs_inflow, rad_efficiency, dt);
        }else{
            BondiRate = getBondiRate(BHmass, GasMass, Mvir, Vvir, cs, gamma, parameter);
            InflowRate = getInflowRate(cs_inflow);
            EddRate = getEddRate(BHmass);

            if(InflowRate > BondiRate){
                minRate = BondiRate;
                type_acc = 1;
            }else{
                minRate = InflowRate;
                type_acc = 2;
            }

            m_acc = minRate / EddRate * rad_efficiency;

            if(minRate < 3.e3 * EddRate){
                if(minRate > EddRate / rad_efficiency)
                {
                    minRate = EddRate / rad_efficiency;
                    type_acc = 3;
                }
            }

            *BHaccretionRate = minRate;

            // printf(" type_acc = %d\t BHaccretionRate = %e\t BHmass = %e\n", type_acc, *BHaccretionRate, BHmass + *BHaccretionMass);
            switch(type_acc)
            {
              case 1:
                  dt_peak = getPeakTime_Bondi(BHmass, GasMass, Mvir, Vvir, cs, gamma, parameter, rad_efficiency, 1., BHaccrete);
                  if(dt_peak < dt) {
                      dt = dt_peak;
                      *BHaccretionRate = 0.;
                  }
                  *BHaccretionMass = getBHaccretionMass_Bondi(BHmass, GasMass, Mvir, Vvir, cs, gamma, parameter, rad_efficiency, 1., dt);
                  break;
              case 2:
                  dt_peak = getPeakTime_Inflow(cs_inflow, rad_efficiency, BHaccrete);
                  if(dt_peak < dt) {
                      dt = dt_peak;
                      *BHaccretionRate = 0.;
                  }
                  *BHaccretionMass = getBHaccretionMass_Inflow(BHmass, cs_inflow, rad_efficiency, dt);
                  break;
              case 3:
                  dt_peak = getPeakTime_EddingtonLimited(BHmass, fEdd, rad_efficiency, (1.-rad_efficiency)*BHaccrete);
                  if(dt_peak < dt) {
                      dt = dt_peak;
                      *BHaccretionRate = 0.;
                  }
                  *BHaccretionMass = getBHaccretionMass_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt);
                  break;
              default:
                break;
            }

        }

        *type_acc_global = type_acc;
        m_acc = getMaccr(*BHaccretionRate, BHmass);
        // if(m_acc < 1.){
        //     *Luminosity = getLuminosity_radIneff_lowMaccr(BHmass, m_acc);
        // }else{
            *Luminosity = getLuminosity_radIneff_log(*BHaccretionRate, m_acc, 20.);
        // }
        if(*Luminosity < 0.){
            printf("%d: maccr = %e\t BHmass = %e\t BHaccretionRate = %e\t Luminosity = %e\n", type_acc, m_acc, BHmass, *BHaccretionRate, *Luminosity);
        }
        // if(*Luminosity*2.e33/3.15e7 > 1.e45){
        //     printf("BondiRate = %e \t InflowRate = %e \t EddRate = %e\n", BondiRate, InflowRate, EddRate);
        //     printf("BHmass = %e\n", BHmass);
        //     printf("BHaccretonRate = %e\t maccr = %e\n", *BHaccretionRate, m_acc);
        //     printf("BHaccretonMass = %e\n", *BHaccretionMass);
        //     printf("Luminosity = %e\n\n", *Luminosity*2.e33/3.15e7);
        // }
    //     *Luminosity = getLuminosity_radIneff_quot(*BHaccretionRate, m_acc, 1.);
    //     *Luminosity = getLuminosity_radEfficient(*BHaccretionRate, rad_efficiency);
  }
  else
  {
      *BHaccretionRate = 0.;
      *BHaccretionMass = 0.;
      *Luminosity = 0.;
  }
}

// Ryu / Inayoshi 2016
void accreteOnBH_Ryu_MCF(double BHmass, double GasMass, double Mvir, double Vvir, double cs, double cs_inflow, double gamma, double mu, double logZ, double rad_efficiency, double parameter, double BHaccrete, double dt, double *BHaccretionRate, double *BHaccretionMass, double *Luminosity, int *type_acc_global)
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
        BondiRate = getBondiRate_MCF(BHmass, Vvir, cs, gamma, mu, logZ);
        InflowRate = getInflowRate(cs_inflow);
        EddRate = getEddRate(BHmass);

        if(InflowRate > BondiRate){
            minRate = BondiRate;
            type_acc = 1;
        }else{
            minRate = InflowRate;
            type_acc = 2;
        }

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
              dt_peak = getPeakTime_Bondi_MCF(BHmass, Vvir, cs, gamma, mu, logZ, rad_efficiency, 1., BHaccrete);
              if(dt_peak < dt) {
                  dt = dt_peak;
                  *BHaccretionRate = 0.;
              }else{
                  *BHaccretionRate = getBHaccretionRate_Bondi_MCF(BHmass, Vvir, cs, gamma, mu, logZ, rad_efficiency, 1., dt);
              }
              *BHaccretionMass = getBHaccretionMass_Bondi_MCF(BHmass, Vvir, cs, gamma, mu, logZ, rad_efficiency, 1., dt);
              break;
          case 2:
              dt_peak = getPeakTime_Inflow(cs_inflow, rad_efficiency, BHaccrete);
              if(dt_peak < dt) {
                  dt = dt_peak;
                  *BHaccretionRate = 0.;
              }
              *BHaccretionMass = getBHaccretionMass_Inflow(BHmass, cs_inflow, rad_efficiency, dt);
              break;
          case 3:
              dt_peak = getPeakTime_EddingtonLimited(BHmass, fEdd, rad_efficiency, (1.-rad_efficiency)*BHaccrete);
              if(dt_peak < dt) {
                  dt = dt_peak;
                  *BHaccretionRate = 0.;
              }else{
                  *BHaccretionRate = getBHaccretionMass_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt);
              }
              *BHaccretionMass = getBHaccretionMass_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt);
              break;
          default:
            break;
        }

    }

    *type_acc_global = type_acc;
    m_acc = getMaccr(*BHaccretionRate, BHmass);
    // if(m_acc < 1.){
    //     *Luminosity = getLuminosity_radIneff_lowMaccr(BHmass, m_acc);
    // }else{
        *Luminosity = getLuminosity_radIneff_log(*BHaccretionRate, m_acc, 20.);
    // }
    if(*Luminosity < 0.){
        printf("%d: maccr = %e\t BHmass = %e\t BHaccretionRate = %e\t Luminosity = %e\n", type_acc, m_acc, BHmass, *BHaccretionRate, *Luminosity);
    }
    // if(*Luminosity*2.e33/3.15e7 > 1.e45){
    //     printf("BondiRate = %e \t InflowRate = %e \t EddRate = %e\n", BondiRate, InflowRate, EddRate);
    //     printf("BHmass = %e\n", BHmass);
    //     printf("BHaccretonRate = %e\t maccr = %e\n", *BHaccretionRate, m_acc);
    //     printf("BHaccretonMass = %e\n", *BHaccretionMass);
    //     printf("Luminosity = %e\n\n", *Luminosity*2.e33/3.15e7);
    // }
//     *Luminosity = getLuminosity_radIneff_quot(*BHaccretionRate, m_acc, 1.);
//     *Luminosity = getLuminosity_radEfficient(*BHaccretionRate, rad_efficiency);
}

void accreteOnBH_Park(double BHmass, double GasMass, double Mvir, double Vvir, double cs, double cs_inflow, double gamma, double mu, double rad_efficiency, double parameter, double BHaccrete, double dt, double *BHaccretionRate, double *BHaccretionMass, double *Luminosity, int *type_acc_global)
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

    while(dt_tmp < dt && tmpBHaccretionMass < BHaccrete && GasMass > 1.e-10 && Mvir > 1.e-10)
    {
        // fprintf(stderr, "dt_tmp = %e\t dt = %e\t tmpBHaccretionMass = %e\t BHaccrete = %e\n", dt_tmp, dt, tmpBHaccretionMass, BHaccrete);

        // determine accretion rate
        density = getDensityAtBondiRadius(BHmass + tmpBHaccretionMass, GasMass, Mvir, Vvir, cs, gamma, mu, parameter);
        accrRegime = density * (BHmass + tmpBHaccretionMass) * massFact;

        // fprintf(stderr, "  check density: BHmass (in h^-1 10^10) = %e\t GasMass / Mvir = %e\t cs = %e cm/s\t Vvir = %e cm/s\n", BHmass + tmpBHaccretionMass, GasMass/Mvir, cs, Vvir);
        // fprintf(stderr, "  regimes:\t subEdd: %e\t densCrit: %e\t densEdd: %e \t superEdd: %e\n", densSubEddRegime, densCritRegime, densEddRegime, densSuperEddRegime);
        // fprintf(stderr, "  regime densities (cm^-3):\t subEdd: %e\t densEdd: %e\t superEdd: %e\n", getSubEddDensity(BHmass, cs, gamma, mu, rad_efficiency), getCriticalDensity(BHmass, cs, gamma, mu, meanPhotEnergy), getEddingtonDensity(BHmass, cs, gamma, mu, meanPhotEnergy, rad_efficiency));
        // fprintf(stderr, "  density = %e cm^-3\t BHmass = %e Msun\n", density, (BHmass + tmpBHaccretionMass) * massFact);
        // fprintf(stderr, "  accrRegime = %e Msun cm^-3\n", accrRegime);

        if(accrRegime > densSuperEddRegime)
        {
            BHmaxRegime = getBHmassBoundary(BHmass + tmpBHaccretionMass, GasMass, Mvir, Vvir, cs, gamma, mu, parameter, 1.e32);
            BHmaxRegime2 = getBHmassBoundary(BHmass + tmpBHaccretionMass, GasMass, Mvir, Vvir, cs, gamma, mu, parameter, densSuperEddRegime - EPSILON);

            if((BHmaxRegime2 > BHmaxRegime) || (BHmaxRegime != BHmaxRegime)) BHmaxRegime = BHmaxRegime2;

            dt_peak = getPeakTime_Bondi(BHmass + tmpBHaccretionMass, GasMass, Mvir, Vvir, cs, gamma, parameter, rad_efficiency, 1., BHaccrete-tmpBHaccretionMass);

            if(dt_peak >= dt - dt_tmp) dt_peak = dt - dt_tmp;

            tmpBHaccretionRate = getBHaccretionRate_Bondi(BHmass + tmpBHaccretionMass, GasMass, Mvir, Vvir, cs, gamma, parameter, rad_efficiency, 1.);
            tmpBHaccretionMass += getBHaccretionMass_Bondi(BHmass + tmpBHaccretionMass, GasMass, Mvir, Vvir, cs, gamma, parameter, rad_efficiency, 1., dt_peak);

            type_acc = 1;
        }

        if(accrRegime > densEddRegime && accrRegime <= densSuperEddRegime)
        {
            BHmaxRegime = getBHmassBoundary(BHmass + tmpBHaccretionMass, GasMass, Mvir, Vvir, cs, gamma, mu, parameter, densSuperEddRegime - EPSILON);
            BHmaxRegime2 = getBHmassBoundary(BHmass + tmpBHaccretionMass, GasMass, Mvir, Vvir, cs, gamma, mu, parameter, densEddRegime - EPSILON);

            if((BHmaxRegime2 > BHmaxRegime) || (BHmaxRegime != BHmaxRegime)) BHmaxRegime = BHmaxRegime2;

            dt_peak = getPeakTime_EddingtonLimited(BHmass + tmpBHaccretionMass, fEdd, rad_efficiency, BHmaxRegime-(BHmass+tmpBHaccretionMass));

            if(dt_peak >= dt - dt_tmp) dt_peak = dt - dt_tmp;

            tmpBHaccretionRate = getBHaccretionRate_EddingtonLimited(BHmass + tmpBHaccretionMass, fEdd, rad_efficiency, dt_peak);
            tmpBHaccretionMass += getBHaccretionMass_EddingtonLimited(BHmass + tmpBHaccretionMass, fEdd, rad_efficiency, dt_peak);

            type_acc = 2;
        }

        if(accrRegime > densSubEddRegime && accrRegime <= densEddRegime)
        {
            BHmaxRegime = getBHmassBoundary(BHmass + tmpBHaccretionMass, GasMass, Mvir, Vvir, cs, gamma, mu, parameter, densEddRegime - EPSILON);
            BHmaxRegime2 = getBHmassBoundary(BHmass + tmpBHaccretionMass, GasMass, Mvir, Vvir, cs, gamma, mu, parameter, densSubEddRegime - EPSILON);

            if((BHmaxRegime2 > BHmaxRegime) || (BHmaxRegime != BHmaxRegime)) BHmaxRegime = BHmaxRegime2;

            lambdaRad = getLambda_rad(BHmass + tmpBHaccretionMass, density, cs, gamma, mu, meanPhotEnergy);

            dt_peak = getPeakTime_Bondi(BHmass + tmpBHaccretionMass, GasMass, Mvir, Vvir, cs, gamma, parameter, rad_efficiency, lambdaRad, BHmaxRegime-(BHmass+tmpBHaccretionMass));

            if(dt_peak >= dt - dt_tmp) dt_peak = dt - dt_tmp;

            tmpBHaccretionRate = getBHaccretionRate_Bondi(BHmass + tmpBHaccretionMass, GasMass, Mvir, Vvir, cs, gamma, parameter, rad_efficiency, lambdaRad);
            tmpBHaccretionMass += getBHaccretionMass_Bondi(BHmass + tmpBHaccretionMass, GasMass, Mvir, Vvir, cs, gamma, parameter, rad_efficiency, lambdaRad, dt_peak);

            type_acc = 3;
        }

        if(accrRegime <= densSubEddRegime)
        {
            BHmaxRegime = getBHmassBoundary(BHmass + tmpBHaccretionMass, GasMass, Mvir, Vvir, cs, gamma, mu, parameter, densSubEddRegime - EPSILON);

            dt_peak = getPeakTime_Bondi(BHmass + tmpBHaccretionMass, GasMass, Mvir, Vvir, cs, gamma, parameter, rad_efficiency, 1., BHaccrete - tmpBHaccretionMass);

            if(dt_peak >= dt - dt_tmp) dt_peak = dt - dt_tmp;

            tmpBHaccretionRate = getBHaccretionRate_Bondi(BHmass + tmpBHaccretionMass, GasMass, Mvir, Vvir, cs, gamma, parameter, rad_efficiency, 1.);
            tmpBHaccretionMass += getBHaccretionMass_Bondi(BHmass + tmpBHaccretionMass, GasMass, Mvir, Vvir, cs, gamma, parameter, rad_efficiency, 1., dt_peak);

            type_acc = 4;
        }

        assert(dt_peak > 0. && "peak time should be larger than 0!!!");

        dt_tmp += dt_peak;

        // printf("  type_acc = %d\t BHaccretionRate = %e\t BHaccretionMass = %e\n\n", type_acc, tmpBHaccretionRate, tmpBHaccretionMass);
    }

    if(dt_tmp >= dt){
        *BHaccretionRate = tmpBHaccretionRate;
    }else{
        *BHaccretionRate = 0.;
    }
    *BHaccretionMass = tmpBHaccretionMass;

    maccr = getMaccr(*BHaccretionRate, BHmass);

    if(type_acc == 4)
    {
        if(*BHaccretionRate == 0.){
            *Luminosity = 0.;
        }else{
            // *Luminosity = getLuminosity_radIneff_quot(*BHaccretionRate, maccr, 10.);
            *Luminosity = getLuminosity_radEfficient(*BHaccretionRate, rad_efficiency);
        }
    }
    else if((type_acc == 3 || type_acc == 2) && (*BHaccretionRate > 0.))
    {
        *Luminosity = getLuminosity_oscillations(BHmass, density, densCritRegime/(BHmass * massFact), cs, gamma, mu, meanPhotEnergy, rad_efficiency);
    }
    else if(type_acc == 1)
    {
        *Luminosity = getLuminosity_radIneff_quot(*BHaccretionRate, maccr, 10.);
    }
    else
    {
        *Luminosity = 0.;
    }

    *type_acc_global = type_acc;
    // printf("\n");
}

void accreteOnBH_Park_MCF(double BHmass, double GasMass, double Mvir, double Vvir, double cs, double cs_inflow, double gamma, double mu, double logZ, double rad_efficiency, double parameter, double BHaccrete, double dt, double *BHaccretionRate, double *BHaccretionMass, double *Luminosity, int *type_acc_global)
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

    while(dt_tmp < dt && tmpBHaccretionMass < BHaccrete && GasMass > 1.e-10)
    {
        // fprintf(stderr, "dt_tmp = %e\t dt = %e\t tmpBHaccretionMass = %e\t BHaccrete = %e\n", dt_tmp, dt, tmpBHaccretionMass, BHaccrete);

        // determine accretion rate
        density = getDensityAtBondiRadius_MCF(BHmass, Vvir, cs, gamma, mu, logZ);
        accrRegime = density * (BHmass + tmpBHaccretionMass) * massFact;

        // fprintf(stderr, "  check density: BHmass (in h^-1 10^10) = %e\t GasMass / Mvir = %e\t cs = %e cm/s\t Vvir = %e cm/s\n", BHmass + tmpBHaccretionMass, GasMass/Mvir, cs, Vvir);
        // fprintf(stderr, "  regimes:\t subEdd: %e\t densCrit: %e\t densEdd: %e \t superEdd: %e\n", densSubEddRegime, densCritRegime, densEddRegime, densSuperEddRegime);
        // fprintf(stderr, "  regime densities (cm^-3):\t subEdd: %e\t densEdd: %e\t superEdd: %e\n", getSubEddDensity(BHmass, cs, gamma, mu, rad_efficiency), getCriticalDensity(BHmass, cs, gamma, mu, meanPhotEnergy), getEddingtonDensity(BHmass, cs, gamma, mu, meanPhotEnergy, rad_efficiency));
        // fprintf(stderr, "  density = %e cm^-3\t BHmass = %e Msun\n", density, (BHmass + tmpBHaccretionMass) * massFact);
        // fprintf(stderr, "  accrRegime = %e Msun cm^-3\n", accrRegime);

        if(accrRegime > densSuperEddRegime)
        {
            BHmaxRegime = getBHmassBoundary(BHmass, GasMass, Mvir, Vvir, cs, gamma, mu, parameter, 1.e32);
            BHmaxRegime2 = getBHmassBoundary(BHmass, GasMass, Mvir, Vvir, cs, gamma, mu, parameter, densSuperEddRegime - EPSILON);

            if(BHmaxRegime2 > BHmaxRegime) BHmaxRegime = BHmaxRegime2;

            dt_peak = getPeakTime_Bondi_MCF(BHmass, Vvir, cs, gamma, mu, logZ, rad_efficiency, 1., BHaccrete-tmpBHaccretionMass);

            if(dt_peak >= dt - dt_tmp) dt_peak = dt - dt_tmp;

            tmpBHaccretionRate = getBHaccretionRate_Bondi_MCF(BHmass, Vvir, cs, gamma, mu, logZ, rad_efficiency, 1., dt_peak);
            tmpBHaccretionMass += getBHaccretionMass_Bondi_MCF(BHmass, Vvir, cs, gamma, mu, logZ, rad_efficiency, 1., dt_peak);

            type_acc = 1;
        }

        if(accrRegime > densEddRegime && accrRegime <= densSuperEddRegime)
        {
            BHmaxRegime = getBHmassBoundary(BHmass, GasMass, Mvir, Vvir, cs, gamma, mu, parameter, densSuperEddRegime - EPSILON);
            BHmaxRegime2 = getBHmassBoundary(BHmass, GasMass, Mvir, Vvir, cs, gamma, mu, parameter, densEddRegime - EPSILON);

            if(BHmaxRegime2 > BHmaxRegime) BHmaxRegime = BHmaxRegime2;

            dt_peak = getPeakTime_EddingtonLimited(BHmass, fEdd, rad_efficiency, BHmaxRegime-(BHmass+tmpBHaccretionMass));

            if(dt_peak >= dt - dt_tmp) dt_peak = dt - dt_tmp;

            tmpBHaccretionRate = getBHaccretionRate_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt_peak);
            tmpBHaccretionMass += getBHaccretionMass_EddingtonLimited(BHmass, fEdd, rad_efficiency, dt_peak);

            type_acc = 2;
        }

        if(accrRegime > densSubEddRegime && accrRegime <= densEddRegime)
        {
            BHmaxRegime = getBHmassBoundary(BHmass, GasMass, Mvir, Vvir, cs, gamma, mu, parameter, densEddRegime - EPSILON);
            BHmaxRegime2 = getBHmassBoundary(BHmass, GasMass, Mvir, Vvir, cs, gamma, mu, parameter, densSubEddRegime - EPSILON);

            if(BHmaxRegime2 > BHmaxRegime) BHmaxRegime = BHmaxRegime2;

            lambdaRad = getLambda_rad(BHmass, density, cs, gamma, mu, meanPhotEnergy);

            dt_peak = getPeakTime_Bondi_MCF(BHmass, Vvir, cs, gamma, mu, logZ, rad_efficiency, lambdaRad, BHmaxRegime-(BHmass+tmpBHaccretionMass));

            if(dt_peak >= dt - dt_tmp) dt_peak = dt - dt_tmp;

            tmpBHaccretionRate = getBHaccretionRate_Bondi_MCF(BHmass, Vvir, cs, gamma, mu, logZ, rad_efficiency, lambdaRad, dt_peak);
            tmpBHaccretionMass += getBHaccretionMass_Bondi_MCF(BHmass, Vvir, cs, gamma, mu, logZ, rad_efficiency, lambdaRad, dt_peak);

            type_acc = 3;
        }

        if(accrRegime <= densSubEddRegime)
        {
            BHmaxRegime = getBHmassBoundary(BHmass, GasMass, Mvir, Vvir, cs, gamma, mu, parameter, densSubEddRegime - EPSILON);

            dt_peak = getPeakTime_Bondi_MCF(BHmass, Vvir, cs, gamma, mu, logZ, rad_efficiency, 1., BHmaxRegime-(BHmass+tmpBHaccretionMass));

            if(dt_peak >= dt - dt_tmp) dt_peak = dt - dt_tmp;

            tmpBHaccretionRate = getBHaccretionRate_Bondi_MCF(BHmass, Vvir, cs, gamma, mu, logZ, rad_efficiency, 1., dt_peak);
            tmpBHaccretionMass += getBHaccretionMass_Bondi_MCF(BHmass, Vvir, cs, gamma, mu, logZ, rad_efficiency, 1., dt_peak);

            type_acc = 4;
        }

        assert(dt_peak > 0. && "peak time should be larger than 0!!!");

        dt_tmp += dt_peak;

        // printf("  type_acc = %d\t BHaccretionRate = %e\t BHaccretionMass = %e\n\n", type_acc, tmpBHaccretionRate, tmpBHaccretionMass);
    }

    if(dt_tmp >= dt){
        *BHaccretionRate = tmpBHaccretionRate;
    }else{
        *BHaccretionRate = 0.;
    }
    *BHaccretionMass = tmpBHaccretionMass;

    maccr = getMaccr(*BHaccretionRate, BHmass);

    if(type_acc == 4)
    {
        if(*BHaccretionRate == 0.){
            *Luminosity = 0.;
        }else{
            // *Luminosity = getLuminosity_radIneff_quot(*BHaccretionRate, maccr, 10.);
            *Luminosity = getLuminosity_radEfficient(*BHaccretionRate, rad_efficiency);
        }
    }
    else if((type_acc == 3 || type_acc == 2) && (*BHaccretionRate > 0.))
    {
        *Luminosity = getLuminosity_oscillations(BHmass, density, densCritRegime/(BHmass * massFact), cs, gamma, mu, meanPhotEnergy, rad_efficiency);
    }
    else if(type_acc == 1)
    {
        *Luminosity = getLuminosity_radIneff_quot(*BHaccretionRate, maccr, 10.);
    }
    else
    {
        *Luminosity = 0.;
    }

    *type_acc_global = type_acc;
    // printf("\n");
}


// void accreteOnBH_Pezulli()
// {
//     LuminosityEdd = ;
//     BHaccretionRateEdd = 16. * LuminosityEdd / (C * C);
//
//     tmp = BHaccretionRateEdd / ;
//     Luminosity = LuminosityEdd * A * (0.985 / (tmp + B) + 0.015 / (tmp + C));
// }
