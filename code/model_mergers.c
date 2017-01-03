#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"



double estimate_merging_time(int sat_halo, int mother_halo, int ngal)
{
  double coulomb, mergtime, SatelliteMass, SatelliteRadius;

  if(sat_halo == mother_halo)
  {
    printf("\t\tSnapNum, Type, IDs, sat radius:\t%i\t%i\t%i\t%i\t--- sat/cent have the same ID\n",
      Gal[ngal].SnapNum, Gal[ngal].Type, sat_halo, mother_halo);
    return -1.0;
  }

  coulomb = log(Halo[mother_halo].Len / ((double) Halo[sat_halo].Len) + 1);

  SatelliteMass = get_virial_mass(sat_halo) + Gal[ngal].StellarMass + Gal[ngal].ColdGas;
  SatelliteRadius = get_virial_radius(mother_halo);

  if(SatelliteMass > 0.0 && coulomb > 0.0)
    mergtime = 2.0 *
    1.17 * SatelliteRadius * SatelliteRadius * get_virial_velocity(mother_halo) / (coulomb * G * SatelliteMass);
  else
    mergtime = -1.0;

  return mergtime;

}

void deal_with_galaxy_merger(int p, int merger_centralgal, int centralgal, double time, double dt, int halonr, int step)
{
  double mi, ma, mass_ratio;

  // calculate mass ratio of merging galaxies
  if(Gal[p].StellarMass + Gal[p].ColdGas <
    Gal[merger_centralgal].StellarMass + Gal[merger_centralgal].ColdGas)
  {
    mi = Gal[p].StellarMass + Gal[p].ColdGas;
    ma = Gal[merger_centralgal].StellarMass + Gal[merger_centralgal].ColdGas;
  }
  else
  {
    mi = Gal[merger_centralgal].StellarMass + Gal[merger_centralgal].ColdGas;
    ma = Gal[p].StellarMass + Gal[p].ColdGas;
  }

  if(ma > 0)
    mass_ratio = mi / ma;
  else
    mass_ratio = 1.0;

  add_galaxies_together(merger_centralgal, p);

  // grow black hole through accretion from cold disk during mergers, a la Kauffmann & Haehnelt (2000)
  if(AGNrecipeOn)
  {
    if(TrackBHgrowthOn != 0)
    {
#ifdef DEBUG
      printf("deal with galaxy merger\t%d: mergTime = %e\t mergTimeInit = %e\t dt = %e\n", p, Gal[p].MergTime, Gal[p].MergTimeInit, dt*(step+1));
#endif
      grow_black_hole_trackBHgrowth(merger_centralgal, p, mass_ratio, time);
    }else if(ContinuousAccretionOn != 0)
    {
      grow_black_hole_continuousAccretion(merger_centralgal, p, mass_ratio, time, dt);
    }else{
      grow_black_hole(merger_centralgal, mass_ratio);
    }
  }

  // starburst recipe similar to Somerville et al. 2001
  collisional_starburst_recipe(mass_ratio, merger_centralgal, centralgal, time, dt, halonr, 0, step);

  if(mass_ratio > 0.1)
		Gal[merger_centralgal].TimeOfLastMinorMerger = time;

  if(mass_ratio > ThreshMajorMerger)
  {
    make_bulge_from_burst(merger_centralgal);
    Gal[merger_centralgal].TimeOfLastMajorMerger = time;
    Gal[p].mergeType = 2;  // mark as major merger
  }
  else
  {
    Gal[p].mergeType = 1;  // mark as minor merger
  }

}


void grow_black_hole(int merger_centralgal, double mass_ratio)
{
  double BHaccrete, metallicity;

  if(Gal[merger_centralgal].ColdGas > 0.0)
  {
    BHaccrete = BlackHoleGrowthRate * mass_ratio /
      (1.0 + pow(280.0 / Gal[merger_centralgal].Vvir, 2.0)) * Gal[merger_centralgal].ColdGas;

    // cannot accrete more gas than is available!
    if(BHaccrete > Gal[merger_centralgal].ColdGas)
      BHaccrete = Gal[merger_centralgal].ColdGas;

    metallicity = get_metallicity(Gal[merger_centralgal].ColdGas, Gal[merger_centralgal].MetalsColdGas);
    Gal[merger_centralgal].BlackHoleMass += BHaccrete;
    Gal[merger_centralgal].ColdGas -= BHaccrete;
    Gal[merger_centralgal].MetalsColdGas -= metallicity * BHaccrete;

    Gal[merger_centralgal].QuasarModeBHaccretionMass += BHaccrete;

    quasar_mode_wind(merger_centralgal, BHaccrete);
  }
}

void grow_black_hole_trackBHgrowth(int merger_centralgal, int p, double mass_ratio, double time)
{
  double BHaccrete, metallicity;

  if(Gal[merger_centralgal].ColdGas > 0.0)
  {
    BHaccrete = BlackHoleGrowthRate * mass_ratio /
      (1.0 + pow(280.0 / Gal[merger_centralgal].Vvir, 2.0)) * Gal[merger_centralgal].ColdGas;

    // cannot accrete more gas than is available!
    if(BHaccrete > Gal[merger_centralgal].ColdGas)
      BHaccrete = Gal[merger_centralgal].ColdGas;

    track_BHgrowth(merger_centralgal, p, BHaccrete, time);

    metallicity = get_metallicity(Gal[merger_centralgal].ColdGas, Gal[merger_centralgal].MetalsColdGas);
    Gal[merger_centralgal].BlackHoleMass += BHaccrete;
    Gal[merger_centralgal].ColdGas -= BHaccrete;
    Gal[merger_centralgal].MetalsColdGas -= metallicity * BHaccrete;

    Gal[merger_centralgal].QuasarModeBHaccretionMass += BHaccrete;

    quasar_mode_wind(merger_centralgal, BHaccrete);
  }
}

void grow_black_hole_continuousAccretion(int merger_centralgal, int p, double mass_ratio, double time, double dt)
{
  double BHaccrete = 0., BHaccretionRate = 0., BHaccretionMass = 0.;
  double Luminosity = 0.;
  double metallicity;
  double currentMvir, currentRvir, currentVvir;

  if(Gal[merger_centralgal].ColdGas > 0.0)
  {
    // 1) Determine BH mass to be accreted (this is model dependent!)
    BHaccrete = 0.;
    if(Gal[merger_centralgal].hasJustMerged == 1)
    {
      BHaccrete = BlackHoleGrowthRate * mass_ratio /
                  (1.0 + pow(280.0 / Gal[merger_centralgal].Vvir, 2.0)) * Gal[merger_centralgal].ColdGas;// * (1. + ZZ[Halo[Gal[merger_centralgal].HaloNr].SnapNum]);
    }
    printf("BHaccrete = %e \t ColdaGasToAccrete = %e\n", BHaccrete, Gal[merger_centralgal].ColdGasToAccrete);
    BHaccrete += Gal[merger_centralgal].ColdGasToAccrete;

    // cannot accrete more gas than is available!
    if(BHaccrete > Gal[merger_centralgal].ColdGas)
      BHaccrete = Gal[merger_centralgal].ColdGas;

    // 2) obtain BH accreted in timestep & luminosity (this is model dependent!)

    currentMvir = Gal[merger_centralgal].Mvir;// + Gal[p].deltaMvir * (1.0 - ((double)step + 1.0) / (double)STEPS);
    currentRvir = get_virial_radius_evolving(merger_centralgal, currentMvir);
    currentVvir = get_virial_velocity_evolving(currentMvir, currentRvir);

    printf("time = %e\tdt = %e\tRvir = %e\tMvir=%e\tVvir = %e\t BHmass = %e \t BHaccrete = %e\n", time, dt, currentRvir, currentMvir, currentVvir, Gal[merger_centralgal].BlackHoleMass, BHaccrete);
    if(BHaccrete > 0. && currentMvir > 0.)
    {
      switch(ContinuousAccretionOn)
      {
        case 1:
          accreteOnBH_EddingtonLimited(Gal[merger_centralgal].BlackHoleMass, BHaccrete, QuasarRadEfficiency, dt, &BHaccretionRate, &BHaccretionMass, &Luminosity);
          break;
        case 2:
          accreteOnBH_EddingtonLimited_redshift(Gal[merger_centralgal].BlackHoleMass, BHaccrete, ZZ[Gal[merger_centralgal].SnapNum], QuasarRadEfficiency, dt, &BHaccretionRate, &BHaccretionMass, &Luminosity);
          break;
        case 3:
          accreteOnBH_Hopkins(Gal[merger_centralgal].BlackHoleMass, BHaccrete, QuasarRadEfficiency, dt, &BHaccretionRate, &BHaccretionMass, &Luminosity);
          break;
        case 4:
          accreteOnBH_Ryu(Gal[merger_centralgal].BlackHoleMass, Gal[merger_centralgal].HotGas, currentMvir, currentVvir*UnitVelocity_in_cm_per_s*0.58, MaxInflowVelocity, GasGamma, QuasarRadEfficiency, GasProfileParameter, BHaccrete, dt, &BHaccretionRate, &BHaccretionMass, &Luminosity);
          break;
        case 5:
          accreteOnBH_Park(Gal[merger_centralgal].BlackHoleMass, Gal[merger_centralgal].HotGas, currentMvir, currentVvir*UnitVelocity_in_cm_per_s*0.58, MaxInflowVelocity, GasGamma, GasMu, QuasarRadEfficiency, GasProfileParameter, BHaccrete, dt, &BHaccretionRate, &BHaccretionMass, &Luminosity);
        default:
          break;
      }
    }
    #if 0
    if(ContinuousAccretionOn == 1  && BHaccrete > 0.)
    {
      // Eddington limited case
      fEdd = 1.;
      BHaccretionMassEdd = getBHaccretionMass_EddingtonLimited(Gal[merger_centralgal].BlackHoleMass, fEdd, QuasarRadEfficiency, dt);
      if(BHaccretionMassEdd - Gal[merger_centralgal].BlackHoleMass > BHaccrete * (1. - QuasarRadEfficiency))
      {
        dt_peak = getPeakTime_EddingtonLimited(Gal[merger_centralgal].BlackHoleMass, fEdd, QuasarRadEfficiency, BHaccrete * (1. - QuasarRadEfficiency));
      }else{
        dt_peak = dt;
      }
      BHaccretionRate = getBHaccretionRate_EddingtonLimited(Gal[merger_centralgal].BlackHoleMass, fEdd, QuasarRadEfficiency, dt_peak);
      BHaccretionMass = getBHaccretionMass_EddingtonLimited(Gal[merger_centralgal].BlackHoleMass, fEdd, QuasarRadEfficiency, dt_peak);
      Luminosity = getLuminosity_radEfficient(BHaccretionRate, QuasarRadEfficiency);
    }
    else if(ContinuousAccretionOn == 2  && BHaccrete > 0.)
    {
      if(ZZ[Gal[merger_centralgal].SnapNum] < 3.0)
      {
        fEdd = 0.3 * pow(0.25*(1 + ZZ[Gal[merger_centralgal].SnapNum]), 0.25);
      }else{
        fEdd = 0.3;
      }
      BHaccretionMassEdd = getBHaccretionMass_EddingtonLimited(Gal[merger_centralgal].BlackHoleMass, fEdd, QuasarRadEfficiency, dt);
      if(BHaccretionMassEdd - Gal[merger_centralgal].BlackHoleMass > BHaccrete * (1. - QuasarRadEfficiency))
      {
        dt_peak = getPeakTime_EddingtonLimited(Gal[merger_centralgal].BlackHoleMass, fEdd, QuasarRadEfficiency, BHaccrete * (1. - QuasarRadEfficiency));
      }else{
        dt_peak = dt;
      }
      BHaccretionRate = getBHaccretionRate_EddingtonLimited(Gal[merger_centralgal].BlackHoleMass, fEdd, QuasarRadEfficiency, dt_peak);
      BHaccretionMass = getBHaccretionMass_EddingtonLimited(Gal[merger_centralgal].BlackHoleMass, fEdd, QuasarRadEfficiency, dt_peak);
      Luminosity = getLuminosity_radEfficient(BHaccretionRate, QuasarRadEfficiency);
    }
    else if(ContinuousAccretionOn == 3 && BHaccrete > 0.)
    {
      fEdd = 1.;
      F = 0.7;
      BHaccretionMassEdd = getBHaccretionMass_EddingtonLimited(Gal[merger_centralgal].BlackHoleMass, fEdd, QuasarRadEfficiency, dt);
      if(Gal[merger_centralgal].BlackHoleMass <= 0.){
        Mpeak = BHaccretionMassEdd;
      }else{
        Mpeak = Gal[merger_centralgal].BlackHoleMass + F * BHaccrete * (1. - QuasarRadEfficiency);
      }

      if(BHaccretionMassEdd <= Mpeak)
      {
        BHaccretionRate = getBHaccretionRate_EddingtonLimited(Gal[merger_centralgal].BlackHoleMass, fEdd, QuasarRadEfficiency, dt);
        BHaccretionMass = getBHaccretionMass_EddingtonLimited(Gal[merger_centralgal].BlackHoleMass, fEdd, QuasarRadEfficiency, dt);
        Luminosity = getLuminosity_radEfficient(BHaccretionRate, QuasarRadEfficiency);
      }else{
        dt_peak = getPeakTime_EddingtonLimited(Gal[merger_centralgal].BlackHoleMass, fEdd, QuasarRadEfficiency, F * BHaccrete * (1. - QuasarRadEfficiency));
        BHaccretionRate = getBHaccretionRate_Hopkins(Gal[merger_centralgal].BlackHoleMass, QuasarRadEfficiency, BHaccrete, F, dt-dt_peak);
        BHaccretionMass = getBHaccretionMass_Hopkins(Gal[merger_centralgal].BlackHoleMass, QuasarRadEfficiency, BHaccrete, F, dt-dt_peak);
        Luminosity = getLuminosity_Hopkins(Gal[merger_centralgal].BlackHoleMass, QuasarRadEfficiency, BHaccrete, F, dt-dt_peak);
        // printf("Mpeak = %e\t dt = %e\t dt_peak = %e\n", Mpeak, dt, dt_peak);
        // printf("BHmass = %e\t BHaccrete = %e BHaccretionMass = %e\n\n", Gal[merger_centralgal].BlackHoleMass, BHaccrete, BHaccretionMass-Gal[merger_centralgal].BlackHoleMass);
      }
    }
    else if(ContinuousAccretionOn == 4  && BHaccrete > 0.)
    {
      // Hirschmann 2014 model
      GasRho = GasMu * PROTONMASS * 9.5e5 * (GasTemperature/1.e4)*(GasTemperature/1.e4) / ((Gal[merger_centralgal].BlackHoleMass+1.e-8*Hubble_h)*(Gal[merger_centralgal].BlackHoleMass+1.e-8*Hubble_h)*1.e12);
      printf("\nrho = %e\t%e\t%e\n", GasRho, GasTemperature, Gal[merger_centralgal].BlackHoleMass*1.e6);
      fEdd = compute_fEdd(Gal[merger_centralgal].BlackHoleMass, GasRho, GasTemperature, GasGamma, GasMu,  fEdd, QuasarRadEfficiency, dt);
      printf("fEdd = %e\n", fEdd);
      double BHaccretionMass_tmp = getBHaccretionMass(Gal[merger_centralgal].BlackHoleMass, GasRho, GasTemperature, GasGamma, GasMu, 1., QuasarRadEfficiency, dt);
      if(fEdd < 1.){
          dt_peak = getPeakTime_Bondi(Gal[merger_centralgal].BlackHoleMass, GasRho, GasTemperature, GasGamma, GasMu, BHaccretionMass_tmp);
      }else{
          dt_peak = getPeakTime_EddingtonLimited(Gal[merger_centralgal].BlackHoleMass, fEdd, QuasarRadEfficiency, BHaccretionMass_tmp);
      }
      if(dt_peak < dt)  dt = dt_peak;

      BHaccretionRate = getBHaccretionRate(Gal[merger_centralgal].BlackHoleMass, GasRho, GasTemperature, GasGamma, GasMu, 1., QuasarRadEfficiency, dt);
      BHaccretionMass = getBHaccretionMass(Gal[merger_centralgal].BlackHoleMass, GasRho, GasTemperature, GasGamma, GasMu, 1., QuasarRadEfficiency, dt);
      Luminosity = getLuminosity(BHaccretionRate, QuasarRadEfficiency, Gal[merger_centralgal].BlackHoleMass, GasRho, GasTemperature, GasGamma, GasMu, fEdd);
    }
    else
    {
      BHaccretionRate = 0.;
      BHaccretionMass = 0.;
      Luminosity = 0.;
    }
#endif
    Gal[merger_centralgal].QSOBHaccretionRate = BHaccretionRate;
    Gal[merger_centralgal].QSOBHaccretionMass = BHaccretionMass;
    Gal[merger_centralgal].QSOLuminosity = Luminosity;

    if(BHaccretionMass > BHaccrete)
    {
      Gal[merger_centralgal].ColdGasToAccrete = 0.;
      BHaccretionMass = BHaccrete;
      Gal[merger_centralgal].QSOBHaccretionMass = BHaccretionMass;
    }else{
      Gal[merger_centralgal].ColdGasToAccrete = BHaccrete - BHaccretionMass;
    }
    BHaccrete = BHaccretionMass;

    if(TrackBHgrowthOn == 1 && Gal[merger_centralgal].hasJustMerged == 1)
    {
        track_BHgrowth(merger_centralgal, p, BHaccrete, time);
    }
    
    if(BHaccrete > 0. && Gal[merger_centralgal].BlackHoleMass < 1.e-8*Hubble_h)
    {
        Gal[merger_centralgal].BlackHoleMass = 1.e-8 * Hubble_h;
    }

    // 3) update values
    metallicity = get_metallicity(Gal[merger_centralgal].ColdGas, Gal[merger_centralgal].MetalsColdGas);
    Gal[merger_centralgal].BlackHoleMass += BHaccrete;
    Gal[merger_centralgal].ColdGas -= BHaccrete;
    Gal[merger_centralgal].MetalsColdGas -= metallicity * BHaccrete;

    Gal[merger_centralgal].QuasarModeBHaccretionMass += BHaccrete;

    quasar_mode_wind(merger_centralgal, BHaccrete);
  }
}


void quasar_mode_wind(int gal, float BHaccrete)
{
  float quasar_energy, cold_gas_energy, hot_gas_energy;

  // work out total energies in quasar wind (eta*m*c^2), cold and hot gas (1/2*m*Vvir^2)
  quasar_energy = QuasarModeEfficiency * 0.1 * BHaccrete * (C / UnitVelocity_in_cm_per_s) * (C / UnitVelocity_in_cm_per_s);
  cold_gas_energy = 0.5 * Gal[gal].ColdGas * Gal[gal].Vvir * Gal[gal].Vvir;
  hot_gas_energy = 0.5 * Gal[gal].HotGas * Gal[gal].Vvir * Gal[gal].Vvir;

  // compare quasar wind and cold gas energies and eject cold
  if(quasar_energy > cold_gas_energy)
  {
    Gal[gal].EjectedMass += Gal[gal].ColdGas;
    Gal[gal].MetalsEjectedMass += Gal[gal].MetalsColdGas;

    Gal[gal].ColdGas = 0.0;
    Gal[gal].MetalsColdGas = 0.0;
  }

  // compare quasar wind and cold+hot gas energies and eject hot
  if(quasar_energy > cold_gas_energy + hot_gas_energy)
  {
    Gal[gal].EjectedMass += Gal[gal].HotGas;
    Gal[gal].MetalsEjectedMass += Gal[gal].MetalsHotGas;

    Gal[gal].HotGas = 0.0;
    Gal[gal].MetalsHotGas = 0.0;
  }
}



void add_galaxies_together(int t, int p)
{
  int step;

  Gal[t].ColdGas += Gal[p].ColdGas;
  Gal[t].MetalsColdGas += Gal[p].MetalsColdGas;

  Gal[t].StellarMass += Gal[p].StellarMass;
  Gal[t].MetalsStellarMass += Gal[p].MetalsStellarMass;

  Gal[t].HotGas += Gal[p].HotGas;
  Gal[t].MetalsHotGas += Gal[p].MetalsHotGas;

  Gal[t].EjectedMass += Gal[p].EjectedMass;
  Gal[t].MetalsEjectedMass += Gal[p].MetalsEjectedMass;

  Gal[t].ICS += Gal[p].ICS;
  Gal[t].MetalsICS += Gal[p].MetalsICS;

  Gal[t].BlackHoleMass += Gal[p].BlackHoleMass;

  // add merger to bulge
	Gal[t].BulgeMass += Gal[p].StellarMass;
	Gal[t].MetalsBulgeMass += Gal[p].MetalsStellarMass;

  for(step = 0; step < STEPS; step++)
  {
    Gal[t].SfrBulge[step] += Gal[p].SfrDisk[step] + Gal[p].SfrBulge[step];
    Gal[t].SfrBulgeColdGas[step] += Gal[p].SfrDiskColdGas[step] + Gal[p].SfrBulgeColdGas[step];
    Gal[t].SfrBulgeColdGasMetals[step] += Gal[p].SfrDiskColdGasMetals[step] + Gal[p].SfrBulgeColdGasMetals[step];
  }
}



void make_bulge_from_burst(int p)
{
  int step;

  // generate bulge
  Gal[p].BulgeMass = Gal[p].StellarMass;
  Gal[p].MetalsBulgeMass = Gal[p].MetalsStellarMass;

  // update the star formation rate
  for(step = 0; step < STEPS; step++)
  {
    Gal[p].SfrBulge[step] += Gal[p].SfrDisk[step];
    Gal[p].SfrBulgeColdGas[step] += Gal[p].SfrDiskColdGas[step];
    Gal[p].SfrBulgeColdGasMetals[step] += Gal[p].SfrDiskColdGasMetals[step];
    Gal[p].SfrDisk[step] = 0.0;
    Gal[p].SfrDiskColdGas[step] = 0.0;
    Gal[p].SfrDiskColdGasMetals[step] = 0.0;
  }
}



void collisional_starburst_recipe(double mass_ratio, int merger_centralgal, int centralgal, double time, double dt, int halonr, int mode, int step)
{
  double stars, reheated_mass, ejected_mass, fac, metallicity, eburst;
  double FracZleaveDiskVal;

  // This is the major and minor merger starburst recipe of Somerville et al. 2001.
  // The coefficients in eburst are taken from TJ Cox's PhD thesis and should be more accurate then previous.

  // the bursting fraction
  if(mode == 1)
    eburst = mass_ratio;
  else
    eburst = 0.56 * pow(mass_ratio, 0.7);

  stars = eburst * Gal[merger_centralgal].ColdGas;
  if(stars < 0.0)
    stars = 0.0;

  // this bursting results in SN feedback on the cold/hot gas
  if(SupernovaRecipeOn == 1)
    reheated_mass = FeedbackReheatingEpsilon * stars;
  else
    reheated_mass = 0.0;

	assert(reheated_mass >= 0.0);

  // can't use more cold gas than is available! so balance SF and feedback
  if((stars + reheated_mass) > Gal[merger_centralgal].ColdGas)
  {
    fac = Gal[merger_centralgal].ColdGas / (stars + reheated_mass);
    stars *= fac;
    reheated_mass *= fac;
  }

  // determine ejection
  if(SupernovaRecipeOn == 1)
  {
    if(Gal[centralgal].Vvir > 0.0)
			ejected_mass =
				(FeedbackEjectionEfficiency * (EtaSNcode * EnergySNcode) / (Gal[centralgal].Vvir * Gal[centralgal].Vvir) -
					FeedbackReheatingEpsilon) * stars;
		else
			ejected_mass = 0.0;

    if(ejected_mass < 0.0)
      ejected_mass = 0.0;
  }
  else
    ejected_mass = 0.0;

  // starbursts add to the bulge
  Gal[merger_centralgal].SfrBulge[step] += stars / dt;
  Gal[merger_centralgal].SfrBulgeColdGas[step] += Gal[merger_centralgal].ColdGas;
  Gal[merger_centralgal].SfrBulgeColdGasMetals[step] += Gal[merger_centralgal].MetalsColdGas;

  metallicity = get_metallicity(Gal[merger_centralgal].ColdGas, Gal[merger_centralgal].MetalsColdGas);
  update_from_star_formation(merger_centralgal, stars, metallicity);

  Gal[merger_centralgal].BulgeMass += (1 - RecycleFraction) * stars;
  Gal[merger_centralgal].MetalsBulgeMass += metallicity * (1 - RecycleFraction) * stars;

  // recompute the metallicity of the cold phase
  metallicity = get_metallicity(Gal[merger_centralgal].ColdGas, Gal[merger_centralgal].MetalsColdGas);

  // update from feedback
  update_from_feedback(merger_centralgal, centralgal, reheated_mass, ejected_mass, metallicity);

  // check for disk instability
  if(DiskInstabilityOn && mode == 0)
    if(mass_ratio < ThreshMajorMerger)
    {
      check_disk_instability(merger_centralgal, centralgal, halonr, time, dt, step);
    }

  // formation of new metals - instantaneous recycling approximation - only SNII
  if(Gal[merger_centralgal].ColdGas > 1e-8 && mass_ratio < ThreshMajorMerger)
  {
    FracZleaveDiskVal = FracZleaveDisk * exp(-1.0 * Gal[centralgal].Mvir / 30.0);  // Krumholz & Dekel 2011 Eq. 22
    Gal[merger_centralgal].MetalsColdGas += Yield * (1.0 - FracZleaveDiskVal) * stars;
    Gal[centralgal].MetalsHotGas += Yield * FracZleaveDiskVal * stars;
    // Gal[centralgal].MetalsEjectedMass += Yield * FracZleaveDiskVal * stars;
  }
  else
    Gal[centralgal].MetalsHotGas += Yield * stars;
    // Gal[centralgal].MetalsEjectedMass += Yield * stars;
}



void disrupt_satellite_to_ICS(int centralgal, int gal)
{
  Gal[centralgal].HotGas += Gal[gal].ColdGas + Gal[gal].HotGas;
  Gal[centralgal].MetalsHotGas += Gal[gal].MetalsColdGas + Gal[gal].MetalsHotGas;

  Gal[centralgal].EjectedMass += Gal[gal].EjectedMass;
  Gal[centralgal].MetalsEjectedMass += Gal[gal].MetalsEjectedMass;

  Gal[centralgal].ICS += Gal[gal].ICS;
  Gal[centralgal].MetalsICS += Gal[gal].MetalsICS;

  Gal[centralgal].ICS += Gal[gal].StellarMass;
  Gal[centralgal].MetalsICS += Gal[gal].MetalsStellarMass;

  // what should we do with the disrupted satellite BH?

  Gal[gal].mergeType = 4;  // mark as disruption to the ICS


}
