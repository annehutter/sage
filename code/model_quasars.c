#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"

void track_BHgrowth(int merger_centralgal, int p, double BHaccrete, double time)
{
    int i,j;
    float tmpQSOBHaccrete, tmpQSOmergeAge, tmpQSOmergeTime, tmpQSOmergeSnap, tmpQSOBH;
    int entryToRemove_index = 0;
    float entryToRemove;
    float tmpi, tmpj;

#ifdef DEBUG
    printf("p = %d\ttime = %e\t Type = %d\tSnapNum = %d\tMergSnap = %d\n", p, time, Gal[p].Type, Gal[p].SnapNum, Gal[p].MergSnap);
#endif

    // if(Gal[p].SnapNum != Gal[p].MergSnap)
    if(Gal[p].Type != 0)
    {
        Gal[merger_centralgal].MergNum++;

        tmpQSOmergeAge = time;
        tmpQSOmergeTime =  Gal[p].MergTimeInit - Gal[p].MergTime;
        if(tmpQSOmergeTime <= 0.)
            tmpQSOmergeTime = Gal[merger_centralgal].dT / STEPS;
        tmpQSOBHaccrete = BHaccrete;
        tmpQSOmergeSnap = Gal[p].MergSnap;
        tmpQSOBH = Gal[merger_centralgal].BlackHoleMass;

        if(TrackBHgrowthOn == 1)
        {
            // find minimum for considering BHaccrete
            entryToRemove = 1.e10;
            entryToRemove_index = MERGER_NUM;
            for(j = 0; j < MERGER_NUM; j++)
            {
                if(Gal[merger_centralgal].QSOBHaccrete[j] < entryToRemove)
                {
                    entryToRemove = Gal[merger_centralgal].QSOBHaccrete[j];
                    entryToRemove_index = j;
                }
            }
        }else if(TrackBHgrowthOn ==2){
            // find maximum for considering time
            float timeFromMerger;
            entryToRemove = 0.;
            entryToRemove_index = 0;
            for(j = 0; j < MERGER_NUM; j++)
            {
                if(Gal[merger_centralgal].QSOmergeAge[j] != 0.)
                    timeFromMerger = Gal[merger_centralgal].QSOmergeAge[j];
                else
                    timeFromMerger = 1.;
                if(timeFromMerger > entryToRemove)
                {
                    entryToRemove = timeFromMerger;
                    entryToRemove_index = j;
                }
            }
#ifdef DEBUG
            printf("entry to remove: %d\n", entryToRemove_index);
#endif
        }else if(TrackBHgrowthOn == 3){
            // find minimum for considering BHaccrete/time
            float BHaccreteOvertimeFromMerger;
            entryToRemove = 1.e10;
            entryToRemove_index = MERGER_NUM;
            for(j = 0; j < MERGER_NUM; j++)
            {
                if(Gal[merger_centralgal].QSOmergeAge[j] != (float)time)
                {
                    BHaccreteOvertimeFromMerger = Gal[merger_centralgal].QSOBHaccrete[j] / (Gal[merger_centralgal].QSOmergeAge[j] - time);
                }else{
                    BHaccreteOvertimeFromMerger = 1.e10;
                }
                if(BHaccreteOvertimeFromMerger < entryToRemove)
                {
                    entryToRemove = BHaccreteOvertimeFromMerger;
                    entryToRemove_index = j;
                }
            }
#ifdef DEBUG
                    printf("entry to remove: %d\n", entryToRemove_index);
#endif
        }

        // write current value into previous minimum value
        Gal[merger_centralgal].QSOBHaccrete[entryToRemove_index] = tmpQSOBHaccrete;
        Gal[merger_centralgal].QSOmergeAge[entryToRemove_index] = tmpQSOmergeAge;
        Gal[merger_centralgal].QSOmergeTime[entryToRemove_index] = tmpQSOmergeTime;
        Gal[merger_centralgal].QSOmergSnap[entryToRemove_index] = tmpQSOmergeSnap;
        Gal[merger_centralgal].QSOBH[entryToRemove_index] = tmpQSOBH;

        if(TrackBHgrowthOn == 1)
        {
            // sort array
            for(i = 0; i < MERGER_NUM; i++)
            {
                for(j = i+1; j < MERGER_NUM; ++j)
                {
                    if(Gal[merger_centralgal].QSOBHaccrete[i]<Gal[merger_centralgal].QSOBHaccrete[j])
                    {
                        tmpQSOBHaccrete = Gal[merger_centralgal].QSOBHaccrete[i];
                        tmpQSOmergeAge = Gal[merger_centralgal].QSOmergeAge[i];
                        tmpQSOmergeTime = Gal[merger_centralgal].QSOmergeTime[i];
                        tmpQSOmergeSnap = Gal[merger_centralgal].QSOmergSnap[i];
                        tmpQSOBH = Gal[merger_centralgal].QSOBH[i];

                        Gal[merger_centralgal].QSOBHaccrete[i] = Gal[merger_centralgal].QSOBHaccrete[j];
                        Gal[merger_centralgal].QSOmergeAge[i] = Gal[merger_centralgal].QSOmergeAge[j];
                        Gal[merger_centralgal].QSOmergeTime[i] = Gal[merger_centralgal].QSOmergeTime[j];
                        Gal[merger_centralgal].QSOmergSnap[i] = Gal[merger_centralgal].QSOmergSnap[j];
                        Gal[merger_centralgal].QSOBH[i] = Gal[merger_centralgal].QSOBH[j];

                        Gal[merger_centralgal].QSOBHaccrete[j] = tmpQSOBHaccrete;
                        Gal[merger_centralgal].QSOmergeAge[j] = tmpQSOmergeAge;
                        Gal[merger_centralgal].QSOmergeTime[j] = tmpQSOmergeTime;
                        Gal[merger_centralgal].QSOmergSnap[j] = tmpQSOmergeSnap;
                        Gal[merger_centralgal].QSOBH[j] = tmpQSOBH;
                    }
                }
            }
        }else if(TrackBHgrowthOn == 2){
            // sort array
            for(i = 0; i < MERGER_NUM; i++)
            {
                if(Gal[merger_centralgal].QSOmergeAge[i] != 0)
                    tmpi = Gal[merger_centralgal].QSOmergeAge[i];
                else
                    tmpi = 1.;
                for(j = i+1; j < MERGER_NUM; ++j)
                {
                    if(Gal[merger_centralgal].QSOmergeAge[j] != 0.)
                        tmpj = Gal[merger_centralgal].QSOmergeAge[j];
                    else
                        tmpj = 1.;
                    if(tmpi > tmpj)
                    {
                        tmpQSOBHaccrete = Gal[merger_centralgal].QSOBHaccrete[i];
                        tmpQSOmergeAge = Gal[merger_centralgal].QSOmergeAge[i];
                        tmpQSOmergeTime = Gal[merger_centralgal].QSOmergeTime[i];
                        tmpQSOmergeSnap = Gal[merger_centralgal].QSOmergSnap[i];
                        tmpQSOBH = Gal[merger_centralgal].QSOBH[i];

                        Gal[merger_centralgal].QSOBHaccrete[i] = Gal[merger_centralgal].QSOBHaccrete[j];
                        Gal[merger_centralgal].QSOmergeAge[i] = Gal[merger_centralgal].QSOmergeAge[j];
                        Gal[merger_centralgal].QSOmergeTime[i] = Gal[merger_centralgal].QSOmergeTime[j];
                        Gal[merger_centralgal].QSOmergSnap[i] = Gal[merger_centralgal].QSOmergSnap[j];
                        Gal[merger_centralgal].QSOBH[i] = Gal[merger_centralgal].QSOBH[j];

                        Gal[merger_centralgal].QSOBHaccrete[j] = tmpQSOBHaccrete;
                        Gal[merger_centralgal].QSOmergeAge[j] = tmpQSOmergeAge;
                        Gal[merger_centralgal].QSOmergeTime[j] = tmpQSOmergeTime;
                        Gal[merger_centralgal].QSOmergSnap[j] = tmpQSOmergeSnap;
                        Gal[merger_centralgal].QSOBH[j] = tmpQSOBH;
                    }
                }
            }
        }else if(TrackBHgrowthOn == 3){
        // sort array
            for(i = 0; i < MERGER_NUM; i++)
            {
                tmpi = Gal[merger_centralgal].QSOBHaccrete[i] / (Gal[merger_centralgal].QSOmergeAge[i] -  (float)time);
                if(Gal[merger_centralgal].QSOmergeAge[i] == (float)time) tmpi = 1.e10;

                for(j = i+1; j < MERGER_NUM; ++j)
                {
                    tmpj = Gal[merger_centralgal].QSOBHaccrete[j] / (Gal[merger_centralgal].QSOmergeAge[j] -  (float)time);
                    if(Gal[merger_centralgal].QSOmergeAge[j] == (float)time) tmpi = 1.e10;

                    if(tmpi < tmpj)
                    {
                        tmpQSOBHaccrete = Gal[merger_centralgal].QSOBHaccrete[i];
                        tmpQSOmergeAge = Gal[merger_centralgal].QSOmergeAge[i];
                        tmpQSOmergeTime = Gal[merger_centralgal].QSOmergeTime[i];
                        tmpQSOmergeSnap = Gal[merger_centralgal].QSOmergSnap[i];
                        tmpQSOBH = Gal[merger_centralgal].QSOBH[i];

                        Gal[merger_centralgal].QSOBHaccrete[i] = Gal[merger_centralgal].QSOBHaccrete[j];
                        Gal[merger_centralgal].QSOmergeAge[i] = Gal[merger_centralgal].QSOmergeAge[j];
                        Gal[merger_centralgal].QSOmergeTime[i] = Gal[merger_centralgal].QSOmergeTime[j];
                        Gal[merger_centralgal].QSOmergSnap[i] = Gal[merger_centralgal].QSOmergSnap[j];
                        Gal[merger_centralgal].QSOBH[i] = Gal[merger_centralgal].QSOBH[j];

                        Gal[merger_centralgal].QSOBHaccrete[j] = tmpQSOBHaccrete;
                        Gal[merger_centralgal].QSOmergeAge[j] = tmpQSOmergeAge;
                        Gal[merger_centralgal].QSOmergeTime[j] = tmpQSOmergeTime;
                        Gal[merger_centralgal].QSOmergSnap[j] = tmpQSOmergeSnap;
                        Gal[merger_centralgal].QSOBH[j] = tmpQSOBH;
                    }
                }
            }
        }

#ifdef DEBUG
        for(j = 0; j < MERGER_NUM; j++)
            printf("BHaccrete = %e\t", Gal[merger_centralgal].QSOBHaccrete[j]);
        printf("\n");
        for(j = 0; j < MERGER_NUM; j++)
            printf("BH = %e\t\t", Gal[merger_centralgal].QSOBH[j]);
        printf("\n");
        for(j = 0; j < MERGER_NUM; j++)
            printf("mergeAge  = %e\t", Gal[merger_centralgal].QSOmergeAge[j]);
        printf("\n");
        for(j = 0; j < MERGER_NUM; j++)
        {
            tmpj = Gal[merger_centralgal].QSOBHaccrete[j] / (Gal[merger_centralgal].QSOmergeAge[j] -  (float)time);
            if(Gal[merger_centralgal].QSOmergeAge[j] == (float)time) tmpi = 1.e10;
            printf("BHaccOvTime  = %e\t", tmpj);
        }
        printf("\n\n");
#endif
    }
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
