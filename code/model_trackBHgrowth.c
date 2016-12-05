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
    int tmpQSOmergeType;
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
        tmpQSOmergeType = Gal[p].Type;

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
        Gal[merger_centralgal].QSOmergeType[entryToRemove_index] = tmpQSOmergeType;

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
                        tmpQSOmergeType = Gal[merger_centralgal].QSOmergeType[i];

                        Gal[merger_centralgal].QSOBHaccrete[i] = Gal[merger_centralgal].QSOBHaccrete[j];
                        Gal[merger_centralgal].QSOmergeAge[i] = Gal[merger_centralgal].QSOmergeAge[j];
                        Gal[merger_centralgal].QSOmergeTime[i] = Gal[merger_centralgal].QSOmergeTime[j];
                        Gal[merger_centralgal].QSOmergSnap[i] = Gal[merger_centralgal].QSOmergSnap[j];
                        Gal[merger_centralgal].QSOBH[i] = Gal[merger_centralgal].QSOBH[j];
                        Gal[merger_centralgal].QSOmergeType[i] = Gal[merger_centralgal].QSOmergeType[j];

                        Gal[merger_centralgal].QSOBHaccrete[j] = tmpQSOBHaccrete;
                        Gal[merger_centralgal].QSOmergeAge[j] = tmpQSOmergeAge;
                        Gal[merger_centralgal].QSOmergeTime[j] = tmpQSOmergeTime;
                        Gal[merger_centralgal].QSOmergSnap[j] = tmpQSOmergeSnap;
                        Gal[merger_centralgal].QSOBH[j] = tmpQSOBH;
                        Gal[merger_centralgal].QSOmergeType[i] = tmpQSOmergeType;
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
                        tmpQSOmergeType = Gal[merger_centralgal].QSOmergeType[i];

                        Gal[merger_centralgal].QSOBHaccrete[i] = Gal[merger_centralgal].QSOBHaccrete[j];
                        Gal[merger_centralgal].QSOmergeAge[i] = Gal[merger_centralgal].QSOmergeAge[j];
                        Gal[merger_centralgal].QSOmergeTime[i] = Gal[merger_centralgal].QSOmergeTime[j];
                        Gal[merger_centralgal].QSOmergSnap[i] = Gal[merger_centralgal].QSOmergSnap[j];
                        Gal[merger_centralgal].QSOBH[i] = Gal[merger_centralgal].QSOBH[j];
                        Gal[merger_centralgal].QSOmergeType[i] = Gal[merger_centralgal].QSOmergeType[j];

                        Gal[merger_centralgal].QSOBHaccrete[j] = tmpQSOBHaccrete;
                        Gal[merger_centralgal].QSOmergeAge[j] = tmpQSOmergeAge;
                        Gal[merger_centralgal].QSOmergeTime[j] = tmpQSOmergeTime;
                        Gal[merger_centralgal].QSOmergSnap[j] = tmpQSOmergeSnap;
                        Gal[merger_centralgal].QSOBH[j] = tmpQSOBH;
                        Gal[merger_centralgal].QSOmergeType[i] = tmpQSOmergeType;
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
                        tmpQSOmergeType = Gal[merger_centralgal].QSOmergeType[i];

                        Gal[merger_centralgal].QSOBHaccrete[i] = Gal[merger_centralgal].QSOBHaccrete[j];
                        Gal[merger_centralgal].QSOmergeAge[i] = Gal[merger_centralgal].QSOmergeAge[j];
                        Gal[merger_centralgal].QSOmergeTime[i] = Gal[merger_centralgal].QSOmergeTime[j];
                        Gal[merger_centralgal].QSOmergSnap[i] = Gal[merger_centralgal].QSOmergSnap[j];
                        Gal[merger_centralgal].QSOBH[i] = Gal[merger_centralgal].QSOBH[j];
                        Gal[merger_centralgal].QSOmergeType[i] = Gal[merger_centralgal].QSOmergeType[j];

                        Gal[merger_centralgal].QSOBHaccrete[j] = tmpQSOBHaccrete;
                        Gal[merger_centralgal].QSOmergeAge[j] = tmpQSOmergeAge;
                        Gal[merger_centralgal].QSOmergeTime[j] = tmpQSOmergeTime;
                        Gal[merger_centralgal].QSOmergSnap[j] = tmpQSOmergeSnap;
                        Gal[merger_centralgal].QSOBH[j] = tmpQSOBH;
                        Gal[merger_centralgal].QSOmergeType[i] = tmpQSOmergeType;
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
