#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include "core_simulation.h"

#define ABORT(sigterm)                                                  \
do {                                                                \
  printf("Error in file: %s\tfunc: %s\tline: %i\n", __FILE__, __FUNCTION__, __LINE__); \
  myexit(sigterm);                                                \
} while(0)

#define  STEPS 10         // Number of integration intervals between two snapshots
#define  MAXGALFAC 1
#define  ALLOCPARAMETER 10.0
#define  MAX_NODE_NAME_LEN 50
#define  ABSOLUTEMAXSNAPS 1000

#define  PI          3.14159265359

#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.3806e-16
#define  GAS_CONST   8.31425e7
#define  C           2.9979e10
#define  PLANCK      6.6262e-27
#define  CM_PER_MPC  3.085678e24
#define  CM_PER_KM   1.e5
#define  PROTONMASS  1.6726e-24
#define  THOMSON_CS  6.652459e-25
#define  HUBBLE      3.2407789e-18   /* in h/sec */

#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7

#define  MERGER_NUM  10    // number of mergers that are stored for each galaxy
#define  TEDD_YR     4.5e9

// This structure contains the properties that are output
struct GALAXY_OUTPUT
{
  int   SnapNum;
  int   Type;

  long long   GalaxyIndex;
  long long   CentralGalaxyIndex;
  int   SAGEHaloIndex;
  int   SAGETreeIndex;
  long long   SimulationHaloIndex;

  int   mergeType;  //0=none; 1=minor merger; 2=major merger; 3=disk instability; 4=disrupt to ICS
  int   mergeIntoID;
  int   mergeIntoSnapNum;
  float dT;

  // (sub)halo properties
  float Pos[3];
  float Vel[3];
  float Spin[3];
  int   Len;
  float Mvir;
  float CentralMvir;
  float Rvir;
  float Vvir;
  float Vmax;
  float VelDisp;

  // baryonic reservoirs
  float ColdGas;
  float StellarMass;
  float BulgeMass;
  float HotGas;
  float EjectedMass;
  float BlackHoleMass;
  float ICS;

  // metals
  float MetalsColdGas;
  float MetalsStellarMass;
  float MetalsBulgeMass;
  float MetalsHotGas;
  float MetalsEjectedMass;
  float MetalsICS;

  // to calculate magnitudes
  float SfrDisk;
  float SfrBulge;
  float SfrDiskZ;
  float SfrBulgeZ;

  // misc
  float DiskScaleRadius;
  float Cooling;
  float Heating;
  float QuasarModeBHaccretionMass;
  float TimeOfLastMajorMerger;
  float TimeOfLastMinorMerger;
  float OutflowRate;

  // infall properties
  float infallMvir;
  float infallVvir;
  float infallVmax;

  // quasar luminosity parameters
  float QSOBHaccrete[MERGER_NUM];
  float QSOmergeAge[MERGER_NUM];
  float QSOmergeTime[MERGER_NUM];
  float QSOBH[MERGER_NUM];

  int MergSnap;
  int QSOmergSnap[MERGER_NUM];
  int MergNum;
  int QSOmergeType[MERGER_NUM];

  float QSOBHaccretionRate;
  float QSOBHaccretionMass;
  float QSOLuminosity;
  float ColdGasToAccrete;
};


// This structure contains the properties used within the code
struct GALAXY
{
  int   SnapNum;
  int   Type;

  int   GalaxyNr;
  int   CentralGal;
  int   HaloNr;
  long long MostBoundID;

  int   mergeType;  //0=none; 1=minor merger; 2=major merger; 3=disk instability; 4=disrupt to ICS
  int   mergeIntoID;
  int   mergeIntoSnapNum;
  float dT;

  // (sub)halo properties
  float Pos[3];
  float Vel[3];
  int   Len;
  float Mvir;
  float deltaMvir;
  float CentralMvir;
  float Rvir;
  float Vvir;
  float Vmax;

  // baryonic reservoirs
  float ColdGas;
  float StellarMass;
  float BulgeMass;
  float HotGas;
  float EjectedMass;
  float BlackHoleMass;
  float ICS;

  // metals
  float MetalsColdGas;
  float MetalsStellarMass;
  float MetalsBulgeMass;
  float MetalsHotGas;
  float MetalsEjectedMass;
  float MetalsICS;

  // to calculate magnitudes
  float SfrDisk[STEPS];
  float SfrBulge[STEPS];
  float SfrDiskColdGas[STEPS];
  float SfrDiskColdGasMetals[STEPS];
  float SfrBulgeColdGas[STEPS];
  float SfrBulgeColdGasMetals[STEPS];

  // misc
  float DiskScaleRadius;
  // float deltaDiskScaleRadius;
  float MergTime;
  float MergTimeInit;
  double Cooling;
  double Heating;
  float r_heat;
  float QuasarModeBHaccretionMass;
  float TimeOfLastMajorMerger;
  float TimeOfLastMinorMerger;
  float OutflowRate;
  float TotalSatelliteBaryons;

  // infall properties
  float infallMvir;
  float infallVvir;
  float infallVmax;

  // quasar luminosity parameters
  float QSOBHaccrete[MERGER_NUM];
  float QSOmergeAge[MERGER_NUM];
  float QSOmergeTime[MERGER_NUM];
  float QSOBH[MERGER_NUM];

  int MergSnap;
  int QSOmergSnap[MERGER_NUM];
  int MergNum;
  int QSOmergeType[MERGER_NUM];

  // quasar accretion model
  int hasJustMerged;
  float QSOBHaccretionRate;
  float QSOBHaccretionMass;
  float QSOLuminosity;
  float ColdGasToAccrete;
}
*Gal, *HaloGal;


// auxiliary halo data
struct halo_aux_data
{
  int DoneFlag;
  int HaloFlag;
  int NGalaxies;
  int FirstGalaxy;
}
*HaloAux;


extern int    FirstFile;    // first and last file for processing
extern int    LastFile;

extern int    Ntrees;      // number of trees in current file
extern int    NumGals;     // Total number of galaxies stored for current tree
extern int    MaxGals;     // Maximum number of galaxies allowed for current tree
extern int    FoF_MaxGals;

extern int    GalaxyCounter;     // unique galaxy ID for main progenitor line in tree

extern int    LastSnapShotNr;

extern char   OutputDir[512];
extern char   FileNameGalaxies[512];
extern char   TreeName[512];
extern char   SimulationDir[512];
extern char   FileWithSnapList[512];

extern int    TotHalos;
extern int    TotGalaxies[ABSOLUTEMAXSNAPS];
extern int    *TreeNgals[ABSOLUTEMAXSNAPS];

extern int    *FirstHaloInSnap;

extern int    *TreeNHalos;
extern int    *TreeFirstHalo;

#ifdef MPI
extern int ThisTask, NTask, nodeNameLen;
extern char *ThisNode;
#endif

extern double Omega;
extern double OmegaLambda;
extern double PartMass;
extern double Hubble_h;
extern double EnergySNcode, EnergySN;
extern double EtaSNcode, EtaSN;

// recipe flags
extern int    ReionizationOn;
extern int    SupernovaRecipeOn;
extern int    DiskInstabilityOn;
extern int    AGNrecipeOn;
extern int    SFprescription;
extern int    TrackBHgrowthOn;
extern int    ContinuousAccretionOn;

// recipe parameters
extern double RecycleFraction;
extern double Yield;
extern double FracZleaveDisk;
extern double ReIncorporationFactor;
extern double ThreshMajorMerger;
extern double BaryonFrac;
extern double SfrEfficiency;
extern double FeedbackReheatingEpsilon;
extern double FeedbackEjectionEfficiency;
extern double RadioModeEfficiency;
extern double QuasarModeEfficiency;
extern double BlackHoleGrowthRate;
extern double QuasarRadEfficiency;
extern double GasTemperature;
extern double GasGamma;
extern double GasMu;
extern double Reionization_z0;
extern double Reionization_zr;
extern double ThresholdSatDisruption;

extern double UnitLength_in_cm,
  UnitTime_in_s,
  UnitVelocity_in_cm_per_s,
  UnitMass_in_g,
  RhoCrit,
  UnitPressure_in_cgs,
  UnitDensity_in_cgs,
  UnitCoolingRate_in_cgs,
  UnitEnergy_in_cgs,
  UnitTime_in_Megayears,
  G,
  Hubble,
  a0, ar;

extern int    ListOutputSnaps[ABSOLUTEMAXSNAPS];

extern double ZZ[ABSOLUTEMAXSNAPS];
extern double AA[ABSOLUTEMAXSNAPS];
extern double *Age;

extern int    MAXSNAPS;
extern int    NOUT;
extern int    Snaplistlen;

extern gsl_rng *random_generator;

extern int TreeID;
extern int FileNum;


#endif  // #ifndef ALLVARS_H
