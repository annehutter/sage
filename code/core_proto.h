#include "core_allvars.h"

size_t myfread(void  *ptr,  size_t  size,  size_t  nmemb,  FILE *stream);
size_t myfwrite(void  *ptr,  size_t  size,  size_t  nmemb,  FILE *stream);
int myfseek(FILE *stream, long offset, int whence);

void construct_galaxies(int halonr, int tree);
void evolve_galaxies(int halonr, int ngal, int tree);
int  join_galaxies_of_progenitors(int halonr, int nstart);
void init(void);
void set_units(void);

void load_tree_table(int filenr);
void load_tree(int filenr, int nr);
void save_galaxies(int filenr, int tree);

void prepare_galaxy_for_output(int filenr, int tree, struct GALAXY *g, struct GALAXY_OUTPUT *o);

void free_galaxies_and_tree(void);
void free_tree_table(void);
void print_allocated(void);

void read_parameter_file(char *fname);
void *mymalloc(size_t n);
void *myrealloc(void *p, size_t n);
void myfree(void *p);
void myexit(int signum);

void finalize_galaxy_file(int filenr);

void starformation_and_feedback(int p, int centralgal, double time, double dt, int halonr, int step);
void add_galaxies_together(int t, int p);
void init_galaxy(int p, int halonr);
double infall_recipe(int centralgal, int ngal, double Zcurr);
void add_infall_to_hot(int centralgal, double infallingGas);
double cooling_recipe(int centralgal, double dt);
void cool_gas_onto_galaxy(int centralgal, double coolingGas);
void reincorporate_gas(int centralgal, double dt);
double estimate_merging_time(int prog, int mother_halo, int ngal);
void deal_with_galaxy_merger(int p, int merger_centralgal, int centralgal, double time, double dt, int halonr, int step);
double dmax(double x, double y);
double do_reionization(int centralgal, double Zcurr);
double do_AGN_heating(double coolingGas, int centralgal, double dt, double x, double rcool);
void collisional_starburst_recipe(double mass_ratio, int merger_centralgal, int centralgal, double time, double dt, int halonr, int mode, int step);
void update_from_star_formation(int p, double stars, double metallicity);
void update_from_feedback(int p, int centralgal, double reheated_mass, double ejected_mass, double metallicity);
void make_bulge_from_burst(int p);
void grow_black_hole_trackBHgrowth(int merger_centralgal, int p, double mass_ratio, double time);
void grow_black_hole(int merger_centralgal, double mass_ratio);
void grow_black_hole_continuousAccretion(int merger_centralgal, int p, double mass_ratio, double time, double dt);
void check_disk_instability(int p, int centralgal, int halonr, double time, double dt, int step);

void track_BHgrowth(int merger_centralgal, int p, double BHaccrete, double time);

void accreteOnBH_EddingtonLimited(double BHmass, double BHaccrete, double rad_efficiency, double dt, double *BHaccretionRate, double *BHaccretionMass, double *Luminosity);
void accreteOnBH_EddingtonLimited_redshift(double BHmass, double BHaccrete, double redshift, double rad_efficiency, double dt, double *BHaccretionRate, double *BHaccretionMass, double *Luminosity);
void accreteOnBH_Hopkins(double BHmass, double BHaccrete, double rad_efficiency, double dt, double *BHaccretionRate, double *BHaccretionMass, double *Luminosity);
void accreteOnBH_Ryu(double BHmass, double GasMass, double Mvir, double Vvir, double cs, double cs_inflow, double gamma, double rad_efficiency, double parameter, double BHaccrete, double dt, double *BHaccretionRate, double *BHaccretionMass, double *Luminosity, int *type_acc_global);
void accreteOnBH_Ryu_MCF(double BHmass, double GasMass, double Mvir, double Vvir, double cs, double cs_inflow, double gamma, double mu, double logZ, double rad_efficiency, double parameter, double BHaccrete, double dt, double *BHaccretionRate, double *BHaccretionMass, double *Luminosity, int *type_acc_global);
void accreteOnBH_Park(double BHmass, double GasMass, double Mvir, double Vvir, double cs, double cs_inflow, double gamma, double mu, double rad_efficiency, double parameter, double BHaccrete, double dt, double *BHaccretionRate, double *BHaccretionMass, double *Luminosity, int *type_acc_global);
void accreteOnBH_Park_MCF(double BHmass, double GasMass, double Mvir, double Vvir, double cs, double cs_inflow, double gamma, double mu, double logZ, double rad_efficiency, double parameter, double BHaccrete, double dt, double *BHaccretionRate, double *BHaccretionMass, double *Luminosity, int *type_acc_global);


double getSubEddDensity(double BHmass, double cs, double gamma, double mu, double rad_efficiency);
double getEddingtonDensity(double BHmass, double cs, double gamma, double mu, double meanPhotEnergy, double rad_efficiency);
double getCriticalDensity(double BHmass, double cs, double gamma, double mu, double meanPhotEnergy);
double getSuperEddDensity(double BHmass, double cs, double gamma, double mu, double meanPhotEnergy, double rad_efficiency);
double getRhoAtBondiRadius(double BHmass, double GasMass, double Mvir, double Vvir, double cs, double gamma, double parameter);
double getDensityAtBondiRadius(double BHmass, double GasMass, double Mvir, double Vvir, double cs, double gamma, double mu, double parameter);
double getRhoAtBondiRadius_MCF(double BHmass, double Vvir, double cs, double gamma, double mu, double logZ);
double getDensityAtBondiRadius_MCF(double BHmass, double Vvir, double cs, double gamma, double mu, double logZ);
double getBondiRate(double BHmass, double GasMass, double Mvir, double Vvir, double cs, double gamma, double parameter);
double getBondiRate_MCF(double BHmass, double Vvir, double cs, double gamma, double mu, double logZ);
double getEddRate(double BHmass);
double getInflowRate(double cs);
double getBondiMass(double BHmass, double GasMass, double Mvir, double parameter);
double getBHaccretionRate_EddingtonLimited(double BHmass, double fEdd, double rad_efficiency, double dt);
double getBHaccretionMass_EddingtonLimited(double BHmass, double fEdd, double rad_efficiency, double dt);
double getPeakTime_EddingtonLimited(double BHmass, double fEdd, double rad_efficiency, double BHaccrete);
double getBHaccretionRate_Bondi(double BHmass, double GasMass, double Mvir, double Vvir, double cs, double gamma, double parameter, double rad_efficiency, double lambdaRad);
double getBHaccretionMass_Bondi(double BHmass, double GasMass, double Mvir, double Vvir, double cs, double gamma, double parameter, double rad_efficiency, double lambdaRad, double dt);
double getPeakTime_Bondi(double BHmass, double GasMass, double Mvir, double Vvir, double cs, double gamma, double parameter, double rad_efficiency, double lambdaRad, double BHaccrete);
double getBHaccretionRate_Bondi_MCF(double BHmass, double Vvir, double cs, double gamma, double mu, double logZ, double rad_efficiency, double lambdaRad, double dt);
double getBHaccretionMass_Bondi_MCF(double BHmass, double Vvir, double cs, double gamma, double mu, double logZ, double rad_efficiency, double lambdaRad, double dt);
double getPeakTime_Bondi_MCF(double BHmass, double Vvir, double cs, double gamma, double mu, double logZ, double rad_efficiency, double lambdaRad, double BHaccrete);
double getBHaccretionRate_Inflow(double cs, double rad_efficiency);
double getBHaccretionMass_Inflow(double BHmass, double cs, double rad_efficiency, double dt);
double getPeakTime_Inflow(double cs, double rad_efficiency, double BHaccrete);
double getMeanPhotEnergy(double spectralIndex);
double getLambda_rad(double BHmass, double density, double cs, double gamma, double mu, double meanPhotEnergy);
double getFduty(double density, double density_crit, double cs, double gamma, double mu);
double getTcycle(double BHmass, double density, double density_crit, double meanPhotEnergy, double rad_efficiency);
double getBHmassBoundary(double BHmass, double GasMass, double Mvir, double Vvir, double cs, double gamma, double mu, double parameter, double BoundaryValue);
double getBHmassBoundary_MCF(double BHmass, double GasMass, double Mvir, double Vvir, double cs, double gamma, double mu, double parameter, double BoundaryValue);
double getLmaxDivLEdd(double BHmass, double density, double density_crit, double cs, double gamma, double mu, double meanPhotEnergy, double rad_efficiency);
double getMaccr(double BHaccretionRate, double BHmass);
double getEddLuminosity(double BHmass);
double getLuminosity_radEfficient(double BHaccretionRate, double rad_efficiency);
double getLuminosity_radIneff_log(double BHaccretionRate, double maccr, double boundary);
double getLuminosity_radIneff_quot(double BHaccretionRate, double maccr, double boundary);
double getLuminosity_radIneff_lowMaccr(double BHmass, double maccr);
double getLuminosity_oscillations(double BHmass, double density, double density_crit, double cs, double gamma, double mu, double meanPhotEnergy, double rad_efficiency);
double getLuminosity_Sadowski(double BHmass, double BHspin, double BHaccretionRate);
double getBHaccretionRate_Hopkins(double BHmass, double rad_efficiency, double BHaccrete, double F, double dt);
double getBHaccretionMass_Hopkins(double BHmass, double rad_efficiency, double BHaccrete, double F, double dt);
double getLuminosity_Hopkins(double BHmass, double rad_efficiency, double BHaccrete, double F, double dt);


void strip_from_satellite(int halonr, int centralgal, int gal);
void disrupt_satellite_to_ICS(int centralgal, int gal);
void quasar_mode_wind(int gal, float BHaccrete);

double get_metallicity(double gas, double metals);
double get_virial_velocity(int halonr);
double get_virial_velocity_evolving(double Mvir, double Rvir);
double get_virial_radius(int halonr);
double get_virial_radius_evolving(int halonr, double Mvir);
double get_virial_mass(int halonr);
double get_disk_radius(int halonr, int p);
double get_disk_radius_evolving(int halonr, double Mvir, double Rvir);

void read_output_snaps(void);
void read_snap_list(void);
void read_cooling_functions(void);
double get_metaldependent_cooling_rate(double logTemp, double logZ);
double get_rate(int tab, double logTemp);

double time_to_present(double z);
double integrand_time_to_present(double a, void *param);

double metallicity_dependent_star_formation(int p);
double Z_dependent_SF(float lower_limit, float upper_limit, float Sigma_c0, float Xi, float gamma);
double integrand_Z_dependent_SF(double q, void *p);
