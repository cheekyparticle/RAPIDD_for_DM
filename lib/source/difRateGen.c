/*#########################################

  _____            _____ _____ _____  _____  
 |  __ \     /\   |  __ \_   _|  __ \|  __ \ 
 | |__) |   /  \  | |__) || | | |  | | |  | |
 |  _  /   / /\ \ |  ___/ | | | |  | | |  | |
 | | \ \  / ____ \| |    _| |_| |__| | |__| |
 |_|  \_\/_/    \_\_|   |_____|_____/|_____/ 
                                             
Based on : arXiv:1802.03174
difRateGen.c 
Calculates differential rates
early versions of this code were developed by M. Peiró
and E. Gerstmayr. 
##########################################*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>

#include "difRateGen.h"
#include "halo.h"
#include "phys_consts.h"
#include "coeffs_eft.h"
#include "effFormFact.h"
#include "atomic_resp.h"

#define SQR(X) ((X)*(X))
#define ABS(X) ((X) > 0 ? (X) : (-(X)))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

#define VERBOSE 0

/*########################################################################################
Combining results from the velocity distribution integrals and the differential cross section
The return value is given in events/(keV*kg*day)
########################################################################################*/

/*#####################################################################################
First we define the spin independent differential rate in the canonical way
#####################################################################################*/

double difrate_SI_standard(double rho, int A, int Z, double fp, double fn, double mchi, double Er){
	double mN = approx_mass_nucleus(A, Z);
	double muN = reduced_mass(mN, mchi);

	double velmin = c*sqrt(1. / (2.*mN*Er*1e-6))*mN*Er*1e-6 / muN;
	double rate = (2.*rho / (mchi*M_PI))*pow((Z*fp + (A - Z)*fn), 2.)*pow(fsi(A, Er), 2.)*halo(velmin, 0)*4.36e+5;
	return rate;
}  

/*#####################################################################################
Below are the functions for the EFT. 
#####################################################################################*/

double difrate_isotope_dEr(int A, int Z, double rhochi, void * input_difcros, int i_coeff, int j_coeff){
	/*Calculate the differential rate dR/dE_R in events/(keV*kg*day) given the NREFT coefficients
	specified by the indexes i_coeff and j_coeff.
	Parameters:
	- A: mass number
	- Z: atomic number
	- rhochi: dark matter density in GeV/cm^3
	- input_difcros: pointer to the differential cross section parameters
	- i_coeff: index of the first NREFT coefficient
	- j_coeff: index of the second NREFT coefficient
	*/
    struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;
	double Er = (val_difcros->Er);
 	double mchi = (val_difcros->mchi);
  	double jchi = (val_difcros->jchi);
	double v_h4 = higgs_vev*higgs_vev*higgs_vev*higgs_vev;
	char * nuclear_framework = (val_difcros->Nucleon);
	double delta = (val_difcros->delta);
	/*
	From natural units to events/(keV*kg*day):
	GeV/cm^3 1/(GeV^5) (km/s)^-1 = 1.698603e+14 (kg day keV)^-1
	GeV/cm^3 1/(GeV^5) km/s = 1.889947e+3 (kg day keV)^-1
	*/
	double conv_factor_v0 = 1.698603e+14;
	double conv_factor_v2 = 1.889947e+3;
	double mTarget = approx_mass_nucleus(A, Z);
	double muN = reduced_mass(mchi, mTarget);
	double Er_GeV = Er * 1.e-6;
	double delta_GeV = delta * 1.e-6;
	double vmin = c * fabs((mTarget * Er_GeV / muN) + delta_GeV) / sqrt(2. * mTarget * Er_GeV);
	double c_p_i, c_n_i, c_p_j, c_n_j;

	char p_name[2] = "";
	char n_name[2] = "";

	// Check the nuclear framework and set the names for basis accordingly (pn: proton/neutron, iso: isospin (+/-))
	if (strncmp(nuclear_framework, "pn", 2) == 0) {
		strcpy(p_name, "p");
		strcpy(n_name, "n");
		c_p_i = Cp(i_coeff);  // c_proton
		c_n_i = Cn(i_coeff);  // c_neutron
		c_p_j = Cp(j_coeff);  // c_proton
		c_n_j = Cn(j_coeff);  // c_neutron
	} else if (strncmp(nuclear_framework, "iso", 3) == 0) {
		strcpy(p_name, "+");
		strcpy(n_name, "-");
		c_p_i = (Cp(i_coeff)+Cn(i_coeff))/2.; // c_plus = (c_proton + c_neutron)/2
		c_n_i = (Cp(i_coeff)-Cn(i_coeff))/2.; // c_minus = (c_proton - c_neutron)/2
		c_p_j = (Cp(j_coeff)+Cn(j_coeff))/2.; // c_plus = (c_proton + c_neutron)/2
		c_n_j = (Cp(j_coeff)-Cn(j_coeff))/2.; // c_minus = (c_proton - c_neutron)/2
	} else {
		fprintf(stderr, "nuclear_framework must start with 'pn' or 'iso' to be valid: %s\n", nuclear_framework);
		exit(1);
	}

	double v0term = (
		c_p_i * c_p_j * FormFact_v0(nuclear_framework, A, Z, i_coeff, j_coeff, p_name, p_name, Er, mchi, jchi) +
		c_n_i * c_n_j * FormFact_v0(nuclear_framework, A, Z, i_coeff, j_coeff, n_name, n_name, Er, mchi, jchi) +
		c_p_i * c_n_j * FormFact_v0(nuclear_framework, A, Z, i_coeff, j_coeff, p_name, n_name, Er, mchi, jchi) +
		c_n_i * c_p_j * FormFact_v0(nuclear_framework, A, Z, i_coeff, j_coeff, n_name, p_name, Er, mchi, jchi))
		* halo(vmin, 0) * conv_factor_v0;
	
	double v2term = (
		c_p_i * c_p_j * FormFact_v2(nuclear_framework, A, Z, i_coeff, j_coeff, p_name, p_name, Er, mchi, jchi) +
		c_n_i * c_n_j * FormFact_v2(nuclear_framework, A, Z, i_coeff, j_coeff, n_name, n_name, Er, mchi, jchi) +
		c_p_i * c_n_j * FormFact_v2(nuclear_framework, A, Z, i_coeff, j_coeff, p_name, n_name, Er, mchi, jchi) +
		c_n_i * c_p_j * FormFact_v2(nuclear_framework, A, Z, i_coeff, j_coeff, n_name, p_name, Er, mchi, jchi))
		* halo(vmin, 2) * conv_factor_v2;

	double coeff_times_v_avg_formfact = v0term + v2term;
	double difrate = rhochi / (2. * M_PI * mchi * v_h4) * coeff_times_v_avg_formfact;
	
	return difrate;
}

double total_difrate_isotope_dEr(int A, int Z, double rhochi, void * input_difcros){

	struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;
	double rate = 0.0;
	int i;
	for ( i = 1; i < 16; i++)
        {
          if (Cp(i) != 0. || Cn(i) != 0.){
			rate += difrate_isotope_dEr(A, Z, rhochi, val_difcros, i, i);
          }
        }
	//INTEREFERENCE TERMS
	// Interferences between different nuclear responses not need both ij and ji terms
	// as they already include both contributions in their definition.
	// For the ones that comes from the same nuclear responses, we need to explicitly add moth ij and ji terms.
	// see eq. 38 and 89 of https://arxiv.org/abs/1308.6288
	if ((Cp(1) != 0. || Cn(1) != 0.) && (Cp(3) != 0. || Cn(3) != 0.)){
		rate += difrate_isotope_dEr(A, Z, rhochi, input_difcros, 1, 3);
	}
	if ((Cp(4) != 0. || Cn(4) != 0.) && (Cp(5) != 0. || Cn(5) != 0.)){
		rate += difrate_isotope_dEr(A, Z, rhochi, input_difcros, 4, 5);
	}
	if ((Cp(4) != 0. || Cn(4) != 0.) && (Cp(6) != 0. || Cn(6) != 0.)){
		rate += (difrate_isotope_dEr(A, Z, rhochi, input_difcros, 4, 6) + difrate_isotope_dEr(A, Z, rhochi, input_difcros, 6, 4));
	}
	if ((Cp(8) != 0. || Cn(8) != 0.) && (Cp(9) != 0. || Cn(9) != 0.)){
		rate += difrate_isotope_dEr(A, Z, rhochi, input_difcros, 8, 9);
	}
	if ((Cp(11) != 0. || Cn(11) != 0.) && (Cp(12) != 0. || Cn(12) != 0.)){
		rate += difrate_isotope_dEr(A, Z, rhochi, input_difcros, 11, 12);
	}
	if ((Cp(11) != 0. || Cn(11) != 0.) && (Cp(15) != 0. || Cn(15) != 0.)){
		rate += difrate_isotope_dEr(A, Z, rhochi, input_difcros, 11, 15);
	}
	if ((Cp(12) != 0. || Cn(12) != 0.) && (Cp(15) != 0. || Cn(15) != 0.)){
		rate += difrate_isotope_dEr(A, Z, rhochi, input_difcros, 12, 15);
	}

    return fabs(rate);
}

int nuclear_model_supports_target(char* NR_framework, char* target){
	/*Checks that the nuclear form factors of NR_framework are tabulated for every isotope
	of target used by difrate_dER (isotopes with zero abundance are skipped). Each tabulated
	isotope has the M-M response, so the check evaluates it at a small recoil energy.
	Without this check difrate_dER silently returns 0 for a missing isotope.
	Parameters:
	- NR_framework: nuclear model, e.g. "iso_GCN5082", "pn_Fitz"
	- target: "Xe", "Ge", "Ar", "F" or "CaWO4"
	Output: 1 if all the isotopes are available, 0 otherwise (also for an unknown target).*/
	const int *A_list;
	const double *abundance;
	int num_isos;
	if (strncmp(target, "Xe", 10) == 0) { A_list = atomic_numbers_Xe; abundance = isotopes_Xe; num_isos = 7; }
	else if (strncmp(target, "Ge", 10) == 0) { A_list = atomic_numbers_Ge; abundance = isotopes_Ge; num_isos = 5; }
	else if (strncmp(target, "Ar", 10) == 0) { A_list = atomic_numbers_Ar; abundance = isotopes_Ar; num_isos = 1; }
	else if (strncmp(target, "F", 10) == 0) { A_list = atomic_numbers_F; abundance = isotopes_F; num_isos = 1; }
	else if (strncmp(target, "CaWO4", 10) == 0) { A_list = atomic_numbers_CaWO4; abundance = isotopes_CaWO4; num_isos = 2; }
	else return 0;
	char * nucleon = (strncmp(NR_framework, "pn", 2) == 0) ? "p" : "+";
	for (int l = 0; l < num_isos; l++) {
		if (abundance[l] == 0.) continue;
		if (FF(NR_framework, A_list[l], "M", "M", nucleon, nucleon, 0.1) == 0.) return 0;
	}
	return 1;
}

double difrate_dER(double rhochi, void * input_difcros, double logenergy, char* model){
	struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;
	char * target = (val_difcros->target);
	double isotopes[10];
	double znumarr[10];
	int atomic_numbers[10]; 
	int num_isos, znum;
	if (strncmp(target, "Xe", 10) == 0){
		num_isos = 7;
		znum = 54;
		for (int l = 0; l < num_isos; l++){
			znumarr[l]=znum;
			isotopes[l]=isotopes_Xe[l];	
			atomic_numbers[l]=atomic_numbers_Xe[l];	
		}
	}
	if (strncmp(target, "Ge", 10) == 0){
		num_isos = 5;
		znum = 32; 
		for (int l = 0; l < num_isos; l++){
			znumarr[l]=znum;
			isotopes[l]=isotopes_Ge[l];	
			atomic_numbers[l]=atomic_numbers_Ge[l];	
		}
	}
	if (strncmp(target, "Ar", 10) == 0){
		num_isos = 1;
		znum = 18; 
		for (int l = 0; l < num_isos; l++){
			znumarr[l]=znum;
			isotopes[l]=isotopes_Ar[l];	
			atomic_numbers[l]=atomic_numbers_Ar[l];	
		}
	}
	if (strncmp(target, "F", 10) == 0){
		num_isos = 1;
		znum = 9; 
		for (int l = 0; l < num_isos; l++){
			znumarr[l]=znum;
			isotopes[l]=isotopes_F[l];	
			atomic_numbers[l]=atomic_numbers_F[l];	
		}
	}
	if (strncmp(target, "CaWO4", 10) == 0){
		num_isos = 2;
		for (int l = 0; l < num_isos; l++){
			znumarr[l]=Z_numbers_CaWO4[l];
			isotopes[l]=isotopes_CaWO4[l];	
			atomic_numbers[l]=atomic_numbers_CaWO4[l];
			printf("Z %i A %i number %f\n", Z_numbers_CaWO4[l],atomic_numbers[l],isotopes[l]);

		}
	}
	double counts=0.0;
	double energy = pow(10., logenergy);
	val_difcros->Er = energy;
	for (int l = 0; l< num_isos; l++){
		//printf("check %i, %lf , %i, %lf\n", l, isotopes[l], atomic_numbers[l],total_difrate_isotope_dEr(atomic_numbers[l],znum, rhochi, val_difcros));
		if (strncmp(model, "Light_Med", 10)==0)
		{
			double med_mass = give_med_mass();
			double prefact = 1./(2*approx_mass_nucleus(atomic_numbers[l],znumarr[l])*energy*1.e-6 + med_mass*med_mass);
			counts +=  prefact*prefact*isotopes[l]*total_difrate_isotope_dEr(atomic_numbers[l],znumarr[l], rhochi, val_difcros);
		}
		else if (strncmp(model, "ElecMom", 10)==0)
		{	
			double C11p = Cp(11);
			set_coeffs();
			set_any_Ncoeff(C11p, 11, "p");
			double med_mass = 0.0;
			double prefact = 1./(2*approx_mass_nucleus(atomic_numbers[l],znumarr[l])*energy*1.e-6 + med_mass*med_mass);
			counts +=  prefact*prefact*isotopes[l]*total_difrate_isotope_dEr(atomic_numbers[l],znumarr[l], rhochi, val_difcros);
		}
		else if (strncmp(model, "MagMom", 10)==0){
			double C1p, C4p, C4n, C5p, C6p, C6n;
			C1p = Cp(1);
			C4p = Cp(4);
			C4n = Cn(4);
			C5p = Cp(5);
			C6p = Cp(6);
			C6n = Cn(6);
			set_coeffs();
			set_any_Ncoeff(C1p, 1, "p");
			set_any_Ncoeff(C4p, 4, "p");
			set_any_Ncoeff(C4n, 4, "n");
			counts +=  isotopes[l]*total_difrate_isotope_dEr(atomic_numbers[l],znumarr[l], rhochi, val_difcros);
			
			set_coeffs();
			set_any_Ncoeff(C5p, 5, "p");
			set_any_Ncoeff(C6p, 6, "p");
			set_any_Ncoeff(C6n, 6, "n");
			double med_mass = 0.0;
			double prefact = 1./(2*approx_mass_nucleus(atomic_numbers[l],znumarr[l])*energy*1.e-6 + med_mass*med_mass);
			counts +=  prefact*prefact*isotopes[l]*total_difrate_isotope_dEr(atomic_numbers[l],znumarr[l], rhochi, val_difcros);
		}
		else {
			counts +=  isotopes[l]*total_difrate_isotope_dEr(atomic_numbers[l],znumarr[l], rhochi, val_difcros);
		}
	}
	return counts;
}

double difrate_dER_python(double rhochi, double mass, double logenergy, char* model, char* Target, char*ISO_switch, double delta){
  struct difcros_params struct_difcros_T= {Target, 1., mass, 0.5, 220, ISO_switch, 0.0, delta};
  
  double rate = difrate_dER(rhochi, &struct_difcros_T, logenergy, model);
  return rate;
}

double difrate_dER_python_2(double rhochi, double mass, double energy, char* model, char* Target, char*ISO_switch){
  struct difcros_params struct_difcros_T= {Target, 1., mass, 0.5, 220, ISO_switch, 0.0};
  
  double rate = difrate_dER_2(rhochi, &struct_difcros_T, energy, model);
  return rate;
}

double difrate_dER_2(double rhochi, void * input_difcros, double energy, char* model){

	struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;
	char * target = (val_difcros->target);
	double isotopes[10];
	double znumarr[10];
	int atomic_numbers[10]; 
	int num_isos, znum;
	if (strncmp(target, "Xe", 10) == 0){
		num_isos = 7;
		znum = 54;
		for (int l = 0; l < num_isos; l++){
			znumarr[l]=znum;
			isotopes[l]=isotopes_Xe[l];	
			atomic_numbers[l]=atomic_numbers_Xe[l];	
		}
	}
	if (strncmp(target, "Ge", 10) == 0){
		num_isos = 5;
		znum = 32; 
		for (int l = 0; l < num_isos; l++){
			znumarr[l]=znum;
			isotopes[l]=isotopes_Ge[l];	
			atomic_numbers[l]=atomic_numbers_Ge[l];	
		}
	}
	if (strncmp(target, "Ar", 10) == 0){
		num_isos = 1;
		znum = 18; 
		for (int l = 0; l < num_isos; l++){
			znumarr[l]=znum;
			isotopes[l]=isotopes_Ar[l];	
			atomic_numbers[l]=atomic_numbers_Ar[l];	
		}
	}
	if (strncmp(target, "F", 10) == 0){
		num_isos = 1;
		znum = 9; 
		for (int l = 0; l < num_isos; l++){
			znumarr[l]=znum;
			isotopes[l]=isotopes_F[l];	
			atomic_numbers[l]=atomic_numbers_F[l];	
		}
	}
	if (strncmp(target, "CaWO4", 10) == 0){
		num_isos = 2;
		for (int l = 0; l < num_isos; l++){
			znumarr[l]=Z_numbers_CaWO4[l];
			isotopes[l]=isotopes_CaWO4[l];	
			atomic_numbers[l]=atomic_numbers_CaWO4[l];
			printf("Z %i A %i number %f\n", Z_numbers_CaWO4[l],atomic_numbers[l],isotopes[l]);

		}
	}
	double counts=0.0;
	//double energy = pow(10., logenergy);
	val_difcros->Er = energy;
	for (int l = 0; l< num_isos; l++){
		//printf("check %i, %lf , %i, %lf\n", l, isotopes[l], atomic_numbers[l],total_difrate_isotope_dEr(atomic_numbers[l],znum, rhochi, val_difcros));
		if (strncmp(model, "Light_Med", 10)==0)
		{
			double med_mass = give_med_mass();
			double prefact = 1./(2*approx_mass_nucleus(atomic_numbers[l],znumarr[l])*energy*1.e-6 + med_mass*med_mass);
			counts +=  prefact*prefact*isotopes[l]*total_difrate_isotope_dEr(atomic_numbers[l],znumarr[l], rhochi, val_difcros);
		}
		else if (strncmp(model, "ElecMom", 10)==0)
		{	
			double C11p = Cp(11);
			set_coeffs();
			set_any_Ncoeff(C11p, 11, "p");
			double med_mass = 0.0;
			double prefact = 1./(2*approx_mass_nucleus(atomic_numbers[l],znumarr[l])*energy*1.e-6 + med_mass*med_mass);
			counts +=  prefact*prefact*isotopes[l]*total_difrate_isotope_dEr(atomic_numbers[l],znumarr[l], rhochi, val_difcros);
		}
		else if (strncmp(model, "MagMom", 10)==0){
			double C1p, C4p, C4n, C5p, C6p, C6n;
			C1p = Cp(1);
			C4p = Cp(4);
			C4n = Cn(4);
			C5p = Cp(5);
			C6p = Cp(6);
			C6n = Cn(6);
			set_coeffs();
			set_any_Ncoeff(C1p, 1, "p");
			set_any_Ncoeff(C4p, 4, "p");
			set_any_Ncoeff(C4n, 4, "n");
			counts +=  isotopes[l]*total_difrate_isotope_dEr(atomic_numbers[l],znumarr[l], rhochi, val_difcros);
			
			set_coeffs();
			set_any_Ncoeff(C5p, 5, "p");
			set_any_Ncoeff(C6p, 6, "p");
			set_any_Ncoeff(C6n, 6, "n");
			double med_mass = 0.0;
			double prefact = 1./(2*approx_mass_nucleus(atomic_numbers[l],znumarr[l])*energy*1.e-6 + med_mass*med_mass);
			counts +=  prefact*prefact*isotopes[l]*total_difrate_isotope_dEr(atomic_numbers[l],znumarr[l], rhochi, val_difcros);
		}
		else {
			counts +=  isotopes[l]*total_difrate_isotope_dEr(atomic_numbers[l],znumarr[l], rhochi, val_difcros);
		}
	}
	return counts;
}

/*########################################################################################
DM-ELECTRON SCATTERING
Ionization rate of Xenon from DM-electron scattering for the NREFT operators O_1, O_3 ... O_15
of arXiv:1912.08204 (contact and long-range couplings, atomic responses W_1 ... W_4) and
XENON10 / XENON1T (S2-only) detector response. The SI/SD rate and the detector response are
a port of the legacy MadDM Fortran (dm_response_direct_e.f and electron_recoil_signal.f).
Please cite: arXiv:1912.08204, arXiv:1509.01598
########################################################################################*/

#define N_ELECTR_MAX_XENON10 90   // max number of electrons for XENON10
#define S2_MAX_XENON10 203        // number of S2 values for XENON10
#define N_BINS_XENON10 7
#define S2_BIN_MAX_XENON1T 63     // number of S2 bins for XENON1T
#define LEN_E_RESP_XENON1T 41     // number of energies of the XENON1T S2 response
#define N_MASS_XENON1T 15         // number of DM masses of the XENON1T results table

static const int bins_Xenon10[N_BINS_XENON10 + 1] = {14, 41, 68, 95, 122, 149, 176, 203};
static const int n_obs_Xenon10[N_BINS_XENON10] = {126, 60, 12, 3, 2, 0, 2};

/* S2 values corresponding to the XENON1T S2 bin index (iS2 = 0 -> S2 = 90.7964 ...),
https://github.com/XENON1T/s2only_data_release/blob/master/s2_binning_info.csv */
static const double S2_val_Xenon1T[S2_BIN_MAX_XENON1T] = {
	90.7964, 92.4105, 94.0533, 95.7253, 97.427, 99.159, 100.9218, 102.7158, 104.5418,
	106.4003, 108.2918, 110.2169, 112.1762, 114.1704, 116.2, 118.2657, 120.3681,
	122.5079, 124.6857, 126.9022, 129.1582, 131.4543, 133.7911, 136.1695, 138.5902,
	141.054, 143.5615, 146.1136, 148.711, 151.3547, 154.0453, 156.7838, 159.571,
	162.4077, 165.2948, 168.2332, 171.2239, 174.2678, 177.3658, 180.5188, 183.7279,
	186.994, 190.3182, 193.7015, 197.145, 200.6496, 204.2166, 207.847, 211.5419,
	215.3025, 219.1299, 223.0254, 226.9901, 231.0254, 235.1323, 239.3123, 243.5665,
	247.8964, 252.3033, 256.7885, 261.3535, 265.9996, 270.7282};

static double ce_q(int i, double q, double m_e){
	/*DM-electron coupling c_i(q) = c_i^s + (q_ref/q)^2 c_i^l of arXiv:1912.08204, with
	q_ref = alpha m_e. The long-range part is added only if c_i^l is not zero, so that
	contact interactions can be evaluated at q = 0.
	Inputs: i (operator index), q (exchanged momentum, GeV), m_e (GeV).*/
	double ci = Ce(i);
	if (Ce_LR(i) != 0.) ci += SQR(alpha_em * m_e / q) * Ce_LR(i);
	return ci;
}

static void dm_responses_electron(double q, double dE, double mchi, double jchi, double m_e,
								  double *R, double *R1_v2){
	/*DM responses R_1 ... R_4 of arXiv:1912.08204 (Eqs. (DM response 1-4), real couplings,
	operators O_1, O_3 ... O_15) without unit factors. v_el^perp is evaluated at k = 0 (the k
	dependence is in the atomic responses W_2 ... W_4) and energy conservation gives (App. C)
	    (v_el^perp)^2 = v^2 + q^2/(4 mu^2) (mchi - m_e)/(mchi + m_e) - dE/mu,
	    v_el^perp . q = dE - q^2/(2 m_e),
	with velocities in units of c and mu the DM-electron reduced mass. R_1 is split as
	R_1 = R[0] + R1_v2 v^2: R[0] is integrated with eta_0 = int f/v and R1_v2 with
	int f v. R_2, R_3 and R_4 do not depend on v at fixed q and dE.
	Inputs: q (GeV), dE = E_B + E_e (GeV), mchi (GeV), jchi (DM spin), m_e (GeV).
	Outputs: R[0..3] (v independent parts of R_1 ... R_4), *R1_v2 (coefficient of v^2 in R_1).*/
	double c[OPERATOR_NUMBER];
	for (int i = 1; i < OPERATOR_NUMBER; i++) c[i] = ce_q(i, q, m_e);
	double jj = jchi * (jchi + 1.);
	double mu = mchi * m_e / (mchi + m_e);
	double q2 = SQR(q / m_e);                          // (q/m_e)^2
	double qv = (dE - q * q / (2. * m_e)) / m_e;       // (q/m_e) . v_el^perp
	double vperp2_const = q * q / (4. * mu * mu) * (mchi - m_e) / (mchi + m_e) - dE / mu; // (v_el^perp)^2 - v^2

	// coefficients of (v_el^perp)^2 and of ((q/m_e) . v_el^perp)^2 in R_1: they coincide with R_3 and R_4
	double A = SQR(c[3]) / 4. * q2 + SQR(c[7]) / 4.
			   + jj / 12. * (4. * SQR(c[8]) + 2. * SQR(c[12])
							 + (4. * SQR(c[5]) + SQR(c[13]) + SQR(c[14]) - 2. * c[12] * c[15]) * q2
							 + SQR(c[15]) * q2 * q2);
	double B = -SQR(c[3]) / 4.
			   + jj / 12. * (-4. * SQR(c[5]) - SQR(c[15]) * q2 + 2. * c[12] * c[15] + 2. * c[13] * c[14]);

	// R_1 without the v_el^perp terms (c_1 and c_4 first, as in the legacy SI/SD expression)
	double R1 = SQR(c[1]) + jj / 4. * SQR(c[4]);
	R1 += SQR(c[10]) / 4. * q2
		  + jj / 12. * (SQR(c[6]) * q2 * q2 + (2. * SQR(c[9]) + 4. * SQR(c[11]) + 2. * c[4] * c[6]) * q2);
	R[0] = R1 + A * vperp2_const + B * qv * qv;
	*R1_v2 = A;

	// R_2: the (q/m_e)^-2 terms are evaluated only if their couplings are not zero
	double R2_inv_q2 = SQR(c[7]) / 2. + jj / 6. * (4. * SQR(c[8]) + 2. * SQR(c[12]));
	double R2 = -jj / 6. * SQR(c[13] + c[14]);
	if (R2_inv_q2 != 0.) R2 -= R2_inv_q2 / q2;
	R[1] = qv * R2;
	R[2] = A;
	R[3] = B;
}

double dm_response_electron(int resp, double q, double dE, double mchi, double jchi, double m_e, double *R_v2){
	/*DM response R_resp of arXiv:1912.08204 built from the DM-electron NREFT coefficients
	Ce(i) (contact) and Ce_LR(i) (long-range), times the legacy MadDM normalisation
	gev_to_pb * pb_to_cm2. This is the only place where the unit factors are applied.
	For contact c_1, c_4 only, R_1 = (Ce(1)^2 + jchi(jchi+1)/4 Ce(4)^2) * gev_to_pb * pb_to_cm2
	at any q (legacy SI/SD expression, used by MadDM at q = 0).
	Excluded operators: O_2 and O_16 (as in arXiv:1912.08204) and O_17 ... O_20 (spin-1 DM,
	not in the DM responses of the paper).
	Inputs: resp (1..4), q (exchanged momentum, GeV), dE (deposited energy E_B + E_e, GeV),
	mchi (GeV), jchi (DM spin), m_e (GeV), R_v2 (can be NULL).
	Output: v independent part of R_resp; *R_v2 is the coefficient of v^2 (v in units of c),
	not zero only for resp = 1. 0 for a wrong resp.*/
	if (resp < 1 || resp > N_ATOMIC_RESP) {
		fprintf(stderr, "dm_response_electron: resp must be in [1, %d]\n", N_ATOMIC_RESP);
		if (R_v2 != NULL) *R_v2 = 0.;
		return 0.;
	}
	double R[N_ATOMIC_RESP], R1_v2;
	dm_responses_electron(q, dE, mchi, jchi, m_e, R, &R1_v2);
	if (R_v2 != NULL) *R_v2 = (resp == 1) ? R1_v2 * gev_to_pb * pb_to_cm2 : 0.;
	return R[resp - 1] * gev_to_pb * pb_to_cm2;
}

static void electron_grids(double m_e, double *k, double *E, double *q){
	/*Fills the log-spaced grids of the outgoing electron momentum k (GeV), its kinetic
	energy E (GeV), both of size GRIDSIZE_K+1, and of the exchanged momentum q (GeV),
	of size GRIDSIZE_Q+1. The last point of each grid lies past the tabulated range and
	is only used for the differentials (as in the legacy Fortran).*/
	int i;
	for (i = 0; i < GRIDSIZE_K + 1; i++) {
		k[i] = exp(log(K_MIN_GEV) + (log(K_MAX_GEV) - log(K_MIN_GEV)) * (double)i / (double)(GRIDSIZE_K - 1));
		E[i] = sqrt(m_e * m_e + k[i] * k[i]) - m_e;
	}
	for (i = 0; i < GRIDSIZE_Q + 1; i++) {
		q[i] = exp(log(Q_MIN_GEV) + (log(Q_MAX_GEV) - log(Q_MIN_GEV)) * (double)i / (double)(GRIDSIZE_Q - 1));
	}
}

#ifdef ER_VALIDATE_LEGACY_HALO
/* Debug only: analytic SHM eta averaged over the day bins with the legacy Earth velocity,
as in the legacy Fortran (get_dRdlogE + get_vEarth). Used to validate the port. */
static double legacy_vMP = 238., legacy_vesc = 544., legacy_vEarth_mod = 29.8;
static double legacy_vSun[3] = {11.1, 12.2, 7.3};
static double legacy_day_min = 0., legacy_day_max = 365.;
static int legacy_day_bins = 10;

void set_legacy_halo_params(double vMP, double vesc, double vEarth_mod, double vSun_r, double vSun_phi,
							double vSun_theta, double day_min, double day_max, int day_bins){
	/*Sets the parameters of the legacy (debug) halo, as in the maddm_card.*/
	legacy_vMP = vMP;
	legacy_vesc = vesc;
	legacy_vEarth_mod = vEarth_mod;
	legacy_vSun[0] = vSun_r;
	legacy_vSun[1] = vSun_phi;
	legacy_vSun[2] = vSun_theta;
	legacy_day_min = day_min;
	legacy_day_max = day_max;
	legacy_day_bins = day_bins;
}

static double legacy_vEarth(double days){
	/*Legacy get_vEarth: Earth speed in the galactic frame (km/s), days since March 22, 2018.*/
	double omega = 0.0172;
	double cosdt = cos(omega * days);
	double sindt = sin(omega * days);
	double v1 = legacy_vSun[0] + legacy_vEarth_mod * (0.9941 * cosdt - 0.0504 * sindt);
	double v2 = legacy_vMP + legacy_vSun[1] + legacy_vEarth_mod * (0.1088 * cosdt + 0.4946 * sindt);
	double v3 = legacy_vSun[2] + legacy_vEarth_mod * (0.0042 * cosdt - 0.8677 * sindt);
	return sqrt(v1 * v1 + v2 * v2 + v3 * v3);
}

static double legacy_eta(double vmin){
	/*Legacy eta (s/cm) at vmin (cm/s), averaged over the day bins.*/
	double v0 = legacy_vMP * 1.e5;
	double vesc = legacy_vesc * 1.e5;
	double kNorm = (v0 * v0 * v0) * M_PI * (sqrt(M_PI) * erf(vesc / v0) - 2. * (vesc / v0) * exp(-SQR(vesc / v0)));
	double dday = (legacy_day_max - legacy_day_min) / (double)legacy_day_bins;
	double eta_avg = 0.;
	for (int id = 0; id < legacy_day_bins; id++) {
		double daymid = legacy_day_min + ((double)id + 0.5) * dday;
		double vE = legacy_vEarth(daymid) * 1.e5;
		double eta;
		if (vmin <= (vesc - vE)) {
			eta = (v0 * v0 * M_PI) / (2. * vE * kNorm) *
				((-4.) * exp(-SQR(vesc / v0)) * vE + sqrt(M_PI) * v0 * (erf((vmin + vE) / v0) - erf((vmin - vE) / v0)));
		}
		else if ((vesc - vE) <= vmin && vmin <= (vesc + vE)) {
			eta = (v0 * v0 * M_PI) / (2. * vE * kNorm) *
				((-2.) * exp(-SQR(vesc / v0)) * (vesc - vmin + vE) + sqrt(M_PI) * v0 * (erf(vesc / v0) - erf((vmin - vE) / v0)));
		}
		else {
			eta = 0.;
		}
		eta_avg += eta * dday / (legacy_day_max - legacy_day_min);
	}
	return eta_avg;
}
#endif

/* Sign of the tabulated atomic response W_2: -1 (default) corrects the sign of the vector
form factor of arXiv:1912.08204 found in arXiv:2405.04855 (Liang et al.), +1 uses W_2 as in
arXiv:1912.08204 (Catena et al.). R_2 is the same in the two works. */
static double electron_w2_sign = -1.;

void set_electron_w2_sign(int sign){
	/*Sets the sign of the atomic response W_2 in the DM-electron rate.
	Inputs: sign (-1: arXiv:2405.04855, default; +1: arXiv:1912.08204).*/
	if (sign != 1 && sign != -1) {
		fprintf(stderr, "set_electron_w2_sign: sign must be +1 or -1\n");
		return;
	}
	electron_w2_sign = (double)sign;
}

int electron_w2_sign_value(){
	/*Returns the sign of W_2 used in the DM-electron rate (+1 or -1).*/
	return (int)electron_w2_sign;
}

static int coeff_e_set(int i){
	/*Returns 1 if the contact or the long-range DM-electron coefficient i is not zero.*/
	return Ce(i) != 0. || Ce_LR(i) != 0.;
}

static void electron_responses_needed(int resp_mask, int *needed, int *needs_v2){
	/*Selects the terms R_j W_j to compute: needed[j-1] = 1 if bit j-1 of resp_mask is set
	and R_j is not identically zero for the coefficients set (R_1 W_1 is always computed
	when selected). *needs_v2 = 1 if the v^2 part of R_1 (integrated with halo(vmin, 2))
	is needed.*/
	int vperp = coeff_e_set(3) || coeff_e_set(5) || coeff_e_set(7) || coeff_e_set(8) || coeff_e_set(12)
				|| coeff_e_set(13) || coeff_e_set(14) || coeff_e_set(15);
	needed[0] = (resp_mask & 1) != 0;
	needed[1] = (resp_mask & 2) != 0 && (coeff_e_set(7) || coeff_e_set(8) || coeff_e_set(12)
										 || coeff_e_set(13) || coeff_e_set(14));
	needed[2] = (resp_mask & 4) != 0 && vperp;
	needed[3] = (resp_mask & 8) != 0 && (coeff_e_set(3) || coeff_e_set(5) || coeff_e_set(15)
										 || (coeff_e_set(12) && coeff_e_set(15))
										 || (coeff_e_set(13) && coeff_e_set(14)));
	*needs_v2 = needed[0] && vperp;
}

static double difrate_electron_dEr_w(int shell, int ik, const double *E, const double *q, double rhochi,
									 double mchi, double jchi, double m_e, const int *needed,
									 gsl_interp_accel *ga, gsl_spline *gs, gsl_interp_accel *ga2, gsl_spline *gs2){
	/*Kernel of difrate_electron_dEr with the grids and the halo splines already built:
	dR/dlnE = n_chi/(128 pi mchi^2 m_e^2) int dq q sum_j [R_j eta_0 + delta_j1 R_1^(v2) eta_2] W_j
	(arXiv:1912.08204), with eta_0 = halo(vmin, 0) (spline ga, gs) and eta_2 = halo(vmin, 2)
	(spline ga2, gs2, can be NULL if the v^2 part of R_1 is not needed). Only the terms with
	needed[j-1] = 1 are summed. Rectangular sum over the q grid, as in the legacy Fortran.
	Output: dR/dE in events/(kg day GeV).*/
	double nDM = rhochi / mchi; // cm^-3
	double const_integral = nDM / (128. * M_PI * SQR(mchi * GeV_to_kg) * SQR(m_e * GeV_to_kg));
	double Eb = shell_binding_energy(shell);
	double conv = GeV_to_J * SQR(m_to_cm) / c_light_cms;
	double dRdlogE = 0.;
	int iq, j;
	for (iq = 0; iq < GRIDSIZE_Q; iq++) {
		double dq = q[iq + 1] - q[iq];
		double vmin = ((E[ik] + Eb) / q[iq] + q[iq] / (2. * mchi)) * c_light_cms; // cm/s
#ifdef ER_VALIDATE_LEGACY_HALO
		(void)ga;
		(void)gs;
		(void)ga2;
		double eta = legacy_eta(vmin);
		double eta2 = 0.;
#else
		double eta = halo_w(vmin * 1.e-5, ga, gs) * 1.e-5; // (km/s)^-1 -> s/cm
		// int f v d^3v in km/s -> int f/v (v/c)^2 d^3v in s/cm
		double eta2 = (gs2 != NULL) ? halo_w(vmin * 1.e-5, ga2, gs2) / SQR(c) * 1.e-5 : 0.;
#endif
		double R[N_ATOMIC_RESP], R1_v2;
		dm_responses_electron(q[iq], E[ik] + Eb, mchi, jchi, m_e, R, &R1_v2);
		double pref = const_integral * (dq * conv) * (q[iq] * conv);
		for (j = 0; j < N_ATOMIC_RESP; j++) {
			if (!needed[j]) continue;
			double W = atomic_response(j + 1, shell, ik, iq);
			if (j == 1) W *= electron_w2_sign;
			double ioniz_amplitude = W * (R[j] * gev_to_pb * pb_to_cm2) * SQR(GeV_to_kg);
			dRdlogE += pref * eta * ioniz_amplitude / (m_Xe_avg_au * atomic_mass_unit_kg);
			if (j == 0 && gs2 != NULL) {
				double ioniz_amplitude_v2 = W * (R1_v2 * gev_to_pb * pb_to_cm2) * SQR(GeV_to_kg);
				dRdlogE += pref * eta2 * ioniz_amplitude_v2 / (m_Xe_avg_au * atomic_mass_unit_kg);
			}
		}
	}
	return dRdlogE * day_to_sec / E[ik];
}

static gsl_spline * electron_halo_spline(int i, gsl_interp_accel **ga){
	/*Builds once the cubic spline of the halo integral halo(vmin, i) read by read_halo.*/
	double v[1000], f[1000];
	int n = halo_table(i, v, f);
	*ga = gsl_interp_accel_alloc();
	gsl_spline *gs = gsl_spline_alloc(gsl_interp_cspline, n);
	gsl_spline_init(gs, v, f, n);
	return gs;
}

struct electron_rate_setup { int needed[N_ATOMIC_RESP]; gsl_interp_accel *ga, *ga2; gsl_spline *gs, *gs2; };

static int electron_rate_setup_init(const char *caller, int resp_mask, struct electron_rate_setup *st){
	/*Checks the inputs needed by the terms R_j W_j selected by resp_mask and the coefficients
	set, and builds the halo splines.
	Output: 0 on success, 1 if a needed atomic response W_j has not been read, 2 if the halo
	table has no halo(vmin, 2) column while the v^2 part of R_1 is needed, 3 for a wrong resp_mask.*/
	int needs_v2, j;
	st->ga = st->ga2 = NULL;
	st->gs = st->gs2 = NULL;
	if (resp_mask <= 0 || resp_mask > 15) {
		fprintf(stderr, "%s: resp_mask must be in [1, 15] (bit j-1 selects R_j W_j)\n", caller);
		return 3;
	}
	electron_responses_needed(resp_mask, st->needed, &needs_v2);
	for (j = 0; j < N_ATOMIC_RESP; j++) {
		if (st->needed[j] && !atomic_responses_loaded(j + 1)) {
			fprintf(stderr, "%s: the atomic response W_%d is needed by the coefficients set: "
					"call read_atomic_responses first\n", caller, j + 1);
			return 1;
		}
	}
#ifdef ER_VALIDATE_LEGACY_HALO
	if (needs_v2) {
		fprintf(stderr, "%s: the legacy (debug) halo has only eta_0, the v^2 part of R_1 cannot be computed\n", caller);
		return 2;
	}
#else
	if (needs_v2 && halo_table_power() < 2) {
		fprintf(stderr, "%s: the operators set need halo(vmin, 2): read a halo table with i >= 2\n", caller);
		return 2;
	}
#endif
	st->gs = electron_halo_spline(0, &st->ga);
	if (needs_v2) st->gs2 = electron_halo_spline(2, &st->ga2);
	return 0;
}

static void electron_rate_setup_free(struct electron_rate_setup *st){
	/*Frees the halo splines built by electron_rate_setup_init.*/
	if (st->gs != NULL) gsl_spline_free(st->gs);
	if (st->ga != NULL) gsl_interp_accel_free(st->ga);
	if (st->gs2 != NULL) gsl_spline_free(st->gs2);
	if (st->ga2 != NULL) gsl_interp_accel_free(st->ga2);
}

double difrate_electron_dEr(int shell, int ik, double rhochi, double mchi, double jchi, double m_e){
	/*Differential ionization rate dR/dE for one Xenon shell at the outgoing electron
	energy E(ik), using the coefficients Ce(i), Ce_LR(i), the atomic responses W_j (read
	with read_atomic_responses) and the halo read with read_halo (eta_0 = halo(vmin, 0),
	and halo(vmin, 2) for the operators depending on v_el^perp).
	Parameters:
	- shell: Xenon shell index (0-based, see xe_shells in atomic_resp.c)
	- ik: index of the k' grid (0..GRIDSIZE_K-1)
	- rhochi: dark matter density in GeV/cm^3
	- mchi: dark matter mass in GeV
	- jchi: dark matter spin
	- m_e: electron mass in GeV
	Output: dR/dE in events/(kg day GeV), 0 if the inputs are missing.*/
	double k[GRIDSIZE_K + 1], E[GRIDSIZE_K + 1], q[GRIDSIZE_Q + 1];
	struct electron_rate_setup st;
	if (electron_rate_setup_init("difrate_electron_dEr", 15, &st) != 0) {
		electron_rate_setup_free(&st);
		return 0.;
	}
	electron_grids(m_e, k, E, q);
	double rate = difrate_electron_dEr_w(shell, ik, E, q, rhochi, mchi, jchi, m_e, st.needed, st.ga, st.gs, st.ga2, st.gs2);
	electron_rate_setup_free(&st);
	return rate;
}

int difrate_electron_dEr_resp(double rhochi, double mchi, double jchi, double m_e, int resp_mask,
							  double *k_out, double *E_out, double *rate_shell_out){
	/*Fills in a single call the k' and E grids and dR/dE of all the Xenon shells, keeping
	only the terms R_j W_j selected by resp_mask (bit j-1 set: R_j W_j included; 15: all, 1: only
	R_1 W_1, 2: only R_2 W_2, ...). The single contributions can be negative.
	Parameters:
	- rhochi, mchi, jchi, m_e: see difrate_electron_dEr
	- resp_mask: selection of the terms R_j W_j (1..15)
	- k_out: outgoing electron momenta in GeV [GRIDSIZE_K+1]
	- E_out: outgoing electron energies in GeV [GRIDSIZE_K+1]
	- rate_shell_out: dR/dE in events/(kg day GeV) [N_SHELL_XE*GRIDSIZE_K], shell-major
	Output: 0 on success, 1 if a needed atomic response W_j has not been read, 2 if the halo
	table has no halo(vmin, 2) column while needed, 3 for a wrong resp_mask.*/
	struct electron_rate_setup st;
	int err = electron_rate_setup_init("difrate_electron_dEr_resp", resp_mask, &st);
	if (err != 0) {
		electron_rate_setup_free(&st);
		return err;
	}
	double q[GRIDSIZE_Q + 1];
	electron_grids(m_e, k_out, E_out, q);
	for (int s = 0; s < N_SHELL_XE; s++) {
		for (int ik = 0; ik < GRIDSIZE_K; ik++) {
			rate_shell_out[s * GRIDSIZE_K + ik] = difrate_electron_dEr_w(s, ik, E_out, q, rhochi, mchi, jchi, m_e,
																		 st.needed, st.ga, st.gs, st.ga2, st.gs2);
		}
	}
	electron_rate_setup_free(&st);
	return 0;
}

int difrate_electron_dEr_all(double rhochi, double mchi, double jchi, double m_e,
							 double *k_out, double *E_out, double *rate_shell_out){
	/*Fills in a single call the k' and E grids and dR/dE of all the Xenon shells (all the
	terms R_j W_j, see difrate_electron_dEr_resp).
	Parameters:
	- rhochi, mchi, jchi, m_e: see difrate_electron_dEr
	- k_out: outgoing electron momenta in GeV [GRIDSIZE_K+1]
	- E_out: outgoing electron energies in GeV [GRIDSIZE_K+1]
	- rate_shell_out: dR/dE in events/(kg day GeV) [N_SHELL_XE*GRIDSIZE_K], shell-major
	Output: 0 on success, 1 if a needed atomic response W_j has not been read (W_1 is always
	needed), 2 if the halo table has no halo(vmin, 2) column while needed.*/
	return difrate_electron_dEr_resp(rhochi, mchi, jchi, m_e, 15, k_out, E_out, rate_shell_out);
}

static double binomial_prob(int k, int n, double p){
	/*Binomial probability of k successes in n trials with success probability p.*/
	if (k > n) return 0.;
	double num = 1., den = 1.;
	for (int i = 1; i <= k; i++) {
		num *= (double)(n - i + 1);
		den *= (double)i;
	}
	return num / den * pow(p, k) * pow(1. - p, n - k);
}

static double gauss_prob(int val, int mean, double sigma){
	/*Gaussian probability density at val with given mean and sigma.*/
	return 1. / (sigma * sqrt(2. * M_PI)) * exp(-0.5 * SQR((double)(val - mean) / sigma));
}

static double lin_interp(const double *x, const double *y, int n, double xi){
	/*Linear interpolation of y(x) at xi (x increasing), 0 outside the range of x.*/
	double xmin = x[0], xmax = x[0];
	for (int j = 1; j < n; j++) {
		xmin = MIN(xmin, x[j]);
		xmax = MAX(xmax, x[j]);
	}
	if (xi < xmin || xi > xmax) return 0.;
	for (int j = 0; j < n - 1; j++) {
		if (xi >= x[j] && xi <= x[j + 1]) {
			double slope = (y[j + 1] - y[j]) / (x[j + 1] - x[j]);
			double intercept = y[j] - slope * x[j];
			return slope * xi + intercept;
		}
	}
	return 0.;
}

static int read_column(const char *path, double *out, int n){
	/*Reads the first n numbers (one per line) of the file at path. Returns 0 on success.*/
	FILE *f = fopen(path, "r");
	if (f == NULL) {
		fprintf(stderr, "cannot open %s\n", path);
		return 1;
	}
	for (int i = 0; i < n; i++) {
		if (fscanf(f, "%lf", &out[i]) != 1) {
			fprintf(stderr, "%s is incomplete\n", path);
			fclose(f);
			return 1;
		}
	}
	fclose(f);
	return 0;
}

int xenon10_bins(int *bins_out){
	/*Fills bins_out [N_BINS_XENON10+1] with the XENON10 S2 bin edges. Returns the number of edges.*/
	for (int i = 0; i < N_BINS_XENON10 + 1; i++) bins_out[i] = bins_Xenon10[i];
	return N_BINS_XENON10 + 1;
}

int difrate_dS2_Xenon10(const double *E, const double *rate_shell, const char *exp_dir,
						double *S2_out, double *dRdS2_out, double *dRdS2_shell_out, double *signal_out, int *n_obs_out){
	/*XENON10 S2 spectrum from the per-shell dR/dE (arXiv:1703.00910): binomial number of
	electrons from n_1 = floor(E/W) primary plus n_2 secondary quanta, gaussian S2 per
	electron, trigger efficiency (exp_dir/Xenon10/Xenon10_TriggerEfficiency.txt), flat
	efficiency and exposure. As in the legacy Fortran, S2 = 1 is left without efficiency
	and exposure, and the shells with zero rate are skipped.
	Parameters:
	- E: outgoing electron energies in GeV [GRIDSIZE_K+1]
	- rate_shell: dR/dE in events/(kg day GeV) [N_SHELL_XE*GRIDSIZE_K]
	- exp_dir: path of RAPIDD_for_DM/lib/experiments
	- S2_out: S2 values [S2_MAX_XENON10]
	- dRdS2_out: dN/dS2 in events summed over the shells [S2_MAX_XENON10]
	- dRdS2_shell_out: per-shell dN/dS2 [N_SHELL_XE*S2_MAX_XENON10] (can be NULL)
	- signal_out: events in the 7 S2 bins; n_obs_out: observed events in the 7 S2 bins
	Output: number of S2 values, or -1 on error.*/
	const double W = 13.8e-9;               // mean energy to produce a quantum (GeV)
	const double f_e = 1. / 1.2;            // fraction of quanta observed as electrons
	const int g2 = 27;                      // secondary scintillation gain factor
	const double sigma_S2 = 6.7;
	const double exposure = 15.;            // kg day
	const double flat_efficiency = 0.92;
	double trigger_eff[S2_MAX_XENON10];
	char path[1024];
	snprintf(path, sizeof(path), "%s/Xenon10/Xenon10_TriggerEfficiency.txt", exp_dir);
	if (read_column(path, trigger_eff, S2_MAX_XENON10) != 0) return -1;

	int S2, ne, iE, s, i;
	for (S2 = 1; S2 <= S2_MAX_XENON10; S2++) {
		S2_out[S2 - 1] = S2;
		dRdS2_out[S2 - 1] = 0.;
	}
	if (dRdS2_shell_out != NULL) {
		for (i = 0; i < N_SHELL_XE * S2_MAX_XENON10; i++) dRdS2_shell_out[i] = 0.;
	}

	for (s = 0; s < N_SHELL_XE; s++) {
		const double *rate = &rate_shell[s * GRIDSIZE_K];
		double tot_rate = 0.;
		for (iE = 0; iE < GRIDSIZE_K; iE++) tot_rate += rate[iE] * (E[iE + 1] - E[iE]);
		if (tot_rate == 0.) continue;

		int n_tot[GRIDSIZE_K];
		for (iE = 0; iE < GRIDSIZE_K; iE++) n_tot[iE] = (int)floor(E[iE] / W) + shell_n_secondary(s);

		double dRdne[N_ELECTR_MAX_XENON10 + 1];
		for (ne = 1; ne <= N_ELECTR_MAX_XENON10; ne++) {
			dRdne[ne] = 0.;
			for (iE = 0; iE < GRIDSIZE_K; iE++) {
				dRdne[ne] += (E[iE + 1] - E[iE]) * rate[iE] * binomial_prob(ne - 1, n_tot[iE], f_e);
			}
		}

		for (S2 = 1; S2 <= S2_MAX_XENON10; S2++) {
			double dRdS2 = 0.;
			for (ne = 1; ne <= N_ELECTR_MAX_XENON10; ne++) {
				dRdS2 += dRdne[ne] * gauss_prob(S2, ne * g2, sqrt((double)ne) * sigma_S2);
			}
			if (S2 >= 2) dRdS2 *= exposure * trigger_eff[S2 - 2] * flat_efficiency;
			if (dRdS2_shell_out != NULL) dRdS2_shell_out[s * S2_MAX_XENON10 + S2 - 1] = dRdS2;
			dRdS2_out[S2 - 1] += dRdS2;
		}
	}

	for (i = 0; i < N_BINS_XENON10; i++) {
		signal_out[i] = 0.;
		n_obs_out[i] = n_obs_Xenon10[i];
	}
	for (S2 = 1; S2 <= S2_MAX_XENON10; S2++) {
		for (i = 0; i < N_BINS_XENON10; i++) {
			if (bins_Xenon10[i] <= S2 && S2 < bins_Xenon10[i + 1]) {
				signal_out[i] += dRdS2_out[S2 - 1];
				break;
			}
		}
	}
	return S2_MAX_XENON10;
}

int difrate_dS2_Xenon1T(const double *E, const double *rate_shell, double mchi, const char *exp_dir,
						double *S2_out, double *dRdS2_out, double *signal_out, double *n_bkg_out, int *n_obs_out){
	/*XENON1T S2-only spectrum from the total dR/dE, using the S2 response
	(exp_dir/Xenon1T_S2only/s2_response_er.csv, which includes efficiencies and cuts,
	https://github.com/XENON1T/s2only_data_release) and the raw exposure. The signal is
	summed over the mass dependent ROI of exp_dir/Xenon1T_S2only/5d_results_dmelectron.dat.
	Parameters:
	- E: outgoing electron energies in GeV [GRIDSIZE_K+1]
	- rate_shell: dR/dE in events/(kg day GeV) [N_SHELL_XE*GRIDSIZE_K]
	- mchi: dark matter mass in GeV
	- exp_dir: path of RAPIDD_for_DM/lib/experiments
	- S2_out: S2 values of the bins [S2_BIN_MAX_XENON1T]; dRdS2_out: events per bin
	- signal_out, n_bkg_out, n_obs_out: signal, background and observed events in the ROI
	  (-1 if mchi is outside the mass range of the results table)
	Output: number of S2 bins, or -1 on error.*/
	const double exposure = 356770.; // kg day (0.97678 tonne years, XENON1T data release)
	double E_resp[LEN_E_RESP_XENON1T], E_lo[LEN_E_RESP_XENON1T], E_hi[LEN_E_RESP_XENON1T];
	static double s2_resp[LEN_E_RESP_XENON1T][S2_BIN_MAX_XENON1T];
	char path[1024];
	int iE, iS2, ik;

	snprintf(path, sizeof(path), "%s/Xenon1T_S2only/s2_response_er.csv", exp_dir);
	FILE *f = fopen(path, "r");
	if (f == NULL) {
		fprintf(stderr, "cannot open %s\n", path);
		return -1;
	}
	char *line = NULL;
	size_t cap = 0;
	if (getline(&line, &cap, f) < 0) { // header
		fclose(f);
		free(line);
		return -1;
	}
	for (iE = 0; iE < LEN_E_RESP_XENON1T; iE++) {
		if (getline(&line, &cap, f) < 0) {
			fprintf(stderr, "%s is incomplete\n", path);
			fclose(f);
			free(line);
			return -1;
		}
		char *p = line, *end;
		double vals[3 + S2_BIN_MAX_XENON1T];
		for (int j = 0; j < 3 + S2_BIN_MAX_XENON1T; j++) {
			vals[j] = strtod(p, &end);
			p = (*end == ',') ? end + 1 : end;
		}
		E_resp[iE] = vals[0] / 1.e6; // keV -> GeV
		E_lo[iE] = vals[1] / 1.e6;
		E_hi[iE] = vals[2] / 1.e6;
		for (iS2 = 0; iS2 < S2_BIN_MAX_XENON1T; iS2++) s2_resp[iE][iS2] = vals[3 + iS2];
	}
	fclose(f);
	free(line);

	// total dR/dE; if the last energies differ by less than 5%, consider them as the same
	double x[GRIDSIZE_K], rate[GRIDSIZE_K];
	for (ik = 0; ik < GRIDSIZE_K; ik++) {
		x[ik] = E[ik];
		rate[ik] = 0.;
		for (int s = 0; s < N_SHELL_XE; s++) rate[ik] += rate_shell[s * GRIDSIZE_K + ik];
	}
	if ((E_resp[LEN_E_RESP_XENON1T - 1] - x[GRIDSIZE_K - 1]) / x[GRIDSIZE_K - 1] <= 0.05) {
		x[GRIDSIZE_K - 1] = E_resp[LEN_E_RESP_XENON1T - 1];
	}
	double rate_interp[LEN_E_RESP_XENON1T];
	for (iE = 0; iE < LEN_E_RESP_XENON1T; iE++) rate_interp[iE] = lin_interp(x, rate, GRIDSIZE_K, E_resp[iE]);

	for (iS2 = 0; iS2 < S2_BIN_MAX_XENON1T; iS2++) {
		S2_out[iS2] = S2_val_Xenon1T[iS2];
		dRdS2_out[iS2] = 0.;
		for (iE = 0; iE < LEN_E_RESP_XENON1T; iE++) {
			dRdS2_out[iS2] += rate_interp[iE] * (E_hi[iE] - E_lo[iE]) * s2_resp[iE][iS2] * exposure;
		}
	}

	// mass dependent ROI, background and observed events
	double m_tab[N_MASS_XENON1T], roi_min_tab[N_MASS_XENON1T], roi_max_tab[N_MASS_XENON1T];
	double n_obs_tab[N_MASS_XENON1T], n_bkg_tab[N_MASS_XENON1T];
	snprintf(path, sizeof(path), "%s/Xenon1T_S2only/5d_results_dmelectron.dat", exp_dir);
	f = fopen(path, "r");
	if (f == NULL) {
		fprintf(stderr, "cannot open %s\n", path);
		return -1;
	}
	char buf[512];
	int n = 0;
	while (n < N_MASS_XENON1T && fgets(buf, sizeof(buf), f) != NULL) {
		if (buf[0] == '#') continue;
		if (sscanf(buf, "%lf %lf %lf %lf %lf", &m_tab[n], &roi_min_tab[n], &roi_max_tab[n], &n_obs_tab[n], &n_bkg_tab[n]) == 5) n++;
	}
	fclose(f);
	if (n != N_MASS_XENON1T) {
		fprintf(stderr, "%s is incomplete\n", path);
		return -1;
	}

	if (mchi < m_tab[0] || mchi > m_tab[N_MASS_XENON1T - 1]) {
		printf("DM mass outside XENON1T bounds: (%.2f : %.1f) GeV\n", m_tab[0], m_tab[N_MASS_XENON1T - 1]);
		*signal_out = -1.;
		*n_bkg_out = -1.;
		*n_obs_out = -1;
		return S2_BIN_MAX_XENON1T;
	}
	int roi_min = (int)round(lin_interp(m_tab, roi_min_tab, N_MASS_XENON1T, mchi));
	int roi_max = (int)round(lin_interp(m_tab, roi_max_tab, N_MASS_XENON1T, mchi));
	*n_obs_out = (int)round(lin_interp(m_tab, n_obs_tab, N_MASS_XENON1T, mchi));
	*n_bkg_out = lin_interp(m_tab, n_bkg_tab, N_MASS_XENON1T, mchi);
	*signal_out = 0.;
	for (iS2 = roi_min; iS2 <= roi_max; iS2++) *signal_out += dRdS2_out[iS2 - 1]; // ROI is 1-based
	return S2_BIN_MAX_XENON1T;
}

int difrate_dS2(const char *experiment, const double *E, const double *rate_shell, double mchi, const char *exp_dir,
				double *S2_out, double *dRdS2_out, double *dRdS2_shell_out, double *signal_out, double *n_bkg_out, int *n_obs_out){
	/*Detector response in S2 for the DM-electron scattering.
	Parameters:
	- experiment: "Xenon10" or "Xenon1T"
	- E, rate_shell: grids and per-shell dR/dE from difrate_electron_dEr_all
	- mchi: dark matter mass in GeV (used by Xenon1T)
	- exp_dir: path of RAPIDD_for_DM/lib/experiments
	- S2_out, dRdS2_out: S2 values and dN/dS2 (203 values for Xenon10, 63 for Xenon1T)
	- dRdS2_shell_out: per-shell dN/dS2 (Xenon10 only, can be NULL)
	- signal_out, n_obs_out: 7 bins for Xenon10, 1 value (ROI) for Xenon1T
	- n_bkg_out: background in the ROI (Xenon1T only)
	Output: number of S2 values, or -1 on error.*/
	if (strncmp(experiment, "Xenon10", 8) == 0) {
		*n_bkg_out = 0.;
		return difrate_dS2_Xenon10(E, rate_shell, exp_dir, S2_out, dRdS2_out, dRdS2_shell_out, signal_out, n_obs_out);
	}
	if (strncmp(experiment, "Xenon1T", 8) == 0) {
		return difrate_dS2_Xenon1T(E, rate_shell, mchi, exp_dir, S2_out, dRdS2_out, signal_out, n_bkg_out, n_obs_out);
	}
	fprintf(stderr, "difrate_dS2: experiment must be Xenon10 or Xenon1T: %s\n", experiment);
	return -1;
}
