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
	specified by the indexes i_coeff and j_coeff.*/

    struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;
	double Er = (val_difcros->Er);
 	double mchi = (val_difcros->mchi);
  	double jchi = (val_difcros->jchi);
	double v_h4 = higgs_vev*higgs_vev*higgs_vev*higgs_vev;
	char * nuclear_framework = (val_difcros->Nucleon);
	/*
	From natural units to events/(keV*kg*day):
	GeV/cm^3 1/(GeV^5) (km/s)^-1 = 1.698603e+14 (kg day keV)^-1
	GeV/cm^3 1/(GeV^5) km/s = 1.889947e+3 (kg day keV)^-1
	*/
	double conv_factor_v0 = 1.698603e+14;
	double conv_factor_v2 = 1.889947e+3;
	double mTarget = approx_mass_nucleus(A, Z);
	double muN = reduced_mass(mchi, mTarget);
	double vmin = c*sqrt((mTarget*Er*1.e-6) / 2.) / muN;
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
	for ( i = 1; i < 16; i++) // FIXME what about Cp(0) et al?
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
	if ((Cp(1) != 0. || Cn(1) != 0.) & (Cp(3) != 0. || Cn(3) != 0.)){
		rate += difrate_isotope_dEr(A, Z, rhochi, input_difcros, 1, 3);
	}
	if ((Cp(4) != 0. || Cn(4) != 0.) & (Cp(5) != 0. || Cn(5) != 0.)){
		rate += difrate_isotope_dEr(A, Z, rhochi, input_difcros, 4, 5);
	}
	if ((Cp(4) != 0. || Cn(4) != 0.) & (Cp(6) != 0. || Cn(6) != 0.)){
		rate += (difrate_isotope_dEr(A, Z, rhochi, input_difcros, 4, 6) + difrate_isotope_dEr(A, Z, rhochi, input_difcros, 6, 4));
	}
	if ((Cp(8) != 0. || Cn(8) != 0.) & (Cp(9) != 0. || Cn(9) != 0.)){
		rate += difrate_isotope_dEr(A, Z, rhochi, input_difcros, 8, 9);
	}
	if ((Cp(11) != 0. || Cn(11) != 0.) & (Cp(12) != 0. || Cn(12) != 0.)){
		rate += difrate_isotope_dEr(A, Z, rhochi, input_difcros, 11, 12);
	}
	if ((Cp(11) != 0. || Cn(11) != 0.) & (Cp(15) != 0. || Cn(15) != 0.)){
		rate += difrate_isotope_dEr(A, Z, rhochi, input_difcros, 11, 15);
	}
	if ((Cp(12) != 0. || Cn(12) != 0.) & (Cp(15) != 0. || Cn(15) != 0.)){
		rate += difrate_isotope_dEr(A, Z, rhochi, input_difcros, 12, 15);
	}

    return fabs(rate);
}

double difrate_dER(double rhochi, void * input_difcros, double logenergy, char* model){
	struct difcros_params * val_difcros = (struct difcros_params *)input_difcros;
	char * target = (val_difcros->target);
	double isotopes[10];
	double znumarr[10];
	int atomic_numbers[10]; 
	int num_isos, znum;
	char result[256];
	if (strncmp(target, "Xe", 10) == 0){
		num_isos = 7;
		znum = 74;
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

double difrate_dER_python(double rhochi, double mass, double logenergy, char* model, char* Target, char*ISO_switch){
  struct difcros_params struct_difcros_T= {Target, 1., mass, 0.5, 220, ISO_switch, 0.0};
  
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
	char result[256];
	if (strncmp(target, "Xe", 10) == 0){
		num_isos = 7;
		znum = 74;
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
