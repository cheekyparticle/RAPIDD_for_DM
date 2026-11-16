/*#########################################

  _____            _____ _____ _____  _____  
 |  __ \     /\   |  __ \_   _|  __ \|  __ \ 
 | |__) |   /  \  | |__) || | | |  | | |  | |
 |  _  /   / /\ \ |  ___/ | | | |  | | |  | |
 | | \ \  / ____ \| |    _| |_| |__| | |__| |
 |_|  \_\/_/    \_\_|   |_____|_____/|_____/ 
                                             
Based on : arXiv:1802.03174
Physical constants
##########################################*/
#include <math.h>
#include<stdio.h>
#include<stdlib.h>

#include "phys_consts.h"


/////Frequently used equations///////////

double reduced_mass(double mass_1, double mass_2){

	return  mass_1*mass_2 / (mass_1 + mass_2);

}

/* Atomic masses in u of the isotopes of the targets: AME2020 (Wang et al., Chinese Phys. C 45, 030003 (2021)),
as tabulated by CIAAW (ciaaw.org/atomic-masses.htm). */
static const struct { int A; int Z; double mass_u; } atomic_masses[] = {
	{16, 8, 15.994914630}, {19, 9, 18.998403224}, {40, 18, 39.962383681}, {40, 20, 39.962590588},
	{70, 32, 69.924249656}, {72, 32, 71.922078859}, {73, 32, 72.923462608}, {74, 32, 73.921177394},
	{76, 32, 75.921401565},
	{124, 54, 123.905894239}, {126, 54, 125.904280981}, {128, 54, 127.903531222}, {129, 54, 128.904780103},
	{130, 54, 129.903509426}, {131, 54, 130.905072012}, {132, 54, 131.904144154}, {134, 54, 133.905394991},
	{136, 54, 135.907213743},
	{0, 0, 0.}
};

double nucleus_mass(int A, int Z){
	/*Mass of the nucleus (A, Z) in GeV: atomic mass (AME2020) minus the mass of the Z electrons (their
	binding energy, < 1e-6 of the mass, is neglected). It is the m_T of q^2 = 2 m_T E_R, of the DM-nucleus
	reduced mass and of vmin. The sum of the masses of the free nucleons, Z m_p + (A - Z) m_n, is ~0.9%
	larger (nuclear binding energy).
	For an isotope that is not tabulated the atomic mass is approximated by A u (error < 0.1%).
	Inputs: A (mass number), Z (atomic number).*/
	double mass_u = (double)A;
	for (int i = 0; atomic_masses[i].A != 0; i++) {
		if (atomic_masses[i].A == A && atomic_masses[i].Z == Z) {
			mass_u = atomic_masses[i].mass_u;
			break;
		}
	}
	return mass_u * atomic_mass_unit_GeV - Z * melectron;
}

void isotope_mass_fractions(const int *A, int Z, const double *abundance, int n, double *mass_fraction){
	/*Mass fractions of the isotopes of an element from their mole (number) fractions:
	w_i = f_i m_i / m_avg, with m_avg = sum_j f_j m_j / sum_j f_j over the n listed isotopes.
	The differential rate of each isotope is per kg of that isotope, so the rate per kg of the element
	is sum_i w_i R_i. sum_i w_i = sum_i f_i: isotopes that are not listed (no form factors) do not
	contribute to the rate, but their mass fraction is not given to the others.
	Inputs: A[n] (mass numbers), Z, abundance[n] (mole fractions). Output: mass_fraction[n].*/
	double sum_f = 0., sum_fm = 0.;
	for (int l = 0; l < n; l++) {
		sum_f += abundance[l];
		sum_fm += abundance[l] * nucleus_mass(A[l], Z);
	}
	for (int l = 0; l < n; l++) {
		mass_fraction[l] = (sum_fm > 0.) ? abundance[l] * nucleus_mass(A[l], Z) * sum_f / sum_fm : 0.;
	}
}

////////Function to calculate couplings f (f_n, f_p) for SI interaction by using the cross section sigma//////////////////////////
double sigma_to_f_11(double sigma, double mchi, double mnucleon){

	double mun = reduced_mass(mchi, mnucleon); //reduced mass DM particle and nucleon (mproton or mneutron)
	return sqrt(sigma*M_PI) / (2.*mun);

}

////////Function to calculate couplings a (a_n, a_p) for SD interaction by using the cross section sigma//////////////////////////
double sigma_to_a_44(double sigma, double mchi, double mnucleon){

	double mun = reduced_mass(mchi, mnucleon); //reduced mass DM particle and nucleon (mproton or mneutron)
	return sqrt(sigma*M_PI) / (sqrt(24.)*mun*G_F);

}


int countlines(char *filename){
  // count the number of lines in the file called filename
  FILE *fp = fopen(filename,"r");
  int ch=0;
  int lines=0;

  if (fp == NULL){
  return 0;
  }


  while(!feof(fp))
  {
    ch = fgetc(fp);
    if(ch == '\n')
    {
      lines++;
    }
  }

  fclose(fp);
	//printf("attempt %i\n", lines );
  return lines;
}
