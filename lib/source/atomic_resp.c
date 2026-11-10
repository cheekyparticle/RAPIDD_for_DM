/*########################################################################################

  _____            _____ _____ _____  _____
 |  __ \     /\   |  __ \_   _|  __ \|  __ \
 | |__) |   /  \  | |__) || | | |  | | |  | |
 |  _  /   / /\ \ |  ___/ | | | |  | | |  | |
 | | \ \  / ____ \| |    _| |_| |__| | |__| |
 |_|  \_\/_/    \_\_|   |_____|_____/|_____/

Based on : arXiv:1802.03174
file : atomic_resp.c
Reads the tabulated atomic responses W_1 ... W_4(k', q) of Xenon and Argon for the DM-electron
scattering (arXiv:1912.08204) stored in RAPIDD_for_DM/lib/atomic_responses.

########################################################################################*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "atomic_resp.h"
#include "phys_consts.h"

/* Xenon shells in the order used by the legacy MadDM code (ls order of the files).
Binding energies in GeV; n_secondary is the number of secondary quanta of Table II of
arXiv:1703.00910 (used for XENON10), e_deexcitation (GeV) the energy of the de-excitation photon
from which they come: the hole of 4d, 4p, 4s is filled by an electron of 5p, 4d, 4p, and
n_secondary = floor(e_deexcitation / 13.8 eV). It is used by the "electron" energy convention of the
XENON1T response (primary electron plus de-excitation photon, as XENON1T and XENONnT). No secondary
quanta are given for the other shells (5p, 5s: e_deexcitation < 13.8 eV; inner shells: not in Table
II, conservatively none). */
#define XE_EB_4D 75.58972072252894e-9
#define XE_EB_4P 163.49493390058458e-9
#define XE_EB_4S 213.7805688618793e-9
#define XE_EB_5P 12.44330433672413e-9
const struct xe_shell xe_shells[N_SHELL_XE] = {
	{"1s", 33317.56055656222e-9, 0, 0.},
	{"2p", 4837.706588171414e-9, 0, 0.},
	{"2s", 5149.213639792182e-9, 0, 0.},
	{"3d", 710.7303605534998e-9, 0, 0.},
	{"3p", 958.4299495823894e-9, 0, 0.},
	{"3s", 1093.2351842564003e-9, 0, 0.},
	{"4d", XE_EB_4D, 4, XE_EB_4D - XE_EB_5P},
	{"4p", XE_EB_4P, 6, XE_EB_4P - XE_EB_4D},
	{"4s", XE_EB_4S, 3, XE_EB_4S - XE_EB_4P},
	{"5p", XE_EB_5P, 0, 0.},
	{"5s", 25.69862365041479e-9, 0, 0.},
};

/* Argon shells (all the occupied orbitals, ls order of the files). Binding energies in GeV: the
Roothaan-Hartree-Fock orbital energies of the DarkART data files (data/Ar_<shell>.txt, in Hartree,
from Bunge et al. 1993) times the Hartree energy 27.211386245988 eV, as for Xenon. No secondary
quanta or de-excitation photons (no detector response for Argon). */
const struct xe_shell ar_shells[N_SHELL_AR] = {
	{"1s", 3227.5520194104365e-9, 0, 0.},
	{"2p", 260.4528038435693e-9, 0, 0.},
	{"2s", 335.3028374537735e-9, 0, 0.},
	{"3p", 16.082364653558844e-9, 0, 0.},
	{"3s", 34.758518644085264e-9, 0, 0.},
};

/* Target data: name (file prefix), shells, mean atomic mass (u) */
static const char * const target_names[N_ELECTRON_TARGETS] = {"Xe", "Ar"};
static const struct xe_shell * const target_shells[N_ELECTRON_TARGETS] = {xe_shells, ar_shells};
static const int target_n_shell[N_ELECTRON_TARGETS] = {N_SHELL_XE, N_SHELL_AR};

static double W[N_ELECTRON_TARGETS][N_ATOMIC_RESP][N_SHELL_MAX][GRIDSIZE_K][GRIDSIZE_Q];
static int W_loaded[N_ELECTRON_TARGETS][N_ATOMIC_RESP] = {{0, 0, 0, 0}, {0, 0, 0, 0}};

static int valid_target(int target){
	return target >= 0 && target < N_ELECTRON_TARGETS;
}

int electron_target_index(const char *target){
	/*Returns the index of the target ("Xe": ELECTRON_TARGET_XE, "Ar": ELECTRON_TARGET_AR), -1 if unknown.*/
	for (int t = 0; t < N_ELECTRON_TARGETS; t++) {
		if (strcmp(target, target_names[t]) == 0) return t;
	}
	return -1;
}

const char * electron_target_name(int target){
	/*Returns the name ("Xe", "Ar") of the target index, NULL if unknown.*/
	return valid_target(target) ? target_names[target] : NULL;
}

double electron_target_mass_au(int target){
	/*Returns the mean atomic mass (u) of the target, 0 if unknown.*/
	switch (target) {
		case ELECTRON_TARGET_XE: return m_Xe_avg_au;
		case ELECTRON_TARGET_AR: return m_Ar_avg_au;
		default: return 0.;
	}
}

static void skip_line(FILE *f){
	int ch;
	while ((ch = fgetc(f)) != EOF && ch != '\n') {}
}

int read_atomic_responses_target(const char* dir, int response_index, int target){
	/*Reads the shell files <target>_<shell>.txt (Xe: 11 shells, Ar: 5 shells) of the atomic
	response W_<response_index> from dir (e.g. .../lib/atomic_responses/atomic_response_1) into
	the global table of the target. Each file has 3 comment lines (//) followed by GRIDSIZE_K
	rows (k') of GRIDSIZE_Q columns (q).
	Inputs: dir (directory path), response_index (1..4), target (ELECTRON_TARGET_XE, ELECTRON_TARGET_AR).
	Output: 0 on success, 1 for a wrong response index or target, 2 if a file is missing,
	3 if a file is incomplete.*/
	if (response_index < 1 || response_index > N_ATOMIC_RESP) {
		fprintf(stderr, "read_atomic_responses: response index must be in [1, %d]\n", N_ATOMIC_RESP);
		return 1;
	}
	if (!valid_target(target)) {
		fprintf(stderr, "read_atomic_responses: unknown target %d\n", target);
		return 1;
	}
	int r = response_index - 1;
	W_loaded[target][r] = 0;
	char path[1024];
	for (int s = 0; s < target_n_shell[target]; s++) {
		snprintf(path, sizeof(path), "%s/%s_%s.txt", dir, target_names[target], target_shells[target][s].name);
		FILE *f = fopen(path, "r");
		if (f == NULL) {
			fprintf(stderr, "read_atomic_responses: cannot open %s\n", path);
			return 2;
		}
		for (int l = 0; l < 3; l++) skip_line(f);
		for (int ik = 0; ik < GRIDSIZE_K; ik++) {
			for (int iq = 0; iq < GRIDSIZE_Q; iq++) {
				if (fscanf(f, "%lf", &W[target][r][s][ik][iq]) != 1) {
					fprintf(stderr, "read_atomic_responses: %s is incomplete\n", path);
					fclose(f);
					return 3;
				}
			}
		}
		fclose(f);
	}
	W_loaded[target][r] = 1;
	return 0;
}

int read_atomic_responses(const char* dir, int response_index){
	/*Reads the 11 Xenon shell files Xe_<shell>.txt of the atomic response W_<response_index>
	(see read_atomic_responses_target).*/
	return read_atomic_responses_target(dir, response_index, ELECTRON_TARGET_XE);
}

double atomic_response_target(int target, int resp, int shell, int ik, int iq){
	/*Returns W_resp(k'_ik, q_iq) of the selected shell of the target (0-based shell, ik, iq; resp = 1..4).*/
	return W[target][resp - 1][shell][ik][iq];
}

double atomic_response(int resp, int shell, int ik, int iq){
	/*Returns W_resp(k'_ik, q_iq) of the selected Xenon shell (0-based shell, ik, iq; resp = 1..4).*/
	return W[ELECTRON_TARGET_XE][resp - 1][shell][ik][iq];
}

int atomic_responses_loaded_target(int target, int resp){
	/*Returns 1 if W_resp (resp = 1..4) of the target has been read, 0 otherwise.*/
	if (!valid_target(target) || resp < 1 || resp > N_ATOMIC_RESP) return 0;
	return W_loaded[target][resp - 1];
}

int atomic_responses_loaded(int resp){
	/*Returns 1 if W_resp (resp = 1..4) of Xenon has been read, 0 otherwise.*/
	return atomic_responses_loaded_target(ELECTRON_TARGET_XE, resp);
}

int n_shell_target(int target){
	/*Returns the number of shells of the target, 0 if unknown.*/
	return valid_target(target) ? target_n_shell[target] : 0;
}

double shell_binding_energy_target(int target, int shell){
	/*Returns the binding energy (GeV) of the selected shell (0-based) of the target.*/
	return target_shells[target][shell].binding_energy;
}

const char * shell_name_target(int target, int shell){
	/*Returns the name ("1s", "2p", ...) of the selected shell (0-based) of the target.*/
	return target_shells[target][shell].name;
}

int n_shell_xe(){
	/*Returns the number of Xenon shells.*/
	return N_SHELL_XE;
}

double shell_binding_energy(int shell){
	/*Returns the binding energy (GeV) of the selected shell (0-based).*/
	return xe_shells[shell].binding_energy;
}

int shell_n_secondary(int shell){
	/*Returns the number of secondary quanta (XENON10) of the selected shell (0-based).*/
	return xe_shells[shell].n_secondary;
}

double shell_deexcitation_energy(int shell){
	/*Returns the energy (GeV) of the de-excitation photon of the selected shell (0-based), 0 if none.*/
	return xe_shells[shell].e_deexcitation;
}

const char * shell_name(int shell){
	/*Returns the name ("1s", "2p", ...) of the selected shell (0-based).*/
	return xe_shells[shell].name;
}
