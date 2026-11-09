/*########################################################################################

  _____            _____ _____ _____  _____
 |  __ \     /\   |  __ \_   _|  __ \|  __ \
 | |__) |   /  \  | |__) || | | |  | | |  | |
 |  _  /   / /\ \ |  ___/ | | | |  | | |  | |
 | | \ \  / ____ \| |    _| |_| |__| | |__| |
 |_|  \_\/_/    \_\_|   |_____|_____/|_____/

Based on : arXiv:1802.03174
file : atomic_resp.c
Reads the tabulated atomic responses W_1 ... W_4(k', q) of Xenon for the DM-electron
scattering (arXiv:1912.08204) stored in RAPIDD_for_DM/lib/atomic_responses.

########################################################################################*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "atomic_resp.h"

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

static double W[N_ATOMIC_RESP][N_SHELL_XE][GRIDSIZE_K][GRIDSIZE_Q];
static int W_loaded[N_ATOMIC_RESP] = {0, 0, 0, 0};

static void skip_line(FILE *f){
	int ch;
	while ((ch = fgetc(f)) != EOF && ch != '\n') {}
}

int read_atomic_responses(const char* dir, int response_index){
	/*Reads the 11 Xenon shell files Xe_<shell>.txt of the atomic response W_<response_index>
	from dir (e.g. .../lib/atomic_responses/atomic_response_1) into the global table.
	Each file has 3 comment lines (//) followed by GRIDSIZE_K rows (k') of GRIDSIZE_Q
	columns (q).
	Inputs: dir (directory path), response_index (1..4).
	Output: 0 on success, 1 for a wrong response index, 2 if a file is missing,
	3 if a file is incomplete.*/
	if (response_index < 1 || response_index > N_ATOMIC_RESP) {
		fprintf(stderr, "read_atomic_responses: response index must be in [1, %d]\n", N_ATOMIC_RESP);
		return 1;
	}
	int r = response_index - 1;
	W_loaded[r] = 0;
	char path[1024];
	for (int s = 0; s < N_SHELL_XE; s++) {
		snprintf(path, sizeof(path), "%s/Xe_%s.txt", dir, xe_shells[s].name);
		FILE *f = fopen(path, "r");
		if (f == NULL) {
			fprintf(stderr, "read_atomic_responses: cannot open %s\n", path);
			return 2;
		}
		for (int l = 0; l < 3; l++) skip_line(f);
		for (int ik = 0; ik < GRIDSIZE_K; ik++) {
			for (int iq = 0; iq < GRIDSIZE_Q; iq++) {
				if (fscanf(f, "%lf", &W[r][s][ik][iq]) != 1) {
					fprintf(stderr, "read_atomic_responses: %s is incomplete\n", path);
					fclose(f);
					return 3;
				}
			}
		}
		fclose(f);
	}
	W_loaded[r] = 1;
	return 0;
}

double atomic_response(int resp, int shell, int ik, int iq){
	/*Returns W_resp(k'_ik, q_iq) of the selected shell (0-based shell, ik, iq; resp = 1..4).*/
	return W[resp - 1][shell][ik][iq];
}

int atomic_responses_loaded(int resp){
	/*Returns 1 if W_resp (resp = 1..4) has been read, 0 otherwise.*/
	if (resp < 1 || resp > N_ATOMIC_RESP) return 0;
	return W_loaded[resp - 1];
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
