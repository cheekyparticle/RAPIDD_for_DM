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
Binding energies in GeV; n_secondary is the number of secondary quanta of
Table II of arXiv:1703.00910 (used for XENON10). */
const struct xe_shell xe_shells[N_SHELL_XE] = {
	{"1s", 33317.56055656222e-9, 0},
	{"2p", 4837.706588171414e-9, 0},
	{"2s", 5149.213639792182e-9, 0},
	{"3d", 710.7303605534998e-9, 0},
	{"3p", 958.4299495823894e-9, 0},
	{"3s", 1093.2351842564003e-9, 0},
	{"4d", 75.58972072252894e-9, 4},
	{"4p", 163.49493390058458e-9, 6},
	{"4s", 213.7805688618793e-9, 3},
	{"5p", 12.44330433672413e-9, 0},
	{"5s", 25.69862365041479e-9, 0},
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

const char * shell_name(int shell){
	/*Returns the name ("1s", "2p", ...) of the selected shell (0-based).*/
	return xe_shells[shell].name;
}
