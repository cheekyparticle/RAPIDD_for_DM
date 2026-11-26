/*Crystal response W_1(q, Delta E) of semiconductor targets for the DM-electron rate (arXiv:2105.02233,
Eq. R_crystal_2D: the rate of the atomic case with the unit cell as target and the deposited energy Delta E).
Silicon: tables lib/crystal_responses/Si_W1.txt (unscreened, Im eps), Si_W1_screened.txt (RPA, Im[-1/eps]) and
Si_W1_nolfe.txt (Im eps without local field effects), from the QCDark2 dielectric functions (C. Dreyer, R. Essig,
M. Fernandez-Serra, M. Hott, A. Singal, arXiv:2603.12326, github.com/meganhott/QCDark2, CC0), see
lib/crystal_responses/README.md. Crystal rate formalism: R. Catena, T. Emken, M. Matas, N. A. Spaldin,
E. Urdshals, arXiv:2105.02233.*/

#include <stdio.h>
#include <string.h>
#include "crystal_resp.h"
#include "phys_consts.h"

static const char * const crystal_names[N_CRYSTAL_TARGETS] = {"Si"};
// unit cell mass (u): 2 atoms, standard atomic weight of Si
static const double crystal_cell_au[N_CRYSTAL_TARGETS] = {2. * 28.0855};
static const char * const crystal_files[N_CRYSTAL_TABLES] = {"W1", "W1_screened", "W1_nolfe"};

static double Wc[N_CRYSTAL_TARGETS][N_CRYSTAL_TABLES][CRYSTAL_NE][CRYSTAL_NQ];
static int Wc_loaded[N_CRYSTAL_TARGETS] = {0};

static int valid_crystal(int target){
	return target >= 0 && target < N_CRYSTAL_TARGETS;
}

int crystal_target_index(const char *target){
	/*Returns the index of the crystal target ("Si": CRYSTAL_TARGET_SI), -1 if unknown.*/
	for (int t = 0; t < N_CRYSTAL_TARGETS; t++) {
		if (strcmp(target, crystal_names[t]) == 0) return t;
	}
	return -1;
}

const char * crystal_target_name(int target){
	/*Returns the name of the crystal target, NULL if unknown.*/
	return valid_crystal(target) ? crystal_names[target] : NULL;
}

double crystal_cell_mass_au(int target){
	/*Returns the mass (u) of the unit cell of the crystal target, 0 if unknown.*/
	return valid_crystal(target) ? crystal_cell_au[target] : 0.;
}

static void skip_comments(FILE *f){
	/*Skips the header lines starting with "//".*/
	int ch;
	while ((ch = fgetc(f)) == '/') {
		while ((ch = fgetc(f)) != EOF && ch != '\n') {}
	}
	if (ch != EOF) ungetc(ch, f);
}

int read_crystal_responses(const char *dir, int target){
	/*Reads <dir>/<target>_W1.txt, <target>_W1_screened.txt and <target>_W1_nolfe.txt: header lines starting with "//",
	then CRYSTAL_NE rows (Delta E) of CRYSTAL_NQ columns (q).
	Output: 0 on success, 1 for an unknown target, 2 if a file is missing, 3 if a file is incomplete.*/
	if (!valid_crystal(target)) {
		fprintf(stderr, "read_crystal_responses: unknown target %d\n", target);
		return 1;
	}
	Wc_loaded[target] = 0;
	char path[1024];
	for (int s = 0; s < N_CRYSTAL_TABLES; s++) {
		snprintf(path, sizeof(path), "%s/%s_%s.txt", dir, crystal_names[target], crystal_files[s]);
		FILE *f = fopen(path, "r");
		if (f == NULL) {
			fprintf(stderr, "read_crystal_responses: cannot open %s\n", path);
			return 2;
		}
		skip_comments(f);
		for (int iE = 0; iE < CRYSTAL_NE; iE++) {
			for (int iq = 0; iq < CRYSTAL_NQ; iq++) {
				if (fscanf(f, "%lf", &Wc[target][s][iE][iq]) != 1) {
					fprintf(stderr, "read_crystal_responses: %s is incomplete\n", path);
					fclose(f);
					return 3;
				}
			}
		}
		fclose(f);
	}
	Wc_loaded[target] = 1;
	return 0;
}

int crystal_responses_loaded(int target){
	/*Returns 1 if the tables of the crystal target have been read.*/
	return valid_crystal(target) ? Wc_loaded[target] : 0;
}

double crystal_response(int target, int table, int iE, int iq){
	/*W_1(q_iq, Delta E_iE) of the crystal target from the table CRYSTAL_UNSCREENED, CRYSTAL_SCREENED or
	CRYSTAL_NOLFE.*/
	return Wc[target][table][iE][iq];
}

double crystal_q_gev(int iq, double m_e){
	/*q (GeV) of the point iq of the grid of the crystal tables.*/
	return (CRYSTAL_Q0_AME + CRYSTAL_DQ_AME * (double)iq) * alpha_em * m_e;
}

double crystal_dE_gev(int iE){
	/*Delta E (GeV) of the point iE of the grid of the crystal tables.*/
	return CRYSTAL_DE_EV * (double)iE * 1.e-9;
}
