/*Crystal responses of the DM-electron rate in semiconductors (crystal_resp.c)*/

#pragma once

#define N_CRYSTAL_TARGETS 1      // targets of the crystal DM-electron rate: 0 = Si
#define CRYSTAL_TARGET_SI 0
// grid of the tables (QCDark2 Si_comp.h5): q = CRYSTAL_Q0_AME + j CRYSTAL_DQ_AME (alpha m_e), j < CRYSTAL_NQ,
// Delta E = i CRYSTAL_DE_EV (eV), i < CRYSTAL_NE
#define CRYSTAL_NQ 1251
#define CRYSTAL_NE 501
#define CRYSTAL_Q0_AME 0.01
#define CRYSTAL_DQ_AME 0.02
#define CRYSTAL_DE_EV 0.1
#define CRYSTAL_UNSCREENED 0     // W_1 from Im eps (with local field effects, LFE)
#define CRYSTAL_SCREENED 1       // W_1 from Im[-1/eps] (RPA screening, with LFE)
#define CRYSTAL_NOLFE 2          // W_1 from Im eps without LFE (independent-particle response, spin couplings)
#define N_CRYSTAL_TABLES 3

int crystal_target_index(const char *target);
const char * crystal_target_name(int target);
double crystal_cell_mass_au(int target);
int read_crystal_responses(const char *dir, int target);
int crystal_responses_loaded(int target);
double crystal_response(int target, int table, int iE, int iq);
double crystal_q_gev(int iq, double m_e);
double crystal_dE_gev(int iE);
