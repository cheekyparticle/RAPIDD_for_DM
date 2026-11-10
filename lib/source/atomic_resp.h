/*Please include all definitions of the functions and variables used in atomic_resp.c
  that you can use in other parts of the code (needed for correct linking)*/

#pragma once

#define N_ATOMIC_RESP 4   // W_1 ... W_4 of arXiv:1912.08204
#define N_SHELL_XE 11     // Xenon shells
#define N_SHELL_AR 5      // Argon shells
#define N_SHELL_MAX 11    // largest number of shells of a target
#define N_ELECTRON_TARGETS 2 // targets of the DM-electron rate: 0 = Xe, 1 = Ar
#define ELECTRON_TARGET_XE 0
#define ELECTRON_TARGET_AR 1
#define GRIDSIZE_K 100    // outgoing electron momentum k' grid
#define GRIDSIZE_Q 200    // exchanged momentum q grid

// k' and q ranges of the tabulated atomic responses, in GeV
#define K_MIN_GEV 1.0e-7
#define K_MAX_GEV 1.0e-4
#define Q_MIN_GEV 1.0e-6
#define Q_MAX_GEV 1.0e-3

struct xe_shell { const char * name; double binding_energy; int n_secondary; double e_deexcitation; };
extern const struct xe_shell xe_shells[N_SHELL_XE];
extern const struct xe_shell ar_shells[N_SHELL_AR];

int read_atomic_responses(const char* dir, int response_index);
double atomic_response(int resp, int shell, int ik, int iq);
int atomic_responses_loaded(int resp);
int n_shell_xe();
double shell_binding_energy(int shell);
int shell_n_secondary(int shell);
double shell_deexcitation_energy(int shell);
const char * shell_name(int shell);

// same as above for the target (ELECTRON_TARGET_XE, ELECTRON_TARGET_AR); the functions without
// target refer to Xenon
int electron_target_index(const char *target);
const char * electron_target_name(int target);
double electron_target_mass_au(int target);
int read_atomic_responses_target(const char* dir, int response_index, int target);
double atomic_response_target(int target, int resp, int shell, int ik, int iq);
int atomic_responses_loaded_target(int target, int resp);
int n_shell_target(int target);
double shell_binding_energy_target(int target, int shell);
const char * shell_name_target(int target, int shell);
