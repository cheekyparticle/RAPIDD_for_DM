/*Please include all definitions of the functions and variables used in atomic_resp.c
  that you can use in other parts of the code (needed for correct linking)*/

#pragma once

#define N_ATOMIC_RESP 4   // W_1 ... W_4 of arXiv:1912.08204
#define N_SHELL_XE 11     // Xenon shells
#define GRIDSIZE_K 100    // outgoing electron momentum k' grid
#define GRIDSIZE_Q 200    // exchanged momentum q grid

// k' and q ranges of the tabulated atomic responses, in GeV
#define K_MIN_GEV 1.0e-7
#define K_MAX_GEV 1.0e-4
#define Q_MIN_GEV 1.0e-6
#define Q_MAX_GEV 1.0e-3

struct xe_shell { const char * name; double binding_energy; int n_secondary; };
extern const struct xe_shell xe_shells[N_SHELL_XE];

int read_atomic_responses(const char* dir, int response_index);
double atomic_response(int resp, int shell, int ik, int iq);
int atomic_responses_loaded(int resp);
int n_shell_xe();
double shell_binding_energy(int shell);
int shell_n_secondary(int shell);
const char * shell_name(int shell);
