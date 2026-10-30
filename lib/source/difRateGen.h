#pragma once
#include <gsl/gsl_spline.h>

struct difcros_params { char * target; double Er; double mchi; double jchi; double v; char * Nucleon; double exposure; double delta; };

double difrate_SI_standard(double rho, int A, int Z, double fp, double fn, double mchi, double Er);
double total_difrate_isotope_dEr(int A, int Z, double rhochi, void * input_difcros);
double difrate_dER(double rhochi, void * input_difcros, double logenergy, char* model);
int nuclear_model_supports_target(char* NR_framework, char* target);
double difrate_dER_2(double rhochi, void * input_difcros, double energy, char* model);
double difrate_isotope_dEr_w(int A, int Z, double rhochi, void * input_difcros, int F_i, int F_j, gsl_interp_accel *ga, gsl_spline * gs, gsl_interp_accel *ga2, gsl_spline * gs2);
double difrate_isotope_v2_dEr_w(int A, int Z, double rhochi, void * input_difcros, int F_i, int F_j, gsl_interp_accel *ga, gsl_spline * gs);
double difrate_isotope_v0_dEr_w(int A, int Z, double rhochi, void * input_difcros, int F_i, int F_j, gsl_interp_accel *ga, gsl_spline * gs);

double difrate_dER_w(double rhochi, void * input_difcros, double logenergy, char* model, gsl_interp_accel *ga, gsl_spline * gs, gsl_interp_accel *ga2, gsl_spline * gs2);
double total_difrate_isotope_dEr_w(int A, int Z, double rhochi, void * input_difcros, gsl_interp_accel *ga, gsl_spline * gs, gsl_interp_accel *ga2, gsl_spline * gs2);
double difrate_dER_python(double rhochi, double mass, double logenergy, char* model, char* Target, char*ISO_switch, double delta);
double difrate_dER_python_2(double rhochi, double mass, double energy, char* model, char* Target, char*ISO_switch);

// DM-electron scattering
double dm_response_electron(int resp, double q, double mchi, double jchi);
double difrate_electron_dEr(int shell, int ik, double rhochi, double mchi, double jchi, double m_e);
int difrate_electron_dEr_all(double rhochi, double mchi, double jchi, double m_e, double *k_out, double *E_out, double *rate_shell_out);
int xenon10_bins(int *bins_out);
int difrate_dS2_Xenon10(const double *E, const double *rate_shell, const char *exp_dir,
						double *S2_out, double *dRdS2_out, double *dRdS2_shell_out, double *signal_out, int *n_obs_out);
int difrate_dS2_Xenon1T(const double *E, const double *rate_shell, double mchi, const char *exp_dir,
						double *S2_out, double *dRdS2_out, double *signal_out, double *n_bkg_out, int *n_obs_out);
int difrate_dS2(const char *experiment, const double *E, const double *rate_shell, double mchi, const char *exp_dir,
				double *S2_out, double *dRdS2_out, double *dRdS2_shell_out, double *signal_out, double *n_bkg_out, int *n_obs_out);
#ifdef ER_VALIDATE_LEGACY_HALO
void set_legacy_halo_params(double vMP, double vesc, double vEarth_mod, double vSun_r, double vSun_phi,
							double vSun_theta, double day_min, double day_max, int day_bins);
#endif
