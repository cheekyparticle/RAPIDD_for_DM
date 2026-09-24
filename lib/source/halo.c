/*########################################################################################

  _____            _____ _____ _____  _____
 |  __ \     /\   |  __ \_   _|  __ \|  __ \
 | |__) |   /  \  | |__) || | | |  | | |  | |
 |  _  /   / /\ \ |  ___/ | | | |  | | |  | |
 | | \ \  / ____ \| |    _| |_| |__| | |__| |
 |_|  \_\/_/    \_\_|   |_____|_____/|_____/

Based on : arXiv:1802.03174
file : halo.c
Produces the DM halo integral for different halo profiles and stores them in the
RAPIDD_for_DM/lib/halo_table folder.
It uses 2D integrations to change from galactic to Earth (lab) coordinates.

Early versions of this code were developed by M. Peiro
and E. Gerstmayr.

########################################################################################*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h> // Contains gamma and beta functions.
#include <unistd.h>           // for access check
#include <time.h>             // for elapsed time
#include <sys/stat.h>         // for mkdir
#include "halo.h"

#define SQR(X) ((X) * (X))
#define ABS(X) ((X) > 0 ? (X) : (-(X)))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "../cubature.h"

#define VERBOSE 0
#define length 200 // Number of divisions between 0 and (vesc+ve) to perform the interpolation.
#define power 10   // Maximum predefined power of the velocity in the halo integral

#if defined(PCUBATURE)
#  define cubature pcubature
#else
#  define cubature hcubature
#endif

float vel[length];
float fv[power][length];

#define Narray 50
double pQ[Narray + 1];
double pQa[Narray + 1];


/*########################################################################################
STANDARD HALO MODEL (SHM)
Analytical integrals for i = 0, 1, 2, where i-1 is the power of v in the integrand.
The integrals are normalized to 1.
Closed-form primitives of the boosted, truncated Maxwellian speed density
    f1(v) = (v/(sqrt(pi) ve v0 Nesc)) * [exp(-(v-ve)^2/v0^2) - exp(-(v+ve)^2/v0^2)]   v < vesc-ve
    f1(v) = (v/(sqrt(pi) ve v0 Nesc)) * [exp(-(v-ve)^2/v0^2) - exp(-vesc^2/v0^2)]     vesc-ve<v<vesc+ve
(Savage, Freese & Gondolo, astro-ph/0607121; McCabe, arXiv:1005.0579)
########################################################################################*/

static double shm_p1(double v, double ve, double v0){
	/*Primitive of the lower-branch exp term, used for i=1.*/
	return -0.5 * v0 * v0 * exp(-SQR(v - ve) / (v0 * v0))
	       + 0.5 * ve * v0 * sqrt(M_PI) * gsl_sf_erf((v - ve) / v0);
}
static double shm_q1(double v, double ve, double v0){
	/*Primitive of the upper-branch (v+ve) exp term, used for i=1.*/
	return -0.5 * v0 * v0 * exp(-SQR(v + ve) / (v0 * v0))
	       - 0.5 * ve * v0 * sqrt(M_PI) * gsl_sf_erf((v + ve) / v0);
}
static double shm_p2(double v, double ve, double v0){
	/*Primitive of the lower-branch exp term, used for i=2.*/
	return -0.5 * v0 * v0 * (v + ve) * exp(-SQR(v - ve) / (v0 * v0))
	       + sqrt(M_PI) * v0 * (0.25 * v0 * v0 + 0.5 * ve * ve) * gsl_sf_erf((v - ve) / v0);
}
static double shm_q2(double v, double ve, double v0){
	/*Primitive of the upper-branch (v+ve) exp term, used for i=2.*/
	return -0.5 * v0 * v0 * (v - ve) * exp(-SQR(v + ve) / (v0 * v0))
	       + sqrt(M_PI) * v0 * (0.25 * v0 * v0 + 0.5 * ve * ve) * gsl_sf_erf((v + ve) / v0);
}

double shm_halo(double vmin, double vesc, double v0, double ve, int i){
	/*Fully analytic SHM halo integral (hard Theta cutoff, beta=0) for i=0,1,2.
	Inputs: vmin, vesc, v0, ve, i (velocity power index; only 0/1/2 implemented).
	Output: normalized halo integral (0 if vmin>vesc+ve or i not implemented).*/
	if (vmin > vesc + ve) return 0.;
	if (i < 0 || i > 2) {
		printf("shm_halo: only i = 0, 1, 2 are implemented analytically\n");
		return 0.;
	}

	double xesc   = vesc / v0;
	double Nesc   = gsl_sf_erf(xesc) - (2. / sqrt(M_PI)) * xesc * exp(-xesc * xesc);
	double expesc = exp(-xesc * xesc);
	double vlo = vesc - ve, vhi = vesc + ve;

	/* i = 0: eta(vmin), same closed form as the pre-existing beta-smoothed formula
	   at beta = 0 (McCabe, arXiv:1005.0579) */
	if (i == 0) {
		double xmin = vmin / v0;
		double xe   = ve / v0;
		double norm3d = pow(M_PI, 1.5) * pow(v0, 3.) * Nesc; /* == norm_shm(vesc,v0,0.0) */

		if ((xe + xmin) < xesc) {
			return pow(M_PI, 1.5) * pow(v0, 2.) / (2. * xe * norm3d)
			       * (gsl_sf_erf(xmin + xe) - gsl_sf_erf(xmin - xe) - 4. * xe / sqrt(M_PI) * expesc);
		}
		if (xmin > fabs(xesc - xe) && (xe + xesc) > xmin) {
			return pow(M_PI, 1.5) * pow(v0, 2.) / (2. * xe * norm3d)
			       * (gsl_sf_erf(xesc) + gsl_sf_erf(xe - xmin) - 2. / sqrt(M_PI) * expesc * (xesc + xe - xmin));
		}
		if (xe > (xmin + xesc)) return 1. / (xe * v0);
		return 0.;
	}

	/* i = 1, 2: int f1(v) dv and int v*f1(v) dv, N = 1/Nesc normalizes f1 to a pdf */
	double N    = 1. / Nesc;
	double pref = N / (sqrt(M_PI) * ve * v0);
	double K;

	if (i == 1) {
		if (vmin > vlo) {
			K = (shm_p1(vhi, ve, v0) - expesc * vhi * vhi / 2.)
			  - (shm_p1(vmin, ve, v0) - expesc * vmin * vmin / 2.);
		} else {
			K = (shm_p1(vlo, ve, v0) - shm_q1(vlo, ve, v0) - shm_p1(vmin, ve, v0) + shm_q1(vmin, ve, v0))
			  + (shm_p1(vhi, ve, v0) - expesc * vhi * vhi / 2. - shm_p1(vlo, ve, v0) + expesc * vlo * vlo / 2.);
		}
	} else { /* i == 2 */
		if (vmin > vlo) {
			K = (shm_p2(vhi, ve, v0) - expesc * vhi * vhi * vhi / 3.)
			  - (shm_p2(vmin, ve, v0) - expesc * vmin * vmin * vmin / 3.);
		} else {
			K = (shm_p2(vlo, ve, v0) - shm_q2(vlo, ve, v0) - shm_p2(vmin, ve, v0) + shm_q2(vmin, ve, v0))
			  + (shm_p2(vhi, ve, v0) - expesc * vhi * vhi * vhi / 3. - shm_p2(vlo, ve, v0) + expesc * vlo * vlo * vlo / 3.);
		}
	}

	return pref * K;
}

/*########################################################################################
STANDARD HALO MODEL (SHM)
Cutoff and smooth distributions controlled by beta. See: 1005.0579 for more info.
Numerical integrals for i = 0, 1, 2, where i-1 is the power of v in the integrand.
########################################################################################*/

double norm_shm(double vesc, double v0, double beta){
	/*Analytic normalization of the (smooth-cutoff) SHM velocity distribution.
	Inputs: vesc (escape velocity), v0 (velocity dispersion), beta (smoothing parameter).
	Output: normalization constant.*/
	double xesc = vesc / v0;
	return pow(M_PI, 1.5) * pow(v0, 3.) * (gsl_sf_erf(xesc) - 4 / sqrt(M_PI)
	       * exp(-xesc * xesc) * (xesc / 2. + beta * pow(xesc, 3.) / 3.));
}

struct Nshm_params { double vesc_Nshm; double v0_Nshm; double beta_Nshm; };

double N_shm(double v, void *p){
	/*Integrand used to numerically evaluate the SHM normalization (see norm_shm_num).
	Inputs: v (velocity), p (pointer to Nshm_params).
	Output: integrand value at v.*/
	struct Nshm_params *params = (struct Nshm_params *)p;
	double vesc_Nshm = (params->vesc_Nshm);
	double v0_Nshm   = (params->v0_Nshm);
	double beta_Nshm = (params->beta_Nshm);
	return 4. * M_PI * v * v * (exp(-(v * v) / (v0_Nshm * v0_Nshm))
	       - beta_Nshm * exp(-vesc_Nshm * vesc_Nshm / (v0_Nshm * v0_Nshm)));
}

double norm_shm_num(double vesc, double v0, double beta){
	/*Numerical normalization of the SHM distribution via GSL QAGS integration.
	Inputs: vesc, v0, beta.
	Output: normalization constant.*/
	double result, error;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

	gsl_function F;
	struct Nshm_params params = {vesc, v0, beta};

	F.function = &N_shm;
	F.params = &params;

	gsl_integration_qags(&F, 0, vesc, 1e-8, 1e-8, 1000, w, &result, &error);
	gsl_integration_workspace_free(w);

	return result;
}

struct shm_params {
	double vesc_param;
	double v0_param;
	double ve_param;
	double beta_param;
	int i_param;
};

int f_shm(unsigned ndim, const double *k, void *p, unsigned fdim, double *fval){
	/*Cubature integrand for the SHM halo integral with a hard (sharp) velocity cutoff.
	Inputs: ndim, k (integration variables: k[0]=cos(theta), k[1]=v), p (shm_params), fdim.
	Output: fval[0] set to 0 below the cutoff angle, else the integrand value; return 0 on success.*/
	struct shm_params *params = (struct shm_params *)p;
	double vesc_param = (params->vesc_param);
	double v0_param   = (params->v0_param);
	double ve_param   = (params->ve_param);
	int i_param       = (params->i_param);
	double v = k[1];
	double cos_theta = k[0];
	double cos_theta_min = (v * v + ve_param * ve_param - vesc_param * vesc_param)
	                        / (2. * v * ve_param);

	if (cos_theta < cos_theta_min) {
		fval[0] = 0.;
		return 0;
	}

	fval[0] = 2. * M_PI * v * pow(v, i_param * 1.0) *
	          exp(-(v * v + ve_param * ve_param - 2. * v * ve_param * cos_theta)
	              / (v0_param * v0_param));

	return 0;
}

int f_shm_ann(unsigned ndim, const double *k, void *p, unsigned fdim, double *fval){
	/*Cubature integrand for the SHM halo integral, alternate ("annual") smooth-cutoff
	form used in arXiv:1112.0524v2 (added 12/02/2016).
	Inputs: ndim, k (integration variables), p (shm_params), fdim.
	Output: fval[0] set to the integrand value; return 0 on success.*/
	struct shm_params *params = (struct shm_params *)p;
	double vesc_param = (params->vesc_param);
	double v0_param   = (params->v0_param);
	double ve_param   = (params->ve_param);
	double beta_param = (params->beta_param);
	int i_param       = (params->i_param);

	fval[0] = 2. * M_PI * k[1] * pow(k[1], i_param * 1.0) *
	          (exp((-3 / 2) * (k[1] * k[1] + ve_param * ve_param - 2. * k[1] * ve_param * k[0]
	           + beta_param * (-vesc_param) * vesc_param) / (v0_param * v0_param)));

	return 0;
}

int f_shm_beta(unsigned ndim, const double *k, void *p, unsigned fdim, double *fval){
	/*Cubature integrand for the SHM halo integral (smooth-cutoff, angle-velocity space).
	Inputs: ndim, k (integration variables: k[0]=cos(theta), k[1]=v), p (shm_params), fdim.
	Output: fval[0] set to the integrand value; return 0 on success.*/
	struct shm_params *params = (struct shm_params *)p;
	double vesc_param = (params->vesc_param);
	double v0_param   = (params->v0_param);
	double ve_param   = (params->ve_param);
	double beta_param = (params->beta_param);
	int i_param       = (params->i_param);

	fval[0] = 2. * M_PI * k[1] * pow(k[1], i_param * 1.0) *
	          (exp(-(k[1] * k[1] + ve_param * ve_param - 2. * k[1] * ve_param * k[0]) / (v0_param * v0_param))
	           - beta_param * exp(-vesc_param * vesc_param / (v0_param * v0_param)));

	/* Edited by Andrew originally 2.*M_PI*k[1] * pow(k[1], i_param*1.0)
	(exp(-(k[1] * k[1] + ve_param*ve_param + 2.*k[1] * ve_param*k[0]) / (v0_param*v0_param))
	 - beta_param*exp(-vesc_param*vesc_param/(v0_param*v0_param))) */

	return 0;
}

double shm_halo_beta(double vmin, double vesc, double v0, double ve, double beta, int i){
	/*SHM halo integral (smooth-cutoff form) via 2D cubature over angle and velocity.
	Inputs: vmin, vesc, v0, ve, beta, i (velocity power).
	Output: value of the halo integral (0 if vmin > vesc+ve).*/
	double val, err;
	struct shm_params params = {vesc, v0, ve, beta, i};

	// limits of integration
	double xl[2] = {-1, vmin};
	double xu[2] = {1, vesc + ve};

	if (vmin > vesc + ve) return 0.;
	else {
		hcubature(1, &f_shm_beta, &params,
		          2, xl, xu,
		          0, 0, 1e-13, ERROR_INDIVIDUAL, &val, &err);
		return val/norm_shm_num(vesc, v0, beta);;
	}
}

struct shm_radial_params { double vesc; double v0; double ve; int i; };

double f_shm_radial(double v, void *p){
	/*Radial integrand for the SHM halo integral, after analytic
	integration over the angular variable (see shm_halo_numeric).
	Inputs: v (velocity), p (shm_radial_params).
	Output: integrand value at v.*/
	struct shm_radial_params *params = (struct shm_radial_params *)p;
	double vesc = params->vesc;
	double v0   = params->v0;
	double ve   = params->ve;
	int i       = params->i;

	/* avoid division by zero near v=0; true integrand -> 0 there anyway */
	if (v < 1e-6) return 0.0;

	double u_min_raw = (v * v + ve * ve - vesc * vesc) / (2. * v * ve);
	double u_min = (u_min_raw > -1.) ? u_min_raw : -1.;

	double exp_lower = exp(-(v - ve) * (v - ve) / (v0 * v0));
	double exp_upper_arg = v * v + ve * ve - 2. * v * ve * u_min;
	double exp_upper = exp(-exp_upper_arg / (v0 * v0));

	return 2. * M_PI * pow(v, i + 1.0) * (v0 * v0 / (2. * v * ve)) * (exp_lower - exp_upper);
}

double shm_halo_numeric(double vmin, double vesc, double v0, double ve, int i){
	/*
	Inputs: vmin, vesc, v0, ve, i (velocity power).
	Output: normalized value of the halo integral (0 if vmin > vesc+ve).*/
	if (vmin > vesc + ve) return 0.;

	double result, error;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

	gsl_function F;
	struct shm_radial_params params = {vesc, v0, ve, i};
	F.function = &f_shm_radial;
	F.params = &params;

	gsl_integration_qags(&F, vmin, vesc + ve, 1e-10, 1e-10, 1000, w, &result, &error);
	gsl_integration_workspace_free(w);

	return result/norm_shm_num(vesc, v0, 0.0);
}

/*########################################################################################
SHM + LMC BOOSTED COMPONENT (SHM_wLMC)
Mixture f = (1-w) f_SHM + w f_LMC, Eq. (1) of arXiv:2609.04175:
  f_LMC ~ exp(-|v-vb|^2/sigb^2) * Theta(vcut - |v-vb|), boosted by vb in the Galactic frame.
Same conventions as shm_halo(): ve = lab speed, i = velocity power (i=0 gives eta(vmin)).
w is a FRACTION in [0,1] (0.6% -> w = 0.006). cosb = cos(angle between vb and v_lab).
########################################################################################*/

double lmc_lab_bulk(double vb, double cosb, double ve){
	/*Speed of the LMC bulk velocity in the lab frame, |vb - v_lab|.
	Inputs: vb (Galactic-frame bulk speed), cosb (cos of angle between vb and v_lab), ve (lab speed).
	Output: lab-frame bulk speed.*/
	return sqrt(vb * vb + ve * ve - 2. * vb * ve * cosb);
}

double norm_lmc(double sigb, double vcut){
	/*Analytic normalization of the truncated Gaussian LMC component (int d^3v f = 1).
	Inputs: sigb (Gaussian width), vcut (hard cutoff around vb).
	Output: normalization constant.*/
	double x = vcut / sigb;
	double sigb3 = sigb * sigb * sigb;
	return pow(M_PI, 1.5) * sigb3 * (gsl_sf_erf(x) - 2. * x / sqrt(M_PI) * exp(-x * x));
}

struct lmc_radial_params { double uc; double sigb; double vcut; int i; };

double f_lmc_radial(double u, void *p){
	/*Radial integrand of the LMC halo integral after analytic angular integration
	(same structure as f_shm_radial). Inputs: u (lab speed), p (lmc_radial_params).
	Output: integrand value at u.*/
	struct lmc_radial_params *par = (struct lmc_radial_params *)p;
	double uc = par->uc, sigb = par->sigb, vcut = par->vcut;
	int i = par->i;

	if (u < 1e-6) return 0.0;
	/* isotropic limit (uc -> 0), avoids cancellation in the general formula */
	if (uc < 1e-3) return (u < vcut) ? 4. * M_PI * pow(u, i + 1.0) * exp(-u * u / (sigb * sigb)) : 0.;

	double c_min_raw = (u * u + uc * uc - vcut * vcut) / (2. * u * uc);
	double c_min = (c_min_raw > -1.) ? c_min_raw : -1.;
	double e_low = exp(-(u - uc) * (u - uc) / (sigb * sigb));
	double e_up  = exp(-(u * u + uc * uc - 2. * u * uc * c_min) / (sigb * sigb));

	return 2. * M_PI * pow(u, i + 1.0) * (sigb * sigb / (2. * u * uc)) * (e_low - e_up);
}

double lmc_halo(double vmin, double vb, double cosb, double sigb, double vcut, double ve, int i){
	/*Un-normalized halo integral of the LMC component (divide by norm_lmc).
	Inputs: vmin, vb, cosb, sigb, vcut, ve, i (velocity power).
	Output: value of the halo integral (0 if vmin is above the lab-frame endpoint).*/
	double uc = lmc_lab_bulk(vb, cosb, ve);
	double u_lo = MAX(MAX(vmin, 0.), uc - vcut);
	double u_hi = uc + vcut;
	if (u_lo >= u_hi) return 0.;

	double result, error;
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
	gsl_function F;
	struct lmc_radial_params params = {uc, sigb, vcut, i};
	F.function = &f_lmc_radial;
	F.params = &params;
	gsl_integration_qags(&F, u_lo, u_hi, 1e-10, 1e-10, 1000, w, &result, &error);
	gsl_integration_workspace_free(w);

	return result;
}

double shm_wlmc_halo(double vmin, double vesc, double v0, double ve, double w,
                     double vb, double cosb, double sigb, double vcut, int i){
	/*Normalized SHM + LMC mixture halo integral (hard-cutoff SHM, beta = 0).
	Inputs: vmin, vesc, v0, ve (lab speed), w (LMC fraction in [0,1]),
	vb, cosb, sigb, vcut (LMC parameters), i (velocity power).
	Output: (1-w) * SHM + w * LMC, each normalized to unity.
	Paper fiducial: vesc=544, v0=238, vb=570, cosb=-0.71, sigb=100, vcut=200.*/
	if (w < 0. || w > 1.) { printf("shm_wlmc_halo: w must be in [0,1]\n"); return 0.; }
	double h_shm = shm_halo(vmin, vesc, v0, ve, i);
	if (w == 0.) return h_shm;
	double h_lmc = lmc_halo(vmin, vb, cosb, sigb, vcut, ve, i) / norm_lmc(sigb, vcut);
	return (1. - w) * h_shm + w * h_lmc;
}

/*########################################################################################
N-BODY EXTRACTED HALO (LISANTI-TYPE)
Extracted phase-space densities from N-body simulations. See: 1010.4300 for more info.
In the following the integer "i" will control the power in velocity included in
calculation. For i=0 (usual calculation); i=1 (operator proportional to v);
i=2 (operator proportional v^2) and so on.
Use equation 15 in 1208.6426.
########################################################################################*/

struct Nk_params { double vesc_Nk; double v0_Nk; double k_Nk; };

double Nk_intd(double v, void *p){
	/*Integrand for the normalization of the Lisanti-type (N-body-extracted)
	velocity distribution.
	Inputs: v (velocity), p (Nk_params).
	Output: integrand value at v.*/
	struct Nk_params *params = (struct Nk_params *)p;
	double vesc_Nk = (params->vesc_Nk);
	double v0_Nk   = (params->v0_Nk);
	double k_Nk    = (params->k_Nk);
	double x = pow(v / v0_Nk, 1);
	double k_threshold = 0.1; // threshold for k_Nk to avoid division by zero or negative powers
	if (k_Nk >= k_threshold)
		return 4 * M_PI * v * v * pow(exp(-(v * v) / (k_Nk * v0_Nk * v0_Nk))
		       - exp(-(vesc_Nk * vesc_Nk) / (v0_Nk * v0_Nk * k_Nk)), k_Nk);
	else
		return v * v * exp(-x);
}

double Nk_uncert(double vesc, double v0, double k){
	/*Numerical normalization of the Lisanti-type velocity distribution via GSL QAGS.
	Inputs: vesc, v0, k (shape parameter).
	Output: normalization constant.*/
	double result, error;

	gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

	gsl_function F;
	struct Nk_params params = {vesc, v0, k};

	F.function = &Nk_intd;
	F.params = &params;

	gsl_integration_qags(&F, 0, vesc, 1e-6, 1e-6, 1000, w, &result, &error);
	// Fast integration --> see GSL documentation for more info
	// size_t neval; gsl_integration_qng (&F, 0, vesc, 1e-1, 1e-1, &result, &error, &neval);
	gsl_integration_workspace_free(w);

	return result;
}

struct lisanti_params { double vesc_lis; double v0_lis; double ve_lis; double k_lis; int i_lis; };

int fint(unsigned ndim, const double *k, void *p, unsigned fdim, double *fval){
	/*Cubature integrand for the Lisanti-type halo integral (angle-velocity space).
	Inputs: ndim, k (integration variables: k[0]=cos(theta), k[1]=v), p (lisanti_params), fdim.
	Output: fval[0] set to the integrand value (0 if the bracket is non-positive); return 0 on success.*/
	struct lisanti_params *params = (struct lisanti_params *)p;
	double vesc_lis = (params->vesc_lis);
	double v0_lis   = (params->v0_lis);
	double ve_lis   = (params->ve_lis);
	double k_lis    = (params->k_lis);
	int i_lis       = (params->i_lis);

	double A = exp(-(k[1] * k[1] + ve_lis * ve_lis + 2. * k[1] * ve_lis * k[0]) / (k_lis * v0_lis * v0_lis));
	double B = exp(-(vesc_lis * vesc_lis) / (k_lis * v0_lis * v0_lis));
	if (A - B > 0.0) {
		fval[0] = 2 * M_PI * k[1] * pow(k[1], (i_lis * 1.0)) * pow(A - B, k_lis);
	} else {
		fval[0] = 0;
	}

	return 0;
}

double lisanti_halo(double vmin, double vesc, double v0, double ve, double k, int i){
	/*Lisanti-type halo integral (integration of the uncertainties distribution)
	via 2D cubature over angle and velocity.
	Normalization is done by Nk_uncert().
	Inputs: vmin, vesc, v0, ve, k (shape parameter), i (velocity power).
	Output: value of the halo integral.*/
	if (vmin > vesc + ve) return 0.;
	
	double val, err;
	struct lisanti_params params = {vesc, v0, ve, k, i};

	// limits of integration
	double xl[2] = {-1, vmin};
	double xu[2] = {1, vesc + ve};

	hcubature(1, &fint, &params,
	          2, xl, xu,
	          0, 0, 1e-6, ERROR_INDIVIDUAL, &val, &err);
	return val/Nk_uncert(vesc, v0, k);;
}

/*########################################################################################
VELOCITY BINNING AND TABULATION
########################################################################################*/

int velocity(double vesc, double ve){
	/*Fills the global vel[] array with 'length' velocity bins from 0 to vesc+ve,
	with increasing point density toward the tail (high-v end).
	Inputs: vesc, ve.
	Output: always 0 (updates the global vel[] array).*/
	int i;
	double vmax = vesc + ve;
	double p = 2.0; /* p>1: higher p = more density near vmax. p=1 recovers linear. */
	double t;
	for (i = 0; i < length; i++) {
		t = (double)i / (length - 1);
		vel[i] = vmax * (1.0 - pow(1.0 - t, p));
	}
	return 0;
}

int fv_v(char *profile, double vesc, double v0, double beta, double ve, double k, int i,
         double w, double vb, double cosb, double sigb, double vcut){
	/*Fills the global fv[][] array (normalized halo integral values, f(v)*v^(i+1) bins)
	for every velocity bin and every velocity power up to i, for the selected profile.
	Inputs: profile (see valid_profiles array), vesc, v0, beta, ve, k, i (max velocity power),
	w, vb, cosb, sigb, vcut.
	Output: always 0 (updates the global vel[] and fv[][] arrays).*/

	velocity(vesc, ve);
	/* LMC tail extends beyond vesc+ve: extend the grid so halo() does not truncate it */
	if (strncmp(profile, "SHM_wLMC", 10) == 0 && w > 0.) {
		double vmax = MAX(vesc + ve, lmc_lab_bulk(vb, cosb, ve) + vcut);
		int jj;
		for (jj = 0; jj < length; jj++) vel[jj] = vmax * jj / length;
	}
	int m;
	int j;
	for (m = 0; m < i + 1; m++) {
		for (j = 0; j < length; j++) {
				fv[m][j] = halo_f(profile, vel[j], vesc, v0, beta, ve, k, m, w, vb, cosb, sigb, vcut);
			}
		}
	return 0;
}

/*########################################################################################
INTERPOLATION
########################################################################################*/

double halo(double vmin, int i){
	/*Cubic-spline interpolation of the precomputed halo integral fv[i][] over vel[],
	evaluated at vmin.
	Inputs: vmin, i (velocity power index into fv).
	Output: interpolated halo value, or 0 if vmin is outside the tabulated range.*/
	int m;
	double x[length], y[length];
	double result;

	for (m = 0; m < length; m++) {
		x[m] = vel[m];
		y[m] = fv[i][m];
	}

	if (vmin >= x[0] && vmin <= x[length - 1]) {
		gsl_interp_accel *fofv = gsl_interp_accel_alloc();
		gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, length);
		gsl_spline_init(spline, x, y, length);
		result = gsl_spline_eval(spline, vmin, fofv);
		gsl_spline_free(spline);
		gsl_interp_accel_free(fofv);

		return result;
	}
	else { // f(v) data out of the range
		return 0.;
	}
}

double halo_w(double vmin, gsl_interp_accel *ga, gsl_spline *gs){
	/*Evaluates an already-built GSL spline/accelerator pair at vmin (range-checked
	against the global vel[] array). Useful when the spline is built once and reused.
	Inputs: vmin, ga (GSL interpolation accelerator), gs (GSL spline).
	Output: interpolated value, or 0 if vmin is outside the tabulated range.*/
	if (vmin >= vel[0] && vmin <= vel[length - 1]) {
		return gsl_spline_eval(gs, vmin, ga);
	}
	else { // f(v) data out of the range
		return 0.;
	}
}

double halo_f(char *profile, double vmin, double vesc, double v0, double beta, double ve, double k, int i,
    		  double w, double vb, double cosb, double sigb, double vcut){
	/*Direct (non-interpolated, non-tabulated) evaluation of the normalized halo
	integral for the selected profile at a single vmin.
	Inputs: profile (see valid_profiles array), vmin, vesc, v0, beta, ve, k, i (velocity power).
	Output: normalized halo integral value, or 0.0 if the profile is not implemented.*/
	if (strncmp(profile, "SHM", 10) == 0) {
		return shm_halo(vmin, vesc, v0, ve, i);
	}
	else if (strncmp(profile, "SHM_beta", 10) == 0) {
		return shm_halo_beta(vmin, vesc, v0, ve, beta, i);
	}
	else if (strncmp(profile, "SHM_numeric", 10) == 0) {
		return shm_halo_numeric(vmin, vesc, v0, ve, i);
	}
	else if (strncmp(profile, "SHM_wLMC", 10) == 0) {
		return shm_wlmc_halo(vmin, vesc, v0, ve, w, vb, cosb, sigb, vcut, i);
	}
	if (strncmp(profile, "Lisanti", 10) == 0) {
		return lisanti_halo(vmin, vesc, v0, ve, k, i);
	}
	else {
		printf("The selected profile is not implemented.\n");
		return 0.0;
	}
}

/*########################################################################################
TABLE I/O
Modified functions for writing and reading halo tables to/from disk. (30-07-2015)
########################################################################################*/

int write_fv_v(char *path, char *profile, double vesc, double v0, double beta, double ve, double k, int i,
			   double w, double vb, double cosb, double sigb, double vcut){
	/*Computes fv[][] for the selected profile (via fv_v) and writes the velocity
	bins and halo values to a text table file.
	Inputs: path (output file path), profile, vesc, v0, beta, ve, k, i (max velocity power),
	w, vb, cosb, sigb, vcut.
	Output: always 0 (writes the table to disk).*/
	fv_v(profile, vesc, v0, beta, ve, k, i, w, vb, cosb, sigb, vcut);
	FILE *table;
	table = fopen(path, "w+");
	fprintf(table, "%d %d \r\n", length, i);
	int j, m;
	for (j = 0; j < length; j++) {
		fprintf(table, "%.5E", vel[j]);
		for (m = 0; m <= i; m++) {
			fprintf(table, " %.5E", fv[m][j]);
		}
		fprintf(table, " \r\n");
	}
	fclose(table);

	return 0;
}

void define_and_write_halo(char *path, char *profile, double vesc, double v0, double beta, double ve, double k, int i, 
						   double w, double vb, double cosb, double sigb, double vcut){
	/*Convenience wrapper that computes the normalization and writes the halo table
	to file for the selected profile, with basic input validation and logging.
	Inputs: path (output file path), profile, vesc, v0, beta, ve, k, i (max velocity power).
	Output: none (writes the table to disk).*/
	if (i > power - 1) {
		printf("Please select a power in velocity up to %d or change the definition of power in source/halo.c\n", power - 1);
	}
	else {
		printf("Calculating halo integrals for %s up to order %d in v...\n", profile, i);
		write_fv_v(path, profile, vesc, v0, beta, ve, k, i, w, vb, cosb, sigb, vcut);
		printf("Done!\n");
	}
}

void read_halo(char *path){
	/*Reads a precomputed halo table (as written by write_fv_v) from file into the
	global vel[] and fv[][] arrays. Supports i = 0, 1, or 2 velocity powers.
	Inputs: path (input file path).
	Output: none (updates the global vel[] and fv[][] arrays).*/
	FILE *table;
	int i;
	int j;
	int check_length;

	table = fopen(path, "r");

	fscanf(table, "%d %d", &check_length, &i);

	for (j = 0; j < length; j++) {
		if (i == 0) fscanf(table, "%E %E", &vel[j], &fv[0][j]);
		if (i == 1) fscanf(table, "%E %E %E", &vel[j], &fv[0][j], &fv[1][j]);
		if (i == 2) fscanf(table, "%E %E %E %E", &vel[j], &fv[0][j], &fv[1][j], &fv[2][j]);
	}

	fclose(table);
}

/*########################################################################################
ANNUAL MODULATION / TIME-DEPENDENT LAB-FRAME VELOCITY
########################################################################################*/

double lab_frame_speed(double v0_lsr, const double v_pec[3], double v_earth_avg, double t0, int t){
	/*Computes the lab-frame speed |v_lab| = |v0_lsr_vec + v_pec + v_earth(t)| for a
	given day. Reference: Eq. 10-11 of arXiv:2105.00599.
	Inputs: v0_lsr (local standard of rest speed, phi component only, km/s; Table 1: 238),
	v_pec[3] (solar peculiar velocity vector: vr, vphi, heta, km/s; Table 1: 11.1, 12.2, 7.3),
	v_earth_avg (<|v_earth|>, km/s; Table 1: 29.8),
	t0 (reference day offset from March 22 2018, i.e. the day corresponding to t=0),
	t (day index).
	Output: lab-frame speed at day t.*/
	double delta_t = t - t0;
	double v_earth[3];
	earth_velocity_vector(delta_t, v_earth_avg, v_earth);

	double vx = 0.0    + v_pec[0] + v_earth[0];
	double vy = v0_lsr + v_pec[1] + v_earth[1];
	double vz = 0.0    + v_pec[2] + v_earth[2];

	return sqrt(vx * vx + vy * vy + vz * vz);
}

double lab_frame_speed_annual_avg(double v0_lsr, const double v_pec[3]){
	/*Annual-average lab-frame speed, using the fixed Earth-velocity vector evaluated
	at March 9 (recommended for analyses not targeting annual modulation).
	Reference: Eq. 12 of arXiv:2105.00599.
	Inputs: v0_lsr, v_pec[3] (see lab_frame_speed).
	Output: annual-average lab-frame speed.*/
	double v_earth_mar9[3] = {29.2, -0.1, 5.9}; // Eq. 12

	double vx = 0.0    + v_pec[0] + v_earth_mar9[0];
	double vy = v0_lsr + v_pec[1] + v_earth_mar9[1];
	double vz = 0.0    + v_pec[2] + v_earth_mar9[2];

	return sqrt(vx * vx + vy * vy + vz * vz);
}

void earth_velocity_vector(double delta_t, double v_earth_avg, double v_out[3]){
	/*Computes the vector Earth velocity relative to the Sun in the galactic frame,
	components (vr, vphi, heta): r points radially inward, phi points in the
	direction of the Milky Way's rotation. Reference: Eq. 11 of arXiv:2105.00599.
	Inputs: delta_t (days since March 22, 2018, arbitrary reference date),
	v_earth_avg (<|v_earth|> = 29.8 km/s, Table 1).
	Output: v_out[3] filled with the Earth velocity vector.*/
	double omega = 0.0172; // rad/day
	v_out[0] = v_earth_avg * (0.9941 * cos(omega * delta_t) - 0.0504 * sin(omega * delta_t));
	v_out[1] = v_earth_avg * (0.1088 * cos(omega * delta_t) + 0.4946 * sin(omega * delta_t));
	v_out[2] = v_earth_avg * (0.0042 * cos(omega * delta_t) - 0.8677 * sin(omega * delta_t));
}

void define_and_write_halo_time(char *profile, double vesc, double v0, double beta,
                                 double v0_lsr, const double v_pec[3],
                                 double k, int i, double t0, double T,
								 double w, double vb, double cosb, double sigb, double vcut){
	/*Computes and writes halo tables for a time series of days T, using the
	time-dependent lab-frame Earth velocity (for annual-modulation studies).
	Inputs: profile, vesc, v0, beta, v0_lsr, v_pec[3], k, i (max velocity power),
	t0 (reference day offset), T (number of days to compute), w, vb, cosb, sigb, vcut.
	Output: none (writes one table file per day under halo_table/).*/
	if (i > power - 1) {
		printf("Please select a power in velocity up to %d or change the definition of power in source/halo.c\n", power - 1);
	}
	else {
		printf("Calculating halo integrals for %s up to order %d in v...\n", profile, i);

		clock_t start, end;
		double cpu_time_used;
		start = clock();

		int t;
		for (t = 0; t < T; t++) {
			char path[32];
			snprintf(path, sizeof(char) * 32, "halo_table/halo_table_%i.dat", t);
			printf("Calculating for t=%i...\n", t);
			double ve_t = lab_frame_speed(v0_lsr, v_pec, 29.8, t0, t); // v_earth_avg = 29.8 km/s, Table 1
			write_fv_v(path, profile, vesc, v0, beta, ve_t, k, i, w, vb, cosb, sigb, vcut);
		}

		printf("Done!\n");
		end = clock();
		cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
		printf("Elapsed time: %f seconds\n", cpu_time_used);
	}
}