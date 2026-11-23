#pragma once

static size_t MEGA=0;
double fsi(int A, double Er);
double f(double Er, int A, int Z, double a0, double a1, double a2, double a3, double a4, 
         double a5, double a6, double a7, double a8, double a9, double a10);
double f_Men(double Er, int A, int Z, double a0, double a1, double a2, double a3, double a4,
			 double a5);

int read_nuclear_responses(const char *dir);
int nuclear_responses_loaded(void);
double FF_iso_GCN5082_MEN(int A, int Z, char* Symbol1, char* Symbol2, char* tau1, char* tau2, double Er);
double FF(char* NR_framework, int A, int Z, char * Symbol1, char * Symbol2, char * N1, char * N2, double Er);

double FormFact_v0(char* NR_framework,int A,int Z,int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);
double FormFact_v2(char* NR_framework,int A,int Z, int i, int j, char * N1, char * N2, double Er, double mchi, double jchi);