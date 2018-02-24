#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define NX 			200 					//number of segments along logx-axis, number of points=NX+1
#define MAXITER     10 						//number of iterations, 10 should be fine
	
#define mp 			1.672621777e-24 		// proton mass in grams
#define me 			9.10938356e-28		    // electron mass in grams				
#define Kb 		    1.3806485e-16 			// Boltzmann constant in erg*K^{-1}
#define qel 		4.80320425e-10 			// electron charge in cgs statcoulombs
#define re 			2.8179403267e-13		// classical electron radius in cm
#define parsec	    3.0856775807e18			// 1 parsec in cm
#define clight 	    299792.458 				// speed of light in km/s
#define gamma_th    5.0/3.0 				// gass adiabatic constant
#define gamma_cr    4.0/3.0 				// cosmic ray adiabatic constant

double SynchF(double x);

void FirstSegment(double y1, double y2, double x0, double x1, double x2, double *I);

void OtherSegment(double y1, double y2, double x1, double x2, double *I);

void PowerSegment(double y1, double y2, double x1, double x2, double *I);

double FindRprecLimit(double Ma, double MS, double zeta);

double Luminosity(double freq, double Bn, double R, double *f0, double *X, double limit, double Chi_esc, double Vshock, double rho0);

void nldsa(double MS, double Radius, double rho0, double T0, double B, double psi, double zeta, double Chi_esc, double *gamma_eff,double *Rtot_out, double *p_max, double *e_max, double *B2, double *izlaz_L, double *alfa);