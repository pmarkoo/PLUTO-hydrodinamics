#include "blasi.h"


// aproksimacija iz rada Fouka & Ouichaoui, RAA, 2013
// maximalna relativna greska +/- 0.26 %
double SynchF(double x)
{
	double a1_1 = -0.97947838884478688;
	double a1_2 = -0.83333239129525072;
	double a1_3 = +0.15541796026816246;

	double a2_1 = -4.69247165562628882e-2;
	double a2_2 = -0.70055018056462881;
	double a2_3 = +1.03876297841949544e-2;

	double F1 = 2.1495282415344786367;
	double F2 = sqrt(M_PI*0.5);

	double delta1, H1, delta2, H2, F;

	H1 = a1_1 * x + a1_2 * sqrt(x) + a1_3 * pow(x, 1./3.);
	delta1 = exp(H1);

	H2 = a2_1 * x + a2_2 * sqrt(x) + a2_3 * pow(x, 1./3.);
	delta2 = 1.0 - exp(H2);

	F = F1 * pow(x, 1./3.) * delta1 + F2 * exp(-x) * sqrt(x) * delta2;

	return F;
}


void FirstSegment(double y1, double y2, double x0, double x1, double x2, double *I)
{
	double a, b;
	b = log(y1/y2)/log(x2/x1);
    a = y1 * pow(x1,b);
    *I = a/pow(0.5*(x0+x1),b);
}


void OtherSegment(double y1, double y2, double x1, double x2, double *I)
{
	double a, b;
	b = log(y1/y2)/log(x2/x1);
    a = y1 * pow(x1,b);
    *I = a/pow(0.5*(x1+x2),b);
}


void PowerSegment(double y1, double y2, double x1, double x2, double *I)
{
	double a, b, F;
	double EPS = 1.e-6;

	/* if (y2 == 0) {
        fprintf(stderr, "TwoSegments()::y2 should be != 0! Aborting...\n");
        exit(EXIT_FAILURE); 
    } 
    else if (x1 == 0) {
    	fprintf(stderr, "TwoSegments()::x1 should be != 0! Aborting...\n");
        exit(EXIT_FAILURE); 
     } */

	// calculate parameters in parabola f(x) = a/x^b 
    b = log(y1/y2)/log(x2/x1);
    a = y1 * pow(x1,b);


    if (fabs(b-1.0) > EPS) {
    	*I = a * (pow(x2,1-b) - pow(x1,1-b))/(1-b);
	} else {
    	*I = a * log(x2/x1);
    }
}


//funkcija koja racuna luminoznost iz spektra elektrona, fe = fp/Kep
double Luminosity(double freq, double Bn, double R, double *f0, double *X, double limit, double Chi_esc, double Vshock, double rho0)
{
	double mp_me = 1836.15267; //proton over electron mass: mp/me
	double GHz = 1.e+9;
	int i;
	double Kep = 0.01;
	double re_factor = 0.2; //renormalization factor, Berezhko and Volk 2004
	double k_CD = 0.9; //position of CD normalized to shock radius Rcd = k_CD * Rshock
	int electron_limit;
	double Const = 3 * qel * Bn * mp_me * mp_me * mp_me / (4 * M_PI * mp * clight * 1.e5);
	double lambda = (3*Chi_esc/0.05) * (Vshock/(Bn/sqrt(4*M_PI * rho0 * mp)/1.e+5)) * (Vshock/clight) * R; //jedinice = [pc]
	double FL[NX+1]; 
	double nu_krit;
	double x, rez, sum=0.0, L, L1;

	nu_krit = X[0] * X[0] * Const;
	x = freq * GHz / nu_krit;	
	// Zirakashvili and Aharonian (2007)
	FL[0] = X[0] * X[0] * Kep * f0[0] * SynchF(x) * (1 + 0.523*pow(X[0]/limit, 2.25))*(1 + 0.523*pow(X[0]/limit, 2.25)) * exp(-X[0] * X[0]/(limit*limit));

	for (i=1; i <= NX; i++) 
	{ 
		nu_krit = X[i] * X[i] * Const;
		x = freq * GHz / nu_krit;
		// Zirakashvili and Aharonian (2007)
		FL[i] = X[i] * X[i] * Kep * f0[i] * SynchF(x) * (1 + 0.523*pow(X[i]/limit, 2.25))*(1 + 0.523*pow(X[i]/limit, 2.25)) * exp(-X[i] * X[i]/(limit*limit));
		rez = 0.5 * (FL[i-1] + FL[i]) * (X[i] - X[i-1]); // obicno trapezno pravilo
		sum += rez;
	}

	//L = 2.94193e-21 * re_factor * (4.0/3.0) * M_PI * (R * parsec) * (R * parsec) * (R * parsec) * Bn * sum;
	// * parsec * parsec * parsec
	L = 3.70206e-20 * re_factor  * Bn * sum * lambda * (R*R - 2*R*lambda + 2*lambda*lambda - 
		 exp((-R + k_CD*R)/lambda)*(k_CD*k_CD*R*R - 2*k_CD*R*lambda + 2*lambda*lambda)) * parsec * parsec * parsec;

	//L = sum;
	return L;
}


double FindRprecLimit(double Ma, double MS, double zeta)
{
	//solution of simple non-standard equation by 2-step iteration
	double x0, F, G, R; //R = Rprec
	double LeftR, RightR;
	x0 = sqrt(1 + 2*Ma/(1-zeta));
	x0 = pow(sqrt(x0) + 2*Ma/(1-zeta), 0.4);
	LeftR = 1; 
	RightR = x0;
	//dobili smo asimptotu, sledi dodatna popravka
	// resavamo slozenu nelinearnu jednacinu binarnom pretragom
	// da bi konacno dobili Ms,1min iz cega sledi Rprec,max
	while ((RightR-LeftR)>0.01)
	{
		R = 0.5*(RightR+LeftR);
		F = MS*pow(R,(-gamma_th-1)*0.5)/sqrt( 1 + zeta*(gamma_th-1)*MS*MS*(1-pow(R,-gamma_th))/Ma);
		G = 2/(2 - (1-zeta)*(pow(R,2.5) - sqrt(R))/Ma);
		if (F > G)
		{
			LeftR = R;
		} else {
			RightR = R;
		}
	}
	return R;
}


void nldsa(double MS, double Radius, double rho0, double T0, double B, double psi, double zeta, double Chi_esc, double *gamma_eff_out, double *Rtot_out, double *p_max, double *e_max, double *B2, double *izlaz_L, double *alfa)
{
	FILE *fp0, *fp1, *fp2;
	int i, j, k, iter;
	double a,b,c, MIN_P, MIN_Pe, MIN_e; //coefficient for quadratic equation
	//double B = 5.e-6; // Magnetic field strength in Gauss (5 microG), typical in ISM
	//double T0 = 1.e+4; // ISM temperature in K
	//double rho0 = 0.1;  //ISM density rho_0/proton_mass in cm^{-3}
	//double psi = 3.5; //injection parameter pinj = psi * pth
	//double zeta = 0.4; //Caprioli parameter
	//double Chi_esc = 0.1; //free escape boundary parameter
	double MS12, B1; //compression imediately upstream (squared) MS12 = Ms,1^2
	//double MS = 2000; //Shock Mach number about Cs = 10 km/s 
	double f0[NX+1], U[NX+1], F[NX+1], X[NX+1], P_CR[NX+1], P[NX+1];
	double f0s[NX]; //sredine f0[p/2]
	double sum, izraz, h, rez, U1, U2, rez1, rez2, rezU, m1, m2;
	double expo=(-1.0) * (gamma_th + 1.0);
	double Rsub, Rprec, Rtot,  Rsub1, Rsub2;
	double T1, T2;
	double k1, k2, k3, k4, dX;
	double gamma_eff;
	double ngas1;
	//double D0 = 3.13404e+16/B; //diffusion coefficient D(p) = D0*p^2/(1+p^2) - Bohm diffusion, kao Ferrand (2010)
	double Rsnr; // SNR radius in cm
	double eff = 4.e-3; // default eta = efficiency of injection if constant
    double Vshock = MS * sqrt(gamma_th * Kb * T0/mp)/1.e+5; // shock velocity in km/s
    double MA = Vshock/(B/sqrt(4*M_PI * rho0 * mp)/1.e+5); //Alfvenic Mach number
    double pth, pinj; // p_{th} thermal momentum in mp*c, p_{inj} injection momentum
    double XBEG, XEND, XEND1; // lower and upper limit for momentum, log scale
    double Alfven_damping; // term for heating due to the damping of Alfven waves in the precursor, Eq. 24 (Ferrand memo 1)
    double Rprec_MAX; //maksimalno Rsub, za Rprec=1
    double Pw1;
    double Rprec_MIN=0.001; //avoid Rprec=1 which gives nan values
    double LeftR, RightR; // leva i desna granica pretrage Rprec, primenjena binarna pretraga slozenosti logN

	Rsnr = Radius * parsec;
    // 1.0*qel in next equation is proton charge
    // Bohm diffusion 
    //XEND = log10(Chi_esc*Vshock*1.e+5*Rsnr/(mp*clight*clight*clight*1.e+15/3/qel/B)); //D0 = mc^3/(3*e*B), assumed Bohm diffusion

  	//XEND1 = log10(3.0 * (Vshock/(clight*clight*1.e+5)) * 1.0 * qel * B * Rsnr/(8.0 * mp * clight * 1.e+5)); //Bell 2013 i Boki
    // ovo daje jako slicno Bohm-ovoj difuziji, razlikuju se jer Bohm uzima Chi_esc=1/10, Boki tu uzima 1/8

   // printf("Vshock = %f  Valfen = %f  MA = %f \n", Vshock, B/sqrt(4*M_PI * rho0 * mp)/1.e+5, MA);
	//fp0 = fopen("U.txt","w");
	//fp1 = fopen("f.txt", "w");
	//fp2 = fopen("rh.txt","w");

	LeftR = 1 + Rprec_MIN;
	
	if (zeta < 1.0)
	{	
		RightR = FindRprecLimit(MA, MS, zeta);
	} else {
		RightR = pow(MS*MS/(1 + (gamma_th-1)*MS*MS*(1-pow(MS, -2*gamma_th/(gamma_th+1)))/MA), 1/(gamma_th+1));
		RightR = pow(MS*MS/(1 + (gamma_th-1)*MS*MS*(1-pow(RightR, -gamma_th))/MA), 1/(gamma_th+1));
	}

	//printf("Vshock = %f Valfen= %f MA = %f  MAlimit = %f \n", Vshock, (B/sqrt(4*M_PI * rho0 * mp)/1.e+5), MA, RightR);

	// Bisection method for root-finding U(pmax)=0
	// method is guaranteed to converge to a root of f if f is a continuous function
	// on the interval [a, b] and f(a) and f(b) have opposite signs
    while ((RightR - LeftR) > 0.004)
    //while (Rprec < RightR)
	{
	//Koraci++;
	Rprec = 0.5*(LeftR + RightR);

	//koristimo Rtot = Rsub * Rprec;
	//Blasi et al. (2005), Eq. 15
	//Rtot = pow(MS, 2/(gamma_th+1)) * pow(((gamma_th+1) * pow(Rsub, gamma_th) - (gamma_th-1) * pow(Rsub, gamma_th+1))*0.5, 1/(gamma_th+1));
	MS12 = MS * MS * pow(Rprec, -gamma_th-1)/(1 + zeta*(gamma_th-1)*MS*MS*(1 - pow(Rprec, -gamma_th))/MA);
	Pw1 = (1 - zeta) * pow(Rprec, 2.5) * (1 - 1/(Rprec*Rprec))/(4*MA);

	//sada ide resenje kvadratne koje daje Rsub ...... kao pozitivno resenje
	if (zeta < 1) {
		a = 2*(gamma_th-2)*MS12*Pw1;
		b = -(2 + (gamma_th-1 + 2*gamma_th*Pw1)*MS12);
		c = (gamma_th+1)*MS12;
		Rsub = (-b - sqrt(b*b-4*a*c))/(2*a);
	} else {
		Rsub = (gamma_th+1)*MS12/((gamma_th-1)*MS12 + 2);
	}

	//eficiency calculation (optional), commenting it leaves eff = const.
	//occurence of multiple solutions is dramatically reduced by using this recipe
	//Blasi et al. (2005), Eq. 27
	eff = 4.0 * (Rsub - 1.0) * psi * psi * psi * exp(-psi*psi)/(3.0 * sqrt(M_PI));

	//Rsub = (gamma_th + 1)*MS12*MS12/((gamma_th - 1)*MS12*MS12 + 2);
	Rtot = Rsub * Rprec;
	gamma_eff = (MS*MS*(Rtot+1) - 2*Rtot)/(MS*MS*(Rtot-1));

	ngas1 = rho0 * Rtot / Rsub;

	//B1 = B * sqrt(0.5*(1 - zeta)*MA*(1 - 1/(Rtot*Rtot))*pow(Rtot,1.5)  + 1.0)/sqrt(11);
	B1 = B * sqrt(0.5*(1 - zeta)*MA*(1 - 1/(Rprec*Rprec))*pow(Rprec,1.5)  + 1.0);
	XEND = log10(Chi_esc*Vshock*1.e+5*Rsnr/(mp*clight*clight*clight*1.e+15/3/qel/B1)); //D0 = mc^3/(3*e*B), assumed Bohm diffusion
	//XEND = 5.0;

	//downstream temperature Blasi memo 2, Eq. (65) + (66)

	T1 = T0 * pow(Rprec, gamma_th-1)*(1 + zeta*(gamma_th-1)*MS*MS*(1 - pow(Rprec, -gamma_th))/MA);
	T2 = T1 * ((gamma_th+1)*Rsub - (gamma_th-1)*(1-(Rsub-1)*(Rsub-1)*(Rsub-1)*gamma_th*Pw1*MS12)) / ((gamma_th + 1 - (gamma_th-1)*Rsub)*Rsub);
	//T2 = T0 * pow(Rprec, gamma_th-1) * (gamma_th + 1 - (gamma_th-1)/Rsub) / (gamma_th + 1 - (gamma_th-1)*Rsub);
	
	pth = sqrt(2 * mp * Kb * T2)/(mp * clight * 1.e+5); // p_{th} thermal momentum in mp*c
    pinj = psi * pth; // p_{inj} injection momentum
   // pinj =1.0e-2;
    XBEG = log10(pinj); // lower limit for momentum, log scale

   // printf("XBEG = %f  T1 = %lf    T2 = %lf \n", XBEG, T1, T2);

	h = (XEND-XBEG)/NX; //step in log(p) scale
	

	//pocetna inicijalizacija
	for (i=0; i <= NX; i++)
	{
		 f0[i] = 0.0;
		 U[i]  = 1/Rprec;
		 X[i]  = pow(10, XBEG + i*h);
	}


	for (iter=0; iter < MAXITER; iter++)
	{
		sum = 0.0;
		U[0] = 1/Rprec;

	    //racuna se sukcesivno f0[p], U[p]
	    F[0] = (3 * Rtot * U[0])/((Rtot * U[0] - 1) * X[0]);

	    //dobijeno kada se primeni limf0(p) kada p -> pinj
	    f0[0] = 3*Rsub/(Rsub-1)*(eff * ngas1/(4 * M_PI * pinj * pinj * pinj));

		for (i=0; i < NX; i++)
		{
			//pravi se integrand za   \int_{pinj}^{p}
			F[i+1] = (3 * Rtot * U[i+1])/((Rtot * U[i+1] - 1) * X[i+1]);
			//integracija provlacenjem parabole
			PowerSegment(F[i], F[i+1], X[i], X[i+1], &rez1);
			f0[i+1] = (eff * ngas1/(4 * M_PI * pinj * pinj * pinj)) * ((3*Rsub)/(Rtot*U[i+1] - 1)) * exp(-(sum+rez1));

			// racunanje medju-tacaka f0s(p+dp/2), neophodno za RK4 i RK2			
			//f0s[i] = pow(2*pinj/(X[i]+X[i+1]), 3*Rsub/(Rsub-1))*(eff * ngas1/(4 * M_PI * pinj * pinj * pinj))*3*Rsub/(Rsub-1);
			if (i>0) {
				if (i==1) {
					FirstSegment(f0[1], f0[2], X[0], X[1], X[2], &rez);
					f0s[0] = rez;
					OtherSegment(f0[1], f0[2], X[1], X[2], &rez);
					f0s[1] = rez;
				} else {
					OtherSegment(f0[i], f0[i+1], X[i], X[i+1], &rez);
					f0s[i] = rez;
				}
			}
			sum += rez1;
		}

		/*
		//cosmic ray pressure, normalized to rho_0*u_0^2, P^{CR,int}_CR/rho_0*u_0^2
		//!!! obratiti paznju na ovu normalizaciju kasnije, mozda nesto treba menjati
		P[NX] = gamma_th * MS * MS * (Rprec/ngas1) * pow(clight/Vshock, 2.0) * 4 * M_PI * pow(X[NX], 4)*f0[NX]/sqrt(1 + pow(X[NX], 2))/3;
		P_CR[NX] = 0.0; //limes je zaista = 0 kada p -> pmax

		for (i=NX-1; i >= 0; i--)
		{
			P[i] = gamma_th * MS * MS * (Rprec/ngas1) * pow(clight/Vshock, 2.0) * 4 * M_PI * pow(X[i], 4)*f0[i]/sqrt(1 + pow(X[i], 2))/3;
			PowerSegment(P[i], P[i+1], X[i], X[i+1], &rez1);
			P_CR[i] = P_CR[i+1] + rez1;
			//ovo cisto radi probe, delimo sa Pth,0
			//P_CR[i] = P_CR[i] * gamma_th * MS * MS;
		}
		*/

    	for (i=0; i < NX; i++)
    	{
    		dX = X[i+1]-X[i];
    		//Cauchy-jevi problemi dU/dp = f(U,p), U(0) = 1/Rprec


			// ---- Euler 1st order method for ODEs, solving U[p] O(h^2)----
			// ---- Vrednost U[i+1] ----		
			/*
			rez = 4 * M_PI * pow(X[i],4.0)/(3 * sqrt(1 + X[i]*X[i])) * (1/rho0) * pow(clight/Vshock, 2.0) * f0[i];
    		k1 = dX * rez/(1.0 - pow(U[i], expo)/pow(MS,2.0));
    		U[i+1] = U[i] + k1;
    		*/


			// !!!!!!!!!! ---- RUNGE-KUTTA 2th order method for ODEs, 2nd version solving U[p] (RK2) ----
    		// ovoj verziji ne trebaju h/2 tacke tj. niz f0s[p], preciznost ista kao za gornji RK2 O(h^3)
    		/*
    		rez = 4 * M_PI * pow(X[i],4.0)/(3 * sqrt(1 + X[i]*X[i])) * (1/rho0) * pow(clight/Vshock, 2.0) * f0[i];
    		k1 = dX * rez/(1.0 - pow(U[i], expo)/pow(MS,2.0));

    		rez = 4 * M_PI * pow(X[i+1],4.0)/(3 * sqrt(1 + X[i+1]*X[i+1])) * (1/rho0) * pow(clight/Vshock, 2.0) * f0[i+1];
    		k2 = dX * rez/(1.0 - pow(U[i]+k1, expo)/pow(MS,2.0));

    		U[i+1] = U[i] + 0.5*(k1 + k2);
    		*/


    		// ---- RUNGE-KUTTA 3rd order method for ODEs, solving U[p] (RK3) ----	
    		/*
    		rez = 4 * M_PI * pow(X[i],4.0)/(3 * sqrt(1 + X[i]*X[i])) * (1/rho0) * pow(clight/Vshock, 2.0) * f0[i];
    		k1 = dX * rez/(1.0 - pow(U[i], expo)/pow(MS,2.0));
    		rez = 4 * M_PI * pow(X[i]+dX/2,4.0)/(3 * sqrt(1 + (X[i]+dX/2)*(X[i]+dX/2))) * (1/rho0) * pow(clight/Vshock, 2.0) * f0s[i];
    		k2 = dX * rez/(1.0 - pow(U[i]+k1/2, exo)/pow(MS,2.0));
    		rez = 4 * M_PI * pow(X[i+1],4.0)/(3 * sqrt(1 + X[i+1]*X[i+1])) * (1/rho0) * pow(clight/Vshock, 2.0) * f0[i+1];
    		k3 = dX * rez/(1.0 - pow(U[i]-k1+2*k2, expo)/pow(MS,2.0));

    		U[i+1] = U[i] + (k1 + 4*k2 + k3)/6.0;
    		*/

    		
			// ---- RUNGE-KUTTA 4th order method for ODEs, solving U[p] (RK4) ----	
			
		Alfven_damping = 1.0 + zeta*(gamma_th-1.0)*MS*MS/MA;
    		rez = 4 * M_PI * pow(X[i],4.0)/(3 * sqrt(1 + X[i]*X[i])) * (1/rho0) * pow(clight/Vshock, 2.0) * f0[i];
    		k1 = dX * rez/(1.0 - pow(U[i], expo)*Alfven_damping/(MS*MS) - (1-zeta)*(U[i]*U[i]+3)/(8*MA*pow(U[i], 2.5)));

    		rez = 4 * M_PI * pow(X[i]+dX/2,4.0)/(3 * sqrt(1 + (X[i]+dX/2)*(X[i]+dX/2))) * (1/rho0) * pow(clight/Vshock, 2.0) * f0s[i];
    		k2 = dX * rez/(1.0 - pow(U[i]+k1/2, expo)*Alfven_damping/(MS*MS) - (1-zeta)*((U[i]+k1/2)*(U[i]+k1/2)+3)/(8*MA*pow(U[i]+k1/2, 2.5)));
  
    		k3 = dX * rez/(1.0 - pow(U[i]+k2/2, expo)*Alfven_damping/(MS*MS) - (1-zeta)*((U[i]+k2/2)*(U[i]+k2/2)+3)/(8*MA*pow(U[i]+k2/2, 2.5)));

    		rez = 4 * M_PI * pow(X[i+1],4.0)/(3 * sqrt(1 + X[i+1]*X[i+1])) * (1/rho0) * pow(clight/Vshock, 2.0) * f0[i+1];
    		k4 = dX * rez/(1.0 - pow(U[i]+k3, expo)*Alfven_damping/(MS*MS) - (1-zeta)*((U[i]+k3/2)*(U[i]+k3/2)+3)/(8*MA*pow(U[i]+k3/2, 2.5)));

    		U[i+1] = U[i] + k1/6 + k2/3 + k3/3 + k4/6;
			
    		//Runge-Kutta 4 bez Alfen wave heating-a
    		/*
    		rez = 4 * M_PI * pow(X[i],4.0)/(3 * sqrt(1 + X[i]*X[i])) * (1/rho0) * pow(clight/Vshock, 2.0) * f0[i];
    		k1 = dX * rez/(1.0 - pow(U[i], expo)/pow(MS,2.0));

    		rez = 4 * M_PI * pow(X[i]+dX/2,4.0)/(3 * sqrt(1 + (X[i]+dX/2)*(X[i]+dX/2))) * (1/rho0) * pow(clight/Vshock, 2.0) * f0s[i];
    		k2 = dX * rez/(1.0 - pow(U[i]+k1/2, expo)/pow(MS,2.0));
    		k3 = dX * rez/(1.0 - pow(U[i]+k2/2, expo)/pow(MS,2.0));

    		rez = 4 * M_PI * pow(X[i+1],4.0)/(3 * sqrt(1 + X[i+1]*X[i+1])) * (1/rho0) * pow(clight/Vshock, 2.0) * f0[i+1];
    		k4 = dX * rez/(1.0 - pow(U[i]+k3, expo)/pow(MS,2.0));

    		U[i+1] = U[i] + k1/6 + k2/3 + k3/3 + k4/6;
    		*/
			
    	} 
    	
    	

    }
	
	//fprintf(fp0, "%f  	 %f     %.9f     %.9f 	%.9f \n", log10(Rtot), U[NX], Rsub, Rprec, Rtot);

    // very simple implementation of binary search for Rprec leading to U(pmax)=1

    if (U[NX] < 1.0) {
    	RightR = Rprec;
    } 
    else 
    {
    	LeftR = Rprec;
    }


	//printf("Rsub = %.9f  Rprec = %f  Rtot = %f  U[pmax] = %f \n", Rsub, Rprec, Rtot, U[NX]);
	//fprintf(fp0, "%f  	 %f      %.9f 	%.9f \n", Rprec, U[NX], Rsub, Rtot);

	}
	//printf("Broj koraka = %d \n", Koraci);
	//printf("Rsub = %.9f  Rprec = %f  Rtot = %f   eff = %f  U[pmax] = %f  gamma_eff = %f \n", Rsub, Rprec, Rtot, eff, U[NX], gamma_eff);

    /*for (k=0; k <= NX; k++)	{
    	//fprintf(fp0, "%f   %f \n", log10(X[k]), log10(X[k]*X[k]*X[k]*X[k]*f0[k]));
    	fprintf(fp1, "%f   %f \n", log10(X[k]), 3*Rtot/(Rtot*U[k] - 1.0));
    } */
	///for (k=0; k <= NX; k++) fprintf(fp0, "%f  	 %f    %f\n", log10(X[k]), log10(f0[k]), log10(f0s[k]));
	//for (k=0; k <= NX; k++) fprintf(fp0, "%f  	 %f \n", log10(X[k]), log10(U[k]));
	//for (k=0; k <= NX; k++) fprintf(fp0, "%f  	 %f \n", X[k], U[k]);
	//for (k=0; k <= NX; k++)	fprintf(fp1, "%f   %f \n", log10(X[k]), log10(X[k]*X[k]*X[k]*X[k]*f0[k]));

    //for (k=0; k <= NX; k++) fprintf(fp0, "%f      %f  \n", log10(X[k]), log10(U[k]));

	//printf("Rsub= %f   Rtot = %f  Rprec = %f eff = %.10f gamma_eff = %f log(T2) = %f U(pmax) = %f \n", Rsub, Rtot, Rprec, eff, gamma_eff, log10(T2), U[NX]);
	//printf("f0[0] = %f \n", f0[0]);
	//printf("Slope of log-log = %f \n", (log10(F[NX]) - log10(F[0]))/(log10(X[NX]) - log10(X[0])));
    //printf("Koraka = %d \n", Koraci);
	//printf("Rsub = %.9f  Rprec = %f  Rtot = %f  U[pmax] = %f \n", Rsub, Rprec, Rtot, U[NX]);
	//fclose(fp0);	
	//fclose(fp1);
	//fclose(fp2);

    //maximum momentum for electrons calculation
    // Morlino, Amato and Blasi (2009)

    MIN_P = 1.5 * sqrt(me*me*me*clight*clight*clight*clight*1.e20/(qel*B1*re)) * (Vshock/clight) * U[0] * sqrt((1 - 1/(Rtot * U[0]))/(1 + sqrt(11)*Rtot*U[0]));
    MIN_P = fabs(X[0]*mp*clight*1.e+5 - MIN_P);
    MIN_e = X[0];

    for (i=1; i <= NX; i++)
	{
	MIN_Pe = 1.5 * sqrt(me*me*me*clight*clight*clight*clight*1.e20/(qel*B1*re)) * (Vshock/clight) * U[i] * sqrt((1 - 1/(Rtot * U[i]))/(1 + sqrt(11)*Rtot*U[i]));
    MIN_Pe = fabs(X[i]*mp*clight*1.e+5 - MIN_Pe);

    if (MIN_Pe < MIN_P) 
    	{
    		MIN_e = X[i];
    		MIN_P = MIN_Pe;
    	}
	}
	//-------------------------------------------

	*e_max = MIN_e; /// menjano!!!
	*gamma_eff_out = gamma_eff;
	*Rtot_out = Rtot;
	*p_max = XEND;
	//quasi-linear model from Caprioli et al. 2009 + 2010
	//*B2 = B * sqrt(0.5*(1 - zeta)*MA*(1 - 1/(Rtot*Rtot))*pow(Rtot,1.5)  + 1.0);

	//izmenjen prema Liju
	*B2 = sqrt(1.0/3.0 + 2.0 * Rsub * Rsub/3.0) * B1;

	//double pom = B * sqrt(0.5*(1 - zeta)*MA*(1 - 1/(Rtot*Rtot))*pow(Rtot,1.5)  + 1.0);
	double pom = 0.5 * sqrt(1.0/3.0 + 2.0 * Rsub * Rsub/3.0) * B1; // = 0.5 * B2 = Bnorm, normal to LOS

	*izlaz_L = (Luminosity(1.0, pom, Radius, f0, X, MIN_e, Chi_esc, Vshock, rho0)/8.644677688e+45)/1.0e-23; //Luminoznost u Jy (Jansky) na 8.5 kpc menjano !!!


	double x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,xs,ys,slope=0.0;
        x1=0.5;
        x2=0.8;
        x3=1.4;
        x4=3.0;
        x5=5.0;

    y1=(Luminosity(x1, pom, Radius, f0, X, MIN_e, Chi_esc, Vshock, rho0)/8.644677688e+45)/1.0e-23;
	y2=(Luminosity(x2, pom, Radius, f0, X, MIN_e, Chi_esc, Vshock, rho0)/8.644677688e+45)/1.0e-23;
	y3=(Luminosity(x3, pom, Radius, f0, X, MIN_e, Chi_esc, Vshock, rho0)/8.644677688e+45)/1.0e-23;
	y4=(Luminosity(x4, pom, Radius, f0, X, MIN_e, Chi_esc, Vshock, rho0)/8.644677688e+45)/1.0e-23;
	y5=(Luminosity(x5, pom, Radius, f0, X, MIN_e, Chi_esc, Vshock, rho0)/8.644677688e+45)/1.0e-23;

	x1=log10(x1);
	x2=log10(x2);
	x3=log10(x3);
	x4=log10(x4);
	x5=log10(x5);

	y1=log10(y1);
	y2=log10(y2);
	y3=log10(y3);
	y4=log10(y4);
	y5=log10(y5);

	xs=(x1+x2+x3+x4+x5)/5.0;
	ys=(y1+y2+y3+y4+y5)/5.0;

	slope = (x1-xs)*(y1-ys)+(x2-xs)*(y2-ys)+(x3-xs)*(y3-ys)+(x4-xs)*(y4-ys)+(x5-xs)*(y5-ys);
	slope = slope/((x1-xs)*(x1-xs)+(x2-xs)*(x2-xs)+(x3-xs)*(x3-xs)+(x4-xs)*(x4-xs)+(x5-xs)*(x5-xs));
	
	*alfa = slope;


}
