

#include <math.h>
#include <float.h>
#include "pluto.h"

const double lowest_double = -DBL_MAX;
const double sn_mu = 1.4;  //Berezhko & Volk (2004)
const double Wind_dens = 9.0; //from Mathematica, at 0.5 pc, in cm^{-3}

#define X1_BEG_SLOOP(k,j,i)  for (i = IBEG-1; i--; ) KTOT_LOOP(k) JTOT_LOOP(j)
#define X2_BEG_SLOOP(k,j,i)  for (j = JBEG-1; j--; ) KTOT_LOOP(k) ITOT_LOOP(i)

#define RESET       "\033[0m"
#define BLACK       "\033[30m"      /* Black */
#define RED         "\033[31m"      /* Red */
#define GREEN       "\033[32m"      /* Green */
#define YELLOW      "\033[33m"      /* Yellow */
#define BLUE        "\033[34m"      /* Blue */
#define MAGENTA     "\033[35m"      /* Magenta */
#define CYAN        "\033[36m"      /* Cyan */
#define WHITE       "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */


/* ********************************************************************* */
void Init (double *us, double x1, double x2, double x3)
/*
 *
 *
 *
 *********************************************************************** */
{
  int i;
  double SoundSpeed, P0, R_0;

  static int first_call = 1;

  int seed;
  double rand_num;

  if (first_call) {
    first_call = 0;
    seed = 1000;             /* choose a seed value */
    srand(seed);            /* initialize the randomizer */
  }

/* ------------------------------------------------------------------
         Set units
   ------------------------------------------------------------------ */

  g_gamma =5.0/3.0;
  

/* Other units, for reference only, defined as in PLUTO UG:  */
  double unit_pressure = UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY;
  double unit_time     = UNIT_LENGTH / UNIT_VELOCITY; //approx default 977.8131 years
  double unit_b        = UNIT_VELOCITY * sqrt(4.0 * CONST_PI * UNIT_DENSITY);
  double unit_mass = UNIT_DENSITY * UNIT_LENGTH * UNIT_LENGTH * UNIT_LENGTH/1.9891e33;
  g_ctime = g_inputParam[t0]/(  (UNIT_LENGTH / UNIT_VELOCITY)/(86400.0*365.2422) );

  SoundSpeed = sqrt(g_inputParam[Ts] * g_gamma * CONST_kB/(sn_mu * CONST_mp))/UNIT_VELOCITY;
  g_SoundSpeed = SoundSpeed;


/* ------------------------------------------------------------------
         Derive initial conditions
   ------------------------------------------------------------------ */

  double r     = x1;
  double theta = x2;
  double phi   = x3;

 
  r = D_EXPAND(x1*x1, + x2*x2, + x3*x3);
  r = sqrt(r); 

  /* Define velocities ---------------------------------------------- */

  double v_r     = 0.0;
  double v_theta = 0.0;
  double v_phi   = 0.0;

  /* Define pressure and density ------------------------------------ */ 

  double x_init = g_inputParam[xc];  //parameter r_s/r_i obtained numericaly in Mathematica


  double r_s = g_inputParam[rs]/UNIT_LENGTH; //default 0.25567: scaling radius
  double r_i = r_s/x_init;  //default 0.5: initial radius
  double rho_s = g_inputParam[RHOs]/UNIT_DENSITY; //from Mathematica

  //double rho_zone  = g_inputParam[RHOamb] + Wind_dens * (r_i/r) * (r_i/r);
   double rho_zone  = g_inputParam[RHOamb];

  V_limit = (4.0/3.0) * g_inputParam[Vsh]/x_init;
  g_vshock = V_limit;

  /* The current shock radius initialization */
  g_time = 1.e-6;

  first_jump = 1;

  //pressure of ISM with given ambient temperature Ts
  double pres_zone = 2*g_inputParam[RHOamb]*UNIT_DENSITY*CONST_kB*g_inputParam[Ts]/(unit_pressure*sn_mu*CONST_amu); 

  if (r <= r_i)
  {
    v_r = (r/r_s)*g_inputParam[Vsh]*1e8/UNIT_VELOCITY;
    rho_zone = rho_s*exp(1.0-(r/r_s));
    rand_num = (2.0*((double)rand()/((double)(RAND_MAX)+(double)(1)))-1.0)*0.05; // gives 5% irregularities
    rho_zone = (1.0+rand_num) * rho_zone;   /* Randomize the ejecta materal only*/
    //1% thermal and 99% kinetic energy = Etot
    pres_zone = (1./(2.*99.*10)) * (g_gamma - 1.) *  rho_zone * v_r * v_r; 
  }

  

  //double pres_zone = 2*rho_zone*UNIT_DENSITY*CONST_kB*g_inputParam[Ts]/(unit_pressure*sn_mu*CONST_amu);   //p = 2nkT, T = const everywhere
    //p = 2nkT, T = const everywhere

  P0 = 2*g_inputParam[RHOamb]*UNIT_DENSITY*CONST_kB*g_inputParam[Ts]/(unit_pressure*sn_mu*CONST_amu);

  /* radial magnetic field */
  double b_r     =  5.e-6/unit_b;
  double b_theta =  0.0;
  double b_phi   =  0.0;


/* ------------------------------------------------------------------
         Assign initial conditions
   ------------------------------------------------------------------ */


  g_eff = g_gamma;
  g_rshock = r_i;
  us[TRC] = g_gamma;
  us[RHO] = rho_zone;
  us[PRS] = pres_zone;
  
    EXPAND(us[VX1] = v_r*x1/r; , 
           us[VX2] = v_r*x2/r; , 
           us[VX3] = v_r*x3/r;)


  #if PHYSICS == MHD || PHYSICS == RMHD

    us[BX1] = b_r;
    us[BX2] = b_theta;
    us[BX3] = b_phi;

    #ifdef STAGGERED_MHD

      us[AX1] = 0.0;
      us[AX2] = 0.0;
      us[AX3] = 0.0;

    #endif

  #endif



  /*******************************************************************************/


}
/* ********************************************************************* */


 #define SHOCK_STRENGTH  1.0
 #define V_STRENGTH      0.1 //5% around forward shock position



void Analysis (const Data *d, Grid *grid)
/*
 *
 *
 *********************************************************************** */
{
  

}

#ifndef EPS_PSHOCK_FLATTEN
 #define EPS_PSHOCK_FLATTEN  5.0
#endif

#ifndef EPS_PSHOCK_ENTROPY
 #define EPS_PSHOCK_ENTROPY  0.05
#endif




void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid)
/*
 *
 *********************************************************************** */
{
  
double y0 = 1.85959;
double A1 = -0.23871;
double t1 = 9.4306;
double A2 = -0.34789;
double t2 = 151.97819;




  int i, j, k, nv, rank;
  double *x1, *x2, *x3;
  double r, r0, cs;
  double Vwind = 1.0, rho, vr;

  x1 = grid[IDIR].xgc;
  x2 = grid[JDIR].xgc;
  x3 = grid[KDIR].xgc;

   //MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 
  g_eff  = A1*exp(-g_time/t1) + A2*exp(-g_time/t2) + y0;


  if (side == 0){

    TOT_LOOP(k,j,i){

    	//printf("%d       %f \n", rank, g_rshock);
       r = sqrt(x1[i]*x1[i] + x2[j]*x2[j] + x3[k]*x3[k]);
        if ((r>=g_rshock) && (r<=1.02*g_rshock)){
          d->Vc[TRC][k][j][i] = g_eff;
        }
        }

    }



}

