#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cstdio>
#include "PatchPluto.H" 
#include "LevelPluto.H"
#include "BoxIterator.H"
#include "AMRIO.H"
#include "SPMD.H"
#include "LoHiSide.H"
#include "CH_Timer.H"
#include "blasi.h"

#include "NamespaceHeader.H"

// Constructor - set up some defaults
LevelPluto::LevelPluto()
{
  m_dx           = 0.0;
  m_dl_min       = 1.e30;
  m_refineCoarse = 0;
  m_patchPluto   = NULL;
  m_isDefined    = false;
}

// Destructor - free up storage
LevelPluto::~LevelPluto()
{
  if (m_patchPluto != NULL)
    {
      delete m_patchPluto;
    }
}

extern double g_rshock;

void LevelPluto::getPrimitiveVars_MP (Data_Arr U, Data *d, Grid *grid)
/*!
 * - Recover primitive variables from the input conservative array \c U
 * - Set physical boundary conditions and convert
 *
 * \date June 25, 2015
 *********************************************************************** */
{
  int i,j,k;
  int dir, err;
  int nx, ny, nz;
  int lft_side[3] = {0,0,0}, rgt_side[3]={0,0,0};
  static unsigned char ***flagEntr;
  RBox cbox, *box;

  nx = grid[IDIR].np_tot;
  ny = grid[JDIR].np_tot;
  nz = grid[KDIR].np_tot;

/* ------------------------------------------------------- 
     Check whether the patch touches a physical boundary
   ------------------------------------------------------- */

  for (dir = 0; dir < DIMENSIONS; dir++){
    lft_side[dir] = (grid[dir - IDIR].lbound != 0);
    rgt_side[dir] = (grid[dir - IDIR].rbound != 0);
  }

/* ---------------------------------------------------
    Extract the portion of the domain where U 
    is defined (i.e. NOT in the physical boundary).
   --------------------------------------------------- */

  cbox.ib = 0; cbox.ie = nx - 1;
  cbox.jb = 0; cbox.je = ny - 1;
  cbox.kb = 0; cbox.ke = nz - 1;

/* -------------------------------------------------
    Exclude physical boundaries since the 
    conservative vector U is not yet defined.    
   ------------------------------------------------- */

  D_EXPAND(if (lft_side[IDIR]) cbox.ib = IBEG;  ,
           if (lft_side[JDIR]) cbox.jb = JBEG;  ,
           if (lft_side[KDIR]) cbox.kb = KBEG;)

  D_EXPAND(if (rgt_side[IDIR]) cbox.ie = IEND;  ,
           if (rgt_side[JDIR]) cbox.je = JEND;  ,
           if (rgt_side[KDIR]) cbox.ke = KEND;)

/* ----------------------------------------------------------
    Convert conservative variables into primitive variables.
    Normally this operation is performed by using total
    energy density.
    However, when the ENTROPY_SWITCH is enabled, we force
    conversion to be done from the entropy by artificially
    setting flagEntr.
   ---------------------------------------------------------- */

#if ENTROPY_SWITCH
  if (flagEntr == NULL) {
    flagEntr = ARRAY_3D(NX3_MAX, NX2_MAX, NX1_MAX, unsigned char);
    for (k = 0; k < NX3_MAX; k++){
    for (j = 0; j < NX2_MAX; j++){
    for (i = 0; i < NX1_MAX; i++){
      flagEntr[k][j][i] = 0;
      flagEntr[k][j][i] |= FLAG_ENTROPY;
    }}}
  }
/*
  BOX_LOOP(&cbox,k,j,i){
    flagEntr[k][j][i]  = 0;
    flagEntr[k][j][i] |= FLAG_ENTROPY;
  }
*/
  ConsToPrim3D(U, d->Vc, flagEntr, &cbox);
#else
  ConsToPrim3D(U, d->Vc, d->flag, &cbox);
#endif

 Boundary (d, ALL_DIR, grid);

/* --------------------------------------------------------------
    Convert primitive variables to conservative in the ghost 
    zones.
   -------------------------------------------------------------- */
/*
#if INTERNAL_BOUNDARY == YES
  box = GetRBox(TOT, CENTER);
  PrimToCons3D(d->Vc, U, box);
#else
  if (lft_side[IDIR]) {
    box = GetRBox(X1_BEG, CENTER);
    PrimToCons3D(d->Vc, U, box);
  }

  if (lft_side[JDIR]) {
    box = GetRBox(X2_BEG, CENTER);
    PrimToCons3D(d->Vc, U, box);
  }

  if (lft_side[KDIR]) {
    box = GetRBox(X3_BEG, CENTER);
    PrimToCons3D(d->Vc, U, box);
  }

  if (rgt_side[IDIR]) {
    box = GetRBox(X1_END, CENTER);
    PrimToCons3D(d->Vc, U, box);
  }

  if (rgt_side[JDIR]) {
    box = GetRBox(X2_END, CENTER);
    PrimToCons3D(d->Vc, U, box);
  }

  if (rgt_side[KDIR]) {
    box = GetRBox(X3_END, CENTER);
    PrimToCons3D(d->Vc, U, box);
  }
#endif    */
}


// Define the object so that time stepping can begin
void LevelPluto::define(const DisjointBoxLayout&  a_thisDisjointBoxLayout,
                        const DisjointBoxLayout&  a_coarserDisjointBoxLayout,
                        const ProblemDomain&      a_domain,
                        const int&                a_refineCoarse,
                        const int&                a_level,
                        const Real&               a_dx,
                        const PatchPluto*         a_patchPlutoFactory,
                        const bool&               a_hasCoarser,
                        const bool&               a_hasFiner)
{
  CH_TIME("LevelPluto::define");

  // Sanity checks
  CH_assert(a_refineCoarse > 0);
  CH_assert(a_dx > 0.0);

  // Make a copy of the current grids
  m_grids  = a_thisDisjointBoxLayout;

  // Cache data
  m_dx = a_dx;
  m_level = a_level;
  m_domain = a_domain;
  m_refineCoarse = a_refineCoarse;
  m_hasCoarser = a_hasCoarser;
  m_hasFiner = a_hasFiner;

  // Remove old patch integrator (if any), create a new one, and initialize
  if (m_patchPluto != NULL)
    {
      delete m_patchPluto;
    }

 // Determing the number of ghost cells necessary here
  m_numGhost = GetNghost();

  m_patchPluto = a_patchPlutoFactory->new_patchPluto();
  m_patchPluto->define(m_domain,m_dx,m_level,m_numGhost);
 
  // Set the grid for the entire level
  setGridLevel();

  // Get the number of conserved variable and face centered fluxes
  m_numCons   = m_patchPluto->numConserved();
  m_numFluxes = m_patchPluto->numFluxes();

  m_exchangeCopier.exchangeDefine(a_thisDisjointBoxLayout,
                                  m_numGhost*IntVect::Unit);

  // Setup an interval corresponding to the conserved variables
  Interval UInterval(0,m_numCons-1);

 #if (TIME_STEPPING == RK2)
  // Create temporary storage with a layer of "m_numGhost" ghost cells
  // for the flags passing from predictor to corrector (RK2 only)
  m_Flags.define(m_grids,1,m_numGhost*IntVect::Unit);
  m_Utmp.define(m_grids,m_numCons,m_numGhost*IntVect::Unit);
 #endif

  // Create temporary storage with a layer of "m_numGhost" ghost cells
  {
    CH_TIME("setup::Udefine");
    m_U.define(m_grids,m_numCons,m_numGhost*IntVect::Unit);
  }

  // Set up the interpolator if there is a coarser level
  if (m_hasCoarser)
    {
      m_patcher.define(a_thisDisjointBoxLayout,
                       a_coarserDisjointBoxLayout,
                       m_numCons,
                       coarsen(a_domain,a_refineCoarse),
                       a_refineCoarse,
                       m_dx,
                       m_numGhost);
    }

  // Everything is defined
  m_isDefined = true;
}

// Advance the solution by "a_dt" by using an unsplit method.
// "a_finerFluxRegister" is the flux register with the next finer level.
// "a_coarseFluxRegister" is flux register with the next coarser level.
// If source terms do not exist, "a_S" should be null constructed and not
// defined (i.e. its define() should not be called).




Real LevelPluto::step(LevelData<FArrayBox>&       a_U,
                      LevelData<FArrayBox>        a_flux[CH_SPACEDIM],
                      LevelFluxRegister&          a_finerFluxRegister,
                      LevelFluxRegister&          a_coarserFluxRegister,
                      LevelData<FArrayBox>&       a_split_tags,
                      const LevelData<FArrayBox>& a_UCoarseOld,
                      const Real&                 a_TCoarseOld,
                      const LevelData<FArrayBox>& a_UCoarseNew,
                      const Real&                 a_TCoarseNew,
                      const Real&                 a_time,
                      const Real&                 a_dt,
                      const Real&                 a_cfl)
{
  double EPS_PSHOCK_FLATTEN = 5.0;
  static double ***pt;
  double divv, gradp, pt_min;
  double *x1, *x2, *x3;
  int i,j,k, nv;
  double Rshock=0.0, radius, RRshock, sum;
  int world_size, world_rank;
  double ***UU[NVAR];
  static Data d;
  double vel, Velocity, VVelocity;

  double dpx1, *dx1, *dV1, ***vx1, pt_min1, dvx1;
  double dpx2, *dx2, *dV2, ***vx2, pt_min2, dvx2;
  double dpx3, *dx3, *dV3, ***vx3, pt_min3, dvx3;

  double time_c = (UNIT_LENGTH / UNIT_VELOCITY)/(86400.0*365.2422);
  double unit_pressure = UNIT_DENSITY * UNIT_VELOCITY * UNIT_VELOCITY;
  double SoundSpeed = sqrt(g_inputParam[Ts] * g_gamma * CONST_kB/(1.4 * CONST_mp))/UNIT_VELOCITY;
  double Vshock, Mach;

  double izlaz_gamma, izlaz_Rtot, izlaz_pmax, izlaz_emax, izlaz_B2, izlaz_L, izlaz_alfa;



  CH_TIMERS("LevelPluto::step");

  CH_TIMER("LevelPluto::step::setup"   ,timeSetup);
  CH_TIMER("LevelPluto::step::update"  ,timeUpdate);
  CH_TIMER("LevelPluto::step::reflux"  ,timeReflux);
  CH_TIMER("LevelPluto::step::conclude",timeConclude);

  // Make sure everything is defined
  CH_assert(m_isDefined);

  CH_START(timeSetup);

  // Clear flux registers with next finer level
  if (m_hasFiner && (g_intStage == 1))
    {
      a_finerFluxRegister.setToZero();
    }

  // Setup an interval corresponding to the conserved variables
  Interval UInterval(0,m_numCons-1);

  {
    CH_TIME("setup::localU");
    for (DataIterator dit = m_U.dataIterator(); dit.ok(); ++dit)
    {
      m_U[dit].setVal(0.0); // Gets rid of denormalized crap.
      m_U[dit].copy(a_U[dit]);
    }

    m_U.exchange(m_exchangeCopier);
  }

  // Fill m_U's ghost cells using fillInterp
  if (m_hasCoarser)
    {
      // Fraction "a_time" falls between the old and the new coarse times
      Real alpha = (a_time - a_TCoarseOld) / (a_TCoarseNew - a_TCoarseOld);

    // Truncate the fraction to the range [0,1] to remove floating-point
    // subtraction roundoff effects
      Real eps = 0.04 * a_dt / m_refineCoarse;

      if (Abs(alpha) < eps)     alpha = 0.0;
      if (Abs(1.0-alpha) < eps) alpha = 1.0;

      // Current time before old coarse time
      if (alpha < 0.0)
        {
          MayDay::Error( "LevelPluto::step: alpha < 0.0");
        }

      // Current time after new coarse time
      if (alpha > 1.0)
        {
          MayDay::Error( "LevelPluto::step: alpha > 1.0");
        }

      // Interpolate ghost cells from next coarser level using both space
      // and time interpolation
      m_patcher.fillInterp(m_U,
                           a_UCoarseOld,
                           a_UCoarseNew,
                           alpha,
                           0,0,m_numCons);
    }

  // Potentially used in boundary conditions
  m_patchPluto->setCurrentTime(a_time);

  // Use to restrict maximum wave speed away from zero
  Real maxWaveSpeed = 1.e-12;
  Real minDtCool    = 1.e38;

  // The grid structure
  Grid *grid;
  static Time_Step Dts;
  Real inv_dt;
  
  #ifdef GLM_MHD
   glm_ch = g_coeff_dl_min*m_dx/(a_dt + 1.e-16)*a_cfl;
//   glm_ch = g_coeff_dl_min/(a_dt + 1.e-16)*a_cfl; /* If subcycling is turned off */
   glm_ch = MIN(glm_ch,glm_ch_max*g_coeff_dl_min);
  #endif

  CH_STOP(timeSetup);
  g_level_dx = m_dx;


  // Beginning of loop through patches/grids.
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit){
    CH_START(timeUpdate);

    // The current box
    Box curBox = m_grids.get(dit());

    // The current grid of conserved variables
    FArrayBox& curU = m_U[dit];

    // The current grid of volumes
    #if GEOMETRY != CARTESIAN
     const FArrayBox& curdV = m_dV[dit()];
    #else
     const FArrayBox  curdV;
    #endif
 
   #ifdef SKIP_SPLIT_CELLS
    // The current grid of split/unsplit tags
    FArrayBox& split_tags = a_split_tags[dit];
   #else
    FArrayBox split_tags;
   #endif

   #if (TIME_STEPPING == RK2)
    // The current storage for flags (RK2 only)
    BaseFab<unsigned char>& flags = m_Flags[dit];
    // Local temporary storage for conserved variables
    FArrayBox& curUtmp = m_Utmp[dit];

   #else
    BaseFab<unsigned char> flags;
    FArrayBox curUtmp;
   #endif

    // The fluxes computed for this grid - used for refluxing and returning
    // other face centered quantities
    FluxBox flux;

    // Set the current box for the patch integrator
    m_patchPluto->setCurrentBox(curBox);

    Real minDtCoolGrid;
    
    grid = m_structs_grid[dit].getGrid();

    x1 = grid[IDIR].x;
    x2 = grid[JDIR].x;
    x3 = grid[KDIR].x;
 
    IBEG = grid[IDIR].lbeg; IEND = grid[IDIR].lend;
    JBEG = grid[JDIR].lbeg; JEND = grid[JDIR].lend;
    KBEG = grid[KDIR].lbeg; KEND = grid[KDIR].lend;

    NX1 = grid[IDIR].np_int;
    NX2 = grid[JDIR].np_int;
    NX3 = grid[KDIR].np_int;

    NX1_TOT = grid[IDIR].np_tot;
    NX2_TOT = grid[JDIR].np_tot;
    NX3_TOT = grid[KDIR].np_tot;

    dV1 = grid[IDIR].dV;
    dV2 = grid[JDIR].dV; 
    dV3 = grid[KDIR].dV; 

    
    SetRBox();  /* RBox structures must be redefined for each patch */
        
    g_dt   = a_dt;
    g_time = a_time;
    g_maxRiemannIter = 0;
    PLM_CoefficientsSet (grid);  /* -- these may be needed by
                                       shock flattening algorithms */
    #if RECONSTRUCTION == PARABOLIC
     PPM_CoefficientsSet (grid);  
    #endif
    
    // reset time step coefficients 
    if (Dts.cmax == NULL) Dts.cmax = ARRAY_1D(NMAX_POINT, double);
    int id;
    Dts.inv_dta = 1.e-18;
    Dts.inv_dtp = 1.e-18;
    Dts.dt_cool = 1.e18;
    Dts.cfl     = a_cfl;
    Where(-1, grid); /* -- store grid for subsequent calls -- */

    //============================================================================= 
  { //----------krecemo kao block, da ne bi mucili memoriju
  int in;
  int nxf, nyf, nzf, indf;
  int nxb, nyb, nzb;

  double ***UU[NVAR];
  double *inv_dl, dl2, cylr;
  static Data d;
#ifdef SKIP_SPLIT_CELLS
  double ***splitcells;
#endif
#if (PARABOLIC_FLUX & EXPLICIT)
  static double **dcoeff;
#endif   
  Index indx;
  static State_1D state;

#if TIME_STEPPING == RK2 
  double wflux = 0.5;
#else
  double wflux = 1.;
#endif
  RBox *rbox = GetRBox(DOM, CENTER);
  

  if (NX1_TOT > NMAX_POINT || NX2_TOT > NMAX_POINT || NX3_TOT > NMAX_POINT){
    print ("! advanceStep(): need to re-allocate matrix\n");
    QUIT_PLUTO(1);
  }

#if GEOMETRY != CARTESIAN
  for (nv = 0; nv < NVAR; nv++) curU.divide(curdV,0,nv);
  #if CHOMBO_CONS_AM == YES
    #if ROTATING_FRAME == YES
      Box curBox = curU.box();
      for(BoxIterator bit(curBox); bit.ok(); ++bit) {
        const IntVect& iv = bit();
        curU(iv,iMPHI) /= curdV(iv,1);
        curU(iv,iMPHI) -= curU(iv,RHO)*curdV(iv,1)*g_OmegaZ;
      }
    #else
      curU.divide(curdV,1,iMPHI);
    #endif
  #endif
#else
  if (g_stretch_fact != 1.) curU /= g_stretch_fact; 
#endif

  for (nv = 0; nv < NVAR; nv++){
    UU[nv] = ArrayMap(NX3_TOT, NX2_TOT, NX1_TOT, curU.dataPtr(nv));
  }
#ifdef SKIP_SPLIT_CELLS
  splitcells = ArrayBoxMap(KBEG, KEND, JBEG, JEND, IBEG, IEND,
                           split_tags.dataPtr(0));
#endif
#if (TIME_STEPPING == RK2)
  d.flag = ArrayCharMap(NX3_TOT, NX2_TOT, NX1_TOT,flags.dataPtr(0));
#endif
#if RESISTIVITY != NO
  if (d.J == NULL) d.J = ARRAY_4D(3,NX3_MAX, NX2_MAX, NX1_MAX, double);
#endif


  if (state.flux == NULL){
    MakeState (&state);
    nxf = nyf = nzf = 1;
    D_EXPAND(nxf = NMAX_POINT;  ,
             nyf = NMAX_POINT;  ,
             nzf = NMAX_POINT;)

    d.Vc = ARRAY_4D(NVAR, nzf, nyf, nxf, double);
#if (TIME_STEPPING != RK2)
    d.flag = ARRAY_3D(nzf, nyf, nxf, unsigned char);
#endif 
#if (PARABOLIC_FLUX & EXPLICIT)
    dcoeff = ARRAY_2D(NMAX_POINT, NVAR, double);
#endif
  }

  if (g_intStage == 1) TOT_LOOP(k,j,i) d.flag[k][j][i] = 0;      

   getPrimitiveVars_MP (UU, &d, grid);

    //============================================================================= 

    for (k = KOFFSET; k < NX3_TOT-KOFFSET; k++){ 
    for (j = JOFFSET; j < NX2_TOT-JOFFSET; j++){
    for (i = IOFFSET; i < NX1_TOT-IOFFSET; i++){ 

    radius = sqrt(x1[i]*x1[i] + x2[j]*x2[j] + x3[k]*x3[k]);

    //#if GEOMETRY == CARTESIAN
     D_EXPAND(  dvx1 = d.Vc[VX1][k][j][i + 1] - d.Vc[VX1][k][j][i - 1];   ,
                 dvx2 = d.Vc[VX2][k][j + 1][i] - d.Vc[VX2][k][j - 1][i];   ,
                 dvx3 = d.Vc[VX2][k + 1][j][i] - d.Vc[VX2][k - 1][j][i]; )

    divv = D_EXPAND(dvx1/dV1[i], + dvx2/dV2[j], + dvx3/dV3[k]);

    if (divv < 0.0){

            pt_min = d.Vc[PRS][k][j][i];
      D_EXPAND(pt_min1 = MIN(d.Vc[PRS][k][j][i+1], d.Vc[PRS][k][j][i-1]); ,
               pt_min2 = MIN(d.Vc[PRS][k][j+1][i], d.Vc[PRS][k][j-1][i]);  ,
               pt_min3 = MIN(d.Vc[PRS][k+1][j][i], d.Vc[PRS][k-1][j][i]); )

      D_EXPAND(pt_min = MIN(pt_min, pt_min1);  ,
               pt_min = MIN(pt_min, pt_min2);  ,
               pt_min = MIN(pt_min, pt_min3);)
      
      D_EXPAND(dpx1 = fabs(d.Vc[PRS][k][j][i+1] - d.Vc[PRS][k][j][i-1]);  ,  
               dpx2 = fabs(d.Vc[PRS][k][j+1][i] - d.Vc[PRS][k][j-1][i]);  , 
               dpx3 = fabs(d.Vc[PRS][k+1][j][i] - d.Vc[PRS][k-1][j][i]);)   
                
      gradp = D_EXPAND(dpx1, + dpx2, + dpx3);

       if (gradp > EPS_PSHOCK_FLATTEN*pt_min)
                   Rshock = MAX(radius, Rshock);  

      }

    
      }}}

    /* -------------------------------------------------
               Free memory 
   ------------------------------------------------- */
  
  	for (nv = 0; nv < NVAR; nv++) FreeArrayMap(UU[nv]);

  	#ifdef SKIP_SPLIT_CELLS
   		FreeArrayBoxMap (splitcells, KBEG, KEND, JBEG, JEND, IBEG, IEND);
  	#endif
 
  	#if (TIME_STEPPING == RK2)
   		FreeArrayCharMap(d.flag);
  	#endif


    } // kraj bloka u kome se kao izlaz dobija Rshock, nema potrebe za free

    // Take one step
    m_patchPluto->advanceStep (curU, curUtmp, curdV, split_tags, flags, flux,
                               &Dts, curBox, grid);
 
    inv_dt = Dts.inv_dta + 2.0*Dts.inv_dtp;
    maxWaveSpeed = Max(maxWaveSpeed, inv_dt); // Now the inverse of the timestep

    minDtCool = Min(minDtCool, Dts.dt_cool/a_cfl);

    CH_STOP(timeUpdate);

    CH_START(timeReflux);

    // Do flux register updates
    for (int idir = 0; idir < SpaceDim; idir++) {
    // Increment coarse flux register between this level and the next
    // finer level - this level is the next coarser level with respect
    // to the next finer level
      if (m_hasFiner) {
        a_finerFluxRegister.incrementCoarse(flux[idir],a_dt,dit(),
                                            UInterval, UInterval,idir);
      }

      // Increment fine flux registers between this level and the next
      // coarser level - this level is the next finer level with respect
      // to the next coarser level
       if (m_hasCoarser) {
         a_coarserFluxRegister.incrementFine(flux[idir],a_dt,dit(),
                                             UInterval, UInterval,idir);
       }
    }

    CH_STOP(timeReflux);


  }
  // ================== End of loop through patches/grids ==================
   //MPI_Comm_rank(Chombo_MPI::comm, &world_rank);
   //MPI_Comm_size(Chombo_MPI::comm, &world_size);
   //MPI_Barrier(Chombo_MPI::comm); //dodao sam posle greske of "free memory"



   int maxLevel = 0;
   double ctime;
   maxLevel = atoi(ParamFileGet("Levels",1));
   int new_maxLevel;

   int result1 = MPI_Allreduce(&Rshock, &RRshock, 1, MPI_DOUBLE, MPI_MAX, Chombo_MPI::comm);
   int result2 = MPI_Allreduce(&g_time, &ctime, 1, MPI_DOUBLE, MPI_MAX, Chombo_MPI::comm);

   if ((result1 != MPI_SUCCESS) || (result2 != MPI_SUCCESS))
   { //bark!!!
      MayDay::Error("sorry, but I had a communcation error in radius calculation");
   }
   //new_maxLevel = (int) (maxLevel + 1.0 - 3.0*(ctime/10.0)); //copy from TagCells.cpp --ALWAYS!!!


   new_maxLevel = (int) log2(pow(2,maxLevel+1)/(g_rshock/0.5));
   new_maxLevel = MIN(new_maxLevel,maxLevel);


   if (m_level == new_maxLevel) g_rshock = RRshock;




   CH_START(timeConclude);

   {
    CH_TIME("conclude::copyU");
    // Now that we have completed the updates of all the patches, we copy the
    // contents of temporary storage, U, into the permanent storage, a_U.
    for(DataIterator dit = m_U.dataIterator(); dit.ok(); ++dit){
      a_U[dit].copy(m_U[dit]);
    }
   }

  // Find the minimum of dt's over this level
  Real local_dtNew = 1. / maxWaveSpeed;
  local_dtNew = Min(local_dtNew,minDtCool);
  Real dtNew;

  {
    CH_TIME("conclude::getDt");
 #ifdef CH_MPI
  #if (TIME_STEPPING == RK2) && (COOLING == NO)
  if (g_intStage == 1) {
  #endif
   int result = MPI_Allreduce(&local_dtNew, &dtNew, 1, MPI_CH_REAL,
                                  MPI_MIN, Chombo_MPI::comm);
   if(result != MPI_SUCCESS){ //bark!!!
      MayDay::Error("sorry, but I had a communcation error on new dt");
   }
  #if (TIME_STEPPING == RK2) && (COOLING == NO)
  } else {
   dtNew = local_dtNew;
  }
  #endif
 #else
   dtNew = local_dtNew;
 #endif
  }

  CH_STOP(timeConclude);

  // Return the maximum stable time step
  return dtNew;
}

void LevelPluto::setGridLevel()
{

 CH_TIME("LevelPluto::setGrid");

 m_structs_grid.define(m_grids);

 #if GEOMETRY != CARTESIAN
  m_dV.define(m_grids,CHOMBO_NDV,m_numGhost*IntVect::Unit);
 #endif

 Real dlMinLoc = 1.e30;

 for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
   {
      // The current box
      Box curBox = m_grids.get(dit());
      struct GRID* grid = m_structs_grid[dit].getGrid();

      #if GEOMETRY != CARTESIAN 
       FArrayBox& curdV = m_dV[dit()];    
      #else
       FArrayBox  curdV; 
      #endif
       
      m_patchPluto->setGrid(curBox, grid, curdV);           
       
      for (int idir = 0; idir < SpaceDim; idir++) dlMinLoc = Min(dlMinLoc,grid[idir].dl_min);   
   }

#if (GEOMETRY == CARTESIAN) || (GEOMETRY == CYLINDRICAL)

   D_EXPAND(m_dl_min = m_dx; ,
            m_dl_min = MIN(m_dl_min,m_dx*g_x2stretch); ,
            m_dl_min = MIN(m_dl_min,m_dx*g_x3stretch); )
#else

 #ifdef CH_MPI
  Real dlMin;
  int result = MPI_Allreduce(&dlMinLoc, &dlMin, 1, MPI_CH_REAL,
                             MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){ //bark!!!
   MayDay::Error("sorry, but I had a communcation error on dlMin");
  }
  m_dl_min = dlMin;
 #else
  m_dl_min = dlMinLoc;
 #endif

#endif

}

#if GEOMETRY != CARTESIAN
const LevelData<FArrayBox>& LevelPluto::getdV() const
{
  return m_dV;
}
#endif

Real LevelPluto::getDlMin()
{

 CH_TIME("LevelPluto::getDlMin");

 return m_dl_min / m_dx;
// return m_dl_min; /* If subcycling is turned off */


}

#include "NamespaceFooter.H"
