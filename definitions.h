#define  PHYSICS                 HD
#define  DIMENSIONS              3
#define  COMPONENTS              3
#define  GEOMETRY                CARTESIAN
#define  BODY_FORCE              NO
#define  COOLING                 NO
#define  RECONSTRUCTION          LINEAR
#define  TIME_STEPPING           RK2
#define  DIMENSIONAL_SPLITTING   NO
#define  NTRACER                 1
#define  USER_DEF_PARAMETERS     9

/* -- physics dependent declarations -- */

#define  EOS                     IDEAL
#define  ENTROPY_SWITCH          NO
#define  THERMAL_CONDUCTION      NO
#define  VISCOSITY               NO
#define  ROTATING_FRAME          NO

/* -- user-defined parameters (labels) -- */

#define  Mej                     0
#define  Ek                      1
#define  t0                      2
#define  Ts                      3
#define  RHOamb                  4
#define  Vsh                     5
#define  xc                      6
#define  RHOs                    7
#define  rs                      8

/* [Beg] user-defined constants (do not change this line) */

#define  CHOMBO_REF_VAR          -1
#define  UNIT_DENSITY            (1.4*CONST_mp)
#define  UNIT_VELOCITY           1.0e8
#define  UNIT_LENGTH             CONST_pc

/* [End] user-defined constants (do not change this line) */

/* -- supplementary constants (user editable) -- */ 

#define  INITIAL_SMOOTHING   NO
#define  WARNING_MESSAGES    YES
#define  PRINT_TO_FILE       YES
#define  INTERNAL_BOUNDARY   YES
#define  SHOCK_FLATTENING    MULTID
#define  CHAR_LIMITING       YES
#define  LIMITER             DEFAULT
