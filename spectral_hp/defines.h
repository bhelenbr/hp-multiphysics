/* DEFAULTS DON'T CHANGE THESE HERE */
/* OVERRIDE AT BOTTOM oF FILE */

/* SIMULATION PV3VIEWER */
#define SIMULATION 

/* KOVASZNAY TEST CYLINDER SPHERE FREESTREAM TWOLAYER CONVDIFF UNSTEADY_DROP DROP TWOLAYERNOINERTIA */
#define FREESTREAM

/* CIRCLE SIN COS NACA */
#define CIRCLE

/* 1 2 3 */
#define NV 3

/* 1 2 3 */
#define MXSTEP 2

/* CONSERV  NONCONSERV */
#define CONSERV

/* AXISYMMETRIC  NO_AXISYMMETRIC */
#define NO_AXISYMMETRIC

/* INERTIALESS NO_INERTIALESS */
#define NO_INERTIALESS

/* TIMEACCURATE NO_TIMEACCURATE */
#define NO_TIMEACCURATE  

/* OLDRECONNECT NEWRECONNECT */
#define OLDRECONNECT 

/* FINE_ERROR NO_FINE_ERROR */
#define NO_FINE_ERROR

/* TWO_LEVEL NO_TWO_LEVEL */
#define NO_TWO_LEVEL

/* DEFORM  NO_DEFORM */
#define DEFORM

/* FOURTH NO_FOURTH */
#define NO_FOURTH

/* GEOMETRIC NO_GEOMETRIC */
#define GEOMETRIC

/* DEBUG NO_DEBUG */
#define NO_DEBUG 

/* ENERGY NO_ENERGY */
#define NO_ENERGY

/* AUTOMATIC SWITCHES FOR PARTICULAR CASES */
#ifdef CIRCLE
#define RADIUS 0.5;
#define CENTERX 0.0;
#define CENTERY 0.0;
#endif

#define AMP 0.01;
#define LAM 1.0;

#ifdef COS
#define AMP 0.375;
#define LAM 1.0;
#endif

#if (defined(SPHERE) || defined(KOVASZNAY) || defined(CYLINDER) || defined(CONVDIFF) || defined(TEST) || defined(FREESTREAM))
#undef DEFORM
#define NO_DEFORM
#endif

#if (defined(SPHERE) || defined(DROP) || defined(UNSTEADY_DROP))
#undef NO_AXISYMMETRIC
#define AXISYMMETRIC
#endif

#ifdef UNSTEADY_DROP
#define AMP 0.0
#endif

#ifdef TWOLAYER
#undef CIRCLE
#define SIN
#endif

#ifdef TWOLAYERNOINERTIA
#define TWOLAYER
#define INERTIALESS
#undef CIRCLE
#define SIN
#endif

#ifdef CONVDIFF
#define NV 1
#define CASE 0
#define AMP 1.0
#define LAM 0.0
#if (CASE == 2)
#undef NO_TIMEACCURATE
#define TIMEACCURATE
#endif
#endif

#if (defined(DROP) || defined(UNSTEADY_DROP))
#undef NO_ENERGY
#define ENERGY
#endif










