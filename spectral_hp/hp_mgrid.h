/*
 *  hp_mgrid.h
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */
#include"spectral_hp.h"
#include"surface.h"


#define MXLG2P 5
#define NSTAGE 5
#define MXSTEP 2
#if (MXSTEP == 1)
#define MXSTEPM1 1
#else
#define MXSTEPM1 MXSTEP-1
#endif

#define CONSERV
#define NO_AXISYMMETRIC

#ifdef AXISYMMETRIC
#define RAD(I,J) crd[0][I][J]
#define RAD1D(I) crd[0][0][I]
#else
#define RAD(I,J) 1
#define RAD1D(I) 1
#endif

/* THESE THINGS ARE SHARED BY ALL MESHES OF THE SAME BLOCK */
struct hp_mgrid_glbls {

   /* SOLUTION STORAGE ON FIRST ENTRY TO NSTAGE */
   struct vsi ug0;

   /* RESIDUAL STORAGE */
   struct vsi res;

   /* VISCOUS FORCE RESIDUAL STORAGE */
   struct vsi vf;  
   
   /* RESIDUAL STORAGE FOR ENTRY TO MULTIGRID */
   struct vsi res0;

   /* MATRIX PRECONDITIONER  */
   FLT (*vprcn)[NV][NV];
   FLT (*sprcn)[NV][NV];
   FLT (*tprcn)[NV][NV];
   
   /* STABILIZATION */
   FLT *tau,*delt;

   /* UNSTEADY SOURCE TERMS (NEEDED ON FINE MESH ONLY) FOR BACKWARDS DIFFERENCE */
   struct vsi ugbd[MXSTEPM1]; // BACKWARDS DIFFERENCE FLOW INFO
   FLT (*vrtxbd[MXSTEPM1])[ND]; // BACKWARDS DIFFERENCE MESH INFO (TO CALCULATE MESH VELOCITY)
   struct bistruct *binfobd[MXSTEPM1][MAXSB];  /* BACKWARDS CURVED BDRY INFORMATION (FINE MESH ONLY) */
   FLT ***dugdt[NV];  // UNSTEADY SOURCE FOR FLOW (ONLY NEEDED ON FINEST MESH)
   struct bistruct *dbinfodt[MAXSB]; // UNSTEADY CURVED SIDE VELOCITY (ONLY NEEDED ON FINEST MESH)
   /* MESH DVRTDT IS NEEDED ON EACH MESH FOR NONLINEAR TERM IN NAVIER-STOKES */
      
   /* PHYSICAL CONSTANTS */
   FLT rho, rhoi, mu, nu;

   /* INITIALIZATION AND BOUNDARY CONDITION FUNCTION */
   FLT (*func)(int n, FLT x, FLT y);
};

class hp_mgrid : public spectral_hp {
   private:
      /* THINGS SHARED BY ALL HP_MGRIDS (STATIC) */
      static const FLT alpha[NSTAGE+1]; // MULTISTAGE TIME STEP CONSTANTS (IMAGINARY)
      static const FLT beta[NSTAGE+1]; // MULTISTAGE TIME STEP CONSTANTS (REAL)
      static FLT **cv00,**cv01,**cv10,**cv11; // LOCAL WORK ARRAYS
      static FLT **e00,**e01,**e10,**e11; // LOCAL WORK ARRAYS
      static FLT **(mvel[ND]); // for local mesh velocity info
      static int extrap; // NUMBER OF STEPS IN BD SCHEME
      static FLT g, dti, time, bd[MXSTEP+1]; // GRAVITY, INVERSE TIME STEP, TIME, BACKWARDS DIFFERENCE CONSTANTS
      static FLT fadd, cfl[MXLG2P];   // ITERATION PARAMETERS  
      static FLT adis; // DISSIPATION CONSTANT
      static int charyes;  // USE CHARACTERISTIC FAR-FIELD B.C'S
      static FLT trncerr, invbdryerr, vlngth_tol, adapt_tol;  //   ADAPTATION CONSTANTS  
      static int changed; // FLAG TO TELL WHEN MESH HAS CHANGED FOR PV3
      static struct vsi ugstr[MXSTEPM1]; // STORAGE FOR UNSTEADY ADAPTATION BD FLOW INFO
      static FLT (*vrtxstr[MXSTEPM1])[ND]; // STORAGE FOR UNSTEADY ADAPTATION MESH BD INFO
      static struct bistruct *binfostr[MXSTEPM1][MAXSB]; // STORAGE FOR UNSTEADY ADAPTATION BOUNDARY BD INFO
      static FLT **bdwk[MXSTEPM1][NV]; // WORK FOR ADAPTATION
      static int size;
  
      /* TELLS WHICH P WE ARE ON FOR P MULTIGRID */
      int log2p;
      
      /* THINGS SHARED BY HP_MGRIDS IN SAME BLOCK */
      struct hp_mgrid_glbls *gbl;
      
      /* THINGS NEEDED ON EACH HP_MGRID MESH FOR MGRID */
      FLT (*dvrtdt)[ND]; // BACKWARDS DIFFERENCE MESH INFO (TO CALCULATE MESH VELOCITY)
      FLT (*vug_frst)[NV]; // SOLUTION ON FIRST ENTRY TO COARSE MESH
      struct vsi dres[MXLG2P]; // DRIVING TERM FOR MULTIGRID
      
      /* SURFACE BOUNDARY CONDITION STUFF */
      class surface *srf;
      
      /* SOME PRIVATE UTILITY FUNCTIONS */
      void restouht_bdry(int tind); // USED IN MINVRT
      void gbl_alloc(struct hp_mgrid_glbls *store);
      
   private:
      bool isfrst; // FLAG TO SET ON FIRST ENTRY TO COARSE MESH

   public:
      void allocate(int mgrid, struct hp_mgrid_glbls *store);
      void inline loadbasis(class hpbasis& bas) { 
         b = bas;
         log2p = 0;
         while ((b.p-1)>>log2p > 0) ++log2p;
      }
      /* SET UP HIGHER-ORDER BACKWARD-DIFFERENCE CONSTANTS */
      static void setbd(int nsteps);

      /* CALCULATE TIMESTEP */
      void tstep1();
      void tstep_mp();
      void tstep2();
      
      /* DETERMINE SOLUTION RESIDUAL */
      void rsdl(int stage, int mgrid);
      
      /* PRINT SOLUTION RESIDUAL */
      void maxres();

      /* INVERT MASS MATRIX (4 STEP PROCESS WITH COMMUNICATION IN BETWEEN EACH STEP) */
      void minvrt1();
      void minvrt2();
      void minvrt3(int mode);
      inline void minvrt3_mp(int mode) {bdry_ssnd(mode);}
      void minvrt4();
      void minvrt_test_bgn(FLT (*func)(int, FLT, FLT));
      void minvrt_test_end();
      void minvrt_test_tstep();

      /* BOUNDARY CONDITION ROUTINES */
      void setinflow();
      void addbflux(int mgrid);
      void bdry_vsnd();
      void bdry_mp();
      void bdry_vrcvandzero();
      void bdry_ssnd(int mode);
      void bdry_srcvandzero(int mode);
      void bdrycheck1();
      void bdrycheck2();
      void drag(int bdry_id);

      /* COUPLED SURFACE BOUNDARY ROUTINES */ 
      void surfvrttoug();
      void surfugtovrt1();
      void surfugtovrt2();
      /* SETUP SURFACE 1D SPRING CONSTANTS */
      void setksprg1d();
      void surfksrc1d();
      void surfrsdl(int bnum, int mgrid);
      void surfinvrt1(int bnum);
      void surfinvrt2(int bnum);
      void surfdt1(int bnum);
      void surfdt2(int bnum);
      void surfnstage1(int bnum);
      void surfnstage2(int bnum, int stage);
      void surfgetfres(int bnum);
      void surfgetcchng(int bnum);
      void surfmaxres();
      
      /* PARTS FOR 5 STEP UPDATE */
      void nstage1();
      void nstage2(int stage);
      
      /* MGRID TRANSFER */
      void getfres();
      void getcchng();

      /* FOR FINEST MESH ONLY ADVANCE TIME SOLUTION */
      void tadvance();
      void getfdvrtdt();  // TO TRANSFER MESH TIME DERIVATIVE TO COARSE MESHES */
      
      /* FUNCTIONS FOR ADAPTION */      
      void length1();
      void length_mp();
      void length2(); 
      void outlength(char *name, FILETYPE type);
      void inlength(char *name);
      void adapt(class hp_mgrid& bgn, char *adaptfile);
      
      friend class block;
      friend class blocks;
      
#ifdef PV3
      void pvstruc(int& knode, int& kequiv, int& kcel1, int& kcel2, int& kcel3, int& kcel4, int& knptet, int &kptet,int& knblock,int &blocks,int &kphedra, int& ksurf,int& knsurf,int& hint);
      void pvcell(int &kn, int &kpoffset, int cel1[][4], int cel2[][5], int cel3[][6], int cel4[][8], int nptet[][8], int ptet[]);
      void pvgrid(int &kn, float (*xyz)[3]);
      void pvsurface(int snum, int &offset, int nsurf[][3], int scon[], int scel[][4], char tsurf[][20]);
      void pvvect(int *key, float v[][3]);
      void flotov(int &kn, struct vsi flo,int nvar, float *v);
      void logflotov(int &kn, struct vsi flo,int nvar, float *v);
#endif
};

