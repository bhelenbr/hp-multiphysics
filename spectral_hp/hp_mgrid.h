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
#define MXSTEP 3

#define HP_MGRID_MP (COMX_MASK +COMY_MASK)

/* THESE THINGS ARE SHARED BY ALL MESHES OF THE SAME BLOCK */
struct hp_mgrid_glbls {

/*	SOLUTION STORAGE ON FIRST ENTRY TO NSTAGE */
   struct vsi ug0;

/*	RESIDUAL STORAGE */
   struct vsi res;

/*	VISCOUS FORCE RESIDUAL STORAGE */
   struct vsi vf;  
   
/*	RESIDUAL STORAGE FOR ENTRY TO MULTIGRID */
   struct vsi res0;

/*	PRECONDITIONER  */
   FLT *gam,*dtstar;
   FLT *vdiagv, *vdiagp;
   FLT *sdiagv, *sdiagp;
   
/* STABILIZATION */
   FLT *tau,*delt;

/*	UNSTEADY SOURCE TERMS (NEEDED ON FINE MESH ONLY) FOR BACKWARDS DIFFERENCE */
   struct vsi ugbd[MXSTEP-1]; // BACKWARDS DIFFERENCE FLOW INFO
   FLT (*vrtxbd[MXSTEP-1])[ND]; // BACKWARDS DIFFERENCE MESH INFO (TO CALCULATE MESH VELOCITY)
   struct bistruct *binfobd[MXSTEP-1][MAXSB];  /* BACKWARDS CURVED BDRY INFORMATION (FINE MESH ONLY) */
   FLT ***dugdt[NV];  // UNSTEADY SOURCE FOR FLOW (ONLY NEEDED ON FINEST MESH)
   struct bistruct *dbinfodt[MAXSB]; // UNSTEADY CURVED SIDE VELOCITY (ONLY NEEDED ON FINEST MESH)
/*	MESH DVRTDT IS NEEDED ON EACH MESH FOR NONLINEAR TERM IN NAVIER-STOKES */
      
/*	PHYSICAL CONSTANTS */
   FLT rho, rhoi, mu, nu;

/*	INITIALIZATION AND BOUNDARY CONDITION FUNCTION */
   FLT (*func)(int n, FLT x, FLT y);
};

class hp_mgrid : public spectral_hp {
   protected:
/*		THINGS SHARED BY ALL HP_MGRIDS (STATIC) */
      static const FLT alpha[NSTAGE+1] = {0.25, 1./6., .375, .5, 1.0, 1.0}; // MULTISTAGE TIME STEP CONSTANTS (IMAGINARY)
      static const FLT beta[NSTAGE+1] = {1.0, 0.0, 5./9., 0.0, 4./9., 1.0}; // MULTISTAGE TIME STEP CONSTANTS (REAL)
      static FLT **cv00,**cv01,**cv10,**cv11; // LOCAL WORK ARRAYS
      static FLT **e00,**e01,**e10,**e11; // LOCAL WORK ARRAYS
      static int nstep; // NUMBER OF STEPS IN BD SCHEME
      static FLT g, dti, time, bd[MXSTEP+1]; // GRAVITY, INVERSE TIME STEP, TIME, BACKWARDS DIFFERENCE CONSTANTS
      static FLT fadd, cfl[MXLG2P];   // ITERATION PARAMETERS  
      static FLT adis; // DISSIPATION CONSTANT
      static int charyes;  // USE CHARACTERISTIC FAR-FIELD B.C'S
      static FLT trncerr, tol;  //	ADAPTATION CONSTANTS  
      static class hp_mgrid hpstr; // STORAGE FOR ADAPTATION 
      static struct vsi ugstr[MXSTEP-1]; // STORAGE FOR UNSTEADY ADAPTATION BD FLOW INFO
      static FLT (*vrtxstr[MXSTEP-1])[ND]; // STORAGE FOR UNSTEADY ADAPTATION MESH BD INFO
      static struct bistruct *binfostr[MXSTEP-1][MAXSB]; // STORAGE FOR UNSTEADY ADAPTATION BOUNDARY BD INFO
      static FLT **bdwk[MXSTEP-1][NV]; // WORK FOR ADAPTATION
      static int size;
  
/*		TELLS WHICH P WE ARE ON FOR P MULTIGRID */
      int log2p;
      
/*		THINGS SHARED BY HP_MGRIDS IN SAME BLOCK */
      struct hp_mgrid_glbls *gbl;
      
/*		THINGS NEEDED ON EACH HP_MGRID MESH FOR MGRID */
      FLT (*dvrtdt)[ND]; // BACKWARDS DIFFERENCE MESH INFO (TO CALCULATE MESH VELOCITY)
      FLT (*vug_frst)[NV]; // SOLUTION ON FIRST ENTRY TO COARSE MESH
      struct vsi dres[MXLG2P]; // DRIVING TERM FOR MULTIGRID
      bool isfrst; // FLAG TO SET ON FIRST ENTRY TO COARSE MESH
      
/*    SURFACE BOUNDARY CONDITION STUFF */
      class surface *srf;
      
/*		SOME PRIVATE UTILITY FUNCTIONS */
      void restouht_bdry(int tind); // USED IN MINVRT
      void gbl_alloc(struct hp_mgrid_glbls *store);

   public:
      void allocate(int mgrid, struct hp_mgrid_glbls *store);
      void inline loadbasis(class hpbasis& bas) { 
         b = bas;
         log2p = 0;
         while ((b.p-1)>>log2p > 0) ++log2p;
      }

/*		CALCULATE TIMESTEP */
      void tstep1();
      void tstep_mp();
      void tstep2();
      
/*		DETERMINE SOLUTION RESIDUAL */
      void rsdl(int stage, int mgrid);
      
/*		PRINT SOLUTION RESIDUAL */
      void maxres();

/*		INVERT MASS MATRIX (4 STEP PROCESS WITH COMMUNICATION IN BETWEEN EACH STEP) */
      void minvrt1();
      void minvrt2();
      void minvrt3(int mode);
      void minvrt4();

/*		BOUNDARY CONDITION ROUTINES */
      void setinflow();
      void addbflux(int mgrid);
      void bdry_vsnd();
      void bdry_mp();
      void bdry_vrcvandzero();
      void bdry_ssnd(int mode);
      void bdry_srcvandzero(int mode);

/*		COUPLED SURFACE BOUNDARY ROUTINES */ 
      void surfvrttoug();
      void surfugtovrt1();
      void surfugtovrt2();
/*		SETUP SURFACE 1D SPRING CONSTANTS */
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
      
/*		PARTS FOR 5 STEP UPDATE */
      void nstage1();
      void nstage2(int stage);
      
/*    MGRID TRANSFER */
      void getfres();
      void getcchng();
      int setfine(class hp_mgrid& tgt);
      int setcoarse(class hp_mgrid& tgt); 

/*		FOR FINEST MESH ONLY ADVANCE TIME SOLUTION */
      void tadvance();
      void getfdvrtdt();  // TO TRANSFER MESH TIME DERIVATIVE TO COARSE MESHES */
      
/*		FUNCTIONS FOR ADAPTION */      
      void length1();
      void length_mp();
      void length2(); 
      void outlength(char *name, FILETYPE type);
      void inlength(char *name);
      void adapt(class hp_mgrid& bgn, FLT tolerance);
      
      friend class block;
      friend class blocks;
};


