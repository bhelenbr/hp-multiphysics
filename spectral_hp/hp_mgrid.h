/*
 *  hp_mgrid.h
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "spectral_hp.h"
#include "surface.h"

#define MXLG2P 5
#define NSTAGE 5
#define NSTEP 2

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
   struct vsi ugbd[NSTEP-1]; // BACKWARDS DIFFERENCE FLOW INFO
   FLT (*vrtxbd[NSTEP-1])[ND]; // BACKWARDS DIFFERENCE MESH INFO (TO CALCULATE MESH VELOCITY)
   struct bistruct *binfobd[NSTEP-1][MAXSB];  /* BACKWARDS CURVED BDRY INFORMATION (FINE MESH ONLY) */
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
      static const FLT alpha[NSTAGE+1] = {0.25, 1./6., .375, .5, 1.0, 1.0};
      static const FLT beta[NSTAGE+1] = {1.0, 0.0, 5./9., 0.0, 4./9., 1.0};
      static FLT **cv00,**cv01,**cv10,**cv11;
      static FLT **e00,**e01,**e10,**e11;
      static FLT g, dti, time, bd[NSTEP+1];
      static FLT fadd, cfl[MXLG2P];   // ITERATION PARAMETERS  
      static FLT adis;
      static int charyes;  // USE CHARACTERISTIC FAR-FIELD B.C'S
      static FLT trncerr, tol;  //	ADAPTATION CONSTANTS  
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
      
/*		MGRID MESH POINTERS */
      class hp_mgrid *cmesh;
      class hp_mgrid *fmesh;
      
/*		SOME PRIVATE UTILITY FUNCTIONS */
      void restouht_bdry(int tind); // USED IN MINVRT
      void gbl_alloc(struct hp_mgrid_glbls *store);

   public:
      void allocate(int mgrid, struct hp_mgrid_glbls *store);
      static inline void setstatics(FLT dtiin, FLT timein, FLT gin) {
         dti = dtiin;
         time = timein;
         g = gin;
         bd[0] = 1.5*dti;
         bd[1] = -2.*dti;
         bd[2] = 0.5*dti;
         bd[3] = 0.0;
      }
      void inline loadbasis(class hpbasis& bas) { 
         b = bas;
         log2p = 0;
         while ((b.p-1)>>log2p > 0) ++log2p;
      }

/*		CREATE SOURCE (FOR UNSTEADY) */
      void allocate_source();
      void dt_source(spectral_hp un0, spectral_hp un1, spectral_hp un2);

/*		CALCULATE TIMESTEP */
      void tstep1();
      void tstep2();
      
/*		DETERMINE SOLUTION RESIDUAL */
      void rsdl(int stage, int mgrid);
      void rsdlp1(int stage, int mgrid);
      
/*		MEASURE SOLUTION RESIDUAL */
      void maxres(FLT mx[NV]);

/*		INVERT MASS MATRIX (4 STEP PROCESS WITH COMMUNICATION IN BETWEEN EACH STEP) */
      void minvrt1();
      void minvrt2();
      void minvrt3(int mode);
      void minvrt4();

/*		BOUNDARY CONDITION ROUTINES */
      void setinflow();
      void addbflux(int mgrid);
      void bdry_rcvandzero(int mode);
      void bdry_snd(int mode);

/*		COUPLED SURFACE BOUNDARY ROUTINES */ 
      void surfvrttoug();
      void surfugtovrt1();
      void surfugtovrt2();
/*		SETUP SURFACE 1D SPRING CONSTANTS */
      void setksprg1d();
      void surfrsdl(int bnum, int mgrid);
      void surfinvrt1(int bnum);
      void surfinvrt2(int bnum);
      void surfdt1(int bnum);
      void surfdt2(int bnum);
      void surfnstage1(int bnum);
      void surfnstage2(int bnum, int stage);
      void surfgetfres(int bnum);
      void surfgetcchng(int bnum);
      
/*		PARTS FOR 5 STEP UPDATE */
      void nstage1();
      void nstage2(int stage);
      
/*    MGRID TRANSFER */
      void getfres();
      void getcchng();
      int setfine(class hp_mgrid& tgt);
      int setcoarse(class hp_mgrid& tgt); 

/*		SETUP MESH DENSITY FUNCTION FOR ADAPTION */      
      void density1();
      void density2();
      void outdensity(char *name);
      
/*		FOR FINEST MESH ONLY ADVANCE TIME SOLUTION */
      void tadvance();
      
      friend class block;
      friend class blocks;
};


