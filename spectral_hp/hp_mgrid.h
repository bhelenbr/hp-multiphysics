/*
 *  hp_mgrid.h
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "spectral_hp.h"

#define MXLG2P 5
#define NSTAGE 5

#define HP_MGRID_MP (PRDX_MASK +PRDY_MASK +COMX_MASK +COMY_MASK)

/* THESE THINGS ARE SHARED BY ALL MESHES OF THE SAME BLOCK */
struct hp_mgrid_glbls {

/*	SOLUTION STORAGE ON FIRST ENTRY TO NSTAGE */
   FLT (*vug0)[NV];
   FLT (*sug0)[NV];
   FLT (*iug0)[NV];

/*	RESIDUAL STORAGE */
   FLT (*vres)[NV];
   FLT (*sres)[NV];
   FLT (*ires)[NV];

/*	VISCOUS FORCE RESIDUAL STORAGE */
   FLT (*vvf)[NV];
   FLT (*svf)[NV];
   FLT (*ivf)[NV];   
   
/*	RESIDUAL STORAGE FOR ENTRY TO MULTIGRID */
   FLT (*vres0)[NV];
   FLT (*sres0)[NV];
   FLT (*ires0)[NV];

/*	PRECONDITIONER  */
   FLT *gam,*dtstar;
   FLT *vdiagv, *vdiagp;
   FLT *sdiagv, *sdiagp;
   
/* STABILIZATION */
   FLT *tau,*delt,adis;

/*	UNSTEADY SOURCE TERMS (NEEDED ON FINE MESH ONLY) */
   FLT (*ug0)[NV], (*ug1)[NV], (*ug2)[NV], *jcb1, *jcb2;
   FLT ***dudt[ND], ***cdjdt;

/*	ITERATION PARAMETERS */
   FLT fadd, flowcfl[MXLG2P];  
      
/*	PHYSICAL CONSTANTS */
   FLT rho, rhoi, mu, nu, sigma, g;

/*	INITIALIZATION AND BOUNDARY CONDITION FUNCTION */
   FLT (*func)(int n, FLT x, FLT y);
   
/*	OTHER CONSTANTS */
   FLT charyes;  // USE CHARACTERISTIC FAR-FIELD B.C'S
   
};

class hp_mgrid : public spectral_hp {
   protected:
/*		THINGS SHARED BY ALL HP_MGRIDS (STATIC) */
      static const FLT alpha[NSTAGE+1] = {0.25, 1./6., .375, .5, 1.0, 1.0};
      static const FLT beta[NSTAGE+1] = {1.0, 0.0, 5./9., 0.0, 4./9., 1.0};
      static FLT **cv00,**cv01,**cv10,**cv11;
      static FLT **e00,**e01,**e10,**e11;
      static FLT dt0, dt1, dt2, dt3;
      static int size;

/*		TELLS WHICH P WE ARE ON FOR P MULTIGRID */
      int log2p;
      
/*		THINGS SHARED BY HP_MGRIDS IN SAME BLOCK */
      struct hp_mgrid_glbls gbl;
      
/*		THINGS NEEDED ON EACH HP_MGRID MESH FOR MGRID */
      FLT (*vug_frst)[NV];
      FLT (*vdres[MXLG2P])[NV];
      FLT (*sdres[MXLG2P])[NV];
      FLT (*idres[MXLG2P])[NV];
      bool isfrst;
      
/*		MGRID MESH POINTERS */
      class hp_mgrid *cmesh;
      class hp_mgrid *fmesh;
      
/*		SOME PRIVATE UTILITY FUNCTIONS */
      void restouht_bdry(int tind); // USED IN MINVRT
      
   public:
      void allocate(struct hp_mgrid_glbls& ginit, int mgrid);
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

/*		SURFACE BOUNDARY ROUTINES */      
      void surfrsdl(int bnum, int mgrid);
      
/*		PARTS FOR 5 STEP UPDATE */
      void nstage1();
      void nstage2(int stage);
      
/*    MGRID TRANSFER */
      void getfres();
      void getcchng();
      int setfine(class hp_mgrid& tgt);
      int setcoarse(class hp_mgrid& tgt);  
};


