/*
 *  hp_mgrid.h
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "spectral_hp.h"

/* THESE THINGS ARE SHARED BY ALL MESHES OF THE SAME BLOCK */
struct hp_mgrid_glbls {

/*	RESIDUAL STORAGE */
   FLT (*gf)[NV];
   FLT (*gf0)[NV];

/*	PRECONDITIONER/STABILIZATION  */
   FLT *gam,*tau,*delt,*dtstar;
   FLT *vdiagv, *vdiagp;
   FLT *sdiagv, *sdiagp;

/*	UNSTEADY SOURCE TERMS */
   FLT (*ug0)[NV], (*ug1)[NV], (*ug2)[NV], *jcb1, *jcb2;
   FLT ***dudt[ND], ***cdjdt;
   
/*	PHYSICAL CONSTANTS */
   FLT rho, mu, nu, sigma;
}

class hp_mgrid : public spectral_hp {
   protected:
/*		THINGS SHARED BY ALL HP_MGRIDS (STATIC) */
      static const double alpha[NSTAGE+1] = {0.25, 1./6., .375, .5, 1.0, 1.0};
      static const double beta[NSTAGE+1] = {1.0, 0.0, 5./9., 0.0, 4./9., 1.0};
      static FLT **cv00,**cv01,**cv10,**cv11;
      static FLT **e00,**e01,**e10,**e11;
      static dt0, dt1, dt2, dt3;
      static int size;
         
/*		THINGS SHARED BY HP_MGRIDS IN SAME BLOCK */
      struct hp_mgrid_glbls gbl;
      
/*		THINGS NEEDED ON EACH HP_MGRID MESH FOR MGRID */
      FLT (*vug_frst)[NV];
      FLT (*sug_frst)[NV];
      FLT (*iug_frst)[NV];
      bool isfrst;
      
/*		MGRID MESH POINTERS */
      class hp_mgrid *cmesh;
      class hp_mgrid *fmesh;

   public:
      hp_mgrid() : size(0) {};
      void allocate(struct hp_mgrid_glbls);

/*		CREATE SOURCE (FOR UNSTEADY) */
      allocate_source();
      dt_source(spectral_hp un0, spectral_hp un1, spectral_hp un2);

/*		DETERMINE SOLUTION RESIDUAL */
      void rsdl();
      void rsdlp1();
      void rsdl_mp();
      void update(int lvl);

/*		CALCULATE TIMESTEP */
      void vddt();
      void vddt_mp();
      void vddti();
      
/*    MGRID TRANSFER */
      void mg_getfres();
      void mg_getcchng();
      int setfine(class hp_mgrid& tgt);
      int setcoarse(class hp_mgrid& tgt);  
      
/*		COMMUNICATION BOUNDARIES */
      void send(int MASK, FLT *base,int bgn,int end, int stride);
      void rcv(int MASK, FLT *base,int bgn,int end, int stride);
};


