/*
 *  surface.h
 *  spectral_hp
 *
 *  Created by helenbrk on Tue Nov 06 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */
#define ND 2
#define MXLG2P 5

struct surface_glbls {
/*	SOLUTION STORAGE ON FIRST ENTRY TO NSTAGE */
   FLT (*vug0)[ND];
   FLT (*sug0)[ND];
/*	AVERAGE PROPERTIES */
   FLT sigma;
   FLT rho2;
   FLT mu2;
/*	RESIDUALS */
   FLT (*vres)[ND];
   FLT (*sres)[ND];
   FLT (*vres0)[ND];
   FLT (*sres0)[ND];
/*	PSEUDO TIME ITERATION */
   FLT (*vdt)[ND][ND];
   FLT (*sdt)[ND][ND];
   FLT *normc;
   FLT *meshc;
}; 

class surface {
   protected:
      FLT *ksprg;
      FLT (*vug)[ND];
      FLT (*sug)[ND];
      FLT (*vug_frst)[ND];
      FLT (*vdres[MXLG2P])[ND];
      FLT (*sdres[MXLG2P])[ND];
      struct surface_glbls *gbl;
      
/*		THINGS USED BY ALL SURFACES */
      static FLT fadd[ND];
      static FLT cfl[MXLG2P][ND];
      
/*		FINE MESH GLBL ALLOCATION */
      void gbl_alloc(int maxside, int p, struct surface_glbls *store);

   public:
      void alloc(int maxside, int log2p, int mgrid, int fmesh, struct surface_glbls *store);
      friend class hp_mgrid;
      friend class block;
      friend class blocks;
};