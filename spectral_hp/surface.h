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
/*	FLAG TO MARK DUPLICATE INTERFACES */
   int first;
/*	SOLUTION STORAGE ON FIRST ENTRY TO NSTAGE */
   FLT (*vug0)[ND];
   FLT (*sug0)[ND];
/*	AVERAGE PROPERTIES */
   FLT drho;
   FLT rhoav;
   FLT muav;
/*	ITERATIVE CONSTANTS */
   FLT fadd[ND];
   FLT cfl[MXLG2P][ND];
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
      struct surface_glbls gbl;
      
/*		FINE MESH GLBL ALLOCATION */
      void gbl_alloc(int maxside, int p, struct surface_glbls& store);

   public:
      void alloc(int maxside, int log2p, int mgrid, int fmesh, struct surface_glbls& store);
      friend class hp_mgrid;
};