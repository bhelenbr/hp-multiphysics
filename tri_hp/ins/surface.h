/*
 *  surface.h
 *  tri_hp
 *
 *  Created by helenbrk on Tue Nov 06 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */
#define ND 2
#define MXLG2P 5

struct surface_glbls {
   /* SOLUTION STORAGE ON FIRST ENTRY TO NSTAGE */
   FLT (*vug0)[ND];
   FLT (*sug0)[ND];
   /* AVERAGE PROPERTIES */
   FLT sigma;
   FLT rho2;
   FLT mu2;
   /* RESIDUALS */
   FLT (*vres)[ND];
   FLT (*sres)[ND];
   FLT (*vres0)[ND];
   FLT (*sres0)[ND];
   /* PSEUDO TIME ITERATION */
   FLT (*vdt)[ND][ND];
   FLT (*sdt)[ND][ND];
   FLT *dtfnrm;
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
      
      /* THINGS USED BY ALL SURFACES */
      static FLT fadd[ND];
      static FLT cfl[MXLG2P][ND];
      
      /* FINE MESH GLBL ALLOCATION */
      void gbl_alloc(int maxside, int p, struct surface_glbls *store);

   public:
      void alloc(int maxside, int log2p, int mgrid, int fmesh, struct surface_glbls *store);
      /* COUPLED SURFACE BOUNDARY ROUTINES */ 
      void surfvrttoug();
      void surfugtovrt1();
      void surfugtovrt2();
      /* SETUP SURFACE 1D SPRING CONSTANTS */
      void setksprg1d();
      void surfksrc1d();
      void surfrsdl(int mgrid);
      void surfinvrt1(int bnum);
      void surfinvrt2(int bnum);
      void surfdt1(int bnum);
      void surfdt2(int bnum);
      void surfnstage1(int bnum);
      void surfnstage2(int bnum, int stage);
      void surfgetfres(int bnum);
      void surfgetcchng(int bnum);
      void surfmaxres();
};

