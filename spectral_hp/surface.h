/*
 *  surface.h
 *  spectral_hp
 *
 *  Created by helenbrk on Tue Nov 06 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "hp_mgrid.h"

struct surface_glbls {
   int first;
   FLT drho;
   FLT rhoav;
   FLT fadd0;
   FLT fadd1;
   FLT (*vres)[ND];
   FLT (*sres)[ND];
   FLT (*vres0)[ND];
   FLT (*sres0)[ND];
   FLT (*vdt)[ND][ND];
   FLT (*sdt)[ND][ND];
}; 

class surface {
   public:
      FLT *ksprg;
      FLT (*vdres[MXLG2P])[ND];
      FLT (*sdres[MXLG2P])[ND];
      FLT *sfct;
      struct surface_glbls gbl;
      
   public:
      void init(int sind, int lg2p);
      friend class hg_mgrid;
};