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
   FLT (*res)[ND];
   FLT fadd0;
   FLT fadd1;
   FLT (*res0)[ND];
}; 

class surface {
   public:
      FLT *ksprg;
      FLT (*dres[MXLG2P])[ND];
      FLT *sfct;
      struct surface_glbls gbl;
      
   public:
      void init(int sind, int lg2p);
      friend class hg_mgrid;
};