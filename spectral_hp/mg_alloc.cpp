/*
 *  mg_alloc.cpp
 *  planar++
 *
 *  Created by helenbrk on Tue Oct 30 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include"hp_mgrid.h"
#include<utilities.h>

/*	STATIC WORK VARIABLES USED BY ALL HP_MGRID OBJECTS */
FLT **hp_mgrid::cv00,**hp_mgrid::cv01,**hp_mgrid::cv10,**hp_mgrid::cv11;
FLT **hp_mgrid::e00,**hp_mgrid::e01,**hp_mgrid::e10,**hp_mgrid::e11;	
int hp_mgrid::size;

/*	STATIC VARIABLES USED BY ALL HP_MGRID OBJECTS */
const FLT hp_mgrid::alpha[NSTAGE+1];
const FLT hp_mgrid::beta[NSTAGE+1];
FLT hp_mgrid::dt0=0.0, hp_mgrid::dt1=0.0, hp_mgrid::dt2=0.0, hp_mgrid::dt3=0.0;


void hp_mgrid::allocate(struct hp_mgrid_glbls& ginit, int mgrid) {
   
   if (spectral_hp::size == 0 or mesh::initialized == 0) {
      printf("must initialize mesh/spectral_hp first\n");
      exit(1);
   }
   
   if (size == 0) {
      mat_alloc(cv00,b.gpx,b.gpn,FLT);
      mat_alloc(cv01,b.gpx,b.gpn,FLT);
      mat_alloc(cv10,b.gpx,b.gpn,FLT);
      mat_alloc(cv11,b.gpx,b.gpn,FLT);
      	
      mat_alloc(e00,b.gpx,b.gpn,FLT);
      mat_alloc(e01,b.gpx,b.gpn,FLT);
      mat_alloc(e10,b.gpx,b.gpn,FLT);
      mat_alloc(e11,b.gpx,b.gpn,FLT);
      size = b.p;
   }
   else {
      if (size < b.p) {
         printf("allocate from largest too smallest\n");
         exit(1);
      }
   }

/*	INITIALIZE GLOBAL STRUCTURE */   
   gbl = ginit;
   
/*	THINGS NEEDED FOR EACH MGRID LEVEL BUT NOT FINEST */
   if (mgrid) {
      log2p = 0;
      while ((b.p-1)>>log2p > 0) ++log2p;
      
      if (log2p > MXLG2P) {
         printf("make MXLG2P bigger %d %d\n",log2p,MXLG2P);
         exit(1);
      }

/*		NEED TO STORE INITIAL RESIDUAL ON EACH COARSE LEVEL */
      vdres[log2p] = (FLT (*)[NV]) xmalloc(NV*maxvst*sizeof(FLT));
      sdres[log2p] = (FLT (*)[NV]) xmalloc(NV*maxvst*b.sm*sizeof(FLT));
      idres[log2p] = (FLT (*)[NV]) xmalloc(NV*maxvst*b.im*sizeof(FLT));
         
/*		THINGS NEEDED ONLY FOR COARSE MESHES (NOT ON COARSE P BUT FINE MESH LEVELS) */
      if (p0 == 1) {
         vug_frst = (FLT (*)[NV]) xmalloc(NV*maxvst*sizeof(FLT));
      }
   }
   
   return;
}

void hp_mgrid::maxres(FLT mx[NV]) {
   int i,n;
   
   for(i=0;i<nvrtx;++i)
      for(n=0;n<NV;++n)
         mx[n] = MAX(mx[n],fabs(gbl.vres[i][n]));

   for(i=0;i<nside*b.sm;++i)
      for(n=0;n<NV;++n)
         mx[n] = MAX(mx[n],fabs(gbl.sres[i][n]));

   for(i=0;i<ntri*b.im;++i)
      for(n=0;n<NV;++n)
         mx[n] = MAX(mx[n],fabs(gbl.ires[i][n]));
         
   return;
}
   
