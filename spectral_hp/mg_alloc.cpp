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

void hp_mgrid::gbl_alloc(int mx, int p, struct hp_mgrid_glbls& store) {
   int sm, im;
   
   sm = p-1;
   im = (p-2)*(p-1)/2;

/*	SOLUTION STORAGE ON FIRST ENTRY TO NSTAGE */
   store.vug0 = (FLT (*)[NV]) xmalloc(NV*mx*sizeof(FLT));
   store.sug0 = (FLT (*)[NV]) xmalloc(NV*mx*sm*sizeof(FLT));
   store.iug0 = (FLT (*)[NV]) xmalloc(NV*mx*im*sizeof(FLT));

/*	RESIDUAL STORAGE */
   store.vres = (FLT (*)[NV]) xmalloc(NV*mx*sizeof(FLT));
   store.sres = (FLT (*)[NV]) xmalloc(NV*mx*sm*sizeof(FLT));
   store.ires = (FLT (*)[NV]) xmalloc(NV*mx*im*sizeof(FLT));

/*	VISCOUS FORCE RESIDUAL STORAGE */
   store.vvf = (FLT (*)[NV]) xmalloc(NV*mx*sizeof(FLT));
   store.svf = (FLT (*)[NV]) xmalloc(NV*mx*sm*sizeof(FLT));
   store.ivf = (FLT (*)[NV]) xmalloc(NV*mx*im*sizeof(FLT));
   
/*	RESIDUAL STORAGE FOR ENTRY TO MULTIGRID */
   store.vres0 = (FLT (*)[NV]) xmalloc(NV*mx*sizeof(FLT));
   p = p>>1;
   sm = p-1;
   im = (p-2)*(p-1)/2;
   if (sm > 0) store.sres0 = (FLT (*)[NV]) xmalloc(NV*mx*sm*sizeof(FLT));
   if (im > 0) store.ires0 = (FLT (*)[NV]) xmalloc(NV*mx*im*sizeof(FLT));

/*	PRECONDITIONER  */
   vect_alloc(store.gam,mx,FLT);
   vect_alloc(store.dtstar,mx,FLT);
   vect_alloc(store.vdiagv,mx,FLT);
   vect_alloc(store.vdiagp,mx,FLT);
   vect_alloc(store.sdiagv,mx,FLT);
   vect_alloc(store.sdiagp,mx,FLT);
   
/* STABILIZATION */
   vect_alloc(store.tau,mx,FLT);
   vect_alloc(store.delt,mx,FLT);

#ifdef SKIP
/*	UNSTEADY SOURCE TERMS (NEEDED ON FINE MESH ONLY) */
   FLT (*ug0)[NV], (*ug1)[NV], (*ug2)[NV], *jcb1, *jcb2;
   FLT ***dudt[ND], ***cdjdt;
#endif

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
   
