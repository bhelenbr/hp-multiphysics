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
FLT hp_mgrid::fadd, hp_mgrid::cfl[MXLG2P];   // ITERATION PARAMETERS  
FLT hp_mgrid::adis; // STABILIZATION
int hp_mgrid::charyes;  // USE CHARACTERISTIC FAR-FIELD B.C'S
FLT hp_mgrid::trncerr, hp_mgrid::tol;
class hp_mgrid hp_mgrid::hpstr; // STORAGE FOR ADAPTATION 
struct vsi hp_mgrid::ugstr[MXSTEP-1]; // STORAGE FOR UNSTEADY ADAPTATION BD FLOW INFO
FLT (*hp_mgrid::vrtxstr[MXSTEP-1])[ND]; // STORAGE FOR UNSTEADY ADAPTATION MESH BD INFO
struct bistruct *hp_mgrid::binfostr[MXSTEP-1][MAXSB]; // STORAGE FOR UNSTEADY ADAPTATION BOUNDARY BD INFO
FLT **hp_mgrid::bdwk[MXSTEP-1][NV]; // WORK FOR ADAPTATION
int hp_mgrid::size;


/*	STATIC VARIABLES USED BY ALL HP_MGRID OBJECTS */
const FLT hp_mgrid::alpha[NSTAGE+1];
const FLT hp_mgrid::beta[NSTAGE+1];
int hp_mgrid::nstep=1;
FLT hp_mgrid::bd[MXSTEP+1];
FLT hp_mgrid::dti=0.0, hp_mgrid::time=0.0, hp_mgrid::g=0.0;


void hp_mgrid::allocate(int mgrid, struct hp_mgrid_glbls *store) {
   int i,j,n,onfmesh;
   
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
      
/*		ALLOCATE UNSTEADY ADAPTATION STORAGE */
      for(i=0;i<MXSTEP-1;++i) {
         ugstr[i].v = (FLT (*)[NV]) xmalloc(NV*maxvst*sizeof(FLT));
         ugstr[i].s = (FLT (*)[NV]) xmalloc(NV*maxvst*b.sm*sizeof(FLT));
         ugstr[i].i = (FLT (*)[NV]) xmalloc(NV*maxvst*b.im*sizeof(FLT));
         
         vrtxstr[i] = (FLT (*)[ND]) xmalloc(ND*maxvst*sizeof(FLT));
         for(j=0;j<nsbd;++j)
            if (sbdry[j].type&CURV_MASK) 
               binfostr[i][j] = new struct bistruct[maxsbel+1 +maxsbel*sm0];
               
         for(n=0;n<NV;++n)
            mat_alloc(bdwk[i][n],b.gpx,b.gpn,FLT);
      }      
   }
   else {
      if (size < b.p) {
         printf("allocate from largest to smallest\n");
         exit(1);
      }
   }

/*	ON FINEST MESH ALLOCATE GLOBAL STORAGE */   
   if (!mgrid) gbl_alloc(store);

/*	COPY GLOBAL POINTER */   
   gbl = store;

/* FIGURE OUT WHERE WE ARE */   
   if (mgrid && p0 == 1) onfmesh = 0;
   else onfmesh = 1;
   log2p = 0;
   while ((b.p-1)>>log2p > 0) ++log2p;
   if (log2p > MXLG2P) {
      printf("make MXLG2P bigger %d %d\n",log2p,MXLG2P);
      exit(1);
   }
   
/*	THINGS NEEDED 1 FOR EACH PHYSICAL GRID */
   if ((onfmesh && b.p == p0) || !onfmesh)
   	dvrtdt = (FLT (*)[ND]) xmalloc(ND*maxvst*sizeof(FLT));

/*	THINGS NEEDED FOR EACH MGRID LEVEL BUT NOT FINEST */
   if (mgrid) {

/*		NEED TO STORE INITIAL RESIDUAL ON EACH COARSE LEVEL */
      dres[log2p].v = (FLT (*)[NV]) xmalloc(NV*maxvst*sizeof(FLT));
      dres[log2p].s = (FLT (*)[NV]) xmalloc(NV*maxvst*b.sm*sizeof(FLT));
      dres[log2p].i = (FLT (*)[NV]) xmalloc(NV*maxvst*b.im*sizeof(FLT));
         
/*		THINGS NEEDED ONLY FOR COARSE MESHES (NOT ON COARSE P LEVELS) */
      if (!onfmesh) {
         vug_frst = (FLT (*)[NV]) xmalloc(NV*maxvst*sizeof(FLT));
      }
   }
   
   return;
}

void hp_mgrid::gbl_alloc(struct hp_mgrid_glbls *store) {
   int i, j, pn, smn, imn;
   
/*	SOLUTION STORAGE ON FIRST ENTRY TO NSTAGE */
   store->ug0.v = (FLT (*)[NV]) xmalloc(NV*maxvst*sizeof(FLT));
   store->ug0.s = (FLT (*)[NV]) xmalloc(NV*maxvst*b.sm*sizeof(FLT));
   store->ug0.i = (FLT (*)[NV]) xmalloc(NV*maxvst*b.im*sizeof(FLT));

/*	RESIDUAL STORAGE */
   store->res.v = (FLT (*)[NV]) xmalloc(NV*maxvst*sizeof(FLT));
   store->res.s = (FLT (*)[NV]) xmalloc(NV*maxvst*b.sm*sizeof(FLT));
   store->res.i = (FLT (*)[NV]) xmalloc(NV*maxvst*b.im*sizeof(FLT));

/*	VISCOUS FORCE RESIDUAL STORAGE */
   store->vf.v = (FLT (*)[NV]) xmalloc(NV*maxvst*sizeof(FLT));
   store->vf.s = (FLT (*)[NV]) xmalloc(NV*maxvst*b.sm*sizeof(FLT));
   store->vf.i = (FLT (*)[NV]) xmalloc(NV*maxvst*b.im*sizeof(FLT));
   
/*	RESIDUAL STORAGE FOR ENTRY TO MULTIGRID NEXT COARSER MESH */
   store->res0.v = (FLT (*)[NV]) xmalloc(NV*maxvst*sizeof(FLT));
   pn = b.p>>1;
   smn = pn-1;
   imn = (pn-2)*(pn-1)/2;
   if (smn > 0) store->res0.s = (FLT (*)[NV]) xmalloc(NV*maxvst*smn*sizeof(FLT));
   if (imn > 0) store->res0.i = (FLT (*)[NV]) xmalloc(NV*maxvst*imn*sizeof(FLT));

/*	PRECONDITIONER  */
   vect_alloc(store->gam,maxvst,FLT);
   vect_alloc(store->dtstar,maxvst,FLT);
   vect_alloc(store->vdiagv,maxvst,FLT);
   vect_alloc(store->vdiagp,maxvst,FLT);
   vect_alloc(store->sdiagv,maxvst,FLT);
   vect_alloc(store->sdiagp,maxvst,FLT);
   
/* STABILIZATION */
   vect_alloc(store->tau,maxvst,FLT);
   vect_alloc(store->delt,maxvst,FLT);

/*	UNSTEADY SOURCE TERMS (NEEDED ON FINE MESH ONLY) */
   for(i=0;i<MXSTEP-1;++i) {
      store->ugbd[i].v = (FLT (*)[NV]) xmalloc(NV*maxvst*sizeof(FLT));
      store->ugbd[i].s = (FLT (*)[NV]) xmalloc(NV*maxvst*b.sm*sizeof(FLT));
      store->ugbd[i].i = (FLT (*)[NV]) xmalloc(NV*maxvst*b.im*sizeof(FLT));
      
      store->vrtxbd[i] = (FLT (*)[ND]) xmalloc(ND*maxvst*sizeof(FLT));
      for(j=0;j<nsbd;++j)
         if (sbdry[j].type&CURV_MASK) 
            store->binfobd[i][j] = new struct bistruct[maxsbel+1 +maxsbel*sm0];
   }
   for(i=0;i<NV;++i)
      tens_alloc(store->dugdt[i],maxvst,b.gpx,b.gpn,FLT);  // UNSTEADY SOURCE FOR FLOW
   
   for(j=0;j<nsbd;++j)
      if (sbdry[j].type&CURV_MASK) 
         store->dbinfodt[j] = new struct bistruct[maxsbel+1 +maxsbel*sm0];

   return;
}

void hp_mgrid::maxres(FLT mx[NV]) {
   int i,n;
   
   for(i=0;i<nvrtx;++i)
      for(n=0;n<NV;++n)
         mx[n] = MAX(mx[n],fabs(gbl->res.v[i][n]));

   for(i=0;i<nside*b.sm;++i)
      for(n=0;n<NV;++n)
         mx[n] = MAX(mx[n],fabs(gbl->res.s[i][n]));

   for(i=0;i<ntri*b.im;++i)
      for(n=0;n<NV;++n)
         mx[n] = MAX(mx[n],fabs(gbl->res.i[i][n]));
         
   return;
}
   
