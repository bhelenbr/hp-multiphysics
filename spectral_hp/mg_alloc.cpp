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

/* STATIC WORK VARIABLES USED BY ALL HP_MGRID OBJECTS */
FLT **hp_mgrid::cv00,**hp_mgrid::cv01,**hp_mgrid::cv10,**hp_mgrid::cv11;
FLT **hp_mgrid::e00,**hp_mgrid::e01,**hp_mgrid::e10,**hp_mgrid::e11; 
FLT **(hp_mgrid::mvel[ND]); // for local mesh velocity info
FLT hp_mgrid::fadd, hp_mgrid::cfl[MXLG2P];   // ITERATION PARAMETERS  
FLT hp_mgrid::adis; // STABILIZATION
int hp_mgrid::charyes;  // USE CHARACTERISTIC FAR-FIELD B.C'S
FLT hp_mgrid::trncerr, hp_mgrid::invbdryerr, hp_mgrid::vlngth_tol, hp_mgrid::adapt_tol;
int hp_mgrid::changed = 1; //FLAG FOR PV3 TO INDICATE STRUCTURE CHANGED
struct vsi hp_mgrid::ugwk[TMADAPT]; // STORAGE FOR UNSTEADY ADAPTATION BD FLOW INFO
FLT (*hp_mgrid::vrtxwk[TMADAPT])[ND]; // STORAGE FOR UNSTEADY ADAPTATION MESH BD INFO
struct bistruct *hp_mgrid::binfowk[TMADAPT][MAXSB]; // STORAGE FOR UNSTEADY ADAPTATION BOUNDARY BD INFO
FLT **hp_mgrid::bdwk[TMADAPT][NV]; // WORK FOR ADAPTATION
int hp_mgrid::size;


/* STATIC VARIABLES USED BY ALL HP_MGRID OBJECTS */
const FLT hp_mgrid::alpha[NSTAGE+1] = {0.25, 1./6., .375, .5, 1.0, 1.0};
const FLT hp_mgrid::beta[NSTAGE+1] = {1.0, 0.0, 5./9., 0.0, 4./9., 1.0};
//const FLT hp_mgrid::alpha[NSTAGE+1] = {1.0, 1.0};  // 1 JACOBI SWEEP
//const FLT hp_mgrid::beta[NSTAGE+1] = {1.0, 1.0};  // 1 JACOBI SWEEP
FLT hp_mgrid::dti=0.0, hp_mgrid::time=0.0, hp_mgrid::g=0.0;


void hp_mgrid::allocate(int mgrid, struct hp_mgrid_glbls *store) {
   int i,j,n,onfmesh;
   
#if (defined(DROP) || defined(UNSTEADY_DROP))
   r_mesh::fixy_mask = 0xff -SYMM_MASK -OUTF_MASK -PRDX_MASK -INFL_MASK;
#else
   r_mesh::fixy_mask = 0xff -SYMM_MASK -OUTF_MASK -PRDX_MASK;
#endif
   
   if (spectral_hp::size == 0 || mesh::initialized == 0) {
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
      
      for(n=0;n<ND;++n)
         mat_alloc(mvel[n],b.gpx,b.gpn,FLT);
      
      /* ALLOCATE UNSTEADY ADAPTATION STORAGE */
      for(i=0;i<TMADAPT;++i) {
         ugwk[i].v = (FLT (*)[NV]) xmalloc(NV*maxvst*sizeof(FLT));
         ugwk[i].s = (FLT (*)[NV]) xmalloc(NV*maxvst*b.sm*sizeof(FLT));
         ugwk[i].i = (FLT (*)[NV]) xmalloc(NV*maxvst*b.im*sizeof(FLT));
         
         vrtxwk[i] = (FLT (*)[ND]) xmalloc(ND*maxvst*sizeof(FLT));
         for(j=0;j<nsbd;++j)
            binfowk[i][j] = new struct bistruct[maxsbel+1 +maxsbel*sm0];
               
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

   /* ON FINEST MESH ALLOCATE GLOBAL STORAGE */   
   if (!mgrid) gbl_alloc(store);

   /* COPY GLOBAL POINTER */   
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
   
   /* THINGS NEEDED 1 FOR EACH PHYSICAL GRID */
   if ((onfmesh && b.p == p0) || !onfmesh)
      dvrtdt = (FLT (*)[ND]) xmalloc(ND*maxvst*sizeof(FLT));

   /* THINGS NEEDED FOR EACH MGRID LEVEL BUT NOT FINEST */
   if (mgrid) {

      /* NEED TO STORE INITIAL RESIDUAL ON EACH COARSE LEVEL */
      dres[log2p].v = (FLT (*)[NV]) xmalloc(NV*maxvst*sizeof(FLT));
      dres[log2p].s = (FLT (*)[NV]) xmalloc(NV*maxvst*b.sm*sizeof(FLT));
      dres[log2p].i = (FLT (*)[NV]) xmalloc(NV*maxvst*b.im*sizeof(FLT));
         
      /* THINGS NEEDED ONLY FOR COARSE MESHES (NOT ON COARSE P LEVELS) */
      if (!onfmesh) {
         vug_frst = (FLT (*)[NV]) xmalloc(NV*maxvst*sizeof(FLT));
      }
   }
   
   return;
}

void hp_mgrid::gbl_alloc(struct hp_mgrid_glbls *store) {
   int i, j, pn, smn, imn;
   
   /* SOLUTION STORAGE ON FIRST ENTRY TO NSTAGE */
   store->ug0.v = (FLT (*)[NV]) xmalloc(NV*maxvst*sizeof(FLT));
   store->ug0.s = (FLT (*)[NV]) xmalloc(NV*maxvst*b.sm*sizeof(FLT));
   store->ug0.i = (FLT (*)[NV]) xmalloc(NV*maxvst*b.im*sizeof(FLT));

   /* RESIDUAL STORAGE */
   store->res.v = (FLT (*)[NV]) xmalloc(NV*maxvst*sizeof(FLT));
   store->res.s = (FLT (*)[NV]) xmalloc(NV*maxvst*b.sm*sizeof(FLT));
   store->res.i = (FLT (*)[NV]) xmalloc(NV*maxvst*b.im*sizeof(FLT));

   /* VISCOUS FORCE RESIDUAL STORAGE */
   store->vf.v = (FLT (*)[NV]) xmalloc(NV*maxvst*sizeof(FLT));
   store->vf.s = (FLT (*)[NV]) xmalloc(NV*maxvst*b.sm*sizeof(FLT));
   store->vf.i = (FLT (*)[NV]) xmalloc(NV*maxvst*b.im*sizeof(FLT));
   
   /* RESIDUAL STORAGE FOR ENTRY TO MULTIGRID NEXT COARSER MESH */
   store->res0.v = (FLT (*)[NV]) xmalloc(NV*maxvst*sizeof(FLT));
   pn = b.p>>1;
   smn = pn-1;
   imn = (pn-2)*(pn-1)/2;
   if (smn > 0) store->res0.s = (FLT (*)[NV]) xmalloc(NV*maxvst*smn*sizeof(FLT));
   if (imn > 0) store->res0.i = (FLT (*)[NV]) xmalloc(NV*maxvst*imn*sizeof(FLT));

   /* PRECONDITIONER  */
   store->vprcn = (FLT (*)[NV][NV]) xmalloc(NV*NV*maxvst*sizeof(FLT));
   store->sprcn = (FLT (*)[NV][NV]) xmalloc(NV*NV*maxvst*sizeof(FLT));
   store->tprcn = (FLT (*)[NV][NV]) xmalloc(NV*NV*maxvst*sizeof(FLT));
   
/* STABILIZATION */
   vect_alloc(store->tau,maxvst,FLT);
   vect_alloc(store->delt,maxvst,FLT);

   /* UNSTEADY SOURCE TERMS (NEEDED ON FINE MESH ONLY) */
   for(i=0;i<TMSTORE;++i) {
      store->ugbd[i].v = (FLT (*)[NV]) xmalloc(NV*maxvst*sizeof(FLT));
      store->ugbd[i].s = (FLT (*)[NV]) xmalloc(NV*maxvst*b.sm*sizeof(FLT));
      store->ugbd[i].i = (FLT (*)[NV]) xmalloc(NV*maxvst*b.im*sizeof(FLT));
      
      store->vrtxbd[i] = (FLT (*)[ND]) xmalloc(ND*maxvst*sizeof(FLT));
      for(j=0;j<nsbd;++j)
         store->binfobd[i][j] = new struct bistruct[maxsbel+1 +maxsbel*sm0];
   }
   for(i=0;i<NV;++i)
      tens_alloc(store->dugdt[i],maxvst,b.gpx,b.gpn,FLT);  // UNSTEADY SOURCE FOR FLOW
   
   for(j=0;j<nsbd;++j)
      store->dbinfodt[j] = new struct bistruct[maxsbel+1 +maxsbel*sm0];

   return;
}

FLT hp_mgrid::maxres(FLT *mxr) {
   int i,n;
   FLT emax = 0.0;
   
   for(n=0;n<NV;++n)
      mxr[n] = 0.0;
   
   for(i=0;i<nvrtx;++i)
      for(n=0;n<NV;++n)
         mxr[n] = MAX(mxr[n],fabs(gbl->res.v[i][n]));

   for(n=0;n<NV;++n)
      emax = MAX(emax,mxr[n]);
      
   return(emax);
}   
   

