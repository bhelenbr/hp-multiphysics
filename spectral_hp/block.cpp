#include"blocks.h"
#include"utilities.h"

void block::meshinit(int n, char *filename, FILETYPE filetype = easymesh, FLT grwfac = 1) {
   int i;
   
   ngrid = n;
   grd = new class hp_mgrid[ngrid];
   grd[0].in_mesh(filename,filetype,grwfac);
   for(i = 1; i< ngrid; ++i) {
      grd[i].coarsen(grd[i-1]);
/*    grd[i].smooth_cofa(2); */
      grd[i].setfine(grd[i-1]);
      grd[i-1].setcoarse(grd[i]);
   }

   return;
}

void block::hpinit(class hpbasis *bin, int lg2p) {
   int maxvst,i;
   
   hpbase = bin;
   lg2pmax = lg2p;

/*	INITIALIZE SPECTRAL_HP'S */
   grd[0].spectral_hp::allocate(hpbase[lg2pmax]);
   for(i=1;i<ngrid;++i)
      grd[i].spectral_hp::allocate(hpbase[0]);
      

/*	ALLOCATE GLOBAL STORAGE */
   maxvst = grd[0].max();
/*	SOLUTION STORAGE ON FIRST ENTRY TO NSTAGE */
   gbl.vug0 = (FLT (*)[NV]) xmalloc(NV*maxvst*sizeof(FLT));
   gbl.sug0 = (FLT (*)[NV]) xmalloc(NV*maxvst*hpbase[lg2pmax].sm*sizeof(FLT));
   gbl.iug0 = (FLT (*)[NV]) xmalloc(NV*maxvst*hpbase[lg2pmax].im*sizeof(FLT));

/*	RESIDUAL STORAGE */
   gbl.vres = (FLT (*)[NV]) xmalloc(NV*maxvst*sizeof(FLT));
   gbl.sres = (FLT (*)[NV]) xmalloc(NV*maxvst*hpbase[lg2pmax].sm*sizeof(FLT));
   gbl.ires = (FLT (*)[NV]) xmalloc(NV*maxvst*hpbase[lg2pmax].im*sizeof(FLT));

/*	VISCOUS FORCE RESIDUAL STORAGE */
   gbl.vvf = (FLT (*)[NV]) xmalloc(NV*maxvst*sizeof(FLT));
   gbl.svf = (FLT (*)[NV]) xmalloc(NV*maxvst*hpbase[lg2pmax].sm*sizeof(FLT));
   gbl.ivf = (FLT (*)[NV]) xmalloc(NV*maxvst*hpbase[lg2pmax].im*sizeof(FLT));
   
/*	RESIDUAL STORAGE FOR ENTRY TO MULTIGRID */
   gbl.vres0 = (FLT (*)[NV]) xmalloc(NV*maxvst*sizeof(FLT));
   if (lg2pmax > 1) gbl.sres0 = (FLT (*)[NV]) xmalloc(NV*maxvst*hpbase[lg2pmax-1].sm*sizeof(FLT));
   if (lg2pmax > 1) gbl.ires0 = (FLT (*)[NV]) xmalloc(NV*maxvst*hpbase[lg2pmax-1].im*sizeof(FLT));

/*	PRECONDITIONER  */
   vect_alloc(gbl.gam,maxvst,FLT);
   vect_alloc(gbl.dtstar,maxvst,FLT);
   vect_alloc(gbl.vdiagv,maxvst,FLT);
   vect_alloc(gbl.vdiagp,maxvst,FLT);
   vect_alloc(gbl.sdiagv,maxvst,FLT);
   vect_alloc(gbl.sdiagp,maxvst,FLT);
   
/* STABILIZATION */
   vect_alloc(gbl.tau,maxvst,FLT);
   vect_alloc(gbl.delt,maxvst,FLT);

#ifdef SKIP
/*	UNSTEADY SOURCE TERMS (NEEDED ON FINE MESH ONLY) */
   FLT (*ug0)[NV], (*ug1)[NV], (*ug2)[NV], *jcb1, *jcb2;
   FLT ***dudt[ND], ***cdjdt;
#endif

/*	INITIALIZE HP_MGRID STORAGE NECESSARY FOR EACH MESH */ 
	grd[0].allocate(gbl,0);  // 0 DENOTES FINEST LEVEL
   for(i = lg2pmax -1; i >= 0; --i) {
      grd[0].loadbasis(hpbase[i]);
      grd[0].allocate(gbl,1); // 1 DENOTES MGRID LEVEL
   }
   for(i=1;i<ngrid;++i)
      grd[i].allocate(gbl,1);
   
   return;
}


void block::reconnect() {
   int i;
   
   for(i = 1; i< ngrid; ++i) {
      grd[i].coarsen(grd[i-1]);
/*    grd[i].smooth_cofa(2); */
      grd[i].setfine(grd[i-1]);
      grd[i-1].setcoarse(grd[i]);
   }

   return;
}


