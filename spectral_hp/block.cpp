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
   int i;
   
   hpbase = bin;
   lg2pmax = lg2p;
   
/*	INITIALIZE SPECTRAL_HP'S */
   grd[0].spectral_hp::allocate(hpbase[lg2pmax]);
   grd[0].init_comm_buf(NV*(hpbase[lg2pmax].sm +2));
   for(i=1;i<ngrid;++i) {
      grd[i].spectral_hp::allocate(hpbase[0]);
      grd[i].init_comm_buf(NV*2);
   }

/*	INITIALIZE HP_MGRID STORAGE NECESSARY FOR EACH MESH */ 
	grd[0].allocate(0,gbl);  // 0 DENOTES FINEST LEVEL ALLOCATES GLOBAL STORAGE
   for(i = lg2pmax -1; i >= 0; --i) {
      grd[0].loadbasis(hpbase[i]);
      grd[0].allocate(1,gbl); // 1 DENOTES MGRID LEVEL
   }
   for(i=1;i<ngrid;++i)
      grd[i].allocate(1,gbl);
   
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


