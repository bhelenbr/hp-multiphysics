#include"blocks.h"
#include"utilities.h"

void block::init(int n, char *filename, FILETYPE filetype = easymesh, FLT grwfac = 1) {
   int i, maxvst;
   
   ngrid = n;
   grd = new class r_mesh[ngrid];
   
   grd[0].in_mesh(filename,filetype,grwfac);   
   grd[0].init_comm_buf(8);
   
/* WORK VARIABLES FOR MGRID */
   maxvst = grd[0].max();
   rglbl.work = (FLT (*)[ND]) xmalloc(maxvst*ND*sizeof(FLT));
   rglbl.res = (FLT (*)[ND]) xmalloc(maxvst*ND*sizeof(FLT));
   rglbl.diag = new FLT[maxvst];
   rglbl.fadd = 0.75;
   rglbl.vnn = 0.5;
   grd[0].init(0,&rglbl);
   
   for(i = 1; i< ngrid; ++i) {
      grd[i].coarsen(grd[i-1]);
      grd[i].init_comm_buf(8);
/*    grd[i].smooth_cofa(2); */
      grd[i].setfine(grd[i-1]);
      grd[i-1].setcoarse(grd[i]);
      grd[i].init(1,&rglbl);
   }

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