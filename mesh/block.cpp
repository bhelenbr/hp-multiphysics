#include"blocks.h"
#include"utilities.h"

void block::init(int n, char *filename, FILETYPE filetype = easymesh, FLT grwfac = 1) {
   int i;
   
   ngrid = n;
   grd = new class r_mesh[ngrid];
   
   grd[0].in_mesh(filename,filetype,grwfac);   
   grd[0].init_comm_buf(8);
   
   /* WORK VARIABLES FOR MGRID */
   r_mesh::fadd = 0.75;
   r_mesh::vnn = 0.5;
   grd[0].allocate(0,&rglbl);
   
   for(i = 1; i< ngrid; ++i) {
      grd[i].coarsen(1.6,grd[i-1]);
      grd[i].init_comm_buf(8);
      /* grd[i].smooth_cofa(2); */
      grd[i].setfine(grd[i-1]);
      grd[i-1].setcoarse(grd[i]);
      grd[i].allocate(1,&rglbl);
   }

   return;
}

void block::reconnect() {
   int i;
   
   for(i = 1; i< ngrid; ++i) {
      grd[i].coarsen(1.6,grd[i-1]);
      grd[i].setbcinfo();
      /* grd[i].smooth_cofa(2); */
      grd[i].setfine(grd[i-1]);
      grd[i-1].setcoarse(grd[i]);
   }

   return;
}

void block::coarsenchk(const char *fname) {
   int i;
   char name[100];

   for(i = 1; i< ngrid; ++i) {
      number_str(name,fname,i,1);
      grd[i].mesh::setbcinfo();
      grd[i].checkintegrity();
      grd[i].out_mesh(name);
      grd[i].setbcinfo();
   }
   return;
}