#include "mesh.h"
#include "boundary.h"

void mesh::coarsen_substructured(const class mesh &zx,int p) {
   int i,sind,sign,p2;

   copy(zx);
   
   p2 = p*p;
   ntri = zx.ntri/p2;
   nside = (ntri*3 +zx.sbdry[0]->nel/p)/2;
   nvrtx = zx.nvrtx -nside*(p-1) -ntri*((p-1)*(p-2))/2;
   
   for(i=0;i<zx.ntri;i+=p2) {
      td[i/p2].vrtx[0] = zx.td[i+2*p-2].vrtx[2];
      td[i/p2].vrtx[1] = zx.td[i].vrtx[0];
      td[i/p2].vrtx[2] = zx.td[i+p2-1].vrtx[1];
      sind = (zx.td[i].vrtx[1] -nvrtx)/(p-1);
      sign = ((zx.td[i].vrtx[1] -nvrtx)%(p-1) == 0 ? 1 : -1);
      td[i/p2].side[0] = sind; 
      td[i/p2].sign[0] = sign;
      sd[sind].vrtx[(1-sign)/2] = td[i/p2].vrtx[1];
      sd[sind].vrtx[(1+sign)/2] = td[i/p2].vrtx[2];
      
      sind = (zx.td[i+p2-1].vrtx[2] -nvrtx)/(p-1);
      sign = ((zx.td[i+p2-1].vrtx[2] -nvrtx)%(p-1) == 0 ? 1 : -1);
      td[i/p2].side[1] = sind;
      td[i/p2].sign[1] = sign;
      sd[sind].vrtx[(1-sign)/2] = td[i/p2].vrtx[2];
      sd[sind].vrtx[(1+sign)/2] = td[i/p2].vrtx[0];
      
      sind = (zx.td[i+2*p-2].vrtx[0] -nvrtx)/(p-1);
      sign = ((zx.td[i+2*p-2].vrtx[0] -nvrtx)%(p-1) == 0 ? 1 : -1);
      td[i/p2].side[2] = sind;
      td[i/p2].sign[2] = sign;
      sd[sind].vrtx[(1-sign)/2] = td[i/p2].vrtx[0];
      sd[sind].vrtx[(1+sign)/2] = td[i/p2].vrtx[1];
   }
   sbdry[0]->nel = zx.sbdry[0]->nel/p;
   i = 0;
   if (zx.sd[zx.sbdry[0]->el[0]].vrtx[0] < nvrtx) ++i;
   for(;i<zx.sbdry[0]->nel;i+=p) {
      sind = (zx.sd[zx.sbdry[0]->el[i]].vrtx[0] -nvrtx)/(p-1);
      sbdry[0]->el[i/p] = sind;
   }
   
   return;
}

void mesh::shift(FLT s[ND]) {
   int n;
   
   for(int i=0;i<nvrtx;++i)
      for(n=0;n<ND;++n)
         vrtx[i][1] += s[n];

   return;
}

void mesh::scale(FLT s[ND]) {
   int n;
   for(int i=0;i<nvrtx;++i)
      for(n=0;n<ND;++n)
         vrtx[i][1] *= s[n];

   return;
}