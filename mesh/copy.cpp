#include "mesh.h"
#include "boundary.h"
#include "utilities.h"
#include<assert.h>

void mesh::copy(const mesh& tgt) {
   int i,n;
      
   if (!initialized) {
      allocate_duplicate(1.0,tgt);
   }
   else {
      /* CHECK IF BIG ENOUGH */
      if (tgt.nside > maxvst) {
         *log << "mesh is too big to copy" << std::endl;
         exit(1);
      }
   }
   
   log = tgt.log;
   
   /* COPY VERTEX INFO OVER */
   nvrtx = tgt.nvrtx;
   for(i=0;i<nvrtx;++i)
      for(n=0;n<ND;++n)
         vrtx[i][n] = tgt.vrtx[i][n];

   for(i=0;i<nvrtx;++i)
      vlngth[i] = tgt.vlngth[i];

   for(i=0;i<nvrtx;++i)
      vd[i].info = tgt.vd[i].info;

   for(i=0;i<nvrtx;++i)
      vd[i].nnbor = tgt.vd[i].nnbor;
      
   for(i=0;i<nvrtx;++i)
      vd[i].tri = tgt.vd[i].tri;
      
   /* COPY VERTEX BOUNDARY INFO */
   for(i=0;i<nvbd;++i)
      vbdry[i]->copy(*tgt.vbdry[i]);
         
/* COPY SIDE INFORMATION */
   nside = tgt.nside;
   for(i=0;i<nside;++i)
      for(n=0;n<2;++n)
         sd[i].vrtx[n] = tgt.sd[i].vrtx[n];

   for(i=0;i<nside;++i)
      for(n=0;n<2;++n)
         sd[i].tri[n] = tgt.sd[i].tri[n];
         
   for(i=0;i<nside;++i)
      sd[i].info = tgt.sd[i].info;
      
   /* COPY SIDE BOUNDARY INFO */
   for(i=0;i<nsbd;++i)
      sbdry[i]->copy(*tgt.sbdry[i]);
   
   /* COPY ELEMENT DATA */
   ntri = tgt.ntri;
   for(i=0;i<ntri;++i)
      for(n=0;n<3;++n)
         td[i].vrtx[n] = tgt.td[i].vrtx[n];
 
   for(i=0;i<ntri;++i)
      for(n=0;n<3;++n)
         td[i].tri[n] = tgt.td[i].tri[n];
   
   for(i=0;i<ntri;++i) {
      for(n=0;n<3;++n) {
         td[i].side[n] = tgt.td[i].side[n];
         td[i].sign[n] = tgt.td[i].sign[n];
      }
   }
   
   for(i=0;i<ntri;++i)
      td[i].info = tgt.td[i].info;
      
   qtree.copy(tgt.qtree);
   qtree.change_vptr(vrtx);
   
   return;  
}

