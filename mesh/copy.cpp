#include"mesh.h"
#include"utilities.h"
#include<assert.h>


void mesh::copy(const mesh& tgt) {
   /* static */int i,j,n;
      
   if (!initialized) {
      allocate(tgt.maxvst);
      bdryalloc(tgt.maxsbel);
      initialized = 1;
   }
   else {
/*		CHECK IF BIG ENOUGH */
      if (tgt.nside > maxvst) {
         printf("mesh is too big to copy\n");
         exit(1);
      }
   }
   
/*	COPY VERTEX INFO OVER */
   nvrtx = tgt.nvrtx;
   for(i=0;i<nvrtx;++i)
      for(n=0;n<ND;++n)
         vrtx[i][n] = tgt.vrtx[i][n];

   for(i=0;i<nvrtx;++i)
      vlngth[i] = tgt.vlngth[i];

   for(i=0;i<nvrtx;++i)
      vinfo[i] = tgt.vinfo[i];

   for(i=0;i<nvrtx;++i)
      nnbor[i] = tgt.nnbor[i];
      
   for(i=0;i<nvrtx;++i)
      vtri[i] = tgt.vtri[i];
      
/*	tgt VERTEX BOUNDARY INFO */
   nvbd = tgt.nvbd;
   for(i=0;i<nvbd;++i) {
      vbdry[i].type = tgt.vbdry[i].type;
      vbdry[i].num = tgt.vbdry[i].num;
      for(j=0;j<vbdry[i].num;++j) 
         vbdry[i].el[j] = tgt.vbdry[i].el[j];
   }
   
/* tgt SIDE INFORMATION */
   nside = tgt.nside;
   for(i=0;i<nside;++i)
      for(n=0;n<2;++n)
         svrtx[i][n] = tgt.svrtx[i][n];

   for(i=0;i<nside;++i)
      for(n=0;n<2;++n)
         stri[i][n] = tgt.stri[i][n];
         
   for(i=0;i<nside;++i)
      sinfo[i] = tgt.sinfo[i];
      
/*	tgt SIDE BOUNDARY INFO */
   nsbd = tgt.nsbd;
   for(i=0;i<nsbd;++i) {
      sbdry[i].type = tgt.sbdry[i].type;
      sbdry[i].num = tgt.sbdry[i].num;
      if (sbdry[i].num > maxsbel) {
         printf("too many side elements to copy\n");
         exit(1);
      }
      for(j=0;j<sbdry[i].num;++j) 
         sbdry[i].el[j] = tgt.sbdry[i].el[j];
   }
   
/*	tgt ELEMENT DATA */
   ntri = tgt.ntri;
   for(i=0;i<ntri;++i)
      for(n=0;n<3;++n)
         tvrtx[i][n] = tgt.tvrtx[i][n];
 
   for(i=0;i<ntri;++i)
      for(n=0;n<3;++n)
         ttri[i][n] = tgt.ttri[i][n];
   
   for(i=0;i<ntri;++i) {
      for(n=0;n<3;++n) {
         tside[i].side[n] = tgt.tside[i].side[n];
         tside[i].sign[n] = tgt.tside[i].sign[n];
      }
   }
   
   for(i=0;i<ntri;++i)
      tinfo[i] = tgt.tinfo[i];
      
   qtree.copy(tgt.qtree);
   qtree.change_vptr(vrtx);
   
   return;  
}

