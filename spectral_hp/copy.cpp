/*
 *  copy.cpp
 *  planar++
 *
 *  Created by helenbrk on Thu Oct 25 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include"spectral_hp.h"
#include"utilities.h"

void spectral_hp::copy(const spectral_hp& tgt) {
   int i,n,j;
   
/*	COPY MESH INFORMATION */
   this->mesh::copy(tgt);
      
/*	SHALLOW COPY BASIS */   
   b = tgt.b;
   
   p0 = tgt.p0;
   sm0 = tgt.sm0;
   im0 = tgt.im0;
   
   if (size == 0) {
      ug.v = (FLT (*)[NV]) xmalloc(NV*maxvst*sizeof(FLT));
      ug.s = (FLT (*)[NV]) xmalloc(NV*maxvst*sm0*sizeof(FLT));
      ug.i = (FLT (*)[NV]) xmalloc(NV*maxvst*im0*sizeof(FLT));
      size = maxvst;
      
/*	ALLOCATE STORAGE FOR BOUNDARIES */
      for(i=0;i<nsbd;++i)
         binfo[i] = new struct bistruct[maxsbel+1 +maxsbel*sm0];
   }
   else if (size < tgt.maxvst) {
      printf("spectral_hp object is not big enough for tgt\n");
      exit(1);
   }
      
   for(i=0;i<nvrtx;++i)
      for(n=0;n<NV;++n)
         ug.v[i][n] = tgt.ug.v[i][n];
         
   for(i=0;i<nside*sm0;++i)
      for(n=0;n<NV;++n)
         ug.s[i][n] = tgt.ug.s[i][n];
         
   for(i=0;i<ntri*im0;++i)
      for(n=0;n<NV;++n)
         ug.i[i][n] = tgt.ug.i[i][n];
         
   for(i=0;i<nsbd;++i)
      for(j=0;j<sbdry[i].num*(1+sm0) +1;++j)
         binfo[i][j] = tgt.binfo[i][j];
         
   return;
}
         
   


