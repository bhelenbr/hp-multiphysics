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

spectral_hp& spectral_hp::operator=(const spectral_hp& tgt) {
   int i,n,j;
   
/*	COPY MESH INFORMATION */
   this->copy(tgt);
      
/*	SHALLOW COPY BASIS */   
   b = tgt.b;
   
   p0 = tgt.p0;
   sm0 = tgt.sm0;
   im0 = tgt.im0;
   
   if (size == 0) {
      vug = (FLT (*)[NV]) xmalloc(NV*maxvst*sizeof(FLT));
      sug = (FLT (*)[NV]) xmalloc(NV*maxvst*sm0*sizeof(FLT));
      iug = (FLT (*)[NV]) xmalloc(NV*maxvst*im0*sizeof(FLT));
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
         vug[i][n] = tgt.vug[i][n];
         
   for(i=0;i<nside*sm0;++i)
      for(n=0;n<NV;++n)
         sug[i][n] = tgt.sug[i][n];
         
   for(i=0;i<ntri*im0;++i)
      for(n=0;n<NV;++n)
         iug[i][n] = tgt.iug[i][n];
         
   for(i=0;i<nsbd;++i)
      for(j=0;j<sbdry[i].num*(1+sm0) +1;++j)
         binfo[i][j] = tgt.binfo[i][j];
         
   return(*this);
}
         
   


