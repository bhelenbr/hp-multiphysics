/*
 *  allocate.cpp
 *  planar++
 *
 *  Created by helenbrk on Thu Oct 11 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "spectral_hp.h"
#include<utilities.h>

/* STATIC WORK ARRAYS */
FLT **spectral_hp::u[NV],**spectral_hp::du[NV][ND], **spectral_hp::res[NV];
FLT **spectral_hp::crd[ND],**spectral_hp::dcrd[ND][ND], **spectral_hp::cjcb;
FLT *spectral_hp::uht[NV], *spectral_hp::lf[NV], *spectral_hp::lf1[NV];
int spectral_hp::pmax = 0;
   
void spectral_hp::allocate(class hpbasis& bas) {
   int i,m,n;
   
   if (initialized == 0) {
      printf("must load mesh before allocating\n");
      exit(1);
   }

/*	SHALLOW COPY BASIS */   
   b = bas;
   
   p0 = b.p;
   sm0 = b.sm;
   im0 = b.im;
      
   /*	LOCAL STORAGE/WORK */
   if (p0 > pmax) {
      if (pmax != 0) {
         printf("allocate from largest to smallest\n");
         exit(1);
      }
      for(n=0;n<NV;++n) {
         mat_alloc(u[n],b.gpx,b.gpn,FLT);
         mat_alloc(res[n],b.gpx,b.gpn,FLT);
         for(m=0;m<ND;++m)
            mat_alloc(du[n][m],b.gpx,b.gpn,FLT);
            
         vect_alloc(uht[n],b.tm,FLT);
         vect_alloc(lf[n],b.tm,FLT);
         vect_alloc(lf1[n],b.tm,FLT);
      }
      
      for(n=0;n<ND;++n) {
         mat_alloc(crd[n],b.gpx,b.gpn,FLT);
         for(m=0;m<ND;++m)
            mat_alloc(dcrd[n][m],b.gpx,b.gpn,FLT);
      }
      mat_alloc(cjcb,b.gpx,b.gpn,FLT);
      pmax = b.p;
   }
   
   ug.v = (FLT (*)[NV]) xmalloc(NV*maxvst*sizeof(FLT));
   ug.s = (FLT (*)[NV]) xmalloc(NV*maxvst*sm0*sizeof(FLT));
   ug.i = (FLT (*)[NV]) xmalloc(NV*maxvst*im0*sizeof(FLT));
   
   size = maxvst;
   
/*	ALLOCATE STORAGE FOR BOUNDARIES */
   for(i=0;i<nsbd;++i)
      binfo[i] = new struct bistruct[maxsbel+1 +maxsbel*sm0];
      
   setbcinfo();

   return;
}

void spectral_hp::setbcinfo() {
   static int i,j,sind;
   
/*	SET UP VRTX BC INFORMATION FOR OUTPUT */
   for(i=0;i<nvrtx;++i)
      vinfo[i] = -1;

   for(i=0;i<nvbd;++i)
      for(j=0;j<vbdry[i].num;++j)
         vinfo[vbdry[i].el[j]] = vbdry[i].type;

/*	SET UP SIDE BC INFORMATION FOR CURVED SIDES OUTPUT */
   for(i=0;i<nside;++i)
      sinfo[i] = -1;
   
   for(i=0;i<ntri;++i)
      tinfo[i] = -1;
   
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type&CURV_MASK) {
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            sinfo[sind] = 0;
            tinfo[stri[sind][0]] = 0;
         }
      }
   }
   
   return;
}
   
   
            

      
      

 
   

