/*
 *  refineby2.cpp
 *  mesh
 *
 *  Created by Brian Helenbrook on Tue Jun 11 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */
#include "mesh.h"
#include<utilities.h>
#include<assert.h>

void mesh::refineby2(const class mesh& inmesh) {
   int i,j,k,n,sind,tind,v0,v1,count,snum,vnear,err,initialsidenumber;
   FLT xpt[ND];
   
   /* INPUT MESH MUST HAVE GROWTH FACTOR OF 4 */
   /* BECAUSE OF INTWK USAGE */
    if (!initialized) {
      /* VERTEX STORAGE ALLOCATION */
      maxvst =  MAX(inmesh.maxvst,10);
      allocate(maxvst);
      nsbd = inmesh.nsbd;
      for(i=0;i<nsbd;++i) {
         getnewsideobject(i,inmesh.sbdry[i]->idnty());
         sbdry[i]->alloc(MAX(2*inmesh.sbdry[i]->mxsz(),10));
      }
      nvbd = inmesh.nvbd;
      for(i=0;i<nvbd;++i) {
         getnewvrtxobject(i,inmesh.vbdry[i]->idnty());
      }
      qtree.allocate(vrtx,maxvst);
      initialized = 1;
   }
   
   this->copy(inmesh);
   
   /* CALCULATE LOCATION OF NEW INTERIOR POINTS */
   count = nvrtx;
   for(sind=0;sind<nside;++sind) {
   
      if (stri[sind][1] < 0) continue;
      
      v0 = svrtx[sind][0];
      v1 = svrtx[sind][1];
      
      /* MIDPOINT */
      for(n=0;n<ND;++n)
         xpt[n] = 0.5*(vrtx[v0][n] +vrtx[v1][n]);
               
      /* INSERT POINT */
      for(n=0;n<ND;++n)
         vrtx[count][n] = xpt[n];
      ++count;
   }
   
   /*	INSERT INTERIOR POINTS */
   for(i=nvrtx;i<count;++i) { 
      qtree.addpt(nvrtx);
      qtree.nearpt(nvrtx,vnear);
      tind = findtri(vrtx[nvrtx],vnear);
      err = insert(tind,nvrtx,0.0);
      ++nvrtx;
   }

   
   /* INSERT BOUNDARY POINTS */
   for(i=0;i<nsbd;++i) {
      initialsidenumber = sbdry[i]->nsd();
      for(j=0;j<initialsidenumber;++j) {
         sind = sbdry[i]->sd(j);
         sbdry[i]->mvpttobdry(j,0.5,vrtx[nvrtx]);
         
         tind = stri[sind][0];
         snum = -2;
         for(k=0;k<3;++k) {
            if (tside[tind].side[k] == sind) {
               snum = k;
               break;
            }
         }
         assert(snum > -2);
         qtree.addpt(nvrtx);
         bdry_insert(tind,snum,nvrtx);
         ++nvrtx;
      }
   }
      
   for (i=0;i<nsbd;++i)
      sbdry[i]->reorder();
      
   cnt_nbor();
   
   return;
}

