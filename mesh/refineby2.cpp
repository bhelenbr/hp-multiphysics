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
   int i,j,k,sind,tind,v0,v1,count,bid,snum,vnear,err,initialsidenumber;
   FLT xpt,ypt;
   
   /* INPUT MESH MUST HAVE GROWTH FACTOR OF 4 */
   /* BECAUSE OF INTWK USAGE */
   
   if (!initialized) {
      /* VERTEX STORAGE ALLOCATION */
      maxvst =  MAX((int) (inmesh.maxvst),10);
      maxsbel = MAX(inmesh.maxsbel,10);
      allocate(maxvst);
      bdryalloc(maxsbel);
      nsbd = inmesh.nsbd;
      nvbd = inmesh.nvbd;
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
      xpt = 0.5*(vrtx[v0][0] +vrtx[v1][0]);
      ypt = 0.5*(vrtx[v0][1] +vrtx[v1][1]);
               
      /* INSERT POINT */
      vrtx[count][0] = xpt;
      vrtx[count++][1] = ypt;
   }
   
   /*	INSERT INTERIOR POINTS */
   for(i=nvrtx;i<count;++i) { 
      qtree.addpt(nvrtx);
      qtree.nearpt(nvrtx,vnear);
      tind = findtri(vrtx[nvrtx][0],vrtx[nvrtx][1],vnear);
      err = insert(tind,nvrtx,0.0);
      ++nvrtx;
   }

   
   /* INSERT BOUNDARY POINTS */
   for(i=0;i<nsbd;++i) {
      initialsidenumber = sbdry[i].num;
      for(j=0;j<initialsidenumber;++j) {
         sind = sbdry[i].el[j];
         v0 = svrtx[sind][0];
         v1 = svrtx[sind][1];
         
         /* MIDPOINT */
         xpt = 0.5*(vrtx[v0][0] +vrtx[v1][0]);
         ypt = 0.5*(vrtx[v0][1] +vrtx[v1][1]);
         
         bid = sbdry[(-stri[sind][1]>>16) -1].type;
         if (bid&CURV_MASK) mvpttobdry(bid,xpt,ypt);
                  
         /* INSERT POINT */
         vrtx[nvrtx][0] = xpt;
         vrtx[nvrtx][1] = ypt;
         
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
      bdrysidereorder(i);
      
   cnt_nbor();
   
   return;
}
   
   

   
   
         

   