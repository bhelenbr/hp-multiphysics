/*
 *  r_length.cpp
 *  mesh
 *
 *  Created by Brian Helenbrook on Wed Jan 09 2002.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include"r_mesh.h"

void r_mesh::perturb() {
   int i,j,sind,v0;
   
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type  == EULR_MASK) {
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            vrtx[v0][0] -= 1.9;
         }
      }
       if (sbdry[i].type  == FSRF_MASK) {
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            vrtx[v0][0] -= 0.0;
            vrtx[v0][1] -= 0.5;
         }
      }     
      
   }
   
   return;
}


void r_mesh::length1() {
   int i,j,v0,sind,bnum,count;
   class mesh *tgt;
   
   /* SET VLNGTH HERE */
//      vlngth[i] = 0.125 +0.0001*(vrtx[i][0] +vrtx[i][1]);

   /* SEND COMMUNICATIONS TO ADJACENT MESHES */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & (COMY_MASK +IFCE_MASK)) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         /* SEND VERTEX INFO */
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            tgt->sbuff[bnum][count++] = vlngth[v0];
         }
         v0 = svrtx[sind][1];
         tgt->sbuff[bnum][count++] = vlngth[v0];
      }
   }

   return;
   
}

void r_mesh::length_mp() {
   int i,j,v0,sind,bnum,count;
   class mesh *tgt;
   
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & (COMY_MASK +IFCE_MASK)) {
         count = 0;
         /* RECV VERTEX INFO */
         for(j=sbdry[i].num-1;j>=0;--j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][1];
            vlngth[v0] = 0.5*(vlngth[v0] +sbuff[i][count++]);
         }
         v0 = svrtx[sind][0];
         vlngth[v0] = 0.5*(vlngth[v0] +sbuff[i][count++]);
      }
   }
   
   /* SEND COMMUNICATIONS TO ADJACENT MESHES */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & COMX_MASK) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         /* SEND VERTEX INFO */
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            tgt->sbuff[bnum][count++] = vlngth[v0];
         }
         v0 = svrtx[sind][1];
         tgt->sbuff[bnum][count++] = vlngth[v0];
      }
   }
   return;
}

void r_mesh::length2() {
   int i,j,v0,sind,count;

   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & COMX_MASK) {
         count = 0;
         /* RECV VERTEX INFO */
         for(j=sbdry[i].num-1;j>=0;--j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][1];
            vlngth[v0] = 0.5*(vlngth[v0] +sbuff[i][count++]);
         }
         v0 = svrtx[sind][0];
         vlngth[v0] = 0.5*(vlngth[v0] +sbuff[i][count++]);
      }
   }
   
   return;
}
