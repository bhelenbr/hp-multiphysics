/*
 *  gtol.cpp
 *  planar++
 *
 *  Created by helenbrk on Sun Oct 14 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "spectral_hp.h"
#include<assert.h>

void spectral_hp::ugtouht(int tind) {
    int i,k,m,n,indx,cnt;
    int sign, msgn;

   /* THIS IS FOR FLOW VARIABLES ON ANY MESH */
   /* VERTICES */   
   for (i=0; i<3; ++i) {
      indx = tvrtx[tind][i];
      for(n=0; n<NV; ++n)
         uht[n][i] = ug.v[indx][n];
   }
   

   /* SIDES */
   cnt = 3;
   for(i=0;i<3;++i) {
      indx = tside[tind].side[i]*sm0;
      sign = tside[tind].sign[i];
      msgn = 1;
      for (m = 0; m < b->sm; ++m) {
         for(n=0; n<NV; ++n)
            uht[n][cnt] = msgn*ug.s[indx +m][n];
         msgn *= sign;
         ++cnt;
      }
   }
 
   /* INTERIORS */   
   if (b->im > 0) {   
      indx = tind*im0;
      for(m = 1; m < b->sm; ++m) {
         for(k = 0; k < b->sm-m; ++k) {
            for(n=0; n<NV; ++n)
               uht[n][cnt] = ug.i[indx][n];
            ++cnt; ++indx;
         }
         indx += sm0 -b->sm;
      }
   }
   
   return;
}

void spectral_hp::ugtouht(int tind, struct vsi ug) {
    int i,k,m,n,indx,cnt;
    int sign, msgn;

   /* THIS IS FOR FLOW VARIABLES ON ANY MESH */
   /* VERTICES */   
   for (i=0; i<3; ++i) {
      indx = tvrtx[tind][i];
      for(n=0; n<NV; ++n)
         uht[n][i] = ug.v[indx][n];
   }
   

   /* SIDES */
   /* SIDES */
   cnt = 3;
   for(i=0;i<3;++i) {
      indx = tside[tind].side[i]*sm0;
      sign = tside[tind].sign[i];
      msgn = 1;
      for (m = 0; m < b->sm; ++m) {
         for(n=0; n<NV; ++n)
            uht[n][cnt] = msgn*ug.s[indx +m][n];
         msgn *= sign;
         ++cnt;
      }
   }

   /* INTERIORS */   
   if (b->im > 0) {   
      indx = tind*im0;
      for(m = 1; m < b->sm; ++m) {
         for(k = 0; k < b->sm-m; ++k) {
            for(n=0; n<NV; ++n)
               uht[n][cnt] = ug.i[indx][n];
            ++cnt; ++indx;
         }
         indx += sm0 -b->sm;
      }
   }
   
   return;
}

void spectral_hp::ugtouht_bdry(int tind) {
    int i,m,n,indx,cnt;
    int sign, msgn;
   
   for (i=0; i<3; ++i) {
      indx = tvrtx[tind][i];
      for(n=0; n<NV; ++n)
         uht[n][i] = ug.v[indx][n];
   }

   /* SIDES */
   cnt = 3;
   for(i=0;i<3;++i) {
      indx = tside[tind].side[i]*sm0;
      sign = tside[tind].sign[i];
      msgn = 1;
      for (m = 0; m < b->sm; ++m) {
         for(n=0; n<NV; ++n)
            uht[n][cnt] = msgn*ug.s[indx +m][n];
         msgn *= sign;
         ++cnt;
      }
   }
   
   return;
}

void spectral_hp::ugtouht_bdry(int tind,struct vsi ug) {
    int i,m,n,indx,cnt;
    int sign, msgn;
   
   for (i=0; i<3; ++i) {
      indx = tvrtx[tind][i];
      for(n=0; n<NV; ++n)
         uht[n][i] = ug.v[indx][n];
   }

   /* SIDES */
   cnt = 3;
   for(i=0;i<3;++i) {
      indx = tside[tind].side[i]*sm0;
      sign = tside[tind].sign[i];
      msgn = 1;
      for (m = 0; m < b->sm; ++m) {
         for(n=0; n<NV; ++n)
            uht[n][cnt] = msgn*ug.s[indx +m][n];
         msgn *= sign;
         ++cnt;
      }
   }
      
   return;
}


void spectral_hp::ugtouht1d(int sind) {
   int m,n,indx,v0,v1;
   
   v0 = svrtx[sind][0];
   v1 = svrtx[sind][1];
   for(n=0;n<NV;++n) {
      uht[n][0] = ug.v[v0][n];
      uht[n][1] = ug.v[v1][n];
   }
    
   indx = sind*sm0;
   for(m=0;m<b->sm;++m)
    for(n=0;n<NV;++n) 
      uht[n][m+2] = ug.s[indx+m][n];
         
   return;
}

void spectral_hp::ugtouht1d(int sind, struct vsi ug) {
   int m,n,indx,v0,v1;
   
   v0 = svrtx[sind][0];
   v1 = svrtx[sind][1];
   for(n=0;n<NV;++n) {
      uht[n][0] = ug.v[v0][n];
      uht[n][1] = ug.v[v1][n];
   }
    
   indx = sind*sm0;
   for(m=0;m<b->sm;++m)
    for(n=0;n<NV;++n) 
      uht[n][m+2] = ug.s[indx+m][n];
         
   return;
}

void spectral_hp::crdtocht(int tind) {
   int i,m,n,cnt,bnum,sind,indx;
   
   /* VERTICES */   
   for (i=0; i<3; ++i) {
      indx = tvrtx[tind][i];
      for(n=0; n<ND; ++n)
         cht[n][i] = vrtx[indx][n];
   }

   if (b->sm == 0) return;
   
   /* SIDES */
   cnt = 3;
   for (i=0; i<3;++i) {   
      sind = tside[tind].side[i];
      if (sinfo[sind] < 0) {
         for(m=0;m<b->sm;++m) {
            for(n=0;n<ND;++n)
               cht[n][cnt] = 0.0;
            ++cnt;
         }
      }
      else {
         bnum = (-stri[sind][1]>>16) -1;
         indx = (-stri[sind][1]&0xFFFF)*sm0;
         assert(bnum > -1 && bnum < nsbd);
         assert(indx > -1 && indx < sbdry[bnum].num*sm0);
         
         for (m = 0; m < b->sm; ++m) {
            for(n=0; n<ND; ++n)
               cht[n][cnt] = binfo[bnum][indx+m].curv[n];
            ++cnt;
         }
      }
   }
   
   return;
}

void spectral_hp::crdtocht(int tind, FLT (*vrtx)[ND], struct bistruct **binfo) {
   int i,m,n,cnt,bnum,sind,indx;
   
   /* VERTICES */   
   for (i=0; i<3; ++i) {
      indx = tvrtx[tind][i];
      for(n=0; n<ND; ++n)
         cht[n][i] = vrtx[indx][n];
   }

   if (b->sm == 0) return;
   
   /* SIDES */
   cnt = 3;
   for (i=0; i<3;++i) {   
      sind = tside[tind].side[i];
      if (sinfo[sind] < 0) {
         for(m=0;m<b->sm;++m) {
            for(n=0;n<ND;++n)
               cht[n][cnt] = 0.0;
            ++cnt;
         }
      }
      else {
         bnum = (-stri[sind][1]>>16) -1;
         indx = (-stri[sind][1]&0xFFFF)*sm0;
         assert(bnum > -1 && bnum < nsbd);
         assert(indx > -1 && indx < sbdry[bnum].num*sm0);
         
         for (m = 0; m < b->sm; ++m) {
            for(n=0; n<ND; ++n)
               cht[n][cnt] = binfo[bnum][indx+m].curv[n];
            ++cnt;
         }
      }
   }
   
   return;
}

void spectral_hp::crdtocht1d(int sind) {
   int m,n,bnum,indx,v0,v1;
   
   v0 = svrtx[sind][0];
   v1 = svrtx[sind][1];
   for(n=0;n<ND;++n) {
      cht[n][0] = vrtx[v0][n];
      cht[n][1] = vrtx[v1][n];
   }
   
   bnum = (-stri[sind][1]>>16) -1;
   indx = (-stri[sind][1]&0xFFFF)*sm0;
   assert(bnum > -1 && bnum < nsbd);
   assert(indx > -1 && indx <= sbdry[bnum].num*sm0);
   
   for(m=0;m<b->sm;++m)
    for(n=0;n<ND;++n) 
         cht[n][m+2] = binfo[bnum][indx+m].curv[n];
         
   return;
}

void spectral_hp::crdtocht1d(int sind,FLT (*vrtx)[ND], struct bistruct **binfo) {
   int m,n,bnum,indx,v0,v1;
   
   v0 = svrtx[sind][0];
   v1 = svrtx[sind][1];
   for(n=0;n<ND;++n) {
      cht[n][0] = vrtx[v0][n];
      cht[n][1] = vrtx[v1][n];
   }
   
   bnum = (-stri[sind][1]>>16) -1;
   indx = (-stri[sind][1]&0xFFFF)*sm0;
   assert(bnum > -1 && bnum < nsbd);
   assert(indx > -1 && indx <= sbdry[bnum].num*sm0);
   
   for(m=0;m<b->sm;++m)
    for(n=0;n<ND;++n) 
         cht[n][m+2] = binfo[bnum][indx+m].curv[n];
         
   return;
}

void spectral_hp::lftog(int tind, struct vsi g) {
   int i,m,n,indx,gindx,sgn,msgn;
   
   /* VERTEX MODES */
   for (m=0; m<3; ++m) {
      gindx = tvrtx[tind][m];
      for(n=0;n<NV;++n)
         g.v[gindx][n] += lf[n][m];
   }


   if (b->p > 1) {
      /* SIDE MODES */
      indx = 3;
      for(i=0;i<3;++i) {
         gindx = tside[tind].side[i]*b->sm;
         sgn = tside[tind].sign[i];
         msgn = 1;
         for (m = 0; m < b->sm; ++m) {
            for(n=0;n<NV;++n)
               g.s[gindx][n] += msgn*lf[n][indx];
            msgn *= sgn;
            ++gindx;
            ++indx;
         }
      }

      gindx = tind*b->im;
      indx = b->bm;
      for(m=0;m<b->im;++m) {
         for(n=0;n<NV;++n)
            g.i[gindx][n] += lf[n][indx];
         ++gindx;
         ++indx;
      }
   }
   
   return;
}
