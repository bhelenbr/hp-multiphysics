/*
 *  curvinit.cpp
 *  planar++
 *
 *  Created by helenbrk on Tue Oct 16 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include"spectral_hp.h"
#include<myblas.h>

void spectral_hp::curvinit(int MASK) {
   int i,j,m,n,bind,indx,sind,v0,v1,info,typ;
   FLT x,y;
   char uplo[] = "U";
   
   for(bind=0;bind<nvbd;++bind) {
      if (!(vbdry[bind].type&MVPT_MASK)) continue;
      
      v0 = vbdry[bind].el[0];
      x = vrtx[v0][0];
      y = vrtx[v0][1];
      mvpttobdry(vbdry[bind].type,x,y);
      vrtx[v0][0] = x;
      vrtx[v0][1] = y;
   }

   /* MAKE SURE VERTICES ARE ON SURFACE */
   for(bind=0;bind<nsbd;++bind) {
      typ = sbdry[bind].type;

      if (!(typ&CURV_MASK) || !(typ&MASK)) continue;

      /* MESS WITH END VERTICES? */
      for(j=1;j<sbdry[bind].num;++j) {
         sind = sbdry[bind].el[j];
         v0 = svrtx[sind][0];
         x = vrtx[v0][0];
         y = vrtx[v0][1];
         mvpttobdry(sbdry[bind].type,x,y);
         vrtx[v0][0] = x;
         vrtx[v0][1] = y;
      }
//      v0 = svrtx[sind][1];
//      x = vrtx[v0][0];
//      y = vrtx[v0][1];
//      mvpttobdry(sbdry[bind].type,x,y);
//      vrtx[v0][0] = x;
//      vrtx[v0][1] = y;
   }

   if (b->p == 1) return;
      
   /*****************************/
   /* SET UP HIGHER ORDER MODES */
   /*****************************/
   for(bind=0;bind<nsbd;++bind) {
      typ = sbdry[bind].type;

      if (!(typ&MASK)) continue;

      if (!(typ&CURV_MASK)) {
         for(j=0;j<sbdry[bind].num*sm0;++j)
            for(n=0;n<ND;++n)
               binfo[bind][j].curv[n] = 0.0;
         continue;
      }
      

      for(j=0;j<sbdry[bind].num;++j) {
         sind = sbdry[bind].el[j];

         v0 = svrtx[sind][0];
         v1 = svrtx[sind][1];
         
         for(n=0;n<ND;++n) 
            b->proj1d(vrtx[v0][n],vrtx[v1][n],crd[n][0]);
      
         for(i=0;i<b->gpx;++i) {
            x = crd[0][0][i];
            y = crd[1][0][i];
            mvpttobdry(typ,x,y);
            
            crd[0][0][i] -= x;
            crd[1][0][i] -= y;
         }
         
         indx = j*sm0;
         for(n=0;n<ND;++n) {
            b->intgrt1d(cf[n],crd[n][0]);
            PBTRS(uplo,b->sm,b->sbwth,1,&b->sdiag1d(0,0),b->sbwth+1,&cf[n][2],b->sm,info);
         
            for(m=0;m<b->sm;++m)
               binfo[bind][indx+m].curv[n] = -cf[n][m+2];
         }
      }
   }

   return;
}

