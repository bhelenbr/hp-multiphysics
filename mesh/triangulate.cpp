/*
 *  triangulate.cpp
 *  mblock
 *
 *  Created by helenbrk on Mon Aug 13 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#include "mesh.h"
#include "utilities.h"
#include<float.h>
#include<assert.h>

#define MAXGOOD 100

template void mesh<2>::triangulate(int nsd);
template void mesh<3>::triangulate(int nsd);

template<int ND> void mesh<ND>::triangulate(int nsd) {
   int i,j,n,vcnt,sck,dirck,stest,nv,vtry;
   int sd;
   int bgn,end;
   int sind,dir;
   int sind1,sindprev;
   int nsidebefore,ntest;
   int ngood,good[MAXGOOD];
   int minv,maxv,itemp;
   int v[2],v2,v3;
   FLT xmid[2],xcen[2];
   FLT xm1dx1,xm2dx2;
   FLT hmin,height;
   FLT temp,ang[MAXGOOD],ds1,ds2;
   FLT xmax[ND],xmin[ND],xmax1[ND],xmin1[ND];
   FLT dx1[ND],dx2[ND],dx3[ND],dx4[ND];
   FLT det,s,t;

   /* CREATE VERTEX LIST */
   /* STORE IN NNBOR SINCE THIS IS UNUSED RIGHT NOW */
   nv = 0;
   for(i=0;i<nsd;++i) {
      sind = abs(intwk1[i]) -1;
      dir = (1 -SIGN(intwk1[i]))/2;
      stri[sind][dir] = -1;
      intwk2[nv++] = svrtx[sind][dir];
      assert(nv < maxvst-1);
   }

   /* SETUP SIDE POINTER INFO */
   for(i=0;i<nv;++i)
      vinfo[intwk2[i]] = -1;
      
   for(i=0;i<nsd;++i) {
      sind = abs(intwk1[i]) -1;
      v[0] = svrtx[sind][0];
      v[1] = svrtx[sind][1];
      if (v[1] > v[0]) {
         minv = v[0];
         maxv = v[1];
      }
      else {
         minv = v[1];
         maxv = v[0];
      }
      sind1 = vinfo[minv];
      while (sind1 >= 0) {
         sindprev = sind1;
         sind1 = sinfo[sind1];
      }
      sinfo[sind] = -1;
      if (vinfo[minv] < 0)
         vinfo[minv] = sind;
      else 
         sinfo[sindprev] = sind;
   }
   
   bgn = 0;
   end = nsd;
   ntest = end;
   while(bgn < end) {
      nsidebefore = nside;
      for(sd=bgn;sd<end;++sd) {
         sind = abs(intwk1[sd]) -1;
         dir =  (1 -SIGN(intwk1[sd]))/2;
         
         if (stri[sind][dir] > -1) continue; // SIDE HAS ALREADY BEEN MATCHED 
         
         v[0] = svrtx[sind][dir];
         v[1] = svrtx[sind][1 -dir];

         /* SEARCH FOR GOOD POINTS */
         for(n=0;n<ND;++n) {
            dx2[n] = vrtx[v[1]][n] -vrtx[v[0]][n];
            xmid[n] = 0.5*(vrtx[v[1]][n] -vrtx[v[0]][n]);
         }
         hmin = 1.0e99;
         
         /* FIND NODES WHICH MAKE POSITIVE TRIANGLE WITH SIDE */
         for(i=0;i<nv;++i) {
            vtry = intwk2[i];
            if (vtry == v[0] || vtry == v[1]) continue;
      
            for(n=0;n<ND;++n)
               dx1[n] = vrtx[v[0]][n] -vrtx[vtry][n];
            det        = dx1[0]*dx2[1] -dx1[1]*dx2[0];
            if (det <= 0.0) continue;
            
            /* CIRCUMCENTER IS AT INTERSECTION OF NORMAL TO SIDES THROUGH MIDPOINT */
            det = 1./det;
            xm1dx1       = -0.5*(dx1[0]*dx1[0] +dx1[1]*dx1[1]);
            xm2dx2       =  0.5*(dx2[0]*dx2[0] +dx2[1]*dx2[1]);
            xcen[0] = det*(xm1dx1*dx2[1] -xm2dx2*dx1[1]);
            xcen[1] = det*(xm2dx2*dx1[0] -xm1dx1*dx2[0]);
                  
            /* FIND TRIANGLE FOR WHICH THE HEIGHT OF THE CIRCUMCENTER */
            /* ABOVE THE EDGE MID-POINT IS MINIMIZED (MINIMIZES RADIUS) */
            height = dx2[0]*(xcen[1] -xmid[1]) -dx2[1]*(xcen[0] -xmid[0]);

            if (height > hmin) continue;
            
            /* CHECK FOR INTERSECTION OF TWO CREATED SIDES */
            /* WITH ALL OTHER BOUNDARY SIDES */
            for(vcnt=0;vcnt<2;++vcnt) {
               minv = MIN(vtry,v[vcnt]);
               maxv = MAX(vtry,v[vcnt]);
               /* LOOK THROUGH ALL SIDES CONNECTED TO MINV FOR DUPLICATE */
               /* IF DUPLICATE THEN SIDE IS OK - NO NEED TO CHECK */
               sind1 = vinfo[minv];
               while (sind1 >= 0) {
                  if (maxv == svrtx[sind1][0] || maxv == svrtx[sind1][1]) goto next_vrt;
                  sind1 = sinfo[sind1];
               }
               
               for(n=0;n<ND;++n) {
                  xmin[n]      = MIN(vrtx[v[0]][n],vrtx[v[1]][n]);
                  xmax[n]      = MAX(vrtx[v[0]][n],vrtx[v[1]][n]);
               }
               
               for(sck=0;sck<ntest;++sck) {
                  stest = abs(intwk1[sck])-1;
                  if (stest == sind) continue;
                  dirck =  (1 -SIGN(intwk1[sck]))/2;
                  v2 = svrtx[stest][dirck];
                  v3 = svrtx[stest][1-dirck];
                  if (v2 == minv || v3 == minv) {
                     /* SPECIAL TEST FOR CONNECTED SIDES */
                     /* CAN ONLY FAIL IF CONVEX BOUNDARY (LOOKING FROM OUTSIDE) */
                     if (area(v2,v3,v[vcnt]) > 0.0 && area(v2,v3,maxv) <= 0.0) goto vtry_failed;
                     continue;
                  }
                  if (v2 == maxv || v3 == maxv) {
                     /* SPECIAL TEST FOR CONNECTED SIDES */
                     /* CAN ONLY FAIL IF CONVEX BOUNDARY (LOOKING FROM OUTSIDE) */
                     if (area(v2,v3,v[vcnt]) > 0.0 && area(v2,v3,minv) <= 0.0) goto vtry_failed;
                     continue;
                  }
                  for(n=0;n<ND;++n) {
                     xmin1[n]      = MIN(vrtx[v2][n],vrtx[v3][n]);
                     xmax1[n]      = MAX(vrtx[v2][n],vrtx[v3][n]);
                  }
                  for(n=0;n<ND;++n)
                     if (xmax[n] < xmin1[n] || xmin[n] > xmax1[n]) goto next_bdry_side;
      
                  for(n=0;n<ND;++n) {
                     dx1[n] = vrtx[maxv][n] -vrtx[minv][n];
                     dx3[n] = vrtx[v3][n]-vrtx[v2][n];
                     dx4[n] = vrtx[minv][n]-vrtx[v2][n];
                  }
                  
                  det = -dx1[0]*dx3[1] +dx1[1]*dx3[0];
                  if (det < EPSILON*100.0*(fabs(xmax[0])+fabs(xmax[1]))) continue;
                  
                  det = 1./det;
                  s = det*(dx4[0]*dx3[1] -dx3[0]*dx4[1]);
                  t = det*(-dx1[0]*dx4[1] +dx4[0]*dx1[1]);
                  
                  if (s < 0.0 || s > 1.0) continue;
                  if (t < 0.0 || t > 1.0) continue;
                  
                  goto vtry_failed;
                  
next_bdry_side:   continue;
               }
next_vrt:      continue;
            }
      
            /* CHECK IF DEGENERATE */
            if (height > hmin-200.0*EPSILON*sqrt(xcen[0]*xcen[0]+xcen[1]*xcen[1])) {
               good[ngood++] = vtry;
               assert(ngood < MAXGOOD);
               continue;
            }
            
            ngood = 0;
            good[ngood++] = vtry;
            hmin = height+100.0*EPSILON*sqrt(xcen[0]*xcen[0]+xcen[1]*xcen[1]);
vtry_failed:continue;
         }
         
         if (ngood > 1) {
            /* ORDER COCIRCULAR POINTS */   
            /* CALCULATE SIDE ANGLE */
            ds2 = 1./sqrt(dx2[0]*dx2[0] +dx2[1]*dx2[1]);
            for(i=0;i<ngood;++i) {
               vtry = good[i];
               dx1[0] = vrtx[v[0]][0] -vrtx[vtry][0];
               dx1[1] = vrtx[v[0]][1] -vrtx[vtry][1];
               ds1 = 1./sqrt(dx1[0]*dx1[0] +dx1[1]*dx1[1]);
               ang[i] = -(dx2[0]*dx1[0]  +dx2[1]*dx1[1])*ds2*ds1;
            }
         
            /* ORDER POINTS BY ANGLE */      
            for(i=0;i<ngood-1;++i) {
               for(j=i+1;j<ngood;++j) {
      
                  /* TO ELIMINATE POSSIBILITY OF REPEATED VERTICES IN intwk2 */
                  if (good[i] == good[j]) {
                     good[j] = good[ngood-1];
                     --ngood;
                  }
      
                  /* ORDER BY ANGLE */
                  if(ang[i] > ang[j]) {
                     temp = ang[i];
                     ang[i] = ang[j];
                     ang[j] = temp;
                     itemp = good[i];
                     good[i] = good[j];
                     good[j] = itemp;
                  }
               }
               // *log << "degenerate case" << v[0] << ' ' << v[1] << std::endl;
            }
         }
         addtri(v[0],v[1],good[0],sind,dir);
         /* ADD ANY DEGENERATE TRIANGLES */           
         for(i=1;i<ngood;++i)
            addtri(good[i-1],v[1],good[i],-1,-1);
      }
     
      bgn = end;
      end += nside -nsidebefore;
      for(i=nsidebefore;i<nside;++i)
         intwk1[bgn+i-nsidebefore] = -(i + 1);
   }
   
   for(i=0;i<maxvst;++i) {
      intwk1[i] = -1;
      intwk2[i] = -1;
   }
      
   return;

}


template void mesh<2>::addtri(int v0,int v1, int v2, int sind, int dir);
template void mesh<3>::addtri(int v0,int v1, int v2, int sind, int dir);

template<int ND> void mesh<ND>::addtri(int v0,int v1, int v2, int sind, int dir) {
   int i,j,k,end,sind1,tind;
   int minv,maxv,order,sindprev,temp;
      
   /* ADD NEW TRIANGLE */
   tvrtx[ntri][0] = v0;
   tvrtx[ntri][1] = v1;
   tvrtx[ntri][2] = v2;
               
   vtri[v0] = ntri;
   vtri[v1] = ntri;
   vtri[v2] = ntri;
   
   end = 3;
   if (sind > -1) {
      /* SIDE 2 INFO IS KNOWN ALREADY */
      tside[ntri].side[2] = sind;
      tside[ntri].sign[2] = 1 -2*dir;
      stri[sind][dir] = ntri;
      tind = stri[sind][1-dir];
      ttri[ntri][2] = tind;
      if (tind > -1) {
         for(i=0;i<3;++i) {
            if (tside[tind].side[i] == sind) {
               ttri[tind][i] = ntri;
               break;
            }
         }
      }
      end = 2;
   }

   /* LOOP THROUGH SIDES */
   for(k=0;k<end;++k) {
      if (v2 > v1) {
         minv = v1;
         maxv = v2;
         order = 1;
      }
      else {
         minv = v2;
         maxv = v1;
         order = 0;
      }
      sind1 = vinfo[minv];
      while (sind1 >= 0) {
         if (maxv == svrtx[sind1][order]) {
            /* SIDE IN SAME DIRECTION */
            if (stri[sind1][0] >= 0) {
               *log << "1:side already matched?" << sind1 << ' ' << v1 << ' ' << v2 << std::endl;
               out_mesh("error",ftype::tecplot);
               out_mesh("error",ftype::grid);
               exit(1);
            }
            stri[sind1][0] = ntri;
            tside[ntri].side[k] = sind1;
            tside[ntri].sign[k] = 1;
            tind = stri[sind1][1];
            ttri[ntri][k] = tind;
            if (tind > -1) {
               for(j=0;j<3;++j) {
                  if (tside[tind].side[j] == sind1) {
                     ttri[tind][j] = ntri;
                     break;
                  }
               }
            }
            goto NEXTTRISIDE;
         }
         else if(maxv == svrtx[sind1][1-order]) {
            /* SIDE IN OPPOSITE DIRECTION */
            if (stri[sind1][1] >= 0) {
               *log << "2:side already matched?" << sind1 << ' ' << v1 << ' ' << v2 << std::endl;
               out_mesh("error",ftype::tecplot);
               out_mesh("error",ftype::grid);
               exit(1);
            }
            stri[sind1][1] = ntri;
            tside[ntri].side[k] = sind1;
            tside[ntri].sign[k] = -1;
            tind = stri[sind1][0];
            ttri[ntri][k] = tind;
            if (tind > -1) {
               for(j=0;j<3;++j) {
                  if (tside[tind].side[j] == sind1) {
                     ttri[tind][j] = ntri;
                     break;
                  }
               }
            }
            goto NEXTTRISIDE;     
         }
         sindprev = sind1;
         sind1 = sinfo[sind1];
      }
      /* NEW SIDE */
      svrtx[nside][0] = v1;
      svrtx[nside][1] = v2;
      stri[nside][0] = ntri;
      stri[nside][1] = -1;
      tside[ntri].side[k] = nside;
      tside[ntri].sign[k] = 1;
      sinfo[nside] = -1;
      if (vinfo[minv] < 0)
         vinfo[minv] = nside;
      else 
         sinfo[sindprev] = nside;
      ++nside;
      assert(nside < maxvst -1);
         
NEXTTRISIDE:
      temp = v0;
      v0 = v1;
      v1 = v2;
      v2 = temp;
   }
   ++ntri;
   
   assert(ntri < maxvst -1);
      
   return;
}

