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

/* USES VINFO/SINFO TO STORE SIDE LOOKUP */
/* USES NNBOR TO STORE LIST OF VERTICES */
/* USES TINFO TO STORE NEW SIDE LIST */ 
/* THIS OBVIOUSLY NEEDS TO BE CLEANED UP A BIT */
/* BUT IT DOESN'T INTERFERE WITH ANYTHING */

void mesh::triangulate(int **sidelst, int *nsdloop, int nloop, int cknbor = 1) {
   int i,j,nv;
   int lp, sd;
   int sind,dir,sindnxt,dirnxt;
   int sind1,sindprev;
   int nsidebefore,newside;
   int ngood,goodvrt[MAXGOOD];
   int minv,maxv;
   int v1,v2,v[4];
   
/*	CREATE VERTEX LIST */
/*	STORE IN NNBOR SINCE THIS IS OBVIOUSLY UNUSED RIGHT NOW */
   nv = 0;
   for(i=0;i<nloop;++i) {
      for(j=0;j<nsdloop[i];++j) {
         sind = abs(sidelst[i][j]) -1;
         dir = (1 -SIGN(sidelst[i][j]))/2;
         stri[sind][dir] = -1;
         nnbor[nv++] = svrtx[sind][dir];
         assert(nv < maxvst-1);
      }
   }

/*	SETUP SIDE POINTER INFO */
   for(i=0;i<nv;++i)
      vinfo[nnbor[i]] = -1;
      
   for(i=0;i<nloop;++i) {
      for(j=0;j<nsdloop[i];++j) {
         sind = abs(sidelst[i][j]) -1;
         v1 = svrtx[sind][0];
         v2 = svrtx[sind][1];
         if (v2 > v1) {
            minv = v1;
            maxv = v2;
         }
         else {
            minv = v2;
            maxv = v1;
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
   }
   
/*	PUT PERIODICITY ENDSIDE ON SIDELST */
   for(lp=0;lp<nloop;++lp) 
      sidelst[lp][nsdloop[lp]] = sidelst[lp][0];
         
   do {
      newside = 0;
      
      for(lp=0;lp<nloop;++lp) {
         sind = abs(sidelst[lp][nsdloop[lp]-1]) -1;
         dir =  (1 -SIGN(sidelst[lp][nsdloop[lp]-1]))/2;
         v[1] = svrtx[sind][dir];
         for(sd=0;sd<nsdloop[lp];++sd) {
            if (cknbor) {
               v[0] = v[1];
               sind = abs(sidelst[lp][sd+1]) -1;
               dir = (1 -SIGN(sidelst[lp][sd+1]))/2;
               v[3] = svrtx[sind][1 -dir];
            }
            sind = abs(sidelst[lp][sd]) -1;
            dir =  (1 -SIGN(sidelst[lp][sd]))/2;
            v[1] = svrtx[sind][dir];
            v[2] = svrtx[sind][1 -dir];

            if (stri[sind][dir] > -1) { // SIDE HAS ALREADY BEEN MATCHED 
               sind = sindnxt;
               dir = dirnxt;
               continue; 
            }
            findpt(nnbor,nv,v,cknbor,goodvrt,ngood);
            
            nsidebefore = nside;
            addtri(v[1],v[2],goodvrt[0],sind,dir);
 /*			ADD ANY DEGENERATE TRIANGLES */           
            for(i=1;i<ngood;++i)
               addtri(goodvrt[i-1],v[2],goodvrt[i],-1,-1);

/*				STORE NEWSIDES IN TINFO SINCE THIS IS ALSO UNUSED RIGHT NOW */ 
            assert(newside +nside - nsidebefore < maxvst);
            for(i=nsidebefore;i<nside;++i)
               tinfo[newside++] = i;
         }
      }
     
      cknbor = 0;
      nloop = 1;
      nsdloop[0] = 0;
      for(i=0;i<newside;++i) {
         if (stri[tinfo[i]][1] < 0) {
            sidelst[0][nsdloop[0]++] = -(tinfo[i] + 1);
            assert(nsdloop[0] < maxsbel*nsbd -1);
         }
      }
         
      nv = 0;
      for(i=0;i<nsdloop[0];++i) {
         sind = -sidelst[0][i] -1;
         nnbor[nv++] = svrtx[sind][1];
      }

   } while(newside > 0);
   
   return;
}
            
            


void mesh::findpt(int *nnbor,int nv,int *v,int chkadj,int good[], int &ngood) {
   int i,j,k,ncnvx,cnvx[2],vtry,itemp;
   FLT dx1,dy1,dx2,dy2,dx2a,dy2a;
   FLT xmid,ymid,area,alpha,beta,xcen,ycen;
   FLT hmin,height;
   FLT temp,ang[MAXGOOD],ds1,ds2;
   
   dx2 = vrtx[v[2]][0] -vrtx[v[1]][0];
   dy2 = vrtx[v[2]][1] -vrtx[v[1]][1];
   xmid = 0.5*(vrtx[v[2]][0] +vrtx[v[1]][0]);
   ymid = 0.5*(vrtx[v[2]][1] +vrtx[v[1]][1]);
   
   hmin = 1.0e99;

   ncnvx = 0;
   if (chkadj) {
/*		CHECK WHETHER ADJACENT SIDES ARE CONCAVE OR CONVEX */
      dx1 = vrtx[v[1]][0] -vrtx[v[0]][0];
      dy1 = vrtx[v[1]][1] -vrtx[v[0]][1]; 
      area        = dx1*dy2 -dy1*dx2;
      if (area > 0.0)
         cnvx[ncnvx++] = 0;

      dx1 = vrtx[v[1]][0] -vrtx[v[3]][0];
      dy1 = vrtx[v[1]][1] -vrtx[v[3]][1]; 
      area        = dx1*dy2 -dy1*dx2;
      if (area > 0.0)
         cnvx[ncnvx++] = 2;
   }
   
/*	FIND NODES WHICH MAKE POSITIVE TRIANGLE WITH SIDE */
   for(i=0;i<nv;++i) {
      vtry = nnbor[i];
      dx1 = vrtx[v[1]][0] -vrtx[vtry][0];
      dy1 = vrtx[v[1]][1] -vrtx[vtry][1]; 
      area        = dx1*dy2 -dy1*dx2;
      if (area < FLT_EPSILON) continue;
      
/*		CIRCUMCENTER IS AT INTERSECTION OF NORMAL TO SIDES THROUGH MIDPOINT */
      area = 1./area;
      alpha       = dx2*xmid +dy2*ymid;
      beta        = .5*(dx1*(vrtx[vtry][0] +vrtx[v[1]][0]) +dy1*(vrtx[vtry][1] +vrtx[v[1]][1]));
      xcen = area*(beta*dy2 -alpha*dy1);
      ycen = area*(alpha*dx1 -beta*dx2);

/*		FIND TRIANGLE FOR WHICH THE HEIGHT OF THE CIRCUMCENTER */
/*		ABOVE THE EDGE MID-POINT IS MINIMIZED */
      height = dx2*(ycen -ymid) -dy2*(xcen -xmid);

      if (height > hmin +10.*FLT_EPSILON) continue;
      
/*		CHECK FOR INTERSECTION WITH CONVEX SIDES */
      for(k=0;k<ncnvx;++k) {
         j = cnvx[k];
         dx1 = vrtx[v[j]][0]-vrtx[vtry][0];
         dy1 = vrtx[v[j]][1]-vrtx[vtry][1];
         dx2a = vrtx[v[j+1]][0]-vrtx[v[j]][0];
         dy2a = vrtx[v[j+1]][1]-vrtx[v[j]][1];
   	
         area = dx1*dy2a -dy1*dx2a;
         if (area < 0.0) goto NEXT;
      }

/*		CHECK IF DEGENERATE */
      if (height > hmin-10.*FLT_EPSILON) {
         good[ngood++] = vtry;
         assert(ngood < MAXGOOD);
         continue;
      }
      
      ngood = 0;
      good[ngood++] = vtry;
      hmin = height;
NEXT: continue;
   }
   
   if (ngood > 1) {
/*		ORDER COCIRCULAR POINTS */   
/*		CALCULATE SIDE ANGLE */
      ds2 = 1./sqrt(dx2*dx2 +dy2*dy2);
      for(i=0;i<ngood;++i) {
         vtry = good[i];
         dx1 = vrtx[v[1]][0] -vrtx[vtry][0];
         dy1 = vrtx[v[1]][1] -vrtx[vtry][1];
         ds1 = 1./sqrt(dx1*dx1 +dy1*dy1);
         ang[i] = -(dx2*dx1  +dy2*dy1)*ds2*ds1;
      }
   
/*		ORDER POINTS BY ANGLE */      
      for(i=0;i<ngood-1;++i) {
         for(j=i+1;j<ngood;++j) {

/*				TO ELIMINATE POSSIBILITY OF REPEATED VERTICES IN NNBOR */
            if (good[i] == good[j]) {
               good[j] = good[ngood-1];
               --ngood;
            }

/*				ORDER BY ANGLE */
            if(ang[i] > ang[j]) {
               temp = ang[i];
               ang[i] = ang[j];
               ang[j] = temp;
               itemp = good[i];
               good[i] = good[j];
               good[j] = itemp;
            }
         }
      }
   }
   
   return;
}

void mesh::addtri(int v0,int v1, int v2, int sind, int dir) {
   int i,j,k,end,sind1,tind;
   int minv,maxv,order,sindprev,temp;
      
/*	ADD NEW TRIANGLE */
   tvrtx[ntri][0] = v0;
   tvrtx[ntri][1] = v1;
   tvrtx[ntri][2] = v2;
               
   vtri[v0] = ntri;
   vtri[v1] = ntri;
   vtri[v2] = ntri;
   
   end = 3;
   if (sind > -1) {
/*		SIDE 2 INFO IS KNOWN ALREADY */
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

/*	LOOP THROUGH SIDES */
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
/*				SIDE IN SAME DIRECTION */
            if (stri[sind1][0] >= 0) {
               printf("1:side already matched?%d %d %d\n",sind1,v1,v2);
               out_mesh("error");
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
/*							SIDE IN OPPOSITE DIRECTION */
            if (stri[sind1][1] >= 0) {
               printf("2:side already matched? %d %d %d\n",sind1,v1,v2);
               out_mesh("error");
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
/*		NEW SIDE */
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
