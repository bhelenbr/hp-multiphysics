/*
 *  srebay.cpp
 *  mblock
 *
 *  Created by helenbrk on Wed Sep 12 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */
#include"mesh.h"
#include"utilities.h"
#include<assert.h>
#include<math.h>

#define REBAY

/* MEMORY USAGE */
/* INTWK2 - LIST OF SIDES TO BE REFINED */
/* INTWK3 - BACK REFERENCE FROM SIDE INTO LIST */
/* vinfo - TRIANGLE ACCEPTANCE INDICATOR >= 0 */

/* THIS KEEPS TRACK OF HOW MANY SIDES NEED REFINEMENT */
/* SIDES ARE STORED IN INTWK2 WITH BACK REFERENCE IN INTWK3 */
int nslst;

/* THESE ARE FROM INSERT & TELL WHICH SIDES/TRIS WERE DELETED */
extern int ntdel, tdel[MAXLST+1];
extern int nsdel, sdel[MAXLST+1];


void mesh::rebay(FLT tolsize) {
   int i,j,tind,sind,v0,v1,v2,vnear,nsnew,ntnew,snum,bid,bdrycnt,intrcnt,err;
   FLT xpt,ypt,wt[3];
   FLT dx,dy,p,q,s1sq,s2sq,rad1,rad2,rs;
   FLT xmid,ymid,densty,cirrad,arg,dist,rn1,rn2,xdif,ydif,rsign;
   
   /* SET UP FLTWK */
   fltwkreb();
         
   bdrycnt = 0;
   intrcnt = 0;
   
   /* USE VINFO[TIND] = 3,2,1,0 TO FIND BOUNDARY REFINE TRI'S */
   for(i=0;i<ntri;++i)
      vinfo[i] = 0;
   
   /* CLASSIFY SIDES AS ACCEPTED OR UNACCEPTED */
   nslst = 0;
   for(i=0;i<nside;++i) {
      if (fltwk[i] < tolsize) {
         vinfo[stri[i][0]] += 1;
         vinfo[MAX(stri[i][1],-1)] += 1;
         putinlst(i);
      }
   }
   
   vinfo[-1] = 0;

   /* BEGIN REFINEMENT ALGORITHM */
   while (nslst > 0) {
      sind = -1;
      for(i=0;i<nslst;++i) {
         if (vinfo[stri[intwk2[i]][0]] +vinfo[MAX(stri[intwk2[i]][1],-1)] == 6) continue;
         sind = intwk2[i];
         break;
      }
      assert(sind != -1);
            

      tind = stri[sind][0];
      snum = -2;
      for(j=0;j<3;++j) {
         if (tside[tind].side[j] == sind) {
            snum = j;
            break;
         }
      }
      assert(snum > -2);
         
      
      v0 = tvrtx[tind][(snum+1)%3];
      v1 = tvrtx[tind][(snum+2)%3];
      v2 = tvrtx[tind][snum];
      
      if (nvrtx > maxvst -2) {
         printf("too many vertices\n");
         exit(1);
      }
      
      
      /* CHECK IF IT IS A BOUNDARY EDGE */
      if (stri[sind][1] < 0) {

         /* MIDPOINT */
         xpt = 0.5*(vrtx[v0][0] +vrtx[v1][0]);
         ypt = 0.5*(vrtx[v0][1] +vrtx[v1][1]);
         
         bid = sbdry[(-stri[sind][1]>>16) -1].type;
         if (bid&CURV_MASK) mvpttobdry(bid,xpt,ypt);
                  
         /* INSERT POINT */
         vrtx[nvrtx][0] = xpt;
         vrtx[nvrtx][1] = ypt;
         qtree.addpt(nvrtx);
         vlngth[nvrtx] = 0.5*(vlngth[v0] +vlngth[v1]);
         
         bdry_insert(tind,snum,nvrtx);
         ++bdrycnt;
         
         nsnew = nsdel +2;
         ntnew = ntdel +1;
      }
      else {
      
#ifdef REBAY
         /* THIS IS TRIANGLE BASED REFINEMENT */
         /*	FIND ACCEPTED EDGE ON TRIANGLE & BUILD OFF OF IT */
         if (vinfo[tind] < 3) {
            for (snum=0;snum<3;++snum) 
               if (intwk3[snum] == -1) break;
         }
         else {
            tind = stri[sind][1];
            for (snum=0;snum<3;++snum) 
               if (intwk3[snum] == -1) break;
         }
         
         v0 = tvrtx[tind][(snum+1)%3];
         v1 = tvrtx[tind][(snum+2)%3];
         v2 = tvrtx[tind][snum];
         
         /* USE REBAY'S ALGORITHM FOR INSERT POINT */
         dx = vrtx[v0][0] -vrtx[v1][0];
         dy = vrtx[v0][1] -vrtx[v1][1];
         p = 0.5*sqrt(dx*dx +dy*dy);
         dx = vrtx[v2][0] -vrtx[v1][0];
         dy = vrtx[v2][1] -vrtx[v1][1];
         s1sq = 0.25*(dx*dx +dy*dy);
         dx = vrtx[v2][0] -vrtx[v0][0];
         dy = vrtx[v2][1] -vrtx[v0][1];
         s2sq = 0.25*(dx*dx +dy*dy);
         tcenter(tind,xpt,ypt);
         if (p*p > s1sq +s2sq) goto INSRT;
         
         xmid = .5*(vrtx[v0][0] +vrtx[v1][0]);
         ymid = .5*(vrtx[v0][1] +vrtx[v1][1]);
         dx = xpt -xmid;
         dy = ypt -ymid;
         q = sqrt(dx*dx +dy*dy);
         if (q < p) goto INSRT;

         densty    = (vlngth[v0]  +vlngth[v1])/sqrt(3.0);
         rad1      = MAX(densty,p);
         rad2      = .5*(p*p  +q*q)/q;
         cirrad    = MIN(rad1,rad2);
         arg       = fabs(cirrad*cirrad  -p*p);
         dist      = cirrad  +sqrt(arg);
         dx        = vrtx[v1][0] -vrtx[v0][0];
         dy        = vrtx[v1][1] -vrtx[v0][1];
         rs        = 1./sqrt(dx*dx  +dy*dy);
         rn1       = dy*rs;
         rn2       = -dx*rs;
         xdif      = vrtx[v2][0]  -xmid;
         ydif      = vrtx[v2][1]  -ymid;
         rsign     = 1.;
         if (xdif*rn1 +ydif*rn2 < 0.) rsign = -1.;
         xpt       = xmid +rsign*dist*rn1;
         ypt       = ymid +rsign*dist*rn2;
#else
         /* MIDPOINT RULE (VERY SIMPLE) */
         xpt = 0.5*(vrtx[v0][0] +vrtx[v1][0]);
         ypt = 0.5*(vrtx[v0][1] +vrtx[v1][1]);
#endif

INSRT:
         /* INSERT POINT */
         vrtx[nvrtx][0] = xpt;
         vrtx[nvrtx][1] = ypt;
         
         qtree.addpt(nvrtx);
         qtree.nearpt(nvrtx,vnear);
         tind = findtri(xpt,ypt,vnear);
         getwgts(wt);
         err = insert(tind,nvrtx,10.0);
         if (!err) {
            vlngth[nvrtx] = 0.0;
            for(i=0;i<3;++i) 
               vlngth[nvrtx] += wt[i]*vlngth[tvrtx[tind][i]];
            nsnew = nsdel +3;
            ntnew = ntdel +2;
            ++intrcnt;
         }
         else {
            /* REBAY ALGORITHM IS TRYING TO INSERT POINT OUTSIDE OF DOMAIN OR CLOSE TO BOUNDARY SKIP SIDE FOR NOW? */
            printf("#Warning: skipping insert\n");
            qtree.dltpt(nvrtx);
            --nvrtx;
            nsdel = 1;
            sdel[0] = sind;
            nsnew = 0;
            ntnew = 2;
            tdel[0] = stri[sind][0];
            tdel[1] = stri[sind][1];
         }
      }
      ++nvrtx;

      for(i=0;i<nsdel;++i) 
         if (intwk3[sdel[i]] > -1) tkoutlst(sdel[i]);
         
      for(i=0;i<nsnew;++i) {
         sind = sdel[i];
         sinfo[sind] = -2; // MARK SIDE AS TOUCHED (OR SHOULD THIS BE DONE IN INSERT?)
         intwk3[sind] = -1;
         fltwkreb(sind);
         if (fltwk[sind] < tolsize) putinlst(sind);
      }
         
      /* RECONSTRUCT AFFECTED TRIS vinfo ARRAY */
      for(i=0;i<ntnew;++i) {
         tind = tdel[i];
         vinfo[tind] = 0;
         for(j=0;j<3;++j)
            if (intwk3[tside[tind].side[j]] > -1) ++vinfo[tind];
      }
   }
   
   for (i=0;i<nsbd;++i)
      bdrysidereorder(i);
      
   cnt_nbor();
   
   printf("#Rebay finished: new interior points %d, new boundary points %d\n",intrcnt,bdrycnt);

   return;
}
   
void mesh::putinlst(int sind) {
   int i, temp, top, bot, mid;
      
   /* CREATE ORDERED LIST OF SIDES BY RATIO FLTWORK (LENGTH/DENSITY) */
   bot = 0;
   if (nslst > 0) {
      top = 0;
      bot = nslst-1;
      if (fltwk[sind] < fltwk[intwk2[top]]) {
         bot = 0;
      }
      else if (fltwk[sind] > fltwk[intwk2[bot]]) {
         bot = nslst;
      }
      else {
         while(top < bot-1) {
            mid = top + (bot -top)/2;
            if (fltwk[sind] > fltwk[intwk2[mid]])
               top = mid;
            else
               bot = mid;
         }
      }
      for(i=nslst-1;i>=bot;--i) {
         temp = intwk2[i];
         intwk2[i+1] = temp;
         intwk3[temp] = i+1;
      }
   }
   intwk2[bot]= sind;
   intwk3[sind] = bot;
   ++nslst;

   assert(nslst < maxvst -1);
   
   return;
}

void mesh::tkoutlst(int sind) {
   int bgn,temp,i;
   
   bgn = intwk3[sind];
   for(i=bgn+1;i<nslst;++i) {
      temp = intwk2[i];
      intwk2[i-1] = temp;
      intwk3[temp] = i-1;
   }
   intwk3[sind] = -1;
   --nslst;
   intwk2[nslst] = -1;
   
   return;
}

void mesh::fltwkreb(int i) {
   FLT dif, av;
   
   /* CALCULATE SIDE LENGTH RATIO FOR YABER */
   /* HAS TO BE A CONTINUOUS FUNCTION SO COMMUNICATION BDRY'S ARE COARSENED PROPERLY */
   /* OTHERWISE 2 BDRY SIDES CAN HAVE SAME EXACT FLTWK (NOT SURE WHICH TO DO FIRST) */
   /*(THIS ESSENTIALLY TAKES THE MINUMUM) */
   dif = 0.5*(vlngth[svrtx[i][0]] -vlngth[svrtx[i][1]]);
   av = 0.5*(vlngth[svrtx[i][0]] +vlngth[svrtx[i][1]]);
   fltwk[i] = (av +(dif*dif/(0.1*av +fabs(dif))))/distance(svrtx[i][0],svrtx[i][1]);
   
   return;
}

void mesh::fltwkreb() {
   int i;
   
   for(i=0;i<nside;++i)
      fltwkreb(i);
   return;
}
      
   
