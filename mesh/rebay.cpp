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

/* nslst KEEPS TRACK OF HOW MANY SIDES NEED REFINEMENT */
/* SIDES ARE STORED IN INTWK2 WITH BACK REFERENCE IN INTWK3 */

/* THESE TELL WHICH SIDES/TRIS WERE DELETED */
/* ntdel, tdel[MAXLST+1]; */
/* nsdel, sdel[MAXLST+1]; */


void mesh::rebay(FLT tolsize) {
   int i,j,n,tind,sind,v0,v1,v2,vnear,nsnew,ntnew,snum,bdrycnt,intrcnt,err;
   int bnum,bel;
   FLT xpt[2],wt[3];
   FLT dx[2],p,q,s1sq,s2sq,rad1,rad2,rs;
   FLT xmid[2],densty,cirrad,arg,dist,rn[2],xdif[2],rsign;
   
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
         *log << "too many vertices" << std::endl;
         exit(1);
      }
      
      
      /* CHECK IF IT IS A BOUNDARY EDGE */
      if (stri[sind][1] < 0) {
         bnum = (-stri[sind][1]>>16) -1;
         bel = -stri[sind][1]&0xFFFF;
         sbdry[bnum]->mvpttobdry(bel,0.5,vrtx[nvrtx]);
                  
         /* INSERT POINT */
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
         p = 0.0;
         for(n=0;n<ND;++n) {
            dx[n] = vrtx[v0][n] -vrtx[v1][n];
            p += pow(dx[n],2);
         }
         p = 0.5*sqrt(p);
         
         s1sq = 0.0;
         for(n=0;n<ND;++n) {
            dx[n] = vrtx[v2][n] -vrtx[v1][n];
            s1sq += pow(dx[n],2);
         }
         s1sq *= 0.25;

         s2sq = 0.0;
         for(n=0;n<ND;++n) {
            dx[n] = vrtx[v2][n] -vrtx[v0][n];
            s2sq += pow(dx[n],2);
         }
         s2sq *= 0.25;         
         tcenter(tind,xpt);
         if (p*p > s1sq +s2sq) goto INSRT;
         
         q = 0.0;
         for(n=0;n<ND;++n) {
            xmid[n] = .5*(vrtx[v0][n] +vrtx[v1][n]);
            dx[n] = xpt[n] -xmid[n];
            q += pow(dx[n],2);
         }
         q = sqrt(q);
         if (q < p) goto INSRT;

         densty = (vlngth[v0]  +vlngth[v1])/sqrt(3.0);
         rad1 = MAX(densty,p);
         rad2 = .5*(p*p  +q*q)/q;
         cirrad = MIN(rad1,rad2);
         arg = fabs(cirrad*cirrad  -p*p);
         dist = cirrad  +sqrt(arg);
         rs = 0.0;
         for(n=0;n<ND;++n) {
            dx[n] = vrtx[v1][n] -vrtx[v0][n];
            rs += pow(dx[n],2);
         }
         rs = 1./sqrt(rs);
         rn[0] =  dx[1]*rs;
         rn[1] = -dx[0]*rs;
         for(n=0;n<ND;++n)
            xdif[n] = vrtx[v2][n]  -xmid[n];
         rsign = 1.;
         if (xdif[0]*rn[0] +xdif[1]*rn[1] < 0.) rsign = -1.;
         for(n=0;n<ND;++n)
            xpt[n] = xmid[n] +rsign*dist*rn[n];
#else
         /* MIDPOINT RULE (VERY SIMPLE) */
         for(n=0;n<ND;++n)
            xpt[n] = 0.5*(vrtx[v0][n] +vrtx[v1][n]);
#endif

INSRT:
         /* INSERT POINT */
         for(n=0;n<ND;++n)
            vrtx[nvrtx][n] = xpt[n];
         
         qtree.addpt(nvrtx);
         qtree.nearpt(nvrtx,vnear);
         tind = findtri(xpt,vnear);
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
            *log << "#Warning: skipping insert" << std::endl;
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
      sbdry[i]->reorder();
   
   bdrylabel();
      
   cnt_nbor();
   
   *log << "#Rebay finished: new interior points " << intrcnt << " new boundary points " << bdrycnt << std::endl;

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
      
   
