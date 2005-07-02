/*
 *  srebay.cpp
 *  mblock
 *
 *  Created by helenbrk on Wed Sep 12 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */
#include "mesh.h"
#include "boundary.h"
#include <utilities.h>
#include <assert.h>
#include <math.h>

#define VERBOSE

#define REBAY

/* MEMORY USAGE */
/* SUMMARY OF VARIABLE USAGE
i1wk = used in findtri for triangle searches
nslst is the numer of sides needing refinement 
i2wk = ordered list of sides requiring modification
i3wk = pointer from side into i2wk list
sdel[] = sides added/deleted during local operation
tdel[] = triangles affected by local operation
td[].info = output from yaber left intact - triangles changed can be determined from sd[].info
sd[].info = output from yaber left intact - -3 deleted sides, -2 touched sides, -1 untouched sides
vd[].info = +n sides of triangle needing coarsening
*/



/* i2wk - LIST OF SIDES TO BE REFINED */
/* i3wk - BACK REFERENCE FROM SIDE INTO LIST */
/* vd.info - TRIANGLE ACCEPTANCE INDICATOR >= 0 */

/* nslst KEEPS TRACK OF HOW MANY SIDES NEED REFINEMENT */
/* SIDES ARE STORED IN i2wk WITH BACK REFERENCE IN i3wk */


int nslst;

void mesh::rebay(FLT tolsize) {
   int i,j,n,tind,sind,v0,v1,v2,vnear,nsnew,ntnew,snum,bdrycnt,intrcnt,err;
   int nsdel,sdel[maxlst];
   int ntdel,tdel[maxlst];
   int bnum,bel;
   TinyVector<FLT,2> xpt,dx,rn,xmid,xdif;
   TinyVector<FLT,3> wt;
   FLT norm;
   FLT p,q,s1sq,s2sq,rad1,rad2,rs;
   FLT densty,cirrad,arg,dist,rsign;
      
   /* SET UP FLTWK */
   fltwkreb();
         
   bdrycnt = 0;
   intrcnt = 0;
   
   /* USE VINFO[TIND] = 3,2,1,0 TO FIND BOUNDARY REFINE TRI'S */
   for(i=0;i<ntri;++i)
      vd(i).info = 0;
   
   /* CLASSIFY SIDES AS ACCEPTED OR UNACCEPTED */
   nslst = 0;
   for(i=0;i<nside;++i) {
      if (fscr1(i) < tolsize) {
         vd(sd(i).tri(0)).info += 1;
         vd(MAX(sd(i).tri(1),-1)).info += 1;
         putinlst(i);
      }
   }
   
   vd(-1).info = 0;

   /* BEGIN REFINEMENT ALGORITHM */
   while (nslst > 0) {
      sind = -2;
      for(i=0;i<nslst;++i) {
         if (vd(sd(i2wk(i)).tri(0)).info +vd(MAX(sd(i2wk(i)).tri(1),-1)).info == 6) continue;
         sind = i2wk(i);
         break;
      }

      tind = sd(sind).tri(0);
      for(j=0;j<3;++j) {
         if (td(tind).side(j) == sind) {
            snum = j;
            break;
         }
      }         
      
      v0 = td(tind).vrtx((snum+1)%3);
      v1 = td(tind).vrtx((snum+2)%3);
      v2 = td(tind).vrtx(snum);
      
      if (nvrtx > maxvst -2) {
         *log << "too many vertices" << std::endl;
         exit(1);
      }
      
      
      /* CHECK IF IT IS A BOUNDARY EDGE */
      if (sd(sind).tri(1) < 0) {
         bnum = (-sd(sind).tri(1)>>16) -1;
         bel = -sd(sind).tri(1)&0xFFFF;
         sbdry(bnum)->mvpttobdry(bel,0.5,vrtx(nvrtx));
                  
         /* INSERT POINT */
         qtree.addpt(nvrtx);
         vlngth(nvrtx) = 0.5*(vlngth(v0) +vlngth(v1));
#ifdef VERBOSE
         *log << "Inserting boundary side at " << vrtx(nvrtx)(0) << ' ' << vrtx(nvrtx)(1) << std::endl;
#endif
         bdry_insert(tind,snum,nvrtx,ntdel,tdel,nsdel,sdel);
         ++bdrycnt;
         
         nsnew = nsdel +2;
         ntnew = ntdel +1;
      }
      else {

#ifdef REBAY
         /* THIS IS TRIANGLE BASED REFINEMENT */
         /*	FIND ACCEPTED EDGE ON TRIANGLE & BUILD OFF OF IT */
         if (vd(tind).info < 3) {
            for (snum=0;snum<3;++snum) 
               if (i3wk(td(tind).side(snum)) == -1) break;
         }
         else {
            tind = sd(sind).tri(1);
            for (snum=0;snum<3;++snum) 
               if (i3wk(td(tind).side(snum)) == -1) break;
         }
         
         v0 = td(tind).vrtx((snum+1)%3);
         v1 = td(tind).vrtx((snum+2)%3);
         v2 = td(tind).vrtx(snum);
         
         /* USE REBAY'S ALGORITHM FOR INSERT POINT */
         p = 0.0;
         for(n=0;n<ND;++n) {
            dx(n) = vrtx(v0)(n) -vrtx(v1)(n);
            p += pow(dx(n),2);
         }
         p = 0.5*sqrt(p);
         
         s1sq = 0.0;
         for(n=0;n<ND;++n) {
            dx(n) = vrtx(v2)(n) -vrtx(v1)(n);
            s1sq += pow(dx(n),2);
         }
         s1sq *= 0.25;

         s2sq = 0.0;
         for(n=0;n<ND;++n) {
            dx[n] = vrtx(v2)(n) -vrtx(v0)(n);
            s2sq += pow(dx(n),2);
         }
         s2sq *= 0.25;         
         tcenter(tind,xpt);
         if (p*p > s1sq +s2sq) goto INSRT;
         
         q = 0.0;
         for(n=0;n<ND;++n) {
            xmid(n) = .5*(vrtx(v0)(n) +vrtx(v1)(n));
            dx(n) = xpt(n) -xmid(n);
            q += pow(dx(n),2);
         }
         q = sqrt(q);
         if (q < p) goto INSRT;

         densty = (vlngth(v0)  +vlngth(v1))/sqrt(3.0);
         rad1 = MAX(densty,p);
         rad2 = .5*(p*p  +q*q)/q;
         cirrad = MIN(rad1,rad2);
         arg = fabs(cirrad*cirrad  -p*p);
         dist = cirrad  +sqrt(arg);
         rs = 0.0;
         for(n=0;n<ND;++n) {
            dx(n) = vrtx(v1)(n) -vrtx(v0)(n);
            rs += pow(dx(n),2);
         }
         rs = 1./sqrt(rs);
         rn[0] =  dx[1]*rs;
         rn[1] = -dx[0]*rs;
         for(n=0;n<ND;++n)
            xdif[n] = vrtx(v2)(n)  -xmid[n];
         rsign = 1.;
         if (xdif[0]*rn[0] +xdif[1]*rn[1] < 0.) rsign = -1.;
         for(n=0;n<ND;++n)
            xpt(n) = xmid(n) +rsign*dist*rn(n);
#else
         /* MIDPOINT RULE (VERY SIMPLE) */
         for(n=0;n<ND;++n)
            xpt(n) = 0.5*(vrtx(v0)(n) +vrtx(v1)(n));
#endif

INSRT:
         /* INSERT POINT */
         for(n=0;n<ND;++n)
            vrtx(nvrtx)(n) = xpt(n);
            
#ifdef VERBOSE
         *log << "Inserting interior side ";
         for(n=0;n<ND;++n)
            *log << vrtx(nvrtx)(n) << ' ';
         *log << std::endl;
#endif
         
         dist = qtree.nearpt(vrtx(nvrtx).data(),vnear);
         norm = 0.0;
         for (n=0;n<ND;++n)
            norm += fabs(vrtx(nvrtx)(n));
         if (dist < 100.0*EPSILON*norm) {
            *log << "#Point to close to insert " << dist << ' ';
            for(n=0;n<ND;++n)
               *log << vrtx(nvrtx)(n) << ' ';
            *log << std::endl;
            --nvrtx;
            nsdel = 1;
            sdel[0] = sind;
            nsnew = 0;
            ntnew = 2;
            tdel[0] = sd(sind).tri(0);
            tdel[1] = sd(sind).tri(1);
            goto IEND;
         }
            
         qtree.addpt(nvrtx);
         tind = findtri(xpt,vnear);
         getwgts(wt);
         vlngth(nvrtx) = 0.0;
         for(i=0;i<3;++i)
            vlngth(nvrtx) += wt(i)*vlngth(td(tind).vrtx(i));

         err = insert(tind,nvrtx,0.001,ntdel,tdel,nsdel,sdel);
         if (!err) {
            nsnew = nsdel +3;
            ntnew = ntdel +2;
            ++intrcnt;
         }
         else {
            *log << "#Warning: Makes Bad Triangle ";
            for(n=0;n<ND;++n)
               *log << vrtx(v0)(n) << ' ';
            for(n=0;n<ND;++n)
               *log << vrtx(v1)(n) << ' ';
            *log << std::endl; 
            qtree.dltpt(nvrtx);
            --nvrtx;
            nsdel = 1;
            sdel[0] = sind;
            nsnew = 0;
            ntnew = 2;
            tdel[0] = sd(sind).tri(0);
            tdel[1] = sd(sind).tri(1);
         }
      }
IEND: ++nvrtx;
         
      for(i=0;i<nsdel;++i) 
         if (i3wk(sdel[i]) > -1) tkoutlst(sdel[i]);
         
      if (i3wk(sind) > -1) tkoutlst(sind);
                  
      for(i=0;i<nsnew;++i) {
         sind = sdel[i];
         sd(sind).info = -2; // MARK SIDE AS TOUCHED (OR SHOULD THIS BE DONE IN INSERT?)
         i3wk(sind) = -1;
         fltwkreb(sind);
         if (fscr1(sind) < tolsize) putinlst(sind);
      }
         
      /* RECONSTRUCT AFFECTED TRIS vinfo ARRAY */
      for(i=0;i<ntnew;++i) {
         tind = tdel[i];
         vd(tind).info = 0;
         for(j=0;j<3;++j)
            if (i3wk(td(tind).side(j)) > -1) ++vd(tind).info;
      }
   }
   
   for (i=0;i<nsbd;++i) {
      sbdry(i)->reorder();
      sbdry(i)->setupcoordinates();
   }
   
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
      if (fscr1(sind) < fscr1(i2wk(top))) {
         bot = 0;
      }
      else if (fscr1(sind) > fscr1(i2wk(bot))) {
         bot = nslst;
      }
      else {
         while(top < bot-1) {
            mid = top + (bot -top)/2;
            if (fscr1(sind) > fscr1(i2wk(mid)))
               top = mid;
            else
               bot = mid;
         }
      }
      for(i=nslst-1;i>=bot;--i) {
         temp = i2wk(i);
         i2wk(i+1) = temp;
         i3wk(temp) = i+1;
      }
   }
   i2wk(bot)= sind;
   i3wk(sind) = bot;
   ++nslst;

   assert(nslst < maxvst -1);
   
   return;
}

void mesh::tkoutlst(int sind) {
   int bgn,temp,i;
   
   bgn = i3wk(sind);
   for(i=bgn+1;i<nslst;++i) {
      temp = i2wk(i);
      i2wk(i-1) = temp;
      i3wk(temp) = i-1;
   }
   i3wk(sind) = -1;
   --nslst;
   i2wk(nslst) = -1;
   
   return;
}

void mesh::fltwkreb(int i) {
   FLT dif, av;
   
   /* CALCULATE SIDE LENGTH RATIO FOR YABER */
   /* HAS TO BE A CONTINUOUS FUNCTION SO COMMUNICATION BDRY'S ARE COARSENED PROPERLY */
   /* OTHERWISE 2 BDRY SIDES CAN HAVE SAME EXACT FLTWK (NOT SURE WHICH TO DO FIRST) */
   /*(THIS ESSENTIALLY TAKES THE MINUMUM) */
   dif = 0.5*(vlngth(sd(i).vrtx(0)) -vlngth(sd(i).vrtx(1)));
   av = 0.5*(vlngth(sd(i).vrtx(0)) +vlngth(sd(i).vrtx(1)));
   fscr1(i) = (av -(dif*dif/(0.01*av +fabs(dif))))/distance(sd(i).vrtx(0),sd(i).vrtx(1));
   //fscr1(i) = MIN(vlngth(sd(i).vrtx(0)),vlngth(sd(i).vrtx(1)))/distance(sd(i).vrtx(0),sd(i).vrtx(1));

   return;
}

void mesh::fltwkreb() {
   int i;
   
   for(i=0;i<nside;++i)
      fltwkreb(i);
   return;
}

void mesh::trebay(FLT tolsize) {
   int i,j,n,tind,tfind,sind,v0,v1,v2,vnear,nsnew,ntnew,snum,bdrycnt,intrcnt,err;
   TinyVector<FLT,ND> xpt,dx,xmid,xdif,rn;
   TinyVector<FLT,3> wt;
   FLT norm,p,q,s1sq,s2sq,rad1,rad2,rs;
   FLT densty,cirrad,arg,dist,rsign;
   int ntdel,tdel[maxlst];
   int nsdel,sdel[maxlst];
   
   /* SET UP FLTWK */
   fltwkreb();
         
   bdrycnt = 0;
   intrcnt = 0;
   
   /* USE VINFO[TIND] = 3,2,1,0 TO FIND BOUNDARY REFINE TRI'S */
   for(i=0;i<ntri;++i)
      vd(i).info = 0;
   
   /* CLASSIFY TRIANGLES AS ACCEPTED OR UNACCEPTED */
   nslst = 0;
   for(i=0;i<ntri;++i) {
      for(j=0;j<3;++j) {
         if (fscr1(td(i).side(j)) < tolsize) {
            putinlst(i);
            break;
         }
      }
   }
   
   vd(-1).info = 0;

   /* BEGIN REFINEMENT ALGORITHM */
   while (nslst > 0) {
      for(i=0;i<nslst;++i) {
         for(j=0;j<3;++j) {
            tind = td(i2wk(i)).tri(j);
            if (tind < 0 || i3wk(tind) == -1) {
               snum = j;
               tind = i2wk(i);
               goto TFOUND;
            }
         }
      }
            
      TFOUND:
      
      if (nvrtx > maxvst -2) {
         *log << "too many vertices" << std::endl;
         exit(1);
      }
      
      /* THIS IS TRIANGLE BASED REFINEMENT */
      /*	FIND ACCEPTED EDGE ON TRIANGLE & BUILD OFF OF IT */
      v0 = td(tind).vrtx((snum+1)%3);
      v1 = td(tind).vrtx((snum+2)%3);
      v2 = td(tind).vrtx(snum);
      
      /* USE REBAY'S ALGORITHM FOR INSERT POINT */
      p = 0.0;
      for(n=0;n<ND;++n) {
         dx(n) = vrtx(v0)(n) -vrtx(v1)(n);
         p += pow(dx[n],2);
      }
      p = 0.5*sqrt(p);
      
      s1sq = 0.0;
      for(n=0;n<ND;++n) {
         dx(n) = vrtx(v2)(n) -vrtx(v1)(n);
         s1sq += pow(dx(n),2);
      }
      s1sq *= 0.25;

      s2sq = 0.0;
      for(n=0;n<ND;++n) {
         dx[n] = vrtx(v2)(n) -vrtx(v0)(n);
         s2sq += pow(dx(n),2);
      }
      s2sq *= 0.25;         
      tcenter(tind,xpt);
      if (p*p > s1sq +s2sq) goto INSRT;
      
      q = 0.0;
      for(n=0;n<ND;++n) {
         xmid[n] = .5*(vrtx(v0)(n) +vrtx(v1)(n));
         dx[n] = xpt(n) -xmid(n);
         q += pow(dx(n),2);
      }
      q = sqrt(q);
      if (q < p) goto INSRT;

      densty = (vlngth(v0)  +vlngth(v1))/sqrt(3.0);
      rad1 = MAX(densty,p);
      rad2 = .5*(p*p  +q*q)/q;
      cirrad = MIN(rad1,rad2);
      arg = fabs(cirrad*cirrad  -p*p);
      dist = cirrad  +sqrt(arg);
      rs = 0.0;
      for(n=0;n<ND;++n) {
         dx(n) = vrtx(v1)(n) -vrtx(v0)(n);
         rs += pow(dx(n),2);
      }
      rs = 1./sqrt(rs);
      rn(0) =  dx(1)*rs;
      rn(1) = -dx(0)*rs;
      for(n=0;n<ND;++n)
         xdif(n) = vrtx(v2)(n)  -xmid(n);
      rsign = 1.;
      if (xdif(0)*rn(0) +xdif(1)*rn(1) < 0.) rsign = -1.;
      for(n=0;n<ND;++n)
         xpt(n) = xmid(n) +rsign*dist*rn(n);

INSRT:
      /* INSERT POINT */
      for(n=0;n<ND;++n)
         vrtx(nvrtx)(n) = xpt(n);
         
#ifdef VERBOSE
      *log << "Inserting interior side " << intrcnt << ' ';
      for(n=0;n<ND;++n)
         *log << vrtx(nvrtx)(n) << ' ';
      *log << std::endl;
#endif
      
      dist = qtree.nearpt(vrtx(nvrtx).data(),vnear);
      norm = 0.0;
      for (n=0;n<ND;++n)
         norm += fabs(vrtx(nvrtx)(n));
      if (dist < 100.0*EPSILON*norm) {
         *log << "#Point to close to insert " << dist << ' ';
         for(n=0;n<ND;++n)
            *log << vrtx(vnear)(n) << ' ';
         *log << std::endl;
         tkoutlst(tind);
         continue;
      }
         
      tfind = findtri(xpt,vnear);
      if (tfind < 0) {
         *log << "#Warning: Trying to insert outd domain ";
         for(n=0;n<ND;++n)
            *log << vrtx(v0)(n) << ' '; 
         for(n=0;n<ND;++n)
            *log << vrtx(v1)(n) << ' ';
         *log << std::endl;
         tkoutlst(tind);
         continue;
      }

      getwgts(wt);
      vlngth(nvrtx) = 0.0;
      for(i=0;i<3;++i)
         vlngth(nvrtx) += wt(i)*vlngth(td(tfind).vrtx(i));

      err = insert(tfind,nvrtx,0.00001,ntdel,tdel,nsdel,sdel);
      if (!err) {
         qtree.addpt(nvrtx);
         nsnew = nsdel +3;
         ntnew = ntdel +2;
         ++intrcnt;
      }
      else {
         *log << "#Warning: Makes Bad Triangle ";
         for(n=0;n<ND;++n)
            *log << vrtx(v0)(n) << ' ';
         for(n=0;n<ND;++n)
            *log << vrtx(v1)(n) << ' ';
         *log << std::endl;
         tkoutlst(tind);
         continue;
      }
      ++nvrtx;
         
      for(i=0;i<ntdel;++i) 
         if (i3wk(tdel[i]) > -1) tkoutlst(tdel[i]);
      
      if (i3wk(tind) > -1) tkoutlst(tind);
      
      for(i=0;i<nsnew;++i) {
         sind = sdel[i];
         sd(sind).info = -2; // MARK SIDE AS TOUCHED (OR SHOULD THIS BE DONE IN INSERT?)
         fltwkreb(sind);
      }
         
      for(i=0;i<ntnew;++i) {
         tind = tdel[i];
         for(j=0;j<3;++j) {
            if (fscr1(td(tind).side(j)) < tolsize) {
               putinlst(tind);
               break;
            }
         }
      }
//      checkintegrity();
   }
   
   for (i=0;i<nsbd;++i) {
      sbdry(i)->reorder();
      sbdry(i)->setupcoordinates();
   }
   
   bdrylabel();
      
   cnt_nbor();
      
   *log << "#Rebay finished: new interior points " << intrcnt << " new boundary points " << bdrycnt << std::endl;

   return;
}


/* STEPS */
// INITIALIZE FLAGS //
// COARSEN BOUNDARIES //
// COARSEN INTERIOR //
// REFINE BOUNDARIES //
// REFINE INTERIOR //




      
   
