/*
 *  yabers.cpp
 *  mblock
 *
 *  Created by helenbrk on Fri Sep 14 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#include "mesh.h"
#include "boundary.h"
#include <utilities.h>
#include <assert.h>

/* THIS IS SUPPOSED TO DO THE REVERSE OF THE REBAY ROUTINE I HOPE */
/* THUS THE NAME YABER -> REBAY */

extern int nslst;  // (THIS NEEDS TO BE FIXED)

/* SUMMARY OF VARIABLE USAGE
i1wk = not used!
nslst is the numer of sides needing coarsening 
i2wk = ordered list of sides requiring modification
i3wk = pointer from side into i2wk list
sdel[] = sides deleted during local operation
tdel[] = triangles affected by local operation
td[].info = -3 deleted tri, -1 no refinement needed, +n sides of triangle needing coarsening
td[].info = i < ntri stores index of original tri or tinfo = -2 touched, tinfo = -1 for untouched, i > ntri stores movment location
sd[].info = during process: -3 deleted sides, -2 touched sides, -1 untouched sides
sd[].info = after clean-up: i < nside -2, -1 or index of original, i > nside movement location
vd[].info = during process: -3 deleted vrtx, -2 special vrtx, -1 untouched vertex
vd[].info = after cleanup: i < nvrtx: > 0 original location, i > nvrtx movement location
*/

void mesh::yaber(FLT tolsize, int yes_swap, FLT swaptol) {
   int i,j,tind,sind,nfail,v0,v1,cnt;
   int ntdel,tdel[maxlst];
   int nsdel,sdel[maxlst];

   /* SET UP FLTWK */
   fltwkyab();
   
   cnt = 0;
   
   /* NEED TO INITIALIZE TO ZERO TO KEEP TRACK OF DELETED TRIS (-3) */
   /* ALSO TO DETERMINE TRI'S ON BOUNDARY OF COARSENING REGION */
   for(i=0;i<ntri;++i)
      td(i).info = -1;

   /* KEEPS TRACK OF DELETED SIDES = -3, TOUCHED SIDES =-2, UNTOUCHED SIDES =-1 */
   for(i=0;i<nside;++i)
      sd(i).info = -1;
      
   /* VINFO TO KEEP TRACK OF UNTOUCHED (-1), SPECIAL VERTICES (-2), DELETED VERTICES (-3) */
   for(i=0;i<nvrtx;++i)
      vd(i).info = -1;
      
   /* MARK BEGINNING/END OF SIDE GROUPS & SPECIAL VERTEX POINTS */
   /* THESE SHOULD NOT BE DELETED */
   for(i=0;i<nsbd;++i) {
      v0 = sd(sbdry(i)->el(0)).vrtx(0);
      v1 = sd(sbdry(i)->el(sbdry(i)->nel-1)).vrtx(1);
      vd(v0).info = -2;
      vd(v1).info = -2;
   }
   
   for(i=0;i<nvbd;++i)
      vd(vbdry(i)->v0).info = -2;
      
   /* SWAP SIDES AND MARK SWAPPED SIDES TOUCHED */
   if (yes_swap) swap(swaptol);
   
   /* CLASSIFY SIDES AS ACCEPTED OR UNACCEPTED */
   nslst = 0;
   for(i=0;i<nside;++i) {
      if (fscr1(i) > tolsize) {
         td(sd(i).tri(0)).info += 1;
         td(MAX(sd(i).tri(1),-1)).info += 1;
         putinlst(i);
      }
   }

   td(-1).info = -1;
   
   /* BEGIN COARSENING ALGORITHM */
   while (nslst > 0) {
      sind = -1;
      for(i=nslst-1;i>=0;--i) {  // START WITH LARGEST SIDE TO DENSITY RATIO
         if (td(sd(i2wk(i)).tri(0)).info +td(MAX(sd(i2wk(i)).tri(1),-1)).info == 4) continue;
         sind = i2wk(i);
         break;
      }
      assert(sind != -1);
   
      /* COLLAPSE EDGE */
      nfail = collapse(sind,ntdel,tdel,nsdel,sdel);
      ++cnt;
            
      for(i=0;i<nsdel;++i) {
         if (i3wk(sdel[i]) > -1) tkoutlst(sdel[i]);
      }
         
      if (nfail < 0) {
         *log << "#Warning: side collapse failed" << sind << std::endl;
         /* MARK SIDE AS ACCEPTED AND MOVE ON */
         for(i=0;i<2;++i) {
            tind = sd(sind).tri(i);
            if (tind < 0) continue;
            if (td(tind).info > -1) --td(tind).info;
         }
         continue;
      }
         
      /* RECONSTRUCT AFFECTED TINFO ARRAY */
      for(i=0;i<ntdel;++i) {
         tind = tdel[i];
         td(tind).info = -1;
         for(j=0;j<3;++j) {
            sind = td(tind).side(j);
            if (i3wk(sind) > -1) tkoutlst(sind);
            fltwkyab(sind);
            if (fscr1(sind) > tolsize) {
               putinlst(sind);
               ++td(tind).info;
            }
         }
      }
   }
   
   /* DELETE LEFTOVER VERTICES */
   /* VINFO > NVRTX STORES VRTX MOVEMENT HISTORY */
   for(i=0;i<nvrtx;++i) 
      if (vd(i).info == -3) 
         dltvrtx(i);
   
   /* FIX BOUNDARY CONDITION POINTERS */
   for(i=0;i<nvbd;++i)
      if (vbdry(i)->v0 >= nvrtx) 
         vbdry(i)->v0 = vd(vbdry(i)->v0).info;  
   
   /* DELETE SIDES FROM BOUNDARY CONDITIONS */
   for(i=0;i<nsbd;++i)
      for(j=sbdry(i)->nel-1;j>=0;--j) 
         if (sd(sbdry(i)->el(j)).info == -3) 
            sbdry(i)->el(j) = sbdry(i)->el(--sbdry(i)->nel);
                        
   /* CLEAN UP SIDES */
   /* SINFO WILL END UP STORING -1 UNTOUCHED, -2 TOUCHED, or INITIAL INDEX OF UNTOUCHED SIDE */
   /* SINFO > NSIDE WILL STORE MOVEMENT LOCATION */  /* TEMPORARY HAVEN"T TESTED THIS */
   for(i=0;i<nside;++i) 
      if (sd(i).info == -3) dltsd(i);
         
   /* FIX BOUNDARY CONDITION POINTERS */
   for(i=0;i<nsbd;++i)
      for(j=0;j<sbdry(i)->nel;++j) 
         if (sbdry(i)->el(j) >= nside) 
            sbdry(i)->el(j) = sd(sbdry(i)->el(j)).info; 

   for (i=0;i<nsbd;++i) {
      sbdry(i)->reorder();
      sbdry(i)->setupcoordinates();
   }
      
   bdrylabel();
   
   /* CLEAN UP DELETED TRIS */
   /* TINFO < NTRI STORES INDEX OF ORIGINAL TRI ( > 0), TINFO = 0 -> UNMOVED */
   /* TINFO > NTRI STORES TRI MOVEMENT HISTORY */
   for(i=0;i<ntri;++i)
      if (td(i).info == -3)  dlttri(i);
      
   *log << "#Yaber finished: " << cnt << " sides coarsened" << std::endl;

   return;
}

void mesh::checkintegrity() {
   int i,j,sind,dir;
   
   for(i=0;i<ntri;++i) {
      if (td(i).info < 0) continue;
      
      if (area(i) < 0.0) *log << "negative area" << i << std::endl;
      
      for(j=0;j<3;++j) {
         sind = td(i).side(j);
         dir = -(td(i).sign(j) -1)/2;
         
         if (sd(sind).info == -3) {
            *log << "references deleted side" <<  i << sind << std::endl;
            for(i=0;i<nside;++i)
               sd(i).info += 2;
            output("error");
            exit(1);
         }

         if (sd(sind).vrtx(dir) != td(i).vrtx((j+1)%3) && sd(sind).vrtx(1-dir) != td(i).vrtx((j+2)%3)) {
            *log << "failed vrtx check tind" << i << "sind" << sind << std::endl;
            for(i=0;i<nside;++i)
               sd(i).info += 2;
            output("error"); 
            exit(1);
         }    
         
         if (sd(sind).tri(dir) != i) {
            *log << "failed side check tind" << i << "sind" << sind << std::endl;
            for(i=0;i<nside;++i)
               sd(i).info += 2;
            output("error"); 
            exit(1);
         }
         
         if (td(i).tri(j) != sd(sind).tri(1-dir)) {
            *log << "failed ttri check tind" << i << "sind" << sind << std::endl;
            for(i=0;i<nside;++i)
               sd(i).info += 2;
            output("error"); 
            exit(1);
         }
         
         if (td(i).tri(j) > 0) {
            if(td(td(i).tri(j)).info < 0) {
               *log << "references deleted tri" << std::endl;
               for(i=0;i<nside;++i)
                  sd(i).info += 2;
               output("error"); 
               exit(1);
            }
         }
      }
   }
   
   return;
}

void mesh::checkintwk() const {
   int i;
   
   for(i=0;i<maxvst;++i)
      if (i1wk(i) != -1) *log << "failed intwk1 check" << i << i1wk(i) << std::endl;
      
   for(i=0;i<maxvst;++i)
      if (i2wk(i) != -1) *log << "failed intwk2 check" << i << i2wk(i) << std::endl;
   
   for(i=0;i<maxvst;++i)
      if (i3wk(i) != -1) *log << "failed intwk3 check" << i << i3wk(i) << std::endl;
   
   return;
}

void mesh::fltwkyab(int i) {
   FLT dif,av;
   
   /* CALCULATE SIDE LENGTH RATIO FOR YABER */
   /* HAS TO BE A CONTINUOUS FUNCTION SO COMMUNICATION BDRY'S ARE COARSENED PROPERLY */
   /* OTHERWISE 2 BDRY SIDES CAN HAVE SAME EXACT FLTWK (NOT SURE WHICH TO DO FIRST) */
   /*(THIS ESSENTIALLY TAKES THE MINUMUM) */
   dif = 0.5*(vlngth(sd(i).vrtx(0)) -vlngth(sd(i).vrtx(1)));
   av = 0.5*(vlngth(sd(i).vrtx(0)) +vlngth(sd(i).vrtx(1)));
   fscr1(i) = (av -(dif*dif/(0.1*av +fabs(dif))))/distance(sd(i).vrtx(0),sd(i).vrtx(1));
   
   return;
}

void mesh::fltwkyab() {
   int i;
   
   for(i=0;i<nside;++i)
      fltwkyab(i);
   return;
}
      
      
      
   
   

