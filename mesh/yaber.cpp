/*
 *  yabers.cpp
 *  mblock
 *
 *  Created by helenbrk on Fri Sep 14 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#include"mesh.h"
#include"utilities.h"
#include<assert.h>

/* THIS IS SUPPOSED TO DO THE REVERSE OF THE REBAY ROUTINE I HOPE */
/* THUS THE NAME YABER -> REBAY */

/* THESE TELL WHICH SIDES/TRIS WERE DELETED */
/* ntdel, tdel[MAXLST+1]; */
/* nsdel, sdel[MAXLST+1]; */

/* nslst IS THE NUMER OF SIDES NEEDING COARSENING */
/* STORED IN INTWK2 WITH BACK REFERENCE IN INTWK3 */

void mesh::yaber(FLT tolsize, int yes_swap, FLT swaptol) {
   int i,j,tind,sind,nfail,v0,v1,cnt;
   
   /* SET UP FLTWK */
   fltwkyab();
   
   cnt = 0;
   
   /* NEED TO INITIALIZE TO ZERO TO KEEP TRACK OF DELETED TRIS (-1) */
   /* ALSO TO DETERMINE TRI'S ON BOUNDARY OF COARSENING REGION */
   for(i=0;i<ntri;++i)
      tinfo[i] = 0;

   /* KEEPS TRACK OF DELETED SIDES = -3, TOUCHED SIDES =-2, UNTOUCHED SIDES =-1 */
   for(i=0;i<nside;++i)
      sinfo[i] = -1;
      
   /* SWAP SIDES AND MARK SWAPPED SIDES TOUCHED */
   if (yes_swap) swap(swaptol);
      
   /* VINFO TO KEEP TRACK OF SPECIAL VERTICES (1) DELETED VERTICES (-1) */
   for(i=0;i<nvrtx;++i)
      vinfo[i] = 0;
      
   /* MARK BEGINNING/END OF SIDE GROUPS & SPECIAL VERTEX POINTS */
   /* THESE SHOULD NOT BE DELETED */
   for(i=0;i<nsbd;++i) {
      v0 = svrtx[sbdry[i]->sd(0)][0];
      v1 = svrtx[sbdry[i]->sd(sbdry[i]->nsd()-1)][1];
      if (v0 != v1) {
         vinfo[v0] = 1;
         vinfo[v1] = 1;
      }
   }
   
   for(i=0;i<nvbd;++i)
      vinfo[vbdry[i]->v()] = 1;
   
   /* CLASSIFY SIDES AS ACCEPTED OR UNACCEPTED */
   nslst = 0;
   for(i=0;i<nside;++i) {
      if (fltwk[i] > tolsize) {
         tinfo[stri[i][0]] += 1;
         tinfo[MAX(stri[i][1],-1)] += 1;
         putinlst(i);
      }
   }

   tinfo[-1] = 0;
   
   /* BEGIN COARSENING ALGORITHM */
   while (nslst > 0) {
      sind = -1;
      for(i=nslst-1;i>=0;--i) {  // START WITH LARGEST SIDE TO DENSITY RATIO
         if (tinfo[stri[intwk2[i]][0]] +tinfo[MAX(stri[intwk2[i]][1],-1)] == 6) continue;
         sind = intwk2[i];
         break;
      }
      assert(sind != -1);
      
      /* COLLAPSE EDGE */
      nfail = collapse(sind);
      ++cnt;
      
      for(i=0;i<nsdel;++i) 
         if (intwk3[sdel[i]] > -1) tkoutlst(sdel[i]);
         
      if (nfail) {
         *log << "#Warning: side collapse failed" << sind << std::endl;
         /* MARK SIDE AS ACCEPTED AND MOVE ON */
         for(i=0;i<2;++i) {
            tind = stri[sind][i];
            if (tind < 0) continue;
            if (tinfo[tind] > 0) --tinfo[tind];
         }
         continue;
      }
         
      /* RECONSTRUCT AFFECTED TINFO ARRAY */
      for(i=0;i<ntdel;++i) {
         tind = tdel[i];
         tinfo[tind] = 0;
         for(j=0;j<3;++j) {
            sind = tside[tind].side[j];
            if (intwk3[sind] > -1) tkoutlst(sind);
            fltwkyab(sind);
            if (fltwk[sind] > tolsize) {
               putinlst(sind);
               ++tinfo[tind];
            }
         }
      }
   }
   
   /* DELETE LEFTOVER VERTICES */
   /* VINFO > NVRTX STORES VRTX MOVEMENT HISTORY */
   for(i=nvrtx-1;i>=0;--i) {
      if (vinfo[i] < 0) {
         vinfo[nvrtx-1] = i;
         dltvrtx(i);
      }
   }
   
   /* FIX BOUNDARY CONDITION POINTERS */
   for(i=0;i<nvbd;++i)
      while (vbdry[i]->v() >= nvrtx) 
         vbdry[i]->v() = vinfo[vbdry[i]->v()];  
   
   /* DELETE SIDES FROM BOUNDARY CONDITIONS */
   for(i=0;i<nsbd;++i)
      for(j=sbdry[i]->nsd()-1;j>=0;--j) 
         if (sinfo[sbdry[i]->sd(j)] == -3) 
            sbdry[i]->sd(j) = sbdry[i]->sd(--sbdry[i]->nsd());
                        
   /* CLEAN UP SIDES */
   /* SINFO WILL END UP STORING -1 UNTOUCHED, -2 TOUCHED, or INITIAL INDEX OF UNTOUCHED SIDE */
   /* SINFO > NSIDE WILL STORE MOVEMENT HISTORY */
   for(i=nside-1;i>=0;--i) {
      if (sinfo[i] == -3) {
         if (sinfo[nside-1] >= -1)
            sinfo[i] = MAX(nside-1,sinfo[nside-1]);
         else 
            sinfo[i] = -2;
         sinfo[nside-1] = i;
         dltside(i);
      }
   }
   
   /* FIX BOUNDARY CONDITION POINTERS */
   for(i=0;i<nsbd;++i)
      for(j=0;j<sbdry[i]->nsd();++j) 
         while (sbdry[i]->sd(j) >= nside) 
            sbdry[i]->sd(j) = sinfo[sbdry[i]->sd(j)]; 

   for (i=0;i<nsbd;++i)
      sbdry[i]->reorder();
      
   bdrylabel();
   
   /* CLEAN UP DELETED TRIS */
   /* TINFO < NTRI STORES INDEX OF ORIGINAL TRI ( > 0), TINFO = 0 -> UNMOVED */
   /* TINFO > NTRI STORES TRI MOVEMENT HISTORY */
   for(i=0;i<ntri;++i)
      assert(tinfo[i] <= 0);  
      
   for(i=ntri-1;i>=0;--i) {
      if (tinfo[i] < 0) {  // DELETED TRI
         tinfo[i] = MAX(ntri-1,tinfo[ntri-1]);  // TINFO STORES ORIGINAL TRI INDEX
         tinfo[ntri-1] = i; // RECORD MOVEMENT HISTORY FOR TRI'S > NTRI
         dlttri(i);
      }
   }
   
   *log << "#Yaber finished: " << cnt << " sides coarsened" << std::endl;

   return;
}

void mesh::checkintegrity() const {
   int i,j,sind,dir;
   
   for(i=0;i<ntri;++i) {
      if (tinfo[i] < 0) continue;
      
      if (area(i) < 0.0) *log << "negative area" << i << std::endl;
      
      for(j=0;j<3;++j) {
         sind = tside[i].side[j];
         dir = -(tside[i].sign[j] -1)/2;
         
         if (sinfo[sind] == -3) {
            *log << "references deleted side" <<  i << sind << std::endl;
            for(i=0;i<nside;++i)
               sinfo[i] += 2;
            out_mesh("error");
            exit(1);;
         }

         if (svrtx[sind][dir] != tvrtx[i][(j+1)%3] && svrtx[sind][1-dir] != tvrtx[i][(j+2)%3]) {
            *log << "failed vrtx check tind" << i << "sind" << sind << std::endl;
            for(i=0;i<nside;++i)
               sinfo[i] += 2;
            out_mesh("error"); 
            exit(1);;
         }    
         
         if (stri[sind][dir] != i) {
            *log << "failed side check tind" << i << "sind" << sind << std::endl;
            for(i=0;i<nside;++i)
               sinfo[i] += 2;
            out_mesh("error"); 
            exit(1);;
         }
         
         if (ttri[i][j] != stri[sind][1-dir]) {
            *log << "failed ttri check tind" << i << "sind" << sind << std::endl;
            for(i=0;i<nside;++i)
               sinfo[i] += 2;
            out_mesh("error"); 
            exit(1);;
         }
         
         if (ttri[i][j] > 0) {
            if(tinfo[ttri[i][j]] < 0) {
               *log << "references deleted tri" << std::endl;
               for(i=0;i<nside;++i)
                  sinfo[i] += 2;
               out_mesh("error"); 
               exit(1);;
            }
         }
      }
   }
   
   return;
}

void mesh::checkintwk() const {
   int i;
   
   for(i=0;i<maxvst;++i)
      if (intwk1[i] != -1) *log << "failed intwk1 check" << i << intwk1[i] << std::endl;
      
   for(i=0;i<maxvst;++i)
      if (intwk2[i] != -1) *log << "failed intwk2 check" << i << intwk2[i] << std::endl;
   
   for(i=0;i<maxvst;++i)
      if (intwk3[i] != -1) *log << "failed intwk3 check" << i << intwk3[i] << std::endl;
   
   return;
}

void mesh::fltwkyab(int i) {
   FLT dif,av;
   
   /* CALCULATE SIDE LENGTH RATIO FOR YABER */
   /* HAS TO BE A CONTINUOUS FUNCTION SO COMMUNICATION BDRY'S ARE COARSENED PROPERLY */
   /* OTHERWISE 2 BDRY SIDES CAN HAVE SAME EXACT FLTWK (NOT SURE WHICH TO DO FIRST) */
   /*(THIS ESSENTIALLY TAKES THE MINUMUM) */
   dif = 0.5*(vlngth[svrtx[i][0]] -vlngth[svrtx[i][1]]);
   av = 0.5*(vlngth[svrtx[i][0]] +vlngth[svrtx[i][1]]);
   fltwk[i] = (av -(dif*dif/(0.1*av +fabs(dif))))/distance(svrtx[i][0],svrtx[i][1]);
   
   return;
}

void mesh::fltwkyab() {
   int i;
   
   for(i=0;i<nside;++i)
      fltwkyab(i);
   return;
}

