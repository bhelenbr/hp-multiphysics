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

/*	THIS IS SUPPOSED TO DO THE REVERSE OF THE REBAY ROUTINE I HOPE */
/* THUS THE NAME YABER -> REBAY */

/*	THESE TELL WHICH SIDES/TRIS WERE DELETED */
extern int ntdel, tdel[MAXLST+1];
extern int nsdel, sdel[MAXLST+1];

/*	THIS IS THE NUMER OF SIDES NEEDING COARSENING */
/*	STORED IN INTWK2 WITH BACK REFERENCE IN INTWK3 */
extern int nslst;

void putinlst(int sind);
void tkoutlst(int sind);

void mesh::yaber(FLT tolsize) {
   int i,j,tind,sind,nfail,v0,v1;
   
/*	NEED TO INITIALIZE TO ZERO TO KEEP TRACK OF DELETED TRIS (-1) */
/* ALSO TO DETERMINE TRI'S ON BOUNDARY OF COARSENING REGION */
   for(i=0;i<ntri;++i)
      tinfo[i] = 0;

/* KEEPS TRACK OF DELETED SIDES = -3, TOUCHED SIDES =-2, UNTOUCHED SIDES =-1 */
   for(i=0;i<nside;++i)
      sinfo[i] = -1;
      
/* VINFO TO KEEP TRACK OF SPECIAL VERTICES (1) DELETED VERTICES (-1) */
   for(i=0;i<nvrtx;++i)
      vinfo[i] = 0;
      
/*	MARK BEGINNING/END OF SIDE GROUPS & SPECIAL VERTEX POINTS */
/*	THESE SHOULD NOT BE DELETED */
   for(i=0;i<nsbd;++i) {
      v0 = svrtx[sbdry[i].el[0]][0];
      v1 = svrtx[sbdry[i].el[sbdry[i].num-1]][1];
      if (v0 != v1) {
         vinfo[v0] = 1;
         vinfo[v1] = 1;
      }
   }
   
   for(i=0;i<nvbd;++i)
      for(j=0;j<vbdry[i].num;++j)
         vinfo[vbdry[i].el[j]] = 1;
 
   
/*	CLASSIFY SIDES AS ACCEPTED OR UNACCEPTED */
   nslst = 0;
   for(i=0;i<nside;++i) {
      if (fltwk[i] > tolsize) {
         tinfo[stri[i][0]] += 1;
         tinfo[MAX(stri[i][1],-1)] += 1;
         putinlst(i);
      }
   }
   
   tinfo[-1] = 0;
   
/*	BEGIN COARSENING ALGORITHM */
   while (nslst > 0) {
      sind = -1;
      for(i=nslst-1;i>=0;--i) {  // START WITH LARGEST SIDE TO DENSITY RATIO
         if (tinfo[stri[intwk2[i]][0]] +tinfo[MAX(stri[intwk2[i]][1],-1)] == 6) continue;
         sind = intwk2[i];
         break;
      }
      assert(sind != -1);
      
/*		COLLAPSE EDGE */
      nfail = collapse(sind);
      printf("collapsing side %d: %d\n",sind,nfail);

      for(i=0;i<nsdel;++i) 
         if (intwk3[sdel[i]] > -1) tkoutlst(sdel[i]);
         
      if (nfail) {
/*			MARK SIDE AS ACCEPTED AND MOVE ON */
         for(i=0;i<2;++i) {
            tind = stri[sind][i];
            if (tind < 0) continue;
            if (tinfo[tind] > 0) --tinfo[tind];
         }
         continue;
      }
         
/*		RECONSTRUCT AFFECTED TINFO ARRAY */
      for(i=0;i<ntdel;++i) {
         tind = tdel[i];
         tinfo[tind] = 0;
         for(j=0;j<3;++j) {
            sind = tside[tind].side[j];
            if (intwk3[sind] > -1) tkoutlst(sind);
            fltwk[sind] = (vlngth[svrtx[sind][0]] +vlngth[svrtx[sind][1]])/
               (2.*distance(svrtx[sind][0],svrtx[sind][1]));
            if (fltwk[sind] > tolsize) {
               putinlst(sind);
               ++tinfo[tind];
            }
         }
      }
   }
      
/*	DELETE LEFTOVER VERTICES */
/* VINFO > NVRTX STORES VRTX MOVEMENT HISTORY */
   for(i=nvrtx-1;i>=0;--i) {
      if (vinfo[i] < 0) {
         vinfo[nvrtx-1] = i;
         dltvrtx(i);
      }
   }
   
/*	FIX BOUNDARY CONDITION POINTERS */
   for(i=0;i<nvbd;++i)
      for(j=0;j<vbdry[i].num;++j) 
         while (vbdry[i].el[j] >= nvrtx) 
            vbdry[i].el[j] = vinfo[vbdry[i].el[j]];  
   
         
/*	DELETE SIDES FROM BOUNDARY CONDITIONS */
   for(i=0;i<nsbd;++i)
      for(j=sbdry[i].num-1;j>=0;--j) 
         if (sinfo[sbdry[i].el[j]] == -3) 
            sbdry[i].el[j] = sbdry[i].el[--sbdry[i].num];
                        
/*	CLEAN UP SIDES */
/* SINFO WILL END UP STORING -1 UNTOUCHED, -2 TOUCHED, or INITIAL INDEX OF UNTOUCHED SIDE */
/*	SINFO > NSIDE WILL STORE MOVEMENT HISTORY */
   for(i=nside-1;i>=0;--i) {
      if (sinfo[i] == -3) {
         if (sinfo[nside-1] >= -1)
            sinfo[i] = MAX(nside-1,sinfo[nside-1]);
         else 
            sinfo[i] = -2;
         sinfo[nside-1] = i;
         fltwk[i] = fltwk[nside-1];  // THIS ALLOWS US TO CALL REBAY WITH OUT RECALCULATING 
         dltside(i);
      }
   }
   
/*	FIX BOUNDARY CONDITION POINTERS */
   for(i=0;i<nsbd;++i)
      for(j=0;j<sbdry[i].num;++j) 
         while (sbdry[i].el[j] >= nside) 
            sbdry[i].el[j] = sinfo[sbdry[i].el[j]]; 
            
   for (i=0;i<nsbd;++i)
      bdrysidereorder(i);
   
/*	CLEAN UP DELETED TRIS */
/* TINFO < NTRI STORES INDEX OF ORIGINAL TRI: TINFO = 0 -> UNMOVED*/
/* TINFO > NTRI STORES TRI MOVEMENT HISTORY */
   for(i=0;i<ntri;++i)
      assert(tinfo[i] <= 0);
      
   for(i=ntri-1;i>=0;--i) {
      if (tinfo[i] < 0) {
         tinfo[i] = MAX(ntri-1,tinfo[ntri-1]);
         tinfo[ntri-1] = i;
         dlttri(i);
      }
   }
   
  return;
}