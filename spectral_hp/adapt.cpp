/*
 *  adapt.cpp
 *  planar++
 *
 *  Created by helenbrk on Tue Oct 23 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "spectral_hp.h"
#include<myblas.h>
#include"utilities.h"
#include<assert.h>

/* THIS IS USED IN THE MVPTTOBDRY FUNCTION */
extern class spectral_hp *tgt;

void spectral_hp::adapt(class spectral_hp& bgn, FLT tolerance) {
   int i,j,m,n,v0,v1,sind,stgt,info,indx,indx1,indx2,nvrt0,tind,touchd,snum;
   FLT x,y,psi,upt[NV];
   char uplo[] = "U";
   
/*	COPY EVERYTHING TO BEGIN */
   bgn = *this;

/*	SET TARGET POINTER: USED IN EXTERNAL FUNCTION MVPTTOBDRY PROVIDED TO MESH */
   tgt = &bgn;

/*	REDUCE maxsrch (KEEP TRIANGLE SEARCHES LOCAL) */
   mesh::maxsrch = 15;
   
/* BEGIN ADAPTION PROCEDURE */
   swap();

/*	CALCULATE SIDE LENGTH RATIO */
   for(i=0;i<nside;++i)
      fltwk[i] = (vlngth[svrtx[i][0]] +vlngth[svrtx[i][1]])/
         (2.*distance(svrtx[i][0],svrtx[i][1]));

/* COARSEN */ 
   nvrt0 = nvrtx;
   yaber(1.0/tolerance);
      
/*	MOVE KEPT VERTEX VALUES TO NEW POSITIONS */
   for(i=nvrt0-1;i>=nvrtx;--i) {
      if (vinfo[i] < 0) continue;
      v0 = vinfo[i];
      for(n=0;n<NV;++n)
         vug[v0][n] = vug[i][n];
   }
 
/*	REFINE */
   nvrt0 = nvrtx;
   rebay(tolerance);

/* MARK BOUNDARY VERTICES */
   for(i=nvrt0;i<nvrtx;++i)
      vinfo[i] = -1;

   for(i=0;i<nsbd;++i)
      for(j=0;j<sbdry[i].num;++j)
         vinfo[svrtx[sbdry[i].el[j]][0]] = 0;
      
/*	ASSIGN NEW VALUES */
   for(i=nvrt0;i<nvrtx;++i)
      if (vinfo[i] < 0) bgn.ptprobe(vrtx[i][0], vrtx[i][1], vug[i]);

/* ASSIGN NEW BOUNDARY VERTEX VALUES */
   for(i=0;i<nsbd;++i) {
      for(j=0;j<sbdry[i].num;++j) {
         v0 = svrtx[sbdry[i].el[j]][0];
         if (v0 >= nvrt0) {
            bgn.ptprobe1d(sbdry[i].type,vrtx[v0][0],vrtx[v0][1],vug[v0]);
         }
      }
   }
         
      
/* ASSIGN NEW SIDE VALUES */
   indx = 0;
   indx1 = 0;
   for(sind=0;sind<nside;++sind) {
      switch (sinfo[sind]) {
         case(-1): 
/*				UNTOUCHED/UNMOVED */            
            break;

         case(-2):
/*				TOUCHED */
            if (stri[sind][1] < 0) break; // DO BOUNDARY SIDES SEPARATELY
            

            v0 = svrtx[sind][0];
            v1 = svrtx[sind][1];

            for(n=0;n<ND;++n)
               b.proj1d(vrtx[v0][n],vrtx[v1][n],crd[n][0]);
 
            for(n=0;n<NV;++n)
               b.proj1d(vug[v0][n],vug[v1][n],res[n][0]);
      
            for(i=0;i<b.gpx;++i) {
               bgn.ptprobe(crd[0][0][i],crd[1][0][i],upt);
               for(n=0;n<NV;++n)
                  res[n][0][i] -= upt[n];
            }
                  
            for(n=0;n<NV;++n)
               b.intgrt1d(res[n][0],lf[n]);
         
            for(n=0;n<NV;++n) {
               PBTRS(uplo,b.sm,b.sbwth,1,b.sdiag1d[0],b.sbwth+1,&lf[n][2],b.sm,info);
               for(m=0;m<b.sm;++m) 
                  sug[indx+m][n] = -lf[n][2+m];
            }
            
            break;

         default:
/*				SIDE INTACT BUT MOVED FROM IT'S ORIGINAL LOCATION */
            assert(sinfo[sind] > -1);
            
            indx2 = sinfo[sind]*bgn.sm0;
            for(m=0;m<b.sm;++m)
               for(n=0;n<NV;++n)
                  sug[indx+m][n] = bgn.sug[indx2+m][n];
                  
            break;
      }     
      indx1 += bgn.sm0;
      indx += sm0;
   }

   
/*	UPDATE BOUNDARY SIDES */
   for(i=0;i<nsbd;++i) {
      indx = 0;
      for(j=0;j<sbdry[i].num;++j) {
         sind = sbdry[i].el[j];
         
         switch(sinfo[sind]) {
            case(-1): // UNTOUCHED
               indx1 = ((-bgn.stri[sind][1])%maxsbel)*bgn.sm0;
               for(m=0;m<b.sm;++m)
                  binfo[i][indx+m] = bgn.binfo[i][indx1+m];
               break;
               
            case(-2): // TOUCHED
               v0 = svrtx[sind][0];
               v1 = svrtx[sind][1];
               
               for(n=0;n<ND;++n)
                  b.proj1d(vrtx[v0][n],vrtx[v1][n],crd[n][0]);
                  
               for(n=0;n<NV;++n)
                  b.proj1d(vug[v0][n],vug[v1][n],res[n][0]);
            
               for(m=0;m<b.gpx;++m) {
                  x = crd[0][0][m];
                  y = crd[1][0][m];

/*						MOVE PT TO BOUNDRY */
                  psi = bgn.bdry_locate(sbdry[i].type,x,y,stgt);
                  crd[0][0][m] -= x;
                  crd[1][0][m] -= y;

/*						CALCULATE VALUE OF SOLUTION AT POINT */
                  bgn.ugtouht1d(stgt);  
                  b.ptprobe1d(NV,uht,upt,psi,lf[0]);
                 
                  for(n=0;n<NV;++n)
                  	res[n][0][m] -= upt[n]; 
               }
               
               if (sbdry[i].type&CURV_MASK) {
                  for(n=0;n<ND;++n) {
                     b.intgrt1d(crd[n][0],lf[n]);
                     PBTRS(uplo,b.sm,b.sbwth,1,b.sdiag1d[0],b.sbwth+1,&lf[n][2],b.sm,info);
                  
                     for(m=0;m<b.sm;++m)
                        binfo[i][indx+m].curv[n] = -lf[n][m+2];
                  }
               }
               else {
                  for(n=0;n<ND;++n)
                     for(m=0;m<b.sm;++m)
                        binfo[i][indx+m].curv[n] = 0.0;
               }

               indx1 = sind*sm0;
               for(n=0;n<NV;++n) {
                  b.intgrt1d(res[n][0],lf[n]);
         
                  PBTRS(uplo,b.sm,b.sbwth,1,b.sdiag1d[0],b.sbwth+1,&lf[n][2],b.sm,info);
                  for(m=0;m<b.sm;++m) 
                     sug[indx1+m][n] = -lf[n][2+m];
               }
               break;
               
            default:
/*					SIDE INTACT BUT MOVED FROM IT'S ORIGINAL LOCATION */
               assert(sinfo[sind] > -1);
               indx1 = ((-bgn.stri[sinfo[sind]][1])%maxsbel)*bgn.sm0;
               for(m=0;m<b.sm;++m)
                  binfo[i][indx+m] = bgn.binfo[i][indx1+m];
               break;
         }
         indx += sm0;
      }
   }
   
   if (b.im <= 0) {
      setbcinfo();
      return;
   }

/*	FIGURE OUT WHICH TRIANGLES HAVE BEEN TOUCHED */
   for(tind=0;tind<ntri;++tind) {
      touchd = 0;
      for(snum=0;snum<3;++snum)
         touchd = MIN(touchd,sinfo[tside[tind].side[snum]]);
         
      if (touchd == -2) intwk2[tind] = -2;
      else intwk2[i] = tinfo[tind];
   }
   
/* SET UP BDRY CONDITION INFO */
   setbcinfo();
      
/*	RESET INTERIOR VALUES */
   indx = 0;
   for(tind=0;tind<ntri;++tind) {
      switch (intwk2[tind]) {
         case(0): // UNTOUCHED & UNMOVED
            break;
            
         case(-2): // TOUCHED
            ugtouht_bdry(tind);
            for(n=0;n<NV;++n)
               b.proj_bdry(uht[n],u[n]);
               
            if (tinfo[tind] < 0) {
               for(n=0;n<ND;++n)
                  b.proj(vrtx[tvrtx[tind][0]][n],vrtx[tvrtx[tind][1]][n],vrtx[tvrtx[tind][2]][n],crd[n]);
            }
            else {
               crdtouht(tind);
               for(n=0;n<ND;++n)
                  b.proj_bdry(uht[n],crd[n]);
            }
               
            for (i=0; i < b.gpx; ++i ) {
               for (j=0; j < b.gpn; ++j ) {
                  bgn.ptprobe(crd[0][i][j], crd[1][i][j], upt);
                  for(n=0;n<NV;++n)
                     u[n][i][j] -= upt[n];
               }
            }
                           
            for(n=0;n<NV;++n) {
               b.intgrt(u[n],lf[n]);
               PBTRS(uplo,b.im,b.ibwth,1,b.idiag[0],b.ibwth+1,&lf[n][b.bm],b.im,info);
               for(i=0;i<b.im;++i)
                  iug[indx+i][n] = -lf[n][b.bm+i];
            }
            
            break;

         default:  // MOVED BUT NOT TOUCHED
            
            indx1 = intwk2[tind]*bgn.im0;
            
            for(i=0;i<b.im;++i)
               for(n=0;n<NV;++n)
                  iug[indx+i][n] = bgn.iug[indx1+m][n];
            
            break;
      }
      indx += im0;
   }
   
/*	RESET intwk2 */
   for(i=0;i<ntri;++i)
      intwk2[i] = -1;
      
/*	RESTORE maxsrch in findtri */
   mesh::maxsrch = 3*MAXLST/4;
            
   return;
}
   
   
   
   
   

   

