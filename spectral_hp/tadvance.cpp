/*
 *  tadvance.cpp
 *  spectral_hp
 *
 *  Created by helenbrk on Wed Dec 05 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "hp_mgrid.h"


void hp_mgrid::tadvance() {
   int i,j,n,tind,step;
   FLT temp;
   
/*********************************************************/	
/* CALCULATE TIME DERIVATIVE SOURCE TERM FOR FINEST MESH */
/*********************************************************/ 
/*	FIRST TERM */ 
	for(tind=0;tind<ntri;++tind) {
      if (tinfo[tind] > -1) {
         crdtouht(tind);
         for(n=0;n<ND;++n)
            b.proj_bdry(uht[n], crd[n], dcrd[n][0], dcrd[n][1]);
      }
      else {
         for(n=0;n<ND;++n)
            b.proj(vrtx[tvrtx[tind][0]][n],vrtx[tvrtx[tind][1]][n],vrtx[tvrtx[tind][2]][n],crd[n]);
            
         for(i=0;i<b.gpx;++i) {
            for(j=0;j<b.gpn;++j) {
               for(n=0;n<ND;++n) {
                  dcrd[n][0][i][j] = 0.5*(vrtx[tvrtx[tind][1]][n] -vrtx[tvrtx[tind][0]][n]);
                  dcrd[n][1][i][j] = 0.5*(vrtx[tvrtx[tind][2]][n] -vrtx[tvrtx[tind][0]][n]);
               }
            }
         }
      }
		ugtouht(tind);
      for(n=0;n<ND;++n)
         b.proj(uht[n],u[n]);
                  
      for(i=0;i<b.gpx;++i) {
         for(j=0;j<b.gpn;++j) {	
            cjcb[i][j] = dcrd[0][0][i][j]*dcrd[1][1][i][j] -dcrd[1][0][i][j]*dcrd[0][1][i][j];
            for(n=0;n<ND;++n)
               gbl->dugdt[n][tind][i][j]  = bd[1]*u[n][i][j]*cjcb[i][j];
            gbl->dugdt[ND][tind][i][j] = bd[1]*gbl->rho*cjcb[i][j];
         }				
      }
   }
   
/*	NOW DO ADDITIONAL TERMS FOR HIGHER-ORDER BD */
   for(step=0;step<NSTEP-1;++step) {
      for(tind=0;tind<ntri;++tind) {
         if (tinfo[tind] > -1) {
            crdtouht(tind,gbl->vrtxbd[step],gbl->binfobd[step]);
            for(n=0;n<ND;++n)
               b.proj_bdry(uht[n], crd[n], dcrd[n][0], dcrd[n][1]);
         }
         else {
            for(n=0;n<ND;++n)
               b.proj(gbl->vrtxbd[step][tvrtx[tind][0]][n],gbl->vrtxbd[step][tvrtx[tind][1]][n],gbl->vrtxbd[step][tvrtx[tind][2]][n],crd[n]);
               
            for(i=0;i<b.gpx;++i) {
               for(j=0;j<b.gpn;++j) {
                  for(n=0;n<ND;++n) {
                     dcrd[n][0][i][j] = 0.5*(gbl->vrtxbd[step][tvrtx[tind][1]][n] -gbl->vrtxbd[step][tvrtx[tind][0]][n]);
                     dcrd[n][1][i][j] = 0.5*(gbl->vrtxbd[step][tvrtx[tind][2]][n] -gbl->vrtxbd[step][tvrtx[tind][0]][n]);
                  }
               }
            }
         }
         ugtouht(tind,gbl->ugbd[step]);
         for(n=0;n<ND;++n)
            b.proj(uht[n],u[n]);
                     
         for(i=0;i<b.gpx;++i) {
            for(j=0;j<b.gpn;++j) {	
               cjcb[i][j] = dcrd[0][0][i][j]*dcrd[1][1][i][j] -dcrd[1][0][i][j]*dcrd[0][1][i][j];
               for(n=0;n<ND;++n)
                  gbl->dugdt[n][tind][i][j]  += bd[step+2]*u[0][i][j]*cjcb[i][j];
               gbl->dugdt[ND][tind][i][j] += bd[step+2]*gbl->rho*cjcb[i][j];
            }				
         }
      }
   }
   
/*	SHIFT BACKWARDS DIFFERENCE STORAGE */
   for(i=0;i<nvrtx;++i)
      for(step=NSTEP-2;step>=1;--step)
         for(n=0;n<ND;++n)
            gbl->ugbd[step].v[i][n] = gbl->ugbd[step-1].v[i][n];

   for(i=0;i<nside*b.sm;++i)
      for(step=NSTEP-2;step>=1;--step)
         for(n=0;n<ND;++n)
            gbl->ugbd[step].s[i][n] = gbl->ugbd[step-1].s[i][n];            

   for(i=0;i<ntri*b.im;++i)
      for(step=NSTEP-2;step>=1;--step)
         for(n=0;n<ND;++n)
            gbl->ugbd[step].i[i][n] = gbl->ugbd[step-1].i[i][n];    
   
/*	SHIFT & EXTRAPOLATE N+1 VALUE */
   for(i=0;i<nvrtx;++i) {
      for(n=0;n<NV;++n) {
         temp = ug.v[i][n] -gbl->ugbd[0].v[i][n];
         gbl->ugbd[0].v[i][n] = ug.v[i][n];
//TEMPORARY         ug.v[i][n] += temp;
      }
   }

   for(i=0;i<nside*b.sm;++i) {
      for(n=0;n<NV;++n) {
         temp = ug.s[i][n] -gbl->ugbd[0].s[i][n];
         gbl->ugbd[0].s[i][n] = ug.s[i][n];
//TEMPORARY         ug.s[i][n] += temp;
      }
   } 

   for(i=0;i<ntri*b.im;++i) {
      for(n=0;n<NV;++n) {
         temp = ug.i[i][n] -gbl->ugbd[0].i[i][n];
         gbl->ugbd[0].i[i][n] = ug.i[i][n];
//TEMPORARY         ug.i[i][n] += temp;
      }
   }   

/*	DO SIMILAR THING FOR MESH VELOCITY TERMS */
/*	CALCULATE MESH VELOCITY SOURCE TERM */
   for(i=0;i<nvrtx;++i) {
      for(n=0;n<ND;++n) {
         dvrtdt[i][n] = bd[1]*vrtx[i][n];
         for(step=0;step<NSTEP-1;++step)
            dvrtdt[i][n] += bd[step+2]*gbl->vrtxbd[step][i][n];
      }
   }
   
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type&CURV_MASK) {
         for(j=0;j<sbdry[i].num;++j) {
            for(n=0;n<ND;++n) {
               gbl->dbinfodt[i][j].curv[n] = bd[1]*binfo[i][j].curv[n];
               for(step=0;step<NSTEP-1;++step)
                  gbl->dbinfodt[i][j].curv[n] += bd[step+2]*gbl->binfobd[step][i][j].curv[n];
            }
         }
      }
   }
   
/* SHIFT BD MESH INFORMATION */
   for(i=0;i<nvrtx;++i)
      for(step=NSTEP-2;step>=1;--step)
         for(n=0;n<ND;++n)
            gbl->vrtxbd[step][i][n] = gbl->vrtxbd[step-1][i][n];
            
   for(i=0;i<nsbd;++i)
      if (sbdry[i].type&CURV_MASK)
         for(j=0;j<sbdry[i].num;++j) 
            for(step=NSTEP-2;step>=1;--step)
                  for(n=0;n<ND;++n)
                     gbl->binfobd[step][i][j].curv[n] = gbl->binfobd[step-1][i][j].curv[n];

/*	SHIFT & EXTRAPOLATE N+1 VALUE */                     
   for(i=0;i<nvrtx;++i) {
      for(n=0;n<ND;++n) {
         temp = vrtx[i][n] -gbl->vrtxbd[0][i][n];
         gbl->vrtxbd[0][i][n] = vrtx[i][n];
//         vrtx[i][n] += temp;
      }
   }
            
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type&CURV_MASK) {
         for(j=0;j<sbdry[i].num;++j) {
            for(n=0;n<ND;++n) {
               temp = binfo[i][j].curv[n] -gbl->binfobd[0][i][j].curv[n];
               gbl->binfobd[0][i][j].curv[n] = binfo[i][j].curv[n];
//               binfo[i][j].curv[n] += temp;
            }
         }
      }
   }

/*	UPDATE UNSTEADY INFLOW VARIABLES */
   setinflow();
         
   return;
}