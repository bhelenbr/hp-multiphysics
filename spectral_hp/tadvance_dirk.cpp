/*
 *  tadvance.cpp
 *  spectral_hp
 *
 *  Created by helenbrk on Wed Dec 05 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "hp_mgrid.h"

#include <stdlib.h> 
// TEMPOR

#define GRK3 0.43586652150845899941601945
#define C2RK3 (2.-9*GRK3+6.*GRK3*GRK3)/(3*(1-4*GRK3+2*GRK3*GRK3))
#define B2RK3 -3*(1-4*GRK3+2*GRK3*GRK3)*(1-4*GRK3+2*GRK3*GRK3)/(4*(-1+6*GRK3-9*GRK3*GRK3+3*GRK3*GRK3*GRK3))
FLT hp_mgrid::bd[1];
int hp_mgrid::extrap=0;
/* THIS IS THE STANDARD FORM */
//const FLT hp_mgrid::adirk[TMSCHEME][TMSCHEME] = {{GRK3,0.0,0.0},{C2RK3-GRK3,GRK3,0.0},{1-B2RK3-GRK3,B2RK3,GRK3}} 
const FLT tdirk[TMSCHEME][TMSCHEME] = {{GRK3,0.0,0.0},{C2RK3-GRK3,GRK3,0.0},{1-B2RK3-GRK3,B2RK3,GRK3}};
//const FLT hp_mgrid::bdirk[TMSCHEME] = {1-B2RK3-GRK3,B2RK3,GRK3};
//const FLT hp_mgrid::cdirk[TMSCHEME] = {GRK3,C2RK3,1.0};
/* THIS IS THE INCREMENTAL FORM */
const FLT hp_mgrid::adirk[TMSCHEME][TMSCHEME] = { {GRK3,0.0,0.0}, {C2RK3-GRK3,GRK3,0.0}, {1.-B2RK3-C2RK3,B2RK3,GRK3} };
const FLT hp_mgrid::cdirk[TMSCHEME] = {GRK3,C2RK3-GRK3,1.0-C2RK3};

void hp_mgrid::tadvance(int stage) {
   int i,j,n,s,tind;
   
   printf("Checking it out\n");
   
   FLT utest[4];
   FLT  temptime[4];
   temptime[0] = 0.0;
   utest[0] = 0.0;
   for(i=1;i<4;++i) {
      temptime[i] = temptime[i-1] +cdirk[i-1]/dti;
      utest[i] = 0.5*pow(temptime[i],2);
      printf("%e: %e\n",temptime[i], utest[i]);
   }
   
   FLT k[4],wt[4];
   wt[0] = 0.0;
   k[0] = (utest[1]-wt[0])/adirk[0][0];
   wt[1] = wt[0] +adirk[1][0]*k[0];
   k[1] = (utest[2]-wt[1])/adirk[1][1];
   wt[2] = wt[1] +adirk[2][0]*k[0]+adirk[2][1]*k[1];
   k[2] = (utest[3]-wt[2])/adirk[2][2];
   printf("wt: %e %e %e\n",wt[0],wt[1],wt[2]);
   printf("k: %e %e %e\n",k[0]*dti,k[1]*dti,k[2]*dti);
   
   wt[0] = 0.0;
   k[0] = (utest[1]-wt[0])/tdirk[0][0];
   wt[1] = wt[0] + tdirk[1][0]*k[0];
   k[1] = (utest[2]-wt[1])/tdirk[1][1];
   wt[2] = wt[0] + tdirk[2][0]*k[0] +tdirk[2][1]*k[1];
   k[2] = (utest[3]-wt[2])/tdirk[2][2];
   printf("wt: %e %e %e\n",wt[0],wt[1],wt[2]);
   printf("k: %e %e %e\n",k[0]*dti,k[1]*dti,k[2]*dti);   
   
   k[0] = temptime[1];
   k[1] = temptime[2];
   k[2] = temptime[3];
   
   printf("%e %e\n",0.0+(tdirk[2][0]*k[0]+tdirk[2][1]*k[1]+tdirk[2][2]*k[2])/dti,0.5*temptime[3]*temptime[3]);
   printf("%e\n",k[2]*dti);
   
   
   exit(1);
   
   
   
   
   
   if (stage == 0) {
      /* STORE W0 */
      for(i=0;i<nvrtx;++i)
         for(n=0;n<NV;++n)
            gbl->ugbd[0].v[i][n] = ug.v[i][n];
      
      for(i=0;i<nside*b.sm;++i)
         for(n=0;n<NV;++n)
            gbl->ugbd[0].s[i][n] = ug.s[i][n];
      
      for(i=0;i<ntri*b.im;++i)
         for(n=0;n<NV;++n)
            gbl->ugbd[0].i[i][n] = ug.i[i][n];
      
      /* SAME FOR MESH INFORMATION */
      for(i=0;i<nvrtx;++i)
         for(n=0;n<ND;++n)
            gbl->vrtxbd[0][i][n] = vrtx[i][n];
      
      for(i=0;i<nsbd;++i)
         if (sbdry[i].type&CURV_MASK)
            for(j=0;j<sbdry[i].num*b.sm;++j) 
               for(n=0;n<ND;++n)
                  gbl->binfobd[0][i][j].curv[n] = binfo[i][j].curv[n];
      
   }
   else {
      /* CALCULATE K TERM */
      for(i=0;i<nvrtx;++i)
         for(n=0;n<NV;++n)
            gbl->ugbd[stage].v[i][n] = (ug.v[i][n] -gbl->ugbd[0].v[i][n])/adirk[stage-1][stage-1];
   
      for(i=0;i<nside*b.sm;++i)
         for(n=0;n<NV;++n)
            gbl->ugbd[stage].s[i][n] = (ug.s[i][n] -gbl->ugbd[0].s[i][n])/adirk[stage-1][stage-1];
      
      for(i=0;i<ntri*b.im;++i)
         for(n=0;n<NV;++n)
            gbl->ugbd[stage].i[i][n] = (ug.i[i][n] -gbl->ugbd[0].i[i][n])/adirk[stage-1][stage-1];
      
      for(i=0;i<nvrtx;++i)
         for(n=0;n<ND;++n)
            gbl->vrtxbd[stage][i][n] = (vrtx[i][n]-gbl->vrtxbd[0][i][n])/adirk[stage-1][stage-1];
      
      for(i=0;i<nsbd;++i)
         if (sbdry[i].type&CURV_MASK)
            for(j=0;j<sbdry[i].num*b.sm;++j) 
               for(n=0;n<ND;++n)
                  gbl->binfobd[stage][i][j].curv[n] = (binfo[i][j].curv[n]-gbl->binfobd[0][i][j].curv[n])/adirk[stage-1][stage-1];
      
      /* RECALCULATE TILDE W */
      for (s=0;s<stage;++s) {         
         for(i=0;i<nvrtx;++i)
            for(n=0;n<NV;++n)
               gbl->ugbd[0].v[i][n] += adirk[stage][s]*gbl->ugbd[s+1].v[i][n];
         
         for(i=0;i<nside*b.sm;++i)
            for(n=0;n<NV;++n)
               gbl->ugbd[0].s[i][n] += adirk[stage][s]*gbl->ugbd[s+1].s[i][n];
         
         for(i=0;i<ntri*b.im;++i)
            for(n=0;n<NV;++n)
               gbl->ugbd[0].i[i][n] += adirk[stage][s]*gbl->ugbd[s+1].i[i][n];
         
         for(i=0;i<nvrtx;++i)
            for(n=0;n<ND;++n)
               gbl->vrtxbd[0][i][n] += adirk[stage][s]*gbl->vrtxbd[s+1][i][n];
         
         for(i=0;i<nsbd;++i)
            if (sbdry[i].type&CURV_MASK)
               for(j=0;j<sbdry[i].num*b.sm;++j) 
                  for(n=0;n<ND;++n)
                     gbl->binfobd[0][i][j].curv[n] += adirk[stage][s]*gbl->binfobd[s+1][i][j].curv[n];
      }
   }

   /*********************************************************/   
   /* CALCULATE TIME DERIVATIVE SOURCE TERM FOR FINEST MESH */
   /*********************************************************/    
   for(tind=0;tind<ntri;++tind) {
      if (tinfo[tind] > -1) {
         crdtocht(tind,gbl->vrtxbd[0],gbl->binfobd[0]);
         for(n=0;n<ND;++n)
            b.proj_bdry(cht[n], crd[n], dcrd[n][0], dcrd[n][1]);
      }
      else {
         for(n=0;n<ND;++n)
            b.proj(gbl->vrtxbd[0][tvrtx[tind][0]][n],gbl->vrtxbd[0][tvrtx[tind][1]][n],gbl->vrtxbd[0][tvrtx[tind][2]][n],crd[n]);
            
         for(i=0;i<b.gpx;++i) {
            for(j=0;j<b.gpn;++j) {
               for(n=0;n<ND;++n) {
                  dcrd[n][0][i][j] = 0.5*(gbl->vrtxbd[0][tvrtx[tind][1]][n] -gbl->vrtxbd[0][tvrtx[tind][0]][n]);
                  dcrd[n][1][i][j] = 0.5*(gbl->vrtxbd[0][tvrtx[tind][2]][n] -gbl->vrtxbd[0][tvrtx[tind][0]][n]);
               }
            }
         }
      }
      ugtouht(tind,gbl->ugbd[0]);
      for(n=0;n<ND;++n)
         b.proj(uht[n],u[n]);
                  
      for(i=0;i<b.gpx;++i) {
         for(j=0;j<b.gpn;++j) {   
            cjcb[i][j] = -bd[0]*gbl->rho*RAD(i,j)*(dcrd[0][0][i][j]*dcrd[1][1][i][j] -dcrd[1][0][i][j]*dcrd[0][1][i][j]);
            for(n=0;n<ND;++n)
               gbl->dugdt[n][tind][i][j]  = u[n][i][j]*cjcb[i][j];
            gbl->dugdt[ND][tind][i][j] = cjcb[i][j];
         }            
      }
   }
   
   /* DO SIMILAR THING FOR MESH VELOCITY TERMS */
   /* CALCULATE MESH VELOCITY SOURCE TERM */
   for(i=0;i<nvrtx;++i) {
      for(n=0;n<ND;++n) {
         dvrtdt[i][n] = -bd[0]*gbl->vrtxbd[0][i][n];
      }
   }
   
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type&CURV_MASK) {
         for(j=0;j<sbdry[i].num*b.sm;++j) {
            for(n=0;n<ND;++n) {
               gbl->dbinfodt[i][j].curv[n] = -bd[0]*gbl->binfobd[0][i][j].curv[n];
            }
         }
      }
   }

   /* MOVE UNSTEADY ANALYTICALLY DEFINED SURFACE POINTS TO NEW LOCATION */
   /* ALSO ELIMINATES ERROR FOR NEW ADAPTATION POINTS ON ANALYTICALLY DEFINED SURFACE */
   /* UPDATE UNSTEADY INFLOW VARIABLES */
   curvinit(EULR_MASK+INFL_MASK);
   setinflow();

   /* MOVE EXTRAPOLATED VERTEX INFO FOR FREE SURFACES TO UKNOWN VECTOR */
   surfvrttoug();

   return;
}

void hp_mgrid::getfdvrtdt() {
   int i,j,n,tind;
   class hp_mgrid *fmesh;
   
   fmesh = static_cast<class hp_mgrid *>(fmpt);
   
   /* CALCULATE MESH VELOCITY SOURCE TERM ON COARSE MESHES */
   /* TO CALCULATE VUG ON COARSE MESH */
   for(i=0;i<nvrtx;++i) {
      tind = fine[i].tri;

      for(n=0;n<ND;++n)
         dvrtdt[i][n] = 0.0;
         
      for(j=0;j<3;++j) {
         for(n=0;n<ND;++n)
            dvrtdt[i][n] += fine[i].wt[j]*fmesh->dvrtdt[fmesh->tvrtx[tind][j]][n];
      }
   }
   
   return;
}

void hp_mgrid::setbd(int stage) {

   bd[0] = dti/adirk[stage][stage];
   
   return;
}

