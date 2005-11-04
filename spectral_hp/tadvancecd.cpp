/*
 *  tadvance.cpp
 *  tri_hp
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
   /* FIRST TERM */ 
   for(tind=0;tind<ntri;++tind) {
      if (td(tind).info > -1) {
         crdtocht(tind);
         for(n=0;n<ND;++n)
            basis::tri(log2p).proj_bdry(&cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), dcrd[n][1],MXGP);
      }
      else {
         for(n=0;n<ND;++n)
            basis::tri(log2p).proj(vrtx(td(tind).vrtx(0))(n),vrtx(td(tind).vrtx(1))(n),vrtx(td(tind).vrtx(2))(n),&crd(n)(0,0),MXGP);
            
         for(i=0;i<basis::tri(log2p).gpx;++i) {
            for(j=0;j<basis::tri(log2p).gpn;++j) {
               for(n=0;n<ND;++n) {
                  dcrd(n,0)(i,j) = 0.5*(vrtx(td(tind).vrtx(1))(n) -vrtx(td(tind).vrtx(0))(n));
                  dcrd(n,1)(i,j) = 0.5*(vrtx(td(tind).vrtx(2))(n) -vrtx(td(tind).vrtx(0))(n));
               }
            }
         }
      }
      ugtouht(tind);
      for(n=0;n<NV;++n)
         basis::tri(log2p).proj(&uht(n)(0),&u(n)(0,0),MXGP);
                  
      for(i=0;i<basis::tri(log2p).gpx;++i) {
         for(j=0;j<basis::tri(log2p).gpn;++j) {   
            cjcb(i,j) = bd[1]*RAD(i,j)*(dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j));
            for(n=0;n<NV;++n)
               hp_gbl->dugdt[n][tind][i][j]  = u(n)(i,j)*cjcb(i,j);
         }            
      }
   }
   
   /* NOW DO ADDITIONAL TERMS FOR HIGHER-ORDER BD */
   for(step=0;step<MXSTEP-1;++step) {
      for(tind=0;tind<ntri;++tind) {
         if (td(tind).info > -1) {
            crdtocht(tind,hp_gbl->vrtxbd[step],hp_gbl->binfobd[step]);
            for(n=0;n<ND;++n)
               basis::tri(log2p).proj_bdry(&cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), dcrd[n][1],MXGP);
         }
         else {
            for(n=0;n<ND;++n)
               basis::tri(log2p).proj(hp_gbl->vrtxbd[step][td(tind).vrtx(0)][n],hp_gbl->vrtxbd[step][td(tind).vrtx(1)][n],hp_gbl->vrtxbd[step][td(tind).vrtx(2)][n],&crd(n)(0,0),MXGP);
               
            for(i=0;i<basis::tri(log2p).gpx;++i) {
               for(j=0;j<basis::tri(log2p).gpn;++j) {
                  for(n=0;n<ND;++n) {
                     dcrd(n,0)(i,j) = 0.5*(hp_gbl->vrtxbd[step][td(tind).vrtx(1)][n] -hp_gbl->vrtxbd[step][td(tind).vrtx(0)][n]);
                     dcrd(n,1)(i,j) = 0.5*(hp_gbl->vrtxbd[step][td(tind).vrtx(2)][n] -hp_gbl->vrtxbd[step][td(tind).vrtx(0)][n]);
                  }
               }
            }
         }
         ugtouht(tind,hp_gbl->ugbd[step]);
         for(n=0;n<NV;++n)
            basis::tri(log2p).proj(&uht(n)(0),&u(n)(0,0),MXGP);
                     
         for(i=0;i<basis::tri(log2p).gpx;++i) {
            for(j=0;j<basis::tri(log2p).gpn;++j) {   
               cjcb(i,j) = bd[step+2]*RAD(i,j)*(dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j));
               for(n=0;n<NV;++n)
                  hp_gbl->dugdt[n][tind][i][j]  += u(n)(i,j)*cjcb(i,j);
            }            
         }
      }
   }
   
   /* SHIFT BACKWARDS DIFFERENCE STORAGE */
   for(i=0;i<nvrtx;++i)
      for(step=MXSTEP-2;step>=1;--step)
         for(n=0;n<NV;++n)
            hp_gbl->ugbd[step].v(i,n) = hp_gbl->ugbd[step-1].v(i,n);

   for(i=0;i<nside*basis::tri(log2p).sm;++i)
      for(step=MXSTEP-2;step>=1;--step)
         for(n=0;n<NV;++n)
            hp_gbl->ugbd[step].s(i)(n) = hp_gbl->ugbd[step-1].s(i)(n);            

   for(i=0;i<ntri*basis::tri(log2p).im;++i)
      for(step=MXSTEP-2;step>=1;--step)
         for(n=0;n<NV;++n)
            hp_gbl->ugbd[step].i(i)(n) = hp_gbl->ugbd[step-1].i(i)(n);    
   
   /* SHIFT & EXTRAPOLATE N+1 VALUE */
   for(i=0;i<nvrtx;++i) {
      for(n=0;n<NV;++n) {
         temp = ug.v(i,n) -hp_gbl->ugbd[0].v(i,n);
         hp_gbl->ugbd[0].v(i,n) = ug.v(i,n);
         ug.v(i,n) += extrap*temp;
      }
   }

   for(i=0;i<nside*basis::tri(log2p).sm;++i) {
      for(n=0;n<NV;++n) {
         temp = ug.s(i)(n) -hp_gbl->ugbd[0].s(i)(n);
         hp_gbl->ugbd[0].s(i)(n) = ug.s(i)(n);
         ug.s(i)(n) += extrap*temp;
      }
   } 

   for(i=0;i<ntri*basis::tri(log2p).im;++i) {
      for(n=0;n<NV;++n) {
         temp = ug.i(i)(n) -hp_gbl->ugbd[0].i(i)(n);
         hp_gbl->ugbd[0].i(i)(n) = ug.i(i)(n);
         ug.i(i)(n) += extrap*temp;
      }
   }   

   /* DO SIMILAR THING FOR MESH VELOCITY TERMS */
   /* CALCULATE MESH VELOCITY SOURCE TERM */
   for(i=0;i<nvrtx;++i) {
      for(n=0;n<ND;++n) {
         dvrtdt[i][n] = bd[1]*vrtx(i)(n);
         for(step=0;step<MXSTEP-1;++step)
            dvrtdt[i][n] += bd[step+2]*hp_gbl->vrtxbd[step][i][n];
      }
   }
   
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type&CURV_MASK) {
         for(j=0;j<sbdry(i)->nel*basis::tri(log2p).sm;++j) {
            for(n=0;n<ND;++n) {
               hp_gbl->dbinfodt[i][j].curv[n] = bd[1]*binfo[i][j].curv[n];
               for(step=0;step<MXSTEP-1;++step)
                  hp_gbl->dbinfodt[i][j].curv[n] += bd[step+2]*hp_gbl->binfobd[step][i][j].curv[n];
            }
         }
      }
   }
   
   /* SHIFT BD MESH INFORMATION */
   for(i=0;i<nvrtx;++i)
      for(step=MXSTEP-2;step>=1;--step)
         for(n=0;n<ND;++n)
            hp_gbl->vrtxbd[step][i][n] = hp_gbl->vrtxbd[step-1][i][n];
            
   for(i=0;i<nsbd;++i)
      if (sbdry[i].type&CURV_MASK)
         for(j=0;j<sbdry(i)->nel*basis::tri(log2p).sm;++j) 
            for(step=MXSTEP-2;step>=1;--step)
                  for(n=0;n<ND;++n)
                     hp_gbl->binfobd[step][i][j].curv[n] = hp_gbl->binfobd[step-1][i][j].curv[n];

   /* SHIFT & EXTRAPOLATE N+1 VALUE */                     
   for(i=0;i<nvrtx;++i) {
      for(n=0;n<ND;++n) {
         temp = vrtx(i)(n) -hp_gbl->vrtxbd[0][i][n];
         hp_gbl->vrtxbd[0][i][n] = vrtx(i)(n);
         vrtx(i)(n) += extrap*temp;
      }
   }
            
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type&CURV_MASK) {
         for(j=0;j<sbdry(i)->nel*basis::tri(log2p).sm;++j) {
            for(n=0;n<ND;++n) {
               temp = binfo[i][j].curv[n] -hp_gbl->binfobd[0][i][j].curv[n];
               hp_gbl->binfobd[0][i][j].curv[n] = binfo[i][j].curv[n];
               binfo[i][j].curv[n] += extrap*temp;
            }
         }
      }
   }
   
   /* MOVE UNSTEADY ANALYTICALLY DEFINED SURFACE POINTS TO NEW LOCATION */
   /* ALSO ELIMINATES ERROR FOR NEW ADAPTATION POINTS ON ANALYTICALLY DEFINED SURFACE */
   curvinit(EULR_MASK+INFL_MASK);
   
   /* UPDATE UNSTEADY INFLOW VARIABLES */
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
            dvrtdt[i][n] += fine[i].wt[j]*fmesh->dvrtdt[fmesh->td(tind).vrtx(j)][n];
      }
   }
   
   return;
}
