/*
 *  tadvance.cpp
 *  spectral_hp
 *
 *  Created by helenbrk on Wed Dec 05 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "hp_mgrid.h"

int hp_mgrid::extrap=0;
FLT hp_mgrid::bd[TMSCHEME+1];

#ifdef DROP
extern FLT lam;
#endif

void hp_mgrid::unsteady_sources(int mgrid) {
   int i,j,n,tind,step;
   class hp_mgrid *fmesh;
   
#ifdef DROP
   lam = 1.0e99;
#endif

   if (!mgrid) {
   /* ON FINEST LEVEL MOVE BD INFO TO WORK */
      for(step=0;step<TMSCHEME-1;++step) {  
         for(i=0;i<nvrtx;++i)
            for(n=0;n<ND;++n)
               ugwk[step].v[i][n] = gbl->ugbd[step].v[i][n];

         for(i=0;i<nside*b.sm;++i)
            for(n=0;n<ND;++n)
               ugwk[step].s[i][n] = gbl->ugbd[step].s[i][n];            

         for(i=0;i<ntri*b.im;++i)
            for(n=0;n<ND;++n)
               ugwk[step].i[i][n] = gbl->ugbd[step].i[i][n]; 
               
         for(i=0;i<nvrtx;++i)
            for(n=0;n<ND;++n)
               vrtxwk[step][i][n] = gbl->vrtxbd[step][i][n];
               
         for(i=0;i<nsbd;++i)
            if (sbdry[i].type&CURV_MASK)
               for(j=0;j<sbdry[i].num*b.sm;++j) 
                     for(n=0;n<ND;++n)
                        binfowk[step][i][j].curv[n] = gbl->binfobd[step][i][j].curv[n];
      }
                     
                     
      /* CALCULATE MESH VELOCITY SOURCE TERM */
      for(i=0;i<nvrtx;++i) {
         for(n=0;n<ND;++n) {
            dvrtdt[i][n] = bd[1]*vrtx[i][n];
            for(step=0;step<TMSCHEME-1;++step)
               dvrtdt[i][n] += bd[step+2]*gbl->vrtxbd[step][i][n];
         }
      }
      
      for(i=0;i<nsbd;++i) {
         if (sbdry[i].type&CURV_MASK) {
            for(j=0;j<sbdry[i].num*b.sm;++j) {
               for(n=0;n<ND;++n) {
                  gbl->dbinfodt[i][j].curv[n] = bd[1]*binfo[i][j].curv[n];
                  for(step=0;step<TMSCHEME-1;++step)
                     gbl->dbinfodt[i][j].curv[n] += bd[step+2]*gbl->binfobd[step][i][j].curv[n];
               }
            }
         }
      }
      
      /* CALCULATE 1D MESH SOURCE ON FINE MESH ONLY */
      surfvrttoug();
      setksprg1d();
      surfksrc1d();
   }
   else if (p0 == 1) {
      /* MOVE WORK INFO TO COARSER MESHES */
      fmesh = static_cast<class hp_mgrid *>(fmpt);
      
      /* CALCULATE UNSTEADY SOURCE TERM ON COARSE MESHES */
      for(i=0;i<nvrtx;++i) {
         tind = fine[i].tri;

         for(n=0;n<ND;++n)
            vrtx[i][n] = 0.0;
            
         for(n=0;n<ND;++n)
            dvrtdt[i][n] = 0.0;
         
         for(n=0;n<NV;++n)
            ug.v[i][n] = 0.0;
            
         for(j=0;j<3;++j) {
            for(n=0;n<ND;++n)
               vrtx[i][n] += fine[i].wt[j]*fmesh->vrtx[fmesh->tvrtx[tind][j]][n];
               
            for(n=0;n<ND;++n)
               dvrtdt[i][n] += fine[i].wt[j]*fmesh->dvrtdt[fmesh->tvrtx[tind][j]][n];
               
            for(n=0;n<NV;++n)
               ug.v[i][n] += fine[i].wt[j]*fmesh->ug.v[fmesh->tvrtx[tind][j]][n];
         }
      }
      
      for(step=0;step<TMSCHEME-1;++step) {   
         /* CALCULATE UNSTEADY SOURCE TERM ON COARSE MESHES */
         for(i=0;i<nvrtx;++i) {
            tind = fine[i].tri;

            for(n=0;n<ND;++n)
               vrtx_frst[i][n] = 0.0;
            
            for(n=0;n<NV;++n)
               vug_frst[i][n] = 0.0;
               
            for(j=0;j<3;++j) {
               for(n=0;n<ND;++n)
                  vrtx_frst[i][n] += fine[i].wt[j]*vrtxwk[step][fmesh->tvrtx[tind][j]][n];
                  
               for(n=0;n<NV;++n)
                  vug_frst[i][n] += fine[i].wt[j]*ugwk[step].v[fmesh->tvrtx[tind][j]][n];
            }
         }
         
         for(i=0;i<nvrtx;++i) {
            for(n=0;n<ND;++n)
               vrtxwk[step][i][n] = vrtx_frst[i][n];
            
            for(n=0;n<NV;++n)
               ugwk[step].v[i][n] = vug_frst[i][n];
         }
      }

      /* SET SPRING CONSTANTS FOR COARSER MESHES */
      setksprg1d();
   }
   
   
   /* CALCULATE SOURCE TERMS AT GAUSS POINTS */
   for(tind=0;tind<ntri;++tind) {
      if (tinfo[tind] > -1) {
         crdtocht(tind);
         for(n=0;n<ND;++n)
            b.proj_bdry(cht[n], crd[n], dcrd[n][0], dcrd[n][1]);
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
            cjcb[i][j] = bd[1]*gbl->rho*RAD(i,j)*(dcrd[0][0][i][j]*dcrd[1][1][i][j] -dcrd[1][0][i][j]*dcrd[0][1][i][j]);
            for(n=0;n<ND;++n) {
               dugdt[log2p][n][tind][i][j]  = u[n][i][j]*cjcb[i][j];
            }
            dugdt[log2p][ND][tind][i][j] = cjcb[i][j];
         }            
      }
   }

   /* NOW DO ADDITIONAL TERMS FOR HIGHER-ORDER BD */
   for(step=0;step<TMSCHEME-1;++step) {
      for(tind=0;tind<ntri;++tind) {
         if (tinfo[tind] > -1) {
            crdtocht(tind,vrtxwk[step],binfowk[step]);
            for(n=0;n<ND;++n)
               b.proj_bdry(cht[n], crd[n], dcrd[n][0], dcrd[n][1]);
         }
         else {
            for(n=0;n<ND;++n)
               b.proj(vrtxwk[step][tvrtx[tind][0]][n],vrtxwk[step][tvrtx[tind][1]][n],vrtxwk[step][tvrtx[tind][2]][n],crd[n]);
               
            for(i=0;i<b.gpx;++i) {
               for(j=0;j<b.gpn;++j) {
                  for(n=0;n<ND;++n) {
                     dcrd[n][0][i][j] = 0.5*(vrtxwk[step][tvrtx[tind][1]][n] -vrtxwk[step][tvrtx[tind][0]][n]);
                     dcrd[n][1][i][j] = 0.5*(vrtxwk[step][tvrtx[tind][2]][n] -vrtxwk[step][tvrtx[tind][0]][n]);
                  }
               }
            }
         }
         ugtouht(tind,ugwk[step]);
         for(n=0;n<ND;++n)
            b.proj(uht[n],u[n]);
                     
         for(i=0;i<b.gpx;++i) {
            for(j=0;j<b.gpn;++j) {   
               cjcb[i][j] = bd[step+2]*gbl->rho*RAD(i,j)*(dcrd[0][0][i][j]*dcrd[1][1][i][j] -dcrd[1][0][i][j]*dcrd[0][1][i][j]);
               for(n=0;n<ND;++n) {
                  dugdt[log2p][n][tind][i][j]  += u[n][i][j]*cjcb[i][j];
               }
               dugdt[log2p][ND][tind][i][j] += cjcb[i][j];
            }            
         }
      }
   }

   return;
}

void hp_mgrid::shift(int stage) {
   int i,j,n,step;
   FLT temp;
   
   /* SHIFT BACKWARDS DIFFERENCE STORAGE */
   for(i=0;i<nvrtx;++i)
      for(step=TMSCHEME-2;step>=1;--step)
         for(n=0;n<ND;++n)
            gbl->ugbd[step].v[i][n] = gbl->ugbd[step-1].v[i][n];

   for(i=0;i<nside*b.sm;++i)
      for(step=TMSCHEME-2;step>=1;--step)
         for(n=0;n<ND;++n)
            gbl->ugbd[step].s[i][n] = gbl->ugbd[step-1].s[i][n];            

   for(i=0;i<ntri*b.im;++i)
      for(step=TMSCHEME-2;step>=1;--step)
         for(n=0;n<ND;++n)
            gbl->ugbd[step].i[i][n] = gbl->ugbd[step-1].i[i][n];    
   
   /* SHIFT & EXTRAPOLATE N+1 VALUE */
   for(i=0;i<nvrtx;++i) {
      for(n=0;n<NV;++n) {
         temp = ug.v[i][n] -gbl->ugbd[0].v[i][n];
         gbl->ugbd[0].v[i][n] = ug.v[i][n];
         ug.v[i][n] += extrap*temp;
      }
   }

   for(i=0;i<nside*b.sm;++i) {
      for(n=0;n<NV;++n) {
         temp = ug.s[i][n] -gbl->ugbd[0].s[i][n];
         gbl->ugbd[0].s[i][n] = ug.s[i][n];
         ug.s[i][n] += extrap*temp;
      }
   } 

   for(i=0;i<ntri*b.im;++i) {
      for(n=0;n<NV;++n) {
         temp = ug.i[i][n] -gbl->ugbd[0].i[i][n];
         gbl->ugbd[0].i[i][n] = ug.i[i][n];
         ug.i[i][n] += extrap*temp;
      }
   }   
   
   /* SHIFT BD MESH INFORMATION */
   for(i=0;i<nvrtx;++i)
      for(step=TMSCHEME-2;step>=1;--step)
         for(n=0;n<ND;++n)
            gbl->vrtxbd[step][i][n] = gbl->vrtxbd[step-1][i][n];
            
   for(i=0;i<nsbd;++i)
      if (sbdry[i].type&CURV_MASK)
         for(j=0;j<sbdry[i].num*b.sm;++j) 
            for(step=TMSCHEME-2;step>=1;--step)
                  for(n=0;n<ND;++n)
                     gbl->binfobd[step][i][j].curv[n] = gbl->binfobd[step-1][i][j].curv[n];

   /* SHIFT & EXTRAPOLATE N+1 VALUE */                     
   for(i=0;i<nvrtx;++i) {
      for(n=0;n<ND;++n) {
         temp = vrtx[i][n] -gbl->vrtxbd[0][i][n];
         gbl->vrtxbd[0][i][n] = vrtx[i][n];
         vrtx[i][n] += extrap*temp;
      }
   }
            
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type&CURV_MASK) {
         for(j=0;j<sbdry[i].num*b.sm;++j) {
            for(n=0;n<ND;++n) {
               temp = binfo[i][j].curv[n] -gbl->binfobd[0][i][j].curv[n];
               gbl->binfobd[0][i][j].curv[n] = binfo[i][j].curv[n];
               binfo[i][j].curv[n] += extrap*temp;
            }
         }
      }
   }
   
   /* MOVE UNSTEADY ANALYTICALLY DEFINED SURFACE POINTS TO NEW LOCATION */
   /* ALSO ELIMINATES ERROR FOR NEW ADAPTATION POINTS ON ANALYTICALLY DEFINED SURFACE */
   curvinit(EULR_MASK+INFL_MASK);
   
   /* MOVE EXTRAPOLATED VERTEX INFO FOR FREE SURFACES TO UKNOWN VECTOR */
   surfvrttoug();

   return;
}


void hp_mgrid::setbd(int nsteps) {
   int i;
   
   for(i=0;i<TMSCHEME+1;++i)
      hp_mgrid::bd[i] = 0.0;
   
   switch(nsteps) {
      case(1):
         hp_mgrid::bd[0] =  hp_mgrid::dti;
         hp_mgrid::bd[1] = -hp_mgrid::dti;
         break;
      case(2):
         hp_mgrid::bd[0] =  1.5*hp_mgrid::dti;
         hp_mgrid::bd[1] = -2.0*hp_mgrid::dti;
         hp_mgrid::bd[2] =  0.5*hp_mgrid::dti;
         hp_mgrid::extrap = 1;
         break;
      case(3):
         hp_mgrid::bd[0] = 11./6*hp_mgrid::dti;
         hp_mgrid::bd[1] = -3.*hp_mgrid::dti;
         hp_mgrid::bd[2] = 1.5*hp_mgrid::dti;
         hp_mgrid::bd[3] = -1./3.*hp_mgrid::dti;
         hp_mgrid::extrap = 1;
         break;
   }

   
   return;
}
