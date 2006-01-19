/*
 *  tadvance.cpp
 *  tri_hp
 *
 *  Created by helenbrk on Wed Dec 05 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp.h"
#include "hp_boundary.h"

/* DIRK SCHEMES */
block::ctrl tri_hp::tadvance(bool coarse,int excpt,Array<mesh::transfer,1> &fv_to_ct,Array<mesh::transfer,1> &cv_to_ft, tri_hp *fmesh) {
   int i,j,n,s,tind;
   static int last_r_mesh;
   block::ctrl state;
   
   /* DO STUFF FOR DEFORMABLE MESH FIRST */
   if (excpt == 0) last_r_mesh = 1;
   
   if (excpt < last_r_mesh) {
      if (log2p == log2pmax && sim::substep == 0 && (mmovement == coupled_deformable || mmovement == uncoupled_deformable)) {
         state = r_mesh::tadvance(coarse,excpt,fv_to_ct,cv_to_ft,fmesh);
         if (state != block::stop) {
            last_r_mesh = excpt+2;
            return(state);
         }
      }
      mp_phase=-1;
      return(block::advance);
   } 
   
   excpt -= last_r_mesh;

   /* MAKE SURE SOLUTION IS CONTINUOUS */
   switch(excpt) {
      case(0): {
         ++mp_phase;
         switch(mp_phase%3) {
            case(0):
               vc0load(mp_phase/3,ug.v.data());
               return(block::stay);
            case(1):
               vmsgpass(mp_phase/3);
               return(block::stay);
            case(2):
               return(static_cast<block::ctrl>(vc0wait_rcv(mp_phase/3,ug.v.data())));
         }
      }
      case(1): {
         mp_phase = -1;
         return(block::advance);
      }
      case(2): {
         if (!sm0) return(block::advance);
         
         ++mp_phase;
         switch(mp_phase%3) {
            case(0):
               sc0load(mp_phase/3,ug.s.data(),0,sm0-1,ug.s.extent(secondDim));
               return(block::stay);
            case(1):
               smsgpass(mp_phase/3);
               return(block::stay);
            case(2):
               return(static_cast<block::ctrl>(sc0wait_rcv(mp_phase/3,ug.s.data(),0,sm0-1,ug.s.extent(secondDim))));
         }
      }
   }

   int stage = sim::substep +sim::esdirk;

   if (!coarse) {
      if (stage > 0) {

         /* BACK CALCULATE K TERM */
         ugbd(stage+1).v(Range(0,nvrtx-1),Range::all()) = (ug.v(Range(0,nvrtx-1),Range::all()) -ugbd(1).v(Range(0,nvrtx-1),Range::all()))*sim::adirk[stage-1][stage-1];
         if (basis::tri(log2p).sm) {
            ugbd(stage+1).s(Range(0,nside-1),Range::all(),Range::all()) = (ug.s(Range(0,nside-1),Range::all(),Range::all()) -ugbd(1).s(Range(0,nside-1),Range::all(),Range::all()))*sim::adirk[stage-1][stage-1];
            if (basis::tri(log2p).im) {
               ugbd(stage+1).i(Range(0,ntri-1),Range::all(),Range::all()) = (ug.i(Range(0,ntri-1),Range::all(),Range::all()) -ugbd(1).i(Range(0,ntri-1),Range::all(),Range::all()))*sim::adirk[stage-1][stage-1];
            }
         }
         for(i=0;i<nvrtx;++i)
            for(n=0;n<ND;++n)
               vrtxbd(stage+1)(i)(n) = (vrtx(i)(n)-vrtxbd(1)(i)(n))*sim::adirk[stage-1][stage-1];
      }
      
      if (sim::substep == 0) {
         /* STORE TILDE W */
         ugbd(1).v(Range(0,nvrtx-1),Range::all()) = ug.v(Range(0,nvrtx-1),Range::all());
         if (basis::tri(log2p).sm) {
            ugbd(1).s(Range(0,nside-1),Range::all(),Range::all()) = ug.s(Range(0,nside-1),Range::all(),Range::all());
            if (basis::tri(log2p).im) {
               ugbd(1).i(Range(0,ntri-1),Range::all(),Range::all()) = ug.i(Range(0,ntri-1),Range::all(),Range::all());
            }
         }

         /* SAME FOR MESH INFORMATION */
         for(i=0;i<nvrtx;++i)
            for(n=0;n<ND;++n)
               vrtxbd(1)(i)(n) = vrtx(i)(n);
      }
         
      /* UPDATE TILDE W */
      for (s=0;s<stage;++s) {         
         ugbd(1).v(Range(0,nvrtx-1),Range::all()) += sim::adirk[stage][s]*ugbd(s+2).v(Range(0,nvrtx-1),Range::all());
         if (basis::tri(log2p).sm) {
            ugbd(1).s(Range(0,nside-1),Range::all(),Range::all()) += sim::adirk[stage][s]*ugbd(s+2).s(Range(0,nside-1),Range::all(),Range::all());
            if (basis::tri(log2p).im) {
               ugbd(1).i(Range(0,ntri-1),Range::all(),Range::all()) += sim::adirk[stage][s]*ugbd(s+2).i(Range(0,ntri-1),Range::all(),Range::all());
            }
         }
         for(i=0;i<nvrtx;++i)
            for(n=0;n<ND;++n)
               vrtxbd(1)(i)(n) += sim::adirk[stage][s]*vrtxbd(s+2)(i)(n);
      }
      
      for(i=0;i<nsbd;++i)
         hp_sbdry(i)->tadvance(0);
   
   }
   else if (p0 == 1) {
      
      /* CALCULATE UNSTEADY SOURCE TERMS ON COARSE MESHES */
      for(i=0;i<nvrtx;++i) {
         tind = cv_to_ft(i).tri;

         ugbd(1).v(i,Range::all()) = 0.0;

         for(n=0;n<ND;++n)
            vrtxbd(1)(i)(n) = 0.0;
         
            
         for(j=0;j<3;++j) {
            ugbd(1).v(i,Range::all()) += cv_to_ft(i).wt(j)*fmesh->ugbd(1).v(fmesh->td(tind).vrtx(j),Range::all());
            for(n=0;n<ND;++n)
               vrtxbd(1)(i)(n) += cv_to_ft(i).wt(j)*fmesh->vrtxbd(1)(fmesh->td(tind).vrtx(j))(n);
         }
      }
   }
      
   calculate_unsteady_sources(coarse);
   
   /* EXTRAPOLATE HERE?? */
   for(i=0;i<nsbd;++i)
      hp_sbdry(i)->tadvance(2);
   
   return(block::stop);
}

/* A GENERIC CALCULATION OF SOURCES FOR AN AUTONOMOUS SYSTEM IN STANDARD FORM */
/* WILL NEED TO BE OVERRIDDEN FOR SPECIAL CASES */
void tri_hp::calculate_unsteady_sources(bool coarse) {
   int i,j,n,tind;
   
   for (int log2p=0;log2p<=log2pmax;++log2p) {
      for(tind=0;tind<ntri;++tind) {
         if (td(tind).info > -1) {
            crdtocht(tind,1);
            for(n=0;n<ND;++n)
               basis::tri(log2p).proj_bdry(&cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), &dcrd(n,1)(0,0),MXGP);
         }
         else {
            for(n=0;n<ND;++n)
               basis::tri(log2p).proj(vrtxbd(1)(td(tind).vrtx(0))(n),vrtxbd(1)(td(tind).vrtx(1))(n),vrtxbd(1)(td(tind).vrtx(2))(n),&crd(n)(0,0),MXGP);

            for(i=0;i<basis::tri(log2p).gpx;++i) {
               for(j=0;j<basis::tri(log2p).gpn;++j) {
                  for(n=0;n<ND;++n) {
                     dcrd(n,0)(i,j) = 0.5*(vrtxbd(1)(td(tind).vrtx(1))(n) -vrtxbd(1)(td(tind).vrtx(0))(n));
                     dcrd(n,1)(i,j) = 0.5*(vrtxbd(1)(td(tind).vrtx(2))(n) -vrtxbd(1)(td(tind).vrtx(0))(n));
                  }
               }
            }
         }
         
         
         ugtouht(tind,1);
         for(n=0;n<NV;++n)
            basis::tri(log2p).proj(&uht(n)(0),&u(n)(0,0),MXGP);
                     
         for(i=0;i<basis::tri(log2p).gpx;++i) {
            for(j=0;j<basis::tri(log2p).gpn;++j) {   
               cjcb(i,j) = -sim::bd[0]*RAD(i,j)*(dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j));
               for(n=0;n<NV;++n)
                  dugdt(log2p,tind,n)(i,j) = u(n)(i,j)*cjcb(i,j);
               for(n=0;n<ND;++n)
                  dxdt(log2p,tind,n)(i,j) = -sim::bd[0]*crd(n)(i,j);
            }            
         }
      }
   }
   
   for(i=0;i<nsbd;++i)
      hp_sbdry(i)->tadvance(1);
   
   return;
}

#ifdef BACKDIFF

void hp_mgrid::unsteady_sources(int mgrid) {
   int i,j,n,tind,step;
   class hp_mgrid *fmesh;

   if (!mgrid) {
   /* ON FINEST LEVEL MOVE BD INFO TO WORK */
      for(step=0;step<TMSCHEME-1;++step) { 
         ugwk(step).v(Range(0,nvrtx-1),Range::all()) = ugbd(step).v(Range(0,nvrtx-1),Range::all());
         ugwk(step).s(Range(0,nside-1),Range::all(),Range::all()) = ugbd(step).s(Range(0,nside-1),Range::all(),Range::all());
         ugwk(step).i(Range(0,ntri-1),Range::all(),Range::all()) = ugbd(step).i(Range(0,ntri-1),Range::all(),Range::all());

         for(i=0;i<nvrtx;++i)
            for(n=0;n<ND;++n)
               vrtxwk[step][i][n] = hp_gbl->vrtxbd[step][i][n];
               
         for(i=0;i<nsbd;++i)
            if (sbdry[i].type&CURV_MASK)
               for(j=0;j<sbdry(i)->nel*basis::tri(log2p).sm;++j) 
                     for(n=0;n<ND;++n)
                        binfowk[step][i][j].curv[n] = hp_gbl->binfobd[step][i][j].curv[n];
      }
                     
                     
      /* CALCULATE MESH VELOCITY SOURCE TERM */
      for(i=0;i<nvrtx;++i) {
         for(n=0;n<ND;++n) {
            dvrtdt[i][n] = bd[1]*vrtx(i)(n);
            for(step=0;step<TMSCHEME-1;++step)
               dvrtdt[i][n] += bd[step+2]*hp_gbl->vrtxbd[step][i][n];
         }
      }
      
      for(i=0;i<nsbd;++i) {
         if (sbdry[i].type&CURV_MASK) {
            for(j=0;j<sbdry(i)->nel*basis::tri(log2p).sm;++j) {
               for(n=0;n<ND;++n) {
                  hp_gbl->dbinfodt[i][j].curv[n] = bd[1]*binfo[i][j].curv[n];
                  for(step=0;step<TMSCHEME-1;++step)
                     hp_gbl->dbinfodt[i][j].curv[n] += bd[step+2]*hp_gbl->binfobd[step][i][j].curv[n];
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
            vrtx(i)(n) = 0.0;
            
         for(n=0;n<ND;++n)
            dvrtdt[i][n] = 0.0;
         
         for(n=0;n<NV;++n)
            ug.v(i,n) = 0.0;
            
         for(j=0;j<3;++j) {
            for(n=0;n<ND;++n)
               vrtx(i)(n) += fine[i].wt[j]*fmesh->vrtx(fmesh->td(tind).vrtx(j))(n);
               
            for(n=0;n<ND;++n)
               dvrtdt[i][n] += fine[i].wt[j]*fmesh->dvrtdt[fmesh->td(tind).vrtx(j)][n];
               
            for(n=0;n<NV;++n)
               ug.v(i,n) += fine[i].wt[j]*fmesh->ug.v(fmesh->td(tind).vrtx(j))(n);
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
                  vrtx_frst[i][n] += fine[i].wt[j]*vrtxwk[step][fmesh->td(tind).vrtx(j)][n];
                  
               for(n=0;n<NV;++n)
                  vug_frst[i][n] += fine[i].wt[j]*ugwk[step].v(fmesh->td(tind).vrtx(j))(n);
            }
         }
         
         for(i=0;i<nvrtx;++i) {
            for(n=0;n<ND;++n)
               vrtxwk[step][i][n] = vrtx_frst[i][n];
            
            for(n=0;n<NV;++n)
               ugwk[step].v(i,n) = vug_frst[i][n];
         }
      }

      /* SET SPRING CONSTANTS FOR COARSER MESHES */
      setksprg1d();
   }
   
   
   /* CALCULATE SOURCE TERMS AT GAUSS POINTS */
   for(tind=0;tind<ntri;++tind) {
      if (td(tind).info > -1) {
         crdtocht(tind);
         for(n=0;n<ND;++n)
            basis::tri(log2p).proj_bdry(&cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), &dcrd(n,1)(0,0),MXGP);
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
      for(n=0;n<ND;++n)
         basis::tri(log2p).proj(&uht(n)(0),&u(n)(0,0),MXGP);
                  
      for(i=0;i<basis::tri(log2p).gpx;++i) {
         for(j=0;j<basis::tri(log2p).gpn;++j) {   
            cjcb(i,j) = bd[1]*hp_gbl->rho*RAD(i,j)*(dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j));
            for(n=0;n<ND;++n) {
               dugdt[log2p][n][tind][i][j]  = u(n)(i,j)*cjcb(i,j);
            }
            dugdt[log2p][ND][tind][i][j] = cjcb(i,j);
         }            
      }
   }

   /* NOW DO ADDITIONAL TERMS FOR HIGHER-ORDER BD */
   for(step=0;step<TMSCHEME-1;++step) {
      for(tind=0;tind<ntri;++tind) {
         if (td(tind).info > -1) {
            crdtocht(tind,vrtxwk[step],binfowk[step]);
            for(n=0;n<ND;++n)
               basis::tri(log2p).proj_bdry(&cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), &dcrd(n,1)(0,0),MXGP);
         }
         else {
            for(n=0;n<ND;++n)
               basis::tri(log2p).proj(vrtxwk[step][td(tind).vrtx(0)][n],vrtxwk[step][td(tind).vrtx(1)][n],vrtxwk[step][td(tind).vrtx(2)][n],&crd(n)(0,0),MXGP);
               
            for(i=0;i<basis::tri(log2p).gpx;++i) {
               for(j=0;j<basis::tri(log2p).gpn;++j) {
                  for(n=0;n<ND;++n) {
                     dcrd(n,0)(i,j) = 0.5*(vrtxwk[step][td(tind).vrtx(1)][n] -vrtxwk[step][td(tind).vrtx(0)][n]);
                     dcrd(n,1)(i,j) = 0.5*(vrtxwk[step][td(tind).vrtx(2)][n] -vrtxwk[step][td(tind).vrtx(0)][n]);
                  }
               }
            }
         }
         ugtouht(tind,ugwk[step]);
         for(n=0;n<ND;++n)
            basis::tri(log2p).proj(&uht(n)(0),&u(n)(0,0),MXGP);
                     
         for(i=0;i<basis::tri(log2p).gpx;++i) {
            for(j=0;j<basis::tri(log2p).gpn;++j) {   
               cjcb(i,j) = bd[step+2]*hp_gbl->rho*RAD(i,j)*(dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j));
               for(n=0;n<ND;++n) {
                  dugdt[log2p][n][tind][i][j]  += u(n)(i,j)*cjcb(i,j);
               }
               dugdt[log2p][ND][tind][i][j] += cjcb(i,j);
            }            
         }
      }
   }

   return;
}

void hp_mgrid::shift() {
   int i,j,n,step;
   FLT temp;
   
   /* SHIFT BACKWARDS DIFFERENCE STORAGE */
   for(step=TMSCHEME-2;step>=1;--step) {
      ugbd(step).v(Range(0,nvrtx-1),Range::all()) = ugbd(step-1).v(Range(0,nvrtx-1),Range::all());
      ugbd(step).s(Range(0,nside-1),Range::all(),Range::all()) = ugbd(step-1).s(Range(0,nside-1),Range::all(),Range::all());            
      ugbd(step).i(Range(0,ntri-1),Range::all(),Range::all()) = ugbd(step-1).i(Range(0,ntri-1),Range::all(),Range::all());            
   }
   
   /* SHIFT & EXTRAPOLATE N+1 VALUE */
   for(i=0;i<nvrtx;++i) {
      for(n=0;n<NV;++n) {
         temp = ug.v(i,n) -ugbd(0).v(i,n);
         hp_gbl->ugbd[0].v(i,n) = ug.v(i,n);
         ug.v(i,n) += extrap*temp;
      }
   }

   for(i=0;i<nside;++i) {
      for(m=0;m<basis::tri(log2p).sm;++m) {
         for(n=0;n<NV;++n) {
            temp = ug.s(i,m,n) -hp_gbl->ugbd[0].s(i,m,n);
            hp_gbl->ugbd[0].s(i,m,n) = ug.s(i,m,n);
            ug.s(i,m,n) += extrap*temp;
         }
      }
   } 

   for(i=0;i<ntri;++i) {
      for(m=0;m<basis::tri(log2p).im;++m) {
         for(n=0;n<NV;++n) {
            temp = ug.i(i,m,n) -hp_gbl->ugbd[0].i(i,m,n);
            hp_gbl->ugbd[0].i(i,m,n) = ug.i(i,m,n);
            ug.i(i,m,n) += extrap*temp;
         }
      }
   }   
   
   /* SHIFT BD MESH INFORMATION */
   for(i=0;i<nvrtx;++i)
      for(step=TMSCHEME-2;step>=1;--step)
         for(n=0;n<ND;++n)
            hp_gbl->vrtxbd[step][i][n] = hp_gbl->vrtxbd[step-1][i][n];
            
   for(i=0;i<nsbd;++i)
      if (sbdry[i].type&CURV_MASK)
         for(j=0;j<sbdry(i)->nel*basis::tri(log2p).sm;++j) 
            for(step=TMSCHEME-2;step>=1;--step)
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
   
   /* MOVE EXTRAPOLATED VERTEX INFO FOR FREE SURFACES TO UKNOWN VECTOR */
   surfvrttoug();
   
   setinflow();

   return;
}

#endif


