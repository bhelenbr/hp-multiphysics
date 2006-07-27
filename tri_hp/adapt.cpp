/*
 *  adapt.cpp
 *  planar++
 *
 *  Created by helenbrk on Tue Oct 23 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp.h"
#include "hp_boundary.h"
#include <myblas.h>

 block::ctrl tri_hp::adapt(block::ctrl ctrl_message,FLT tol) {
   block::ctrl state;
   
   if (ctrl_message == block::begin) {
      /* ASSUMES STEPS FOR SETTING UP VLENGTH HAVE BEEN TAKEN ALREADY */
      gbl_ptr->pstr->copy_data(*this);
      treeinit();
      excpt = 0;
   }
   
   if (excpt == 0) {
      state = mesh::adapt(ctrl_message,tol);
      if (state == block::stop) {
         setinfo();
         ++excpt;
      }
      return(state);
   }
   return(block::stop);
}

void tri_hp::updatevdata(int v0) {
   int n,tind,step; 
   FLT r,s;     
      
   gbl_ptr->pstr->findandmvptincurved(vrtx(v0),tind,r,s);
   basis::tri(log2p).ptvalues_rs(r,s);
   
   for(step=0;step<sim::nadapt+1;++step) {
      gbl_ptr->pstr->ugtouht(tind,step);
      basis::tri(log2p).ptprobe(NV,&ugbd(step).v(v0,0),&uht(0)(0),MXTM);
   }
   
   if (gbl_ptr->pstr->td(tind).info > -1) {
      for(step=1;step<sim::nadapt+1;++step) {
         gbl_ptr->pstr->crdtocht(tind,step);
         basis::tri(log2p).ptprobe_bdry(ND,&vrtxbd(step)(v0)(0),&cht(0,0),MXTM);
      }
   }
   else {
      for(step=1;step<sim::nadapt+1;++step) {
         for(n=0;n<ND;++n) 
            vrtxbd(step)(v0)(n) = gbl_ptr->pstr->vrtxbd(step)(gbl_ptr->pstr->td(tind).vrtx(0))(n)*(s +1.)/2.
                                 +gbl_ptr->pstr->vrtxbd(step)(gbl_ptr->pstr->td(tind).vrtx(1))(n)*(-r -s)/2.
                                 +gbl_ptr->pstr->vrtxbd(step)(gbl_ptr->pstr->td(tind).vrtx(2))(n)*(r +1.)/2.;
      }
   }

   return;
}

void tri_hp::updatevdata_bdry(int bnum, int bel, int endpt) {
   int n,sind,sidloc,v0,step;
   FLT psi;
   
   v0 = sd(sbdry(bnum)->el(bel)).vrtx(endpt);
   gbl_ptr->pstr->hp_sbdry(bnum)->findbdrypt(vrtx(v0),sidloc,psi);
   sind = gbl_ptr->pstr->sbdry(bnum)->el(sidloc);

   for(step=0;step<sim::nadapt+1;++step) {
      gbl_ptr->pstr->ugtouht1d(sind,step);
      basis::tri(log2p).ptprobe1d(NV,&ugbd(step).v(v0,0),&uht(0)(0),MXTM);
   }

   if (hp_sbdry(bnum)->is_curved()) {
      for(step=1;step<sim::nadapt+1;++step) {
         gbl_ptr->pstr->crdtocht1d(sind,step);
         basis::tri(log2p).ptprobe1d(ND,&vrtxbd(step)(v0)(0),&cht(0,0),MXTM);
      }
   }
   else {
      for(step=1;step<sim::nadapt+1;++step) {
         for(n=0;n<ND;++n) 
            vrtxbd(step)(v0)(n) = gbl_ptr->pstr->vrtxbd(step)(gbl_ptr->pstr->sd(sind).vrtx(0))(n)*(1. -psi)/2.
                                     +gbl_ptr->pstr->vrtxbd(step)(gbl_ptr->pstr->sd(sind).vrtx(1))(n)*(1. +psi)/2.;
      }
   }
   
   /* FOR INTERNALLY STORED DATA */
   hp_sbdry(bnum)->updatevdata_bdry(bel,endpt,gbl_ptr->pstr->hp_sbdry(bnum));
   
   return;
}

void tri_hp::movevdata(int from, int to) {
   int n,step;
      
   for(step=0;step<sim::nadapt+1;++step) {
      for(n=0;n<NV;++n)
         ugbd(step).v(to,n) = ugbd(step).v(from,n);
         
      for(n=0;n<ND;++n)
         vrtxbd(step)(to)(n) = vrtxbd(step)(from)(n);
   }
   
   return;
}

void tri_hp::movevdata_bdry(int bnum,int bel,int endpt) {
   /* This is just for internal data (if any) */
   hp_sbdry(bnum)->movevdata_bdry(bel,endpt,gbl_ptr->pstr->hp_sbdry(bnum));
}

void tri_hp::updatesdata(int sind) {
   int i,m,n,v0,v1,tind,step,info;
   FLT r,s,upt[NV];
   char uplo[] = "U";
   TinyVector<FLT,2> pt;
   
   if (!sm0) return;
   
   v0 = sd(sind).vrtx(0);
   v1 = sd(sind).vrtx(1);
   
   for(n=0;n<ND;++n)
      basis::tri(log2p).proj1d(vrtx(v0)(n),vrtx(v1)(n),&crd(n)(0,0));

   for(step=0;step<sim::nadapt+1;++step)
      for(n=0;n<NV;++n)
         basis::tri(log2p).proj1d(ugbd(step).v(v0,n),ugbd(step).v(v1,n),&bdwk(step,n)(0,0));
      
   for(i=0;i<basis::tri(log2p).gpx;++i) {
      pt(0) = crd(0)(0,i);
      pt(1) = crd(1)(0,i);
      gbl_ptr->pstr->findinteriorpt(pt,tind,r,s);
      basis::tri(log2p).ptvalues_rs(r,s);
         
      for(step=0;step<sim::nadapt+1;++step) {
         gbl_ptr->pstr->ugtouht(tind,step);
         basis::tri(log2p).ptprobe(NV,upt,&uht(0)(0),MXTM);
         for(n=0;n<NV;++n)   
            bdwk(step,n)(0,i) -= upt[n];
      }
   }            

   for(step=0;step<sim::nadapt+1;++step) {
      for(n=0;n<NV;++n)
         basis::tri(log2p).intgrt1d(&lf(n)(0),&bdwk(step,n)(0,0));

      for(n=0;n<NV;++n) {
         PBTRS(uplo,basis::tri(log2p).sm,basis::tri(log2p).sbwth,1,&basis::tri(log2p).sdiag1d(0,0),basis::tri(log2p).sbwth+1,&lf(n)(2),basis::tri(log2p).sm,info);
         for(m=0;m<basis::tri(log2p).sm;++m) 
            ugbd(step).s(sind,m,n) = -lf(n)(2+m);
      }
   }
   return;
}

void tri_hp::updatesdata_bdry(int bnum,int bel) {
   int m,n,sind,v0,v1,step,stgt,info;
   TinyVector<FLT,2> pt;
   FLT psi;
   FLT upt[NV];
   char uplo[] = "U";
   
   if (!sm0) return;
   
   sind = sbdry(bnum)->el(bel);
   v0 = sd(sind).vrtx(0);
   v1 = sd(sind).vrtx(1);

   for(step=0;step<sim::nadapt+1;++step)
      for(n=0;n<NV;++n)
         basis::tri(log2p).proj1d(ugbd(step).v(v0,n),ugbd(step).v(v1,n),&bdwk(step,n)(0,0));
         
   for(step=0;step<sim::nadapt+1;++step)
      for(n=0;n<ND;++n)
         basis::tri(log2p).proj1d(vrtxbd(step)(v0)(n),vrtxbd(step)(v1)(n),&bdwk(step,n)(1,0));
         
   if (hp_sbdry(bnum)->is_curved()) {

      for(m=0;m<basis::tri(log2p).gpx;++m) {
         pt(0) = bdwk(0,0)(1,m);
         pt(1) = bdwk(0,1)(1,m);
         gbl_ptr->pstr->hp_sbdry(bnum)->findbdrypt(pt,stgt,psi);
         stgt = gbl_ptr->pstr->sbdry(bnum)->el(stgt);


         for(step=0;step<sim::nadapt+1;++step) {
            gbl_ptr->pstr->ugtouht1d(stgt,step);
            basis::tri(log2p).ptprobe1d(NV,upt,&uht(0)(0),MXTM);
            for(n=0;n<NV;++n)   
               bdwk(step,n)(0,m) -= upt[n];
        
            gbl_ptr->pstr->crdtocht1d(stgt,step);
            basis::tri(log2p).ptprobe1d(ND,upt,&cht(0,0),MXTM);
            for(n=0;n<ND;++n)   
               bdwk(step,n)(1,m) -= upt[n];
         }                    
      }     

      for(step=0;step<sim::nadapt+1;++step) {
         for(n=0;n<ND;++n) {
            basis::tri(log2p).intgrt1d(&lf(n)(0),&bdwk(step,n)(1,0));
            PBTRS(uplo,basis::tri(log2p).sm,basis::tri(log2p).sbwth,1,&basis::tri(log2p).sdiag1d(0,0),basis::tri(log2p).sbwth+1,&lf(n)(2),basis::tri(log2p).sm,info);
         
            for(m=0;m<basis::tri(log2p).sm;++m)
               hp_sbdry(bnum)->crdsbd(step,bel,m,n) = -lf(n)(m+2);
         }
      }
   }
   else {
      for(m=0;m<basis::tri(log2p).gpx;++m) {
         pt(0) = bdwk(0,0)(1,m);
         pt(1) = bdwk(0,1)(1,m);
         
         /* FIND PSI */            
         gbl_ptr->pstr->hp_sbdry(bnum)->findbdrypt(pt,stgt,psi);
         stgt = gbl_ptr->pstr->sbdry(bnum)->el(stgt);
         
         /* CALCULATE VALUE OF SOLUTION AT POINT */
         for(step=0;step<sim::nadapt+1;++step) {
            gbl_ptr->pstr->ugtouht1d(stgt,step);
            basis::tri(log2p).ptprobe1d(NV,upt,&uht(0)(0),MXTM);
            for(n=0;n<NV;++n)   
               bdwk(step,n)(0,m) -= upt[n];
         }
      }
   }

   for(step=0;step<sim::nadapt+1;++step) {
      for(n=0;n<NV;++n)
         basis::tri(log2p).intgrt1d(&lf(n)(0),&bdwk(step,n)(0,0));

      for(n=0;n<NV;++n) {
         PBTRS(uplo,basis::tri(log2p).sm,basis::tri(log2p).sbwth,1,&basis::tri(log2p).sdiag1d(0,0),basis::tri(log2p).sbwth+1,&lf(n)(2),basis::tri(log2p).sm,info);
         for(m=0;m<basis::tri(log2p).sm;++m) {
            ugbd(step).s(sind,m,n) = -lf(n)(2+m);
         }
      }
   }
   
   /* UPDATE INTERNAL INFORMATION */
   hp_sbdry(bnum)->updatesdata_bdry(bel,gbl_ptr->pstr->hp_sbdry(bnum));
   
   return;
}

void tri_hp::movesdata(int from, int to) {
   int step;
   
   if (!sm0) return;
   
   for(step=0;step<sim::nadapt+1;++step)
      ugbd(step).s(to,Range::all(),Range::all()) = ugbd(step).s(from,Range::all(),Range::all());

   return;
}

void tri_hp::movesdata_bdry(int bnum,int bel) {
   
   hp_sbdry(bnum)->movesdata_bdry(bel,gbl_ptr->pstr->hp_sbdry(bnum));

   return;
}

void tri_hp::updatetdata(int tind) {
   int i,j,n,ttgt,step,info;
   FLT r,s;
   FLT upt[NV];
   char uplo[] = "U";
   TinyVector<FLT,2> pt;
   
   if (!im0) return;  /* TEMPORARY NEED TO FIX THIS IN MESH SO CAN TURN OFF ENTIRE LOOP */
      
   for(step=0;step<sim::nadapt+1;++step) {
      ugtouht_bdry(tind,step);
      for(n=0;n<NV;++n)
         basis::tri(log2p).proj_bdry(&uht(n)(0),&bdwk(step,n)(0,0),MXGP);
   }
   
   crdtocht(tind);
   for(n=0;n<ND;++n)
      basis::tri(log2p).proj_bdry(&cht(n,0),&crd(n)(0,0),MXGP);
   
   for (i=0; i < basis::tri(log2p).gpx; ++i ) {
      for (j=0; j < basis::tri(log2p).gpn; ++j ) {
         pt(0) = crd(0)(i,j);
         pt(1) = crd(1)(i,j);
         gbl_ptr->pstr->findinteriorpt(pt,ttgt,r,s);
         
         for(step=0;step<sim::nadapt+1;++step) {
            gbl_ptr->pstr->ugtouht(ttgt,step);
            basis::tri(log2p).ptprobe(NV,upt,&uht(0)(0),MXTM);
            for(n=0;n<NV;++n)
               bdwk(step,n)(i,j) -= upt[n];
         }
      }
   }
                  
   for(step=0;step<sim::nadapt+1;++step) {
      for(n=0;n<NV;++n) {
         basis::tri(log2p).intgrt(lf(n).data(),&bdwk(step,n)(0,0),MXGP);
         PBTRS(uplo,basis::tri(log2p).im,basis::tri(log2p).ibwth,1,&basis::tri(log2p).idiag(0,0),basis::tri(log2p).ibwth+1,&lf(n)(basis::tri(log2p).bm),basis::tri(log2p).im,info);
         for(i=0;i<basis::tri(log2p).im;++i)
            ugbd(step).i(tind,i,n) = -lf(n)(basis::tri(log2p).bm+i);
      }
   }
      
   return;
}

void tri_hp::movetdata(int from, int to) {
   int step;
   
   if (!im0) return;  /* TEMPORARY NEED TO FIX THIS IN MESH SO CAN TURN OFF ENTIRE LOOP */

   for(step=0;step<sim::nadapt+1;++step) {
      ugbd(step).i(to,Range::all(),Range::all()) = ugbd(step).i(from,Range::all(),Range::all());
   }
            
   return;
}
      

   

