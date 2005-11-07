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

 block::ctrl tri_hp::adapt(int excpt,FLT tol) {
   
   /* ASSUMES STEPS FOR SETTING UP VLENGTH HAVE BEEN TAKEN ALREADY */
   if (excpt == 0) {
      hp_gbl->pstr->copy_data(*this);
      treeinit();
   }
   
   return(mesh::adapt(excpt,tol));
   //setbcinfo();
}

void tri_hp::updatevdata(int v0) {
   int n,tind,step; 
   FLT r,s;     
      
   hp_gbl->pstr->findandmvptincurved(vrtx(v0),tind,r,s);
   basis::tri(log2p).ptvalues(r,s);
   
   for(step=0;step<sim::nhist+1;++step) {
      hp_gbl->pstr->ugtouht(tind,step);
      basis::tri(log2p).ptprobe(NV,&ugbd(step).v(v0,0),&uht(0)(0),MXTM);
   }
   
   if (hp_gbl->pstr->td(tind).info > -1) {
      for(step=1;step<sim::nhist+1;++step) {
         hp_gbl->pstr->crdtocht(tind,step);
         basis::tri(log2p).ptprobe_bdry(ND,&vrtxbd(step)(v0)(0),&cht(0,0),MXTM);
      }
   }
   else {
      for(step=1;step<sim::nhist+1;++step) {
         for(n=0;n<ND;++n) 
            vrtxbd(step)(v0)(n) = hp_gbl->pstr->vrtxbd(step)(hp_gbl->pstr->td(tind).vrtx(0))(n)*(s +1.)/2.
                                 +hp_gbl->pstr->vrtxbd(step)(hp_gbl->pstr->td(tind).vrtx(1))(n)*(-r -s)/2.
                                 +hp_gbl->pstr->vrtxbd(step)(hp_gbl->pstr->td(tind).vrtx(2))(n)*(r +1.)/2.;
      }
   }

   return;
}

void tri_hp::updatevdata_bdry(int bnum, int bel, int endpt) {
   int n,sind,v0,step;
   FLT psi;
   
   sind = sbdry(bnum)->el(bel);
   v0 = sd(sind).vrtx(endpt);
   hp_gbl->pstr->hp_sbdry(bnum)->findbdrypt(vrtx(v0),sind,psi);

   for(step=0;step<sim::nhist+1;++step) {
      hp_gbl->pstr->ugtouht1d(sind,step);
      basis::tri(log2p).ptprobe1d(NV,&ugbd(step).v(v0,0),&uht(0)(0),MXTM);
   }

   if (hp_sbdry(bnum)->is_curved()) {
      for(step=1;step<sim::nhist+1;++step) {
         hp_gbl->pstr->crdtocht1d(sind,step);
         basis::tri(log2p).ptprobe1d(ND,&vrtxbd(step)(v0)(0),&cht(0,0),MXTM);
      }
   }
   else {
      for(step=1;step<sim::nhist+1;++step) {
         for(n=0;n<ND;++n) 
            vrtxbd(step)(v0)(n) = hp_gbl->pstr->vrtxbd(step)(hp_gbl->pstr->sd(sind).vrtx(0))(n)*(1. -psi)/2.
                                     +hp_gbl->pstr->vrtxbd(step)(hp_gbl->pstr->sd(sind).vrtx(1))(n)*(1. +psi)/2.;
      }
   }
   
   /* FOR INTERNALLY STORED DATA */
   hp_sbdry(bnum)->updatevdata_bdry(bel,endpt,hp_gbl->pstr->hp_sbdry(bnum));
   
   return;
}

void tri_hp::movevdata(int from, int to) {
   int n,step;
      
   for(step=0;step<sim::nhist+1;++step) {
      for(n=0;n<NV;++n)
         ugbd(step).v(to,n) = hp_gbl->pstr->ugbd(step).v(from,n);
         
      for(n=0;n<ND;++n)
         vrtxbd(step)(to)(n) = hp_gbl->pstr->vrtxbd(step)(from)(n);
   }
   
   return;
}

void tri_hp::movevdata_bdry(int bnum,int bel,int endpt) {
   /* This is just for internal data (if any) */
   hp_sbdry(bnum)->movevdata_bdry(bel,endpt,hp_gbl->pstr->hp_sbdry(bnum));
}

void tri_hp::updatesdata(int sind) {
   int i,m,n,v0,v1,tind,step,info;
   FLT r,s,upt[NV];
   char uplo[] = "U";
   TinyVector<FLT,2> pt;
   
   v0 = sd(sind).vrtx(0);
   v1 = sd(sind).vrtx(1);
   
   for(n=0;n<ND;++n)
      basis::tri(log2p).proj1d(vrtx(v0)(n),vrtx(v1)(n),&crd(n)(0,0));

   for(step=0;step<sim::nhist+1;++step)
      for(n=0;n<NV;++n)
         basis::tri(log2p).proj1d(ugbd(step).v(v0,n),ugbd(step).v(v1,n),&bdwk(step,n)(0,0));
      
   for(i=0;i<basis::tri(log2p).gpx;++i) {
      pt(0) = crd(0)(0,i);
      pt(1) = crd(1)(0,i);
      hp_gbl->pstr->findinteriorpt(pt,tind,r,s);
      basis::tri(log2p).ptvalues(r,s);
         
      for(step=0;step<sim::nhist+1;++step) {
         hp_gbl->pstr->ugtouht(tind,step);
         basis::tri(log2p).ptprobe(NV,upt,&uht(0)(0),MXTM);
         for(n=0;n<NV;++n)   
            bdwk(step,n)(0,i) -= upt[n];
      }
   }            

   for(step=0;step<sim::nhist+1;++step) {
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
   
   sind = sbdry(bnum)->el(bel);
   v0 = sd(sind).vrtx(0);
   v1 = sd(sind).vrtx(1);

   for(n=0;n<ND;++n)
      basis::tri(log2p).proj1d(vrtx(v0)(n),vrtx(v1)(n),&crd(n)(0,0));

   for(step=0;step<sim::nhist+1;++step)
      for(n=0;n<NV;++n)
         basis::tri(log2p).proj1d(ugbd(step).v(v0,n),ugbd(step).v(v1,n),&bdwk(step,n)(0,0));
         
   for(step=0;step<sim::nhist+1;++step)
      for(n=0;n<ND;++n)
         basis::tri(log2p).proj1d(vrtxbd(step)(v0)(n),vrtxbd(step)(v1)(n),&bdwk(step,n)(1,0));
         
   if (hp_sbdry(bnum)->is_curved()) {

      for(m=0;m<basis::tri(log2p).gpx;++m) {
         pt(0) = bdwk(0,0)(1,m);
         pt(1) = bdwk(0,1)(1,m);
         hp_gbl->pstr->hp_sbdry(bnum)->findbdrypt(pt,stgt,psi);

         for(step=0;step<sim::nhist+1;++step) {
            hp_gbl->pstr->ugtouht1d(stgt,step);
            basis::tri(log2p).ptprobe1d(NV,upt,&uht(0)(0),MXTM);
            for(n=0;n<NV;++n)   
               bdwk(step,n)(0,m) -= upt[n];
        
            hp_gbl->pstr->crdtocht1d(stgt,step);
            basis::tri(log2p).ptprobe1d(ND,upt,&cht(0,0),MXTM);
            for(n=0;n<ND;++n)   
               bdwk(step,n)(1,m) -= upt[n];
         }                    
      }     

      for(step=0;step<sim::nhist+1;++step) {
         for(n=0;n<ND;++n) {
            basis::tri(log2p).intgrt1d(&bdwk(step,n)(1,0),&lf(n)(0));
            PBTRS(uplo,basis::tri(log2p).sm,basis::tri(log2p).sbwth,1,&basis::tri(log2p).sdiag1d(0,0),basis::tri(log2p).sbwth+1,&lf(n)(2),basis::tri(log2p).sm,info);
         
            for(m=0;m<basis::tri(log2p).sm;++m)
               hp_sbdry(bnum)->crdsbd(bel,m,n,step) = -lf(n)(m+2);
         }
      }
   }
   else {
      for(m=0;m<basis::tri(log2p).gpx;++m) {
         pt(0) = crd(0)(0,m);
         pt(1) = crd(1)(0,m);
         
         /* FIND PSI */            
         hp_gbl->pstr->hp_sbdry(bnum)->findbdrypt(pt,stgt,psi);

         /* CALCULATE VALUE OF SOLUTION AT POINT */
         for(step=0;step<sim::nhist+1;++step) {
            hp_gbl->pstr->ugtouht1d(stgt,step);
            basis::tri(log2p).ptprobe1d(NV,upt,&uht(0)(0),MXTM);
            for(n=0;n<NV;++n)   
               bdwk(step,n)(0,m) -= upt[n];
         }
      }
   }

   for(step=0;step<sim::nhist+1;++step) {
      for(n=0;n<NV;++n)
         basis::tri(log2p).intgrt1d(&lf(n)(0),&bdwk(step,n)(0,0));

      for(n=0;n<NV;++n) {
         PBTRS(uplo,basis::tri(log2p).sm,basis::tri(log2p).sbwth,1,&basis::tri(log2p).sdiag1d(0,0),basis::tri(log2p).sbwth+1,&lf(n)(2),basis::tri(log2p).sm,info);
         for(m=0;m<basis::tri(log2p).sm;++m) 
            ugbd(step).s(sind,m,n) = -lf(n)(2+m);
      }
   }
   
   /* UPDATE INTERNAL INFORMATION */
   hp_sbdry(bnum)->updatesdata_bdry(bel,hp_gbl->pstr->hp_sbdry(bnum));
   
   return;
}

void tri_hp::movesdata(int from, int to) {
   int indx,indx1,step;
   
   indx = sm0*to;
   indx1 = sm0*from;
   
   for(step=0;step<sim::nhist+1;++step)
      ugbd(step).s(to,Range::all(),Range::all()) = hp_gbl->pstr->ugbd(step).s(from,Range::all(),Range::all());
            
   return;
}

void tri_hp::movesdata_bdry(int bnum,int bel) {
   int m,n,sind,tgtel,step;
   
   if (hp_sbdry(bnum)->is_curved()) {
      sind = sbdry(bnum)->el(bel);
      tgtel = hp_gbl->pstr->getbdryel(sind);
      for(step=0;step<sim::nhist;++step) {
         for(m=0;m<sm0;++m) {
            for(n=0;n<ND;++n) {
               hp_sbdry(bnum)->crdsbd(bel,m,n,step) = hp_gbl->pstr->hp_sbdry(bnum)->crdsbd(tgtel,m,n,step);
            }
         }
      }
   }
   
   hp_sbdry(bnum)->movesdata_bdry(bel,hp_gbl->pstr->hp_sbdry(bnum));

   return;
}

void tri_hp::updatetdata(int tind) {
   int i,j,n,ttgt,step,info;
   FLT r,s;
   FLT upt[NV];
   char uplo[] = "U";
   TinyVector<FLT,2> pt;
      
   for(step=0;step<sim::nhist+1;++step) {
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
         hp_gbl->pstr->findinteriorpt(pt,ttgt,r,s);
         
         for(step=0;step<sim::nhist+1;++step) {
            hp_gbl->pstr->ugtouht(ttgt,step);
            basis::tri(log2p).ptprobe(NV,upt,&uht(0)(0),MXTM);
            for(n=0;n<NV;++n)
               bdwk(step,n)(i,j) -= upt[n];
         }
      }
   }
                  
   for(step=0;step<sim::nhist+1;++step) {
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
      
   for(step=0;step<sim::nhist+1;++step) {
      ugbd(step).i(to,Range::all(),Range::all()) = hp_gbl->pstr->ugbd(step).i(from,Range::all(),Range::all());
   }
            
   return;
}
      

   

