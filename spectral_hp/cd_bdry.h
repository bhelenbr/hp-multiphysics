/*
 *  cd_bdry.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_cd.h"

template<class BASE> class dirichlet_cd : public BASE {
   
   dirichlet_cd(tri_hp_cd &xin, side_bdry &bin) : BASE(xin,bin) {mytype = "dirichlet_cd";}
   dirichlet_cd* create(tri_hp_cd& xin, side_bdry &bin) const {return new dirichlet_cd(xin,bin);}
   virtual void vdirichlet() {
      int sind,v0;
      
      for(int j=0;j<BASE::base.nel;++j) {
         sind = base.el(j);
         v0 = x.sd(sind).vrtx(0);
         x.cd_gbl->res.v(v0,0) = 0.0;
      }
      v0 = x.sd(sind).vrtx(1);
      x.cd_gbl->res.v(v0,0) = 0.0;
   }
   
   virtual void sdirichlet(int mode) {
      int sind,v0;

      for(int j=0;j<base.nel;++j) {
         sind = base.el(j);
         x.cd_gbl->res.s(sind,mode,0) = 0.0;
      }
   }
      
   block::ctrl tadvance(int excpt) {
      hp_sgeneric<tri_hp_cd>::tadvance(excpt);
      
      if (excpt == 2) {
         /* UPDATE BOUNDARY CONDITION VALUES */
         for(j=0;j<base.nel;++j) {
            sind = base.el(j);
            v0 = x.sd(sind).vrtx(0);
            x = x.vrtx(v0)(0);      
            y = x.vrtx(v0)(1);
            x.ug.v(v0,0) = (*(x.hp_gbl->func))(0,x,y);
         }
         v0 = x.sd(sind).vrtx(1);
         x = x.vrtx(v0)(0);      
         y = x.vrtx(v0)(1);
         x.ug.v(v0,0) = (*(x.hp_gbl->func))(0,x,y);
         
         /*******************/   
         /* SET SIDE VALUES */
         /*******************/
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = base.el(j);
            v0 = x.sd(sind).vrtx(0);
            v1 = x.sd(sind).vrtx(1);
            
            if (is_curved()) {
               x.crdtocht1d(sind);
               for(n=0;n<ND;++n)
                  basis::tri(log2p).proj1d(&x.cht(n,0),&x.crd(n)(0,0),&x.dcrd(n,0)(0,0));
            }
            else {
               for(n=0;n<ND;++n) {
                  basis::tri(log2p).proj1d(x.vrtx(v0)(n),x.vrtx(v1)(n),&x.crd(n)(0,0));
                  
                  for(k=0;k<basis::tri(log2p).gpx;++k)
                     x.dcrd(n,0)(0,k) = 0.5*(x.vrtx(v1)(n)-x.vrtx(v0)(n));
               }
            }

            if (basis::tri(log2p).sm) {
               for(n=0;n<NV;++n)
                  basis::tri(log2p).proj1d(x.ug.v(v0,n),x.ug.v(v1,n),&x.res(n)(0,0));
         
               for(k=0;k<basis::tri(log2p).gpx; ++k)
                  for(n=0;n<NV;++n)
                     x.res(n)(0,k) -= (*(x.hp_gbl->func))(n,x.crd(0)(0,k),x.crd(1)(0,k));
                     
               for(n=0;n<NV;++n)
                  basis::tri(log2p).intgrt1d(&x.lf(n)(0),&x.res(n)(0,0));
         
               indx = sind*sm0;
               for(n=0;n<NV;++n) {
                  PBTRS(uplo,basis::tri(log2p).sm,basis::tri(log2p).sbwth,1,&basis::tri(log2p).sdiag1d(0,0),basis::tri(log2p).sbwth+1,&x.lf(n)(2),basis::tri(log2p).sm,info);
                  for(m=0;m<basis::tri(log2p).sm;++m) 
                     ug.s(sind,m,n) = -lf(n)(2+m);
               }
            }
         }
      }
      return(block::advance);
   }
};
