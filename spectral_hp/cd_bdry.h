/*
 *  cd_bdry.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_cd.h"
#include "myblas.h"

template<class BASE> class dirichlet_cd : public BASE {
   tri_hp_cd &x;
   
   public:
      dirichlet_cd(tri_hp_cd &xin, side_bdry &bin) : BASE(xin,bin), x(xin) {BASE::mytype = "dirichlet_cd";}
      dirichlet_cd* create(tri_hp_cd& xin, side_bdry &bin) const {return new dirichlet_cd(xin,bin);}
      virtual void vdirichlet() {
         int sind,v0;
         
         for(int j=0;j<BASE::base.nel;++j) {
            sind = BASE::base.el(j);
            v0 = BASE::x.sd(sind).vrtx(0);
            BASE::x.hp_gbl->res.v(v0,0) = 0.0;
         }
         v0 = BASE::x.sd(sind).vrtx(1);
         BASE::x.hp_gbl->res.v(v0,0) = 0.0;
      }
      
      virtual void sdirichlet(int mode) {
         int sind;

         for(int j=0;j<BASE::base.nel;++j) {
            sind = BASE::base.el(j);
            BASE::x.hp_gbl->res.s(sind,mode,0) = 0.0;
         }
      }
         
      block::ctrl tadvance(int excpt) {
         int j,k,m,n,v0,v1,sind,indx,info;
         FLT xpt,ypt;
         char uplo[] = "U";
         
         BASE::tadvance(excpt);
         
         if (excpt == 2) {
            /* UPDATE BOUNDARY CONDITION VALUES */
            for(j=0;j<BASE::base.nel;++j) {
               sind = BASE::base.el(j);
               v0 = BASE::x.sd(sind).vrtx(0);
               xpt = BASE::x.vrtx(v0)(0);      
               ypt = BASE::x.vrtx(v0)(1);
               BASE::x.ug.v(v0,0) = (*(BASE::x.hp_gbl->initfunc))(0,xpt,ypt);
            }
            v0 = BASE::x.sd(sind).vrtx(1);
            xpt = BASE::x.vrtx(v0)(0);      
            ypt = BASE::x.vrtx(v0)(1);
            BASE::x.ug.v(v0,0) = (*(BASE::x.hp_gbl->initfunc))(0,xpt,ypt);
            
            /*******************/   
            /* SET SIDE VALUES */
            /*******************/
            for(j=0;j<BASE::base.nel;++j) {
               sind = BASE::base.el(j);
               v0 = BASE::x.sd(sind).vrtx(0);
               v1 = BASE::x.sd(sind).vrtx(1);
               
               if (BASE::is_curved()) {
                  BASE::x.crdtocht1d(sind);
                  for(n=0;n<mesh::ND;++n)
                     basis::tri(BASE::x.log2p).proj1d(&BASE::x.cht(n,0),&BASE::x.crd(n)(0,0),&BASE::x.dcrd(n,0)(0,0));
               }
               else {
                  for(n=0;n<mesh::ND;++n) {
                     basis::tri(BASE::x.log2p).proj1d(BASE::x.vrtx(v0)(n),BASE::x.vrtx(v1)(n),&BASE::x.crd(n)(0,0));
                     
                     for(k=0;k<basis::tri(BASE::x.log2p).gpx;++k)
                        BASE::x.dcrd(n,0)(0,k) = 0.5*(BASE::x.vrtx(v1)(n)-BASE::x.vrtx(v0)(n));
                  }
               }

               if (basis::tri(BASE::x.log2p).sm) {
                  for(n=0;n<BASE::x.NV;++n)
                     basis::tri(BASE::x.log2p).proj1d(BASE::x.ug.v(v0,n),BASE::x.ug.v(v1,n),&BASE::x.res(n)(0,0));
            
                  for(k=0;k<basis::tri(BASE::x.log2p).gpx; ++k)
                     for(n=0;n<BASE::x.NV;++n)
                        BASE::x.res(n)(0,k) -= (*(BASE::x.hp_gbl->initfunc))(n,BASE::x.crd(0)(0,k),BASE::x.crd(1)(0,k));
                        
                  for(n=0;n<BASE::x.NV;++n)
                     basis::tri(BASE::x.log2p).intgrt1d(&BASE::x.lf(n)(0),&BASE::x.res(n)(0,0));
            
                  indx = sind*BASE::x.sm0;
                  for(n=0;n<BASE::x.NV;++n) {
                     PBTRS(uplo,basis::tri(BASE::x.log2p).sm,basis::tri(BASE::x.log2p).sbwth,1,&basis::tri(BASE::x.log2p).sdiag1d(0,0),basis::tri(BASE::x.log2p).sbwth+1,&BASE::x.lf(n)(2),basis::tri(BASE::x.log2p).sm,info);
                     for(m=0;m<basis::tri(BASE::x.log2p).sm;++m) 
                        BASE::x.ug.s(sind,m,n) = -BASE::x.lf(n)(2+m);
                  }
               }
            }
         }
         return(block::advance);
      }
};
