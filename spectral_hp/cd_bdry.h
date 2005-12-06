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
      void vdirichlet() {
         int sind,v0;
                  
         for(int j=0;j<BASE::base.nel;++j) {
            sind = BASE::base.el(j);
            v0 = x.sd(sind).vrtx(0);
            x.hp_gbl->res.v(v0,0) = 0.0;
         }
         v0 = x.sd(sind).vrtx(1);
         x.hp_gbl->res.v(v0,0) = 0.0;
      }
      
      void sdirichlet(int mode) {
         int sind;
                  
         for(int j=0;j<BASE::base.nel;++j) {
            sind = BASE::base.el(j);
            x.hp_gbl->res.s(sind,mode,0) = 0.0;
         }
      }
         
      block::ctrl tadvance(int excpt) {
         int j,k,m,n,v0,v1,sind,indx,info;
         TinyVector<FLT,mesh::ND> pt;
         char uplo[] = "U";
         
         BASE::tadvance(excpt);
         
         if (excpt == 2) {
            /* UPDATE BOUNDARY CONDITION VALUES */
            for(j=0;j<BASE::base.nel;++j) {
               sind = BASE::base.el(j);
               v0 = x.sd(sind).vrtx(0);
               x.ug.v(v0,0) = (*(x.hp_gbl->initfunc))(0,x.vrtx(v0));
            }
            v0 = x.sd(sind).vrtx(1);
            x.ug.v(v0,0) = (*(x.hp_gbl->initfunc))(0,x.vrtx(v0));
            
            /*******************/   
            /* SET SIDE VALUES */
            /*******************/
            for(j=0;j<BASE::base.nel;++j) {
               sind = BASE::base.el(j);
               v0 = x.sd(sind).vrtx(0);
               v1 = x.sd(sind).vrtx(1);
               
               if (BASE::is_curved()) {
                  x.crdtocht1d(sind);
                  for(n=0;n<mesh::ND;++n)
                     basis::tri(x.log2p).proj1d(&x.cht(n,0),&x.crd(n)(0,0),&x.dcrd(n,0)(0,0));
               }
               else {
                  for(n=0;n<mesh::ND;++n) {
                     basis::tri(x.log2p).proj1d(x.vrtx(v0)(n),x.vrtx(v1)(n),&x.crd(n)(0,0));
                     
                     for(k=0;k<basis::tri(x.log2p).gpx;++k)
                        x.dcrd(n,0)(0,k) = 0.5*(x.vrtx(v1)(n)-x.vrtx(v0)(n));
                  }
               }

               if (basis::tri(x.log2p).sm) {
                  for(n=0;n<x.NV;++n)
                     basis::tri(x.log2p).proj1d(x.ug.v(v0,n),x.ug.v(v1,n),&x.res(n)(0,0));
            
                  for(k=0;k<basis::tri(x.log2p).gpx; ++k) {
                     pt(0) = x.crd(0)(0,k);
                     pt(1) = x.crd(1)(0,k);
                     for(n=0;n<x.NV;++n)
                        x.res(n)(0,k) -= (*(x.hp_gbl->initfunc))(n,pt);
                  }
                  for(n=0;n<x.NV;++n)
                     basis::tri(x.log2p).intgrt1d(&x.lf(n)(0),&x.res(n)(0,0));
            
                  indx = sind*x.sm0;
                  for(n=0;n<x.NV;++n) {
                     PBTRS(uplo,basis::tri(x.log2p).sm,basis::tri(x.log2p).sbwth,1,&basis::tri(x.log2p).sdiag1d(0,0),basis::tri(x.log2p).sbwth+1,&x.lf(n)(2),basis::tri(x.log2p).sm,info);
                     for(m=0;m<basis::tri(x.log2p).sm;++m) 
                        x.ug.s(sind,m,n) = -x.lf(n)(2+m);
                  }
               }
            }
         }
         return(block::advance);
      }
};

template<class BASE> class neumann_cd : public BASE {
   protected:
      tri_hp_cd &x;
      virtual FLT flux(FLT u, TinyVector<FLT,mesh::ND> x, TinyVector<FLT,mesh::ND> mv, TinyVector<FLT,mesh::ND> norm) {return(0.0);}
   
   public:
      neumann_cd(tri_hp_cd &xin, side_bdry &bin) : BASE(xin,bin), x(xin) {BASE::mytype = "neumann_cd";}
      neumann_cd* create(tri_hp_cd& xin, side_bdry &bin) const {return new neumann_cd(xin,bin);}
      void addbflux() {
         int j,k,n,v0,v1,sind;
         TinyVector<FLT,2> pt,mvel,nrm;
   
         for(j=0;j<BASE::base.nel;++j) {
            sind = BASE::base.el(j);
            v0 = x.sd(sind).vrtx(0);
            v1 = x.sd(sind).vrtx(1);
            
            x.crdtocht1d(sind);
            for(n=0;n<mesh::ND;++n)
               basis::tri(x.log2p).proj1d(&x.cht(n,0),&x.crd(n)(0,0),&x.dcrd(n,0)(0,0));
            
            x.crdtocht1d(sind,1);
            for(n=0;n<mesh::ND;++n)
               basis::tri(x.log2p).proj1d(&x.cht(n,0),&x.crd(n)(1,0));
            
            x.ugtouht1d(sind);
            for(n=0;n<x.NV;++n)
               basis::tri(x.log2p).proj1d(&x.uht(n)(0),&x.u(n)(0,0));
      
            for(k=0;k<basis::tri(x.log2p).gpx;++k) {
               pt(0) = x.crd(0)(0,k);
               pt(1) = x.crd(1)(0,k);
               nrm(0) = x.dcrd(1,0)(0,k);
               nrm(1) = -x.dcrd(0,0)(0,k);
               for(n=0;n<mesh::ND;++n)
                  mvel(n) = sim::bd[0]*(x.crd(n)(0,k) -x.crd(n)(1,k));
   
               x.res(0)(0,k) = RAD1D(k)*flux(x.u(0)(0,k),pt,mvel,nrm);
            }
      
            for(n=0;n<x.NV;++n)
               basis::tri(x.log2p).intgrt1d(&x.lf(n)(0),&x.res(n)(0,0));
            
            for(n=0;n<x.NV;++n)
               x.hp_gbl->res.v(v0,n) += x.lf(n)(0);

            for(n=0;n<x.NV;++n)
               x.hp_gbl->res.v(v1,n) += x.lf(n)(1);
            
            for(k=0;k<basis::tri(x.log2p).sm;++k) {
               for(n=0;n<x.NV;++n)
                  x.hp_gbl->res.s(sind,k,n) += x.lf(n)(k+2);
            }
         }
         return;
      }
};

template<class BASE> class char_cd : public neumann_cd<BASE> {
   public:
      FLT flux(FLT u, TinyVector<FLT,mesh::ND> pt, TinyVector<FLT,mesh::ND> mv, TinyVector<FLT,mesh::ND> norm) {
         FLT vel;

         vel =  (neumann_cd<BASE>::x.cd_gbl->ax-mv(0))*norm(0) +(neumann_cd<BASE>::x.cd_gbl->ay -mv(1))*norm(1);      


         if (vel > 0.0)
            return(vel*u);

         return((*neumann_cd<BASE>::x.cd_gbl->initfunc)(0, pt)*vel);
      }
      char_cd(tri_hp_cd &xin, side_bdry &bin) : neumann_cd<BASE>(xin,bin) {BASE::mytype = "char_cd";}
      char_cd* create(tri_hp_cd& xin, side_bdry &bin) const {return new char_cd(xin,bin);}
};

