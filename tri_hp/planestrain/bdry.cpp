#include "bdry_ps.h"
#include <myblas.h>

/*************************************************/
/* SET DIRICHLET BOUNDARY VALUES & FLUXES ********/
/* (THINGS THAT ARE INDEPENDENT OF THE SOLUTION) */
/* BUT NOT INDEPENDENT OF TIME/MESH POSITION    */
/*************************************************/
//void chrctr(FLT rho, FLT gam, double wl[NV], double wr[NV], double norm[ND], double mv[ND]);


using namespace bdry_ps;

block::ctrl neumann::rsdl(block::ctrl ctrl_message) {
   int j,k,n,v0,v1,sind;
   TinyVector<FLT,2> pt,mvel,nrm;
   TinyVector<FLT,3> u,flx;

   if (ctrl_message == block::begin) {
      for(j=0;j<base.nel;++j) {
         sind = base.el(j);
         v0 = x.sd(sind).vrtx(0);
         v1 = x.sd(sind).vrtx(1);
         
         x.crdtocht1d(sind);
         for(n=0;n<mesh::ND;++n)
            basis::tri(x.log2p).proj1d(&x.cht(n,0),&x.crd(n)(0,0),&x.dcrd(n,0)(0,0));
         
         x.ugtouht1d(sind);
         for(n=0;n<x.NV;++n)
            basis::tri(x.log2p).proj1d(&x.uht(n)(0),&x.u(n)(0,0));

         for(k=0;k<basis::tri(x.log2p).gpx;++k) {
            nrm(0) = x.dcrd(1,0)(0,k);
            nrm(1) = -x.dcrd(0,0)(0,k);            
            for(n=0;n<mesh::ND;++n) {
               pt(n) = x.crd(n)(0,k);
            }
            for(n=0;n<x.NV;++n)
               u(n) = x.u(n)(0,k);
            
            flux(u,pt,nrm,flx);
                     
            for(n=0;n<x.NV;++n)
               x.res(n)(0,k) = flx(n);

         }
         
         for(n=0;n<x.NV;++n)
            basis::tri(x.log2p).intgrt1d(&x.lf(n)(0),&x.res(n)(0,0));
                  
         for(n=0;n<x.NV;++n)
            x.gbl_ptr->res.v(v0,n) += x.lf(n)(0);

         for(n=0;n<x.NV;++n)
            x.gbl_ptr->res.v(v1,n) += x.lf(n)(1);
         
         for(k=0;k<basis::tri(x.log2p).sm;++k) {
            for(n=0;n<x.NV;++n)
               x.gbl_ptr->res.s(sind,k,n) += x.lf(n)(k+2);
         }
      }
   }
   
   return(block::stop);
}




block::ctrl dirichlet::tadvance(bool coarse, block::ctrl ctrl_message) {
   int j,k,m,n,v0,v1,sind,indx,info;
   TinyVector<FLT,mesh::ND> pt;
   char uplo[] = "U";
   block::ctrl state;
   
   if (ctrl_message == block::begin) excpt1 = 0;
   
   if (excpt1 == 0) {
      if ((state = hp_side_bdry::tadvance(coarse,ctrl_message)) != block::stop) return(state);
      ++excpt1;
      
      /* UPDATE BOUNDARY CONDITION VALUES */
      for(j=0;j<base.nel;++j) {
         sind = base.el(j);
         v0 = x.sd(sind).vrtx(0);
         for(n=0;n<x.ND;++n)
            x.ug.v(v0,n) = x.gbl_ptr->ibc->f(n,x.vrtx(v0));
      }
      v0 = x.sd(sind).vrtx(1);
      for(n=0;n<x.ND;++n)
         x.ug.v(v0,n) = x.gbl_ptr->ibc->f(n,x.vrtx(v0));
      
      /*******************/   
      /* SET SIDE VALUES */
      /*******************/
      for(j=0;j<base.nel;++j) {
         sind = base.el(j);
         v0 = x.sd(sind).vrtx(0);
         v1 = x.sd(sind).vrtx(1);
         
         if (is_curved()) {
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
            for(n=0;n<x.ND;++n)
               basis::tri(x.log2p).proj1d(x.ug.v(v0,n),x.ug.v(v1,n),&x.res(n)(0,0));
      
            for(k=0;k<basis::tri(x.log2p).gpx; ++k) {
               pt(0) = x.crd(0)(0,k);
               pt(1) = x.crd(1)(0,k);
               for(n=0;n<x.ND;++n)
                  x.res(n)(0,k) -= x.gbl_ptr->ibc->f(n,pt);
            }
            for(n=0;n<x.ND;++n)
               basis::tri(x.log2p).intgrt1d(&x.lf(n)(0),&x.res(n)(0,0));
      
            indx = sind*x.sm0;
            for(n=0;n<x.ND;++n) {
               PBTRS(uplo,basis::tri(x.log2p).sm,basis::tri(x.log2p).sbwth,1,&basis::tri(x.log2p).sdiag1d(0,0),basis::tri(x.log2p).sbwth+1,&x.lf(n)(2),basis::tri(x.log2p).sm,info);
               for(m=0;m<basis::tri(x.log2p).sm;++m) 
                  x.ug.s(sind,m,n) = -x.lf(n)(2+m);
            }
         }
      }
   }
   return(block::stop);
}

block::ctrl friction_wall::rsdl(block::ctrl ctrl_message) {
   int j,k,m,n,sd,v0,v1,sind,tind;
   TinyVector<FLT,2> pt,mvel,nrm;
   TinyVector<FLT,3> u,flx;
   TinyVector<FLT,mesh::ND> stress;
   FLT visc[x.ND][x.ND][x.ND][x.ND];

   if (ctrl_message == block::begin) {
      for(j=0;j<base.nel;++j) {
         sind = base.el(j);
         v0 = x.sd(sind).vrtx(0);
         v1 = x.sd(sind).vrtx(1);
         tind = x.sd(sind).tri(0);      

         for(sd=0;sd<3;++sd)
            if (x.td(tind).side(sd) == sind) break;
         
         x.crdtocht(tind);
         for(m=basis::tri(x.log2p).bm;m<basis::tri(x.log2p).tm;++m)
            for(n=0;n<x.ND;++n)
               x.cht(n,m) = 0.0;
               
         for(n=0;n<x.ND;++n)
            basis::tri(x.log2p).proj_side(sd,&x.cht(n,0), &x.crd(n)(0,0), &x.dcrd(n,0)(0,0), &x.dcrd(n,1)(0,0));

         x.ugtouht(tind);

         for(n=0;n<x.NV;++n)
            basis::tri(x.log2p).proj_side(sd,&x.uht(n)(0),&x.u(n)(0,0),&x.du(n,0)(0,0),&x.du(n,1)(0,0));

         for (k=0;k<basis::tri(x.log2p).gpx;++k) {
            x.cjcb(0,k) = x.gbl_ptr->mu/(x.dcrd(0,0)(0,k)*x.dcrd(1,1)(0,k) -x.dcrd(1,0)(0,k)*x.dcrd(0,1)(0,k));
            
            /* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
            /* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
            visc[0][0][0][0] =  x.cjcb(0,k)*(2.*x.dcrd(1,1)(0,k)*x.dcrd(1,1)(0,k) +x.dcrd(0,1)(0,k)*x.dcrd(0,1)(0,k));
            visc[0][0][1][1] =  x.cjcb(0,k)*(2.*x.dcrd(1,0)(0,k)*x.dcrd(1,0)(0,k) +x.dcrd(0,0)(0,k)*x.dcrd(0,0)(0,k));
            visc[0][0][0][1] = -x.cjcb(0,k)*(2.*x.dcrd(1,1)(0,k)*x.dcrd(1,0)(0,k) +x.dcrd(0,1)(0,k)*x.dcrd(0,0)(0,k));
   #define     viscI0II0II1II0I visc[0][0][0][1]

            visc[1][1][0][0] =  x.cjcb(0,k)*(x.dcrd(1,1)(0,k)*x.dcrd(1,1)(0,k) +2.*x.dcrd(0,1)(0,k)*x.dcrd(0,1)(0,k));
            visc[1][1][1][1] =  x.cjcb(0,k)*(x.dcrd(1,0)(0,k)*x.dcrd(1,0)(0,k) +2.*x.dcrd(0,0)(0,k)*x.dcrd(0,0)(0,k));
            visc[1][1][0][1] = -x.cjcb(0,k)*(x.dcrd(1,1)(0,k)*x.dcrd(1,0)(0,k) +2.*x.dcrd(0,1)(0,k)*x.dcrd(0,0)(0,k));
   #define     viscI1II1II1II0I visc[1][1][0][1]
            
            visc[0][1][0][0] = -x.cjcb(0,k)*x.dcrd(0,1)(0,k)*x.dcrd(1,1)(0,k);
            visc[0][1][1][1] = -x.cjcb(0,k)*x.dcrd(0,0)(0,k)*x.dcrd(1,0)(0,k);
            visc[0][1][0][1] =  x.cjcb(0,k)*x.dcrd(0,1)(0,k)*x.dcrd(1,0)(0,k);
            visc[0][1][1][0] =  x.cjcb(0,k)*x.dcrd(0,0)(0,k)*x.dcrd(1,1)(0,k);

            /* OTHER SYMMETRIES    */            
   #define     viscI1II0II0II0I visc[0][1][0][0]
   #define     viscI1II0II1II1I visc[0][1][1][1]
   #define     viscI1II0II0II1I visc[0][1][1][0]
   #define     viscI1II0II1II0I visc[0][1][0][1]

            stress(0) =   (-x.u(2)(0,k)*x.dcrd(1,0)(0,k) 
                        -viscI0II0II1II0I*x.du(0,0)(0,k) -visc[0][1][1][0]*x.du(1,0)(0,k)
                        -visc[0][0][1][1]*x.du(0,1)(0,k) -visc[0][1][1][1]*x.du(1,1)(0,k));															
            stress(1) =   ( x.u(2)(0,k)*x.dcrd(0,0)(0,k)
                        -viscI1II0II1II0I*x.du(0,0)(0,k) -viscI1II1II1II0I*x.du(1,0)(0,k)
                        -viscI1II0II1II1I*x.du(0,1)(0,k) -visc[1][1][1][1]*x.du(1,1)(0,k));
                        
            nrm(0) = x.dcrd(1,0)(0,k);
            nrm(1) = -x.dcrd(0,0)(0,k);            
            for(n=0;n<mesh::ND;++n) {
               pt(n) = x.crd(n)(0,k);
            }

            for(n=0;n<x.NV;++n)
               u(n) = x.u(n)(0,k);
            
            flux(u,stress,pt,nrm,flx);
                     
            for(n=0;n<x.NV;++n)
               x.res(n)(0,k) = flx(n);

         }

         for(n=0;n<x.NV;++n)
            basis::tri(x.log2p).intgrt1d(&x.lf(n)(0),&x.res(n)(0,0));
                  
         for(n=0;n<x.NV;++n)
            x.gbl_ptr->res.v(v0,n) += x.lf(n)(0);

         for(n=0;n<x.NV;++n)
            x.gbl_ptr->res.v(v1,n) += x.lf(n)(1);
         
         for(k=0;k<basis::tri(x.log2p).sm;++k) {
            for(n=0;n<x.NV;++n)
               x.gbl_ptr->res.s(sind,k,n) += x.lf(n)(k+2);
         }
      }
   }
   
   return(block::stop);
}
