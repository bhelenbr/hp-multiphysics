#include "bdry_ins.h"
#include <myblas.h>

/*************************************************/
/* SET DIRICHLET BOUNDARY VALUES & FLUXES ********/
/* (THINGS THAT ARE INDEPENDENT OF THE SOLUTION) */
/* BUT NOT INDEPENDENT OF TIME/MESH POSITION    */
/*************************************************/
//void chrctr(FLT rho, FLT gam, double wl[NV], double wr[NV], double norm[ND], double mv[ND]);

#ifdef DROP
extern FLT dydt;
#endif

using namespace bdry_ins;

void neumann::addbflux() {
   int j,k,n,v0,v1,sind;
   TinyVector<FLT,2> pt,mvel,nrm;
   TinyVector<FLT,3> u,flx;

   for(j=0;j<base.nel;++j) {
      sind = base.el(j);
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
         nrm(0) = x.dcrd(1,0)(0,k);
         nrm(1) = -x.dcrd(0,0)(0,k);            
         for(n=0;n<mesh::ND;++n) {
            pt(n) = x.crd(n)(0,k);
            mvel(n) = sim::bd[0]*(x.crd(n)(0,k) -x.crd(n)(1,k));
         }
#ifdef DROP
         mvel(1) += dydt;
#endif
         for(n=0;n<x.NV;++n)
            u(n) = x.u(n)(0,k);
         
         flux(u,pt,mvel,nrm,flx);
                  
         for(n=0;n<x.NV;++n)
            x.res(n)(0,k) = RAD1D(k)*flx(n);

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




block::ctrl inflow::tadvance(int excpt) {
   int j,k,m,n,v0,v1,sind,indx,info;
   TinyVector<FLT,mesh::ND> pt;
   char uplo[] = "U";
   
   hp_side_bdry::tadvance(excpt);
   
   if (excpt == 2) {
      /* UPDATE BOUNDARY CONDITION VALUES */
      for(j=0;j<base.nel;++j) {
         sind = base.el(j);
         v0 = x.sd(sind).vrtx(0);
         for(n=0;n<x.ND;++n)
            x.ug.v(v0,n) = x.hp_gbl->ibc->f(n,x.vrtx(v0));
      }
      v0 = x.sd(sind).vrtx(1);
      for(n=0;n<x.ND;++n)
         x.ug.v(v0,n) = x.hp_gbl->ibc->f(n,x.vrtx(v0));
      
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
                  x.res(n)(0,k) -= x.hp_gbl->ibc->f(n,pt);
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
   return(block::advance);
}

#ifdef SKIP
      /* OUTFLOW BOUNDARY CONDITION    */
      if (sbdry[i].type&OUTF_MASK) {
         indx = 0;
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(0);
            v1 = sd(sind).vrtx(1);
            
            if (sbdry[i].type&CURV_MASK) {
               crdtocht1d(sind);
               for(n=0;n<ND;++n)
                  basis::tri(log2p).proj1d(&cht(n,0),&crd(n)(0,0),&dcrd(n,0)(0,0));
               
               crdtocht1d(sind,dvrtdt,hp_gbl->dbinfodt);
               for(n=0;n<ND;++n)
                  basis::tri(log2p).proj1d(&cht(n,0),&crd(n)(1,0));
            }
            else {
               for(n=0;n<ND;++n) {
                  basis::tri(log2p).proj1d(vrtx(v0)(n),vrtx(v1)(n),&crd(n)(0,0));
                  
                  for(k=0;k<basis::tri(log2p).gpx;++k)
                     dcrd(n,0)(0,k) = 0.5*(vrtx(v1)(n)-vrtx(v0)(n));
               
                  basis::tri(log2p).proj1d(dvrtdt[v0][n],dvrtdt[v1][n],&crd(n)(1,0));
               }
            }
            
            ugtouht1d(sind);
            for(n=0;n<NV;++n)
               basis::tri(log2p).proj1d(&uht(n)(0),&u(n)(0,0));
            
            gam = hp_gbl->rhoi*hp_gbl->tprcn[sd(sind).tri(0)][0][0]/hp_gbl->tprcn[sd(sind).tri(0)][NV-1][NV-1];
            for(k=0;k<basis::tri(log2p).gpx;++k) {
               pt(0) = crd(0)(0,k);
               pt(1) = crd(1)(0,k);

               for(n=0;n<NV;++n) {
                  wl[n] = u(n)(0,k);
                  wr[n] = hp_gbl->ibc->f(n,pt);
               }
               nrm[0] = dcrd(1,0)(0,k);
               nrm[1] = -dcrd(0,0)(0,k);

               for(n=0;n<ND;++n)
                  mvel[n] = sim::bd[0]*crd(n)(0,k) +crd(n)(1,k);
#ifdef DROP
               mvel[1] += dydt;
#endif
                  
               if (!charyes)
                  wl[2] = wr[2];
               else 
                  chrctr(hp_gbl->rho,gam,wl,wr,nrm,mvel);
               
               res(2)(0,k) = hp_gbl->rho*RAD1D(k)*((wl[0] -mvel[0])*nrm[0] +(wl[1] -mvel[1])*nrm[1]);
#ifndef INERTIALESS
               res(0)(0,k) = res(2)(0,k)*wl[0] +wl[2]*RAD1D(k)*nrm[0];
               res(1)(0,k) = res(2)(0,k)*wl[1] +wl[2]*RAD1D(k)*nrm[1];
#else
               res(0)(0,k) = wl[2]*RAD1D(k)*nrm[0];
               res(1)(0,k) = wl[2]*RAD1D(k)*nrm[1]; 
#endif
            }
            
            for(n=0;n<NV;++n)
               basis::tri(log2p).intgrt1d(&lf(n)(0),&res(n)(0,0));
            
            for(n=0;n<NV;++n)
               hp_gbl->res.v(v0,n) += lf(n)(0);

            for(n=0;n<NV;++n)
               hp_gbl->res.v(v1,n) += lf(n)(1);
            
            indx1 = sind*basis::tri(log2p).sm;
            indx = 2;
            for(k=0;k<basis::tri(log2p).sm;++k) {
               for(n=0;n<NV;++n)
                  hp_gbl->res.s(indx1)(n) += lf(n)(indx);
               ++indx1;
               ++indx;
            }
         }
      }

void chrctr(FLT rho, FLT gam, double wl[NV], double wr[NV], double norm[ND], double mv[ND]) {
   FLT ul,vl,ur,vr,pl,pr,cl,cr,rhoi;
   FLT u,um,v,c,den,lam0,lam1,lam2,uvp[3],mag;
   
   rhoi = 1./rho;

   /* CHARACTERISTIC FAR-FIELD B.C. */   
   mag = sqrt(norm[0]*norm[0] + norm[1]*norm[1]);
   
   norm[0] /= mag;
   norm[1] /= mag;
   
   ul =  wl[0]*norm[0] +wl[1]*norm[1];
   vl = -wl[0]*norm[1] +wl[1]*norm[0];
   pl =  wl[2];
      
   /* DEPENDENT ON FREESTREAM CONDITIONS */
   ur =  wr[0]*norm[0] +wr[1]*norm[1];
   vr = -wr[0]*norm[1] +wr[1]*norm[0];
   pr =  wr[2];
      
   um = mv[0]*norm[0] +mv[1]*norm[1];
   
   cl = sqrt((ul-.5*um)*(ul-.5*um) +gam);
   cr = sqrt((ur-.5*um)*(ur-.5*um) +gam);
   c = 0.5*(cl+cr);
   u = 0.5*(ul+ur);
   v = 0.5*(vl+vr);
   
   den = 1./(2*c);
   lam0 = u -um;
   lam1 = u-.5*um +c; /* always positive */
   lam2 = u-.5*um -c; /* always negative */
      
   /* PERFORM CHARACTERISTIC SWAP */
   /* BASED ON LINEARIZATION AROUND UL,VL,PL */
   uvp[0] = ((pl-pr)*rhoi +(ul*lam1 -ur*lam2))*den;
   if (lam0 > 0.0)
      uvp[1] = v*((pr-pl)*rhoi +lam2*(ur-ul))*den/(lam0-lam2) +vl;
   else
      uvp[1] = v*((pr-pl)*rhoi +lam1*(ur-ul))*den/(lam0-lam1) +vr;
   uvp[2] = (rho*(ul -ur)*gam - lam2*pl +lam1*pr)*den;
   
   /* CHANGE BACK TO X,Y COORDINATES */
   wl[0] =  uvp[0]*norm[0] -uvp[1]*norm[1];
   wl[1] =  uvp[0]*norm[1] +uvp[1]*norm[0];
   wl[2] =  uvp[2];

   /* SHOULDN'T CHANGE NORM */   
   norm[0] *= mag;
   norm[1] *= mag;
      
   return;
 
}
#endif
