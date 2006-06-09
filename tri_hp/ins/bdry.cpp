#include "bdry_ins.h"
#include <myblas.h>

/*************************************************/
/* SET DIRICHLET BOUNDARY VALUES & FLUXES ********/
/* (THINGS THAT ARE INDEPENDENT OF THE SOLUTION) */
/* BUT NOT INDEPENDENT OF TIME/MESH POSITION    */
/*************************************************/
//void chrctr(FLT rho, FLT gam, double wl[NV], double wr[NV], double norm[ND], double mv[ND]);

using namespace bdry_ins;

void generic::output(std::ostream& fout, tri_hp::filetype typ,int tlvl) {
	int i,m,n,ind,sind,tind,sd;
	FLT visc[mesh::ND][mesh::ND][mesh::ND][mesh::ND];
   TinyVector<FLT,2> drag, norm, mvel;
   TinyVector<FLT,3> flux;
   FLT circ,moment,convect;
   
   switch(typ) {
      case(tri_hp::text): {
         hp_side_bdry::output(fout,typ);
         break;
      }
      case(tri_hp::tecplot): {
         if (!report_flag) return;
         
         drag = 0.0;
         flux = 0.0;
         moment = 0.0;
         circ = 0.0;
         
         for(ind=0; ind < base.nel; ++ind) {
            sind = base.el(ind);
            tind = x.sd(sind).tri(0);      
            
            for(sd=0;sd<3;++sd)
               if (x.td(tind).side(sd) == sind) break;
            assert(sd != 3);
            
            x.crdtocht(tind);
            for(m=basis::tri(x.log2p).bm;m<basis::tri(x.log2p).tm;++m)
               for(n=0;n<mesh::ND;++n)
                  x.cht(n,m) = 0.0;
                  
            for(n=0;n<mesh::ND;++n)
               basis::tri(x.log2p).proj_side(sd,&x.cht(n,0), &x.crd(n)(0,0), &x.dcrd(n,0)(0,0), &x.dcrd(n,1)(0,0));

            x.ugtouht(tind);
            for(n=0;n<x.NV;++n)
               basis::tri(x.log2p).proj_side(sd,&x.uht(n)(0),&x.u(n)(0,0),&x.du(n,0)(0,0),&x.du(n,1)(0,0));

            for (i=0;i<basis::tri(x.log2p).gpx;++i) {
               circ += basis::tri(x.log2p).wtx(i)*sqrt(x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i) +x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i));
               x.cjcb(0,i) = x.ins_gbl->mu*RAD(x.crd(0)(0,i))/(x.dcrd(0,0)(0,i)*x.dcrd(1,1)(0,i) -x.dcrd(1,0)(0,i)*x.dcrd(0,1)(0,i));
               
               /* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
               /* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
               visc[0][0][0][0] =  x.cjcb(0,i)*(2.*x.dcrd(1,1)(0,i)*x.dcrd(1,1)(0,i) +x.dcrd(0,1)(0,i)*x.dcrd(0,1)(0,i));
               visc[0][0][1][1] =  x.cjcb(0,i)*(2.*x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i));
               visc[0][0][0][1] = -x.cjcb(0,i)*(2.*x.dcrd(1,1)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,1)(0,i)*x.dcrd(0,0)(0,i));
#define        viscI0II0II1II0I visc[0][0][0][1]

               visc[1][1][0][0] =  x.cjcb(0,i)*(x.dcrd(1,1)(0,i)*x.dcrd(1,1)(0,i) +2.*x.dcrd(0,1)(0,i)*x.dcrd(0,1)(0,i));
               visc[1][1][1][1] =  x.cjcb(0,i)*(x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i) +2.*x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i));
               visc[1][1][0][1] = -x.cjcb(0,i)*(x.dcrd(1,1)(0,i)*x.dcrd(1,0)(0,i) +2.*x.dcrd(0,1)(0,i)*x.dcrd(0,0)(0,i));
#define        viscI1II1II1II0I visc[1][1][0][1]
               
               visc[0][1][0][0] = -x.cjcb(0,i)*x.dcrd(0,1)(0,i)*x.dcrd(1,1)(0,i);
               visc[0][1][1][1] = -x.cjcb(0,i)*x.dcrd(0,0)(0,i)*x.dcrd(1,0)(0,i);
               visc[0][1][0][1] =  x.cjcb(0,i)*x.dcrd(0,1)(0,i)*x.dcrd(1,0)(0,i);
               visc[0][1][1][0] =  x.cjcb(0,i)*x.dcrd(0,0)(0,i)*x.dcrd(1,1)(0,i);

               /* OTHER SYMMETRIES    */            
#define        viscI1II0II0II0I visc[0][1][0][0]
#define        viscI1II0II1II1I visc[0][1][1][1]
#define        viscI1II0II0II1I visc[0][1][1][0]
#define        viscI1II0II1II0I visc[0][1][0][1]



               drag(0) -=   basis::tri(x.log2p).wtx(i)*(-x.u(2)(0,i)*RAD(x.crd(0)(0,i))*x.dcrd(1,0)(0,i) 
                           -viscI0II0II1II0I*x.du(0,0)(0,i) -visc[0][1][1][0]*x.du(1,0)(0,i)
                           -visc[0][0][1][1]*x.du(0,1)(0,i) -visc[0][1][1][1]*x.du(1,1)(0,i));															
               drag(1) -=   basis::tri(x.log2p).wtx(i)*( x.u(2)(0,i)*RAD(x.crd(0)(0,i))*x.dcrd(0,0)(0,i)
                           -viscI1II0II1II0I*x.du(0,0)(0,i) -viscI1II1II1II0I*x.du(1,0)(0,i)
                           -viscI1II0II1II1I*x.du(0,1)(0,i) -visc[1][1][1][1]*x.du(1,1)(0,i));
                           
                           
               norm(0) = x.dcrd(1,0)(0,i);
               norm(1) = -x.dcrd(0,0)(0,i);            
               for(n=0;n<mesh::ND;++n) {
                  mvel(n) = sim::bd[0]*(x.crd(n)(0,i) -dxdt(x.log2p,ind)(n,i));
#ifdef DROP
                  mvel(n) += tri_hp_ins::mesh_ref_vel(n);
#endif
               }

               
               convect = basis::tri(x.log2p).wtx(i)*RAD(x.crd(0)(0,i))*((x.u(0)(0,i)-mvel(0))*norm(0) +(x.u(1)(0,i)-mvel(1))*norm(1));
               flux(2) -= convect;
               flux(0) -= x.u(0)(0,i)*convect;
               flux(1) -= x.u(1)(0,i)*convect;
            }				
         }
         fout << base.idprefix << " circumference: " << circ << std::endl;
         fout << base.idprefix << " drag: " << drag << std::endl;
         fout << base.idprefix << " flux: " << flux << std::endl; 
      }
   }
   
	return;
}




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
               mvel(n) = sim::bd[0]*(x.crd(n)(0,k) -dxdt(x.log2p,j)(n,k));
#ifdef DROP
               mvel(n) += tri_hp_ins::mesh_ref_vel(n);
#endif
            }
            for(n=0;n<x.NV;++n)
               u(n) = x.u(n)(0,k);
            
            flux(u,pt,mvel,nrm,flx);
                     
            for(n=0;n<x.NV;++n)
               x.res(n)(0,k) = RAD(x.crd(0)(0,k))*flx(n);

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
   }
   return(block::stop);
}




block::ctrl inflow::tadvance(bool coarse, block::ctrl ctrl_message) {
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
   return(block::stop);
}

block::ctrl symmetry::tadvance(bool coarse, block::ctrl ctrl_message) {
   int j,m,v0,sind;
   TinyVector<FLT,mesh::ND> pt;
   block::ctrl state;
   
   if (ctrl_message == block::begin) excpt1 = 0;
   
   if (excpt1 == 0) {
      if ((state = hp_side_bdry::tadvance(coarse,ctrl_message)) != block::stop) return(state);
      ++excpt1;
   
      /* UPDATE BOUNDARY CONDITION VALUES */
      for(j=0;j<base.nel;++j) {
         sind = base.el(j);
         v0 = x.sd(sind).vrtx(0);
         x.ug.v(v0,0) = 0.0;
      }
      v0 = x.sd(sind).vrtx(1);
      x.ug.v(v0,0) = 0.0;

      /*******************/   
      /* SET SIDE VALUES */
      /*******************/
      for(j=0;j<base.nel;++j) {
         sind = base.el(j);
         for(m=0;m<basis::tri(x.log2p).sm;++m) {
            x.ug.s(sind,m,0) = 0.0;
         }
      }
   }
   
   return(block::stop);
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
                  
               if (!charyes)
                  wl[2] = wr[2];
               else 
                  chrctr(hp_gbl->rho,gam,wl,wr,nrm,mvel);
               
               res(2)(0,k) = hp_gbl->rho*RAD(x.crd(0)(0,k))*((wl[0] -mvel[0])*nrm[0] +(wl[1] -mvel[1])*nrm[1]);
#ifndef INERTIALESS
               res(0)(0,k) = res(2)(0,k)*wl[0] +wl[2]*RAD(x.crd(0)(0,k))*nrm[0];
               res(1)(0,k) = res(2)(0,k)*wl[1] +wl[2]*RAD(x.crd(0)(0,k))*nrm[1];
#else
               res(0)(0,k) = wl[2]*RAD(x.crd(0)(0,k))*nrm[0];
               res(1)(0,k) = wl[2]*RAD(x.crd(0)(0,k))*nrm[1]; 
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
