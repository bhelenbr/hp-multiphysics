/*
 *  rsdl.cpp
 *  lvlset++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_lvlset.h"
#include "../hp_boundary.h"
   
block::ctrl tri_hp_lvlset::rsdl(block::ctrl ctrl_message, int stage) {
   int i,j,n,tind;
   FLT fluxx,fluxy;
   const int NV = 4;
   TinyVector<int,3> v;
   TinyMatrix<FLT,ND,ND> ldcrd;
   TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,NV,ND> du;
   int lgpx = basis::tri(log2p).gpx, lgpn = basis::tri(log2p).gpn;
   FLT  rhorbd0, cjcb, cjcbi, oneminusbeta;
   TinyMatrix<TinyMatrix<FLT,ND,ND>,NV-1,NV-1> visc;
   TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,NV-1,NV-1> cv, df;
   TinyVector<FLT,NV> tres;
   block::ctrl state;
   TinyMatrix<FLT,MXGP,MXGP> rho, mu;
   FLT heavy, delt, length, phidw, signphi, deltw;
   TinyVector<FLT,ND> tang,norm;
   TinyVector<TinyMatrix<FLT,MXGP,MXGP>,ND> phivel;

   if (ctrl_message == block::begin) {
      oneminusbeta = 1.0-sim::beta[stage];
      
      gbl_ptr->res.v(Range(0,nvrtx-1),Range::all()) = 0.0;
      gbl_ptr->res_r.v(Range(0,nvrtx-1),Range::all()) *= oneminusbeta;

      if (basis::tri(log2p).sm) {
         gbl_ptr->res.s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all()) = 0.0;
         gbl_ptr->res_r.s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all()) *= oneminusbeta;
         
         if (basis::tri(log2p).im) {
            gbl_ptr->res.i(Range(0,ntri-1),Range(0,basis::tri(log2p).im-1),Range::all()) = 0.0;
            gbl_ptr->res_r.i(Range(0,ntri-1),Range(0,basis::tri(log2p).im-1),Range::all()) *= oneminusbeta;
         }
      }
      excpt = 0;
   }
   
   switch(excpt) {
      case 0: {
         /* THIS IS FOR DEFORMING MESH STUFF */
         if (ctrl_message != block::advance1) {
            state = tri_hp::rsdl(ctrl_message,stage);
            if (state != block::stop) return(state);
            return(block::advance1);
         }
         else {
            ++excpt;
            ctrl_message = block::begin;
         }
      }
      
      case 1: {
         if (ctrl_message != block::advance1) {
            state = mover->rsdl(ctrl_message);
   
            if (state != block::stop) return(state);
            return(block::advance1);
         }
         else 
            ++excpt;
            ctrl_message = block::begin;
      }
      
      case 2: {
         if (ctrl_message != block::advance1) {
            state = block::stop;
            
            for(i=0;i<nsbd;++i)
               state &= hp_sbdry(i)->rsdl(ctrl_message);
   
            if (state != block::stop) return(state);
            return(block::advance1);
         }
         else 
            ++excpt;
      }
      
      case 3: {
      
          for(tind = 0; tind<ntri;++tind) {
            /* LOAD INDICES OF VERTEX POINTS */
            v = td(tind).vrtx;
         
            /* IF TINFO > -1 IT IS CURVED ELEMENT */
            if (td(tind).info > -1) {
               /* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
               crdtocht(tind);
               
               /* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
               for(n=0;n<ND;++n)
                  basis::tri(log2p).proj_bdry(&cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), &dcrd(n,1)(0,0),MXGP);
            }
            else {
               /* PROJECT VERTEX COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
               for(n=0;n<ND;++n)
                  basis::tri(log2p).proj(vrtx(v(0))(n),vrtx(v(1))(n),vrtx(v(2))(n),&crd(n)(0,0),MXGP);

               /* CALCULATE COORDINATE DERIVATIVES A SIMPLE WAY */
               for(n=0;n<ND;++n) {
                  ldcrd(n,0) = 0.5*(vrtx(v(2))(n) -vrtx(v(1))(n));
                  ldcrd(n,1) = 0.5*(vrtx(v(0))(n) -vrtx(v(1))(n));
               }
            }

            /* CALCULATE MESH VELOCITY */
            for(i=0;i<lgpx;++i) {
               for(j=0;j<lgpn;++j) {
                  mvel(0)(i,j) = sim::bd[0]*(crd(0)(i,j) -dxdt(log2p,tind,0)(i,j));
                  mvel(1)(i,j) = sim::bd[0]*(crd(1)(i,j) -dxdt(log2p,tind,1)(i,j));
                }
            }

            /* LOAD SOLUTION COEFFICIENTS FOR THIS ELEMENT */
            /* PROJECT SOLUTION TO GAUSS POINTS WITH DERIVATIVES IF NEEDED FOR VISCOUS TERMS */
            ugtouht(tind);
            if (sim::beta[stage] > 0.0) {
               for(n=0;n<NV-1;++n)
                  basis::tri(log2p).proj(&uht(n)(0),&u(n)(0,0),&du(n,0)(0,0),&du(n,1)(0,0),MXGP);
               basis::tri(log2p).proj(&uht(NV-1)(0),&u(NV-1)(0,0),MXGP);
            }
            else {
               for(n=0;n<ND;++n)
                  basis::tri(log2p).proj(&uht(n)(0),&u(n)(0,0),MXGP);
               basis::tri(log2p).proj(&uht(NV-2)(0),&u(NV-2)(0,0),&du(NV-2,0)(0,0),&du(NV-2,1)(0,0),MXGP);
               basis::tri(log2p).proj(&uht(NV-1)(0),&u(NV-1)(0,0),MXGP);
            }
            
            /* lf IS WHERE I WILL STORE THE ELEMENT RESIDUAL */
            for(n=0;n<NV;++n)
               for(i=0;i<basis::tri(log2p).tm;++i)
                  lf(n)(i) = 0.0;

            if (td(tind).info > -1) {
               /* CURVED ELEMENT */
               /* CONVECTIVE TERMS (IMAGINARY FIRST)*/
               for(i=0;i<lgpx;++i) {
                  for(j=0;j<lgpn;++j) {
                  
                     if (u(2)(i,j) < -gbl_ptr->width) {
                        rho(i,j) = gbl_ptr->rho;
                        mu(i,j)  = gbl_ptr->mu;
                        
                        fluxx = rho(i,j)*RAD(crd(0)(i,j))*(u(0)(i,j) -mvel(0)(i,j));
                        fluxy = rho(i,j)*RAD(crd(0)(i,j))*(u(1)(i,j) -mvel(1)(i,j));
                        
                        /* CONTINUITY EQUATION FLUXES */
                        du(NV-1,0)(i,j) = +dcrd(1,1)(i,j)*fluxx -dcrd(0,1)(i,j)*fluxy;
                        du(NV-1,1)(i,j) = -dcrd(1,0)(i,j)*fluxx +dcrd(0,0)(i,j)*fluxy;
                       
                        /* CONVECTIVE FLUXES */
                        for(n=0;n<NV-1;++n) {
                           cv(n,0)(i,j) = u(n)(i,j)*du(NV-1,0)(i,j);
                           cv(n,1)(i,j) = u(n)(i,j)*du(NV-1,1)(i,j);
                        }
                     }
                     else if (u(2)(i,j) > gbl_ptr->width) {
                        rho(i,j) = gbl_ptr->rho2;
                        mu(i,j)  = gbl_ptr->mu2;
                        
                        fluxx = rho(i,j)*RAD(crd(0)(i,j))*(u(0)(i,j) -mvel(0)(i,j));
                        fluxy = rho(i,j)*RAD(crd(0)(i,j))*(u(1)(i,j) -mvel(1)(i,j));
                        
                        /* CONTINUITY EQUATION FLUXES */
                        du(NV-1,0)(i,j) = +dcrd(1,1)(i,j)*fluxx -dcrd(0,1)(i,j)*fluxy;
                        du(NV-1,1)(i,j) = -dcrd(1,0)(i,j)*fluxx +dcrd(0,0)(i,j)*fluxy;
                       
                        /* CONVECTIVE FLUXES */
                        for(n=0;n<NV-1;++n) {
                           cv(n,0)(i,j) = u(n)(i,j)*du(NV-1,0)(i,j);
                           cv(n,1)(i,j) = u(n)(i,j)*du(NV-1,1)(i,j);
                        }
                        
                     }
                     else {
                        phidw = u(2)(i,j)/gbl_ptr->width;
                        heavy = heavyside(phidw);
                        delt = delta(phidw);
                        cjcb = dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j);
                        tang(0) = (-dcrd(0,1)(i,j)*du(2,0)(i,j) +dcrd(0,0)(i,j)*du(2,1)(i,j))/cjcb;
                        tang(1) = (-dcrd(1,1)(i,j)*du(2,0)(i,j) +dcrd(1,0)(i,j)*du(2,1)(i,j))/cjcb;
                        length = sqrt(tang(0)*tang(0) +tang(1)*tang(1));
                        tang /= length;
                        rho(i,j) = gbl_ptr->rho +(gbl_ptr->rho2 -gbl_ptr->rho)*heavy;
                        mu(i,j) = gbl_ptr->mu +(gbl_ptr->mu2 -gbl_ptr->mu)*heavy;
                        
                       
                        fluxx = rho(i,j)*RAD(crd(0)(i,j))*(u(0)(i,j) -mvel(0)(i,j));
                        fluxy = rho(i,j)*RAD(crd(0)(i,j))*(u(1)(i,j) -mvel(1)(i,j));
                        
                        /* CONTINUITY EQUATION FLUXES */
                        du(NV-1,0)(i,j) = +dcrd(1,1)(i,j)*fluxx -dcrd(0,1)(i,j)*fluxy;
                        du(NV-1,1)(i,j) = -dcrd(1,0)(i,j)*fluxx +dcrd(0,0)(i,j)*fluxy;
                       
                        /* CONVECTIVE FLUXES */
                        for(n=0;n<NV-1;++n) {
                           cv(n,0)(i,j) = u(n)(i,j)*du(NV-1,0)(i,j);
                           cv(n,1)(i,j) = u(n)(i,j)*du(NV-1,1)(i,j);
                        }
                        
                        /* SURFACE TENSION TERMS */
                        cv(0,0)(i,j) += -delt*gbl_ptr->sigma*tang(0)*length*(dcrd(1,1)(i,j)*tang(0) -dcrd(0,1)(i,j)*tang(1));
                        cv(0,1)(i,j) += +delt*gbl_ptr->sigma*tang(0)*length*(dcrd(1,0)(i,j)*tang(0) -dcrd(0,0)(i,j)*tang(1));
                        cv(1,0)(i,j) += -delt*gbl_ptr->sigma*tang(1)*length*(dcrd(1,1)(i,j)*tang(0) -dcrd(0,1)(i,j)*tang(1));
                        cv(1,1)(i,j) += +delt*gbl_ptr->sigma*tang(1)*length*(dcrd(1,0)(i,j)*tang(0) -dcrd(0,0)(i,j)*tang(1));
                     }

                     /* PRESSURE TERMS */
                     /* U-MOMENTUM */
                     cv(0,0)(i,j) += dcrd(1,1)(i,j)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
                     cv(0,1)(i,j) -= dcrd(1,0)(i,j)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
                     /* V-MOMENTUM */
                     cv(1,0)(i,j) -=  dcrd(0,1)(i,j)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
                     cv(1,1)(i,j) +=  dcrd(0,0)(i,j)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
                  }
               }
               for(n=0;n<NV-1;++n)
                  basis::tri(log2p).intgrtrs(&lf(n)(0),&cv(n,0)(0,0),&cv(n,1)(0,0),MXGP);
               basis::tri(log2p).intgrtrs(&lf(NV-1)(0),&du(NV-1,0)(0,0),&du(NV-1,1)(0,0),MXGP);

               /* ASSEMBLE GLOBAL FORCING (IMAGINARY TERMS) */
               lftog(tind,gbl_ptr->res);

               /* NEGATIVE REAL TERMS */
               if (sim::beta[stage] > 0.0) {
                  /* TIME DERIVATIVE TERMS */ 
                  for(i=0;i<lgpx;++i) {
                     for(j=0;j<lgpn;++j) {
                        cjcb = dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j);
                        rhorbd0 = rho(i,j)*sim::bd[0]*RAD(crd(0)(i,j))*cjcb;
                        cjcbi = mu(i,j)*RAD(crd(0)(i,j))/cjcb;
                        
                        /* UNSTEADY TERMS */
                        for(n=0;n<NV-1;++n)
                           res(n)(i,j) = rhorbd0*u(n)(i,j) +dugdt(log2p,tind,n)(i,j);
                        res(NV-1)(i,j) = rhorbd0 +dugdt(log2p,tind,NV-1)(i,j);
#ifdef AXISYMMETRIC
                        res(0)(i,j) -= cjcb*(u(NV-1)(i,j) -2.*mu(i,j)*u(0)(i,j)/crd(0)(i,j));
#endif

#ifdef BODYFORCE
                        res(0)(i,j) -= gbl_ptr->rho*RAD(crd(0)(i,j))*cjcb*sim::body(0);
                        res(1)(i,j) -= gbl_ptr->rho*RAD(crd(0)(i,j))*cjcb*sim::body(1);
#endif

                        /* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
                        /* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
                        visc(0,0)(0,0) = -cjcbi*(2.*dcrd(1,1)(i,j)*dcrd(1,1)(i,j) +dcrd(0,1)(i,j)*dcrd(0,1)(i,j));
                        visc(0,0)(1,1) = -cjcbi*(2.*dcrd(1,0)(i,j)*dcrd(1,0)(i,j) +dcrd(0,0)(i,j)*dcrd(0,0)(i,j));
                        visc(0,0)(0,1) =  cjcbi*(2.*dcrd(1,1)(i,j)*dcrd(1,0)(i,j) +dcrd(0,1)(i,j)*dcrd(0,0)(i,j));
#define                 viscI0II0II1II0I visc(0,0)(0,1)

                        visc(1,1)(0,0) = -cjcbi*(dcrd(1,1)(i,j)*dcrd(1,1)(i,j) +2.*dcrd(0,1)(i,j)*dcrd(0,1)(i,j));
                        visc(1,1)(1,1) = -cjcbi*(dcrd(1,0)(i,j)*dcrd(1,0)(i,j) +2.*dcrd(0,0)(i,j)*dcrd(0,0)(i,j));
                        visc(1,1)(0,1) =  cjcbi*(dcrd(1,1)(i,j)*dcrd(1,0)(i,j) +2.*dcrd(0,1)(i,j)*dcrd(0,0)(i,j));
#define                 viscI1II1II1II0I visc(1,1)(0,1)
                        
                        visc(0,1)(0,0) =  cjcbi*dcrd(0,1)(i,j)*dcrd(1,1)(i,j);
                        visc(0,1)(1,1) =  cjcbi*dcrd(0,0)(i,j)*dcrd(1,0)(i,j);
                        visc(0,1)(0,1) = -cjcbi*dcrd(0,1)(i,j)*dcrd(1,0)(i,j);
                        visc(0,1)(1,0) = -cjcbi*dcrd(0,0)(i,j)*dcrd(1,1)(i,j);

                           /* OTHER SYMMETRIES    */            
#define                 viscI1II0II0II0I visc(0,1)(0,0)
#define                 viscI1II0II1II1I visc(0,1)(1,1)
#define                 viscI1II0II0II1I visc(0,1)(1,0)
#define                 viscI1II0II1II0I visc(0,1)(0,1)                        


                        df(0,0)(i,j) = +visc(0,0)(0,0)*du(0,0)(i,j) +visc(0,1)(0,0)*du(1,0)(i,j)
                                    +visc(0,0)(0,1)*du(0,1)(i,j) +visc(0,1)(0,1)*du(1,1)(i,j);

                        df(0,1)(i,j) = +viscI0II0II1II0I*du(0,0)(i,j) +visc(0,1)(1,0)*du(1,0)(i,j)
                                    +visc(0,0)(1,1)*du(0,1)(i,j) +visc(0,1)(1,1)*du(1,1)(i,j);

                        df(1,0)(i,j) = +viscI1II0II0II0I*du(0,0)(i,j) +visc(1,1)(0,0)*du(1,0)(i,j)
                                    +viscI1II0II0II1I*du(0,1)(i,j) +visc(1,1)(0,1)*du(1,1)(i,j);

                        df(1,1)(i,j) = +viscI1II0II1II0I*du(0,0)(i,j) +viscI1II1II1II0I*du(1,0)(i,j)
                                    +viscI1II0II1II1I*du(0,1)(i,j) +visc(1,1)(1,1)*du(1,1)(i,j);
                                    
                        for(n=0;n<ND;++n) {
                           cv(n,0)(i,j) += df(n,0)(i,j);
                           cv(n,1)(i,j) += df(n,1)(i,j);
                        }
                      }
                  }
                  for(n=0;n<NV;++n)
                     basis::tri(log2p).intgrt(&lf(n)(0),&res(n)(0,0),MXGP);

                  /* CALCULATE RESIDUAL TO GOVERNING EQUATION & STORE IN RES */
                  for(n=0;n<ND;++n) {
                     basis::tri(log2p).derivr(&cv(n,0)(0,0),&res(n)(0,0),MXGP);
                     basis::tri(log2p).derivs(&cv(n,1)(0,0),&res(n)(0,0),MXGP);
                  }
                  basis::tri(log2p).derivr(&du(NV-1,0)(0,0),&res(NV-1)(0,0),MXGP);
                  basis::tri(log2p).derivs(&du(NV-1,1)(0,0),&res(NV-1)(0,0),MXGP);

                  /* THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES */
                  for(i=0;i<lgpx;++i) {
                     for(j=0;j<lgpn;++j) {
                                 
                        tres(0) = gbl_ptr->tau(tind,0)*res(0)(i,j);
                        tres(1) = gbl_ptr->tau(tind,0)*res(1)(i,j);
                        tres(2) = gbl_ptr->tau(tind,2)*res(2)(i,j);
                        tres(NV-1) = gbl_ptr->tau(tind,NV-1)*res(NV-1)(i,j);

                        df(0,0)(i,j) -= (dcrd(1,1)(i,j)*(2*u(0)(i,j)-mvel(0)(i,j))
                                    -dcrd(0,1)(i,j)*(u(1)(i,j)-mvel(1)(i,j)))*tres(0)
                                    -dcrd(0,1)(i,j)*u(0)(i,j)*tres(1)
                                    +dcrd(1,1)(i,j)*tres(NV-1);
                        df(0,1)(i,j) -= (-dcrd(1,0)(i,j)*(2*u(0)(i,j)-mvel(0)(i,j))
                                    +dcrd(0,0)(i,j)*(u(1)(i,j)-mvel(1)(i,j)))*tres(0)
                                    +dcrd(0,0)(i,j)*u(0)(i,j)*tres(1)
                                    -dcrd(1,0)(i,j)*tres(NV-1);
                        df(1,0)(i,j) -= +dcrd(1,1)(i,j)*u(1)(i,j)*tres(0)
                                    +(dcrd(1,1)(i,j)*(u(0)(i,j)-mvel(0)(i,j))
                                    -dcrd(0,1)(i,j)*(2.*u(1)(i,j)-mvel(1)(i,j)))*tres(1)
                                    -dcrd(0,1)(i,j)*tres(NV-1);
                        df(1,1)(i,j) -= -dcrd(1,0)(i,j)*u(1)(i,j)*tres(0)
                              +(-dcrd(1,0)(i,j)*(u(0)(i,j)-mvel(0)(i,j))
                                    +dcrd(0,0)(i,j)*(2.*u(1)(i,j)-mvel(1)(i,j)))*tres(1)
                                    +dcrd(0,0)(i,j)*tres(NV-1);
                                    
                        df(2,0)(i,j) = -(dcrd(1,1)(i,j)*(u(0)(i,j)-mvel(0)(i,j))
                                 -dcrd(0,1)(i,j)*(u(1)(i,j)-mvel(1)(i,j)))*tres(2);
                        df(2,1)(i,j) = -(-dcrd(1,0)(i,j)*(u(0)(i,j)-mvel(0)(i,j))
                                 +dcrd(0,0)(i,j)*(u(1)(i,j)-mvel(1)(i,j)))*tres(2);

                                    
                        du(NV-1,0)(i,j) = -(dcrd(1,1)(i,j)*tres(0) -dcrd(0,1)(i,j)*tres(1));
                        du(NV-1,1)(i,j) = -(-dcrd(1,0)(i,j)*tres(0) +dcrd(0,0)(i,j)*tres(1));
                     }
                  }
                  for(n=0;n<NV-1;++n)
                     basis::tri(log2p).intgrtrs(&lf(n)(0),&df(n,0)(0,0),&df(n,1)(0,0),MXGP);
                  basis::tri(log2p).intgrtrs(&lf(NV-1)(0),&du(NV-1,0)(0,0),&du(NV-1,1)(0,0),MXGP);
                 
                  for(n=0;n<NV;++n)
                     for(i=0;i<basis::tri(log2p).tm;++i)
                        lf(n)(i) *= sim::beta[stage];
                        
                  lftog(tind,gbl_ptr->res_r);
               }
            }
            else {
               /* LINEAR ELEMENT */
               /* CONVECTIVE TERMS (IMAGINARY FIRST)*/
               for(i=0;i<lgpx;++i) {
                  for(j=0;j<lgpn;++j) {
                  
                     /* STUFF FOR LEVEL SET */
                     phidw = u(2)(i,j)/gbl_ptr->width;
                     cjcb = ldcrd(0,0)*ldcrd(1,1) -ldcrd(1,0)*ldcrd(0,1);
                     norm(0) = (+ldcrd(1,1)*du(2,0)(i,j) -ldcrd(1,0)*du(2,1)(i,j))/cjcb;
                     norm(1) = (-ldcrd(0,1)*du(2,0)(i,j) +ldcrd(0,0)*du(2,1)(i,j))/cjcb;
                     length = sqrt(norm(0)*norm(0) +norm(1)*norm(1));
                     norm /= length;
                     if (phidw < -1.0) {
                        rho(i,j) = gbl_ptr->rho;
                        mu(i,j)  = gbl_ptr->mu;
                        res(2)(i,j) = (1.-length)*cjcb;
                        phivel(0)(i,j) = -norm(0);
                        phivel(1)(i,j) = -norm(1);
                        
                        fluxx = rho(i,j)*RAD(crd(0)(i,j))*(u(0)(i,j) -mvel(0)(i,j));
                        fluxy = rho(i,j)*RAD(crd(0)(i,j))*(u(1)(i,j) -mvel(1)(i,j));
                        
                        /* CONTINUITY EQUATION FLUXES */
                        du(NV-1,0)(i,j) = +ldcrd(1,1)*fluxx -ldcrd(0,1)*fluxy;
                        du(NV-1,1)(i,j) = -ldcrd(1,0)*fluxx +ldcrd(0,0)*fluxy;
                       
                        /* CONVECTIVE FLUXES */
                        for(n=0;n<ND;++n) {
                           cv(n,0)(i,j) = u(n)(i,j)*du(NV-1,0)(i,j);
                           cv(n,1)(i,j) = u(n)(i,j)*du(NV-1,1)(i,j);
                        }

                     }
                     else if (phidw > 1.0) {
                        rho(i,j) = gbl_ptr->rho2;
                        mu(i,j)  = gbl_ptr->mu2;
                        res(2)(i,j) = (length-1.)*cjcb;
                        phivel(0)(i,j) = +norm(0);
                        phivel(1)(i,j) = +norm(1);
                        
                        fluxx = rho(i,j)*RAD(crd(0)(i,j))*(u(0)(i,j) -mvel(0)(i,j));
                        fluxy = rho(i,j)*RAD(crd(0)(i,j))*(u(1)(i,j) -mvel(1)(i,j));
                        
                        /* CONTINUITY EQUATION FLUXES */
                        du(NV-1,0)(i,j) = +ldcrd(1,1)*fluxx -ldcrd(0,1)*fluxy;
                        du(NV-1,1)(i,j) = -ldcrd(1,0)*fluxx +ldcrd(0,0)*fluxy;
                       
                        /* CONVECTIVE FLUXES */
                        for(n=0;n<ND;++n) {
                           cv(n,0)(i,j) = u(n)(i,j)*du(NV-1,0)(i,j);
                           cv(n,1)(i,j) = u(n)(i,j)*du(NV-1,1)(i,j);
                        }

                     }
                     else {
                        //signphi = phidw/(fabs(phidw)+0.01);
                        signphi = (phidw > 0.0 ? 1.0 : -1.0);
                        heavy = heavyside(phidw);
                        delt = delta(phidw);
                        rho(i,j) = gbl_ptr->rho +(gbl_ptr->rho2 -gbl_ptr->rho)*heavy;
                        mu(i,j) = gbl_ptr->mu +(gbl_ptr->mu2 -gbl_ptr->mu)*heavy;
                                                
                        deltw = gbl_ptr->width*delt;
                        phivel(0)(i,j) = deltw*(u(0)(i,j)-mvel(0)(i,j)) +(1.0-deltw)*signphi*norm(0);
                        phivel(1)(i,j) = deltw*(u(1)(i,j)-mvel(1)(i,j)) +(1.0-deltw)*signphi*norm(1);
                        
                        res(2)(i,j) = deltw*(sim::bd[0]*cjcb*u(2)(i,j) +dugdt(log2p,tind,2)(i,j)) -(1.0-deltw)*signphi*cjcb
                                    +(phivel(0)(i,j)*norm(0) +phivel(1)(i,j)*norm(1))*length*cjcb;

                        fluxx = rho(i,j)*RAD(crd(0)(i,j))*(u(0)(i,j) -mvel(0)(i,j));
                        fluxy = rho(i,j)*RAD(crd(0)(i,j))*(u(1)(i,j) -mvel(1)(i,j));
                        
                        /* CONTINUITY EQUATION FLUXES */
                        du(NV-1,0)(i,j) = +ldcrd(1,1)*fluxx -ldcrd(0,1)*fluxy;
                        du(NV-1,1)(i,j) = -ldcrd(1,0)*fluxx +ldcrd(0,0)*fluxy;
                       
                        /* CONVECTIVE FLUXES */
                        for(n=0;n<ND;++n) {
                           cv(n,0)(i,j) = u(n)(i,j)*du(NV-1,0)(i,j);
                           cv(n,1)(i,j) = u(n)(i,j)*du(NV-1,1)(i,j);
                        }
                        
                        /* SURFACE TENSION TERMS */
                        cv(0,0)(i,j) += -delt*gbl_ptr->sigma*norm(1)*length*(ldcrd(1,1)*norm(1) +ldcrd(0,1)*norm(0));
                        cv(0,1)(i,j) += +delt*gbl_ptr->sigma*norm(1)*length*(ldcrd(1,0)*norm(1) +ldcrd(0,0)*norm(0));
                        cv(1,0)(i,j) += +delt*gbl_ptr->sigma*norm(0)*length*(ldcrd(1,1)*norm(1) +ldcrd(0,1)*norm(0));
                        cv(1,1)(i,j) += -delt*gbl_ptr->sigma*norm(0)*length*(ldcrd(1,0)*norm(1) +ldcrd(0,0)*norm(0));
                     }
                       
                     /* PRESSURE TERMS */
                     /* U-MOMENTUM */
                     cv(0,0)(i,j) += ldcrd(1,1)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
                     cv(0,1)(i,j) -= ldcrd(1,0)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
                     /* V-MOMENTUM */
                     cv(1,0)(i,j) -=  ldcrd(0,1)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
                     cv(1,1)(i,j) +=  ldcrd(0,0)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
                  }
               }
               for(n=0;n<ND;++n)
                  basis::tri(log2p).intgrtrs(&lf(n)(0),&cv(n,0)(0,0),&cv(n,1)(0,0),MXGP);
               basis::tri(log2p).intgrtrs(&lf(NV-1)(0),&du(NV-1,0)(0,0),&du(NV-1,1)(0,0),MXGP);
               

               /* ASSEMBLE GLOBAL FORCING (IMAGINARY TERMS) */
               lftog(tind,gbl_ptr->res);

               /* NEGATIVE REAL TERMS */
               if (sim::beta[stage] > 0.0) {
                  cjcb = ldcrd(0,0)*ldcrd(1,1) -ldcrd(1,0)*ldcrd(0,1);
                  cjcbi = 1./cjcb;
                  
                  /* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
                  /* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
                  visc(0,0)(0,0) = -cjcbi*(2.*ldcrd(1,1)*ldcrd(1,1) +ldcrd(0,1)*ldcrd(0,1));
                  visc(0,0)(1,1) = -cjcbi*(2.*ldcrd(1,0)*ldcrd(1,0) +ldcrd(0,0)*ldcrd(0,0));
                  visc(0,0)(0,1) =  cjcbi*(2.*ldcrd(1,1)*ldcrd(1,0) +ldcrd(0,1)*ldcrd(0,0));
#define           viscI0II0II1II0I visc(0,0)(0,1)

                  visc(1,1)(0,0) = -cjcbi*(ldcrd(1,1)*ldcrd(1,1) +2.*ldcrd(0,1)*ldcrd(0,1));
                  visc(1,1)(1,1) = -cjcbi*(ldcrd(1,0)*ldcrd(1,0) +2.*ldcrd(0,0)*ldcrd(0,0));
                  visc(1,1)(0,1) =  cjcbi*(ldcrd(1,1)*ldcrd(1,0) +2.*ldcrd(0,1)*ldcrd(0,0));
#define           viscI1II1II1II0I visc(1,1)(0,1)
                  
                  visc(0,1)(0,0) =  cjcbi*ldcrd(0,1)*ldcrd(1,1);
                  visc(0,1)(1,1) =  cjcbi*ldcrd(0,0)*ldcrd(1,0);
                  visc(0,1)(0,1) = -cjcbi*ldcrd(0,1)*ldcrd(1,0);
                  visc(0,1)(1,0) = -cjcbi*ldcrd(0,0)*ldcrd(1,1);

                  /* OTHER SYMMETRIES    */            
#define           viscI1II0II0II0I visc(0,1)(0,0)
#define           viscI1II0II1II1I visc(0,1)(1,1)
#define           viscI1II0II0II1I visc(0,1)(1,0)
#define           viscI1II0II1II0I visc(0,1)(0,1)

                  /* TIME DERIVATIVE TERMS */ 
                  for(i=0;i<lgpx;++i) {
                     for(j=0;j<lgpn;++j) {
                        rhorbd0 = RAD(crd(0)(i,j))*rho(i,j)*sim::bd[0]*cjcb;
                        
                        /* UNSTEADY TERMS */
                        for(n=0;n<ND;++n)
                           res(n)(i,j) = rhorbd0*u(n)(i,j) +dugdt(log2p,tind,n)(i,j);
                        res(NV-1)(i,j) = rhorbd0 +dugdt(log2p,tind,NV-1)(i,j);
            
#ifdef AXISYMMETRIC
                        res(0)(i,j) -= cjcb*(u(2)(i,j) -2.*mu(i,j)*u(0)(i,j)/crd(0)(i,j));
#endif

#ifdef BODYFORCE
                        res(0)(i,j) -= gbl_ptr->rho*RAD(crd(0)(i,j))*cjcb*sim::body(0);
                        res(1)(i,j) -= gbl_ptr->rho*RAD(crd(0)(i,j))*cjcb*sim::body(1);
#endif
                        df(0,0)(i,j) = mu(i,j)*RAD(crd(0)(i,j))*(+visc(0,0)(0,0)*du(0,0)(i,j) +visc(0,1)(0,0)*du(1,0)(i,j)
                                              +visc(0,0)(0,1)*du(0,1)(i,j) +visc(0,1)(0,1)*du(1,1)(i,j));

                        df(0,1)(i,j) = mu(i,j)*RAD(crd(0)(i,j))*(+viscI0II0II1II0I*du(0,0)(i,j) +visc(0,1)(1,0)*du(1,0)(i,j)
                                              +visc(0,0)(1,1)*du(0,1)(i,j) +visc(0,1)(1,1)*du(1,1)(i,j));

                        df(1,0)(i,j) = mu(i,j)*RAD(crd(0)(i,j))*(+viscI1II0II0II0I*du(0,0)(i,j) +visc(1,1)(0,0)*du(1,0)(i,j)
                                              +viscI1II0II0II1I*du(0,1)(i,j) +visc(1,1)(0,1)*du(1,1)(i,j));

                        df(1,1)(i,j) = mu(i,j)*RAD(crd(0)(i,j))*(+viscI1II0II1II0I*du(0,0)(i,j) +viscI1II1II1II0I*du(1,0)(i,j)
                                              +viscI1II0II1II1I*du(0,1)(i,j) +visc(1,1)(1,1)*du(1,1)(i,j));
                        
                        for(n=0;n<ND;++n) {
                           cv(n,0)(i,j) += df(n,0)(i,j);
                           cv(n,1)(i,j) += df(n,1)(i,j);
                        } 
                      }
                  }
                  for(n=0;n<NV;++n)
                     basis::tri(log2p).intgrt(&lf(n)(0),&res(n)(0,0),MXGP);

                  /* CALCULATE RESIDUAL TO GOVERNING EQUATION & STORE IN RES */
                  for(n=0;n<ND;++n) {
                     basis::tri(log2p).derivr(&cv(n,0)(0,0),&res(n)(0,0),MXGP);
                     basis::tri(log2p).derivs(&cv(n,1)(0,0),&res(n)(0,0),MXGP);
                  }
                  basis::tri(log2p).derivr(&du(NV-1,0)(0,0),&res(NV-1)(0,0),MXGP);
                  basis::tri(log2p).derivs(&du(NV-1,1)(0,0),&res(NV-1)(0,0),MXGP);
                  
                  /* THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES */
                  for(i=0;i<lgpx;++i) {
                     for(j=0;j<lgpn;++j) {
                        tres(0) = gbl_ptr->tau(tind,0)*res(0)(i,j);
                        tres(1) = gbl_ptr->tau(tind,0)*res(1)(i,j);
                        tres(2) = gbl_ptr->tau(tind,2)*res(2)(i,j);
                        tres(NV-1) = gbl_ptr->tau(tind,NV-1)*res(NV-1)(i,j);

                        df(0,0)(i,j) -= (ldcrd(1,1)*(2*u(0)(i,j)-mvel(0)(i,j))
                                    -ldcrd(0,1)*(u(1)(i,j)-mvel(1)(i,j)))*tres(0)
                                    -ldcrd(0,1)*u(0)(i,j)*tres(1)
                                    +ldcrd(1,1)*tres(NV-1);
                        df(0,1)(i,j) -= (-ldcrd(1,0)*(2*u(0)(i,j)-mvel(0)(i,j))
                                    +ldcrd(0,0)*(u(1)(i,j)-mvel(1)(i,j)))*tres(0)
                                    +ldcrd(0,0)*u(0)(i,j)*tres(1)
                                    -ldcrd(1,0)*tres(NV-1);
                        df(1,0)(i,j) -= +ldcrd(1,1)*u(1)(i,j)*tres(0)
                                    +(ldcrd(1,1)*(u(0)(i,j)-mvel(0)(i,j))
                                    -ldcrd(0,1)*(2.*u(1)(i,j)-mvel(1)(i,j)))*tres(1)
                                    -ldcrd(0,1)*tres(NV-1);
                        df(1,1)(i,j) -= -ldcrd(1,0)*u(1)(i,j)*tres(0)
                                    +(-ldcrd(1,0)*(u(0)(i,j)-mvel(0)(i,j))
                                    +ldcrd(0,0)*(2.*u(1)(i,j)-mvel(1)(i,j)))*tres(1)
                                    +ldcrd(0,0)*tres(NV-1);
                                    
                        df(2,0)(i,j) = -(ldcrd(1,1)*phivel(0)(i,j) -ldcrd(0,1)*phivel(1)(i,j))*tres(2);
                        df(2,1)(i,j) = -(-ldcrd(1,0)*phivel(0)(i,j) +ldcrd(0,0)*phivel(1)(i,j))*tres(2);

                        du(NV-1,0)(i,j) = -(ldcrd(1,1)*tres(0) -ldcrd(0,1)*tres(1));
                        du(NV-1,1)(i,j) = -(-ldcrd(1,0)*tres(0) +ldcrd(0,0)*tres(1));
                     }
                  }
                  for(n=0;n<NV-1;++n)
                     basis::tri(log2p).intgrtrs(&lf(n)(0),&df(n,0)(0,0),&df(n,1)(0,0),MXGP);
                  basis::tri(log2p).intgrtrs(&lf(NV-1)(0),&du(NV-1,0)(0,0),&du(NV-1,1)(0,0),MXGP);
                  
                  for(n=0;n<NV;++n)
                     for(i=0;i<basis::tri(log2p).tm;++i)
                        lf(n)(i) *= sim::beta[stage];
                        
                  lftog(tind,gbl_ptr->res_r);
               }
            }
         }

         /* ADD IN VISCOUS/DISSIPATIVE FLUX */
         gbl_ptr->res.v(Range(0,nvrtx-1),Range::all()) += gbl_ptr->res_r.v(Range(0,nvrtx-1),Range::all());
         if (basis::tri(log2p).sm) {
            gbl_ptr->res.s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all()) += gbl_ptr->res_r.s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all());        
            if (basis::tri(log2p).im) {
               gbl_ptr->res.i(Range(0,ntri-1),Range(0,basis::tri(log2p).im-1),Range::all()) += gbl_ptr->res_r.i(Range(0,ntri-1),Range(0,basis::tri(log2p).im-1),Range::all());     
            }
         }

      /*********************************************/
         /* MODIFY RESIDUALS ON COARSER MESHES         */
      /*********************************************/   
         if(coarse) {
         /* CALCULATE DRIVING TERM ON FIRST ENTRY TO COARSE MESH */
            if(isfrst) {
               dres(log2p).v(Range(0,nvrtx-1),Range::all()) = fadd*gbl_ptr->res0.v(Range(0,nvrtx-1),Range::all()) -gbl_ptr->res.v(Range(0,nvrtx-1),Range::all());
               if (basis::tri(log2p).sm) dres(log2p).s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all()) = fadd*gbl_ptr->res0.s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all()) -gbl_ptr->res.s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all());     
               if (basis::tri(log2p).im) dres(log2p).i(Range(0,ntri-1),Range(0,basis::tri(log2p).im-1),Range::all()) = fadd*gbl_ptr->res0.i(Range(0,ntri-1),Range(0,basis::tri(log2p).im-1),Range::all()) -gbl_ptr->res.i(Range(0,ntri-1),Range(0,basis::tri(log2p).im-1),Range::all());
               isfrst = false;
            }
            gbl_ptr->res.v(Range(0,nvrtx-1),Range::all()) += dres(log2p).v(Range(0,nvrtx-1),Range::all()); 
            if (basis::tri(log2p).sm) gbl_ptr->res.s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all()) += dres(log2p).s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all());
            if (basis::tri(log2p).im) gbl_ptr->res.i(Range(0,ntri-1),Range(0,basis::tri(log2p).im-1),Range::all()) += dres(log2p).i(Range(0,ntri-1),Range(0,basis::tri(log2p).im-1),Range::all());  
         }
         
         ++excpt;
         
//         for(i=0;i<nvrtx;++i) {
//            printf("rsdl v: %d ",i);
//            for (n=0;n<NV;++n) 
//               printf("%e ",gbl_ptr->res.v(i,n));
//            printf("\n");
//         }
//            
//         for(i=0;i<nside;++i) {
//            for(int m=0;m<basis::tri(log2p).sm;++m) {
//               printf("rsdl s: %d %d ",i,m); 
//               for(n=0;n<NV;++n)
//                  printf("%e ",gbl_ptr->res.s(i,m,n));
//               printf("\n");
//            }
//         }
//
//         for(i=0;i<ntri;++i) {
//            for(int m=0;m<basis::tri(log2p).im;++m) {
//               printf("rsdl i: %d %d ",i,m);
//               for(n=0;n<NV;++n) 
//                  printf("%e %e %e\n",gbl_ptr->res.i(i,m,n));
//               printf("\n");
//            }
//         }

      }
   }

   return(block::stop);
}
