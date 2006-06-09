/*
 *  rsdl.cpp
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_ins.h"
#include "hp_boundary.h"
   
block::ctrl tri_hp_ins::rsdl(block::ctrl ctrl_message, int stage) {
   int i,j,n,tind;
   FLT fluxx,fluxy;
   const int NV = 3;
   TinyVector<int,3> v;
   TinyMatrix<FLT,ND,ND> ldcrd;
   TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,NV,ND> du;
   int lgpx = basis::tri(log2p).gpx, lgpn = basis::tri(log2p).gpn;
   FLT rhobd0 = ins_gbl->rho*sim::bd[0], lmu = ins_gbl->mu, rhorbd0, lrhorbd0, cjcb, cjcbi, oneminusbeta;
   FLT visc[ND][ND][ND][ND], tres[NV];
   FLT cv00[MXGP][MXGP],cv01[MXGP][MXGP],cv10[MXGP][MXGP],cv11[MXGP][MXGP]; // LOCAL WORK ARRAYS
   FLT e00[MXGP][MXGP],e01[MXGP][MXGP],e10[MXGP][MXGP],e11[MXGP][MXGP]; // LOCAL WORK ARRAYS
   block::ctrl state;
   

   
   if (ctrl_message == block::begin) {
      oneminusbeta = 1.0-sim::beta[stage];
      
      hp_gbl->res.v(Range(0,nvrtx-1),Range::all()) = 0.0;
      hp_gbl->res_r.v(Range(0,nvrtx-1),Range::all()) *= oneminusbeta;

      if (basis::tri(log2p).sm) {
         hp_gbl->res.s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all()) = 0.0;
         hp_gbl->res_r.s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all()) *= oneminusbeta;
         
         if (basis::tri(log2p).im) {
            hp_gbl->res.i(Range(0,ntri-1),Range(0,basis::tri(log2p).im-1),Range::all()) = 0.0;
            hp_gbl->res_r.i(Range(0,ntri-1),Range(0,basis::tri(log2p).im-1),Range::all()) *= oneminusbeta;
         }
      }
      excpt = 0;
   }


   
   switch(excpt) {
      case 0: {
#ifdef CTRL_DEBUG
         *sim::log << "step 0 of ins::rsdl: ctrl_message: " << ctrl_message << " excpt: " << excpt << " stage: " << stage << std::endl;
#endif
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
#ifdef CTRL_DEBUG
         *sim::log << "step 1 of ins::rsdl: ctrl_message: " << ctrl_message << " excpt: " << excpt << " stage: " << stage << std::endl;
#endif
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
#ifdef CTRL_DEBUG
         *sim::log << "step 1 of ins::rsdl: ctrl_message: " << ctrl_message << " excpt: " << excpt << " stage: " << stage << std::endl;
#endif
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
#ifdef CTRL_DEBUG
         *sim::log << "step 2 of ins::rsdl: ctrl_message: " << ctrl_message << " excpt: " << excpt << " stage: " << stage << std::endl;
#endif         
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
#ifdef DROP
                  mvel(0)(i,j) += mesh_ref_vel(0);
                  mvel(1)(i,j) += mesh_ref_vel(1);
#endif
               }
            }

            /* LOAD SOLUTION COEFFICIENTS FOR THIS ELEMENT */
            /* PROJECT SOLUTION TO GAUSS POINTS WITH DERIVATIVES IF NEEDED FOR VISCOUS TERMS */
            ugtouht(tind);
            if (sim::beta[stage] > 0.0) {
               basis::tri(log2p).proj(&uht(0)(0),&u(0)(0,0),&du(0,0)(0,0),&du(0,1)(0,0),MXGP);
               basis::tri(log2p).proj(&uht(1)(0),&u(1)(0,0),&du(1,0)(0,0),&du(1,1)(0,0),MXGP);
               basis::tri(log2p).proj(&uht(2)(0),&u(2)(0,0),MXGP);
            }
            else {
               basis::tri(log2p).proj(&uht(0)(0),&u(0)(0,0),MXGP);
               basis::tri(log2p).proj(&uht(1)(0),&u(1)(0,0),MXGP);
               basis::tri(log2p).proj(&uht(2)(0),&u(2)(0,0),MXGP);
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

                     fluxx = ins_gbl->rho*RAD(crd(0)(i,j))*(u(0)(i,j) -mvel(0)(i,j));
                     fluxy = ins_gbl->rho*RAD(crd(0)(i,j))*(u(1)(i,j) -mvel(1)(i,j));
                     
                     /* CONTINUITY EQUATION FLUXES */
                     du(2,0)(i,j) = +dcrd(1,1)(i,j)*fluxx -dcrd(0,1)(i,j)*fluxy;
                     du(2,1)(i,j) = -dcrd(1,0)(i,j)*fluxx +dcrd(0,0)(i,j)*fluxy;
#ifndef INERTIALESS
                     /* U-MOMENTUM */
                     cv00[i][j] =  u(0)(i,j)*du(2,0)(i,j) +dcrd(1,1)(i,j)*RAD(crd(0)(i,j))*u(2)(i,j);
                     cv01[i][j] =  u(0)(i,j)*du(2,1)(i,j) -dcrd(1,0)(i,j)*RAD(crd(0)(i,j))*u(2)(i,j);
                     /* V-MOMENTUM */
                     cv10[i][j] =  u(1)(i,j)*du(2,0)(i,j) -dcrd(0,1)(i,j)*RAD(crd(0)(i,j))*u(2)(i,j);
                     cv11[i][j] =  u(1)(i,j)*du(2,1)(i,j) +dcrd(0,0)(i,j)*RAD(crd(0)(i,j))*u(2)(i,j);
#else
                     cv00[i][j] =  +dcrd(1,1)(i,j)*RAD(crd(0)(i,j))*u(2)(i,j);
                     cv01[i][j] =  -dcrd(1,0)(i,j)*RAD(crd(0)(i,j))*u(2)(i,j);
                     cv10[i][j] =  -dcrd(0,1)(i,j)*RAD(crd(0)(i,j))*u(2)(i,j);
                     cv11[i][j] =  +dcrd(0,0)(i,j)*RAD(crd(0)(i,j))*u(2)(i,j);
#endif
                  }
               }
               basis::tri(log2p).intgrtrs(&lf(0)(0),cv00[0],cv01[0],MXGP);
               basis::tri(log2p).intgrtrs(&lf(1)(0),cv10[0],cv11[0],MXGP);
               basis::tri(log2p).intgrtrs(&lf(2)(0),&du(2,0)(0,0),&du(2,1)(0,0),MXGP);

               /* ASSEMBLE GLOBAL FORCING (IMAGINARY TERMS) */
               lftog(tind,ins_gbl->res);

               /* NEGATIVE REAL TERMS */
               if (sim::beta[stage] > 0.0) {
                  /* TIME DERIVATIVE TERMS */ 
                  for(i=0;i<lgpx;++i) {
                     for(j=0;j<lgpn;++j) {
                        cjcb = dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j);
                        rhorbd0 = rhobd0*RAD(crd(0)(i,j))*cjcb;
                        cjcbi = lmu*RAD(crd(0)(i,j))/cjcb;
                        
                        /* UNSTEADY TERMS */
                        res(0)(i,j) = rhorbd0*u(0)(i,j) +dugdt(log2p,tind,0)(i,j);
                        res(1)(i,j) = rhorbd0*u(1)(i,j) +dugdt(log2p,tind,1)(i,j);
                        res(2)(i,j) = rhorbd0 +dugdt(log2p,tind,2)(i,j);
#ifdef AXISYMMETRIC
                        res(0)(i,j) -= cjcb*(u(2)(i,j) -2.*lmu*u(0)(i,j)/crd(0)(i,j));
#endif
#ifdef SOURCETERM
                        res(0)(i,j) -= RAD(crd(0)(i,j))*cjcb*amp*pow(crd(0)(i,j)-0.5,amp-1.);
#endif

#ifdef BODY
                        res(0)(i,j) -= ins_gbl->rho*RAD(crd(0)(i,j))*cjcb*body[0];
                        res(1)(i,j) -= ins_gbl->rho*RAD(crd(0)(i,j))*cjcb*body[1];
#ifdef INERTIALESS
                        res(0)(i,j) = -ins_gbl->rho*RAD(crd(0)(i,j))*cjcb*body[0];
                        res(1)(i,j) = -ins_gbl->rho*RAD(crd(0)(i,j))*cjcb*body[1];
#endif
#endif              
                                                
                        /* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
                        /* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
                        visc[0][0][0][0] = -cjcbi*(2.*dcrd(1,1)(i,j)*dcrd(1,1)(i,j) +dcrd(0,1)(i,j)*dcrd(0,1)(i,j));
                        visc[0][0][1][1] = -cjcbi*(2.*dcrd(1,0)(i,j)*dcrd(1,0)(i,j) +dcrd(0,0)(i,j)*dcrd(0,0)(i,j));
                        visc[0][0][0][1] =  cjcbi*(2.*dcrd(1,1)(i,j)*dcrd(1,0)(i,j) +dcrd(0,1)(i,j)*dcrd(0,0)(i,j));
#define                 viscI0II0II1II0I visc[0][0][0][1]

                        visc[1][1][0][0] = -cjcbi*(dcrd(1,1)(i,j)*dcrd(1,1)(i,j) +2.*dcrd(0,1)(i,j)*dcrd(0,1)(i,j));
                        visc[1][1][1][1] = -cjcbi*(dcrd(1,0)(i,j)*dcrd(1,0)(i,j) +2.*dcrd(0,0)(i,j)*dcrd(0,0)(i,j));
                        visc[1][1][0][1] =  cjcbi*(dcrd(1,1)(i,j)*dcrd(1,0)(i,j) +2.*dcrd(0,1)(i,j)*dcrd(0,0)(i,j));
#define                 viscI1II1II1II0I visc[1][1][0][1]
                        
                        visc[0][1][0][0] =  cjcbi*dcrd(0,1)(i,j)*dcrd(1,1)(i,j);
                        visc[0][1][1][1] =  cjcbi*dcrd(0,0)(i,j)*dcrd(1,0)(i,j);
                        visc[0][1][0][1] = -cjcbi*dcrd(0,1)(i,j)*dcrd(1,0)(i,j);
                        visc[0][1][1][0] = -cjcbi*dcrd(0,0)(i,j)*dcrd(1,1)(i,j);

                           /* OTHER SYMMETRIES    */            
#define                 viscI1II0II0II0I visc[0][1][0][0]
#define                 viscI1II0II1II1I visc[0][1][1][1]
#define                 viscI1II0II0II1I visc[0][1][1][0]
#define                 viscI1II0II1II0I visc[0][1][0][1]
                        


                        e00[i][j] = +visc[0][0][0][0]*du(0,0)(i,j) +visc[0][1][0][0]*du(1,0)(i,j)
                                    +visc[0][0][0][1]*du(0,1)(i,j) +visc[0][1][0][1]*du(1,1)(i,j);

                        e01[i][j] = +viscI0II0II1II0I*du(0,0)(i,j) +visc[0][1][1][0]*du(1,0)(i,j)
                                    +visc[0][0][1][1]*du(0,1)(i,j) +visc[0][1][1][1]*du(1,1)(i,j);

                        e10[i][j] = +viscI1II0II0II0I*du(0,0)(i,j) +visc[1][1][0][0]*du(1,0)(i,j)
                                    +viscI1II0II0II1I*du(0,1)(i,j) +visc[1][1][0][1]*du(1,1)(i,j);

                        e11[i][j] = +viscI1II0II1II0I*du(0,0)(i,j) +viscI1II1II1II0I*du(1,0)(i,j)
                                    +viscI1II0II1II1I*du(0,1)(i,j) +visc[1][1][1][1]*du(1,1)(i,j);
                                    
                        cv00[i][j] += e00[i][j];
                        cv01[i][j] += e01[i][j];
                        cv10[i][j] += e10[i][j];
                        cv11[i][j] += e11[i][j];
                        
                      }
                  }
                  basis::tri(log2p).intgrt(&lf(0)(0),&res(0)(0,0),MXGP);
                  basis::tri(log2p).intgrt(&lf(1)(0),&res(1)(0,0),MXGP);
                  basis::tri(log2p).intgrt(&lf(2)(0),&res(2)(0,0),MXGP);

                  /* CALCULATE RESIDUAL TO GOVERNING EQUATION & STORE IN RES */
                  basis::tri(log2p).derivr(cv00[0],&res(0)(0,0),MXGP);
                  basis::tri(log2p).derivs(cv01[0],&res(0)(0,0),MXGP);
                  basis::tri(log2p).derivr(cv10[0],&res(1)(0,0),MXGP);
                  basis::tri(log2p).derivs(cv11[0],&res(1)(0,0),MXGP);
                  basis::tri(log2p).derivr(&du(2,0)(0,0),&res(2)(0,0),MXGP);
                  basis::tri(log2p).derivs(&du(2,1)(0,0),&res(2)(0,0),MXGP);
                  


                  /* THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES */
                  for(i=0;i<lgpx;++i) {
                     for(j=0;j<lgpn;++j) {
                                 
                        tres[0] = ins_gbl->tau(tind)*res(0)(i,j);
                        tres[1] = ins_gbl->tau(tind)*res(1)(i,j);
                        tres[2] = ins_gbl->delt(tind)*res(2)(i,j);

#ifndef INERTIALESS
                        e00[i][j] -= (dcrd(1,1)(i,j)*(2*u(0)(i,j)-mvel(0)(i,j))
                                    -dcrd(0,1)(i,j)*(u(1)(i,j)-mvel(1)(i,j)))*tres[0]
                                    -dcrd(0,1)(i,j)*u(0)(i,j)*tres[1]
                                    +dcrd(1,1)(i,j)*tres[2];
                        e01[i][j] -= (-dcrd(1,0)(i,j)*(2*u(0)(i,j)-mvel(0)(i,j))
                                    +dcrd(0,0)(i,j)*(u(1)(i,j)-mvel(1)(i,j)))*tres[0]
                                    +dcrd(0,0)(i,j)*u(0)(i,j)*tres[1]
                                    -dcrd(1,0)(i,j)*tres[2];
                        e10[i][j] -= +dcrd(1,1)(i,j)*u(1)(i,j)*tres[0]
                                    +(dcrd(1,1)(i,j)*(u(0)(i,j)-mvel(0)(i,j))
                                    -dcrd(0,1)(i,j)*(2.*u(1)(i,j)-mvel(1)(i,j)))*tres[1]
                                    -dcrd(0,1)(i,j)*tres[2];
                        e11[i][j] -= -dcrd(1,0)(i,j)*u(1)(i,j)*tres[0]
                                    +(-dcrd(1,0)(i,j)*(u(0)(i,j)-mvel(0)(i,j))
                                    +dcrd(0,0)(i,j)*(2.*u(1)(i,j)-mvel(1)(i,j)))*tres[1]
                                    +dcrd(0,0)(i,j)*tres[2];
#endif
                                    
                        du(2,0)(i,j) = -(dcrd(1,1)(i,j)*tres[0] -dcrd(0,1)(i,j)*tres[1]);
                        du(2,1)(i,j) = -(-dcrd(1,0)(i,j)*tres[0] +dcrd(0,0)(i,j)*tres[1]);
                     }
                  }
                  basis::tri(log2p).intgrtrs(&lf(0)(0),e00[0],e01[0],MXGP);
                  basis::tri(log2p).intgrtrs(&lf(1)(0),e10[0],e11[0],MXGP);
                  basis::tri(log2p).intgrtrs(&lf(2)(0),&du(2,0)(0,0),&du(2,1)(0,0),MXGP);

                  for(n=0;n<NV;++n)
                     for(i=0;i<basis::tri(log2p).tm;++i)
                        lf(n)(i) *= sim::beta[stage];
                        
                  lftog(tind,ins_gbl->res_r);
               }
            }
            else {
               /* LINEAR ELEMENT */
               /* CONVECTIVE TERMS (IMAGINARY FIRST)*/
               for(i=0;i<lgpx;++i) {
                  for(j=0;j<lgpn;++j) {

                     fluxx = ins_gbl->rho*RAD(crd(0)(i,j))*(u(0)(i,j) -mvel(0)(i,j));
                     fluxy = ins_gbl->rho*RAD(crd(0)(i,j))*(u(1)(i,j) -mvel(1)(i,j));
                     
                     /* CONTINUITY EQUATION FLUXES */
                     du(2,0)(i,j) = +ldcrd(1,1)*fluxx -ldcrd(0,1)*fluxy;
                     du(2,1)(i,j) = -ldcrd(1,0)*fluxx +ldcrd(0,0)*fluxy;
#ifndef INERTIALESS
                     /* U-MOMENTUM */
                     cv00[i][j] =  u(0)(i,j)*du(2,0)(i,j) +ldcrd(1,1)*RAD(crd(0)(i,j))*u(2)(i,j);
                     cv01[i][j] =  u(0)(i,j)*du(2,1)(i,j) -ldcrd(1,0)*RAD(crd(0)(i,j))*u(2)(i,j);
                     /* V-MOMENTUM */
                     cv10[i][j] =  u(1)(i,j)*du(2,0)(i,j) -ldcrd(0,1)*RAD(crd(0)(i,j))*u(2)(i,j);
                     cv11[i][j] =  u(1)(i,j)*du(2,1)(i,j) +ldcrd(0,0)*RAD(crd(0)(i,j))*u(2)(i,j);
#else
                     cv00[i][j] =  +ldcrd(1,1)*RAD(crd(0)(i,j))*u(2)(i,j);
                     cv01[i][j] =  -ldcrd(1,0)*RAD(crd(0)(i,j))*u(2)(i,j);
                     cv10[i][j] =  -ldcrd(0,1)*RAD(crd(0)(i,j))*u(2)(i,j);
                     cv11[i][j] =  +ldcrd(0,0)*RAD(crd(0)(i,j))*u(2)(i,j);
#endif
                  }
               }
               basis::tri(log2p).intgrtrs(&lf(0)(0),cv00[0],cv01[0],MXGP);
               basis::tri(log2p).intgrtrs(&lf(1)(0),cv10[0],cv11[0],MXGP);
               basis::tri(log2p).intgrtrs(&lf(2)(0),&du(2,0)(0,0),&du(2,1)(0,0),MXGP);
               

               /* ASSEMBLE GLOBAL FORCING (IMAGINARY TERMS) */
               lftog(tind,ins_gbl->res);

               /* NEGATIVE REAL TERMS */
               if (sim::beta[stage] > 0.0) {
                  cjcb = ldcrd(0,0)*ldcrd(1,1) -ldcrd(1,0)*ldcrd(0,1);
                  cjcbi = lmu/cjcb;
                  lrhorbd0 = rhobd0*cjcb;
                  
                  /* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
                  /* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
                  visc[0][0][0][0] = -cjcbi*(2.*ldcrd(1,1)*ldcrd(1,1) +ldcrd(0,1)*ldcrd(0,1));
                  visc[0][0][1][1] = -cjcbi*(2.*ldcrd(1,0)*ldcrd(1,0) +ldcrd(0,0)*ldcrd(0,0));
                  visc[0][0][0][1] =  cjcbi*(2.*ldcrd(1,1)*ldcrd(1,0) +ldcrd(0,1)*ldcrd(0,0));
#define           viscI0II0II1II0I visc[0][0][0][1]

                  visc[1][1][0][0] = -cjcbi*(ldcrd(1,1)*ldcrd(1,1) +2.*ldcrd(0,1)*ldcrd(0,1));
                  visc[1][1][1][1] = -cjcbi*(ldcrd(1,0)*ldcrd(1,0) +2.*ldcrd(0,0)*ldcrd(0,0));
                  visc[1][1][0][1] =  cjcbi*(ldcrd(1,1)*ldcrd(1,0) +2.*ldcrd(0,1)*ldcrd(0,0));
#define           viscI1II1II1II0I visc[1][1][0][1]
                  
                  visc[0][1][0][0] =  cjcbi*ldcrd(0,1)*ldcrd(1,1);
                  visc[0][1][1][1] =  cjcbi*ldcrd(0,0)*ldcrd(1,0);
                  visc[0][1][0][1] = -cjcbi*ldcrd(0,1)*ldcrd(1,0);
                  visc[0][1][1][0] = -cjcbi*ldcrd(0,0)*ldcrd(1,1);

                  /* OTHER SYMMETRIES    */            
#define           viscI1II0II0II0I visc[0][1][0][0]
#define           viscI1II0II1II1I visc[0][1][1][1]
#define           viscI1II0II0II1I visc[0][1][1][0]
#define           viscI1II0II1II0I visc[0][1][0][1]

                  /* TIME DERIVATIVE TERMS */ 
                  for(i=0;i<lgpx;++i) {
                     for(j=0;j<lgpn;++j) {
                        rhorbd0 = RAD(crd(0)(i,j))*lrhorbd0;
                        
                        /* UNSTEADY TERMS */
                        res(0)(i,j) = rhorbd0*u(0)(i,j) +dugdt(log2p,tind,0)(i,j);
                        res(1)(i,j) = rhorbd0*u(1)(i,j) +dugdt(log2p,tind,1)(i,j);
                        res(2)(i,j) = rhorbd0 +dugdt(log2p,tind,2)(i,j);
                        
#ifdef AXISYMMETRIC
                        res(0)(i,j) -= cjcb*(u(2)(i,j) -2.*lmu*u(0)(i,j)/crd(0)(i,j));
#endif
#ifdef SOURCETERM
                        res(0)(i,j) -= RAD(crd(0)(i,j))*cjcb*amp*pow(crd(0)(i,j)-0.5,amp-1.);
#endif

#ifdef BODY
                        res(0)(i,j) -= ins_gbl->rho*RAD(crd(0)(i,j))*cjcb*body[0];
                        res(1)(i,j) -= ins_gbl->rho*RAD(crd(0)(i,j))*cjcb*body[1];
#ifdef INERTIALESS
                        res(0)(i,j) = -ins_gbl->rho*RAD(crd(0)(i,j))*cjcb*body[0];
                        res(1)(i,j) = -ins_gbl->rho*RAD(crd(0)(i,j))*cjcb*body[1];
#endif
#endif
                        e00[i][j] = RAD(crd(0)(i,j))*(+visc[0][0][0][0]*du(0,0)(i,j) +visc[0][1][0][0]*du(1,0)(i,j)
                                              +visc[0][0][0][1]*du(0,1)(i,j) +visc[0][1][0][1]*du(1,1)(i,j));

                        e01[i][j] = RAD(crd(0)(i,j))*(+viscI0II0II1II0I*du(0,0)(i,j) +visc[0][1][1][0]*du(1,0)(i,j)
                                              +visc[0][0][1][1]*du(0,1)(i,j) +visc[0][1][1][1]*du(1,1)(i,j));

                        e10[i][j] = RAD(crd(0)(i,j))*(+viscI1II0II0II0I*du(0,0)(i,j) +visc[1][1][0][0]*du(1,0)(i,j)
                                              +viscI1II0II0II1I*du(0,1)(i,j) +visc[1][1][0][1]*du(1,1)(i,j));

                        e11[i][j] = RAD(crd(0)(i,j))*(+viscI1II0II1II0I*du(0,0)(i,j) +viscI1II1II1II0I*du(1,0)(i,j)
                                              +viscI1II0II1II1I*du(0,1)(i,j) +visc[1][1][1][1]*du(1,1)(i,j));
                        
                        cv00[i][j] += e00[i][j];
                        cv01[i][j] += e01[i][j];
                        cv10[i][j] += e10[i][j];
                        cv11[i][j] += e11[i][j];
                      }
                  }

                  
                  basis::tri(log2p).intgrt(&lf(0)(0),&res(0)(0,0),MXGP);
                  basis::tri(log2p).intgrt(&lf(1)(0),&res(1)(0,0),MXGP);
                  basis::tri(log2p).intgrt(&lf(2)(0),&res(2)(0,0),MXGP);

                  /* CALCULATE RESIDUAL TO GOVERNING EQUATION & STORE IN RES */
                  basis::tri(log2p).derivr(cv00[0],&res(0)(0,0),MXGP);
                  basis::tri(log2p).derivs(cv01[0],&res(0)(0,0),MXGP);
                  basis::tri(log2p).derivr(cv10[0],&res(1)(0,0),MXGP);
                  basis::tri(log2p).derivs(cv11[0],&res(1)(0,0),MXGP);
                  basis::tri(log2p).derivr(&du(2,0)(0,0),&res(2)(0,0),MXGP);
                  basis::tri(log2p).derivs(&du(2,1)(0,0),&res(2)(0,0),MXGP);
                  
                  /* THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES */
                  for(i=0;i<lgpx;++i) {
                     for(j=0;j<lgpn;++j) {
                        tres[0] = ins_gbl->tau(tind)*res(0)(i,j);
                        tres[1] = ins_gbl->tau(tind)*res(1)(i,j);
                        tres[2] = ins_gbl->delt(tind)*res(2)(i,j);

#ifndef INERTIALESS
                        e00[i][j] -= (ldcrd(1,1)*(2*u(0)(i,j)-mvel(0)(i,j))
                                    -ldcrd(0,1)*(u(1)(i,j)-mvel(1)(i,j)))*tres[0]
                                    -ldcrd(0,1)*u(0)(i,j)*tres[1]
                                    +ldcrd(1,1)*tres[2];
                        e01[i][j] -= (-ldcrd(1,0)*(2*u(0)(i,j)-mvel(0)(i,j))
                                    +ldcrd(0,0)*(u(1)(i,j)-mvel(1)(i,j)))*tres[0]
                                    +ldcrd(0,0)*u(0)(i,j)*tres[1]
                                    -ldcrd(1,0)*tres[2];
                        e10[i][j] -= +ldcrd(1,1)*u(1)(i,j)*tres[0]
                                    +(ldcrd(1,1)*(u(0)(i,j)-mvel(0)(i,j))
                                    -ldcrd(0,1)*(2.*u(1)(i,j)-mvel(1)(i,j)))*tres[1]
                                    -ldcrd(0,1)*tres[2];
                        e11[i][j] -= -ldcrd(1,0)*u(1)(i,j)*tres[0]
                                    +(-ldcrd(1,0)*(u(0)(i,j)-mvel(0)(i,j))
                                    +ldcrd(0,0)*(2.*u(1)(i,j)-mvel(1)(i,j)))*tres[1]
                                    +ldcrd(0,0)*tres[2];
#endif
                                    
                        du(2,0)(i,j) = -(ldcrd(1,1)*tres[0] -ldcrd(0,1)*tres[1]);
                        du(2,1)(i,j) = -(-ldcrd(1,0)*tres[0] +ldcrd(0,0)*tres[1]);
                     }
                  }
                  basis::tri(log2p).intgrtrs(&lf(0)(0),e00[0],e01[0],MXGP);
                  basis::tri(log2p).intgrtrs(&lf(1)(0),e10[0],e11[0],MXGP);
                  basis::tri(log2p).intgrtrs(&lf(2)(0),&du(2,0)(0,0),&du(2,1)(0,0),MXGP);
                  
                  for(n=0;n<NV;++n)
                     for(i=0;i<basis::tri(log2p).tm;++i)
                        lf(n)(i) *= sim::beta[stage];
                        
                  lftog(tind,ins_gbl->res_r);
               }
            }
         }

         /* ADD IN VISCOUS/DISSIPATIVE FLUX */
         hp_gbl->res.v(Range(0,nvrtx-1),Range::all()) += hp_gbl->res_r.v(Range(0,nvrtx-1),Range::all());
         if (basis::tri(log2p).sm) {
            hp_gbl->res.s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all()) += hp_gbl->res_r.s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all());        
            if (basis::tri(log2p).im) {
               hp_gbl->res.i(Range(0,ntri-1),Range(0,basis::tri(log2p).im-1),Range::all()) += hp_gbl->res_r.i(Range(0,ntri-1),Range(0,basis::tri(log2p).im-1),Range::all());     
            }
         }

      /*********************************************/
         /* MODIFY RESIDUALS ON COARSER MESHES         */
      /*********************************************/   
         if(coarse) {
         /* CALCULATE DRIVING TERM ON FIRST ENTRY TO COARSE MESH */
            if(isfrst) {
               dres(log2p).v(Range(0,nvrtx-1),Range::all()) = fadd*hp_gbl->res0.v(Range(0,nvrtx-1),Range::all()) -hp_gbl->res.v(Range(0,nvrtx-1),Range::all());
               if (basis::tri(log2p).sm) dres(log2p).s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all()) = fadd*hp_gbl->res0.s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all()) -hp_gbl->res.s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all());     
               if (basis::tri(log2p).im) dres(log2p).i(Range(0,ntri-1),Range(0,basis::tri(log2p).im-1),Range::all()) = fadd*hp_gbl->res0.i(Range(0,ntri-1),Range(0,basis::tri(log2p).im-1),Range::all()) -hp_gbl->res.i(Range(0,ntri-1),Range(0,basis::tri(log2p).im-1),Range::all());
               isfrst = false;
            }
            hp_gbl->res.v(Range(0,nvrtx-1),Range::all()) += dres(log2p).v(Range(0,nvrtx-1),Range::all()); 
            if (basis::tri(log2p).sm) hp_gbl->res.s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all()) += dres(log2p).s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all());
            if (basis::tri(log2p).im) hp_gbl->res.i(Range(0,ntri-1),Range(0,basis::tri(log2p).im-1),Range::all()) += dres(log2p).i(Range(0,ntri-1),Range(0,basis::tri(log2p).im-1),Range::all());  
            

         }
         
         ++excpt;
      }
   }

   return(block::stop);
}
