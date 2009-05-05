/*
 *  rsdl.cpp
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_ins.h"
#include "../hp_boundary.h"

// #define BODYFORCE
    
void tri_hp_ins::rsdl(int stage) {
    int i,j,n,tind;
    FLT fluxx,fluxy;
    const int NV = 3;
    TinyVector<int,3> v;
    TinyMatrix<FLT,ND,ND> ldcrd;
    TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,NV,ND> du;
    int lgpx = basis::tri(log2p).gpx, lgpn = basis::tri(log2p).gpn;
    FLT rhobd0 = gbl->rho*gbl->bd(0), lmu = gbl->mu, rhorbd0, lrhorbd0, cjcb, cjcbi, oneminusbeta;
    TinyMatrix<TinyMatrix<FLT,ND,ND>,NV-1,NV-1> visc;
    TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,NV-1,NV-1> cv, df;
    TinyVector<FLT,NV> tres;
    
    tri_hp::rsdl(stage);
    
    oneminusbeta = 1.0-gbl->beta(stage);
    
    for(tind = 0; tind<ntri;++tind) {
        /* LOAD INDICES OF VERTEX POINTS */
        v = tri(tind).pnt;
    
        /* IF TINFO > -1 IT IS CURVED ELEMENT */
        if (tri(tind).info > -1) {
            /* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
            crdtocht(tind);

            /* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
            for(n=0;n<ND;++n)
                basis::tri(log2p).proj_bdry(&cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), &dcrd(n,1)(0,0),MXGP);
        }
        else {
            /* PROJECT VERTEX COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
            for(n=0;n<ND;++n)
                basis::tri(log2p).proj(pnts(v(0))(n),pnts(v(1))(n),pnts(v(2))(n),&crd(n)(0,0),MXGP);

            /* CALCULATE COORDINATE DERIVATIVES A SIMPLE WAY */
            for(n=0;n<ND;++n) {
                ldcrd(n,0) = 0.5*(pnts(v(2))(n) -pnts(v(1))(n));
                ldcrd(n,1) = 0.5*(pnts(v(0))(n) -pnts(v(1))(n));
            }
        }

        /* CALCULATE MESH VELOCITY */
        for(i=0;i<lgpx;++i) {
            for(j=0;j<lgpn;++j) {
                mvel(0)(i,j) = gbl->bd(0)*(crd(0)(i,j) -dxdt(log2p,tind,0)(i,j));
                mvel(1)(i,j) = gbl->bd(0)*(crd(1)(i,j) -dxdt(log2p,tind,1)(i,j));
#ifdef DROP
                mvel(0)(i,j) += mesh_ref_vel(0);
                mvel(1)(i,j) += mesh_ref_vel(1);
#endif                        
            }
        }

        /* LOAD SOLUTION COEFFICIENTS FOR THIS ELEMENT */
        /* PROJECT SOLUTION TO GAUSS POINTS WITH DERIVATIVES IF NEEDED FOR VISCOUS TERMS */
        ugtouht(tind);
        if (gbl->beta(stage) > 0.0) {
            for(n=0;n<NV-1;++n)
                basis::tri(log2p).proj(&uht(n)(0),&u(n)(0,0),&du(n,0)(0,0),&du(n,1)(0,0),MXGP);
            basis::tri(log2p).proj(&uht(NV-1)(0),&u(NV-1)(0,0),MXGP);
        }
        else {
            for(n=0;n<NV;++n)
                basis::tri(log2p).proj(&uht(n)(0),&u(n)(0,0),MXGP);
        }
        
        /* lf IS WHERE I WILL STORE THE ELEMENT RESIDUAL */
        for(n=0;n<NV;++n)
            for(i=0;i<basis::tri(log2p).tm;++i)
                lf(n)(i) = 0.0;

        if (tri(tind).info > -1) {
            /* CURVED ELEMENT */
            /* CONVECTIVE TERMS (IMAGINARY FIRST)*/
            for(i=0;i<lgpx;++i) {
                for(j=0;j<lgpn;++j) {

                    fluxx = gbl->rho*RAD(crd(0)(i,j))*(u(0)(i,j) -mvel(0)(i,j));
                    fluxy = gbl->rho*RAD(crd(0)(i,j))*(u(1)(i,j) -mvel(1)(i,j));
                    
                    /* CONTINUITY EQUATION FLUXES */
                    du(NV-1,0)(i,j) = +dcrd(1,1)(i,j)*fluxx -dcrd(0,1)(i,j)*fluxy;
                    du(NV-1,1)(i,j) = -dcrd(1,0)(i,j)*fluxx +dcrd(0,0)(i,j)*fluxy;
                  
#ifndef INERTIALESS
                    /* CONVECTIVE FLUXES */
                    for(n=0;n<NV-1;++n) {
                        cv(n,0)(i,j) = u(n)(i,j)*du(NV-1,0)(i,j);
                        cv(n,1)(i,j) = u(n)(i,j)*du(NV-1,1)(i,j);
                    }
#else
                    for(n=0;n<NV-1;++n) {
                        cv(n,0)(i,j) = 0.0;
                        cv(n,1)(i,j) = 0.0;
                    }
#endif
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
            lftog(tind,gbl->res);

            /* NEGATIVE REAL TERMS */
            if (gbl->beta(stage) > 0.0) {
                /* TIME DERIVATIVE TERMS */ 
                for(i=0;i<lgpx;++i) {
                    for(j=0;j<lgpn;++j) {
                        cjcb = dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j);
                        rhorbd0 = rhobd0*RAD(crd(0)(i,j))*cjcb;
                        cjcbi = lmu*RAD(crd(0)(i,j))/cjcb;
                        
                        /* UNSTEADY TERMS */
                        for(n=0;n<NV-1;++n)
                            res(n)(i,j) = rhorbd0*u(n)(i,j) +dugdt(log2p,tind,n)(i,j);
                        res(NV-1)(i,j) = rhorbd0 +dugdt(log2p,tind,NV-1)(i,j);
						
#ifdef INERTIALESS
                        res(0)(i,j) = 0.0;
                        res(1)(i,j) = 0.0;
#endif
#ifdef AXISYMMETRIC
                        res(0)(i,j) -= cjcb*(u(NV-1)(i,j) -2.*lmu*u(0)(i,j)/crd(0)(i,j));
#endif
#ifdef BODYFORCE
                        res(0)(i,j) -= gbl->rho*RAD(crd(0)(i,j))*cjcb*gbl->body(0);
                        res(1)(i,j) -= gbl->rho*RAD(crd(0)(i,j))*cjcb*gbl->body(1);
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

                        /* OTHER SYMMETRIES     */                
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
                        
                        for(n=0;n<NV-1;++n) {
                            cv(n,0)(i,j) += df(n,0)(i,j);
                            cv(n,1)(i,j) += df(n,1)(i,j);
                        }
                     }
                }
                for(n=0;n<NV;++n)
                    basis::tri(log2p).intgrt(&lf(n)(0),&res(n)(0,0),MXGP);

                /* CALCULATE RESIDUAL TO GOVERNING EQUATION & STORE IN RES */
                for(n=0;n<NV-1;++n) {
                    basis::tri(log2p).derivr(&cv(n,0)(0,0),&res(n)(0,0),MXGP);
                    basis::tri(log2p).derivs(&cv(n,1)(0,0),&res(n)(0,0),MXGP);
                }
                basis::tri(log2p).derivr(&du(NV-1,0)(0,0),&res(NV-1)(0,0),MXGP);
                basis::tri(log2p).derivs(&du(NV-1,1)(0,0),&res(NV-1)(0,0),MXGP);

                /* THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES */
                for(i=0;i<lgpx;++i) {
                    for(j=0;j<lgpn;++j) {
                                    
                        tres(0) = gbl->tau(tind,0)*res(0)(i,j);
                        tres(1) = gbl->tau(tind,0)*res(1)(i,j);
                        tres(NV-1) = gbl->tau(tind,NV-1)*res(NV-1)(i,j);

#ifndef INERTIALESS
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
#endif
                                        
                        du(NV-1,0)(i,j) = -(dcrd(1,1)(i,j)*tres(0) -dcrd(0,1)(i,j)*tres(1));
                        du(NV-1,1)(i,j) = -(-dcrd(1,0)(i,j)*tres(0) +dcrd(0,0)(i,j)*tres(1));
                    }
                }
                for(n=0;n<NV-1;++n)
                    basis::tri(log2p).intgrtrs(&lf(n)(0),&df(n,0)(0,0),&df(n,1)(0,0),MXGP);
                basis::tri(log2p).intgrtrs(&lf(NV-1)(0),&du(NV-1,0)(0,0),&du(NV-1,1)(0,0),MXGP);
              
                for(n=0;n<NV;++n)
                    for(i=0;i<basis::tri(log2p).tm;++i)
                        lf(n)(i) *= gbl->beta(stage);
                        
                lftog(tind,gbl->res_r);
            }
        }
        else {
            /* LINEAR ELEMENT */
            /* CONVECTIVE TERMS (IMAGINARY FIRST)*/
            for(i=0;i<lgpx;++i) {
                for(j=0;j<lgpn;++j) {

                    fluxx = gbl->rho*RAD(crd(0)(i,j))*(u(0)(i,j) -mvel(0)(i,j));
                    fluxy = gbl->rho*RAD(crd(0)(i,j))*(u(1)(i,j) -mvel(1)(i,j));
                    
                    /* CONTINUITY EQUATION FLUXES */
                    du(NV-1,0)(i,j) = +ldcrd(1,1)*fluxx -ldcrd(0,1)*fluxy;
                    du(NV-1,1)(i,j) = -ldcrd(1,0)*fluxx +ldcrd(0,0)*fluxy;
                    
#ifndef INERTIALESS
                    /* CONVECTIVE FLUXES */
                    for(n=0;n<NV-1;++n) {
                        cv(n,0)(i,j) = u(n)(i,j)*du(NV-1,0)(i,j);
                        cv(n,1)(i,j) = u(n)(i,j)*du(NV-1,1)(i,j);
                    }
#else
                    for(n=0;n<NV-1;++n) {
                        cv(n,0)(i,j) = 0.0;
                        cv(n,1)(i,j) = 0.0;
                    }
#endif
                    /* PRESSURE TERMS */
                    /* U-MOMENTUM */
                    cv(0,0)(i,j) += ldcrd(1,1)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
                    cv(0,1)(i,j) -= ldcrd(1,0)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
                    /* V-MOMENTUM */
                    cv(1,0)(i,j) -=  ldcrd(0,1)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
                    cv(1,1)(i,j) +=  ldcrd(0,0)*RAD(crd(0)(i,j))*u(NV-1)(i,j);
                }
            }
            for(n=0;n<NV-1;++n)
                basis::tri(log2p).intgrtrs(&lf(n)(0),&cv(n,0)(0,0),&cv(n,1)(0,0),MXGP);
            basis::tri(log2p).intgrtrs(&lf(NV-1)(0),&du(NV-1,0)(0,0),&du(NV-1,1)(0,0),MXGP);
            

            /* ASSEMBLE GLOBAL FORCING (IMAGINARY TERMS) */
            lftog(tind,gbl->res);

            /* NEGATIVE REAL TERMS */
            if (gbl->beta(stage) > 0.0) {
                cjcb = ldcrd(0,0)*ldcrd(1,1) -ldcrd(1,0)*ldcrd(0,1);
                cjcbi = lmu/cjcb;
                lrhorbd0 = rhobd0*cjcb;
                
                /* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
                /* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
                visc(0,0)(0,0) = -cjcbi*(2.*ldcrd(1,1)*ldcrd(1,1) +ldcrd(0,1)*ldcrd(0,1));
                visc(0,0)(1,1) = -cjcbi*(2.*ldcrd(1,0)*ldcrd(1,0) +ldcrd(0,0)*ldcrd(0,0));
                visc(0,0)(0,1) =  cjcbi*(2.*ldcrd(1,1)*ldcrd(1,0) +ldcrd(0,1)*ldcrd(0,0));
#define         viscI0II0II1II0I visc(0,0)(0,1)

                visc(1,1)(0,0) = -cjcbi*(ldcrd(1,1)*ldcrd(1,1) +2.*ldcrd(0,1)*ldcrd(0,1));
                visc(1,1)(1,1) = -cjcbi*(ldcrd(1,0)*ldcrd(1,0) +2.*ldcrd(0,0)*ldcrd(0,0));
                visc(1,1)(0,1) =  cjcbi*(ldcrd(1,1)*ldcrd(1,0) +2.*ldcrd(0,1)*ldcrd(0,0));
#define         viscI1II1II1II0I visc(1,1)(0,1)
                
                visc(0,1)(0,0) =  cjcbi*ldcrd(0,1)*ldcrd(1,1);
                visc(0,1)(1,1) =  cjcbi*ldcrd(0,0)*ldcrd(1,0);
                visc(0,1)(0,1) = -cjcbi*ldcrd(0,1)*ldcrd(1,0);
                visc(0,1)(1,0) = -cjcbi*ldcrd(0,0)*ldcrd(1,1);

                /* OTHER SYMMETRIES     */                
#define         viscI1II0II0II0I visc(0,1)(0,0)
#define         viscI1II0II1II1I visc(0,1)(1,1)
#define         viscI1II0II0II1I visc(0,1)(1,0)
#define         viscI1II0II1II0I visc(0,1)(0,1)

                /* TIME DERIVATIVE TERMS */ 
                for(i=0;i<lgpx;++i) {
                    for(j=0;j<lgpn;++j) {
                        rhorbd0 = RAD(crd(0)(i,j))*lrhorbd0;
                        
                        /* UNSTEADY TERMS */
                        for(n=0;n<NV-1;++n)
                            res(n)(i,j) = rhorbd0*u(n)(i,j) +dugdt(log2p,tind,n)(i,j);
                        res(NV-1)(i,j) = rhorbd0 +dugdt(log2p,tind,NV-1)(i,j);
						
#ifdef INERTIALESS
                        res(0)(i,j) = 0.0;
                        res(1)(i,j) = 0.0;
#endif
#ifdef AXISYMMETRIC
                        res(0)(i,j) -= cjcb*(u(2)(i,j) -2.*lmu*u(0)(i,j)/crd(0)(i,j));
#endif
#ifdef BODYFORCE
                        res(0)(i,j) -= gbl->rho*RAD(crd(0)(i,j))*cjcb*gbl->body(0);
                        res(1)(i,j) -= gbl->rho*RAD(crd(0)(i,j))*cjcb*gbl->body(1);
#endif        
                        df(0,0)(i,j) = RAD(crd(0)(i,j))*(+visc(0,0)(0,0)*du(0,0)(i,j) +visc(0,1)(0,0)*du(1,0)(i,j)
                                                     +visc(0,0)(0,1)*du(0,1)(i,j) +visc(0,1)(0,1)*du(1,1)(i,j));

                        df(0,1)(i,j) = RAD(crd(0)(i,j))*(+viscI0II0II1II0I*du(0,0)(i,j) +visc(0,1)(1,0)*du(1,0)(i,j)
                                                     +visc(0,0)(1,1)*du(0,1)(i,j) +visc(0,1)(1,1)*du(1,1)(i,j));

                        df(1,0)(i,j) = RAD(crd(0)(i,j))*(+viscI1II0II0II0I*du(0,0)(i,j) +visc(1,1)(0,0)*du(1,0)(i,j)
                                                     +viscI1II0II0II1I*du(0,1)(i,j) +visc(1,1)(0,1)*du(1,1)(i,j));

                        df(1,1)(i,j) = RAD(crd(0)(i,j))*(+viscI1II0II1II0I*du(0,0)(i,j) +viscI1II1II1II0I*du(1,0)(i,j)
                                                     +viscI1II0II1II1I*du(0,1)(i,j) +visc(1,1)(1,1)*du(1,1)(i,j));
                        
                        for(n=0;n<NV-1;++n) {
                            cv(n,0)(i,j) += df(n,0)(i,j);
                            cv(n,1)(i,j) += df(n,1)(i,j);
                        } 
                     }
                }
                for(n=0;n<NV;++n)
                    basis::tri(log2p).intgrt(&lf(n)(0),&res(n)(0,0),MXGP);

                /* CALCULATE RESIDUAL TO GOVERNING EQUATION & STORE IN RES */
                for(n=0;n<NV-1;++n) {
                    basis::tri(log2p).derivr(&cv(n,0)(0,0),&res(n)(0,0),MXGP);
                    basis::tri(log2p).derivs(&cv(n,1)(0,0),&res(n)(0,0),MXGP);
                }
                basis::tri(log2p).derivr(&du(NV-1,0)(0,0),&res(NV-1)(0,0),MXGP);
                basis::tri(log2p).derivs(&du(NV-1,1)(0,0),&res(NV-1)(0,0),MXGP);
                
                /* THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES */
                for(i=0;i<lgpx;++i) {
                    for(j=0;j<lgpn;++j) {
                        tres(0) = gbl->tau(tind,0)*res(0)(i,j);
                        tres(1) = gbl->tau(tind,0)*res(1)(i,j);
                        tres(NV-1) = gbl->tau(tind,NV-1)*res(NV-1)(i,j);

#ifndef INERTIALESS
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
#endif
                                        
                        du(NV-1,0)(i,j) = -(ldcrd(1,1)*tres(0) -ldcrd(0,1)*tres(1));
                        du(NV-1,1)(i,j) = -(-ldcrd(1,0)*tres(0) +ldcrd(0,0)*tres(1));
                    }
                }
                for(n=0;n<NV-1;++n)
                    basis::tri(log2p).intgrtrs(&lf(n)(0),&df(n,0)(0,0),&df(n,1)(0,0),MXGP);
                basis::tri(log2p).intgrtrs(&lf(NV-1)(0),&du(NV-1,0)(0,0),&du(NV-1,1)(0,0),MXGP);
                
                for(n=0;n<NV;++n)
                    for(i=0;i<basis::tri(log2p).tm;++i)
                        lf(n)(i) *= gbl->beta(stage);
                        
                lftog(tind,gbl->res_r);
            }
        }
    }

    /* ADD IN VISCOUS/DISSIPATIVE FLUX */
    gbl->res.v(Range(0,npnt-1),Range::all()) += gbl->res_r.v(Range(0,npnt-1),Range::all());
    if (basis::tri(log2p).sm) {
        gbl->res.s(Range(0,nseg-1),Range(0,basis::tri(log2p).sm-1),Range::all()) += gbl->res_r.s(Range(0,nseg-1),Range(0,basis::tri(log2p).sm-1),Range::all());          
        if (basis::tri(log2p).im) {
            gbl->res.i(Range(0,ntri-1),Range(0,basis::tri(log2p).im-1),Range::all()) += gbl->res_r.i(Range(0,ntri-1),Range(0,basis::tri(log2p).im-1),Range::all());      
        }
    }
    
    
#ifdef DEBUG
    if (coarse_flag) {
    for(i=0;i<npnt;++i) {
        printf("rsdl v: %d ",i);
        for (n=0;n<NV;++n) 
            printf("%e ",gbl->res.v(i,n));
        printf("\n");
    }
        
    for(i=0;i<nseg;++i) {
        for(int m=0;m<basis::tri(log2p).sm;++m) {
            printf("rsdl s: %d %d ",i,m); 
            for(n=0;n<NV;++n)
                printf("%e ",gbl->res.s(i,m,n));
            printf("\n");
        }
    }

    for(i=0;i<ntri;++i) {
        for(int m=0;m<basis::tri(log2p).im;++m) {
            printf("rsdl i: %d %d ",i,m);
            for(n=0;n<NV;++n) 
                printf("%e %e %e\n",gbl->res.i(i,m,n));
            printf("\n");
        }
    }
    }
#endif

    /*********************************************/
    /* MODIFY RESIDUALS ON COARSER MESHES            */
    /*********************************************/    
    if(coarse_flag) {
    /* CALCULATE DRIVING TERM ON FIRST ENTRY TO COARSE MESH */
        if(isfrst) {
            dres(log2p).v(Range(0,npnt-1),Range::all()) = fadd*gbl->res0.v(Range(0,npnt-1),Range::all()) -gbl->res.v(Range(0,npnt-1),Range::all());
            if (basis::tri(log2p).sm) dres(log2p).s(Range(0,nseg-1),Range(0,basis::tri(log2p).sm-1),Range::all()) = fadd*gbl->res0.s(Range(0,nseg-1),Range(0,basis::tri(log2p).sm-1),Range::all()) -gbl->res.s(Range(0,nseg-1),Range(0,basis::tri(log2p).sm-1),Range::all());      
            if (basis::tri(log2p).im) dres(log2p).i(Range(0,ntri-1),Range(0,basis::tri(log2p).im-1),Range::all()) = fadd*gbl->res0.i(Range(0,ntri-1),Range(0,basis::tri(log2p).im-1),Range::all()) -gbl->res.i(Range(0,ntri-1),Range(0,basis::tri(log2p).im-1),Range::all());
            isfrst = false;
        }
        gbl->res.v(Range(0,npnt-1),Range::all()) += dres(log2p).v(Range(0,npnt-1),Range::all()); 
        if (basis::tri(log2p).sm) gbl->res.s(Range(0,nseg-1),Range(0,basis::tri(log2p).sm-1),Range::all()) += dres(log2p).s(Range(0,nseg-1),Range(0,basis::tri(log2p).sm-1),Range::all());
        if (basis::tri(log2p).im) gbl->res.i(Range(0,ntri-1),Range(0,basis::tri(log2p).im-1),Range::all()) += dres(log2p).i(Range(0,ntri-1),Range(0,basis::tri(log2p).im-1),Range::all());  
    }
    else {
        if (stage == gbl->nstage) {
            /* HACK FOR AUXILIARY FLUXES */
            for (i=0;i<nebd;++i)
                hp_ebdry(i)->output(*gbl->log, tri_hp::auxiliary);
        }
    }
        

            
    return;
}
