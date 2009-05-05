/*
 *  rsdl.cpp
 *  swirl++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_swirl.h"
#include "../hp_boundary.h"
    
void tri_hp_swirl::rsdl(int stage) {
    int i,j,n,tind;
    FLT fluxx,fluxy;
    const int NV = 4;
    TinyVector<int,3> v;
    TinyMatrix<FLT,ND,ND> ldcrd;
    TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,NV,ND> du;
    int lgpx = basis::tri(log2p).gpx, lgpn = basis::tri(log2p).gpn;
    FLT rhobd0 = gbl->rho*gbl->bd(0), lmu = gbl->mu, rhorbd0, cjcb, cjcbi, oneminusbeta;
    FLT visc[ND+1][ND+1][ND][ND], tres[NV];
    FLT cv00[MXGP][MXGP],cv01[MXGP][MXGP],cv10[MXGP][MXGP],cv11[MXGP][MXGP],cv20[MXGP][MXGP],cv21[MXGP][MXGP]; // LOCAL WORK ARRAYS
    FLT e00[MXGP][MXGP],e01[MXGP][MXGP],e10[MXGP][MXGP],e11[MXGP][MXGP],e20[MXGP][MXGP],e21[MXGP][MXGP]; // LOCAL WORK ARRAYS

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
             }
        }

        /* LOAD SOLUTION COEFFICIENTS FOR THIS ELEMENT */
        /* PROJECT SOLUTION TO GAUSS POINTS WITH DERIVATIVES IF NEEDED FOR VISCOUS TERMS */
        ugtouht(tind);
        if (gbl->beta(stage) > 0.0) {
            basis::tri(log2p).proj(&uht(0)(0),&u(0)(0,0),&du(0,0)(0,0),&du(0,1)(0,0),MXGP);
            basis::tri(log2p).proj(&uht(1)(0),&u(1)(0,0),&du(1,0)(0,0),&du(1,1)(0,0),MXGP);
            basis::tri(log2p).proj(&uht(2)(0),&u(2)(0,0),&du(2,0)(0,0),&du(2,1)(0,0),MXGP);
            basis::tri(log2p).proj(&uht(3)(0),&u(3)(0,0),MXGP);
        }
        else {
            basis::tri(log2p).proj(&uht(0)(0),&u(0)(0,0),MXGP);
            basis::tri(log2p).proj(&uht(1)(0),&u(1)(0,0),MXGP);
            basis::tri(log2p).proj(&uht(2)(0),&u(2)(0,0),MXGP);
            basis::tri(log2p).proj(&uht(3)(0),&u(3)(0,0),MXGP);
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
                    du(3,0)(i,j) = +dcrd(1,1)(i,j)*fluxx -dcrd(0,1)(i,j)*fluxy;
                    du(3,1)(i,j) = -dcrd(1,0)(i,j)*fluxx +dcrd(0,0)(i,j)*fluxy;

                    /* U-MOMENTUM */
                    cv00[i][j] =  u(0)(i,j)*du(3,0)(i,j) +dcrd(1,1)(i,j)*RAD(crd(0)(i,j))*u(3)(i,j);
                    cv01[i][j] =  u(0)(i,j)*du(3,1)(i,j) -dcrd(1,0)(i,j)*RAD(crd(0)(i,j))*u(3)(i,j);
                    /* V-MOMENTUM */
                    cv10[i][j] =  u(1)(i,j)*du(3,0)(i,j) -dcrd(0,1)(i,j)*RAD(crd(0)(i,j))*u(3)(i,j);
                    cv11[i][j] =  u(1)(i,j)*du(3,1)(i,j) +dcrd(0,0)(i,j)*RAD(crd(0)(i,j))*u(3)(i,j);
                    /* W-MOMENTUM */
                    cv20[i][j] =  u(2)(i,j)*du(3,0)(i,j);
                    cv21[i][j] =  u(2)(i,j)*du(3,1)(i,j);
                }
            }
            basis::tri(log2p).intgrtrs(&lf(0)(0),cv00[0],cv01[0],MXGP);
            basis::tri(log2p).intgrtrs(&lf(1)(0),cv10[0],cv11[0],MXGP);
            basis::tri(log2p).intgrtrs(&lf(2)(0),cv20[0],cv21[0],MXGP);
            basis::tri(log2p).intgrtrs(&lf(3)(0),&du(3,0)(0,0),&du(3,1)(0,0),MXGP);

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
                        res(0)(i,j) = rhorbd0*u(0)(i,j) +dugdt(log2p,tind,0)(i,j);
                        res(1)(i,j) = rhorbd0*u(1)(i,j) +dugdt(log2p,tind,1)(i,j);
                        res(2)(i,j) = rhorbd0*u(2)(i,j) +dugdt(log2p,tind,2)(i,j);
                        res(3)(i,j) = rhorbd0 +dugdt(log2p,tind,3)(i,j);
                        
#ifdef AXISYMMETRIC
						/* SOURCE TERMS */
                        res(0)(i,j) -= cjcb*(u(3)(i,j) -2.*lmu*u(0)(i,j)/crd(0)(i,j) +gbl->rho*u(2)(i,j)*u(2)(i,j));
                        //res(2)(i,j) -= cjcb*(lmu*((dcrd(1,1)(i,j)*du(2,1)(i,j) -dcrd(1,0)(i,j)*du(2,1)(i,j))/cjcb -u(2)(i,j)/crd(0)(i,j)) -gbl->rho*u(0)(i,j)*u(2)(i,j));	
                        //res(2)(i,j) += lmu*(dcrd(1,1)(i,j)*du(2,0)(i,j)-dcrd(1,0)(i,j)*du(2,1)(i,j))-lmu*cjcb*u(2)(i,j)/crd(0)(i,j) +cjcb*gbl->rho*u(0)(i,j)*u(2)(i,j);
                        //res(2)(i,j) += cjcb*gbl->rho*u(0)(i,j)*u(2)(i,j) -cjcb*lmu*((dcrd(1,1)(i,j)*du(2,0)(i,j)-dcrd(1,0)(i,j)*du(2,1)(i,j)/cjcb) -u(2)(i,j)/crd(0)(i,j));
                        res(2)(i,j) += cjcb*gbl->rho*u(0)(i,j)*u(2)(i,j) -lmu*dcrd(1,1)(i,j)*du(2,0)(i,j) +lmu*dcrd(1,0)(i,j)*du(2,1)(i,j) +cjcb*lmu*u(2)(i,j)/crd(0)(i,j);
#else
						res(2)(i,j) -= cjcb*RAD(crd(0)(i,j))*dpdz;
#endif
                        
                        /* Source Terms for Test Function */
                        //res(0)(i,j) -= crd(0)(i,j)*cjcb*(gbl->rho*pow(crd(0)(i,j),3.0) +2*crd(0)(i,j) -3*lmu);
                        
                        //res(1)(i,j) -= crd(0)(i,j)*cjcb*(6*gbl->rho*pow(crd(0)(i,j),2.0)*crd(1)(i,j) +3/crd(0)(i,j)*lmu*crd(1)(i,j));
                        
                        //res(2)(i,j) -= crd(0)(i,j)*cjcb*(3*gbl->rho*pow(crd(0)(i,j),3.0) -3*lmu);
                        

                        /* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH) */
                        /* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S) */
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
                  
                        visc[2][2][0][0] = -cjcbi*(dcrd(0,1)(i,j)*dcrd(0,1)(i,j) +dcrd(1,1)(i,j)*dcrd(1,1)(i,j));
                        visc[2][2][0][1] =  cjcbi*(dcrd(0,0)(i,j)*dcrd(0,1)(i,j) +dcrd(1,0)(i,j)*dcrd(1,1)(i,j));
                        visc[2][2][1][1] = -cjcbi*(dcrd(0,0)(i,j)*dcrd(0,0)(i,j) +dcrd(1,0)(i,j)*dcrd(1,0)(i,j));
#define					viscI2II2II1II0I visc[2][2][0][1]

                        /* OTHER SYMMETRIES     */                
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
                  
                        e20[i][j] = +visc[2][2][0][0]*du(2,0)(i,j) +visc[2][2][0][1]*du(2,1)(i,j);
#ifdef AXISYMMETRIC
						e20[i][j] += lmu*dcrd(1,1)(i,j)*u(2)(i,j);
#endif
                        e21[i][j] = +viscI2II2II1II0I*du(2,0)(i,j) +visc[2][2][1][1]*du(2,1)(i,j);
#ifdef AXISYMMETRIC
						e21[i][j] -= lmu*dcrd(1,0)(i,j)*u(2)(i,j);
#endif
                
                                        
                        cv00[i][j] += e00[i][j];
                        cv01[i][j] += e01[i][j];
                        cv10[i][j] += e10[i][j];
                        cv11[i][j] += e11[i][j];
                        cv20[i][j] += e20[i][j];
                        cv21[i][j] += e21[i][j];
                     }
                }
                basis::tri(log2p).intgrt(&lf(0)(0),&res(0)(0,0),MXGP);
                basis::tri(log2p).intgrt(&lf(1)(0),&res(1)(0,0),MXGP);
                basis::tri(log2p).intgrt(&lf(2)(0),&res(2)(0,0),MXGP);
                basis::tri(log2p).intgrt(&lf(3)(0),&res(3)(0,0),MXGP);

                /* CALCULATE RESIDUAL TO GOVERNING EQUATION & STORE IN RES */
                basis::tri(log2p).derivr(cv00[0],&res(0)(0,0),MXGP);
                basis::tri(log2p).derivs(cv01[0],&res(0)(0,0),MXGP);
                basis::tri(log2p).derivr(cv10[0],&res(1)(0,0),MXGP);
                basis::tri(log2p).derivs(cv11[0],&res(1)(0,0),MXGP);
                basis::tri(log2p).derivr(cv20[0],&res(2)(0,0),MXGP);
                basis::tri(log2p).derivs(cv21[0],&res(2)(0,0),MXGP);
                basis::tri(log2p).derivr(&du(3,0)(0,0),&res(3)(0,0),MXGP);
                basis::tri(log2p).derivs(&du(3,1)(0,0),&res(3)(0,0),MXGP);
                

                /* THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES */
                for(i=0;i<lgpx;++i) {
                    for(j=0;j<lgpn;++j) {
                        tres[0] = gbl->tau(tind,0)*res(0)(i,j);
                        tres[1] = gbl->tau(tind,0)*res(1)(i,j);
                        tres[2] = gbl->tau(tind,0)*res(2)(i,j);
                        tres[3] = gbl->tau(tind,NV-1)*res(3)(i,j);


                        e00[i][j] -= (dcrd(1,1)(i,j)*(2*u(0)(i,j)-mvel(0)(i,j))
                                        -dcrd(0,1)(i,j)*(u(1)(i,j)-mvel(1)(i,j)))*tres[0]
                                        -dcrd(0,1)(i,j)*u(0)(i,j)*tres[1]
                                        +dcrd(1,1)(i,j)*tres[3];
                        e01[i][j] -= (-dcrd(1,0)(i,j)*(2*u(0)(i,j)-mvel(0)(i,j))
                                        +dcrd(0,0)(i,j)*(u(1)(i,j)-mvel(1)(i,j)))*tres[0]
                                        +dcrd(0,0)(i,j)*u(0)(i,j)*tres[1]
                                        -dcrd(1,0)(i,j)*tres[3];
                        e10[i][j] -= +dcrd(1,1)(i,j)*u(1)(i,j)*tres[0]
                                        +(dcrd(1,1)(i,j)*(u(0)(i,j)-mvel(0)(i,j))
                                        -dcrd(0,1)(i,j)*(2.*u(1)(i,j)-mvel(1)(i,j)))*tres[1]
                                        -dcrd(0,1)(i,j)*tres[3];
                        e11[i][j] -= -dcrd(1,0)(i,j)*u(1)(i,j)*tres[0]
                                +(-dcrd(1,0)(i,j)*(u(0)(i,j)-mvel(0)(i,j))
                                +dcrd(0,0)(i,j)*(2.*u(1)(i,j)-mvel(1)(i,j)))*tres[1]
                                +dcrd(0,0)(i,j)*tres[3];
                        e20[i][j] -= +dcrd(1,1)(i,j)*u(2)(i,j)*tres[0]
                                    -dcrd(0,1)(i,j)*u(2)(i,j)*tres[1]
                                    +(dcrd(1,1)(i,j)*(u(0)(i,j)-mvel(0)(i,j))
                                    -dcrd(0,1)(i,j)*(u(1)(i,j)-mvel(1)(i,j)))*tres[2];
                        e21[i][j] -= -dcrd(1,0)(i,j)*u(2)(i,j)*tres[0]
                                    +dcrd(0,0)(i,j)*u(2)(i,j)*tres[1]
                                    +(-dcrd(1,0)(i,j)*(u(0)(i,j)-mvel(0)(i,j))
                                    +dcrd(0,0)(i,j)*(u(1)(i,j)-mvel(1)(i,j)))*tres[2];

                                        
                        du(3,0)(i,j) = -(dcrd(1,1)(i,j)*tres[0] -dcrd(0,1)(i,j)*tres[1]);
                        du(3,1)(i,j) = -(-dcrd(1,0)(i,j)*tres[0] +dcrd(0,0)(i,j)*tres[1]);
                    }
                }
                basis::tri(log2p).intgrtrs(&lf(0)(0),e00[0],e01[0],MXGP);
                basis::tri(log2p).intgrtrs(&lf(1)(0),e10[0],e11[0],MXGP);
                basis::tri(log2p).intgrtrs(&lf(2)(0),e20[0],e21[0],MXGP);
                basis::tri(log2p).intgrtrs(&lf(3)(0),&du(3,0)(0,0),&du(3,1)(0,0),MXGP);

                for(n=0;n<NV;++n)
                    for(i=0;i<basis::tri(log2p).tm;++i)
                        lf(n)(i) *= gbl->beta(stage);
                        
                lftog(tind,gbl->res_r);
            }
        }
        else {
            /* LINEAR ELEMENT */
            for(i=0;i<lgpx;++i) {
                for(j=0;j<lgpn;++j) {

                    fluxx = gbl->rho*RAD(crd(0)(i,j))*(u(0)(i,j) -mvel(0)(i,j));
                    fluxy = gbl->rho*RAD(crd(0)(i,j))*(u(1)(i,j) -mvel(1)(i,j));
                    
                    /* CONTINUITY EQUATION FLUXES */
                    du(3,0)(i,j) = +ldcrd(1,1)*fluxx -ldcrd(0,1)*fluxy;
                    du(3,1)(i,j) = -ldcrd(1,0)*fluxx +ldcrd(0,0)*fluxy;

                    /* U-MOMENTUM */
                    cv00[i][j] =  u(0)(i,j)*du(3,0)(i,j) +ldcrd(1,1)*RAD(crd(0)(i,j))*u(3)(i,j);
                    cv01[i][j] =  u(0)(i,j)*du(3,1)(i,j) -ldcrd(1,0)*RAD(crd(0)(i,j))*u(3)(i,j);
                    /* V-MOMENTUM */
                    cv10[i][j] =  u(1)(i,j)*du(3,0)(i,j) -ldcrd(0,1)*RAD(crd(0)(i,j))*u(3)(i,j);
                    cv11[i][j] =  u(1)(i,j)*du(3,1)(i,j) +ldcrd(0,0)*RAD(crd(0)(i,j))*u(3)(i,j);
                    /* W-MOMENTUM */
                    cv20[i][j] =  u(2)(i,j)*du(3,0)(i,j);
                    cv21[i][j] =  u(2)(i,j)*du(3,1)(i,j);
                }
            }
            basis::tri(log2p).intgrtrs(&lf(0)(0),cv00[0],cv01[0],MXGP);
            basis::tri(log2p).intgrtrs(&lf(1)(0),cv10[0],cv11[0],MXGP);
            basis::tri(log2p).intgrtrs(&lf(2)(0),cv20[0],cv21[0],MXGP);
            basis::tri(log2p).intgrtrs(&lf(3)(0),&du(3,0)(0,0),&du(3,1)(0,0),MXGP);

            /* ASSEMBLE GLOBAL FORCING (IMAGINARY TERMS) */
            lftog(tind,gbl->res);

            /* NEGATIVE REAL TERMS */
            if (gbl->beta(stage) > 0.0) {
                /* TIME DERIVATIVE TERMS */ 
                for(i=0;i<lgpx;++i) {
                    for(j=0;j<lgpn;++j) {
                        cjcb = ldcrd(0,0)*ldcrd(1,1) -ldcrd(1,0)*ldcrd(0,1);
                        rhorbd0 = rhobd0*RAD(crd(0)(i,j))*cjcb;
                        cjcbi = lmu*RAD(crd(0)(i,j))/cjcb;
                        
                        /* UNSTEADY TERMS */
                        res(0)(i,j) = rhorbd0*u(0)(i,j) +dugdt(log2p,tind,0)(i,j);
                        res(1)(i,j) = rhorbd0*u(1)(i,j) +dugdt(log2p,tind,1)(i,j);
                        res(2)(i,j) = rhorbd0*u(2)(i,j) +dugdt(log2p,tind,2)(i,j);
                        res(3)(i,j) = rhorbd0 +dugdt(log2p,tind,3)(i,j);
                        
#ifdef AXISYMMETRIC
                        /* SOURCE TERMS */
                        res(0)(i,j) -= cjcb*(u(3)(i,j) -2.*lmu*u(0)(i,j)/crd(0)(i,j) +gbl->rho*u(2)(i,j)*u(2)(i,j));
                        //res(2)(i,j) -= cjcb*(lmu*((ldcrd(1,1)*du(2,1)(i,j) -ldcrd(1,0)*du(2,1)(i,j))/cjcb -u(2)(i,j)/crd(0)(i,j)) -gbl->rho*u(0)(i,j)*u(2)(i,j));	
                        //res(2)(i,j) += lmu*(ldcrd(1,1)*du(2,0)(i,j)-ldcrd(1,0)*du(2,1)(i,j))-lmu*cjcb*u(2)(i,j)/crd(0)(i,j) +cjcb*gbl->rho*u(0)(i,j)*u(2)(i,j);
                        //res(2)(i,j) += cjcb*gbl->rho*u(0)(i,j)*u(2)(i,j) -cjcb*lmu*((ldcrd(1,1)*du(2,0)(i,j)-ldcrd(1,0)*du(2,1)(i,j)/cjcb) -u(2)(i,j)/crd(0)(i,j));
                        res(2)(i,j) += cjcb*gbl->rho*u(0)(i,j)*u(2)(i,j) -lmu*ldcrd(1,1)*du(2,0)(i,j) +lmu*ldcrd(1,0)*du(2,1)(i,j) +cjcb*lmu*u(2)(i,j)/crd(0)(i,j);
#else
						res(2)(i,j) -= cjcb*RAD(crd(0)(i,j))*dpdz;
#endif                        
                        /* Source Terms for Test Function */
                        //res(0)(i,j) -= crd(0)(i,j)*cjcb*(gbl->rho*pow(crd(0)(i,j),3.0) +2*crd(0)(i,j) -3*lmu);
                        
                        //res(1)(i,j) -= crd(0)(i,j)*cjcb*(6*gbl->rho*pow(crd(0)(i,j),2.0)*crd(1)(i,j) +3/crd(0)(i,j)*lmu*crd(1)(i,j));
                        
                        //res(2)(i,j) -= crd(0)(i,j)*cjcb*(3*gbl->rho*pow(crd(0)(i,j),3.0) -3*lmu);
                        

                        /* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH) */
                        /* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S) */
                        visc[0][0][0][0] = -cjcbi*(2.*ldcrd(1,1)*ldcrd(1,1) +ldcrd(0,1)*ldcrd(0,1));
                        visc[0][0][1][1] = -cjcbi*(2.*ldcrd(1,0)*ldcrd(1,0) +ldcrd(0,0)*ldcrd(0,0));
                        visc[0][0][0][1] =  cjcbi*(2.*ldcrd(1,1)*ldcrd(1,0) +ldcrd(0,1)*ldcrd(0,0));
#define                 viscI0II0II1II0I visc[0][0][0][1]

                        visc[1][1][0][0] = -cjcbi*(ldcrd(1,1)*ldcrd(1,1) +2.*ldcrd(0,1)*ldcrd(0,1));
                        visc[1][1][1][1] = -cjcbi*(ldcrd(1,0)*ldcrd(1,0) +2.*ldcrd(0,0)*ldcrd(0,0));
                        visc[1][1][0][1] =  cjcbi*(ldcrd(1,1)*ldcrd(1,0) +2.*ldcrd(0,1)*ldcrd(0,0));
#define                 viscI1II1II1II0I visc[1][1][0][1]
                        
                        visc[0][1][0][0] =  cjcbi*ldcrd(0,1)*ldcrd(1,1);
                        visc[0][1][1][1] =  cjcbi*ldcrd(0,0)*ldcrd(1,0);
                        visc[0][1][0][1] = -cjcbi*ldcrd(0,1)*ldcrd(1,0);
                        visc[0][1][1][0] = -cjcbi*ldcrd(0,0)*ldcrd(1,1);
                  
                        visc[2][2][0][0] = -cjcbi*(ldcrd(0,1)*ldcrd(0,1) +ldcrd(1,1)*ldcrd(1,1));
                        visc[2][2][0][1] =  cjcbi*(ldcrd(0,0)*ldcrd(0,1) +ldcrd(1,0)*ldcrd(1,1));
                        visc[2][2][1][1] = -cjcbi*(ldcrd(0,0)*ldcrd(0,0) +ldcrd(1,0)*ldcrd(1,0));
#define					viscI2II2II1II0I visc[2][2][0][1]
						
                        /* OTHER SYMMETRIES     */                
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
                  
                        e20[i][j] = +visc[2][2][0][0]*du(2,0)(i,j) +visc[2][2][0][1]*du(2,1)(i,j);
#ifdef AXISYMMETRIC
						e20[i][j] += lmu*ldcrd(1,1)*u(2)(i,j);
#endif
                        e21[i][j] = +viscI2II2II1II0I*du(2,0)(i,j) +visc[2][2][1][1]*du(2,1)(i,j);
#ifdef AXISYMMETRIC
						e20[i][j] -= lmu*ldcrd(1,0)*u(2)(i,j);
#endif                
                                        
                        cv00[i][j] += e00[i][j];
                        cv01[i][j] += e01[i][j];
                        cv10[i][j] += e10[i][j];
                        cv11[i][j] += e11[i][j];
                        cv20[i][j] += e20[i][j];
                        cv21[i][j] += e21[i][j];
                     }
                }
                basis::tri(log2p).intgrt(&lf(0)(0),&res(0)(0,0),MXGP);
                basis::tri(log2p).intgrt(&lf(1)(0),&res(1)(0,0),MXGP);
                basis::tri(log2p).intgrt(&lf(2)(0),&res(2)(0,0),MXGP);
                basis::tri(log2p).intgrt(&lf(3)(0),&res(3)(0,0),MXGP);

                /* CALCULATE RESIDUAL TO GOVERNING EQUATION & STORE IN RES */
                basis::tri(log2p).derivr(cv00[0],&res(0)(0,0),MXGP);
                basis::tri(log2p).derivs(cv01[0],&res(0)(0,0),MXGP);
                basis::tri(log2p).derivr(cv10[0],&res(1)(0,0),MXGP);
                basis::tri(log2p).derivs(cv11[0],&res(1)(0,0),MXGP);
                basis::tri(log2p).derivr(cv20[0],&res(2)(0,0),MXGP);
                basis::tri(log2p).derivs(cv21[0],&res(2)(0,0),MXGP);
                basis::tri(log2p).derivr(&du(3,0)(0,0),&res(3)(0,0),MXGP);
                basis::tri(log2p).derivs(&du(3,1)(0,0),&res(3)(0,0),MXGP);
                

                /* THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES */
                for(i=0;i<lgpx;++i) {
                    for(j=0;j<lgpn;++j) {
                        tres[0] = gbl->tau(tind,0)*res(0)(i,j);
                        tres[1] = gbl->tau(tind,0)*res(1)(i,j);
                        tres[2] = gbl->tau(tind,0)*res(2)(i,j);
                        tres[3] = gbl->tau(tind,NV-1)*res(3)(i,j);


                        e00[i][j] -= (ldcrd(1,1)*(2*u(0)(i,j)-mvel(0)(i,j))
                                        -ldcrd(0,1)*(u(1)(i,j)-mvel(1)(i,j)))*tres[0]
                                        -ldcrd(0,1)*u(0)(i,j)*tres[1]
                                        +ldcrd(1,1)*tres[3];
                        e01[i][j] -= (-ldcrd(1,0)*(2*u(0)(i,j)-mvel(0)(i,j))
                                        +ldcrd(0,0)*(u(1)(i,j)-mvel(1)(i,j)))*tres[0]
                                        +ldcrd(0,0)*u(0)(i,j)*tres[1]
                                        -ldcrd(1,0)*tres[3];
                        e10[i][j] -= +ldcrd(1,1)*u(1)(i,j)*tres[0]
                                        +(ldcrd(1,1)*(u(0)(i,j)-mvel(0)(i,j))
                                        -ldcrd(0,1)*(2.*u(1)(i,j)-mvel(1)(i,j)))*tres[1]
                                        -ldcrd(0,1)*tres[3];
                        e11[i][j] -= -ldcrd(1,0)*u(1)(i,j)*tres[0]
                                +(-ldcrd(1,0)*(u(0)(i,j)-mvel(0)(i,j))
                                +ldcrd(0,0)*(2.*u(1)(i,j)-mvel(1)(i,j)))*tres[1]
                                +ldcrd(0,0)*tres[3];
                        e20[i][j] -= +ldcrd(1,1)*u(2)(i,j)*tres[0]
                                    -ldcrd(0,1)*u(2)(i,j)*tres[1]
                                    +(ldcrd(1,1)*(u(0)(i,j)-mvel(0)(i,j))
                                    -ldcrd(0,1)*(u(1)(i,j)-mvel(1)(i,j)))*tres[2];
                        e21[i][j] -= -ldcrd(1,0)*u(2)(i,j)*tres[0]
                                    +ldcrd(0,0)*u(2)(i,j)*tres[1]
                                    +(-ldcrd(1,0)*(u(0)(i,j)-mvel(0)(i,j))
                                    +ldcrd(0,0)*(u(1)(i,j)-mvel(1)(i,j)))*tres[2];

                                        
                        du(3,0)(i,j) = -(ldcrd(1,1)*tres[0] -ldcrd(0,1)*tres[1]);
                        du(3,1)(i,j) = -(-ldcrd(1,0)*tres[0] +ldcrd(0,0)*tres[1]);
                    }
                }
                basis::tri(log2p).intgrtrs(&lf(0)(0),e00[0],e01[0],MXGP);
                basis::tri(log2p).intgrtrs(&lf(1)(0),e10[0],e11[0],MXGP);
                basis::tri(log2p).intgrtrs(&lf(2)(0),e20[0],e21[0],MXGP);
                basis::tri(log2p).intgrtrs(&lf(3)(0),&du(3,0)(0,0),&du(3,1)(0,0),MXGP);

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
            
    /*********************************************/
    /* MODIFY RESIDUALS ON COARSER MESHES            */
    /*********************************************/    
    if (coarse_flag) {
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

    return;
}
