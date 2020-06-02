//
//  rsdl.cpp
//  tri_hp
//
//  Created by Brian Helenbrook on 5/30/20.
//

#include <stdio.h>
/*
 *  rsdl.cpp
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_explicit.h"
#include "../hp_boundary.h"

void tri_hp_explicit::element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uht,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im){
    FLT visc[ND][ND][ND][ND], tres[NV];
    FLT cv00[MXGP][MXGP],cv01[MXGP][MXGP];
    FLT e00[MXGP][MXGP],e01[MXGP][MXGP];
    const int lgpx = basis::tri(log2p)->gpx(), lgpn = basis::tri(log2p)->gpn();
    TinyVector<TinyMatrix<FLT,MXGP,MXGP>,ND> mvel; // for local mesh velocity info

    /* IF TINFO > -1 IT IS CURVED ELEMENT */
    if (tri(tind).info > -1) {
        /* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
        crdtocht(tind);

        /* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
        for(int n=0;n<ND;++n)
            basis::tri(log2p)->proj_bdry(&cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), &dcrd(n,1)(0,0),MXGP);
    }
    else {
        /* LOAD INDICES OF VERTEX POINTS */
        TinyVector<int,3> v;
        v = tri(tind).pnt;
        
        /* PROJECT VERTEX COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
        for(int n=0;n<ND;++n)
            basis::tri(log2p)->proj(pnts(v(0))(n),pnts(v(1))(n),pnts(v(2))(n),&crd(n)(0,0),MXGP);

        /* CALCULATE COORDINATE DERIVATIVES A SIMPLE WAY */
        for(int i=0;i<basis::tri(log2p)->gpx();++i) {
            for(int j=0;j<basis::tri(log2p)->gpn();++j) {
                for(int n=0;n<ND;++n) {
                    dcrd(n,0)(i,j) = 0.5*(pnts(v(2))(n) -pnts(v(1))(n));
                    dcrd(n,1)(i,j) = 0.5*(pnts(v(0))(n) -pnts(v(1))(n));
                }
            }
        }
    }

    /* CALCULATE MESH VELOCITY */
    for(int i=0;i<lgpx;++i) {
        for(int j=0;j<lgpn;++j) {
            mvel(0)(i,j) = gbl->bd(0)*(crd(0)(i,j) -dxdt(log2p)(tind,0,i,j));
            mvel(1)(i,j) = gbl->bd(0)*(crd(1)(i,j) -dxdt(log2p)(tind,1,i,j));
#ifdef MESH_REF_VEL
            mvel(0)(i,j) += gbl->mesh_ref_vel(0);
            mvel(1)(i,j) += gbl->mesh_ref_vel(1);
#endif
        }
    }

    basis::tri(log2p)->proj(&uht(0)(0),&u(0)(0,0),&du(0,0)(0,0),&du(0,1)(0,0),MXGP);

    /* lf IS WHERE I WILL STORE THE ELEMENT RESIDUAL */
    for(int n=0;n<NV;++n){
        for(int i=0;i<basis::tri(log2p)->tm();++i){
            lf_re(n)(i) = 0.0;
            lf_im(n)(i) = 0.0;
        }
    }

    /* CONVECTION */
    for(int i=0;i<lgpx;++i) {
        for(int j=0;j<lgpn;++j) {

#ifdef CONST_A
            FLT fluxx = gbl->rhocv*RAD(crd(0)(i,j))*(gbl->ax -mvel(0)(i,j))*u(0)(i,j);
            FLT fluxy = gbl->rhocv*RAD(crd(0)(i,j))*(gbl->ay -mvel(1)(i,j))*u(0)(i,j);
#else
            pt(0) = crd(0)(i,j);
            pt(1) = crd(1)(i,j);
            FLT fluxx = gbl->rhocv*RAD(crd(0)(i,j))*(gbl->a->f(0,pt,gbl->time) -mvel(0)(i,j))*u(0)(i,j);
            FLT fluxy = gbl->rhocv*RAD(crd(0)(i,j))*(gbl->a->f(1,pt,gbl->time) -mvel(1)(i,j))*u(0)(i,j);
#endif
            cv00[i][j] = +dcrd(1,1)(i,j)*fluxx -dcrd(0,1)(i,j)*fluxy;
            cv01[i][j] = -dcrd(1,0)(i,j)*fluxx +dcrd(0,0)(i,j)*fluxy;
        }
    }
    basis::tri(log2p)->intgrtrs(&lf_im(0)(0),&cv00[0][0],&cv01[0][0],MXGP);



    /* NEGATIVE REAL TERMS */
    if (gbl->beta(stage) > 0.0) {

        /* TIME DERIVATIVE TERMS */
        for(int i=0;i<lgpx;++i) {
            for(int j=0;j<lgpn;++j) {
                TinyVector<FLT,ND> pt(crd(0)(i,j),crd(1)(i,j));
                cjcb(i,j) = dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j);
                res(0)(i,j) = RAD(crd(0)(i,j))*gbl->sigma*u(0)(i,j)*cjcb(i,j);
                res(0)(i,j) -= RAD(crd(0)(i,j))*cjcb(i,j)*gbl->src->f(0,pt,gbl->time);
            }
        }
        basis::tri(log2p)->intgrt(&lf_re(0)(0),&res(0)(0,0),MXGP);

        /* DIFFUSIVE TERMS  */
        for(int i=0;i<lgpx;++i) {
            for(int j=0;j<lgpn;++j) {

                cjcb(i,j) = gbl->kcond*RAD(crd(0)(i,j))/cjcb(i,j);

                /* DIFFUSION TENSOR (LOTS OF SYMMETRY THOUGH)*/
                /* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
                visc[0][0][0][0] = -cjcb(i,j)*(dcrd(1,1)(i,j)*dcrd(1,1)(i,j) +dcrd(0,1)(i,j)*dcrd(0,1)(i,j));
                visc[0][0][1][1] = -cjcb(i,j)*(dcrd(1,0)(i,j)*dcrd(1,0)(i,j) +dcrd(0,0)(i,j)*dcrd(0,0)(i,j));
                visc[0][0][0][1] =  cjcb(i,j)*(dcrd(1,1)(i,j)*dcrd(1,0)(i,j) +dcrd(0,1)(i,j)*dcrd(0,0)(i,j));
#define          viscI0II0II1II0I visc[0][0][0][1]

                e00[i][j] = +visc[0][0][0][0]*du(0,0)(i,j) +visc[0][0][0][1]*du(0,1)(i,j);
                e01[i][j] = +viscI0II0II1II0I*du(0,0)(i,j) +visc[0][0][1][1]*du(0,1)(i,j);

                cv00[i][j] += e00[i][j];
                cv01[i][j] += e01[i][j];
            }
        }
        basis::tri(log2p)->derivr(cv00[0],&res(0)(0,0),MXGP);
        basis::tri(log2p)->derivs(cv01[0],&res(0)(0,0),MXGP);

        /* THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES */
        for(int i=0;i<lgpx;++i) {
            for(int j=0;j<lgpn;++j) {
                TinyVector<FLT,ND> pt(crd(0)(i,j),crd(1)(i,j));
                tres[0] = gbl->tau(tind)*res(0)(i,j);

#ifdef CONST_A
                e00[i][j] -= (dcrd(1,1)(i,j)*(gbl->ax-mvel(0)(i,j))
                    -dcrd(0,1)(i,j)*(gbl->ay-mvel(1)(i,j)))*tres[0];
                e01[i][j] -= (-dcrd(1,0)(i,j)*(gbl->ax-mvel(0)(i,j))
                    +dcrd(0,0)(i,j)*(gbl->ay-mvel(1)(i,j)))*tres[0];
#else
                e00[i][j] -= (dcrd(1,1)(i,j)*(gbl->a->f(0,pt,gbl->time)-mvel(0)(i,j))
                                            -dcrd(0,1)(i,j)*(gbl->a->f(1,pt,gbl->time)-mvel(1)(i,j)))*tres[0];
                e01[i][j] -= (-dcrd(1,0)(i,j)*(gbl->a->f(0,pt,gbl->time)-mvel(0)(i,j))
                                            +dcrd(0,0)(i,j)*(gbl->a->f(1,pt,gbl->time)-mvel(1)(i,j)))*tres[0];
#endif
            }
        }
        basis::tri(log2p)->intgrtrs(&lf_re(0)(0),e00[0],e01[0],MXGP);

        for(int n=0;n<NV;++n)
            for(int i=0;i<basis::tri(log2p)->tm();++i)
                lf_re(n)(i) *= gbl->beta(stage);

    }
    

    return;
}
