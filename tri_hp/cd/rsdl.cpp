/*
 *  rsdl.cpp
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_cd.h"
#include "../hp_boundary.h"

    
void tri_hp_cd::rsdl(int stage) {
    int i,j,n,tind;
    FLT fluxx,fluxy;
    FLT visc[ND][ND][ND][ND], tres[NV];
    FLT cv00[MXGP][MXGP],cv01[MXGP][MXGP];
    FLT e00[MXGP][MXGP],e01[MXGP][MXGP];
    FLT oneminusbeta;
    TinyVector<FLT,ND> pt;
    TinyVector<int,3> v;
    int lgpx = basis::tri(log2p).gpx, lgpn = basis::tri(log2p).gpn;

    tri_hp::rsdl(stage);
    
    oneminusbeta = 1.0-sim::beta[stage];

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
            for(i=0;i<basis::tri(log2p).gpx;++i) {
                for(j=0;j<basis::tri(log2p).gpn;++j) {
                    for(n=0;n<ND;++n) {
                        dcrd(n,0)(i,j) = 0.5*(pnts(v(2))(n) -pnts(v(1))(n));
                        dcrd(n,1)(i,j) = 0.5*(pnts(v(0))(n) -pnts(v(1))(n));
                    }
                }
            }
        }

        /* CALCULATE MESH VELOCITY */
        for(i=0;i<lgpx;++i) {
            for(j=0;j<lgpn;++j) {
                mvel(0)(i,j) = gbl->bd[0]*(crd(0)(i,j) -dxdt(log2p,tind,0)(i,j));
                mvel(1)(i,j) = gbl->bd[0]*(crd(1)(i,j) -dxdt(log2p,tind,1)(i,j));
            }
        }
        
        ugtouht(tind);
        basis::tri(log2p).proj(&uht(0)(0),&u(0)(0,0),&du(0,0)(0,0),&du(0,1)(0,0),MXGP);

        for(n=0;n<NV;++n)
            for(i=0;i<basis::tri(log2p).tm;++i)
                lf(n)(i) = 0.0;
                
        /* CONVECTION */
        for(i=0;i<basis::tri(log2p).gpx;++i) {
            for(j=0;j<basis::tri(log2p).gpn;++j) {

                fluxx = RAD(crd(0)(i,j))*(gbl->ax -mvel(0)(i,j))*u(0)(i,j);
                fluxy = RAD(crd(0)(i,j))*(gbl->ay -mvel(1)(i,j))*u(0)(i,j);
                
                cv00[i][j] = +dcrd(1,1)(i,j)*fluxx -dcrd(0,1)(i,j)*fluxy;
                cv01[i][j] = -dcrd(1,0)(i,j)*fluxx +dcrd(0,0)(i,j)*fluxy;
            }
        }
        basis::tri(log2p).intgrtrs(&lf(0)(0),&cv00[0][0],&cv01[0][0],MXGP);

        /* ASSEMBLE GLOBAL FORCING (IMAGINARY TERMS) */
        lftog(tind,gbl->res);

        /* NEGATIVE REAL TERMS */
        if (sim::beta[stage] > 0.0) {
                
            /* TIME DERIVATIVE TERMS */
            for(i=0;i<basis::tri(log2p).gpx;++i) {
                for(j=0;j<basis::tri(log2p).gpn;++j) {
                    pt(0) = crd(0)(i,j);
                    pt(1) = crd(1)(i,j);
                    cjcb(i,j) = dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j);
                    res(0)(i,j) = RAD(crd(0)(i,j))*gbl->bd[0]*u(0)(i,j)*cjcb(i,j) +dugdt(log2p,tind,0)(i,j);
                    res(0)(i,j) -= RAD(crd(0)(i,j))*cjcb(i,j)*gbl->src->f(0,pt,gbl->time);
                }
            }            
            basis::tri(log2p).intgrt(&lf(0)(0),&res(0)(0,0),MXGP);

            /* DIFFUSIVE TERMS  */
            for(i=0;i<basis::tri(log2p).gpx;++i) {
                for(j=0;j<basis::tri(log2p).gpn;++j) {

                    cjcb(i,j) = gbl->nu*RAD(crd(0)(i,j))/cjcb(i,j);
                    
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
            basis::tri(log2p).derivr(cv00[0],&res(0)(0,0),MXGP);
            basis::tri(log2p).derivs(cv01[0],&res(0)(0,0),MXGP);
            
            /* THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES */
            for(i=0;i<basis::tri(log2p).gpx;++i) {
                for(j=0;j<basis::tri(log2p).gpn;++j) {

                    tres[0] = gbl->tau(tind)*res(0)(i,j);

                    e00[i][j] -= (dcrd(1,1)(i,j)*(gbl->ax-mvel(0)(i,j))
                                     -dcrd(0,1)(i,j)*(gbl->ay-mvel(1)(i,j)))*tres[0];
                    e01[i][j] -= (-dcrd(1,0)(i,j)*(gbl->ax-mvel(0)(i,j))
                                      +dcrd(0,0)(i,j)*(gbl->ay-mvel(1)(i,j)))*tres[0];
              }
            }
            basis::tri(log2p).intgrtrs(&lf(0)(0),e00[0],e01[0],MXGP); 

            for(n=0;n<NV;++n)
                for(i=0;i<basis::tri(log2p).tm;++i)
                    lf(n)(i) *= sim::beta[stage];
                    
            lftog(tind,gbl->res_r);
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
