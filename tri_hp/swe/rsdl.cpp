/*
 *  rsdl.cpp
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_swe.h"
#include "../hp_boundary.h"
    
block::ctrl tri_hp_swe::rsdl(block::ctrl ctrl_message, int stage) {
    int i,j,n,tind;
    FLT fluxx,fluxy,pres;
    const int NV = 3;
    TinyVector<int,3> v;
    TinyMatrix<FLT,ND,ND> ldcrd;
    TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,NV,ND> du;
    int lgpx = basis::tri(log2p).gpx, lgpn = basis::tri(log2p).gpn;
    FLT cjcb, oneminusbeta;
    TinyMatrix<TinyMatrix<FLT,ND,ND>,NV-1,NV-1> visc;
    TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,NV-1,NV-1> cv, df;
    TinyVector<FLT,ND> vel,pt;
    TinyVector<FLT,NV> tres;
    FLT drag;
    block::ctrl state;
    
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
#ifdef CTRL_DEBUG
            *sim::log << idprefix << " step 0 of ins::rsdl: ctrl_message: " << ctrl_message << " excpt: " << excpt << " stage: " << stage << std::endl;
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
            *sim::log << idprefix << " step 1 of ins::rsdl: ctrl_message: " << ctrl_message << " excpt: " << excpt << " stage: " << stage << std::endl;
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
            *sim::log << idprefix << " step 1 of ins::rsdl: ctrl_message: " << ctrl_message << " excpt: " << excpt << " stage: " << stage << std::endl;
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
            *sim::log << idprefix << " step 2 of ins::rsdl: ctrl_message: " << ctrl_message << " excpt: " << excpt << " stage: " << stage << std::endl;
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
                    }
                }

                /* LOAD SOLUTION COEFFICIENTS FOR THIS ELEMENT */
                /* PROJECT SOLUTION TO GAUSS POINTS WITH DERIVATIVES IF NEEDED FOR VISCOUS TERMS */
                ugtouht(tind);
                for(n=0;n<NV;++n)
                    basis::tri(log2p).proj(&uht(n)(0),&u(n)(0,0),MXGP);
                
                /* lf IS WHERE I WILL STORE THE ELEMENT RESIDUAL */
                for(n=0;n<NV;++n)
                    for(i=0;i<basis::tri(log2p).tm;++i)
                        lf(n)(i) = 0.0;

                if (td(tind).info > -1) {
                    /* CURVED ELEMENT */
                    /* CONVECTIVE TERMS (IMAGINARY FIRST)*/
                    for(i=0;i<lgpx;++i) {
                        for(j=0;j<lgpn;++j) {

                            fluxx = u(0)(i,j) -u(NV-1)(i,j)*mvel(0)(i,j);
                            fluxy = u(1)(i,j) -u(NV-1)(i,j)*mvel(1)(i,j);
                            vel(0) = u(0)(i,j)/u(NV-1)(i,j);
                            vel(1) = u(1)(i,j)/u(NV-1)(i,j);
                            pres = u(NV-1)(i,j)*u(NV-1)(i,j)*sim::g/2.0;
                            
                            /* CONTINUITY EQUATION FLUXES */
                            du(NV-1,0)(i,j) = +dcrd(1,1)(i,j)*fluxx -dcrd(0,1)(i,j)*fluxy;
                            du(NV-1,1)(i,j) = -dcrd(1,0)(i,j)*fluxx +dcrd(0,0)(i,j)*fluxy;
                          
                            /* CONVECTIVE FLUXES */
                            for(n=0;n<NV-1;++n) {
                                cv(n,0)(i,j) = vel(n)*du(NV-1,0)(i,j);
                                cv(n,1)(i,j) = vel(n)*du(NV-1,1)(i,j);
                            }

                            /* PRESSURE TERMS */
                            /* U-MOMENTUM */
                            cv(0,0)(i,j) += dcrd(1,1)(i,j)*pres;
                            cv(0,1)(i,j) -= dcrd(1,0)(i,j)*pres;
                            /* V-MOMENTUM */
                            cv(1,0)(i,j) -=  dcrd(0,1)(i,j)*pres;
                            cv(1,1)(i,j) +=  dcrd(0,0)(i,j)*pres;
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
                                pt(0) = crd(0)(i,j);
                                pt(1) = crd(1)(i,j);
                                
                                drag = gbl_ptr->cd*sqrt(u(0)(i,j)*u(0)(i,j) +u(1)(i,j)*u(1)(i,j))/(u(NV-1)(i,j)*u(NV-1)(i,j));
                                
                                /* UNSTEADY TERMS */
                                res(0)(i,j) = cjcb*(sim::bd[0]*u(0)(i,j) -u(NV-1)(i,j)*sim::g*gbl_ptr->bathy->f(0,pt) -(gbl_ptr->f0 +gbl_ptr->beta*crd(1)(i,j))*u(1)(i,j) +drag*u(0)(i,j)) +dugdt(log2p,tind,0)(i,j);
                                res(1)(i,j) = cjcb*(sim::bd[0]*u(1)(i,j) -u(NV-1)(i,j)*sim::g*gbl_ptr->bathy->f(1,pt) +(gbl_ptr->f0 +gbl_ptr->beta*crd(1)(i,j))*u(0)(i,j) +drag*u(1)(i,j)) +dugdt(log2p,tind,1)(i,j);                                
                                res(NV-1)(i,j) = cjcb*sim::bd[0]*u(NV-1)(i,j) +dugdt(log2p,tind,NV-1)(i,j);
                                
                                /* TEMPORARY TO MAINTAIN FREE-STREAM SOLUTION */
//                                TinyVector<FLT,mesh::ND> xtemp;
//                                xtemp(0) = 0.0;
//                                xtemp(1) = 0.0;
//                                res(0)(i,j) += cjcb*(-gbl_ptr->cd*pow(gbl_ptr->ibc->f(0,xtemp),2));
//                                res(1)(i,j) += cjcb*(-(gbl_ptr->f0 +gbl_ptr->beta*crd(1)(i,j)))*gbl_ptr->ibc->f(0,xtemp);                                
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
                                            
                                tres(0) = -gbl_ptr->tau(tind,0)*res(0)(i,j);
                                tres(1) = -gbl_ptr->tau(tind,1)*res(1)(i,j);
                                tres(NV-1) = -gbl_ptr->tau(tind,NV-1)*res(NV-1)(i,j);
                                
                                vel(0) = u(0)(i,j)/u(NV-1)(i,j);
                                vel(1) = u(1)(i,j)/u(NV-1)(i,j);

                                df(0,0)(i,j) = (dcrd(1,1)(i,j)*(2*vel(0)-mvel(0)(i,j))
                                                -dcrd(0,1)(i,j)*(vel(1)-mvel(1)(i,j)))*tres(0)
                                                -dcrd(0,1)(i,j)*vel(0)*tres(1)
                                                +dcrd(1,1)(i,j)*(sim::g*u(NV-1)(i,j) -vel(0)*vel(0))*tres(NV-1);
                                df(0,1)(i,j) = (-dcrd(1,0)(i,j)*(2*vel(0)-mvel(0)(i,j))
                                                +dcrd(0,0)(i,j)*(vel(1)-mvel(1)(i,j)))*tres(0)
                                                +dcrd(0,0)(i,j)*vel(0)*tres(1)
                                                -dcrd(1,0)(i,j)*(sim::g*u(NV-1)(i,j) -vel(0)*vel(0))*tres(NV-1);
                                df(1,0)(i,j) = +dcrd(1,1)(i,j)*vel(1)*tres(0)
                                                +(dcrd(1,1)(i,j)*(vel(0)-mvel(0)(i,j))
                                                -dcrd(0,1)(i,j)*(2.*vel(1)-mvel(1)(i,j)))*tres(1)
                                                -dcrd(0,1)(i,j)*(sim::g*u(NV-1)(i,j) -vel(1)*vel(1))*tres(NV-1);
                                df(1,1)(i,j) = -dcrd(1,0)(i,j)*vel(1)*tres(0)
                                                +(-dcrd(1,0)(i,j)*(vel(0)-mvel(0)(i,j))
                                                +dcrd(0,0)(i,j)*(2.*vel(1)-mvel(1)(i,j)))*tres(1)
                                                +dcrd(0,0)(i,j)*(sim::g*u(NV-1)(i,j) -vel(1)*vel(1))*tres(NV-1);
                                                
                                du(NV-1,0)(i,j) = (dcrd(1,1)(i,j)*tres(0) -dcrd(0,1)(i,j)*tres(1) -(dcrd(1,1)(i,j)*mvel(0)(i,j) -dcrd(0,1)(i,j)*mvel(1)(i,j))*tres(2));
                                du(NV-1,1)(i,j) = (-dcrd(1,0)(i,j)*tres(0) +dcrd(0,0)(i,j)*tres(1) -(-dcrd(1,0)(i,j)*mvel(0)(i,j) +dcrd(0,0)(i,j)*mvel(1)(i,j))*tres(2));
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

                            fluxx = u(0)(i,j) -u(NV-1)(i,j)*mvel(0)(i,j);
                            fluxy = u(1)(i,j) -u(NV-1)(i,j)*mvel(1)(i,j);
                            vel(0) = u(0)(i,j)/u(NV-1)(i,j);
                            vel(1) = u(1)(i,j)/u(NV-1)(i,j);
                            pres = u(NV-1)(i,j)*u(NV-1)(i,j)*sim::g/2.0;
                            
                            /* CONTINUITY EQUATION FLUXES */
                            du(NV-1,0)(i,j) = +ldcrd(1,1)*fluxx -ldcrd(0,1)*fluxy;
                            du(NV-1,1)(i,j) = -ldcrd(1,0)*fluxx +ldcrd(0,0)*fluxy;
                          
                            /* CONVECTIVE FLUXES */
                            for(n=0;n<NV-1;++n) {
                                cv(n,0)(i,j) = vel(n)*du(NV-1,0)(i,j);
                                cv(n,1)(i,j) = vel(n)*du(NV-1,1)(i,j);
                            }

                            /* PRESSURE TERMS */
                            /* U-MOMENTUM */
                            cv(0,0)(i,j) += ldcrd(1,1)*pres;
                            cv(0,1)(i,j) -= ldcrd(1,0)*pres;
                            /* V-MOMENTUM */
                            cv(1,0)(i,j) -=  ldcrd(0,1)*pres;
                            cv(1,1)(i,j) +=  ldcrd(0,0)*pres;
                        }
                    }
                    for(n=0;n<NV-1;++n)
                        basis::tri(log2p).intgrtrs(&lf(n)(0),&cv(n,0)(0,0),&cv(n,1)(0,0),MXGP);
                    basis::tri(log2p).intgrtrs(&lf(NV-1)(0),&du(NV-1,0)(0,0),&du(NV-1,1)(0,0),MXGP);
                    

                    /* ASSEMBLE GLOBAL FORCING (IMAGINARY TERMS) */
                    lftog(tind,gbl_ptr->res);

                    /* NEGATIVE REAL TERMS */
                    if (sim::beta[stage] > 0.0) {
                        cjcb = ldcrd(0,0)*ldcrd(1,1) -ldcrd(1,0)*ldcrd(0,1);
                        
                        /* TIME DERIVATIVE TERMS */ 
                        for(i=0;i<lgpx;++i) {
                            for(j=0;j<lgpn;++j) {
                                pt(0) = crd(0)(i,j);
                                pt(1) = crd(1)(i,j);

                                drag = gbl_ptr->cd*sqrt(u(0)(i,j)*u(0)(i,j) +u(1)(i,j)*u(1)(i,j))/(u(NV-1)(i,j)*u(NV-1)(i,j));
                                /* UNSTEADY TERMS */
                                res(0)(i,j) = cjcb*(sim::bd[0]*u(0)(i,j) -u(NV-1)(i,j)*sim::g*gbl_ptr->bathy->f(0,pt) -(gbl_ptr->f0 +gbl_ptr->beta*crd(1)(i,j))*u(1)(i,j) +drag*u(0)(i,j)) +dugdt(log2p,tind,0)(i,j);
                                res(1)(i,j) = cjcb*(sim::bd[0]*u(1)(i,j) -u(NV-1)(i,j)*sim::g*gbl_ptr->bathy->f(1,pt) +(gbl_ptr->f0 +gbl_ptr->beta*crd(1)(i,j))*u(0)(i,j) +drag*u(1)(i,j)) +dugdt(log2p,tind,1)(i,j);                                
                                res(NV-1)(i,j) = cjcb*sim::bd[0]*u(NV-1)(i,j) +dugdt(log2p,tind,NV-1)(i,j);
                  
                                
                                /* TEMPORARY TO MAINTAIN FREE-STREAM SOLUTION */
//                                TinyVector<FLT,mesh::ND> xtemp;
//                                xtemp(0) = 0.0;
//                                xtemp(1) = 0.0;
//                                res(0)(i,j) += cjcb*(-gbl_ptr->cd*pow(gbl_ptr->ibc->f(0,xtemp),2));
//                                res(1)(i,j) += cjcb*(-(gbl_ptr->f0 +gbl_ptr->beta*crd(1)(i,j)))*gbl_ptr->ibc->f(0,xtemp);     
                                
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
                                            
                                tres(0) = -gbl_ptr->tau(tind,0)*res(0)(i,j);
                                tres(1) = -gbl_ptr->tau(tind,1)*res(1)(i,j);
                                tres(NV-1) = -gbl_ptr->tau(tind,NV-1)*res(NV-1)(i,j);
                                
                                vel(0) = u(0)(i,j)/u(NV-1)(i,j);
                                vel(1) = u(1)(i,j)/u(NV-1)(i,j);

                                df(0,0)(i,j) = (ldcrd(1,1)*(2*vel(0)-mvel(0)(i,j))
                                                -ldcrd(0,1)*(vel(1)-mvel(1)(i,j)))*tres(0)
                                                -ldcrd(0,1)*vel(0)*tres(1)
                                                +ldcrd(1,1)*(sim::g*u(NV-1)(i,j) -vel(0)*vel(0))*tres(NV-1);
                                df(0,1)(i,j) = (-ldcrd(1,0)*(2*vel(0)-mvel(0)(i,j))
                                                +ldcrd(0,0)*(vel(1)-mvel(1)(i,j)))*tres(0)
                                                +ldcrd(0,0)*vel(0)*tres(1)
                                                -ldcrd(1,0)*(sim::g*u(NV-1)(i,j) -vel(0)*vel(0))*tres(NV-1);
                                df(1,0)(i,j) = +ldcrd(1,1)*vel(1)*tres(0)
                                                +(ldcrd(1,1)*(vel(0)-mvel(0)(i,j))
                                                -ldcrd(0,1)*(2.*vel(1)-mvel(1)(i,j)))*tres(1)
                                                -ldcrd(0,1)*(sim::g*u(NV-1)(i,j) -vel(1)*vel(1))*tres(NV-1);
                                df(1,1)(i,j) = -ldcrd(1,0)*vel(1)*tres(0)
                                                +(-ldcrd(1,0)*(vel(0)-mvel(0)(i,j))
                                                +ldcrd(0,0)*(2.*vel(1)-mvel(1)(i,j)))*tres(1)
                                                +ldcrd(0,0)*(sim::g*u(NV-1)(i,j) -vel(1)*vel(1))*tres(NV-1);
                                                
                                du(NV-1,0)(i,j) = (ldcrd(1,1)*tres(0) -ldcrd(0,1)*tres(1) -(ldcrd(1,1)*mvel(0)(i,j) -ldcrd(0,1)*mvel(1)(i,j))*tres(2));
                                du(NV-1,1)(i,j) = (-ldcrd(1,0)*tres(0) +ldcrd(0,0)*tres(1) -(-ldcrd(1,0)*mvel(0)(i,j) +ldcrd(0,0)*mvel(1)(i,j))*tres(2));
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
            /* MODIFY RESIDUALS ON COARSER MESHES            */
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
            
//            for(i=0;i<nvrtx;++i)
//                printf("rsdl v: %d %e %e %e\n",i,gbl_ptr->res.v(i,0),gbl_ptr->res.v(i,1),gbl_ptr->res.v(i,2));
//                
//            for(i=0;i<nside;++i)
//                for(int m=0;m<basis::tri(log2p).sm;++m)
//                    printf("rsdl s: %d %d %e %e %e\n",i,m,gbl_ptr->res.s(i,m,0),gbl_ptr->res.s(i,m,1),gbl_ptr->res.s(i,m,2));
//
//            for(i=0;i<ntri;++i)
//                for(int m=0;m<basis::tri(log2p).im;++m)
//                    printf("rsdl i: %d %d %e %e %e\n",i,m,gbl_ptr->res.i(i,m,0),gbl_ptr->res.i(i,m,1),gbl_ptr->res.i(i,m,2));

        }
    }

    return(block::stop);
}
