#include "tri_hp_lvlset.h"
#include "../hp_boundary.h"
#include <math.h>

block::ctrl tri_hp_lvlset::setup_preconditioner(block::ctrl ctrl_message) {
    
    if (gbl_ptr->diagonal_preconditioner) {
        /* SET-UP DIAGONAL PRECONDITIONER */
        int tind,i,j,side,v0;
        FLT jcb,jcbphi,h,hmax,q,qmax,lam1,lam2,gam,rho,rhomax,mu,nu,heavy,delt,deltamax,strss;
        TinyVector<int,3> v;
            
        if (ctrl_message == block::begin) excpt = 0;
        if (ctrl_message == block::advance1) ++excpt;    
        
        if (excpt == 3) {
            ++excpt;
            
            /***************************************/
            /** DETERMINE FLOW PSEUDO-TIME STEP ****/
            /***************************************/
            gbl_ptr->vprcn(Range(0,nvrtx-1),Range::all()) = 0.0;
            if (basis::tri(log2p).sm > 0) {
                gbl_ptr->sprcn(Range(0,nside-1),Range::all()) = 0.0;
            }

            for(tind = 0; tind < ntri; ++tind) {
                jcb = 0.25*area(tind);  // area is 2 x triangle area
                v = td(tind).vrtx;
                hmax = 0.0;
                for(j=0;j<3;++j) {
                    h = pow(vrtx(v(j))(0) -vrtx(v((j+1)%3))(0),2.0) + 
                    pow(vrtx(v(j))(1) -vrtx(v((j+1)%3))(1),2.0);
                    hmax = (h > hmax ? h : hmax);
                }
                hmax = sqrt(hmax);
                
                if (!(jcb > 0.0)) {  // THIS CATCHES NAN'S TOO
                    *sim::log << "negative triangle area caught in tstep. Problem triangle is : " << tind << std::endl;
                    *sim::log << "approximate location: " << vrtx(v(0))(0) << ' ' << vrtx(v(0))(1) << std::endl;
                    mesh::output("negative",grid);
                    exit(1);
                }
                h = 4.*jcb/(0.25*(basis::tri(log2p).p +1)*(basis::tri(log2p).p+1)*hmax);
                hmax = hmax/(0.25*(basis::tri(log2p).p +1)*(basis::tri(log2p).p+1));
            
                qmax = 0.0;
                nu = 0.0;
                rhomax = 0.0;
                deltamax = 0.0;
                for(j=0;j<3;++j) {
                    v0 = v(j);
                    q = pow(ug.v(v0,0)-0.5*(sim::bd[0]*(vrtx(v0)(0) -vrtxbd(1)(v0)(0))),2.0) 
                        +pow(ug.v(v0,1)-0.5*(sim::bd[0]*(vrtx(v0)(1) -vrtxbd(1)(v0)(1))),2.0);
                    q += pow(ug.v(v0,2),2.0);
                    qmax = MAX(qmax,q);
                    heavyside_and_delta_if(ug.v(v0,2)/gbl_ptr->width,heavy,delt);
                    rho = gbl_ptr->rho +(gbl_ptr->rho2 -gbl_ptr->rho)*heavy;
                    mu = gbl_ptr->mu +(gbl_ptr->mu2 -gbl_ptr->mu)*heavy;
                    nu = MAX(nu,mu/rho);
                    rhomax = MAX(rhomax,rho); 
                    deltamax = MAX(deltamax,delt);              
                }
                gam = 3.0*qmax +(0.5*hmax*sim::bd[0] +2.*nu/hmax)*(0.5*hmax*sim::bd[0] +2.*nu/hmax);
                if (nu + sim::bd[0] == 0.0) gam = MAX(gam,0.01);
                q = sqrt(qmax);
                lam1 = q + sqrt(qmax +gam);
                
                strss =  4.*gbl_ptr->sigma*deltamax/h; // +fabs(drho*sim::g*nrm(1)/h);
                lam2 = q/h +0.5*sim::bd[0] +strss/(gbl_ptr->rho*lam1);
                delt *= gbl_ptr->width;
                lam2 = lam2*delt +2.0/h*(1.0 -delt);
                
                /* SET UP DISSIPATIVE COEFFICIENTS */
                gbl_ptr->tau(tind,0) = adis*h/(jcb*sqrt(gam));
                gbl_ptr->tau(tind,2) = adis/(jcb*lam2);
                gbl_ptr->tau(tind,NV-1) = qmax*gbl_ptr->tau(tind,0);
                
                /* SET UP DIAGONAL PRECONDITIONER */
                jcbphi = jcb*lam2;
                jcbphi *= RAD((vrtx(v(0))(0) +vrtx(v(1))(0) +vrtx(v(2))(0))/3.);

                // jcb *= 8.*nu*(1./(hmax*hmax) +1./(h*h)) +2*lam1/h +2*sqrt(gam)/hmax +sim::bd[0];
                jcb *= 2.*nu*(1./(hmax*hmax) +1./(h*h)) +3*lam1/h;  // heuristically tuned
                jcb *= RAD((vrtx(v(0))(0) +vrtx(v(1))(0) +vrtx(v(2))(0))/3.);
                
                gbl_ptr->tprcn(tind,0) = rhomax*jcb;    
                gbl_ptr->tprcn(tind,1) = rhomax*jcb;  
                gbl_ptr->tprcn(tind,2) = jcbphi;      
                gbl_ptr->tprcn(tind,3) =  jcb/gam;
                for(i=0;i<3;++i) {
                    gbl_ptr->vprcn(v(i),Range::all())  += gbl_ptr->tprcn(tind,Range::all());
                    if (basis::tri(log2p).sm > 0) {
                        side = td(tind).side(i);
                        gbl_ptr->sprcn(side,Range::all()) += gbl_ptr->tprcn(tind,Range::all());
                    }
                }
            }
        }
        else {
            ctrl_message = tri_hp::setup_preconditioner(ctrl_message);  
        }
    }
//    else {
//        /* SET-UP MATRIX PRECONDITIONER */
//        int tind,i,j,side,v0;
//        FLT jcb,jcbphi,h,hmax,q,qmax,lam1,lam2,gam,ubar,vbar,phibar,rho,rhomax,mu,nu,heavy,delt,deltamax,strss;
//        TinyVector<int,3> v;
//        
//        if (ctrl_message == block::begin) excpt = 0;
//        if (ctrl_message == block::advance1) ++excpt;    
//        
//        if (excpt == 3) {
//            ++excpt;
//                          
//            /***************************************/
//            /** DETERMINE FLOW PSEUDO-TIME STEP ****/
//            /***************************************/
//            gbl_ptr->vprcn_ut(Range(0,nvrtx-1),Range::all(),Range::all()) = 0.0;
//            if (basis::tri(log2p).sm > 0) {
//                gbl_ptr->sprcn_ut(Range(0,nside-1),Range::all(),Range::all()) = 0.0;
//            }
//            gbl_ptr->tprcn_ut(Range(0,ntri-1),Range::all(),Range::all()) = 0.0;
//
//
//            for(tind = 0; tind < ntri; ++tind) {
//                jcb = 0.25*area(tind);  // area is 2 x triangle area
//                v = td(tind).vrtx;
//                hmax = 0.0;
//                for(j=0;j<3;++j) {
//                    h = pow(vrtx(v(j))(0) -vrtx(v((j+1)%3))(0),2.0) + 
//                    pow(vrtx(v(j))(1) -vrtx(v((j+1)%3))(1),2.0);
//                    hmax = (h > hmax ? h : hmax);
//                }
//                hmax = sqrt(hmax);
//                
//                if (!(jcb > 0.0)) {  // THIS CATCHES NAN'S TOO
//                    *sim::log << "negative triangle area caught in tstep. Problem triangle is : " << tind << std::endl;
//                    *sim::log << "approximate location: " << vrtx(v(0))(0) << ' ' << vrtx(v(0))(1) << std::endl;
//                    mesh::output("negative",grid);
//                    exit(1);
//                }
//                h = 4.*jcb/(0.25*(basis::tri(log2p).p +1)*(basis::tri(log2p).p+1)*hmax);
//                hmax = hmax/(0.25*(basis::tri(log2p).p +1)*(basis::tri(log2p).p+1));
//            
//                qmax = 0.0;
//                ubar = 0.0;
//                vbar = 0.0;
//                phibar = 0.0;
//                nu = 0.0;
//                rhomax = 0.0;
//                deltamax = 0.0;
//                for(j=0;j<3;++j) {
//                    v0 = v(j);
//                    q = pow(ug.v(v0,0)-(sim::bd[0]*(vrtx(v0)(0) -vrtxbd(1)(v0)(0))),2.0) 
//                        +pow(ug.v(v0,1)-(sim::bd[0]*(vrtx(v0)(1) -vrtxbd(1)(v0)(1))),2.0);  
//                    ubar += ug.v(v0,0);
//                    vbar += ug.v(v0,1);
//                    phibar += ug.v(v0,2);
//                    qmax = MAX(qmax,q);
//                    heavyside_and_delta_if(ug.v(v0,2)/gbl_ptr->width,heavy,delt);
//                    rho = gbl_ptr->rho +(gbl_ptr->rho2 -gbl_ptr->rho)*heavy;
//                    mu = gbl_ptr->mu +(gbl_ptr->mu2 -gbl_ptr->mu)*heavy;
//                    nu = MAX(nu,mu/rho);
//                    rhomax = MAX(rhomax,rho); 
//                    deltamax = MAX(deltamax,delt);              
//                }
//                
//                phibar /= 3.0;
//                ubar /= 3.0;
//                vbar /= 3.0;
//
//                gam = 3.0*qmax +(0.5*hmax*sim::bd[0] +2.*nu/hmax)*(0.5*hmax*sim::bd[0] +2.*nu/hmax);
//                if (nu + sim::bd[0] == 0.0) gam = MAX(gam,0.1);
//
//                q = sqrt(qmax);
//                lam1 = q + sqrt(qmax +gam);          
//                
//                strss =  4.*gbl_ptr->sigma*deltamax/h; // +fabs(drho*sim::g*nrm(1)/h);
//                lam2 = q/h +0.5*sim::bd[0] +strss/(gbl_ptr->rho*lam1);
//                
//                /* SET UP DISSIPATIVE COEFFICIENTS */
//                gbl_ptr->tau(tind,0) = adis*h/(jcb*sqrt(gam));
//                gbl_ptr->tau(tind,2) = adis/(jcb*lam2);
//                gbl_ptr->tau(tind,NV-1) = qmax*gbl_ptr->tau(tind,0);
//                
//                /* SET UP DIAGONAL PRECONDITIONER */
//                delt *= gbl_ptr->width;
//                jcbphi = jcb*(lam2*delt +1.0*(1.0 -delt));
//                jcbphi *= RAD((vrtx(v(0))(0) +vrtx(v(1))(0) +vrtx(v(2))(0))/3.);
//                
//                /* SET UP DIAGONAL PRECONDITIONER */
//                // jcb *= 8.*nu*(1./(hmax*hmax) +1./(h*h)) +2*lam1/h +2*sqrt(gam)/hmax +sim::bd[0];
//                jcb *= 2.*nu*(1./(hmax*hmax) +1./(h*h)) +3*lam1/h;  // heuristically tuned
//                jcb *= RAD((vrtx(v(0))(0) +vrtx(v(1))(0) +vrtx(v(2))(0))/3.);
//
//                gbl_ptr->tprcn_ut(tind,0,0) = gbl_ptr->rho*jcb;    
//                gbl_ptr->tprcn_ut(tind,1,1) = gbl_ptr->rho*jcb;  
//                gbl_ptr->tprcn_ut(tind,2,2) = gbl_ptr->rho*jcbphi;    
//                gbl_ptr->tprcn_ut(tind,3,3) = jcb/gam;
//                gbl_ptr->tprcn_ut(tind,0,3) = jcb*ubar/gam;
//                gbl_ptr->tprcn_ut(tind,1,3) = jcb*vbar/gam;
//                gbl_ptr->tprcn_ut(tind,2,3) = jcb*phibar/gam;
//                for(i=0;i<3;++i) {
//                    gbl_ptr->vprcn_ut(v(i),Range::all(),Range::all())  += basis::tri(log2p).vdiag*gbl_ptr->tprcn_ut(tind,Range::all(),Range::all());
//                    if (basis::tri(log2p).sm > 0) {
//                        side = td(tind).side(i);
//                        gbl_ptr->sprcn_ut(side,Range::all(),Range::all()) += gbl_ptr->tprcn_ut(tind,Range::all(),Range::all());
//                    }
//                }
//            }
//        }
//        else {
//            ctrl_message = tri_hp::setup_preconditioner(ctrl_message);  
//        }
//    }

    return(ctrl_message);
}