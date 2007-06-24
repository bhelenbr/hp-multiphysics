#include <math.h>

#include "tri_hp_buoyancy.h"
#include "../hp_boundary.h"

block::ctrl tri_hp_buoyancy::setup_preconditioner(block::ctrl ctrl_message) {
    int tind,i,j,side,v0;
    FLT jcb,jcb1,h,hmax,q,qmax,lam1,lam2,gam,rhoav;
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
            rhoav = 0.0;
            for(j=0;j<3;++j) {
                v0 = v(j);
                q = pow(ug.v(v0,0) -0.5*(sim::bd[0]*(vrtx(v0)(0) -vrtxbd(1)(v0)(0))),2.0) 
                    +pow(ug.v(v0,1) -0.5*(sim::bd[0]*(vrtx(v0)(1) -vrtxbd(1)(v0)(1))),2.0);
                qmax = MAX(qmax,q);
                rhoav += gbl_ptr->rhovsT.Eval(ug.v(v0,2)); 
            }
            rhoav /= 3.0;
            FLT nu = gbl_ptr->mu/rhoav;
            FLT alpha = gbl_ptr->kcond/(rhoav*gbl_ptr->cp);
        
            gam = 3.0*qmax +(0.5*hmax*sim::bd[0] +2.*nu/hmax)*(0.5*hmax*sim::bd[0] +2.*nu/hmax);
            if (gbl_ptr->mu + sim::bd[0] == 0.0) gam = MAX(gam,0.01);
            q = sqrt(qmax);
            lam1 = q + sqrt(qmax +gam);
            lam2  = (q +1.5*alpha/h +hmax*sim::bd[0]);

            /* SET UP DISSIPATIVE COEFFICIENTS */
            gbl_ptr->tau(tind,0) = adis*h/(jcb*sqrt(gam));
            gbl_ptr->tau(tind,2)  = adis*h/(jcb*lam2);
            gbl_ptr->tau(tind,NV-1) = qmax*gbl_ptr->tau(tind,0);
                        
            /* SET UP DIAGONAL PRECONDITIONER */
            jcb1 = 2.5*jcb*lam1/h;
            jcb1 *= RAD((vrtx(v(0))(0) +vrtx(v(1))(0) +vrtx(v(2))(0))/3.);

            // jcb *= 8.*nu*(1./(hmax*hmax) +1./(h*h)) +2*lam1/h +2*sqrt(gam)/hmax +sim::bd[0];
            jcb *= 2.*nu*(1./(hmax*hmax) +1./(h*h)) +3*lam1/h;  // heuristically tuned
            jcb *= RAD((vrtx(v(0))(0) +vrtx(v(1))(0) +vrtx(v(2))(0))/3.);

            gbl_ptr->tprcn(tind,0) = rhoav*jcb;    
            gbl_ptr->tprcn(tind,1) = rhoav*jcb;  
            gbl_ptr->tprcn(tind,2) =  rhoav*gbl_ptr->cp*jcb1;      
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

    return(ctrl_message);
}
