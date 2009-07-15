#include "tri_hp_lvlset.h"
#include "../hp_boundary.h"
#include <math.h>

void tri_hp_lvlset::setup_preconditioner() {

	if (gbl->diagonal_preconditioner) {
		/* SET-UP DIAGONAL PRECONDITIONER */
		int tind,i,j,side,v0;
		FLT jcb,jcbphi,h,hmax,q,qmax,lam1,lam2,gam,rho,rhomax,mu,nu,heavy,delt,deltamax,strss;
		TinyVector<int,3> v;

		/***************************************/
		/** DETERMINE FLOW PSEUDO-TIME STEP ****/
		/***************************************/
		gbl->vprcn(Range(0,npnt-1),Range::all()) = 0.0;
		if (basis::tri(log2p).sm > 0) {
			gbl->sprcn(Range(0,nseg-1),Range::all()) = 0.0;
		}

		for(tind = 0; tind < ntri; ++tind) {
			jcb = 0.25*area(tind);  // area is 2 x triangle area
			v = tri(tind).pnt;
			hmax = 0.0;
			for(j=0;j<3;++j) {
				h = pow(pnts(v(j))(0) -pnts(v((j+1)%3))(0),2.0) + 
				pow(pnts(v(j))(1) -pnts(v((j+1)%3))(1),2.0);
				hmax = (h > hmax ? h : hmax);
			}
			hmax = sqrt(hmax);

			if (!(jcb > 0.0)) {  // THIS CATCHES NAN'S TOO
				*gbl->log << "negative triangle area caught in tstep. Problem triangle is : " << tind << std::endl;
				*gbl->log << "approximate location: " << pnts(v(0))(0) << ' ' << pnts(v(0))(1) << std::endl;
				tri_mesh::output("negative",grid);
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
				q = pow(ug.v(v0,0)-0.5*(gbl->bd(0)*(pnts(v0)(0) -vrtxbd(1)(v0)(0))),2.0) 
					+pow(ug.v(v0,1)-0.5*(gbl->bd(0)*(pnts(v0)(1) -vrtxbd(1)(v0)(1))),2.0);
				qmax = MAX(qmax,q);
				heavyside_and_delta_if(ug.v(v0,2)/gbl->width,heavy,delt);
				rho = gbl->rho +(gbl->rho2 -gbl->rho)*heavy;
				mu = gbl->mu +(gbl->mu2 -gbl->mu)*heavy;
				nu = MAX(nu,mu/rho);
				rhomax = MAX(rhomax,rho); 
				deltamax = MAX(deltamax,delt);              
			}
			gam = 3.0*qmax +(0.5*hmax*gbl->bd(0) +2.*nu/hmax)*(0.5*hmax*gbl->bd(0) +2.*nu/hmax);
			if (nu + gbl->bd(0) == 0.0) gam = MAX(gam,0.01);
			q = sqrt(qmax);
			lam1 = q + sqrt(qmax +gam);

			/* SET UP DISSIPATIVE COEFFICIENTS */
			gbl->tau(tind,0) = adis*h/(jcb*sqrt(gam));
			gbl->tau(tind,NV-1) = qmax*gbl->tau(tind,0);

			strss =  4.*gbl->sigma*deltamax/h; // +fabs(drho*gbl->g*nrm(1)/h);
#ifdef LOCALIZED_WITH_DISTANCE_FUNCTION
			delt *= gbl->width;
			gbl->tau(tind,2) = adis/(jcb*(delt*(q/h +0.5*gbl->bd(0)) +(1.0-delt)/h));
			lam2 = q/h +0.5*gbl->bd(0) +strss/(gbl->rho*lam1);
			lam2 = lam2*delt +2.0/h*(1.0 -delt);
#else
			gbl->tau(tind,2) = adis/(jcb*(q/h +0.5*gbl->bd(0)));
			lam2 = q/h +0.5*gbl->bd(0) +strss/(gbl->rho*lam1);
#endif
#ifdef CONSERVATIVE
			lam2 *= rhomax;
#endif            

			/* SET UP DIAGONAL PRECONDITIONER */
			jcbphi = 2.0*jcb*lam2;
			jcbphi *= RAD((pnts(v(0))(0) +pnts(v(1))(0) +pnts(v(2))(0))/3.);

			// jcb *= 8.*nu*(1./(hmax*hmax) +1./(h*h)) +2*lam1/h +2*sqrt(gam)/hmax +gbl->bd(0);
			jcb *= 2.*nu*(1./(hmax*hmax) +1./(h*h)) +3*lam1/h;  // heuristically tuned
			jcb *= RAD((pnts(v(0))(0) +pnts(v(1))(0) +pnts(v(2))(0))/3.);

			/* BRUTE FORCE WAY TO FREEZE FLOW 
			gbl->tprcn(tind,0) = 1.0e99*gbl->rho*jcb;   
			gbl->tprcn(tind,1) = 1.0e99*gbl->rho*jcb;   
			gbl->tprcn(tind,3) =  1.0e99*jcb/gam; */

			gbl->tprcn(tind,0) = rhomax*jcb; 
			gbl->tprcn(tind,1) = rhomax*jcb;  
			gbl->tprcn(tind,2) = jcbphi;      
			gbl->tprcn(tind,3) = jcb/gam;
			for(i=0;i<3;++i) {
				gbl->vprcn(v(i),Range::all())  += gbl->tprcn(tind,Range::all());
				if (basis::tri(log2p).sm > 0) {
					side = tri(tind).seg(i);
					gbl->sprcn(side,Range::all()) += gbl->tprcn(tind,Range::all());
				}
			}
		}
	}
//    else {
//        /* SET-UP MATRIX PRECONDITIONER */
//        int tind,i,j,side,v0;
//        FLT jcb,jcbphi,h,hmax,q,qmax,lam1,lam2,gam,ubar,vbar,phibar,rho,rhomax,mu,nu,heavy,delt,deltamax,strss;
//        TinyVector<int,3> v;
//        
//                      
//        /***************************************/
//        /** DETERMINE FLOW PSEUDO-TIME STEP ****/
//        /***************************************/
//        gbl->vprcn_ut(Range(0,npnt-1),Range::all(),Range::all()) = 0.0;
//        if (basis::tri(log2p).sm > 0) {
//            gbl->sprcn_ut(Range(0,nseg-1),Range::all(),Range::all()) = 0.0;
//        }
//        gbl->tprcn_ut(Range(0,ntri-1),Range::all(),Range::all()) = 0.0;
//
//
//        for(tind = 0; tind < ntri; ++tind) {
//            jcb = 0.25*area(tind);  // area is 2 x triangle area
//            v = tri(tind).pnt;
//            hmax = 0.0;
//            for(j=0;j<3;++j) {
//                h = pow(pnts(v(j))(0) -pnts(v((j+1)%3))(0),2.0) + 
//                pow(pnts(v(j))(1) -pnts(v((j+1)%3))(1),2.0);
//                hmax = (h > hmax ? h : hmax);
//            }
//            hmax = sqrt(hmax);
//            
//            if (!(jcb > 0.0)) {  // THIS CATCHES NAN'S TOO
//                *gbl->log << "negative triangle area caught in tstep. Problem triangle is : " << tind << std::endl;
//                *gbl->log << "approximate location: " << pnts(v(0))(0) << ' ' << pnts(v(0))(1) << std::endl;
//                tri_mesh::output("negative",grid);
//                exit(1);
//            }
//            h = 4.*jcb/(0.25*(basis::tri(log2p).p +1)*(basis::tri(log2p).p+1)*hmax);
//            hmax = hmax/(0.25*(basis::tri(log2p).p +1)*(basis::tri(log2p).p+1));
//        
//            qmax = 0.0;
//            ubar = 0.0;
//            vbar = 0.0;
//            phibar = 0.0;
//            nu = 0.0;
//            rhomax = 0.0;
//            deltamax = 0.0;
//            for(j=0;j<3;++j) {
//                v0 = v(j);
//                q = pow(ug.v(v0,0)-(gbl->bd(0)*(pnts(v0)(0) -vrtxbd(1)(v0)(0))),2.0) 
//                    +pow(ug.v(v0,1)-(gbl->bd(0)*(pnts(v0)(1) -vrtxbd(1)(v0)(1))),2.0);  
//                ubar += ug.v(v0,0);
//                vbar += ug.v(v0,1);
//                phibar += ug.v(v0,2);
//                qmax = MAX(qmax,q);
//                heavyside_and_delta_if(ug.v(v0,2)/gbl->width,heavy,delt);
//                rho = gbl->rho +(gbl->rho2 -gbl->rho)*heavy;
//                mu = gbl->mu +(gbl->mu2 -gbl->mu)*heavy;
//                nu = MAX(nu,mu/rho);
//                rhomax = MAX(rhomax,rho); 
//                deltamax = MAX(deltamax,delt);              
//            }
//            
//            phibar /= 3.0;
//            ubar /= 3.0;
//            vbar /= 3.0;
//
//            gam = 3.0*qmax +(0.5*hmax*gbl->bd(0) +2.*nu/hmax)*(0.5*hmax*gbl->bd(0) +2.*nu/hmax);
//            if (nu + gbl->bd(0) == 0.0) gam = MAX(gam,0.1);
//
//            q = sqrt(qmax);
//            lam1 = q + sqrt(qmax +gam);          
//            
//            strss =  4.*gbl->sigma*deltamax/h; // +fabs(drho*gbl->g*nrm(1)/h);
//            lam2 = q/h +0.5*gbl->bd(0) +strss/(gbl->rho*lam1);
//            
//            /* SET UP DISSIPATIVE COEFFICIENTS */
//            gbl->tau(tind,0) = adis*h/(jcb*sqrt(gam));
//            gbl->tau(tind,2) = adis/(jcb*lam2);
//            gbl->tau(tind,NV-1) = qmax*gbl->tau(tind,0);
//            
//            /* SET UP DIAGONAL PRECONDITIONER */
//            delt *= gbl->width;
//            jcbphi = jcb*(lam2*delt +1.0*(1.0 -delt));
//            jcbphi *= RAD((pnts(v(0))(0) +pnts(v(1))(0) +pnts(v(2))(0))/3.);
//            
//            /* SET UP DIAGONAL PRECONDITIONER */
//            // jcb *= 8.*nu*(1./(hmax*hmax) +1./(h*h)) +2*lam1/h +2*sqrt(gam)/hmax +gbl->bd(0);
//            jcb *= 2.*nu*(1./(hmax*hmax) +1./(h*h)) +3*lam1/h;  // heuristically tuned
//            jcb *= RAD((pnts(v(0))(0) +pnts(v(1))(0) +pnts(v(2))(0))/3.);
//
//            gbl->tprcn_ut(tind,0,0) = gbl->rho*jcb;    
//            gbl->tprcn_ut(tind,1,1) = gbl->rho*jcb;  
//            gbl->tprcn_ut(tind,2,2) = gbl->rho*jcbphi;    
//            gbl->tprcn_ut(tind,3,3) = jcb/gam;
//            gbl->tprcn_ut(tind,0,3) = jcb*ubar/gam;
//            gbl->tprcn_ut(tind,1,3) = jcb*vbar/gam;
//            gbl->tprcn_ut(tind,2,3) = jcb*phibar/gam;
//            for(i=0;i<3;++i) {
//                gbl->vprcn_ut(v(i),Range::all(),Range::all())  += basis::tri(log2p).vdiag*gbl->tprcn_ut(tind,Range::all(),Range::all());
//                if (basis::tri(log2p).sm > 0) {
//                    side = tri(tind).seg(i);
//                    gbl->sprcn_ut(side,Range::all(),Range::all()) += gbl->tprcn_ut(tind,Range::all(),Range::all());
//                }
//            }
//        }
//    }

	tri_hp::setup_preconditioner();

	return;
}
