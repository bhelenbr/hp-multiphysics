#include <math.h>

#include "tet_hp_ins.h"
#include "../hp_boundary.h"

void tet_hp_ins::setup_preconditioner() {
	if (gbl->diagonal_preconditioner) {
		/* SET-UP DIAGONAL PRECONDITIONER */
		int tind,i,j,side,v0,v1,v2,find;
		FLT jcb,h,hmax,q,qmax,lam1,gam,a,amax,amin,dx1,dx2,dy1,dy2,dz1,dz2,cpi,cpj,cpk;
		TinyVector<int,4> v;
		
		FLT nu = gbl->mu/gbl->rho;

		/***************************************/
		/** DETERMINE FLOW PSEUDO-TIME STEP ****/
		/***************************************/
		gbl->vprcn(Range(0,npnt-1),Range::all()) = 0.0;
		if (basis::tet(log2p).em > 0) {
			gbl->eprcn(Range(0,nseg-1),Range::all()) = 0.0;
		}
		if (basis::tet(log2p).fm > 0) {
			gbl->fprcn(Range(0,ntri-1),Range::all()) = 0.0;
		}
		
#ifdef TIMEACCURATE
		gam = 10.0;
		FLT dtstari = 0.0;
#endif

		for(tind = 0; tind < ntet; ++tind) {
			jcb = tet(tind).vol/8.0;   // area is 2 x triangle area
			v = tet(tind).pnt;	
			amin=1000000.0;
			amax = 0.0;
			for(j=0;j<4;++j) { // FIND MAX FACE AREA AND THEN DIVIDE VOLUME BY IT 
				find = tet(tind).tri(j);
				v0 = tri(find).pnt(0);
				v1 = tri(find).pnt(1);
				v2 = tri(find).pnt(2);
				
				dx1 = pnts(v0)(0)-pnts(v1)(0);
				dy1 = pnts(v0)(1)-pnts(v1)(1);
				dz1 = pnts(v0)(2)-pnts(v1)(2);
				dx2 = pnts(v0)(0)-pnts(v2)(0);
				dy2 = pnts(v0)(1)-pnts(v2)(1);
				dz2 = pnts(v0)(2)-pnts(v2)(2);
				cpi = dy1*dz2-dz1*dy2;
				cpj = -dx1*dz2+dz1*dx2;
				cpk = dx1*dy2-dy1*dx2;
				a =	.5*sqrt(cpi*cpi+cpj*cpj+cpk*cpk);
				amax = (a > amax ? a : amax);
				amin = (a < amin ? a : amin);

			}
			
			h = 4.0*jcb/(0.25*(basis::tet(log2p).p+1)*(basis::tet(log2p).p+1)*amax); // 3*8/6=4
			hmax = 4.0*jcb/(0.25*(basis::tet(log2p).p+1)*(basis::tet(log2p).p+1)*amin); // 3*8/6=4
		
			qmax = 0.0;
			for(j=0;j<4;++j) {
				v0 = v(j);
//				cout << "v0 = "<< v0 << " pnts = " << pnts(v0) << " vrtxbd = " << vrtxbd(1)(v0) << endl;
//				cout << "ug.v = " << ug.v(v0,0) << ' ' << ug.v(v0,1) << ' ' << ug.v(v0,2) << endl;
				q = pow(ug.v(v0,0)-0.5*(gbl->bd(0)*(pnts(v0)(0) -vrtxbd(1)(v0)(0))),2.0) 
					+pow(ug.v(v0,1)-0.5*(gbl->bd(0)*(pnts(v0)(1) -vrtxbd(1)(v0)(1))),2.0)  
					+pow(ug.v(v0,2)-0.5*(gbl->bd(0)*(pnts(v0)(2) -vrtxbd(1)(v0)(2))),2.0);
				qmax = MAX(qmax,q);
			}
			cout << "qmax = " << qmax << endl;
			if (!(jcb > 0.0)) { 
				*gbl->log << "negative tetrahedral volume caught in tstep. Problem tet is : " << tind << std::endl;
				*gbl->log << "approximate location: " << pnts(v(0))(0) << ' ' << pnts(v(0))(1)<< ' ' << pnts(v(0))(2) << std::endl;
				tet_mesh::output("negative",grid);
				exit(1);
			}

			if  (!(qmax >= 0.0)) {  // THIS CATCHES NAN'S TOO
				*gbl->log << "flow solution has nan's" << std::endl;
				output("nan",tecplot);
				exit(1);
			}

#ifndef INERTIALESS

#ifndef TIMEACCURATE
			gam = 3.0*qmax +(0.5*hmax*gbl->bd(0) +2.*nu/hmax)*(0.5*hmax*gbl->bd(0) +2.*nu/hmax);
			if (gbl->mu + gbl->bd(0) == 0.0) gam = MAX(gam,0.1);
#endif
			q = sqrt(qmax);
			lam1 = q + sqrt(qmax +gam);

			/* SET UP DISSIPATIVE COEFFICIENTS */
			gbl->tau(tind,0) = adis*h/(jcb*sqrt(gam));
			gbl->tau(tind,NV-1) = qmax*gbl->tau(tind,0);

			/* SET UP DIAGONAL PRECONDITIONER */
			// jcb *= 8.*nu*(1./(hmax*hmax) +1./(h*h)) +2*lam1/h +2*sqrt(gam)/hmax +gbl->bd(0);
			jcb *= 2.*nu*(1./(hmax*hmax) +1./(h*h)) +3*lam1/h;  // heuristically tuned
#else
			gam = pow(2.*nu/hmax,2); 
			lam1 = sqrt(gam);

			/* SET UP DISSIPATIVE COEFFICIENTS */
			gbl->tau(tind,0)  = adis*h/(jcb*sqrt(gam));
			gbl->tau(tind,NV-1) = 0.0;

			jcb *= 8.*nu*(1./(hmax*hmax) +1./(h*h)) +2*lam1/h +2*sqrt(gam)/hmax;
	#endif
	#ifdef TIMEACCURATE
			dtstari = MAX((nu/(h*h) +lam1/h +gbl->bd(0)),dtstari);

		}
		printf("#iterative to physical time step ratio: %f\n",gbl->bd(0)/dtstari);

		for(tind=0;tind<ntet;++tind) {
			v = tet(tind).pnt;
			jcb = 0.125*tet(tind).vol*dtstari;
	#endif

			gbl->iprcn(tind,0) = gbl->rho*jcb;    
			gbl->iprcn(tind,1) = gbl->rho*jcb;
			gbl->iprcn(tind,2) = gbl->rho*jcb;      
			gbl->iprcn(tind,3) =  jcb/gam;
			for(i=0;i<4;++i) 
				gbl->vprcn(v(i),Range::all())  += gbl->iprcn(tind,Range::all());            
			if (basis::tet(log2p).em > 0) {
				for(i=0;i<6;++i){
					side = tet(tind).seg(i);
					gbl->eprcn(side,Range::all()) += gbl->iprcn(tind,Range::all());
				}
			}
			if (basis::tet(log2p).fm > 0) {
				for(i=0;i<4;++i){
					side = tet(tind).tri(i);
					gbl->fprcn(side,Range::all()) += gbl->iprcn(tind,Range::all());
				}
			}
		}
	}
	else {
		cout << "matrix preconditioner being called and doesn't work " << endl;
//        /* SET-UP MATRIX PRECONDITIONER */
//        int tind,i,j,side,v0;
//        FLT jcb,h,hmax,q,qmax,lam1,gam,ubar,vbar;
//        TinyVector<int,4> v;
//        
//        FLT nu = gbl->mu/gbl->rho;
//              
//        /***************************************/
//        /** DETERMINE FLOW PSEUDO-TIME STEP ****/
//        /***************************************/
//        gbl->vprcn_ut(Range(0,npnt-1),Range::all(),Range::all()) = 0.0;
//        if (basis::tet(log2p).em > 0) {
//            gbl->eprcn_ut(Range(0,nseg-1),Range::all(),Range::all()) = 0.0;
//        }
//		if (basis::tet(log2p).fm > 0) {
//            gbl->fprcn_ut(Range(0,ntri-1),Range::all(),Range::all()) = 0.0;
//        }
//        gbl->iprcn_ut(Range(0,ntet-1),Range::all(),Range::all()) = 0.0;
//
//        for(tind = 0; tind < ntet; ++tind) {
//            jcb = tet(tind).vol/8.0; 
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
//            for(j=0;j<3;++j) {
//                v0 = v(j);
//
//                q = pow(ug.v(v0,0)-(gbl->bd(0)*(pnts(v0)(0) -vrtxbd(1)(v0)(0))),2.0) 
//                    +pow(ug.v(v0,1)-(gbl->bd(0)*(pnts(v0)(1) -vrtxbd(1)(v0)(1))),2.0);  
//                ubar += ug.v(v0,0);
//                vbar += ug.v(v0,1);
//                qmax = MAX(qmax,q);
//            }
//            
//            ubar /= 3.0;
//            vbar /= 3.0;
//
//            gam = 3.0*qmax +(0.5*hmax*gbl->bd(0) +2.*nu/hmax)*(0.5*hmax*gbl->bd(0) +2.*nu/hmax);
//            if (gbl->mu + gbl->bd(0) == 0.0) gam = MAX(gam,0.1);
//
//            q = sqrt(qmax);
//            lam1 = q + sqrt(qmax +gam);
//            
//            /* SET UP DISSIPATIVE COEFFICIENTS */
//            gbl->tau(tind,1) = adis*h/(jcb*sqrt(gam));
//            gbl->tau(tind,0) = qmax*gbl->tau(tind,1);
//            
//            /* SET UP DIAGONAL PRECONDITIONER */
//            // jcb *= 8.*nu*(1./(hmax*hmax) +1./(h*h)) +2*lam1/h +2*sqrt(gam)/hmax +gbl->bd(0);
//            jcb *= 2.*nu*(1./(hmax*hmax) +1./(h*h)) +3*lam1/h;  // heuristically tuned
//            jcb *= RAD((pnts(v(0))(0) +pnts(v(1))(0) +pnts(v(2))(0))/3.);
//
//            gbl->tprcn_ut(tind,0,0) = gbl->rho*jcb;    
//            gbl->tprcn_ut(tind,1,1) = gbl->rho*jcb;      
//            gbl->tprcn_ut(tind,2,2) = jcb/gam;
//            gbl->tprcn_ut(tind,0,2) = jcb*ubar/gam;
//            gbl->tprcn_ut(tind,1,2) = jcb*vbar/gam;
//            for(i=0;i<3;++i) {
//                gbl->vprcn_ut(v(i),Range::all(),Range::all())  += basis::tri(log2p).vdiag*gbl->tprcn_ut(tind,Range::all(),Range::all());
//                if (basis::tri(log2p).sm > 0) {
//                    side = tri(tind).seg(i);
//                    gbl->sprcn_ut(side,Range::all(),Range::all()) += gbl->tprcn_ut(tind,Range::all(),Range::all());
//                }
//            }
//        }
	}


	
	tet_hp::setup_preconditioner();
}
