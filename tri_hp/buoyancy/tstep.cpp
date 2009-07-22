#include <math.h>

#include "tri_hp_buoyancy.h"
#include "../hp_boundary.h"

void tri_hp_buoyancy::setup_preconditioner() {
	int tind,i,j,side,v0;
	FLT jcb,jcb1,h,hmax,q,qmax,lam1,lam2,gam,rhoav;
	TinyVector<int,3> v;

	/***************************************/
	/** DETERMINE FLOW PSEUDO-TIME STEP ****/
	/***************************************/
	gbl->vprcn(Range(0,npnt-1),Range::all()) = 0.0;
	if (basis::tri(log2p)->sm() > 0) {
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
		h = 4.*jcb/(0.25*(basis::tri(log2p)->p() +1)*(basis::tri(log2p)->p()+1)*hmax);
		hmax = hmax/(0.25*(basis::tri(log2p)->p() +1)*(basis::tri(log2p)->p()+1));

		qmax = 0.0;
		rhoav = 0.0;
		for(j=0;j<3;++j) {
			v0 = v(j);
			q = pow(ug.v(v0,0) -0.5*(gbl->bd(0)*(pnts(v0)(0) -vrtxbd(1)(v0)(0))),2.0) 
				+pow(ug.v(v0,1) -0.5*(gbl->bd(0)*(pnts(v0)(1) -vrtxbd(1)(v0)(1))),2.0);
			qmax = MAX(qmax,q);
			rhoav += gbl->rhovsT.Eval(ug.v(v0,2)); 
		}
		rhoav /= 3.0;
		FLT nu = gbl->mu/rhoav;
		FLT alpha = gbl->kcond/(rhoav*gbl->cp);

		gam = 3.0*qmax +(0.5*hmax*gbl->bd(0) +2.*nu/hmax)*(0.5*hmax*gbl->bd(0) +2.*nu/hmax);
		if (gbl->mu + gbl->bd(0) == 0.0) gam = MAX(gam,0.01);
		q = sqrt(qmax);
		lam1 = q + sqrt(qmax +gam);
		lam2  = (q +1.5*alpha/h +hmax*gbl->bd(0));

		/* SET UP DISSIPATIVE COEFFICIENTS */
		gbl->tau(tind,0) = adis*h/(jcb*sqrt(gam));
		gbl->tau(tind,2)  = adis*h/(jcb*lam2);
		gbl->tau(tind,NV-1) = qmax*gbl->tau(tind,0);

		/* SET UP DIAGONAL PRECONDITIONER */
		jcb1 = 2.5*jcb*lam1/h;
		jcb1 *= RAD((pnts(v(0))(0) +pnts(v(1))(0) +pnts(v(2))(0))/3.);

		// jcb *= 8.*nu*(1./(hmax*hmax) +1./(h*h)) +2*lam1/h +2*sqrt(gam)/hmax +gbl->bd(0);
		jcb *= 2.*nu*(1./(hmax*hmax) +1./(h*h)) +3*lam1/h;  // heuristically tuned
		jcb *= RAD((pnts(v(0))(0) +pnts(v(1))(0) +pnts(v(2))(0))/3.);

		gbl->tprcn(tind,0) = rhoav*jcb;    
		gbl->tprcn(tind,1) = rhoav*jcb;  
		gbl->tprcn(tind,2) =  rhoav*gbl->cp*jcb1;      
		gbl->tprcn(tind,3) =  jcb/gam;
		for(i=0;i<3;++i) {
			gbl->vprcn(v(i),Range::all())  += gbl->tprcn(tind,Range::all());
			if (basis::tri(log2p)->sm() > 0) {
				side = tri(tind).seg(i);
				gbl->sprcn(side,Range::all()) += gbl->tprcn(tind,Range::all());
			}
		}
	}
	tri_hp::setup_preconditioner();

	return;
}
