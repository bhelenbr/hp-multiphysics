#include <math.h>

#include "tri_hp_thermoelastic.h"
#include "../hp_boundary.h"

/* ELASTICITY PART FOLLOWS planestrain/setup_preconditioner.cpp, TEMPERATURE PART FOLLOWS cd/setup_preconditioner.cpp */
int tri_hp_thermoelastic::setup_preconditioner() {
	int tind,i,j,side,v0;
	FLT jcb,jcbT,h,hmax,hmaxp,lam1,lamT,gam,gami,q,qmax;
	TinyVector<int,3> v;
	TinyVector<FLT,ND> mvel;
	const int T = 2, P = 3;
	int err = 0;
	FLT alpha = (hp_te_gbl->rhocv > 0.0 ? hp_te_gbl->kcond/hp_te_gbl->rhocv : 0.0);

	/***************************************/
	/** DETERMINE PSEUDO-TIME STEP ****/
	/***************************************/
	hp_gbl->vprcn(Range(0,npnt-1),Range::all()) = 0.0;
	if (basis::tri(log2p)->sm() > 0) {
		hp_gbl->sprcn(Range(0,nseg-1),Range::all()) = 0.0;
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
			err = 1;
			break;
		}
		h = 4.*jcb/(0.25*(basis::tri(log2p)->p() +1)*(basis::tri(log2p)->p()+1)*hmax);
		hmaxp = hmax/(0.25*(basis::tri(log2p)->p() +1)*(basis::tri(log2p)->p()+1));

		/* ELASTICITY */
		gami = pow(hmaxp/(2.*hp_ps_gbl->mu),2.0) +hp_ps_gbl->lami*hmaxp*hmaxp/hp_ps_gbl->mu;
		gam = 1./gami;
		lam1 = sqrt(gam);
		hp_ps_gbl->tau(tind)  = adis*h/(jcb*lam1);

		/* TEMPERATURE */
		qmax = 0.0;
		for(j=0;j<3;++j) {
			v0 = v(j);
			mvel(0) = gbl->bd(0)*(pnts(v0)(0) -vrtxbd(1)(v0)(0));
			mvel(1) = gbl->bd(0)*(pnts(v0)(1) -vrtxbd(1)(v0)(1));
#ifdef MESH_REF_VEL
			mvel += hp_gbl->mesh_ref_vel;
#endif
			q = pow(hp_te_gbl->ax -mvel(0),2.0) +pow(hp_te_gbl->ay -mvel(1),2.0);
			qmax = MAX(qmax,q);
		}
		q = sqrt(qmax);
		lamT = (q +1.5*alpha/h +h*gbl->bd(0));
		hp_te_gbl->tau(tind) = (lamT > 0.0 ? adis*h/(jcb*lamT) : 0.0);
		/* PSEUDO-TIME STEP: DIFFUSION ESTIMATED THE SAME WAY AS THE ELASTICITY BELOW */
		/* (cd's 1.5*alpha/h ESTIMATE IS UNSTABLE AT THE DEFAULT CFL FOR DIFFUSION DOMINATED PROBLEMS) */
		jcbT = jcb*(hp_te_gbl->rhocv*(q +h*gbl->bd(0))/h +8.*hp_te_gbl->kcond*(1./(hmaxp*hmaxp) +1./(h*h)));

		jcb *= (8.*hp_ps_gbl->mu*(1./(hmaxp*hmaxp) +1./(h*h)));

		FLT rad = RAD((pnts(v(0))(0) +pnts(v(1))(0) +pnts(v(2))(0))/3.);
		hp_gbl->tprcn(tind,0) = jcb;
		hp_gbl->tprcn(tind,1) = jcb;
		hp_gbl->tprcn(tind,T) = jcbT*rad;
		hp_gbl->tprcn(tind,P) = jcb/gam;
		for(i=0;i<3;++i) {
			hp_gbl->vprcn(v(i),Range::all())  += hp_gbl->tprcn(tind,Range::all());
			if (basis::tri(log2p)->sm() > 0) {
				side = tri(tind).seg(i);
				hp_gbl->sprcn(side,Range::all()) += hp_gbl->tprcn(tind,Range::all());
			}
		}
	}
	return(tri_hp::setup_preconditioner()+err);
}
