#include <math.h>
#include "tri_hp_swirl.h"
#include "../hp_boundary.h"

int tri_hp_swirl::setup_preconditioner() {
	int tind,i,j,side,v0;
	FLT jcb,h,hmax,q,qs,qmax,qsmax,lam1,gam;
	TinyVector<FLT,ND> mvel;
	TinyVector<int,3> v;
    int err = 0;
    
	FLT nu = gbl->mu/gbl->rho;

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
            err = 1;
            break;
		}
		h = 4.*jcb/(0.25*(basis::tri(log2p)->p() +1)*(basis::tri(log2p)->p()+1)*hmax);
		hmax = hmax/(0.25*(basis::tri(log2p)->p() +1)*(basis::tri(log2p)->p()+1));

		qmax = 0.0;
		qsmax = 0.0;
		for(j=0;j<3;++j) {
			v0 = v(j);

			mvel(0) = gbl->bd(0)*(pnts(v0)(0) -vrtxbd(1)(v0)(0));
			mvel(1) = gbl->bd(0)*(pnts(v0)(1) -vrtxbd(1)(v0)(1));
#ifdef MESH_REF_VEL
			mvel += gbl->mesh_ref_vel;
#endif
			
			q = pow(ug.v(v0,0)-0.5*mvel(0),2.0) +pow(ug.v(v0,1)-0.5*mvel(1),2.0);

			qs = q +pow(ug.v(v0,2),2.0); // FIXME?
			qmax = MAX(qmax,q);
			qsmax = MAX(qsmax,qs);
		}
        if  (std::isnan(qmax)) {
            *gbl->log << gbl->idprefix << ' ' << tind << std::endl;
            *gbl->log << "flow solution has nan's " << qmax << std::endl;
            err = 1;
        }
        
		gam = 3.0*qsmax +(0.5*hmax*gbl->bd(0) +2.*nu/hmax)*(0.5*hmax*gbl->bd(0) +2.*nu/hmax);
		if (gbl->mu + gbl->bd(0) == 0.0) gam = MAX(gam,0.01);
		qs = sqrt(qsmax);
		lam1 = qs + sqrt(qsmax +gam);

		/* SET UP DISSIPATIVE COEFFICIENTS */
		gbl->tau(tind,0) = adis*h/(jcb*sqrt(gam));
		gbl->tau(tind,NV-1) = qsmax*gbl->tau(tind,0);

		/* THIS IS FOR SCALARS */
		for (int n=3;n<NV-1;++n) {
			double diff = gbl->D(n)/gbl->rho;
			double dth = qmax +(0.5*hmax*gbl->bd(0) +2.*diff/hmax)*(0.5*hmax*gbl->bd(0) +2.*diff/hmax);
			double rtdth = sqrt(dth);
			gbl->tau(tind,n) = adis*h/(jcb*rtdth);
			gbl->tprcn(tind,n) = gbl->rho*rtdth/h*RAD((pnts(v(0))(0) +pnts(v(1))(0) +pnts(v(2))(0))/3.)*jcb;
		}

		/* SET UP DIAGONAL PRECONDITIONER */
		// jcb *= 8.*nu*(1./(hmax*hmax) +1./(h*h)) +2*lam1/h +2*sqrt(gam)/hmax +gbl->bd(0);
		jcb *= 2.*nu*(1./(hmax*hmax) +1./(h*h)) +3*lam1/h;  // heuristically tuned
		jcb *= RAD((pnts(v(0))(0) +pnts(v(1))(0) +pnts(v(2))(0))/3.);

		gbl->tprcn(tind,0) = gbl->rho*jcb;    
		gbl->tprcn(tind,1) = gbl->rho*jcb;  
		gbl->tprcn(tind,2) = gbl->rho*jcb;      
		gbl->tprcn(tind,NV-1) =  jcb/gam;
		for(i=0;i<3;++i) {
			gbl->vprcn(v(i),Range::all())  += gbl->tprcn(tind,Range::all());
			if (basis::tri(log2p)->sm() > 0) {
				side = tri(tind).seg(i);
				gbl->sprcn(side,Range::all()) += gbl->tprcn(tind,Range::all());
			}
		}
	}
	return(tri_hp::setup_preconditioner()+err);
}
