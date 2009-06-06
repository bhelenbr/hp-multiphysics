#include <math.h>

#include "tri_hp_ps.h"
#include "../hp_boundary.h"

void tri_hp_ps::setup_preconditioner() {
    int tind,i,j,side;
    FLT jcb,h,hmax,lam1,gam,gami;
    TinyVector<int,3> v;

    /***************************************/
    /** DETERMINE PSEUDO-TIME STEP ****/
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

		gami = pow(hmax/(2.*gbl->mu),2.0) +gbl->lami*hmax*hmax/gbl->mu;
		gam = 1./gami; 
		lam1 = sqrt(gam);

		/* SET UP DISSIPATIVE COEFFICIENTS */
		gbl->tau(tind)  = adis*h/(jcb*lam1);
		jcb *= (8.*gbl->mu*(1./(hmax*hmax) +1./(h*h))) ;

		gbl->tprcn(tind,0) = jcb;    
		gbl->tprcn(tind,1) = jcb;      
		gbl->tprcn(tind,2) =  jcb/gam;
		for(i=0;i<3;++i) {
			gbl->vprcn(v(i),Range::all())  += gbl->tprcn(tind,Range::all());
			if (basis::tri(log2p).sm > 0) {
				side = tri(tind).seg(i);
				gbl->sprcn(side,Range::all()) += gbl->tprcn(tind,Range::all());
			}
		}
    }
    tri_hp::setup_preconditioner();

    return;
}
