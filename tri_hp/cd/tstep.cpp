#include "tri_hp_cd.h"
#include <math.h>
#include <utilities.h>
#include "../hp_boundary.h"


void tri_hp_cd::setup_preconditioner() {
    int tind,i,j,side,v0;
    FLT jcb,h,hmax,q,qmax,lam1;
    TinyVector<int,3> v;


    /***************************************/
    /** DETERMINE FLOW PSEUDO-TIME STEP ****/
    /***************************************/
    gbl->vprcn(Range(0,npnt-1),Range::all()) = 0.0;
    if (basis::tri(log2p).sm > 0) {
		gbl->sprcn(Range(0,nseg-1),Range::all()) = 0.0;
    }

#ifdef TIMEACCURATE
    FLT dtstari = 0.0;
#endif

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
		for(j=0;j<3;++j) {
			v0 = v(j);
			q = pow(gbl->ax -(gbl->bd(0)*(pnts(v0)(0) -vrtxbd(1)(v0)(0))),2.0) 
				+pow(gbl->ay -(gbl->bd(0)*(pnts(v0)(1) -vrtxbd(1)(v0)(1))),2.0);
			qmax = MAX(qmax,q);
		}
		q = sqrt(qmax);

		lam1  = (q +1.5*gbl->nu/h +h*gbl->bd(0));

		/* SET UP DISSIPATIVE COEFFICIENTS */
		gbl->tau(tind)  = adis*h/(jcb*lam1);

		jcb *= lam1/h;


		/* SET UP DIAGONAL PRECONDITIONER */
#ifdef TIMEACCURATE
		dtstari = MAX(lam1/h,dtstari);
    }
    printf("#iterative to physical time step ratio: %f\n",gbl->bd(0)/dtstari);

    for(tind=0;tind<ntri;++tind) {
		v = tri(tind).pnt;
		jcb = 0.25*area(tind)*dtstari;
#endif
		jcb *= RAD((pnts(v(0))(0) +pnts(v(1))(0) +pnts(v(2))(0))/3.);
		gbl->tprcn(tind,0) = jcb;    
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
