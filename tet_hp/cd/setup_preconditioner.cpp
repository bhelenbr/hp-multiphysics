#include "tet_hp_cd.h"
#include <math.h>
#include <utilities.h>
#include "../hp_boundary.h"

//#define TIMEACCURATE

void tet_hp_cd::setup_preconditioner() {
	int tind,find,i,j,side,p0,p1,p2,v0;
	FLT jcb,a,h,amax,lam1,q,qmax,amin,hmax,hmin,dtcheck;
	FLT dx1,dy1,dx2,dy2,dz1,dz2,cpi,cpj,cpk;
	TinyVector<int,4> v;

	/***************************************/
	/** DETERMINE FLOW PSEUDO-TIME STEP ****/
	/***************************************/
	gbl->vprcn(Range::all(),Range::all())=0.0;
	if (basis::tet(log2p).em > 0) {
		gbl->eprcn(Range::all(),Range::all())=0.0;
		if (basis::tet(log2p).fm > 0) {
			gbl->fprcn(Range::all(),Range::all())=0.0;
			if (basis::tet(log2p).im > 0) {
				gbl->iprcn(Range::all(),Range::all())=0.0;
			}
		}
	}
	
#ifdef TIMEACCURATE
	FLT dtstari = 0.0;
#endif
//	hmax = 0;
//	hmin = 1000000;
	dtcheck = 0.0;
	for(tind = 0; tind < ntet; ++tind) {
		jcb = 0.125*tet(tind).vol; 
		v = tet(tind).pnt;
		amax = 0.0;
		for(j=0;j<4;++j) { // FIND MAX FACE AREA AND THEN DIVIDE VOLUME BY IT 
			find = tet(tind).tri(j);
			p0 = tri(find).pnt(0);
			p1 = tri(find).pnt(1);
			p2 = tri(find).pnt(2);
	
			dx1 = pnts(p0)(0)-pnts(p1)(0);
			dy1 = pnts(p0)(1)-pnts(p1)(1);
			dz1 = pnts(p0)(2)-pnts(p1)(2);
			dx2 = pnts(p0)(0)-pnts(p2)(0);
			dy2 = pnts(p0)(1)-pnts(p2)(1);
			dz2 = pnts(p0)(2)-pnts(p2)(2);
			cpi = dy1*dz2-dz1*dy2;
			cpj = -dx1*dz2+dz1*dx2;
			cpk = dx1*dy2-dy1*dx2;
			a =	.5*sqrt(cpi*cpi+cpj*cpj+cpk*cpk);
			amax = (a > amax ? a : amax);
		}

		
		h = 4.0*jcb/(0.25*(basis::tet(log2p).p+1)*(basis::tet(log2p).p+1)*amax); // 3*8/6=4
			
		qmax = 0.0;
		for(j=0;j<4;++j) {
			v0 = v(j);
			q = pow(gbl->ax -(gbl->bd(0)*(pnts(v0)(0) -vrtxbd(1)(v0)(0))),2.0) 
				+pow(gbl->ay -(gbl->bd(0)*(pnts(v0)(1) -vrtxbd(1)(v0)(1))),2.0)
				+pow(gbl->az -(gbl->bd(0)*(pnts(v0)(2) -vrtxbd(1)(v0)(2))),2.0);
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
	
	*gbl->log << "#iterative to physical time step ratio: " << gbl->bd(0)/dtstari << ' ' << dtstari << endl;
		
	for(tind=0;tind<ntet;++tind) {
		v = tet(tind).pnt;
		jcb = 0.125*tet(tind).vol*dtstari;
#endif
		/* explicit */
		jcb = 0.125*tet(tind).vol*12000.0;
		
		gbl->iprcn(tind,0) = jcb; 
			
		for(i=0;i<4;++i) 
			gbl->vprcn(v(i),0) += gbl->iprcn(tind,0);
			
		if (basis::tet(log2p).em > 0) {
			for(i=0;i<6;++i){
				side = tet(tind).seg(i);
				gbl->eprcn(side,0) += gbl->iprcn(tind,0);
			}
		
			if (basis::tet(log2p).fm > 0) {
				for(i=0;i<4;++i){
					gbl->fprcn(tet(tind).tri(i),0) += gbl->iprcn(tind,0);
				}
			}
		}	

	}

	tet_hp::setup_preconditioner();
	
	return; 
}
