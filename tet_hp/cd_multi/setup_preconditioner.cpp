#include "tet_hp_cd_multi.h"
#include <math.h>
#include "../hp_boundary.h"

//#define TIMEACCURATE
#define NEWWAY

void tet_hp_cd_multi::setup_preconditioner() {
	int tind,find,i,j,side,p0,p1,p2,v0;
	FLT jcb,a,h,amax,lam1,q,qmax,dtcheck;
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
		
		const FLT rhocv = gbl->rhocv(marks(tind));
		const FLT kcond = gbl->kcond(marks(tind));
		const FLT alpha = kcond/rhocv;

		jcb = 0.125*tet(tind).vol;
		v = tet(tind).pnt;
		
#ifdef NEWWAY
		TinyVector<TinyVector<FLT,3>,3> d, dcrd;
		
		for(int n=0;n<ND;++n) {
			dcrd(n)(0) = 0.5*(pnts(v(3))(n) -pnts(v(2))(n));
			dcrd(n)(1) = 0.5*(pnts(v(1))(n) -pnts(v(2))(n));
			dcrd(n)(2) = 0.5*(pnts(v(0))(n) -pnts(v(2))(n));
		}
		
		d(0)(0) =  dcrd(1)(1)*dcrd(2)(2)-dcrd(1)(2)*dcrd(2)(1);
		d(0)(1) = -dcrd(0)(1)*dcrd(2)(2)+dcrd(0)(2)*dcrd(2)(1);
		d(0)(2) =  dcrd(0)(1)*dcrd(1)(2)-dcrd(0)(2)*dcrd(1)(1);
		d(1)(0) = -dcrd(1)(0)*dcrd(2)(2)+dcrd(1)(2)*dcrd(2)(0);
		d(1)(1) =  dcrd(0)(0)*dcrd(2)(2)-dcrd(0)(2)*dcrd(2)(0);
		d(1)(2) = -dcrd(0)(0)*dcrd(1)(2)+dcrd(0)(2)*dcrd(1)(0);
		d(2)(0) =  dcrd(1)(0)*dcrd(2)(1)-dcrd(1)(1)*dcrd(2)(0);
		d(2)(1) = -dcrd(0)(0)*dcrd(2)(1)+dcrd(0)(1)*dcrd(2)(0);
		d(2)(2) =  dcrd(0)(0)*dcrd(1)(1)-dcrd(0)(1)*dcrd(1)(0);
		
		
		double jcbi = kcond*0.25*(basis::tet(log2p).p+1)*(basis::tet(log2p).p+1)/jcb;
		FLT diff_dti0 = d(0)(0)*d(0)(0)+d(0)(1)*d(0)(1)+d(0)(2)*d(0)(2);
		FLT diff_dti1 =	d(1)(0)*d(1)(0)+d(1)(1)*d(1)(1)+d(1)(2)*d(1)(2);
		FLT diff_dti2 =	d(2)(0)*d(2)(0)+d(2)(1)*d(2)(1)+d(2)(2)*d(2)(2);
		/* This is to make matrix norm = 1 */
		jcb *= gbl->bd(0)*rhocv;
		jcb += 2*jcbi*MAX(MAX(diff_dti0,diff_dti1),diff_dti2);

		/* This is for testing. For p=1, the largest diagonal entry should match */ 
//		Array<FLT,2> K(basis::tet(log2p).tm,basis::tet(log2p).tm);
//		element_jacobian(tind,K);
//		*gbl->log << K << std::endl;
//		for(int k=0;k<basis::tet(log2p).tm;++k) {
//			FLT columnsum = 0; 
//			for(int l=0;l<basis::tet(log2p).tm;++l) {
//				columnsum += fabs(K(k,l));
//			}
//			*gbl->log << columnsum*basis::tet(log2p).vdiag/jcb << std::endl;
//		}
//		sim::finalize(__LINE__, __FILE__, gbl->log);
#else
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

		lam1  = (q +1.5*alpha/h +h*gbl->bd(0));

        /* SET UP DISSIPATIVE COEFFICIENTS */
		gbl->tau(tind)  = adis*h/(jcb*lam1);
		
		jcb *= rhocv*lam1/h;


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
		//jcb = 0.125*tet(tind).vol*12000.0;
#endif
		
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

void tet_hp_cd_multi::connect(multigrid_interface& in) {
	tet_hp_cd_multi &tgt = dynamic_cast<tet_hp_cd_multi&>(in);
	tet_hp_cd::connect(in);
	
	/* To fill in material marks for multigrid mesh */
	/* This just picks the matarial located at the center of the coarse element on the fine grid */
	TinyVector<FLT,ND> center;
	/* Fill in Marks in some way */
	for (int tind=0;tind<ntet;++tind) {
		center = 0.25*(pnts(tet(tind).pnt(0)) +pnts(tet(tind).pnt(1)) +pnts(tet(tind).pnt(2)) +pnts(tet(tind).pnt(3)));
		int ttgt = fcnnct(tet(tind).pnt(0)).tet;
		tgt.findtet(center,ttgt);
		marks(tind) = tgt.marks(ttgt);
	}
}





