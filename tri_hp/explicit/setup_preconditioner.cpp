#include "tri_hp_explicit.h"
#include <math.h>
#include "../hp_boundary.h"

// #define TESTING

int tri_hp_explicit::setup_preconditioner() {
	int side;
    FLT jcb,dtstari = 0.0;
	TinyVector<int,3> v;
	TinyVector<FLT,2> mvel;
	const int ltm = basis::tri(log2p)->tm();
	const int lsm = basis::tri(log2p)->sm();
	const int lbm = basis::tri(log2p)->bm();
	const int lim = basis::tri(log2p)->im();
	Array<FLT,2> mass(ltm,ltm);
	Array<FLT,1> uht(ltm), lf(ltm);
	Range all(0,ltm-1);
	const FLT alpha = hp_cd_gbl->kcond/hp_cd_gbl->rhocv;
    int err = 0;
    
	hp_gbl->vprcn(Range(0,npnt-1),Range::all()) = 0.0;
	hp_explicit_gbl->sprcn2(Range(0,nseg-1),Range::all(),Range::all()) = 0.0;

	
	/* Determine appropriate time step */
	for(int tind = 0; tind < ntri; ++tind) {
		jcb = 0.25*area(tind);  // area is 2 x triangle area
		v = tri(tind).pnt;
		FLT hmax = 0.0;
		for(int j=0;j<3;++j) {
			FLT h = pow(pnts(v(j))(0) -pnts(v((j+1)%3))(0),2.0) +
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
		FLT h = 4.*jcb/(0.25*(basis::tri(log2p)->p() +1)*(basis::tri(log2p)->p()+1)*hmax);
		hmax = hmax/(0.25*(basis::tri(log2p)->p() +1)*(basis::tri(log2p)->p()+1));
		
		FLT qmax = 0.0;
		for(int j=0;j<3;++j) {
			int v0 = v(j);
			
			mvel(0) = gbl->bd(0)*(pnts(v0)(0) -vrtxbd(1)(v0)(0));
			mvel(1) = gbl->bd(0)*(pnts(v0)(1) -vrtxbd(1)(v0)(1));
#ifdef MESH_REF_VEL
			mvel += hp_gbl->mesh_ref_vel;
#endif
			
#ifdef CONST_A
			FLT q = pow(hp_cd_gbl->ax -mvel(0),2.0)
			+pow(hp_cd_gbl->ay -mvel(1),2.0);
#else
			FLT q = pow(hp_cd_gbl->a->f(0,pnts(v0),gbl->time) -mvel(0),2.0)
			+pow(hp_cd_gbl->a->f(1,pnts(v0),gbl->time) -mvel(1)),2.0);
#endif
			
			qmax = MAX(qmax,q);
		}
		FLT q = sqrt(qmax);
		
		FLT lam1  = (q +1.5*alpha/h +h*hp_explicit_gbl->sigma);
		
		/* SET UP DISSIPATIVE COEFFICIENTS */
		hp_cd_gbl->tau(tind)  = adis*h/(jcb*lam1);

		/* DETERMINE MAX TIME STEP */
		dtstari = MAX(lam1/h,dtstari);
	}
	*gbl->log << "#iterative to physical time step ratio: " << gbl->bd(0) << ' ' << dtstari << '\n';

	uht(all) = 0.0;
	/*	SET TIME STEP TO BE 1 */
	for(int tind = 0; tind < ntri; ++tind) {
		jcb = 0.25*area(tind);
		v = tri(tind).pnt;

		/* CALCULATE MIXED BASIS MASS MATRIX ON THIS TRIANGLE */
		if (tri(tind).info > -1) {
			/* IF TINFO > -1 IT IS CURVED ELEMENT */
			/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
			crdtocht(tind);

			/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
			for(int n=0;n<ND;++n)
				basis::tri(log2p)->proj_bdry(&cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), &dcrd(n,1)(0,0),MXGP);

			for(int i=0;i<basis::tri(log2p)->gpx();++i) {
				for(int j=0;j<basis::tri(log2p)->gpn();++j) {
					cjcb(i,j) = RAD(crd(0)(i,j))*(dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j));
				}
			}

			Array<FLT,2>& lmass = hp_explicit_gbl->mass[tind];
			/* JUST GOING TO BE BRUTE FORCE FOR NOW */
			for (int m=0;m<ltm;++m) {
				uht(m) = 1.0;
				basis::tri(log2p)->proj(&uht(0), &u(0)(0,0),MXGP);
				for(int i=0;i<basis::tri(log2p)->gpx();++i) {
					for(int j=0;j<basis::tri(log2p)->gpn();++j) {
						res(0)(i,j) = u(0)(i,j)*cjcb(i,j)*dtstari;
					}
				}

				lf(all) = 0.0;
				basis::tri(log2p)->intgrt(&lf(0),&res(0)(0,0),MXGP);

				lmass(m,all) = lf(all);
				uht(m) = 0.0;
			}

			/* NOW REARRANGE ROWS TO MAKE PSI BASIS */
			/* SIDES FROM VERTICES */
			int indx = 3;
			for (int sind=0;sind<3;++sind) {
				for (int k=0; k <lsm; ++k) {
					for (int i=0; i<2; ++i) {
						lmass((sind+1+i)%3,all) -= basis::tri(log2p)->sfmv(i,k)*lmass(indx,all);
					}
					++indx;
				}
			}

			/* INTERIORS FROM SIDE & VERTEX */
			indx = 0;
			int indx2 = 3;
			for (int i=0; i<3; ++i) {
				for (int k=0;k<lim;++k)
					lmass(i,all) -= basis::tri(log2p)->ifmb(i,k)*lmass(k+lbm,all);

				for (int j=0;j<lsm;++j) {
					for (int k=0;k<lim;++k)
						lmass(indx2,all) -= basis::tri(log2p)->ifmb(indx2,k)*lmass(k+lbm,all);
					++indx2;
				}
			}

#ifdef TESTING
			lmass = lmass/jcb;
			*gbl->log << lmass << std::endl;
			*gbl->log << "vdiag " << basis::tri(log2p)->vdiag() << std::endl;
			for (int m=0;m<lsm;++m)
				*gbl->log << "sdiag " << basis::tri(log2p)->sdiag(m) << std::endl;
			for (int m=0; m<lsm;++m)
				*gbl->log << "vfms " << basis::tri(log2p)->vfms(0,m) << ' ' << basis::tri(log2p)->vfms(1,m) << ' ' << basis::tri(log2p)->vfms(2,m) << std::endl;
			for (int m=0; m<lsm-1;++m)
				for (int m1=m+1;m1 <lsm; ++m1)
					*gbl->log << "sfms " << basis::tri(log2p)->sfms(m,m1,0) << ' ' << basis::tri(log2p)->sfms(m,m1,1) << ' ' << basis::tri(log2p)->sfms(m,m1,2)<< std::endl;
			for(int mb=0;mb<lbm;++mb)
				for(int mi=0;mi<lim;++mi)
					*gbl->log << "bfmi " << basis::tri(log2p)->bfmi(mb,mi) << std::endl;
			for(int mi=0;mi<lim;++mi)
				*gbl->log << "idiag " << basis::tri(log2p)->idiag(mi,0) << std::endl;
            err = 1;
#endif

			indx = 3;
			for(int i=0;i<3;++i) {
				hp_gbl->vprcn(v(i),Range::all())  += lmass(i,0) +lmass(i,1) +lmass(i,2);
				if (basis::tri(log2p)->sm() > 0) {
					side = tri(tind).seg(i);
					for (int m=0;m<lsm;++m) {
						hp_explicit_gbl->sprcn2(side,m,Range::all()) += lmass(indx,indx);
						++indx;
					}
				}
			}
			hp_gbl->tprcn(tind,Range::all()) = gbl->bd(0)*jcb*RAD((pnts(v(0))(0) +pnts(v(1))(0) +pnts(v(2))(0))/3.);
		}
		else {
			hp_gbl->tprcn(tind,Range::all()) = gbl->bd(0)*jcb*RAD((pnts(v(0))(0) +pnts(v(1))(0) +pnts(v(2))(0))/3.);

			for(int i=0;i<3;++i) {
				hp_gbl->vprcn(v(i),Range::all())  += hp_gbl->tprcn(tind,Range::all())*basis::tri(log2p)->vdiag();
				if (basis::tri(log2p)->sm() > 0) {
					side = tri(tind).seg(i);
					for (int m=0;m<lsm;++m) {
						hp_explicit_gbl->sprcn2(side,m,Range::all()) += hp_gbl->tprcn(tind,Range::all())/basis::tri(log2p)->sdiag(m);
					}
				}
			}
		}
	}
	
	err += tri_hp::setup_preconditioner();
	
	hp_gbl->vprcn(Range(0,npnt-1),Range::all()) *= basis::tri(log2p)->vdiag();  // undoes vdiag term done by tri_hp::setup_preconditioner

	if (log2p) {
		sc0load(hp_explicit_gbl->sprcn2.data(),0,lsm*NV-1,lsm*NV);
		smsgpass(boundary::all,0,boundary::symmetric);
		sc0wait_rcv(hp_explicit_gbl->sprcn2.data(),0,lsm*NV-1,lsm*NV);   
		/* INVERT DIAGANOL PRECONDITIONER FOR SIDES */                
		hp_explicit_gbl->sprcn2(Range(0,nseg-1),Range::all(),Range::all()) = 1.0/hp_explicit_gbl->sprcn2(Range(0,nseg-1),Range::all(),Range::all());
	}

	return(err); 
}
