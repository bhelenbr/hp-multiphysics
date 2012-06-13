#include "tri_hp_explicit.h"
#include <math.h>
#include <utilities.h>
#include "../hp_boundary.h"

void tri_hp_explicit::setup_preconditioner() {
	int side;
	FLT jcb,dtstari;
	TinyVector<int,3> v;
	int ltm = basis::tri(log2p)->tm();
	int lsm = basis::tri(log2p)->sm();
	int lbm = basis::tri(log2p)->bm();
	int lim = basis::tri(log2p)->im();
	Array<FLT,2> mass(ltm,ltm);
	Array<FLT,1> uht(ltm), lf(ltm);
	Range all(0,ltm-1);

	gbl->vprcn(Range(0,npnt-1),Range::all()) = 0.0;
	if (basis::tri(log2p)->sm() > 0) {
		gbl->sprcn2(Range(0,nseg-1),Range::all(),Range::all()) = 0.0;
	}
	uht(all) = 0.0;


	/*	SET TIME STEP TO BE 1 */
	for(int tind = 0; tind < ntri; ++tind) {
		jcb = 0.25*area(tind);
		v = tri(tind).pnt;

		/* SET UP DIAGONAL PRECONDITIONER */
		dtstari = jcb*RAD((pnts(v(0))(0) +pnts(v(1))(0) +pnts(v(2))(0))/3.);

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
					cjcb(i,j) = dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j);
				}
			}

			Array<FLT,2>& lmass = gbl->mass[tind];
			/* JUST GOING TO BE BRUTE FORCE FOR NOW */
			for (int m=0;m<ltm;++m) {
				uht(m) = 1.0;
				basis::tri(log2p)->proj(&uht(0), &u(0)(0,0),MXGP);
				for(int i=0;i<basis::tri(log2p)->gpx();++i) {
					for(int j=0;j<basis::tri(log2p)->gpn();++j) {
						res(0)(i,j) = u(0)(i,j)*cjcb(i,j);
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

//			lmass = lmass/jcb;
//			std::cout << lmass << std::endl;
//			
//			std::cout << basis::tri(log2p)->vdiag() << std::endl;
//			
//			std::cout <<  basis::tri(log2p)->sdiag << std::endl;
//			
//			std::cout << basis::tri(log2p)->vfms << std::endl;
//
//			std::cout << basis::tri(log2p)->sfms << std::endl;
//			
//			std::cout << basis::tri(log2p)->bfmi << std::endl;
//			
//			std::cout << basis::tri(log2p)->idiag << std::endl;
//			
//			sim::abort(_);

			indx = 3;
			for(int i=0;i<3;++i) {
				gbl->vprcn(v(i),Range::all())  += lmass(i,0) +lmass(i,1) +lmass(i,2);
				if (basis::tri(log2p)->sm() > 0) {
					side = tri(tind).seg(i);
					for (int m=0;m<lsm;++m) {
						gbl->sprcn2(side,m,Range::all()) += lmass(indx,indx);
						++indx;
					}
				}
			}

			gbl->tprcn(tind,Range::all()) = dtstari;
		}
		else {
			gbl->tprcn(tind,Range::all()) = dtstari;        

			for(int i=0;i<3;++i) {
				gbl->vprcn(v(i),Range::all())  += gbl->tprcn(tind,Range::all())*basis::tri(log2p)->vdiag();
				if (basis::tri(log2p)->sm() > 0) {
					side = tri(tind).seg(i);
					for (int m=0;m<lsm;++m) {
						gbl->sprcn2(side,m,Range::all()) += gbl->tprcn(tind,Range::all())/basis::tri(log2p)->sdiag(m);
					}
				}
			}
		}
	}

	int i,last_phase,mp_phase;

	/* SET UP TSTEP FOR MESH MOVEMENT */
	if (mmovement == coupled_deformable && log2p == 0) {
		r_tri_mesh::setup_preconditioner();    
	}

	/* SET UP TSTEP FOR ACTIVE BOUNDARIES */
	for(i=0;i<nebd;++i)
		hp_ebdry(i)->setup_preconditioner();

	/* SET UP TSTEP FOR HELPER */
	helper->setup_preconditioner();    

	if (gbl->diagonal_preconditioner) {
		for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
			vc0load(mp_phase,gbl->vprcn.data());
			pmsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
			last_phase = true;
			last_phase &= vc0wait_rcv(mp_phase,gbl->vprcn.data());
		}
		/* PREINVERT PRECONDITIONER FOR VERTICES */
		gbl->vprcn(Range(0,npnt-1),Range::all()) = 1.0/gbl->vprcn(Range(0,npnt-1),Range::all());

		if (log2p) {
			sc0load(gbl->sprcn2.data(),0,lsm*NV-1,lsm*NV);
			smsgpass(boundary::all,0,boundary::symmetric);
			sc0wait_rcv(gbl->sprcn2.data(),0,lsm*NV-1,lsm*NV);   
			/* INVERT DIAGANOL PRECONDITIONER FOR SIDES */                
			gbl->sprcn2(Range(0,nseg-1),Range::all(),Range::all()) = 1.0/gbl->sprcn2(Range(0,nseg-1),Range::all(),Range::all());
		}

	}
	else {
		/* NEED STUFF HERE FOR CONTINUITY OF MATRIX PRECONDITIONER */
		for(int stage = 0; stage<NV; ++stage) {
			for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
				vc0load(mp_phase,gbl->vprcn_ut.data() +stage*NV,NV);
				pmsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
				last_phase = true;
				last_phase &= vc0wait_rcv(mp_phase,gbl->vprcn_ut.data()+stage*NV,NV);
			}
			if (log2p) {
				sc0load(gbl->sprcn_ut.data()+stage*NV,0,0,NV);
				smsgpass(boundary::all,0,boundary::symmetric);
				sc0wait_rcv(gbl->sprcn_ut.data()+stage*NV,0,0,NV);
			}
		}

		/* FACTORIZE PRECONDITIONER FOR VERTICES ASSUMES LOWER TRIANGULAR NOTHING  */
//        for(i=0;i<npnt;++i)
//            for(int n=0;n<NV;++n)
//                gbl->vprcn_ut(i,n,n) = 1.0/(basis::tri(log2p)->vdiag()*gbl->vprcn_ut(i,n,n));
//      
//        if (basis::tri(log2p)->sm() > 0) {
//            /* INVERT DIAGANOL PRECONDITIONER FOR SIDES ASSUMES LOWER TRIANGULAR */     
//            for(i=0;i<nseg;++i)
//                for(int n=0;n<NV;++n)
//                    gbl->sprcn_ut(i,n,n)= 1.0/gbl->sprcn_ut(i,n,n);
//        }
	}


	return; 
}
