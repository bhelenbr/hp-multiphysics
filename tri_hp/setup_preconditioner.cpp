#include "tri_hp.h"
#include "hp_boundary.h"
#include <myblas.h>


/* This is for a diagonal or block diagonal preconditioner where the block could be of size NVxNV */
/* Physics classes should call this routine after they have finished setting up their diagonal preconditioner */

void tri_hp::setup_preconditioner() {
	int i,last_phase,mp_phase;
	
#ifndef petsc
	/* SET UP TSTEP FOR MESH MOVEMENT */
	if (mmovement == coupled_deformable && log2p == 0) {
		r_tri_mesh::setup_preconditioner();
	}
#else
	if (mmovement == coupled_deformable) {
		r_tri_mesh::setup_preconditioner();
	}
#endif
	
	/* Do communication for preconditioner */
	if (gbl->diagonal_preconditioner) {
		for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
			vc0load(mp_phase,gbl->vprcn.data());
			pmsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
			last_phase = true;
			last_phase &= vc0wait_rcv(mp_phase,gbl->vprcn.data());
		}
		
		if (log2p) {
			sc0load(gbl->sprcn.data(),0,0,1);
			smsgpass(boundary::all,0,boundary::symmetric);
			sc0wait_rcv(gbl->sprcn.data(),0,0,1);
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
	}
	
	/* Modify according to boundary conditions */
	for(i=0;i<nebd;++i)
		hp_ebdry(i)->setup_preconditioner();
	
	/* SET UP TSTEP FOR ACTIVE BOUNDARIES */
	for(i=0;i<nvbd;++i)
		hp_vbdry(i)->setup_preconditioner();
	
	/* SET UP TSTEP FOR HELPER */
	helper->setup_preconditioner();
	
	/* Invert Preconditioner */
	if (gbl->diagonal_preconditioner) {
		/* PREINVERT PRECONDITIONER FOR VERTICES */
		gbl->vprcn(Range(0,npnt-1),Range::all()) = 1.0/(basis::tri(log2p)->vdiag()*gbl->vprcn(Range(0,npnt-1),Range::all()));
		if (log2p) {
			/* INVERT DIAGANOL PRECONDITIONER FOR SIDES */
			gbl->sprcn(Range(0,nseg-1),Range::all()) = 1.0/gbl->sprcn(Range(0,nseg-1),Range::all());
		}
		
	}
	else {
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
	
#ifdef petsc
	petsc_setup_preconditioner();
#endif
	
	return;
}


/************************************************/
/**********        INVERT MASS MATRIX     **********/
/************************************************/
void tri_hp::minvrt() {
	int i,j,k,m,n,tind,sind,v0,indx,indx1,indx2,sgn,msgn;
	TinyVector<int,3> sign,side;
	Array<FLT,2> tinv(NV,NV);
	Array<FLT,1> temp(NV);
	int last_phase, mp_phase;

	/* LOOP THROUGH SIDES */
	if (basis::tri(log2p)->sm() > 0) {
		indx = 0;
		for(sind = 0; sind<nseg;++sind) {
			/* SUBTRACT SIDE CONTRIBUTIONS TO VERTICES */            
			for (k=0; k <basis::tri(log2p)->sm(); ++k) {
				for (i=0; i<2; ++i) {
					v0 = seg(sind).pnt(i);
					for(n=0;n<NV;++n)
						gbl->res.v(v0,n) -= basis::tri(log2p)->sfmv(i,k)*gbl->res.s(sind,k,n);
				}
				++indx;
			}
		}

		if (basis::tri(log2p)->im() > 0) {
			/* SUBTRACT INTERIORS */
			indx = 0;
			for(tind = 0; tind<ntri;++tind) {
				indx2 = 3;
				for (i=0; i<3; ++i) {
					v0 = tri(tind).pnt(i);
					for (k=0;k<basis::tri(log2p)->im();++k)
						for(n=0;n<NV;++n)
							gbl->res.v(v0,n) -= basis::tri(log2p)->ifmb(i,k)*gbl->res.i(tind,k,n);

					sind = tri(tind).seg(i);
					sgn = tri(tind).sgn(i);
					msgn = 1;
					for (j=0;j<basis::tri(log2p)->sm();++j) {
						for (k=0;k<basis::tri(log2p)->im();++k)
							for(n=0;n<NV;++n)
								gbl->res.s(sind,j,n) -= msgn*basis::tri(log2p)->ifmb(indx2,k)*gbl->res.i(tind,k,n);
						msgn *= sgn;
						++indx2;
					}
				}
				indx += basis::tri(log2p)->im();
			}
		}
	}

	if (gbl->diagonal_preconditioner) {
		gbl->res.v(Range(0,npnt-1),Range::all()) *= gbl->vprcn(Range(0,npnt-1),Range::all());
    } else {
		/* ASSUMES LOWER TRIANGULAR FOR NOW */
		for(i=0;i<npnt;++i) {
			/* forward substitution */
			for(j=0;j<NV;++j){
				FLT lcl = gbl->res.v(i,j);
				for(k=0;k<j;++k){
					lcl -= gbl->vprcn_ut(i,j,k)*gbl->res.v(i,k);
				}
				gbl->res.v(i,j) = lcl/gbl->vprcn_ut(i,j,j);
			}
//			DGETLS(&gbl->vprcn_ut(i,0,0), NV, NV, &gbl->res.v(i,0));
			
			gbl->res.v(i,Range::all()) /= basis::tri(log2p)->vdiag();

		}
	}

	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		vc0load(mp_phase,gbl->res.v.data());
		pmsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
		last_phase = true;
		last_phase &= vc0wait_rcv(mp_phase,gbl->res.v.data());
	}

	/* APPLY VERTEX DIRICHLET B.C.'S */
	for(i=0;i<nebd;++i)
		hp_ebdry(i)->vdirichlet();

	for(i=0;i<nvbd;++i)
		hp_vbdry(i)->vdirichlet();

	if(basis::tri(log2p)->sm() == 0) return;


	/* REMOVE VERTEX CONTRIBUTION FROM SIDE MODES */
	/* SOLVE FOR SIDE MODES */
	/* PART 1 REMOVE VERTEX CONTRIBUTIONS */

	if (gbl->diagonal_preconditioner) {
		for(tind=0;tind<ntri;++tind) {
			for(i=0;i<3;++i) {
				v0 = tri(tind).pnt(i);
				for(n=0;n<NV;++n)
					uht(n)(i) = gbl->res.v(v0,n)*gbl->tprcn(tind,n);
			}

			for(i=0;i<3;++i) {
				sind = tri(tind).seg(i);
				sgn  = tri(tind).sgn(i);
				for(j=0;j<3;++j) {
					indx1 = (i+j)%3;
					msgn = 1;
					for(k=0;k<basis::tri(log2p)->sm();++k) {
						for(n=0;n<NV;++n)
							gbl->res.s(sind,k,n) -= msgn*basis::tri(log2p)->vfms(j,k)*uht(n)(indx1);
						msgn *= sgn;
					}
				}
			}
		}
	}
	else {
		/* THIS IS TO USE A MATRIX PRECONDITIONER */
		for(tind=0;tind<ntri;++tind) {
			for(i=0;i<3;++i) {
				v0 = tri(tind).pnt(i);
				for(n=0;n<NV;++n)
					uht(n)(i) = gbl->res.v(v0,n);
			}

			for(i=0;i<3;++i) {
				indx = tri(tind).seg(i);
				sgn  = tri(tind).sgn(i);
				for(j=0;j<3;++j) {
					indx1 = (i+j)%3;
					msgn = 1;
					for(k=0;k<basis::tri(log2p)->sm();++k) {
						for(n=0;n<NV;++n) {
							for(m=0;m<NV;++m) {
								gbl->res.s(indx,k,n) -= msgn*basis::tri(log2p)->vfms(j,k)*gbl->tprcn_ut(tind,n,m)*uht(m)(indx1);
							}
						}
						msgn *= sgn;
					}
				}
			}
		}
	} 


	for (int mode = 0; mode < basis::tri(log2p)->sm()-1; ++ mode) {
		/* SOLVE FOR SIDE MODE */
		if (gbl->diagonal_preconditioner) {
			gbl->res.s(Range(0,nseg-1),mode,Range::all()) *= gbl->sprcn(Range(0,nseg-1),Range::all())*basis::tri(log2p)->sdiag(mode);
		}
		else {
			for(sind = 0; sind < nseg; ++sind) {
				DGETLS(&gbl->sprcn_ut(sind,0,0), NV, NV, &gbl->res.s(sind,mode,0));
				gbl->res.s(sind,mode,Range::all()) /= basis::tri(log2p)->sdiag(mode);
			}
		}

		sc0load(gbl->res.s.data(),mode,mode,gbl->res.s.extent(secondDim));
		smsgpass(boundary::all,0,boundary::symmetric);
		sc0wait_rcv(gbl->res.s.data(),mode,mode,gbl->res.s.extent(secondDim));


		/* APPLY DIRCHLET B.C.S TO MODE */
		for(i=0;i<nebd;++i)
			hp_ebdry(i)->sdirichlet(mode);


		/* REMOVE MODE FROM HIGHER MODES */
		for(tind=0;tind<ntri;++tind) {

			if (gbl->diagonal_preconditioner) {
				for(i=0;i<3;++i) {
					side(i) = tri(tind).seg(i);
					sign(i) = tri(tind).sgn(i);
					sgn      = (mode % 2 ? sign(i) : 1);
					for(n=0;n<NV;++n)
						uht(n)(i) = sgn*gbl->res.s(side(i),mode,n)*gbl->tprcn(tind,n);
				}

				/* REMOVE MODES J,K FROM MODE I,M */
				for(i=0;i<3;++i) {
					msgn = ( (mode +1) % 2 ? sign(i) : 1);
					for(m=mode+1;m<basis::tri(log2p)->sm();++m) {
						for(j=0;j<3;++j) {
							indx = (i+j)%3;
							for(n=0;n<NV;++n) {
								gbl->res.s(side(i),m,n) -= msgn*basis::tri(log2p)->sfms(mode,m,j)*uht(n)(indx);
							}
						}
						msgn *= sign(i);
					}
				}
			}
			else {
				for(i=0;i<3;++i) {
					side(i) = tri(tind).seg(i);
					sign(i) = tri(tind).sgn(i);
					sgn      = (mode % 2 ? sign(i) : 1);
					for(n=0;n<NV;++n)
						uht(n)(i) = sgn*gbl->res.s(side(i),mode,n);
				}

				/* REMOVE MODES J,K FROM MODE I,M */
				for(i=0;i<3;++i) {
					msgn = ( (mode +1) % 2 ? sign(i) : 1);
					for(m=mode+1;m<basis::tri(log2p)->sm();++m) {
						for(j=0;j<3;++j) {
							indx = (i+j)%3;
							for(n=0;n<NV;++n) {
								gbl->res.s(side(i),m,n) -= msgn*basis::tri(log2p)->sfms(mode,m,j)*gbl->tprcn_ut(tind,n,0)*uht(0)(indx);
								for(k=1;k<NV;++k) {
									gbl->res.s(side(i),m,n) -= msgn*basis::tri(log2p)->sfms(mode,m,j)*gbl->tprcn_ut(tind,n,k)*uht(k)(indx);
								}
							}
						}
						msgn *= sign(i);
					}
				}
			}
		}
	}
	/* SOLVE FOR HIGHEST MODE */
	int mode = basis::tri(log2p)->sm()-1;
	if (gbl->diagonal_preconditioner) {
		gbl->res.s(Range(0,nseg-1),mode,Range::all()) *= gbl->sprcn(Range(0,nseg-1),Range::all())*basis::tri(log2p)->sdiag(mode);
	}
	else {
		for(sind = 0; sind < nseg; ++sind) {
			DGETLS(&gbl->sprcn_ut(sind,0,0), NV, NV, &gbl->res.s(sind,mode,0));
			gbl->res.s(sind,mode,Range::all()) /= basis::tri(log2p)->sdiag(mode);
		}
	}

	sc0load(gbl->res.s.data(),mode,mode,gbl->res.s.extent(secondDim));
	smsgpass(boundary::all,0,boundary::symmetric);
	sc0wait_rcv(gbl->res.s.data(),mode,mode,gbl->res.s.extent(secondDim));

	/* APPLY DIRCHLET B.C.S TO MODE */
	for(i=0;i<nebd;++i)
		hp_ebdry(i)->sdirichlet(mode);


    if (basis::tri(log2p)->im() == 0) return;

    /* SOLVE FOR INTERIOR MODES */
    for(tind = 0; tind < ntri; ++tind) {
		DPBTRSNU2((double *) &basis::tri(log2p)->idiag(0,0),basis::tri(log2p)->ibwth()+1,basis::tri(log2p)->im(),basis::tri(log2p)->ibwth(),&(gbl->res.i(tind,0,0)),NV);
		restouht_bdry(tind);
		if (gbl->diagonal_preconditioner) {
			for(k=0;k<basis::tri(log2p)->im();++k) {
				gbl->res.i(tind,k,Range::all()) /= gbl->tprcn(tind,Range::all());

				for (i=0;i<basis::tri(log2p)->bm();++i)
					for(n=0;n<NV;++n) 
						gbl->res.i(tind,k,n) -= basis::tri(log2p)->bfmi(i,k)*uht(n)(i);
			}
		}
		else {
			for(k=0;k<basis::tri(log2p)->im();++k) {
				/* SUBTRACT BOUNDARY MODES (bfmi is multipled by interior inverse matrix so do this after DPBSLN) */
				for (i=0;i<basis::tri(log2p)->bm();++i) {
					for(n=0;n<NV;++n) {
						for(m=0;m<NV;++m) {
							gbl->res.i(tind,k,n) -= basis::tri(log2p)->bfmi(i,k)*uht(m)(i)*gbl->tprcn_ut(tind,n,m);
						}
					}
				}
				/* INVERT PRECONDITIONER (ASSUMES LOWER TRIANGULAR) */
				DGETLS(&gbl->tprcn_ut(tind,0,0), NV, NV, &gbl->res.i(tind,k,0));
			}
		}
	}

	return;
}

void tri_hp::restouht_bdry(int tind) {
	int i,m,n,indx,cnt;
	int sign, msgn;

	for (i=0; i<3; ++i) {
		indx = tri(tind).pnt(i);
		for(n=0; n<NV; ++n)
			uht(n)(i) = gbl->res.v(indx,n);
	}

	cnt = 3;
	for(i=0;i<3;++i) {
		indx = tri(tind).seg(i);
		sign = tri(tind).sgn(i);
		msgn = 1;
		for (m = 0; m < basis::tri(log2p)->sm(); ++m) {
			for(n=0; n<NV; ++n)
				uht(n)(cnt) = msgn*gbl->res.s(indx,m,n);
			msgn *= sign;
			++cnt;
		}
	}

	return;
}


