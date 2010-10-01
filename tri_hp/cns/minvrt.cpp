#include "tri_hp_cns.h"
#include "../hp_boundary.h"
#include <myblas.h>

/************************************************/
/**********        INVERT MASS MATRIX     **********/
/************************************************/
void tri_hp_cns::minvrt() {
	int i,j,k,m,n,tind,sind,v0,indx,indx1,indx2,sgn,msgn;
	TinyVector<int,3> sign,side;
	Array<FLT,2> tinv(NV,NV);
	Array<FLT,1> temp(NV);
	int last_phase, mp_phase;

	if (gbl->diagonal_preconditioner) {
		tri_hp::minvrt();
		//gbl->res.v(Range(0,npnt-1),Range::all()) *= gbl->vprcn(Range(0,npnt-1),Range::all());
    } 
	else {
		/* only p = 1 for now */
		/* uses Pinv */
		for(i=0;i<npnt;++i) {
			
			/*  LU factorization  */
			int info,ipiv[NV];
			GETRF(NV, NV, &gbl->vprcn_ut(i,0,0), NV, ipiv, info);
			
			if (info != 0) {
				*gbl->log << "DGETRF FAILED IN CNS MINVRT" << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			
			/* SOLVE */
			char trans[] = "T";
			GETRS(trans, NV, 1, &gbl->vprcn_ut(i,0,0), NV, ipiv, &gbl->res.v(i,0), NV, info);
			
			if (info != 0) {
				*gbl->log << "DGETRS FAILED IN CNS MINVRT" << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
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
		hp_vbdry(i)->vdirichlet2d();


	return;
}
