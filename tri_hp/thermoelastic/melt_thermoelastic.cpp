/*
 *  melt_thermoelastic.cpp
 *  tri_hp
 *
 */

#include "melt_thermoelastic.h"

using namespace bdry_thermoelastic;

#ifdef petsc
/* Coupled interface jacobian plus the displacement & continuity flux jacobian */
/* (same as hp_edge_bdry::petsc_jacobian, using the generic element_jacobian) */
void melt::petsc_jacobian() {

	melt_cd::petsc_jacobian();

	const int sm = basis::tri(x.log2p)->sm();
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;

	Array<FLT,2> K(x.NV*(sm+2),vdofs*2 +x.NV*sm);
	Array<int,1> rows(x.NV*(sm+2));
	Array<int,1> cols(vdofs*2 +x.NV*sm);

	elastic_only = true;
	for(int j=0;j<base.nseg;++j) {
		int sind = base.seg(j);

		int rind = 0;
		int cind = 0;
		for(int k=0;k<2;++k) {
			int gindx = vdofs*x.seg(sind).pnt(k);
			for(int n=0;n<x.NV;++n) {
				rows(rind++) = gindx;
				cols(cind++) = gindx++;
			}
			for(int n=x.NV;n<vdofs;++n) {
				cols(cind++) = gindx++;
			}
		}

		/* EDGE MODES */
		if (sm > 0) {
			int gbl_eind = x.npnt*vdofs + sind*sm*x.NV;
			for (int m = 0; m < sm; ++m) {
				for(int n = 0; n < x.NV; ++n) {
					rows(rind++) = gbl_eind;
					cols(cind++) = gbl_eind++;
				}
			}
		}

		hp_edge_bdry::element_jacobian(j,K);
#ifdef MY_SPARSE
		x.J.add_values(x.NV*(sm+2),rows,vdofs*2 +x.NV*sm,cols,K);
#else
		MatSetValuesLocal(x.petsc_J,x.NV*(sm+2),rows.data(),vdofs*2 +x.NV*sm,cols.data(),K.data(),ADD_VALUES);
#endif
	}
	elastic_only = false;
}
#endif
