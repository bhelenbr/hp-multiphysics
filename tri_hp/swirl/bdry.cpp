#include "bdry_swirl.h"

/*************************************************/
/* SET DIRICHLET BOUNDARY VALUES & FLUXES ********/
/* (THINGS THAT ARE INDEPENDENT OF THE SOLUTION) */
/* BUT NOT INDEPENDENT OF TIME/MESH POSITION     */
/*************************************************/

using namespace bdry_swirl;

void symmetry::tadvance() {
	int j,m,v0,sind,indx;
	TinyVector<FLT,tri_mesh::ND> pt;

	hp_edge_bdry::tadvance();

	/* UPDATE BOUNDARY CONDITION VALUES */
	for(j=0;j<base.nseg;++j) {
		sind = base.seg(j);
		v0 = x.seg(sind).pnt(0);
		x.ug.v(v0,0) = 0.0;
		x.ug.v(v0,2) = 0.0;
	}
	v0 = x.seg(sind).pnt(1);
	x.ug.v(v0,0) = 0.0;
	x.ug.v(v0,2) = 0.0;

	/*******************/    
	/* SET SIDE VALUES */
	/*******************/
	for(j=0;j<base.nseg;++j) {
		sind = base.seg(j);
		indx = sind*x.sm0;
		for(m=0;m<basis::tri(x.log2p)->sm();++m) {
			x.ug.s(sind,m,0) = 0.0;
			x.ug.s(sind,m,2) = 0.0;
		}
	}

	return;
}
