#include "bdry_ins.h"
#include <myblas.h>

/*************************************************/
/* SET DIRICHLET BOUNDARY VALUES & FLUXES ********/
/* (THINGS THAT ARE INDEPENDENT OF THE SOLUTION) */
/* BUT NOT INDEPENDENT OF TIME/MESH POSITION     */
/*************************************************/
using namespace bdry_ins;

void generic::output(std::ostream& fout, tet_hp::filetype typ,int tlvl) {
	cout << "warning generic::output in bdry.cpp is being called and is not functioning" << endl;
	return;
}

void neumann::rsdl(int stage) {
	int i,j,k,n,v0,v1,v2,sind,find;
	TinyVector<FLT,3> pt,mvel,nrm;
	Array<FLT,1> u(x.NV),flx(x.NV);
	
	for(j=0;j<base.ntri;++j) {
		find = base.tri(j).gindx;
		v0 = x.tri(find).pnt(0);
		v1 = x.tri(find).pnt(1);
		v2 = x.tri(find).pnt(2);
		
		x.crdtocht2d(find);
		for(n=0;n<tet_mesh::ND;++n)
			basis::tet(x.log2p).proj2d(&x.cht(n)(0),&x.crd2d(n)(0)(0),&x.dcrd2d(n)(0)(0)(0),&x.dcrd2d(n)(1)(0)(0),MXGP);
		
		x.ugtouht2d(find);
		for(n=0;n<x.NV;++n)
			basis::tet(x.log2p).proj2d(&x.uht(n)(0),&x.u2d(n)(0)(0),MXGP);

		for(i=0;i<basis::tet(x.log2p).gpx;++i) {
			for(k=0;k<basis::tet(x.log2p).gpy;++k) {
				/* co-variant vectors: d vec(x)/dxi crossed with d vec(x)/deta */
				nrm(0) = x.dcrd2d(1)(0)(i)(k)*x.dcrd2d(2)(1)(i)(k)-x.dcrd2d(2)(0)(i)(k)*x.dcrd2d(1)(1)(i)(k);
				nrm(1) = -x.dcrd2d(0)(0)(i)(k)*x.dcrd2d(2)(1)(i)(k)+x.dcrd2d(2)(0)(i)(k)*x.dcrd2d(0)(1)(i)(k);
				nrm(2) = x.dcrd2d(0)(0)(i)(k)*x.dcrd2d(1)(1)(i)(k)-x.dcrd2d(1)(0)(i)(k)*x.dcrd2d(0)(1)(i)(k); 
				
				for(n=0;n<tet_mesh::ND;++n) {
					pt(n) = x.crd2d(n)(i)(k);
					mvel(n) = x.gbl->bd(0)*(x.crd2d(n)(i)(k) -dxdt(x.log2p,j)(n)(i)(k));

				}
				
				for(n=0;n<x.NV;++n)
					u(n) = x.u2d(n)(i)(k);
				
				flux(u,pt,mvel,nrm,flx);

				for(n=0;n<x.NV;++n)
					x.res2d(n)(i)(k) = flx(n);

			}
		}

		for(n=0;n<x.NV;++n)
			basis::tet(x.log2p).intgrt2d(&x.lf(n)(0),&x.res2d(n)(0)(0),MXGP);

		for(n=0;n<x.NV;++n)
			x.gbl->res.v(v0,n) += x.lf(n)(0);

		for(n=0;n<x.NV;++n)
			x.gbl->res.v(v1,n) += x.lf(n)(1);
		
		for(n=0;n<x.NV;++n)
			x.gbl->res.v(v2,n) += x.lf(n)(2);
		
		for(i=0;i<3;++i) {
			sind=x.tri(find).seg(i);
			for(k=0;k<basis::tet(x.log2p).em;++k) {
				for(n=0;n<x.NV;++n)
					x.gbl->res.e(sind,k,n) += x.lf(n)(i*basis::tet(x.log2p).em+k+3);
			}
		}
													  
	    for(k=0;k<basis::tet(x.log2p).fm;++k) {
		    for(n=0;n<x.NV;++n)
				x.gbl->res.f(find,k,n) += x.lf(n)(3*basis::tet(x.log2p).em+k+3);
		}
	}
	return;
}

void applied_stress::init(input_map& inmap,void* gbl_in) {
	std::string keyword;
	std::ostringstream nstr;

	neumann::init(inmap,gbl_in);
	
	stress.resize(tet_mesh::ND);
	for(int n=0;n<tet_mesh::ND;++n) {
		nstr.str("");
		nstr << base.idprefix << "_stress" << n << std::flush;
		if (inmap.find(nstr.str()) != inmap.end()) {
			stress(n).init(inmap,nstr.str());
		}
		else {
			*x.gbl->log << "couldn't find stress function " << nstr.str() << '\n';
			exit(1);
		}
	}
	
	return;
}

