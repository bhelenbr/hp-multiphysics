#include "bdry_ins.h"
#include <myblas.h>
#include<blitz/tinyvec-et.h>
/*************************************************/
/* SET DIRICHLET BOUNDARY VALUES & FLUXES ********/
/* (THINGS THAT ARE INDEPENDENT OF THE SOLUTION) */
/* BUT NOT INDEPENDENT OF TIME/MESH POSITION     */
/*************************************************/
using namespace bdry_ins;

void generic::output(std::ostream& fout, tet_hp::filetype typ,int tlvl) {
//	cout << "warning generic::output in bdry.cpp is being called and is not functioning" << endl;
	return;
}

void neumann::rsdl(int stage){
	int sind;
	FLT sgn,msgn;
	
	for(int i=0;i<base.ntri;++i){
		
		int find = base.tri(i).gindx;
		x.ugtouht2d(find);
		
		element_rsdl(find,stage);
		
		for(int n=0;n<x.NV;++n)
			x.gbl->res.v(x.tri(find).pnt(0),n) += x.lf(n)(0);
		
		for(int n=0;n<x.NV;++n)
			x.gbl->res.v(x.tri(find).pnt(1),n) += x.lf(n)(1);
		
		for(int n=0;n<x.NV;++n)
			x.gbl->res.v(x.tri(find).pnt(2),n) += x.lf(n)(2);
		
		int indx = 3;
		for(int j=0;j<3;++j) {
			sind=x.tri(find).seg(j);
			sgn = x.tri(find).sgn(j);
			msgn = 1.0;
			for(int k=0;k<basis::tet(x.log2p).em;++k) {
				for(int n=0;n<x.NV;++n)
					x.gbl->res.e(sind,k,n) += msgn*x.lf(n)(indx);
				msgn *= sgn;
				++indx;
			}
		}
		
	    for(int k=0;k<basis::tet(x.log2p).fm;++k) {
		    for(int n=0;n<x.NV;++n)
				x.gbl->res.f(find,k,n) += x.lf(n)(indx);
			++indx;
		}		
	}
	
	return;
}

void neumann::element_rsdl(int find,int stage) {
	int j,k,n;
	TinyVector<FLT,3> pt,mvel,nrm,vec1,vec2;
	Array<FLT,1> u(x.NV),flx(x.NV);
	
	x.lf = 0.0;
	
	x.crdtocht2d(find);
	for(n=0;n<tet_mesh::ND;++n)
		basis::tet(x.log2p).proj2d(&x.cht(n)(0),&x.crd2d(n)(0)(0),&x.dcrd2d(n)(0)(0)(0),&x.dcrd2d(n)(1)(0)(0),MXGP);
	
	//x.ugtouht2d(find);
	for(n=0;n<x.NV;++n)
		basis::tet(x.log2p).proj2d(&x.uht(n)(0),&x.u2d(n)(0)(0),MXGP);
	
	for(j=0;j<basis::tet(x.log2p).gpx;++j) {
		for(k=0;k<basis::tet(x.log2p).gpy;++k) {
			/* co-variant vectors: d vec(x)/dxi crossed with d vec(x)/deta */
			for(n=0;n<3;++n){
				vec1(n)=x.dcrd2d(n)(0)(j)(k);
				vec2(n)=x.dcrd2d(n)(1)(j)(k);
			}
			nrm=cross(vec1,vec2);
			
			for(n=0;n<tet_mesh::ND;++n) {
				pt(n) = x.crd2d(n)(j)(k);
				mvel(n) = 0.0;//x.gbl->bd(0)*(x.crd2d(n)(j)(k) -dxdt(x.log2p,j)(n)(j)(k));
			}
			
			for(n=0;n<x.NV;++n)
				u(n) = x.u2d(n)(j)(k);
			
			flux(u,pt,mvel,nrm,flx);
			
			for(n=0;n<x.NV;++n)
				x.res2d(n)(j)(k) = flx(n);			
		}
	}
	
	for(n=0;n<x.NV;++n)
		basis::tet(x.log2p).intgrt2d(&x.lf(n)(0),&x.res2d(n)(0)(0),MXGP);
	
	return;
}


//void neumann::rsdl(int stage) {
//	int i,j,k,n,v0,v1,v2,sind,find,indx;
//	FLT sgn,msgn;
//	TinyVector<FLT,3> pt,mvel,nrm,vec1,vec2;
//	Array<FLT,1> u(x.NV),flx(x.NV);
//	
//	for(i=0;i<base.ntri;++i) {
//		find = base.tri(i).gindx;
//		v0 = x.tri(find).pnt(0);
//		v1 = x.tri(find).pnt(1);
//		v2 = x.tri(find).pnt(2);
//		
//		x.crdtocht2d(find);
//		for(n=0;n<tet_mesh::ND;++n)
//			basis::tet(x.log2p).proj2d(&x.cht(n)(0),&x.crd2d(n)(0)(0),&x.dcrd2d(n)(0)(0)(0),&x.dcrd2d(n)(1)(0)(0),MXGP);
//		
//		x.ugtouht2d(find);
//		for(n=0;n<x.NV;++n)
//			basis::tet(x.log2p).proj2d(&x.uht(n)(0),&x.u2d(n)(0)(0),MXGP);
//
//		for(j=0;j<basis::tet(x.log2p).gpx;++j) {
//			for(k=0;k<basis::tet(x.log2p).gpy;++k) {
//				/* co-variant vectors: d vec(x)/dxi crossed with d vec(x)/deta */
//				for(n=0;n<3;++n){
//					vec1(n)=x.dcrd2d(n)(0)(j)(k);
//					vec2(n)=x.dcrd2d(n)(1)(j)(k);
//				}
//				nrm=cross(vec1,vec2);
//
//				for(n=0;n<tet_mesh::ND;++n) {
//					pt(n) = x.crd2d(n)(j)(k);
//					mvel(n) = 0.0;//x.gbl->bd(0)*(x.crd2d(n)(j)(k) -dxdt(x.log2p,j)(n)(j)(k));
//				}
//				
//				for(n=0;n<x.NV;++n)
//					u(n) = x.u2d(n)(j)(k);
//				
//				flux(u,pt,mvel,nrm,flx);
//
//				for(n=0;n<x.NV;++n)
//					x.res2d(n)(j)(k) = flx(n);
//
//			}
//		}
//
//		for(n=0;n<x.NV;++n)
//			basis::tet(x.log2p).intgrt2d(&x.lf(n)(0),&x.res2d(n)(0)(0),MXGP);
//
//		for(n=0;n<x.NV;++n)
//			x.gbl->res.v(v0,n) += x.lf(n)(0);
//
//		for(n=0;n<x.NV;++n)
//			x.gbl->res.v(v1,n) += x.lf(n)(1);
//		
//		for(n=0;n<x.NV;++n)
//			x.gbl->res.v(v2,n) += x.lf(n)(2);
//		
//		indx = 3;
//		for(j=0;j<3;++j) {
//			sind=x.tri(find).seg(j);
//			sgn = x.tri(find).sgn(j);
//			msgn = 1.0;
//			for(k=0;k<basis::tet(x.log2p).em;++k) {
//				for(n=0;n<x.NV;++n)
//					x.gbl->res.e(sind,k,n) += msgn*x.lf(n)(indx);
//				msgn *= sgn;
//				++indx;
//			}
//		}
//													  
//	    for(k=0;k<basis::tet(x.log2p).fm;++k) {
//		    for(n=0;n<x.NV;++n)
//				x.gbl->res.f(find,k,n) += x.lf(n)(indx);
//			++indx;
//		}		
//	}
//	return;
//}

void symmetry::tadvance() {
	int v0,sind,find;
	TinyVector<FLT,tet_mesh::ND> pt;
	
	hp_face_bdry::tadvance();
	
	/* UPDATE BOUNDARY CONDITION VALUES */
	for(int j=0;j<base.npnt;++j) {
		v0 = base.pnt(j).gindx;
		x.ug.v(v0,dir) = 0.0;
	}

	if (basis::tet(x.log2p).em > 0) {
		for(int j=0;j<base.nseg;++j) {
			sind = base.seg(j).gindx;
			x.ug.e(sind,Range(0,basis::tet(x.log2p).em-1),dir) = 0.0;
		}
	}

	if (basis::tet(x.log2p).fm > 0) {
		for(int j=0;j<base.ntri;++j) {
			find = base.tri(j).gindx;
			x.ug.f(find,Range(0,basis::tet(x.log2p).fm-1),dir) = 0.0;
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

