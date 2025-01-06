#include "bdry_cns.h"
//#include <myblas.h>
//#include<blitz/tinyvec-et.h>
/*************************************************/
/* SET DIRICHLET BOUNDARY VALUES & FLUXES ********/
/* (THINGS THAT ARE INDEPENDENT OF THE SOLUTION) */
/* BUT NOT INDEPENDENT OF TIME/MESH POSITION     */
/*************************************************/
using namespace bdry_cns;

/* ---------------------- */
/* face boundary routines */
/* ---------------------- */

void neumann::rsdl(int stage){
	int sind;
	FLT sgn,msgn;
	
	for(int i=0;i<base.ntri;++i){
		
		int find = base.tri(i).gindx;
		x.ugtouht2d(find);
		
		element_rsdl(find,stage);
		
		for(int n=0;n<x.NV;++n)
			x.hp_gbl->res.v(x.tri(find).pnt(0),n) += x.lf(n)(0);
		
		for(int n=0;n<x.NV;++n)
			x.hp_gbl->res.v(x.tri(find).pnt(1),n) += x.lf(n)(1);
		
		for(int n=0;n<x.NV;++n)
			x.hp_gbl->res.v(x.tri(find).pnt(2),n) += x.lf(n)(2);
		
		int indx = 3;
		for(int j=0;j<3;++j) {
			sind=x.tri(find).seg(j);
			sgn = x.tri(find).sgn(j);
			msgn = 1.0;
			for(int k=0;k<basis::tet(x.log2p).em;++k) {
				for(int n=0;n<x.NV;++n)
					x.hp_gbl->res.e(sind,k,n) += msgn*x.lf(n)(indx);
				msgn *= sgn;
				++indx;
			}
		}
		
	    for(int k=0;k<basis::tet(x.log2p).fm;++k) {
		    for(int n=0;n<x.NV;++n)
				x.hp_gbl->res.f(find,k,n) += x.lf(n)(indx);
			++indx;
		}		
	}
	
	return;
}

void neumann::element_rsdl(int find,int stage) {
	int j,k,n;
	TinyVector<FLT,3> pt,mvel,nrm,vec1,vec2;
	Array<FLT,1> u(x.NV),flx(x.NV);
	
	beta_squared = x.hp_cns_gbl->betasquared(x.tri(find).tet(0));
	
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


void neumann::flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm, Array<FLT,1>& flx) {
	Array<FLT,1> v(3);
	FLT pr = ibc->f(0,xpt,x.gbl->time);
	v(0) = u(1)-mv(0);
	v(1) = u(2)-mv(1);
	v(2) = u(3)-mv(2);
	FLT RT = u(4);
	FLT rho = (pr+x.hp_cns_gbl->atm_pressure)/RT;
	
	/* CONTINUITY */
	flx(0) = rho*(v(0)*norm(0)+v(1)*norm(1)+v(2)*norm(2));
	
	/* XYZ MOMENTUM */
	for (int n=1;n<tet_mesh::ND+1;++n)
		flx(n) = flx(0)*v(n-1) + pr*norm(n-1);
	
	/* ENERGY EQUATION */
	FLT h = x.hp_cns_gbl->gamma/(x.hp_cns_gbl->gamma-1.0)*RT +0.5*(v(0)*v(0)+v(1)*v(1)+v(2)*v(2));				
	flx(x.NV-1) = h*flx(0);
	
	return;
}

void inflow::flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm,  Array<FLT,1>& flx) {
	Array<FLT,1> v(3);
	FLT pr = u(0);
	v(0) = ibc->f(1,xpt,x.gbl->time)-mv(0);
	v(1) = ibc->f(2,xpt,x.gbl->time)-mv(1);
	v(2) = ibc->f(3,xpt,x.gbl->time)-mv(2);
	FLT RT = ibc->f(4,xpt,x.gbl->time);
	FLT rho = (pr+x.hp_cns_gbl->atm_pressure)/RT;

	/* CONTINUITY */
	flx(0) = rho*(v(0)*norm(0)+v(1)*norm(1)+v(2)*norm(2));
	
	/* XYZ MOMENTUM */
	for (int n=1;n<tet_mesh::ND+1;++n)
		flx(n) = flx(0)*v(n-1) + pr*norm(n-1);
	
	/* ENERGY EQUATION */
	FLT h = x.hp_cns_gbl->gamma/(x.hp_cns_gbl->gamma-1.0)*RT +0.5*(v(0)*v(0)+v(1)*v(1)+v(2)*v(2));				
	flx(x.NV-1) = h*flx(0);

	return;
}

//void inflow::flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm, Array<FLT,1>& flx) {	
//	
//	
//	TinyVector<FLT,5> fluxtemp,ub;
//	TinyVector<FLT,3> vec1,vec2,vec3,vecr;
//	FLT gogm1 = x.hp_cns_gbl->gamma/(x.hp_cns_gbl->gamma-1.0);
//	
//	FLT mag = sqrt(norm(0)*norm(0)+norm(1)*norm(1)+norm(2)*norm(2));
//	norm /= mag;
//	
//	/* align u to the normal */
//	vec1 = norm;
//	/* find random vector different than norm */
//	vecr(0) = vec1(1),	vecr(1) = vec1(2), vecr(2) = vec1(0);
//	vec2 = cross(norm, vecr);
//	vec3 = cross(norm, vec2);
//	
//	for(int n=0;n<x.NV;++n)
//		ub(n) = ibc->f(n,xpt,x.gbl->time);
//	
//	/* Rotate Coordinate System so that U aligns with normal vector */
//	FLT pr = u(0);
//	FLT ur = ub(1)*vec1(0)+ub(2)*vec1(1)+ub(3)*vec1(2);	
//	FLT vr = ub(1)*vec2(0)+ub(2)*vec2(1)+ub(3)*vec2(2);
//	FLT wr = ub(1)*vec3(0)+ub(2)*vec3(1)+ub(3)*vec3(2);
//	FLT RT = ub(4);
//	
//	FLT rho = (pr+x.hp_cns_gbl->atm_pressure)/RT;
//	
//	fluxtemp(0) = rho*ur;
//	fluxtemp(1) = rho*ur*ur+pr;
//	fluxtemp(2) = rho*ur*vr;
//	fluxtemp(3) = rho*ur*wr;
//	fluxtemp(4) = rho*ur*(gogm1*RT+0.5*(ur*ur+vr*vr+wr*wr));	
//	
//	/* CHANGE BACK TO X,Y,Z COORDINATES */
//	FLT temp2 = 1.0/(1.0-norm(2)*norm(1)-norm(1)*norm(0)-norm(0)*norm(2))/(1.0+norm(2)*norm(1)+norm(1)*norm(0)+norm(0)*norm(2));
//	vec1(0) = norm(0), vec1(1) = temp2*(norm(1)*norm(0)-norm(2)*norm(2)), vec1(2) = temp2*(norm(0)*norm(1)*norm(2)-norm(1)*norm(1)*norm(1)-norm(2)*norm(2)*norm(1)+norm(0)*norm(0)*norm(2));
//	vec2(0) = norm(1), vec2(1) = temp2*(norm(2)*norm(1)-norm(0)*norm(0)), vec2(2) = temp2*(norm(0)*norm(1)*norm(2)-norm(2)*norm(2)*norm(2)-norm(0)*norm(0)*norm(2)+norm(1)*norm(1)*norm(0));
//	vec3(0) = norm(2), vec3(1) = temp2*(norm(2)*norm(0)-norm(1)*norm(1)), vec3(2) = temp2*(norm(0)*norm(1)*norm(2)-norm(0)*norm(0)*norm(0)-norm(1)*norm(1)*norm(0)+norm(2)*norm(2)*norm(1));
//	
//	flx(0) = fluxtemp(0);
//	flx(1) = fluxtemp(1)*vec1(0) + fluxtemp(2)*vec1(1) + fluxtemp(3)*vec1(2);
//	flx(2) = fluxtemp(1)*vec2(0) + fluxtemp(2)*vec2(1) + fluxtemp(3)*vec2(2);
//	flx(3) = fluxtemp(1)*vec3(0) + fluxtemp(2)*vec3(1) + fluxtemp(3)*vec3(2);
//	flx(4) = fluxtemp(4);
//	
//	flx *= mag;
//	
//	return;
//}


void inflow::vdirichlet() {	
	for(int j=0;j<base.npnt;++j) {
		int v0 = base.pnt(j).gindx;
		for (int n=0; n<ndirichlets; ++n) 
			x.hp_gbl->res.v(v0,dirichlets(n)) = 0.0;
		
	}
}

void inflow::edirichlet() {
	//if (basis::tet(x.log2p).em > 0) {
		for(int j=0;j<base.nseg;++j) {
			int sind = base.seg(j).gindx;
			for (int n=0; n<ndirichlets; ++n) 
				x.hp_gbl->res.e(sind,Range::all(),dirichlets(n)) = 0.0;
		}
	//}
}		

void inflow::fdirichlet() {
	if (basis::tet(x.log2p).fm > 0) {
		for(int j=0;j<base.ntri;++j) {
			int find = base.tri(j).gindx;
			for (int n=0; n<ndirichlets; ++n) 
				x.hp_gbl->res.f(find,Range::all(),dirichlets(n)) = 0.0;
		}
	}
}	

//void inflow::apply_sparse_dirichlet(bool compressed_column) {
//	int gind;
//	int em=basis::tet(x.log2p).em;
//	int fm=basis::tet(x.log2p).fm;
//	
//	for(int i=0;i<base.npnt;++i){
//		gind = base.pnt(i).gindx*x.NV;
//		for(int n=0;n<ndirichlets;++n)
//			x.sparse_dirichlet(gind+dirichlets(n),compressed_column);
//	}
//	
//	for(int i=0;i<base.nseg;++i){
//		gind = x.npnt*x.NV+base.seg(i).gindx*em*x.NV;
//		for(int m=0; m<em; ++m)
//			for(int n=0;n<ndirichlets;++n)
//				x.sparse_dirichlet(gind+m*x.NV+dirichlets(n),compressed_column);
//	}
//	
//	for(int i=0;i<base.ntri;++i){
//		gind = x.npnt*x.NV+x.nseg*em*x.NV+base.tri(i).gindx*fm*x.NV;
//		for(int m=0; m<fm; ++m)
//			for(int n=0;n<ndirichlets;++n)
//				x.sparse_dirichlet(gind+m*x.NV+dirichlets(n),compressed_column);
//	}			
//}

//void inflow::modify_boundary_residual() {
//	int j,k,m,n,v0,v1,sind;
//	TinyVector<FLT,tet_mesh::ND> pt;
//	TinyVector<double,MXGP> res1d;
//	TinyVector<double,MXTM> rescoef;
//	
//	FLT ogm1 = 1.0/(x.hp_cns_gbl->gamma-1.0);
//	
//	for(int j=0;j<base.npnt;++j) {
//		int v0 = base.pnt(j).gindx;
//		FLT KE = 0.5*(ibc->f(1, x.pnts(v0), x.gbl->time)*ibc->f(1, x.pnts(v0), x.gbl->time)+ibc->f(2, x.pnts(v0), x.gbl->time)*ibc->f(2, x.pnts(v0), x.gbl->time)+ibc->f(3, x.pnts(v0), x.gbl->time)*ibc->f(3, x.pnts(v0), x.gbl->time));
//		x.hp_gbl->res.v(v0,1) = x.hp_gbl->res.v(v0,0)*ibc->f(1, x.pnts(v0), x.gbl->time);
//		x.hp_gbl->res.v(v0,2) = x.hp_gbl->res.v(v0,0)*ibc->f(2, x.pnts(v0), x.gbl->time);
//		x.hp_gbl->res.v(v0,3) = x.hp_gbl->res.v(v0,0)*ibc->f(3, x.pnts(v0), x.gbl->time);
//		x.hp_gbl->res.v(v0,4) = x.hp_gbl->res.v(v0,0)*(ibc->f(4, x.pnts(v0), x.gbl->time)*ogm1+KE);
//	}
//	
//	if(basis::tet(x.log2p).em) {
//		for(j=0;j<base.nseg;++j) {
//			sind = base.seg(j).gindx;
//			v0 = x.seg(sind).pnt(0);
//			v1 = x.seg(sind).pnt(1);
//			
//			if (is_curved()) {
//				x.crdtocht1d(sind);
//				for(n=0;n<tet_mesh::ND;++n)
//					basis::tet(x.log2p).proj1d(&x.cht(n)(0),&x.crd1d(n)(0));
//			}
//			else {
//				for(n=0;n<tet_mesh::ND;++n) {
//					basis::tet(x.log2p).proj1d(x.pnts(v0)(n),x.pnts(v1)(n),&x.crd1d(n)(0));
//					
//					for(k=0;k<basis::tet(x.log2p).gpx;++k)
//						x.dcrd1d(n)(k) = 0.5*(x.pnts(v1)(n)-x.pnts(v0)(n));
//				}
//			}
//			
//			/* take global coefficients and put into local vector */
//			/*  only need res_rho */
//			for (m=0; m<2; ++m) 
//				rescoef(m) = x.hp_gbl->res.v(x.seg(sind).pnt(m),0);					
//			
//			for (m=0;m<basis::tet(x.log2p).em;++m) 
//				rescoef(m+2) = x.hp_gbl->res.e(sind,m,0);					
//			
//			basis::tet(x.log2p).proj1d(&rescoef(0),&res1d(0));
//			
//			for(n=1;n<x.NV;++n)
//				basis::tet(x.log2p).proj1d(x.hp_gbl->res.v(v0,n),x.hp_gbl->res.v(v1,n),&x.res1d(n)(0));
//			
//			for(k=0;k<basis::tet(x.log2p).gpx; ++k) {
//				pt(0) = x.crd1d(0)(k);
//				pt(1) = x.crd1d(1)(k);
//				
//				FLT KE = 0.5*(ibc->f(1, pt, x.gbl->time)*ibc->f(1, pt, x.gbl->time)+ibc->f(2, pt, x.gbl->time)*ibc->f(2, pt, x.gbl->time)+ibc->f(3, pt, x.gbl->time)*ibc->f(3, pt, x.gbl->time));
//				x.res1d(1)(k) -= res1d(k)*ibc->f(1, pt, x.gbl->time);
//				x.res1d(2)(k) -= res1d(k)*ibc->f(2, pt, x.gbl->time);
//				x.res1d(3)(k) -= res1d(k)*ibc->f(3, pt, x.gbl->time);
//				x.res1d(4)(k) -= res1d(k)*(ibc->f(4, pt, x.gbl->time)*ogm1+KE);
//				
//			}
//			
//			for(n=1;n<x.NV;++n){
//				basis::tet(x.log2p).intgrt1d(&x.lf(n)(0),&x.res1d(n)(0));
//				
//				for(m=0;m<basis::tet(x.log2p).em;++m) 
//					x.hp_gbl->res.e(sind,m,n) = -x.lf(n)(2+m)*basis::tet(x.log2p).diag1d(m);					
//				
//			}								
//		}
//	}
//	
//	if (basis::tet(x.log2p).fm) {
//		*x.gbl->log << "warning face modes not implemented in cns inflow modify boundary residual" << endl;
//	}
//	return;
//}

void adiabatic::flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm,  Array<FLT,1>& flx) {
	Array<FLT,1> v(3);
	FLT pr = u(0);
	v(0) = ibc->f(1,xpt,x.gbl->time);
	v(1) = ibc->f(2,xpt,x.gbl->time);
	v(2) = ibc->f(3,xpt,x.gbl->time);
	FLT RT =u(4);
	FLT rho = (pr+x.hp_cns_gbl->atm_pressure)/RT;
	
	/* CONTINUITY */
	flx(0) = rho*(v(0)*norm(0)+v(1)*norm(1)+v(2)*norm(2));
	
	/* XYZ MOMENTUM */
	for (int n=1;n<tet_mesh::ND+1;++n)
		flx(n) = flx(0)*v(n-1) + pr*norm(n-1);
	
	/* ENERGY EQUATION */
	FLT h = x.hp_cns_gbl->gamma/(x.hp_cns_gbl->gamma-1.0)*RT +0.5*(v(0)*v(0)+v(1)*v(1)+v(2)*v(2));				
	flx(x.NV-1) = h*flx(0);
	
	return;
}

void adiabatic::vdirichlet() {			
	for(int j=0;j<base.npnt;++j) {
		int v0 = base.pnt(j).gindx;		
		for (int n=0; n<ndirichlets; ++n) 
			x.hp_gbl->res.v(v0,dirichlets(n)) = 0.0;		
	}
}

void adiabatic::edirichlet() {
	if (basis::tet(x.log2p).em > 0) {
		for(int j=0;j<base.nseg;++j) {
			int sind = base.seg(j).gindx;
			for (int n=0; n<ndirichlets; ++n) 
				x.hp_gbl->res.e(sind,Range::all(),dirichlets(n)) = 0.0;
		}
	}
}		

void adiabatic::fdirichlet() {
	if (basis::tet(x.log2p).fm > 0) {
		for(int j=0;j<base.ntri;++j) {
			int find = base.tri(j).gindx;
			for (int n=0; n<ndirichlets; ++n) 
				x.hp_gbl->res.f(find,Range::all(),dirichlets(n)) = 0.0;
		}
	}
}		

//void adiabatic::apply_sparse_dirichlet(bool compressed_column) {
//	int gind;
//	int em=basis::tet(x.log2p).em;
//	int fm=basis::tet(x.log2p).fm;
//	
//	for(int i=0;i<base.npnt;++i){
//		gind = base.pnt(i).gindx*x.NV;
//		for(int n=0;n<ndirichlets;++n)
//			x.sparse_dirichlet(gind+dirichlets(n),compressed_column);
//	}
//	
//	for(int i=0;i<base.nseg;++i){
//		gind = x.npnt*x.NV+base.seg(i).gindx*em*x.NV;
//		for(int m=0; m<em; ++m)
//			for(int n=0;n<ndirichlets;++n)
//				x.sparse_dirichlet(gind+m*x.NV+dirichlets(n),compressed_column);
//	}
//	
//	for(int i=0;i<base.ntri;++i){
//		gind = x.npnt*x.NV+x.nseg*em*x.NV+base.tri(i).gindx*fm*x.NV;
//		for(int m=0; m<fm; ++m)
//			for(int n=0;n<ndirichlets;++n)
//				x.sparse_dirichlet(gind+m*x.NV+dirichlets(n),compressed_column);
//	}			
//}

//void adiabatic::modify_boundary_residual() {
//	int j,k,m,n,v0,v1,sind;
//	TinyVector<FLT,tet_mesh::ND> pt;
//	TinyVector<double,MXGP> res1d;
//	TinyVector<double,MXTM> rescoef;
//	
//	for(int j=0;j<base.npnt;++j) {
//		int v0 = base.pnt(j).gindx;
//		x.hp_gbl->res.v(v0,1) = x.hp_gbl->res.v(v0,0)*ibc->f(1, x.pnts(v0), x.gbl->time);
//		x.hp_gbl->res.v(v0,2) = x.hp_gbl->res.v(v0,0)*ibc->f(2, x.pnts(v0), x.gbl->time);
//		x.hp_gbl->res.v(v0,3) = x.hp_gbl->res.v(v0,0)*ibc->f(3, x.pnts(v0), x.gbl->time);
//	}
//	
//	if(basis::tet(x.log2p).em) {
//		for(j=0;j<base.nseg;++j) {
//			sind = base.seg(j).gindx;
//			v0 = x.seg(sind).pnt(0);
//			v1 = x.seg(sind).pnt(1);
//			
//			if (is_curved()) {
//				x.crdtocht1d(sind);
//				for(n=0;n<tet_mesh::ND;++n)
//					basis::tet(x.log2p).proj1d(&x.cht(n)(0),&x.crd1d(n)(0));
//			}
//			else {
//				for(n=0;n<tet_mesh::ND;++n) {
//					basis::tet(x.log2p).proj1d(x.pnts(v0)(n),x.pnts(v1)(n),&x.crd1d(n)(0));
//					
//					for(k=0;k<basis::tet(x.log2p).gpx;++k)
//						x.dcrd1d(n)(k) = 0.5*(x.pnts(v1)(n)-x.pnts(v0)(n));
//				}
//			}
//			
//			/* take global coefficients and put into local vector */
//			/*  only need res_rho */
//			for (m=0; m<2; ++m) 
//				rescoef(m) = x.hp_gbl->res.v(x.seg(sind).pnt(m),0);					
//			
//			for (m=0;m<basis::tet(x.log2p).em;++m) 
//				rescoef(m+2) = x.hp_gbl->res.e(sind,m,0);					
//			
//			basis::tet(x.log2p).proj1d(&rescoef(0),&res1d(0));
//			
//			for(n=1;n<x.NV;++n)
//				basis::tet(x.log2p).proj1d(x.hp_gbl->res.v(v0,n),x.hp_gbl->res.v(v1,n),&x.res1d(n)(0));
//			
//			for(k=0;k<basis::tet(x.log2p).gpx; ++k) {
//				pt(0) = x.crd1d(0)(k);
//				pt(1) = x.crd1d(1)(k);
//				
//				x.res1d(1)(k) -= res1d(k)*ibc->f(1, pt, x.gbl->time);
//				x.res1d(2)(k) -= res1d(k)*ibc->f(2, pt, x.gbl->time);
//				x.res1d(3)(k) -= res1d(k)*ibc->f(3, pt, x.gbl->time);
//				
//			}
//			
//			for(n=1;n<x.NV;++n){
//				basis::tet(x.log2p).intgrt1d(&x.lf(n)(0),&x.res1d(n)(0));
//				
//				for(m=0;m<basis::tet(x.log2p).em;++m) 
//					x.hp_gbl->res.e(sind,m,n) = -x.lf(n)(2+m)*basis::tet(x.log2p).diag1d(m);					
//				
//			}								
//		}
//	}
//	if (basis::tet(x.log2p).fm) {
//		*x.gbl->log << "warning face modes not implemented in cns adiabatic modify boundary residual" << endl;
//	}
//	
//	return;
//}


void characteristic::flux(Array<FLT,1>& pvu, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm, Array<FLT,1>& flx) {	
	
	TinyVector<FLT,5> lambda,Rl,Rr,ul,ur,ub,Roe,fluxtemp,fluxleft, fluxright;
	TinyVector<FLT,3> vec1,vec2,vec3,vecr;
	Array<FLT,2> A(x.NV,x.NV),V(x.NV,x.NV),VINV(x.NV,x.NV),temp(x.NV,x.NV),P(x.NV,x.NV),Pinv(x.NV,x.NV),dpdc(x.NV,x.NV), dcdp(x.NV,x.NV);
	Array<FLT,1> Aeigs(x.NV);
	FLT gam = x.hp_cns_gbl->gamma;
	FLT gm1 = gam-1.0;
	FLT gogm1 = gam/gm1;
	
	FLT mag = sqrt(norm(0)*norm(0)+norm(1)*norm(1)+norm(2)*norm(2));
	norm /= mag;

	/* align u to the normal */
	vec1 = norm;
	/* find random vector different than norm */
    vecr(0) = vec1(1);	vecr(1) = vec1(2); vecr(2) = vec1(0);
	vec2 = cross(norm, vecr);
	vec3 = cross(norm, vec2);
	
	/* Left */
	/* Rotate Coordinate System so that U aligns with normal vector */
	ul(0) = pvu(0);
	ul(1) = pvu(1)*vec1(0)+pvu(2)*vec1(1)+pvu(3)*vec1(2);	
	ul(2) = pvu(1)*vec2(0)+pvu(2)*vec2(1)+pvu(3)*vec2(2);
	ul(3) = pvu(1)*vec3(0)+pvu(2)*vec3(1)+pvu(3)*vec3(2);
	ul(4) = pvu(4);

	/* Roe Variables */
	Rl(0) = sqrt((ul(0)+x.hp_cns_gbl->atm_pressure)/ul(x.NV-1)); // sqrt(rho)
	Rl(1) = Rl(0)*ul(1); // sqrt(rho)*u
	Rl(2) = Rl(0)*ul(2); // sqrt(rho)*v
	Rl(3) = Rl(0)*ul(3); // sqrt(rho)*w
	Rl(4) = Rl(0)*(ul(x.NV-1)/gm1+0.5*(ul(1)*ul(1)+ul(2)*ul(2)+ul(3)*ul(3))); // sqrt(rho)*E
	
	/* Right */
	for(int n=0;n<x.NV;++n)
		ub(n) = ibc->f(n,xpt,x.gbl->time);
	
	/* Rotate Coordinate System */
	ur(0) = ub(0);
	ur(1) = ub(1)*vec1(0)+ub(2)*vec1(1)+ub(3)*vec1(2);	
	ur(2) = ub(1)*vec2(0)+ub(2)*vec2(1)+ub(3)*vec2(2);
	ur(3) = ub(1)*vec3(0)+ub(2)*vec3(1)+ub(3)*vec3(2);
	ur(4) = ub(4);
	
	/* Roe Variables */
	Rr(0) = sqrt((ur(0)+x.hp_cns_gbl->atm_pressure)/ur(x.NV-1));
	Rr(1) = Rr(0)*ur(1);
	Rr(2) = Rr(0)*ur(2);
	Rr(3) = Rr(0)*ur(3);
	Rr(4) = Rr(0)*(ur(x.NV-1)/gm1+0.5*(ur(1)*ur(1)+ur(2)*ur(2)+ur(3)*ur(3)));	
	
	/* Average left and right flux */
	fluxleft(0) = (ul(0)+x.hp_cns_gbl->atm_pressure)*ul(1)/ul(x.NV-1);
	fluxleft(1) = fluxleft(0)*ul(1)+ul(0);
	fluxleft(2) = fluxleft(0)*ul(2);
	fluxleft(3) = fluxleft(0)*ul(3);
	fluxleft(4) = fluxleft(0)*(gogm1*ul(x.NV-1)+0.5*(ul(1)*ul(1)+ul(2)*ul(2)+ul(3)*ul(3)));
	
	fluxright(0) = (ur(0)+x.hp_cns_gbl->atm_pressure)*ur(1)/ur(x.NV-1);
	fluxright(1) = fluxright(0)*ur(1)+ur(0);
	fluxright(2) = fluxright(0)*ur(2);
	fluxright(3) = fluxright(0)*ur(3);
	fluxright(4) = fluxright(0)*(gogm1*ur(x.NV-1)+0.5*(ur(1)*ur(1)+ur(2)*ur(2)+ur(3)*ur(3)));
	
	/* Average Flux */
	fluxtemp = 0.5*(fluxleft+fluxright);
	
	/* Average Roe Variables */
	Roe = 0.5*(Rl+Rr);
	
	/* Calculate u,v,c Variables */
	FLT rho = Roe(0)*Roe(0);
	FLT u = Roe(1)/Roe(0);
	FLT v = Roe(2)/Roe(0);
	FLT w = Roe(3)/Roe(0);
	FLT ke = 0.5*(u*u+v*v+w*w);
	FLT E = Roe(4)/Roe(0);
	FLT rt = gm1*(E-ke);
	//FLT pr = rho*rt;
	FLT c2 = gam*rt;
	FLT c = sqrt(c2);
	
	/* stuff for preconditioner */
	FLT nu = x.hp_cns_gbl->mu/rho;
	FLT cp = gogm1*x.hp_cns_gbl->R;
	FLT alpha = x.hp_cns_gbl->kcond/(rho*cp);
	
	/* make sure this stuff matches what is in tstep */
	FLT h =  mag*2.0/(0.25*(basis::tet(x.log2p).p+1.0)*(basis::tet(x.log2p).p+1.0));
	FLT hdt = 0.5*h*x.gbl->bd(0)/c;
	FLT umag = sqrt(u*u+v*v+w*w);
	FLT M = MAX(umag/c,1.0e-5);
	FLT nuh = 2.0*MAX(4.0*nu/(3.0*h*c),alpha/(h*c));
	
	FLT b2;
	if(M < 0.8 && x.hp_cns_gbl->preconditioner){
	    b2 = MIN(M*M/(1.0-M*M) + hdt*hdt + nuh*nuh, 1.0);
	}
	else {
		b2 = 1.0; // turn off preconditioner
	}

	/* b2 computed on cell in tstep */
	b2 = beta_squared;
	
	//cout << "characteristic BC" << endl;
	//cout << sqrt(b2) << ' ' << sqrt(M*M/(1.0-M*M)) << ' ' << hdt << ' ' << nuh << ' ' << 4.0*nu/(3.0*h*c) << ' ' << alpha/(h*c) << endl;

	/* Inverse of Preconditioner */
	Pinv = 1.0/b2,					 0.0, 0.0, 0.0, 0.0,
		   0.0,						 1.0, 0.0, 0.0, 0.0,
		   0.0,						 0.0, 1.0, 0.0, 0.0,
		   0.0,						 0.0, 0.0, 1.0, 0.0,
		   -(b2-1.0)/(gogm1*rho*b2), 0.0, 0.0, 0.0, 1.0;
	
	/* jacobian of conservative wrt primitive */
	dcdp = 1.0/rt,               0.0,   0.0,   0.0,   -rho/rt,
		   u/rt,                 rho,   0.0,   0.0,   -rho*u/rt,
		   v/rt,                 0.0,   rho,   0.0,   -rho*v/rt,
		   w/rt,                 0.0,   0.0,   rho,   -rho*w/rt,
		   (rt+gm1*ke)/(gm1*rt), rho*u, rho*v, rho*w, -rho*ke/rt;		
	
	temp = 0.0;
	for(int i=0; i<x.NV; ++i)
		for(int j=0; j<x.NV; ++j)
			for(int k=0; k<x.NV; ++k)
				temp(i,j)+=dcdp(i,k)*Pinv(k,j);
	
	Pinv = temp;
	
	// need to fill in Tinv not sure if it is needed
	//	if(hp_cns_gbl->preconditioner == 2) {
	//		temp = 0.0;
	//		for(int i=0; i<NV; ++i)
	//			for(int j=0; j<NV; ++j)
	//				for(int k=0; k<NV; ++k)
	//					temp(i,j)+=Pinv(i,k)*Tinv(k,j);
	//		
	//		Pinv = temp/spectral_radius(Tinv);
	//	}
	
	
	FLT temp1 = sqrt(u*u*(1.0-2.0*b2+b2*b2)+4.0*b2*c2);
	
	V = 0.5*(u*(b2-1.0)+temp1)*rho,			0.5*(u*(b2-1.0)-temp1)*rho,			0.0, 0.0, 0.0,
		1.0,								1.0,								0.0, 0.0, 0.0,
		0.0,								0.0,								1.0, 0.0, 0.0,
		0.0,								0.0,								0.0, 1.0, 0.0,
		0.5*(u*(b2-1.0)*gm1+gm1*temp1)/gam, 0.5*(u*(b2-1.0)*gm1-gm1*temp1)/gam, 0.0, 0.0, 1.0;
	
	Aeigs = 0.5*(u+u*b2+temp1), 0.5*(u+u*b2-temp1),u,u,u;

	for(int i=0; i<x.NV; ++i)
		Aeigs(i) = abs(Aeigs(i));
	
	VINV = 1.0/(rho*temp1),	 0.5*(temp1-u*(b2-1.0))/temp1, 0.0, 0.0, 0.0,
		   -1.0/(rho*temp1), 0.5*(u*(b2-1.0)+temp1)/temp1, 0.0, 0.0, 0.0,
		   0.0,				 0.0,						   1.0, 0.0, 0.0,
		   0.0,				 0.0,						   0.0, 1.0, 0.0,
		   -gm1/(gam*rho),	 0.0,						   0.0, 0.0, 1.0;					
	
	for(int i=0; i < x.NV; ++i)
		for(int j=0; j < x.NV; ++j)
			VINV(i,j) = Aeigs(i)*VINV(i,j);
	
	A = 0.0;
	for(int i=0; i<x.NV; ++i)
		for(int j=0; j<x.NV; ++j)
			for(int k=0; k<x.NV; ++k)
				A(i,j)+=V(i,k)*VINV(k,j);
	
	temp = 0.0;
	for(int i=0; i<x.NV; ++i)
		for(int j=0; j<x.NV; ++j)
			for(int k=0; k<x.NV; ++k)
				temp(i,j)+=Pinv(i,k)*A(k,j);
	
	A = temp;
	
	for(int i = 0; i < x.NV; ++i)
		for(int j = 0; j < x.NV; ++j)
			fluxtemp(i) -= 0.5*A(i,j)*(ur(j)-ul(j));
	
	/* CHANGE BACK TO X,Y,Z COORDINATES */
	FLT temp2 = 1.0/(1.0-norm(2)*norm(1)-norm(1)*norm(0)-norm(0)*norm(2))/(1.0+norm(2)*norm(1)+norm(1)*norm(0)+norm(0)*norm(2));
    vec1(0) = norm(0); vec1(1) = temp2*(norm(1)*norm(0)-norm(2)*norm(2)); vec1(2) = temp2*(norm(0)*norm(1)*norm(2)-norm(1)*norm(1)*norm(1)-norm(2)*norm(2)*norm(1)+norm(0)*norm(0)*norm(2));
    vec2(0) = norm(1); vec2(1) = temp2*(norm(2)*norm(1)-norm(0)*norm(0)); vec2(2) = temp2*(norm(0)*norm(1)*norm(2)-norm(2)*norm(2)*norm(2)-norm(0)*norm(0)*norm(2)+norm(1)*norm(1)*norm(0));
    vec3(0) = norm(2); vec3(1) = temp2*(norm(2)*norm(0)-norm(1)*norm(1)); vec3(2) = temp2*(norm(0)*norm(1)*norm(2)-norm(0)*norm(0)*norm(0)-norm(1)*norm(1)*norm(0)+norm(2)*norm(2)*norm(1));
	
	flx(0) = fluxtemp(0);
	flx(1) = fluxtemp(1)*vec1(0) + fluxtemp(2)*vec1(1) + fluxtemp(3)*vec1(2);
	flx(2) = fluxtemp(1)*vec2(0) + fluxtemp(2)*vec2(1) + fluxtemp(3)*vec2(2);
	flx(3) = fluxtemp(1)*vec3(0) + fluxtemp(2)*vec3(1) + fluxtemp(3)*vec3(2);
	flx(4) = fluxtemp(4);
	
	flx *= mag;
	
	return;
}


void applied_stress::init(input_map& inmap) {
	std::string keyword;
	std::ostringstream nstr;

	neumann::init(inmap);
	
	stress.resize(tet_mesh::ND+1);
	for(int n=0;n<tet_mesh::ND+1;++n) {
		nstr.str("");
		nstr << base.idprefix << "_stress" << n << std::flush;
		if (inmap.find(nstr.str()) != inmap.end()) {
			stress(n).init(inmap,nstr.str());
		}
		else {
			*x.gbl->log << "couldn't find stress function " << nstr.str() << '\n';
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}
	}
	
	return;
}


void applied_stress::flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm, Array<FLT,1>& flx) {
	Array<FLT,1> v(3);
	FLT pr = ibc->f(0,xpt,x.gbl->time);
	v(0) = u(1)-mv(0);
	v(1) = u(2)-mv(1);
	v(2) = u(3)-mv(2);
	FLT RT = u(4);
	FLT rho = (pr+x.hp_cns_gbl->atm_pressure)/RT;
	
	/* CONTINUITY */
	flx(0) = rho*(v(0)*norm(0)+v(1)*norm(1)+v(2)*norm(2));
	
	FLT length = sqrt(norm(0)*norm(0)+norm(1)*norm(1)+norm(2)*norm(2));

	/* XYZ MOMENTUM */
	for (int n=1;n<tet_mesh::ND+1;++n)
		flx(n) = flx(0)*v(n-1) + pr*norm(n-1)-stress(n-1).Eval(xpt,x.gbl->time)*length;
	
	/* ENERGY EQUATION */
	FLT h = x.hp_cns_gbl->gamma/(x.hp_cns_gbl->gamma-1.0)*RT +0.5*(v(0)*v(0)+v(1)*v(1)+v(2)*v(2));				
	flx(x.NV-1) = h*flx(0)-stress(3).Eval(xpt,x.gbl->time)*length;
	
	return;
}


void pure_slip::flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm, Array<FLT,1>& flx) {
	
	//FLT pr = ibc->f(0,xpt,x.gbl->time);
	FLT pr = u(0);

	/* CONTINUITY */
	flx(0) = 0.0;
	
	/* XYZ MOMENTUM */
	for (int n=1;n<tet_mesh::ND+1;++n)
		flx(n) = pr*norm(n-1);
	
	/* ENERGY EQUATION */
	flx(x.NV-1) = 0.0;
	
	return;
}


//void pure_slip::vdirichlet() {			
//	for(int j=0;j<base.npnt;++j) {
//		int v0 = base.pnt(j).gindx;		
//		for (int n=0; n<ndirichlets; ++n) 
//			x.hp_gbl->res.v(v0,dirichlets(n)) = 0.0;		
//	}
//}
//
//void pure_slip::edirichlet() {
//	if (basis::tet(x.log2p).em > 0) {
//		for(int j=0;j<base.nseg;++j) {
//			int sind = base.seg(j).gindx;
//			for (int n=0; n<ndirichlets; ++n) 
//				x.hp_gbl->res.e(sind,Range::all(),dirichlets(n)) = 0.0;
//		}
//	}
//}		
//
//void pure_slip::fdirichlet() {
//	if (basis::tet(x.log2p).fm > 0) {
//		for(int j=0;j<base.ntri;++j) {
//			int find = base.tri(j).gindx;
//			for (int n=0; n<ndirichlets; ++n) 
//				x.hp_gbl->res.f(find,Range::all(),dirichlets(n)) = 0.0;
//		}
//	}
//}


/* ---------------------- */
/* edge boundary routines */
/* ---------------------- */


void neumann_edge::rsdl(int stage){
	int sind = -1;

	for(int i=0;i<base.nseg;++i){
		sind = base.seg(i).gindx;
		
		x.ugtouht1d(sind);
		
		element_rsdl(sind,stage);
		
		for(int n=0;n<x.NV;++n)
			x.hp_gbl->res.v(x.seg(sind).pnt(0),n) += x.lf(n)(0);
		
		for(int n=0;n<x.NV;++n)
			x.hp_gbl->res.v(x.seg(sind).pnt(1),n) += x.lf(n)(1);
		
		int indx = 3;
		
		for(int k=0;k<basis::tet(x.log2p).em;++k) {
			for(int n=0;n<x.NV;++n)
				x.hp_gbl->res.e(sind,k,n) += x.lf(n)(indx);
			++indx;
		}	
	}
	
	return;
}

void neumann_edge::element_rsdl(int eind,int stage) {
	int j,n;
	TinyVector<FLT,3> pt,mvel,nrm,vec1,vec2;
	Array<FLT,1> u(x.NV),flx(x.NV);
		
	x.lf = 0.0;
	
	x.crdtocht1d(eind);
	for(n=0;n<tet_mesh::ND;++n)
		basis::tet(x.log2p).proj1d(&x.cht(n)(0),&x.crd1d(n)(0),&x.dcrd1d(n)(0));
	
	for(n=0;n<x.NV;++n)
		basis::tet(x.log2p).proj1d(&x.uht(n)(0),&x.u1d(n)(0));
	
	for(j=0;j<basis::tet(x.log2p).gpx;++j) {
		
		for(n=0;n<3;++n)
			vec1(n)=x.dcrd1d(n)(j);
		
		vec2(0)=vec1(1);
		vec2(1)=vec1(2);
		vec2(2)=vec1(0);
		nrm=cross(vec1,vec2);
		
		for(n=0;n<tet_mesh::ND;++n) {
			pt(n) = x.crd1d(n)(j);
			mvel(n) = 0.0;//x.gbl->bd(0)*(x.crd2d(n)(j)(k) -dxdt(x.log2p,j)(n)(j)(k));
		}
		
		for(n=0;n<x.NV;++n)
			u(n) = x.u1d(n)(j);
		
		flux(u,pt,mvel,nrm,flx);
		
		for(n=0;n<x.NV;++n)
			x.res1d(n)(j) = flx(n);			
		
	}
	
	for(n=0;n<x.NV;++n)
		basis::tet(x.log2p).intgrt1d(&x.lf(n)(0),&x.res1d(n)(0));
	
	return;
	
	
}


void neumann_edge::flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm, Array<FLT,1>& flx) {
	
	Array<FLT,1> v(3);
	FLT pr = ibc->f(0,xpt,x.gbl->time);
	v(0) = u(1)-mv(0);
	v(1) = u(2)-mv(1);
	v(2) = u(3)-mv(2);
	FLT RT = u(4);
	FLT rho = (pr+x.hp_cns_gbl->atm_pressure)/RT;
	
	/* CONTINUITY */
	flx(0) = rho*(v(0)*norm(0)+v(1)*norm(1)+v(2)*norm(2));
	
	/* XYZ MOMENTUM */
	for (int n=1;n<tet_mesh::ND+1;++n)
		flx(n) = flx(0)*v(n-1) + pr*norm(n-1);
	
	/* ENERGY EQUATION */
	FLT h = x.hp_cns_gbl->gamma/(x.hp_cns_gbl->gamma-1.0)*RT +0.5*(v(0)*v(0)+v(1)*v(1)+v(2)*v(2));				
	flx(x.NV-1) = h*flx(0);
	
	return;
}

void inflow_edge::flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm,  Array<FLT,1>& flx) {
	
	Array<FLT,1> v(3);
	FLT pr = u(0);
	v(0) = ibc->f(1,xpt,x.gbl->time)-mv(0);
	v(1) = ibc->f(2,xpt,x.gbl->time)-mv(1);
	v(2) = ibc->f(3,xpt,x.gbl->time)-mv(2);
	FLT RT = ibc->f(4,xpt,x.gbl->time);
	FLT rho = (pr+x.hp_cns_gbl->atm_pressure)/RT;
	
	/* CONTINUITY */
	flx(0) = rho*(v(0)*norm(0)+v(1)*norm(1)+v(2)*norm(2));
	
	/* XYZ MOMENTUM */
	for (int n=1;n<tet_mesh::ND+1;++n)
		flx(n) = flx(0)*v(n-1) + pr*norm(n-1);
	
	/* ENERGY EQUATION */
	FLT h = x.hp_cns_gbl->gamma/(x.hp_cns_gbl->gamma-1.0)*RT +0.5*(v(0)*v(0)+v(1)*v(1)+v(2)*v(2));				
	flx(x.NV-1) = h*flx(0);
	
	return;
}

void inflow_edge::vdirichlet3d() {	
	int v0=-1,sind=-1;
	
	for(int j=0;j<base.nseg;++j) {
		sind = base.seg(j).gindx;
		v0 = x.seg(sind).pnt(0);
		for(int n=0; n<ndirichlets; ++n) 
			x.hp_gbl->res.v(v0,dirichlets(n)) = 0.0;		
	}
	
	v0 = x.seg(sind).pnt(1);
	for(int n=0; n<ndirichlets; ++n) 
		x.hp_gbl->res.v(v0,dirichlets(n)) = 0.0;
	
}

void inflow_edge::edirichlet3d() {
	//if (basis::tet(x.log2p).em > 0) {
		for(int j=0;j<base.nseg;++j) {
			int sind = base.seg(j).gindx;
			for(int n=0; n<ndirichlets; ++n) 
				x.hp_gbl->res.e(sind,Range::all(),dirichlets(n)) = 0.0;
		}
	//}
}		



void adiabatic_edge::flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm,  Array<FLT,1>& flx) {
	
	Array<FLT,1> v(3);
	FLT pr = u(0)+x.hp_cns_gbl->atm_pressure;
	v(0) = ibc->f(1,xpt,x.gbl->time);
	v(1) = ibc->f(2,xpt,x.gbl->time);
	v(2) = ibc->f(3,xpt,x.gbl->time);
	FLT RT = u(4);
	FLT rho = (pr+x.hp_cns_gbl->atm_pressure)/RT;
	
	/* CONTINUITY */
	flx(0) = rho*(v(0)*norm(0)+v(1)*norm(1)+v(2)*norm(2));
	
	/* XYZ MOMENTUM */
	for (int n=1;n<tet_mesh::ND+1;++n)
		flx(n) = flx(0)*v(n-1) + pr*norm(n-1);
	
	/* ENERGY EQUATION */
	FLT h = x.hp_cns_gbl->gamma/(x.hp_cns_gbl->gamma-1.0)*RT +0.5*(v(0)*v(0)+v(1)*v(1)+v(2)*v(2));				
	flx(x.NV-1) = h*flx(0);
	
	return;
}

void adiabatic_edge::vdirichlet3d() {	
	int v0=-1,sind=-1;
	
	for(int j=0;j<base.nseg;++j) {
		sind = base.seg(j).gindx;
		v0 = x.seg(sind).pnt(0);
		for(int n=0; n<ndirichlets; ++n) 
			x.hp_gbl->res.v(v0,dirichlets(n)) = 0.0;
	}	
	
	v0 = x.seg(sind).pnt(1);
	for(int n=0; n<ndirichlets; ++n) 
		x.hp_gbl->res.v(v0,dirichlets(n)) = 0.0;
}

void adiabatic_edge::edirichlet3d() {
	if (basis::tet(x.log2p).em > 0) {
		for(int j=0;j<base.nseg;++j) {
			int sind = base.seg(j).gindx;
			for(int n=0; n<ndirichlets; ++n) 
				x.hp_gbl->res.e(sind,Range::all(),dirichlets(n)) = 0.0;
		}
	}
}	



//void inflow_pt::tadvance() {
//	
//	hp_vrtx_bdry::tadvance(); 
//	
//	for(int n=1;n<x.NV;++n)
//		x.ug.v(base.pnt,n) = x.hp_gbl->ibc->f(n,x.pnts(base.pnt),x.gbl->time);  
//	
//	return;
//
//}



