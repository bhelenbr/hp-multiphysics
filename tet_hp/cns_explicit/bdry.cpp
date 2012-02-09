#include "bdry_cns_explicit.h"
#include <myblas.h>
#include<blitz/tinyvec-et.h>
/*************************************************/
/* SET DIRICHLET BOUNDARY VALUES & FLUXES ********/
/* (THINGS THAT ARE INDEPENDENT OF THE SOLUTION) */
/* BUT NOT INDEPENDENT OF TIME/MESH POSITION     */
/*************************************************/
using namespace bdry_cns_explicit;

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


void neumann::flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm, Array<FLT,1>& flx) {
	
	/* CONTINUITY */
	flx(0) = ibc->f(0, xpt, x.gbl->time)/u(x.NV-1)*((u(1) -mv(0))*norm(0) +(u(2) -mv(1))*norm(1)+(u(3) -mv(2))*norm(2));
	
	/* X&Y MOMENTUM */
#ifdef INERTIALESS
	for (int n=1;n<tet_mesh::ND+1;++n)
		flx(n) = ibc->f(0, xpt, x.gbl->time)*norm(n-1);
#else
	for (int n=1;n<tet_mesh::ND+1;++n)
		flx(n) = flx(0)*u(n) +ibc->f(0, xpt, x.gbl->time)*norm(n-1);
#endif
	
	/* ENERGY EQUATION */
	double h = x.gbl->gamma/(x.gbl->gamma-1.0)*u(x.NV-1) +0.5*(u(1)*u(1)+u(2)*u(2)+u(3)*u(3));				
	flx(x.NV-1) = h*flx(0);
	
	return;
}

void inflow::flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm,  Array<FLT,1>& flx) {
	
	for(int i = 0; i < x.NV; ++i)
		flx(i) = 0.0;
	return;
	
	
	Array<FLT,1> v(3);
	FLT pr = (x.gbl->gamma-1.0)*(u(4)-0.5/u(0)*(u(1)*u(1)+u(2)*u(2)+u(3)*u(3)));
	v(0) = ibc->f(1,xpt,x.gbl->time)/ibc->f(0,xpt,x.gbl->time)-mv(0);
	v(1) = ibc->f(2,xpt,x.gbl->time)/ibc->f(0,xpt,x.gbl->time)-mv(1);
	v(2) = ibc->f(3,xpt,x.gbl->time)/ibc->f(0,xpt,x.gbl->time)-mv(2);
	FLT RT = (x.gbl->gamma-1.0)*(ibc->f(4,xpt,x.gbl->time)-0.5/ibc->f(0,xpt,x.gbl->time)*(ibc->f(1,xpt,x.gbl->time)*ibc->f(1,xpt,x.gbl->time)+ibc->f(2,xpt,x.gbl->time)*ibc->f(2,xpt,x.gbl->time)+ibc->f(3,xpt,x.gbl->time)*ibc->f(3,xpt,x.gbl->time)))/ibc->f(0,xpt,x.gbl->time);
	FLT rho = pr/RT;
	
	/* CONTINUITY */
	flx(0) = rho*(v(0)*norm(0)+v(1)*norm(1)+v(2)*norm(2));
	
	/* XYZ MOMENTUM */
	for (int n=1;n<tet_mesh::ND+1;++n)
		flx(n) = flx(0)*v(n-1) + pr*norm(n-1);
	
	/* ENERGY EQUATION */
	FLT h = x.gbl->gamma/(x.gbl->gamma-1.0)*RT +0.5*(v(0)*v(0)+v(1)*v(1)+v(2)*v(2));				
	flx(x.NV-1) = h*flx(0);
	
	
	return;
}

void inflow::vdirichlet() {	
	for(int j=0;j<base.npnt;++j) {
		int v0 = base.pnt(j).gindx;
		
		for (int n=0; n<ndirichlets; ++n) 
			x.gbl->res.v(v0,dirichlets(n)) = 0.0;
		
	}
}

void inflow::edirichlet() {
	if (basis::tet(x.log2p).em > 0) {
		for(int j=0;j<base.nseg;++j) {
			int sind = base.seg(j).gindx;
			for (int n=0; n<ndirichlets; ++n) 
				x.gbl->res.e(sind,Range::all(),dirichlets(n)) = 0.0;
		}
	}
}		

void inflow::fdirichlet() {
	if (basis::tet(x.log2p).fm > 0) {
		for(int j=0;j<base.ntri;++j) {
			int find = base.tri(j).gindx;
			for (int n=0; n<ndirichlets; ++n) 
				x.gbl->res.f(find,Range::all(),dirichlets(n)) = 0.0;
		}
	}
}	

void inflow::update(int stage) {
	int j,k,m,n,v0,v1,sind;
	double KE,rho,u,v,w,RT;
	TinyVector<FLT,tet_mesh::ND> pt;
	
	TinyVector<double,MXTM> ucoef;
		
	j = 0;
	do {
		sind = base.seg(j).gindx;
		v0 = x.seg(sind).pnt(0);
		
		u = ibc->f(1, x.pnts(v0), x.gbl->time)/ibc->f(0, x.pnts(v0), x.gbl->time);
		v = ibc->f(2, x.pnts(v0), x.gbl->time)/ibc->f(0, x.pnts(v0), x.gbl->time);
		w = ibc->f(3, x.pnts(v0), x.gbl->time)/ibc->f(0, x.pnts(v0), x.gbl->time);
		KE = 0.5*(u*u+v*v+w*w);
		RT = (ibc->f(4, x.pnts(v0), x.gbl->time)/ibc->f(0, x.pnts(v0), x.gbl->time)-KE)*(x.gbl->gamma-1.0);
		
		rho = x.ug.v(v0,0);
		
		x.ug.v(v0,1) = rho*u;
		x.ug.v(v0,2) = rho*v;
		x.ug.v(v0,3) = rho*w;
		x.ug.v(v0,4) = rho*(RT/(x.gbl->gamma-1.0)+KE);
		
	} while (++j < base.nseg);
	v0 = x.seg(sind).pnt(1);
	
	u = ibc->f(1, x.pnts(v0), x.gbl->time)/ibc->f(0, x.pnts(v0), x.gbl->time);
	v = ibc->f(2, x.pnts(v0), x.gbl->time)/ibc->f(0, x.pnts(v0), x.gbl->time);
	w = ibc->f(3, x.pnts(v0), x.gbl->time)/ibc->f(0, x.pnts(v0), x.gbl->time);
	KE = 0.5*(u*u+v*v+w*w);
	RT = (ibc->f(4, x.pnts(v0), x.gbl->time)/ibc->f(0, x.pnts(v0), x.gbl->time)-KE)*(x.gbl->gamma-1.0);
	
	rho = x.ug.v(v0,0);
	
	x.ug.v(v0,1) = rho*u;
	x.ug.v(v0,2) = rho*v;
	x.ug.v(v0,3) = rho*w;
	x.ug.v(v0,4) = rho*(RT/(x.gbl->gamma-1.0)+KE);
	
	if(basis::tet(x.log2p).p > 2){
		*x.gbl->log << "update in boundary face only works for p=1,2" << endl;
		exit(3);
	}
	if(basis::tet(x.log2p).em){
		for(j=0;j<base.nseg;++j) {
			sind = base.seg(j).gindx;
			v0 = x.seg(sind).pnt(0);
			v1 = x.seg(sind).pnt(1);
			
			if (x.seg(sind).info < 0){
				for(n=0;n<x.ND;++n)
					basis::tet(x.log2p).proj1d(x.pnts(v0)(n),x.pnts(v1)(n),&x.crd1d(n)(0));
			}
			else {
				x.crdtocht1d(sind,0);
				for(n=0;n<x.ND;++n)
					basis::tet(x.log2p).proj1d(&x.cht(n)(0),&x.crd1d(n)(0));
			}
			
			/* take global coefficients and put into local vector */
			/*  only need rho */
			for (m=0; m<2; ++m) 
				ucoef(m) = x.ug.v(x.seg(sind).pnt(m),0);					
			
			for (m=0;m<basis::tet(x.log2p).em;++m) 
				ucoef(m+2) = x.ug.e(sind,m,0);					
			
			basis::tet(x.log2p).proj1d(&ucoef(0),&x.u1d(0)(0));
			
			for(n=1;n<x.NV;++n)
				basis::tet(x.log2p).proj1d(x.ug.v(v0,n),x.ug.v(v1,n),&x.u1d(n)(0));

			
			for(k=0;k<basis::tet(x.log2p).gpx; ++k) {
				pt(0) = x.crd1d(0)(k);
				pt(1) = x.crd1d(1)(k);
				pt(2) = x.crd1d(2)(k);
				
				u = ibc->f(1, pt, x.gbl->time)/ibc->f(0, pt, x.gbl->time);
				v = ibc->f(2, pt, x.gbl->time)/ibc->f(0, pt, x.gbl->time);
				w = ibc->f(3, pt, x.gbl->time)/ibc->f(0, pt, x.gbl->time);
				KE = 0.5*(u*u+v*v+w*w);
				RT = (ibc->f(4, pt, x.gbl->time)/ibc->f(0, pt, x.gbl->time)-KE)*(x.gbl->gamma-1.0);
				
				rho = x.u1d(0)(k);
				
				x.u1d(1)(k) -= rho*u;
				x.u1d(2)(k) -= rho*v;
				x.u1d(3)(k) -= rho*w;
				x.u1d(4)(k) -= rho*(RT/(x.gbl->gamma-1.0)+KE);
				
				
			}
			for(n=1;n<x.NV;++n)
				basis::tet(x.log2p).intgrt1d(&x.lf(n)(0),&x.u1d(n)(0));
			
			for(n=1;n<x.NV;++n) {
				for(m=0;m<basis::tet(x.log2p).em;++m){ 
					x.ug.e(sind,m,n) = -x.lf(n)(2+m)*basis::tet(x.log2p).diag1d(m);
				}
			}			
		}
	}					
}

void adiabatic::flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm,  Array<FLT,1>& flx) {
	
	FLT KE = 0.5*(pow(ibc->f(1, xpt, x.gbl->time)/ibc->f(0, xpt, x.gbl->time),2.0)
				  +pow(ibc->f(2, xpt, x.gbl->time)/ibc->f(0, xpt, x.gbl->time),2.0)
				  +pow(ibc->f(3, xpt, x.gbl->time)/ibc->f(0, xpt, x.gbl->time),2.0));
	
	FLT pr = (ibc->f(4, xpt, x.gbl->time)-ibc->f(0, xpt, x.gbl->time)*KE)*(x.gbl->gamma-1.0);
	
	FLT u1 = u(1)/u(0);
	FLT u2 = u(2)/u(0);
	FLT u3 = u(3)/u(0);
	FLT RT = (x.gbl->gamma-1.0)*(u(4)-0.5*u(0)*(u1*u1+u2*u2+u3*u3))/u(0);
	
	FLT rho = pr/RT;
	
	flx(0) = rho*(u1*norm(0)+u2*norm(1)+u3*norm(2));	

	/* XYZ MOMENTUM */
	for (int n=0;n<tet_mesh::ND;++n)
		flx(n+1) = 0.0;
	
	/* ENERGY EQUATION */
	FLT h = x.gbl->gamma/(x.gbl->gamma-1.0)*RT +0.5*(u1*u1+u2*u2+u3*u3);
	flx(4) = h*flx(0);
	
	return;
}

void adiabatic::vdirichlet() {			
	for(int j=0;j<base.npnt;++j) {
		int v0 = base.pnt(j).gindx;
		
		for (int n=0; n<ndirichlets; ++n) 
			x.gbl->res.v(v0,dirichlets(n)) = 0.0;
		
	}
}

void adiabatic::edirichlet() {
	if (basis::tet(x.log2p).em > 0) {
		for(int j=0;j<base.nseg;++j) {
			int sind = base.seg(j).gindx;
			for (int n=0; n<ndirichlets; ++n) 
				x.gbl->res.e(sind,Range::all(),dirichlets(n)) = 0.0;
		}
	}
}		

void adiabatic::fdirichlet() {
	if (basis::tet(x.log2p).fm > 0) {
		for(int j=0;j<base.ntri;++j) {
			int find = base.tri(j).gindx;
			for (int n=0; n<ndirichlets; ++n) 
				x.gbl->res.f(find,Range::all(),dirichlets(n)) = 0.0;
		}
	}
}		


void characteristic::flux(Array<FLT,1>& cvu, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm, Array<FLT,1>& flx) {	
	
	TinyVector<FLT,5> lambda,Rl,Rr,ul,ur,ub,Roe,fluxtemp;
	TinyVector<FLT,3> vec1,vec2,vec3,vecr;
	Array<FLT,2> A(x.NV,x.NV),V(x.NV,x.NV),VINV(x.NV,x.NV),temp(x.NV,x.NV);
	Array<FLT,1> eigs(x.NV);
	FLT gam = x.gbl->gamma;
	FLT gm1 = gam-1.0;

	FLT mag = sqrt(norm(0)*norm(0)+norm(1)*norm(1)+norm(2)*norm(2));
	norm /= mag;

	/* align u to the normal */
	vec1 = norm;
	/* find random vector different than norm */
	vecr(0) = vec1(1),	vecr(1) = vec1(2), vecr(2) = vec1(0);
	vec2 = cross(norm, vecr);
	vec3 = cross(norm, vec2);

	/* Left */
	/* Rotate Coordinate System so that U aligns with normal vector */
	ul(0) = cvu(0);
	ul(1) = cvu(1)*vec1(0)+cvu(2)*vec1(1)+cvu(3)*vec1(2);	
	ul(2) = cvu(1)*vec2(0)+cvu(2)*vec2(1)+cvu(3)*vec2(2);
	ul(3) = cvu(1)*vec3(0)+cvu(2)*vec3(1)+cvu(3)*vec3(2);
	ul(4) = cvu(4);

	/* Right */
	for(int n=0;n<x.NV;++n)
		ub(n) = ibc->f(n,xpt,x.gbl->time);

	/* Rotate Coordinate System */
	ur(0) = ub(0);
	ur(1) = ub(1)*vec1(0)+ub(2)*vec1(1)+ub(3)*vec1(2);	
	ur(2) = ub(1)*vec2(0)+ub(2)*vec2(1)+ub(3)*vec2(2);
	ur(3) = ub(1)*vec3(0)+ub(2)*vec3(1)+ub(3)*vec3(2);
	ur(4) = ub(4);

	/* Roe Variables from left */
	Rl(0) = sqrt(ul(0)); // sqrt(rho)
	Rl(1) = ul(1)/Rl(0); // sqrt(rho)*u
	Rl(2) = ul(2)/Rl(0); // sqrt(rho)*v
	Rl(3) = ul(3)/Rl(0); // sqrt(rho)*w
	Rl(4) = ul(4)/Rl(0); // sqrt(rho)*E
	
	/* Roe Variables from right */
	Rr(0) = sqrt(ur(0));
	Rr(1) = ur(1)/Rr(0);
	Rr(2) = ur(2)/Rr(0);
	Rr(3) = ur(3)/Rr(0);
	Rr(4) = ur(4)/Rr(0);	
	
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
	FLT c2 = gam*rt;
	FLT c = sqrt(c2);

	/* eigenvectors of df_x/dw */
	V = 0.0, 1.0,    0.0, 1.0,                     1.0,
		0.0, u,      0.0, u+c,                     u-c,
		0.0, 0.0,    1.0, v,                       v,
		1.0, 0.0,    0.0, w,                       w,
		w,   u*u-ke, v,   (gm1*u*c+gm1*ke+c2)/gm1, (-gm1*u*c+gm1*ke+c2)/gm1;
	
	eigs = u,u,u,u+c,u-c;

	VINV = -w*ke*gm1,        u*w*gm1,        v*w*gm1,    gm1*w*w+c2, -w*gm1,
		   c2-gm1*ke,        u*gm1,          v*gm1,      w*gm1,      -gm1,
		   -v*ke*gm1,        u*v*gm1,        gm1*v*v+c2, v*w*gm1,    -v*gm1,
		   0.5*(gm1*ke-u*c), -0.5*(u*gm1-c), -0.5*v*gm1, -0.5*w*gm1, 0.5*gm1,
		   0.5*(gm1*ke+u*c), -0.5*(u*gm1+c), -0.5*v*gm1, -0.5*w*gm1, 0.5*gm1;
	
	VINV /= c2;
	
	for(int i=0; i < x.NV; ++i)
		for(int j=0; j < x.NV; ++j)
			temp(i,j) = eigs(i)*VINV(i,j);
	
	A = 0.0;
	for(int i=0; i<x.NV; ++i)
		for(int j=0; j<x.NV; ++j)
			for(int k=0; k<x.NV; ++k)
				A(i,j)+=V(i,k)*temp(k,j);
	
	fluxtemp = 0.0;

	for(int i = 0; i < x.NV; ++i)
		for(int j = 0; j < x.NV; ++j)
			fluxtemp(i) += 0.5*A(i,j)*(ur(j)+ul(j));
	
	eigs = fabs(u),fabs(u),fabs(u),fabs(u+c),fabs(u-c);

	for(int i=0; i < x.NV; ++i)
		for(int j=0; j < x.NV; ++j)
			temp(i,j) = eigs(i)*VINV(i,j);
	
	A = 0.0;
	for(int i=0; i<x.NV; ++i)
		for(int j=0; j<x.NV; ++j)
			for(int k=0; k<x.NV; ++k)
				A(i,j)+=V(i,k)*temp(k,j);
	
	for(int i = 0; i < x.NV; ++i)
		for(int j = 0; j < x.NV; ++j)
			fluxtemp(i) -= 0.5*A(i,j)*(ur(j)-ul(j));
	
	/* CHANGE BACK TO X,Y,Z COORDINATES */
	FLT temp1 = 1.0/(1.0-norm(2)*norm(1)-norm(1)*norm(0)-norm(0)*norm(2))/(1.0+norm(2)*norm(1)+norm(1)*norm(0)+norm(0)*norm(2));
	vec1(0) = norm(0), vec1(1) = temp1*(norm(1)*norm(0)-norm(2)*norm(2)), vec1(2) = temp1*(norm(0)*norm(1)*norm(2)-norm(1)*norm(1)*norm(1)-norm(2)*norm(2)*norm(1)+norm(0)*norm(0)*norm(2));
	vec2(0) = norm(1), vec2(1) = temp1*(norm(2)*norm(1)-norm(0)*norm(0)), vec2(2) = temp1*(norm(0)*norm(1)*norm(2)-norm(2)*norm(2)*norm(2)-norm(0)*norm(0)*norm(2)+norm(1)*norm(1)*norm(0));
	vec3(0) = norm(2), vec3(1) = temp1*(norm(2)*norm(0)-norm(1)*norm(1)), vec3(2) = temp1*(norm(0)*norm(1)*norm(2)-norm(0)*norm(0)*norm(0)-norm(1)*norm(1)*norm(0)+norm(2)*norm(2)*norm(1));
	
	flx(0) = fluxtemp(0);
	flx(1) = fluxtemp(1)*vec1(0) + fluxtemp(2)*vec1(1) + fluxtemp(3)*vec1(2);
	flx(2) = fluxtemp(1)*vec2(0) + fluxtemp(2)*vec2(1) + fluxtemp(3)*vec2(2);
	flx(3) = fluxtemp(1)*vec3(0) + fluxtemp(2)*vec3(1) + fluxtemp(3)*vec3(2);
	flx(4) = fluxtemp(4);
	
	flx *= mag;
	
	return;
}


void applied_stress::init(input_map& inmap,void* gbl_in) {
	std::string keyword;
	std::ostringstream nstr;

	neumann::init(inmap,gbl_in);
	
	stress.resize(tet_mesh::ND+1);
	for(int n=0;n<tet_mesh::ND+1;++n) {
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


void applied_stress::flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm, Array<FLT,1>& flx) {
	
	
	FLT KE = 0.5*(pow(ibc->f(1, xpt, x.gbl->time)/ibc->f(0, xpt, x.gbl->time),2.0)
				  +pow(ibc->f(2, xpt, x.gbl->time)/ibc->f(0, xpt, x.gbl->time),2.0)
				  +pow(ibc->f(3, xpt, x.gbl->time)/ibc->f(0, xpt, x.gbl->time),2.0));
	
	FLT pr = (ibc->f(4, xpt, x.gbl->time)-ibc->f(0, xpt, x.gbl->time)*KE)*(x.gbl->gamma-1.0);
	
	FLT u1 = u(1)/u(0);
	FLT u2 = u(2)/u(0);
	FLT u3 = u(3)/u(0);
	FLT RT = (x.gbl->gamma-1.0)*(u(4)-0.5*u(0)*(u1*u1+u2*u2+u3*u3))/u(0);

	FLT rho = pr/RT;
	
	flx(0) = rho*(u1*norm(0)+u2*norm(1)+u3*norm(2));
	
	FLT length = sqrt(norm(0)*norm(0)+norm(1)*norm(1)+norm(2)*norm(2));
		
	/* XYZ MOMENTUM */
#ifdef INERTIALESS
	for (int n=0;n<tet_mesh::ND;++n)
		flx(n+1) = -stress(n).Eval(xpt,x.gbl->time)*length+pr*norm(n);
#else
	
	flx(1) = flx(0)*u1+pr*norm(0)-stress(0).Eval(xpt,x.gbl->time)*length;
	flx(2) = flx(0)*u2+pr*norm(1)-stress(1).Eval(xpt,x.gbl->time)*length;
	flx(3) = flx(0)*u3+pr*norm(2)-stress(2).Eval(xpt,x.gbl->time)*length;

#endif
	
	/* ENERGY EQUATION */
	FLT h = x.gbl->gamma/(x.gbl->gamma-1.0)*RT +0.5*(u1*u1+u2*u2+u3*u3);
	flx(4) = h*flx(0)-stress(3).Eval(xpt,x.gbl->time)*length;		
	
	return;
}

void specified_flux::init(input_map& inmap,void* gbl_in) {
	std::string keyword;
	std::ostringstream nstr;
	
	neumann::init(inmap,gbl_in);
	
	stress.resize(x.NV);
	for(int n=0;n<x.NV;++n) {
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


void specified_flux::flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm, Array<FLT,1>& flx) {
	
	
	FLT length = sqrt(norm(0)*norm(0)+norm(1)*norm(1)+norm(2)*norm(2));
	
	for (int n=0; n<x.NV; ++n) 
		flx(n) = stress(n).Eval(xpt,x.gbl->time)*length;
	
	
	return;
}


