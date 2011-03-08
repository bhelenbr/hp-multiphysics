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
	
	/* CONTINUITY */
	flx(0) = ibc->f(0, xpt, x.gbl->time)/u(x.NV-1)*((u(1) -mv(0))*norm(0) +(u(2) -mv(1))*norm(1)+(u(3) -mv(2))*norm(2));
	
	/* EVERYTHING ELSE DOESN'T MATTER */
	for (int n=1;n<x.NV;++n)
		flx(n) = 0.0;
	
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

void adiabatic::flux(Array<FLT,1>& u, TinyVector<FLT,tet_mesh::ND> xpt, TinyVector<FLT,tet_mesh::ND> mv, TinyVector<FLT,tet_mesh::ND> norm,  Array<FLT,1>& flx) {
	
	/* CONTINUITY */
	flx(0) = ibc->f(0, xpt, x.gbl->time)/u(x.NV-1)*((u(1) -mv(0))*norm(0) +(u(2) -mv(1))*norm(1)+(u(3) -mv(2))*norm(2));
	
	/* EVERYTHING ELSE DOESN'T MATTER */
	for (int n=1;n<x.NV;++n)
		flx(n) = 0.0;
	
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
	
	/* CONTINUITY */
	flx(0) = ibc->f(0, xpt, x.gbl->time)/u(x.NV-1)*((u(1) -mv(0))*norm(0) +(u(2) -mv(1))*norm(1)+(u(3) -mv(2))*norm(2));
	
	FLT length = sqrt(norm(0)*norm(0)+norm(1)*norm(1)+norm(2)*norm(2));
	/* XYZ MOMENTUM */
#ifdef INERTIALESS
	for (int n=0;n<tet_mesh::ND;++n)
		flx(n+1) = -stress(n).Eval(xpt,x.gbl->time)*length +ibc->f(0, xpt, x.gbl->time)*norm(n);
#else
	for (int n=0;n<tet_mesh::ND;++n)
		flx(n+1) = flx(0)*u(n+1) -stress(n).Eval(xpt,x.gbl->time)*length +ibc->f(0, xpt, x.gbl->time)*norm(n);
#endif
	
	/* ENERGY EQUATION */
	double h = x.gbl->gamma/(x.gbl->gamma-1.0)*u(x.NV-1) +0.5*(u(1)*u(1)+u(2)*u(2)+u(3)*u(3));
	flx(x.NV-1) = h*flx(0)-stress(3).Eval(xpt,x.gbl->time)*length;		
	
	return;
}



