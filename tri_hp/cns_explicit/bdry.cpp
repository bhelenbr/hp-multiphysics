#include "bdry_cns_explicit.h"
#include <myblas.h>

/*************************************************/
/* SET DIRICHLET BOUNDARY VALUES & FLUXES ********/
/* (THINGS THAT ARE INDEPENDENT OF THE SOLUTION) */
/* BUT NOT INDEPENDENT OF TIME/MESH POSITION     */
/*************************************************/
using namespace bdry_cns_explicit;


void applied_stress::init(input_map& inmap,void* gbl_in) {
	std::string keyword;
	std::ostringstream nstr;

	generic::init(inmap,gbl_in);

	stress.resize(x.NV-1);

	for(int n=0;n<x.NV-1;++n) {
		nstr.str("");
		nstr << base.idprefix << "_stress" << n << std::flush;
		if (inmap.find(nstr.str()) != inmap.end()) {
			stress(n).init(inmap,nstr.str());
		}
		else {
			*x.gbl->log << "couldn't find stress function " << nstr.str() << std::endl;
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}
	}

	return;
}

void characteristic::flux(Array<FLT,1>& cvu, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, Array<FLT,1>& flx) {	

	TinyVector<FLT,4> lambda,Rl,Rr,ub,Roe,fluxtemp;
	Array<FLT,2> A(x.NV,x.NV),V(x.NV,x.NV),VINV(x.NV,x.NV),temp(x.NV,x.NV);
	Array<FLT,1> Aeigs(x.NV);
	
	FLT mag = sqrt(norm(0)*norm(0) + norm(1)*norm(1));
	norm /= mag;

	/* Left */
	/* Rotate Coordinate System */
	FLT ul =  cvu(1)*norm(0) +cvu(2)*norm(1);
	FLT vl = -cvu(1)*norm(1) +cvu(2)*norm(0);
	cvu(1) = ul, cvu(2) = vl;

	/* Roe Variables */
	Rl(0) = sqrt(cvu(0));
	Rl(1) = ul/Rl(0);
	Rl(2) = vl/Rl(0);
	Rl(3) = cvu(3)/Rl(0);	
	
	/* Right */
	for(int n=0;n<x.NV;++n)
		ub(n) = ibc->f(n,xpt,x.gbl->time);

	/* Rotate Coordinate System */
	FLT ur =  ub(1)*norm(0) +ub(2)*norm(1);
	FLT vr = -ub(1)*norm(1) +ub(2)*norm(0);
	ub(1) = ur, ub(2) = vr;

	/* Roe Variables */
	Rr(0) = sqrt(ub(0));
	Rr(1) = ur/Rr(0);
	Rr(2) = vr/Rr(0);
	Rr(3) = ub(3)/Rr(0);	
	
	/* Average Roe Variables */
	Roe = 0.5*(Rl+Rr);

	/* Calculate u,v,c Variables */
	FLT u = Roe(1)/Roe(0);
	FLT v = Roe(2)/Roe(0);
	FLT ke = 0.5*(u*u+v*v);
	FLT E = Roe(3)/Roe(0);
	FLT gam = x.gbl->gamma;
	FLT gm1 = gam-1.0;
	FLT RT = gm1*(E-ke);
	FLT c2 = gam*RT;
	FLT c = sqrt(c2);

//	/* df/dw in u,v,c Variables */
//	A = 0.0,                                       1.0,                                                                         0.0,            0.0,
//	    gam*ke-1.5*u*u-0.5*v*v,                3.0*u-u*gam,                                                               -v*(gam-1.0),  gam-1,
//	    -u*v,                                    v,                                                                          u,             0.0,
//	    u*(-3*gam*ke+2*ke-c2+gam*gam*ke)/(gam-1), -(1.5*u*u-2.5*u*u*gam-0.5*v*v*gam+0.5*v*v-c2+u*u*gam*gam)/(gam-1), -u*(gam-1)*v, u*gam;

	
	/* eigenvectors of df/dw */
	V = 0.0, 1.0,           1.0,              1.0,
	    0.0, u,             u+c,              u-c,
	    1.0, 0.0,           v,                v,
	    v,   0.5*(u*u-v*v), gam*E-gm1*ke+u*c, gam*E-gm1*ke-u*c;
	
	/* eigenvalues of df/dw */
	Aeigs = u,u,u+c,u-c;

		/* inverse of eigenvectors of df/dw */
	VINV = -v*ke*gm1,         v*u*gm1,         c2+gm1*v*v, -v*gm1,
		   -gm1*ke+c2,        u*gm1,           v*gm1,      -gm1,
		   0.5*(-u*c+gm1*ke), -0.5*(-c+u*gm1), -0.5*v*gm1, 0.5*gm1,
		   0.5*(u*c+gm1*ke),  -0.5*(c+u*gm1),  -0.5*v*gm1, 0.5*gm1;
	
	VINV /= c2;
	
	for(int i=0; i < x.NV; ++i)
		for(int j=0; j < x.NV; ++j)
			temp(i,j) = Aeigs(i)*VINV(i,j);
	
	A = 0.0;
	for(int i=0; i<x.NV; ++i)
		for(int j=0; j<x.NV; ++j)
			for(int k=0; k<x.NV; ++k)
				A(i,j)+=V(i,k)*temp(k,j);
	
	fluxtemp = 0.0;
	
	for(int i = 0; i < x.NV; ++i)
		for(int j = 0; j < x.NV; ++j)
			fluxtemp(i) += 0.5*A(i,j)*(ub(j)+cvu(j));
	
//	/* or this way */
//	FLT rho = Roe(0)*Roe(0);
//	FLT pr = rho*RT;
//	
//	fluxtemp(0) = rho*u;
//	fluxtemp(1) = rho*u*u+pr;
//	fluxtemp(2) = rho*u*v;
//	fluxtemp(3) = rho*u*(RT*gam/gm1+ke);
	
	Aeigs = fabs(u),fabs(u),fabs(u+c),fabs(u-c);

	//matrix_absolute_value(A);

	for(int i=0; i < x.NV; ++i)
		for(int j=0; j < x.NV; ++j)
			temp(i,j) = Aeigs(i)*VINV(i,j);
	
	A = 0.0;
	for(int i=0; i<x.NV; ++i)
		for(int j=0; j<x.NV; ++j)
			for(int k=0; k<x.NV; ++k)
				A(i,j)+=V(i,k)*temp(k,j);
	
	for(int i = 0; i < x.NV; ++i)
		for(int j = 0; j < x.NV; ++j)
			fluxtemp(i) -= 0.5*A(i,j)*(ub(j)-cvu(j));

	/* CHANGE BACK TO X,Y COORDINATES */
	flx(0) = fluxtemp(0);
	flx(1) = fluxtemp(1)*norm(0) - fluxtemp(2)*norm(1);
	flx(2) = fluxtemp(1)*norm(1) + fluxtemp(2)*norm(0);
	flx(3) = fluxtemp(3);
	
	flx *= mag;
	
	return;
}
