#include "bdry_cns_explicit.h"
#include <myblas.h>

/*************************************************/
/* SET DIRICHLET BOUNDARY VALUES & FLUXES ********/
/* (THINGS THAT ARE INDEPENDENT OF THE SOLUTION) */
/* BUT NOT INDEPENDENT OF TIME/MESH POSITION     */
/*************************************************/
using namespace bdry_cns_explicit;

void generic::output(std::ostream& fout, tri_hp::filetype typ,int tlvl) {
//	int i,m,n,ind,sind,tind,seg;
//	FLT visc[tri_mesh::ND+1][tri_mesh::ND+1][tri_mesh::ND][tri_mesh::ND];
//	TinyVector<FLT,tri_mesh::ND> norm, mvel;
//	FLT convect,jcb;
//
//	switch(typ) {
//		case(tri_hp::text): case(tri_hp::binary): {
//			hp_edge_bdry::output(fout,typ,tlvl);
//			break;
//		}
//		case(tri_hp::tecplot): {
//			if (!report_flag) return;
//
//			conv_flux = 0.0;
//			diff_flux = 0.0;
//			moment = 0.0;
//			circumference = 0.0;
//			circulation = 0.0;
//#ifdef L2_ERROR
//			FLT l2error = 0.0;
//			TinyVector<FLT,2> xpt;
//#endif
//
//			for(ind=0; ind < base.nseg; ++ind) {
//				sind = base.seg(ind);
//				tind = x.seg(sind).tri(0);        
//
//				for(seg=0;seg<3;++seg)
//					if (x.tri(tind).seg(seg) == sind) break;
//				assert(seg != 3);
//
//				x.crdtocht(tind);
//				for(m=basis::tri(x.log2p)->bm();m<basis::tri(x.log2p)->tm();++m)
//					for(n=0;n<tri_mesh::ND;++n)
//						x.cht(n,m) = 0.0;
//
//				for(n=0;n<tri_mesh::ND;++n)
//					basis::tri(x.log2p)->proj_side(seg,&x.cht(n,0), &x.crd(n)(0,0), &x.dcrd(n,0)(0,0), &x.dcrd(n,1)(0,0));
//
//				x.ugtouht(tind);
//				for(n=0;n<x.NV;++n)
//					basis::tri(x.log2p)->proj_side(seg,&x.uht(n)(0),&x.u(n)(0,0),&x.du(n,0)(0,0),&x.du(n,1)(0,0));
//
//				for (i=0;i<basis::tri(x.log2p)->gpx();++i) {
//					jcb =  basis::tri(x.log2p)->wtx(i)*RAD(x.crd(0)(0,i))*sqrt(x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i) +x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i));
//					circumference += jcb;
//
//					x.cjcb(0,i) = x.gbl->mu*RAD(x.crd(0)(0,i))/(x.dcrd(0,0)(0,i)*x.dcrd(1,1)(0,i) -x.dcrd(1,0)(0,i)*x.dcrd(0,1)(0,i));
//
//					/* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
//					/* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
//					visc[0][0][0][0] =  x.cjcb(0,i)*(2.*x.dcrd(1,1)(0,i)*x.dcrd(1,1)(0,i) +x.dcrd(0,1)(0,i)*x.dcrd(0,1)(0,i));
//					visc[0][0][1][1] =  x.cjcb(0,i)*(2.*x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i));
//					visc[0][0][0][1] = -x.cjcb(0,i)*(2.*x.dcrd(1,1)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,1)(0,i)*x.dcrd(0,0)(0,i));
//#define             viscI0II0II1II0I visc[0][0][0][1]
//
//					visc[1][1][0][0] =  x.cjcb(0,i)*(x.dcrd(1,1)(0,i)*x.dcrd(1,1)(0,i) +2.*x.dcrd(0,1)(0,i)*x.dcrd(0,1)(0,i));
//					visc[1][1][1][1] =  x.cjcb(0,i)*(x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i) +2.*x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i));
//					visc[1][1][0][1] = -x.cjcb(0,i)*(x.dcrd(1,1)(0,i)*x.dcrd(1,0)(0,i) +2.*x.dcrd(0,1)(0,i)*x.dcrd(0,0)(0,i));
//#define             viscI1II1II1II0I visc[1][1][0][1]
//
//					visc[0][1][0][0] = -x.cjcb(0,i)*x.dcrd(0,1)(0,i)*x.dcrd(1,1)(0,i);
//					visc[0][1][1][1] = -x.cjcb(0,i)*x.dcrd(0,0)(0,i)*x.dcrd(1,0)(0,i);
//					visc[0][1][0][1] =  x.cjcb(0,i)*x.dcrd(0,1)(0,i)*x.dcrd(1,0)(0,i);
//					visc[0][1][1][0] =  x.cjcb(0,i)*x.dcrd(0,0)(0,i)*x.dcrd(1,1)(0,i);
//
//					/* OTHER SYMMETRIES     */                
//#define             viscI1II0II0II0I visc[0][1][0][0]
//#define             viscI1II0II1II1I visc[0][1][1][1]
//#define             viscI1II0II0II1I visc[0][1][1][0]
//#define             viscI1II0II1II0I visc[0][1][0][1]
//
//					/* DIFFUSIVE FLUXES ( FOR EXTRA VARIABLES) */
//					visc[2][2][1][0] =  x.cjcb(0,i)*(x.dcrd(1,1)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,1)(0,i)*x.dcrd(0,0)(0,i));
//					visc[2][2][1][1] = -x.cjcb(0,i)*(x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i));
//
//					for (n=tri_mesh::ND;n<x.NV-1;++n) 
//						diff_flux(n) -= x.gbl->D(n)/x.gbl->mu*basis::tri(x.log2p)->wtx(i)*(-visc[2][2][1][0]*x.du(2,0)(0,i) -visc[2][2][1][1]*x.du(2,1)(0,i));
//
//					diff_flux(0) -=    basis::tri(x.log2p)->wtx(i)*(-x.u(2)(0,i)*RAD(x.crd(0)(0,i))*x.dcrd(1,0)(0,i) 
//									-viscI0II0II1II0I*x.du(0,0)(0,i) -visc[0][1][1][0]*x.du(1,0)(0,i)
//									-visc[0][0][1][1]*x.du(0,1)(0,i) -visc[0][1][1][1]*x.du(1,1)(0,i));															
//					diff_flux(1) -=    basis::tri(x.log2p)->wtx(i)*( x.u(2)(0,i)*RAD(x.crd(0)(0,i))*x.dcrd(0,0)(0,i)
//									-viscI1II0II1II0I*x.du(0,0)(0,i) -viscI1II1II1II0I*x.du(1,0)(0,i)
//									-viscI1II0II1II1I*x.du(0,1)(0,i) -visc[1][1][1][1]*x.du(1,1)(0,i));
//
//
//					norm(0) = x.dcrd(1,0)(0,i);
//					norm(1) = -x.dcrd(0,0)(0,i);                
//					for(n=0;n<tri_mesh::ND;++n) {
//						mvel(n) = x.gbl->bd(0)*(x.crd(n)(0,i) -dxdt(x.log2p,ind)(n,i));
//#ifdef DROP
//						mvel(n) += tri_hp_cns::mesh_ref_vel(n);
//#endif
//					}
//
//					circulation += basis::tri(x.log2p)->wtx(i)*(-norm(1)*(x.u(0)(0,i)-mvel(0)) +norm(0)*(x.u(1)(0,i)-mvel(1)));
//
//					convect = basis::tri(x.log2p)->wtx(i)*RAD(x.crd(0)(0,i))*((x.u(0)(0,i)-mvel(0))*norm(0) +(x.u(1)(0,i)-mvel(1))*norm(1));
//					conv_flux(2) -= convect;
//					conv_flux(0) -= x.u(0)(0,i)*convect;
//					conv_flux(1) -= x.u(1)(0,i)*convect;
//
//#ifdef L2_ERROR
//					xpt(0) = x.crd(0)(0,i);
//					xpt(1) = x.crd(1)(0,i);
//					l2error += jcb*l2norm.Eval(xpt,x.gbl->time);
//#endif
//				}				
//			}
//			*x.gbl->log << base.idprefix << " circumference: " << circumference << std::endl;
//			*x.gbl->log << base.idprefix << " viscous/pressure flux: " << diff_flux << std::endl;
//			*x.gbl->log << base.idprefix << " convective flux: " << conv_flux << std::endl;
//			*x.gbl->log << base.idprefix << " circulation: " << circulation << std::endl;
//
//			/* OUTPUT AUXILIARY FLUXES */
//			*x.gbl->log << base.idprefix << "total fluxes: " << total_flux << std::endl;
//#ifdef L2_ERROR
//			*x.gbl->log << base.idprefix << "l2error: " << sqrt(l2error) << std::endl;
//#endif
//			break;
//		}
//
//		case(tri_hp::auxiliary): {
//			if (!report_flag) return;                
//
//			/* AUXILIARY FLUX METHOD */
//			int v0;
//			total_flux = 0.0;
//			for(ind=0; ind < base.nseg; ++ind) {
//				sind = base.seg(ind);
//				v0 = x.seg(sind).pnt(0);
//				total_flux += x.gbl->res.v(v0,Range::all());
//			}
//			v0 = x.seg(sind).pnt(1);
//			total_flux += x.gbl->res.v(v0,Range::all());
//		}
//	}

	return;
}


void neumann::element_rsdl(int eind, int stage) {
	int k,n,sind;
	TinyVector<FLT,2> pt,mvel,nrm;
	Array<FLT,1> u(x.NV),flx(x.NV);
	
	x.lf = 0.0;

	sind = base.seg(eind);
	
	x.crdtocht1d(sind);
	for(n=0;n<tri_mesh::ND;++n)
		basis::tri(x.log2p)->proj1d(&x.cht(n,0),&x.crd(n)(0,0),&x.dcrd(n,0)(0,0));

	for(n=0;n<x.NV;++n)
		basis::tri(x.log2p)->proj1d(&x.uht(n)(0),&x.u(n)(0,0));
	
	for(k=0;k<basis::tri(x.log2p)->gpx();++k) {
		nrm(0) = x.dcrd(1,0)(0,k);
		nrm(1) = -x.dcrd(0,0)(0,k);                
		for(n=0;n<tri_mesh::ND;++n) {
			pt(n) = x.crd(n)(0,k);
			mvel(n) = x.gbl->bd(0)*(x.crd(n)(0,k) -dxdt(x.log2p,eind)(n,k));
#ifdef DROP
			mvel(n) += tri_hp_cns_explicit::mesh_ref_vel(n);
#endif
		}
		
	
		for(n=0;n<x.NV;++n)
			u(n) = x.u(n)(0,k);
		
		flux(u,pt,mvel,nrm,flx);
		
		for(n=0;n<x.NV;++n)
			x.res(n)(0,k) = RAD(x.crd(0)(0,k))*flx(n);
		
	}
	
	for(n=0;n<x.NV;++n)
		basis::tri(x.log2p)->intgrt1d(&x.lf(n)(0),&x.res(n)(0,0));
	
	return;
}

void applied_stress::init(input_map& inmap,void* gbl_in) {
	std::string keyword;
	std::ostringstream nstr;

	neumann::init(inmap,gbl_in);

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
