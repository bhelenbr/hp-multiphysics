#include "bdry_cns.h"
#include <myblas.h>

/*************************************************/
/* SET DIRICHLET BOUNDARY VALUES & FLUXES ********/
/* (THINGS THAT ARE INDEPENDENT OF THE SOLUTION) */
/* BUT NOT INDEPENDENT OF TIME/MESH POSITION     */
/*************************************************/
using namespace bdry_cns;

void generic::output(std::ostream& fout, tri_hp::filetype typ,int tlvl) {
	int i,m,n,ind,sind,tind,seg;
	FLT visc[tri_mesh::ND+1][tri_mesh::ND+1][tri_mesh::ND][tri_mesh::ND];
	TinyVector<FLT,tri_mesh::ND> norm, mvel;
	FLT convect,jcb;

	switch(typ) {
		case(tri_hp::text): case(tri_hp::binary): {
			hp_edge_bdry::output(fout,typ,tlvl);
			break;
		}
		case(tri_hp::tecplot): {
			if (!report_flag) return;

			conv_flux = 0.0;
			diff_flux = 0.0;
			moment = 0.0;
			circumference = 0.0;
			circulation = 0.0;
#ifdef L2_ERROR
			FLT l2error = 0.0;
			TinyVector<FLT,2> xpt;
#endif

			for(ind=0; ind < base.nseg; ++ind) {
				sind = base.seg(ind);
				tind = x.seg(sind).tri(0);        

				for(seg=0;seg<3;++seg)
					if (x.tri(tind).seg(seg) == sind) break;
				assert(seg != 3);

				x.crdtocht(tind);
				for(m=basis::tri(x.log2p)->bm();m<basis::tri(x.log2p)->tm();++m)
					for(n=0;n<tri_mesh::ND;++n)
						x.cht(n,m) = 0.0;

				for(n=0;n<tri_mesh::ND;++n)
					basis::tri(x.log2p)->proj_side(seg,&x.cht(n,0), &x.crd(n)(0,0), &x.dcrd(n,0)(0,0), &x.dcrd(n,1)(0,0));

				x.ugtouht(tind);
				for(n=0;n<x.NV;++n)
					basis::tri(x.log2p)->proj_side(seg,&x.uht(n)(0),&x.u(n)(0,0),&x.du(n,0)(0,0),&x.du(n,1)(0,0));

				for (i=0;i<basis::tri(x.log2p)->gpx();++i) {
					jcb =  basis::tri(x.log2p)->wtx(i)*RAD(x.crd(0)(0,i))*sqrt(x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i) +x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i));
					circumference += jcb;

					x.cjcb(0,i) = x.gbl->mu*RAD(x.crd(0)(0,i))/(x.dcrd(0,0)(0,i)*x.dcrd(1,1)(0,i) -x.dcrd(1,0)(0,i)*x.dcrd(0,1)(0,i));

					/* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
					/* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
					visc[0][0][0][0] =  x.cjcb(0,i)*(2.*x.dcrd(1,1)(0,i)*x.dcrd(1,1)(0,i) +x.dcrd(0,1)(0,i)*x.dcrd(0,1)(0,i));
					visc[0][0][1][1] =  x.cjcb(0,i)*(2.*x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i));
					visc[0][0][0][1] = -x.cjcb(0,i)*(2.*x.dcrd(1,1)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,1)(0,i)*x.dcrd(0,0)(0,i));
#define             viscI0II0II1II0I visc[0][0][0][1]

					visc[1][1][0][0] =  x.cjcb(0,i)*(x.dcrd(1,1)(0,i)*x.dcrd(1,1)(0,i) +2.*x.dcrd(0,1)(0,i)*x.dcrd(0,1)(0,i));
					visc[1][1][1][1] =  x.cjcb(0,i)*(x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i) +2.*x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i));
					visc[1][1][0][1] = -x.cjcb(0,i)*(x.dcrd(1,1)(0,i)*x.dcrd(1,0)(0,i) +2.*x.dcrd(0,1)(0,i)*x.dcrd(0,0)(0,i));
#define             viscI1II1II1II0I visc[1][1][0][1]

					visc[0][1][0][0] = -x.cjcb(0,i)*x.dcrd(0,1)(0,i)*x.dcrd(1,1)(0,i);
					visc[0][1][1][1] = -x.cjcb(0,i)*x.dcrd(0,0)(0,i)*x.dcrd(1,0)(0,i);
					visc[0][1][0][1] =  x.cjcb(0,i)*x.dcrd(0,1)(0,i)*x.dcrd(1,0)(0,i);
					visc[0][1][1][0] =  x.cjcb(0,i)*x.dcrd(0,0)(0,i)*x.dcrd(1,1)(0,i);

					/* OTHER SYMMETRIES     */                
#define             viscI1II0II0II0I visc[0][1][0][0]
#define             viscI1II0II1II1I visc[0][1][1][1]
#define             viscI1II0II0II1I visc[0][1][1][0]
#define             viscI1II0II1II0I visc[0][1][0][1]

					/* DIFFUSIVE FLUXES ( FOR EXTRA VARIABLES) */
					visc[2][2][1][0] =  x.cjcb(0,i)*(x.dcrd(1,1)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,1)(0,i)*x.dcrd(0,0)(0,i));
					visc[2][2][1][1] = -x.cjcb(0,i)*(x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i));

					for (n=tri_mesh::ND;n<x.NV-1;++n) 
						diff_flux(n) -= x.gbl->D(n)/x.gbl->mu*basis::tri(x.log2p)->wtx(i)*(-visc[2][2][1][0]*x.du(2,0)(0,i) -visc[2][2][1][1]*x.du(2,1)(0,i));

					diff_flux(0) -=    basis::tri(x.log2p)->wtx(i)*(-x.u(2)(0,i)*RAD(x.crd(0)(0,i))*x.dcrd(1,0)(0,i) 
									-viscI0II0II1II0I*x.du(0,0)(0,i) -visc[0][1][1][0]*x.du(1,0)(0,i)
									-visc[0][0][1][1]*x.du(0,1)(0,i) -visc[0][1][1][1]*x.du(1,1)(0,i));															
					diff_flux(1) -=    basis::tri(x.log2p)->wtx(i)*( x.u(2)(0,i)*RAD(x.crd(0)(0,i))*x.dcrd(0,0)(0,i)
									-viscI1II0II1II0I*x.du(0,0)(0,i) -viscI1II1II1II0I*x.du(1,0)(0,i)
									-viscI1II0II1II1I*x.du(0,1)(0,i) -visc[1][1][1][1]*x.du(1,1)(0,i));


					norm(0) = x.dcrd(1,0)(0,i);
					norm(1) = -x.dcrd(0,0)(0,i);                
					for(n=0;n<tri_mesh::ND;++n) {
						mvel(n) = x.gbl->bd(0)*(x.crd(n)(0,i) -dxdt(x.log2p,ind)(n,i));
#ifdef DROP
						mvel(n) += tri_hp_cns::mesh_ref_vel(n);
#endif
					}

					circulation += basis::tri(x.log2p)->wtx(i)*(-norm(1)*(x.u(0)(0,i)-mvel(0)) +norm(0)*(x.u(1)(0,i)-mvel(1)));

					convect = basis::tri(x.log2p)->wtx(i)*RAD(x.crd(0)(0,i))*((x.u(0)(0,i)-mvel(0))*norm(0) +(x.u(1)(0,i)-mvel(1))*norm(1));
					conv_flux(2) -= convect;
					conv_flux(0) -= x.u(0)(0,i)*convect;
					conv_flux(1) -= x.u(1)(0,i)*convect;

#ifdef L2_ERROR
					xpt(0) = x.crd(0)(0,i);
					xpt(1) = x.crd(1)(0,i);
					l2error += jcb*l2norm.Eval(xpt,x.gbl->time);
#endif
				}				
			}
			*x.gbl->log << base.idprefix << " circumference: " << circumference << std::endl;
			*x.gbl->log << base.idprefix << " viscous/pressure flux: " << diff_flux << std::endl;
			*x.gbl->log << base.idprefix << " convective flux: " << conv_flux << std::endl;
			*x.gbl->log << base.idprefix << " circulation: " << circulation << std::endl;

			/* OUTPUT AUXILIARY FLUXES */
			*x.gbl->log << base.idprefix << "total fluxes: " << total_flux << std::endl;
#ifdef L2_ERROR
			*x.gbl->log << base.idprefix << "l2error: " << sqrt(l2error) << std::endl;
#endif
			break;
		}

		case(tri_hp::auxiliary): {
			if (!report_flag) return;                

			/* AUXILIARY FLUX METHOD */
			int v0;
			total_flux = 0.0;
			for(ind=0; ind < base.nseg; ++ind) {
				sind = base.seg(ind);
				v0 = x.seg(sind).pnt(0);
				total_flux += x.gbl->res.v(v0,Range::all());
			}
			v0 = x.seg(sind).pnt(1);
			total_flux += x.gbl->res.v(v0,Range::all());
		}
	}

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
			mvel(n) += tri_hp_cns::mesh_ref_vel(n);
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

//void inflow::setvalues(init_bdry_cndtn *ibc) {
//	int j,k,m,n,v0,v1,sind,indx,info;
//	TinyVector<FLT,tri_mesh::ND> pt;
//	char uplo[] = "U";
//
//	/* UPDATE BOUNDARY CONDITION VALUES */
//	for(j=0;j<base.nseg;++j) {
//		sind = base.seg(j);
//		v0 = x.seg(sind).pnt(0);
//		for(n=1;n<x.NV;++n)
//			x.ug.v(v0,n) = ibc->f(n,x.pnts(v0),x.gbl->time);
//    }
//	v0 = x.seg(sind).pnt(1);
//	for(n=1;n<x.NV;++n)
//		x.ug.v(v0,n) = ibc->f(n,x.pnts(v0),x.gbl->time);
//
//	/*******************/    
//	/* SET SIDE VALUES */
//	/*******************/
//	for(j=0;j<base.nseg;++j) {
//		sind = base.seg(j);
//		v0 = x.seg(sind).pnt(0);
//		v1 = x.seg(sind).pnt(1);
//
//		if (is_curved()) {
//			x.crdtocht1d(sind);
//			for(n=0;n<tri_mesh::ND;++n)
//				basis::tri(x.log2p)->proj1d(&x.cht(n,0),&x.crd(n)(0,0),&x.dcrd(n,0)(0,0));
//		}
//		else {
//			for(n=0;n<tri_mesh::ND;++n) {
//				basis::tri(x.log2p)->proj1d(x.pnts(v0)(n),x.pnts(v1)(n),&x.crd(n)(0,0));
//
//				for(k=0;k<basis::tri(x.log2p)->gpx();++k)
//					x.dcrd(n,0)(0,k) = 0.5*(x.pnts(v1)(n)-x.pnts(v0)(n));
//			}
//		}
//
//		if (basis::tri(x.log2p)->sm()) {
//			for(n=1;n<x.NV;++n)
//				basis::tri(x.log2p)->proj1d(x.ug.v(v0,n),x.ug.v(v1,n),&x.res(n)(0,0));
//
//			for(k=0;k<basis::tri(x.log2p)->gpx(); ++k) {
//				pt(0) = x.crd(0)(0,k);
//				pt(1) = x.crd(1)(0,k);
//				for(n=1;n<x.NV;++n)
//					x.res(n)(0,k) -= ibc->f(n,pt,x.gbl->time);
//			}
//			for(n=1;n<x.NV;++n)
//				basis::tri(x.log2p)->intgrt1d(&x.lf(n)(0),&x.res(n)(0,0));
//
//			indx = sind*x.sm0;
//			for(n=1;n<x.NV;++n) {
//				PBTRS(uplo,basis::tri(x.log2p)->sm(),basis::tri(x.log2p)->sbwth(),1,(double *) &basis::tri(x.log2p)->sdiag1d(0,0),basis::tri(x.log2p)->sbwth()+1,&x.lf(n)(2),basis::tri(x.log2p)->sm(),info);
//				for(m=0;m<basis::tri(x.log2p)->sm();++m) 
//					x.ug.s(sind,m,n) = -x.lf(n)(2+m);
//			}
//		}
//	}
//	return;
//}

//void symmetry::tadvance() {
//	int j,m,v0,sind;
//	TinyVector<FLT,tri_mesh::ND> pt;
//
//	hp_edge_bdry::tadvance();
//
//	/* UPDATE BOUNDARY CONDITION VALUES */
//	for(j=0;j<base.nseg;++j) {
//		sind = base.seg(j);
//		v0 = x.seg(sind).pnt(0);
//		x.ug.v(v0,dir) = 0.0;
//	}
//	v0 = x.seg(sind).pnt(1);
//	x.ug.v(v0,dir) = 0.0;
//
//	/*******************/    
//	/* SET SIDE VALUES */
//	/*******************/
//	for(j=0;j<base.nseg;++j) {
//		sind = base.seg(j);
//		for(m=0;m<basis::tri(x.log2p)->sm();++m) {
//			x.ug.s(sind,m,dir) = 0.0;
//		}
//	}
//
//	return;
//}


void characteristic::flux(Array<FLT,1>& pvu, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, Array<FLT,1>& flx) {	
	
	TinyVector<FLT,4> lambda,Rl,Rr,ub,Roe,fluxtemp,cvl,cvr;
	Array<FLT,2> A(x.NV,x.NV),V(x.NV,x.NV),VINV(x.NV,x.NV),temp(x.NV,x.NV);
	Array<FLT,1> Aeigs(x.NV);
	FLT gam = x.gbl->gamma;
	FLT gm1 = gam-1.0;
	
	FLT mag = sqrt(norm(0)*norm(0) + norm(1)*norm(1));
	norm /= mag;
	
	/* Left */
	/* Rotate Coordinate System */
	FLT ul =  pvu(1)*norm(0) +pvu(2)*norm(1);
	FLT vl = -pvu(1)*norm(1) +pvu(2)*norm(0);
	pvu(1) = ul, pvu(2) = vl;
	
	/* Roe Variables */
	Rl(0) = sqrt(pvu(0)/pvu(x.NV-1)); // sqrt(rho)
	Rl(1) = Rl(0)*ul; // sqrt(rho)*u
	Rl(2) = Rl(0)*vl; // sqrt(rho)*v
	Rl(3) = Rl(0)*(pvu(x.NV-1)/gm1+0.5*(ul*ul+vl*vl)); // sqrt(rho)*E
	
	/* Right */
	for(int n=0;n<x.NV;++n)
		ub(n) = ibc->f(n,xpt,x.gbl->time);
	
	/* Rotate Coordinate System */
	FLT ur =  ub(1)*norm(0) +ub(2)*norm(1);
	FLT vr = -ub(1)*norm(1) +ub(2)*norm(0);
	ub(1) = ur, ub(2) = vr;
	
	/* Roe Variables */
	Rr(0) = sqrt(ub(0)/ub(x.NV-1));
	Rr(1) = Rr(0)*ur;
	Rr(2) = Rr(0)*vr;
	Rr(3) = Rr(0)*(ub(x.NV-1)/gm1+0.5*(ur*ur+vr*vr));	
	
	/* Average Roe Variables */
	Roe = 0.5*(Rl+Rr);
	
	/* Calculate u,v,c Variables */
	FLT u = Roe(1)/Roe(0);
	FLT v = Roe(2)/Roe(0);
	FLT ke = 0.5*(u*u+v*v);
	FLT E = Roe(3)/Roe(0);
	FLT rt = gm1*(E-ke);
	FLT c2 = gam*rt;
	FLT c = sqrt(c2);
	
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
	
	/* convert from primitive to conservative variables */
	cvl(0) = pvu(0)/pvu(x.NV-1); // rho
	cvl(1) = cvl(0)*pvu(1); // rho*u
	cvl(2) = cvl(0)*pvu(2); // rho*v
	cvl(3) = cvl(0)*(pvu(x.NV-1)/gm1+0.5*(pvu(1)*pvu(1)+pvu(2)*pvu(2))); // rho*E
	
	cvr(0) = ub(0)/ub(x.NV-1);
	cvr(1) = cvr(0)*ub(1);
	cvr(2) = cvr(0)*ub(2);
	cvr(3) = cvr(0)*(ub(x.NV-1)/gm1+0.5*(ub(1)*ub(1)+ub(2)*ub(2)));
	
	fluxtemp = 0.0;
	
	for(int i = 0; i < x.NV; ++i)
		for(int j = 0; j < x.NV; ++j)
			fluxtemp(i) += 0.5*A(i,j)*(cvr(j)+cvl(j));
		
	/* eigenvalues of df/dw */
	Aeigs = fabs(u),fabs(u),fabs(u+c),fabs(u-c);	
	
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
			fluxtemp(i) -= 0.5*A(i,j)*(cvr(j)-cvl(j));
	
	/* CHANGE BACK TO X,Y COORDINATES */
	flx(0) = fluxtemp(0);
	flx(1) = fluxtemp(1)*norm(0) - fluxtemp(2)*norm(1);
	flx(2) = fluxtemp(1)*norm(1) + fluxtemp(2)*norm(0);
	flx(3) = fluxtemp(3);
	
	flx *= mag;
	
	return;
}

// old way
//void characteristic::flux(Array<FLT,1>& pvu, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, Array<FLT,1>& flx) {	
//	
//	TinyVector<FLT,4> lambda,Rl,Rr,ub,Roe,fluxtemp,fluxleft, fluxright;
//	Array<FLT,2> A(x.NV,x.NV),V(x.NV,x.NV),VINV(x.NV,x.NV),temp(x.NV,x.NV);
//	Array<FLT,1> Aeigs(x.NV);
//	FLT gam = x.gbl->gamma;
//	FLT gm1 = gam-1.0;
//	FLT gogm1 = gam/gm1;
//	
//	FLT mag = sqrt(norm(0)*norm(0) + norm(1)*norm(1));
//	norm /= mag;
//	
//	/* Left */
//	/* Rotate Coordinate System */
//	FLT ul =  pvu(1)*norm(0) +pvu(2)*norm(1);
//	FLT vl = -pvu(1)*norm(1) +pvu(2)*norm(0);
//	pvu(1) = ul, pvu(2) = vl;
//	
//	/* Roe Variables */
//	Rl(0) = sqrt(pvu(0)/pvu(x.NV-1)); // sqrt(rho)
//	Rl(1) = Rl(0)*ul; // sqrt(rho)*u
//	Rl(2) = Rl(0)*vl; // sqrt(rho)*v
//	Rl(3) = Rl(0)*(pvu(x.NV-1)/gm1+0.5*(ul*ul+vl*vl)); // sqrt(rho)*E
//	
//	/* Right */
//	for(int n=0;n<x.NV;++n)
//		ub(n) = ibc->f(n,xpt,x.gbl->time);
//	
//	/* Rotate Coordinate System */
//	FLT ur =  ub(1)*norm(0) +ub(2)*norm(1);
//	FLT vr = -ub(1)*norm(1) +ub(2)*norm(0);
//	ub(1) = ur, ub(2) = vr;
//	
//	/* Roe Variables */
//	Rr(0) = sqrt(ub(0)/ub(x.NV-1));
//	Rr(1) = Rr(0)*ur;
//	Rr(2) = Rr(0)*vr;
//	Rr(3) = Rr(0)*(ub(x.NV-1)/gm1+0.5*(ur*ur+vr*vr));	
//	
//	/* Average Roe Variables */
//	Roe = 0.5*(Rl+Rr);
//	
//	/* Calculate u,v,c Variables */
//	FLT rho = Roe(0)*Roe(0);
//	FLT u = Roe(1)/Roe(0);
//	FLT v = Roe(2)/Roe(0);
//	FLT ke = 0.5*(u*u+v*v);
//	FLT E = Roe(3)/Roe(0);
//	FLT rt = gm1*(E-ke);
//	FLT pr = rho*rt;
//	FLT c2 = gam*rt;
//	FLT c = sqrt(c2);
//	
////	/* eigenvectors of P*df/dw */
////	V = 0.0, 0.0, 1.0, 1.0,
////		0.0, 0.0, c/(gam*pr), -c/(gam*pr),
////		1.0, 0.0, 0.0, 0.0,
////		0.0, 1.0, gm1/(gam*rho), gm1/(gam*rho);
////	
////	/* eigenvalues of P*df/dw  */
////	Aeigs = u,u,u+c,u-c;
////	
////	/* inverse of eigenvectors (P*df/dw)^-1 */		
////	VINV = 0.0,            0.0,           1.0, 0.0,
////		   -gm1/(gam*rho), 0.0,           0.0, 1.0,
////		   0.5,            0.5*gam*pr/c,  0.0, 0.0,
////		   0.5,            -0.5*gam*pr/c, 0.0, 0.0;
////	
////	for(int i=0; i < x.NV; ++i)
////		for(int j=0; j < x.NV; ++j)
////			temp(i,j) = Aeigs(i)*VINV(i,j);
////	
////	A = 0.0;
////	for(int i=0; i<x.NV; ++i)
////		for(int j=0; j<x.NV; ++j)
////			for(int k=0; k<x.NV; ++k)
////				A(i,j)+=V(i,k)*temp(k,j);
//	
//	
//	/* df/dw */
//	A = u/rt,               rho,                       0.0,     -rho*u/rt,
//		u*u/rt+1.0,         2.0*rho*u,                 0.0,     -rho*u*u/rt,
//		u*v/rt,             rho*v,                     rho*u,   -rho*u*v/rt,
//		u*(gogm1*rt+ke)/rt, rho*(gogm1*rt+ke)+rho*u*u, rho*u*v, -rho*u*(gogm1*rt+ke)/rt+rho*u*gogm1;
//	
//	fluxtemp = 0.0;
//	
//	for(int i = 0; i < x.NV; ++i)
//		for(int j = 0; j < x.NV; ++j)
//			fluxtemp(i) += 0.5*A(i,j)*(ub(j)+pvu(j));
//	
//	/* or this way */
//	fluxtemp(0) = rho*u;
//	fluxtemp(1) = rho*u*u+pr;
//	fluxtemp(2) = rho*u*v;
//	fluxtemp(3) = rho*u*(gogm1*rt+ke);
//
////	fluxleft(0) = pvu(0)*pvu(1)/pvu(x.NV-1);
////	fluxleft(1) = fluxleft(0)*pvu(1)+pvu(0);
////	fluxleft(2) = fluxleft(0)*pvu(2);
////	fluxleft(3) = fluxleft(0)*(gogm1*pvu(x.NV-1)+0.5*(pvu(1)*pvu(1)+pvu(2)*pvu(2)));
////
////	fluxright(0) = ub(0)*ub(1)/ub(x.NV-1);
////	fluxright(1) = fluxright(0)*ub(1)+ub(0);
////	fluxright(2) = fluxright(0)*ub(2);
////	fluxright(3) = fluxright(0)*(gogm1*ub(x.NV-1)+0.5*(ub(1)*ub(1)+ub(2)*ub(2)));
////	
////	fluxtemp = 0.5*(fluxleft+fluxright);
//	
////	Aeigs = fabs(u),fabs(u),fabs(u+c),fabs(u-c);
////	
////	
////	for(int i=0; i < x.NV; ++i)
////		for(int j=0; j < x.NV; ++j)
////			temp(i,j) = Aeigs(i)*VINV(i,j);
////	
////	A = 0.0;
////	for(int i=0; i<x.NV; ++i)
////		for(int j=0; j<x.NV; ++j)
////			for(int k=0; k<x.NV; ++k)
////				A(i,j)+=V(i,k)*temp(k,j);
//
//	matrix_absolute_value(A);
//	
//	for(int i = 0; i < x.NV; ++i)
//		for(int j = 0; j < x.NV; ++j)
//			fluxtemp(i) -= 0.5*A(i,j)*(ub(j)-pvu(j));
//	
//	/* CHANGE BACK TO X,Y COORDINATES */
//	flx(0) = fluxtemp(0);
//	flx(1) = fluxtemp(1)*norm(0) - fluxtemp(2)*norm(1);
//	flx(2) = fluxtemp(1)*norm(1) + fluxtemp(2)*norm(0);
//	flx(3) = fluxtemp(3);
//
//	flx *= mag;
//		
//	return;
//}

//void hybrid_slave_pt::update(int stage) {
//
//	if (stage == -1) return;
//
//	base.sndsize() = 4;
//	base.sndtype() = boundary::flt_msg;
//	base.fsndbuf(0) = 0.0;
//	base.fsndbuf(1) = 0.0;
//	base.fsndbuf(2) = 0.0;
//	base.fsndbuf(3) = 0.0;
//
//	base.comm_prepare(boundary::all,0,boundary::symmetric);
//	base.comm_exchange(boundary::all,0,boundary::symmetric);
//	base.comm_wait(boundary::all,0,boundary::symmetric);
//
//	for(int m=0;m<base.nmatches();++m) {
//		for(int i=0;i<4;++i) 
//			base.fsndbuf(i) += base.frcvbuf(m,i);
//	}
//
//	if (base.fsndbuf(0)*base.fsndbuf(2) > 0.0) {
//		*x.gbl->log << "uh-oh opposite characteristics at hybrid point" << std::endl;
//		*x.gbl->log << "local " << base.idprefix << ' ' << base.fsndbuf(0) << "remote " << base.fsndbuf(2) << std::endl;
//	}
//
//	if (base.fsndbuf(0) > 0.0) {
//		x.pnts(base.pnt)(1) = base.fsndbuf(3);
//	}
//	else {
//		x.pnts(base.pnt)(1) = base.fsndbuf(1);
//	}
//}

//void hybrid_pt::rsdl(int stage) {
//	int sind,v0,v1;
//	TinyVector<FLT,2> tang,vel;
//	FLT tangvel;
//
//	if (surfbdry == 0) {
//		sind = x.ebdry(base.ebdry(0))->seg(x.ebdry(base.ebdry(0))->nseg-1);
//		v0 = x.seg(sind).pnt(1);
//		v1 = x.seg(sind).pnt(0);
//	}
//	else {
//		sind = x.ebdry(base.ebdry(1))->seg(0);
//		v0 = x.seg(sind).pnt(0);
//		v1 = x.seg(sind).pnt(1);
//	}
//
//
//	/* TANGENT POINTS INTO DOMAIN ALONG SURFACE */
//	tang(0) =  (x.pnts(v1)(0) -x.pnts(v0)(0));
//	tang(1) =  (x.pnts(v1)(1) -x.pnts(v0)(1));
//
//	vel(0) = 0.5*(x.ug.v(v0,0)-(x.gbl->bd(0)*(x.pnts(v0)(0) -x.vrtxbd(1)(v0)(0))) +
//				  x.ug.v(v1,0)-(x.gbl->bd(0)*(x.pnts(v1)(0) -x.vrtxbd(1)(v1)(0))));
//	vel(1) = 0.5*(x.ug.v(v0,1)-(x.gbl->bd(0)*(x.pnts(v0)(1) -x.vrtxbd(1)(v0)(1))) +
//				  x.ug.v(v1,1)-(x.gbl->bd(0)*(x.pnts(v1)(1) -x.vrtxbd(1)(v1)(1))));
//	tangvel = vel(0)*tang(0)+vel(1)*tang(1);
//
//	if (tangvel > 0.0)
//		fix_norm = 1;
//	else
//		fix_norm = 0;
//
//	surface_outflow_planar::rsdl(stage);
//}



//void hybrid_pt::update(int stage) {
//
//	if (stage == -1) return;
//
//	base.sndsize() = 4;
//	base.sndtype() = boundary::flt_msg;
//	base.fsndbuf(0) = 2*fix_norm-1.0;
//	base.fsndbuf(1) = x.pnts(base.pnt)(1);
//	base.fsndbuf(2) = 0.0;
//	base.fsndbuf(3) = 0.0;
//
//	base.comm_prepare(boundary::all,0,boundary::symmetric);
//	base.comm_exchange(boundary::all,0,boundary::symmetric);
//	base.comm_wait(boundary::all,0,boundary::symmetric);
//
//	for(int m=0;m<base.nmatches();++m) {
//		for(int i=0;i<4;++i) 
//			base.fsndbuf(i) += base.frcvbuf(m,i);
//	}
//
//	if (base.fsndbuf(0)*base.fsndbuf(2) > 0.0) {
//		*x.gbl->log << "uh-oh opposite characteristics at hybrid point" << std::endl;
//		*x.gbl->log << "local "  << base.idprefix << ' ' << base.fsndbuf(0) << "remote " << base.fsndbuf(2) << std::endl;
//	}
//
//	if (fix_norm) {
//		x.pnts(base.pnt)(1) = base.fsndbuf(3);
//	}
//}
