/*
 *  cd_bdry.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 1/13/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

/*
 *  cd_bdry.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "bdry_cd.h"
#include "myblas.h"

using namespace bdry_cd;

void generic::output(std::ostream& fout, tri_hp::filetype typ,int tlvl) {
	int i,m,n,ind,sind,tind,seg;
	FLT visc[tri_mesh::ND];
	TinyVector<FLT,tri_mesh::ND> norm, mvel;
	FLT l,conv,diff,jcb;
	
	switch(typ) {
		case(tri_hp::text): case(tri_hp::binary): {
			hp_edge_bdry::output(fout,typ,tlvl);
			break;
		}
		case(tri_hp::tecplot): {
			if (!report_flag) return;
			
			std::ostringstream fname;
			fname << base.idprefix << '_' << x.gbl->tstep << ".dat";
			
			std::ofstream fout;
			fout.open(fname.str().c_str());
			
			fout << "VARIABLES=\"S\",\"X\",\"Y\",\"CFLUX\",\"DFLUX\"\nTITLE = " << base.idprefix << '\n'<< "ZONE\n";
			
			conv_total = 0.0;
			diff_total = 0.0;
			circumference = 0.0;
			ind = 0; 
			do { 
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
					
					x.cjcb(0,i) = x.gbl->nu/(x.dcrd(0,0)(0,i)*x.dcrd(1,1)(0,i) -x.dcrd(1,0)(0,i)*x.dcrd(0,1)(0,i));
										
					/* DIFFUSIVE FLUXES ( FOR EXTRA VARIABLES) */
					visc[0] =  x.cjcb(0,i)*(x.dcrd(1,1)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,1)(0,i)*x.dcrd(0,0)(0,i));
					visc[1] = -x.cjcb(0,i)*(x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i));
					l = sqrt(x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i)  +x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i));
					norm(0) = x.dcrd(1,0)(0,i)/l;
					norm(1) = -x.dcrd(0,0)(0,i)/l; 
					for(n=0;n<tri_mesh::ND;++n) {
						mvel(n) = x.gbl->bd(0)*(x.crd(n)(0,i) -dxdt(x.log2p,ind)(n,i));
					}					
#ifdef CONST_A
					conv = x.u(0)(0,i)*((x.gbl->ax -mvel(0))*norm(0) +(x.gbl->ay -mvel(1))*norm(1));
#else
					TinyVector<FLT,tri_mesh::ND> pt;
					pt(0) = x.crd(0)(0,i);
					pt(1) = x.crd(1)(0,i);
					conv = x.u(0)(0,i)*((x.gbl->a->f(0,pt,x.gbl->time) -mvel(0))*norm(0) +(x.gbl->a->f(1,pt,x.gbl->time) -mvel(1))*norm(1));
#endif
					conv_total -= basis::tri(x.log2p)->wtx(i)*RAD(x.crd(0)(0,i))*conv*l;
					
					diff = -visc[0]*x.du(0,0)(0,i) -visc[1]*x.du(0,1)(0,i)/l;
					
					diff_total -= basis::tri(x.log2p)->wtx(i)*RAD(x.crd(0)(0,i))*diff*l;
					
					fout << circumference << ' ' << x.crd(0)(0,i) << ' ' << x.crd(1)(0,i) << ' ' << conv << ' ' << diff << std::endl;

				}	
			} while (++ind < base.nseg);
			streamsize oldprecision = (*x.gbl->log).precision(10);
			*x.gbl->log << base.idprefix << " circumference: " << circumference << std::endl;
			*x.gbl->log << base.idprefix << " total diffusive flux: " << diff_total << std::endl;
			*x.gbl->log << base.idprefix << " total convective flux: " << conv_total << std::endl;
			(*x.gbl->log).precision(oldprecision);
			
			fout.close();
			break;
		}
		default: 
			break;
	}
	
	return;
}


void dirichlet::tadvance() {
	int j,k,m,n,v0,v1,sind=-2,indx,info;
	TinyVector<FLT,tri_mesh::ND> pt;
	char uplo[] = "U";

	hp_edge_bdry::tadvance(); 


	/* UPDATE BOUNDARY CONDITION VALUES */
	for(j=0;j<base.nseg;++j) {
		sind = base.seg(j);
		v0 = x.seg(sind).pnt(0);
		x.ug.v(v0,0) = ibc->f(0,x.pnts(v0),x.gbl->time);
	}
	v0 = x.seg(sind).pnt(1);
	x.ug.v(v0,0) = ibc->f(0,x.pnts(v0),x.gbl->time);

	/*******************/    
	/* SET SIDE VALUES */
	/*******************/
	for(j=0;j<base.nseg;++j) {
		sind = base.seg(j);
		v0 = x.seg(sind).pnt(0);
		v1 = x.seg(sind).pnt(1);

		if (is_curved()) {
			x.crdtocht1d(sind);
			for(n=0;n<tri_mesh::ND;++n)
				basis::tri(x.log2p)->proj1d(&x.cht(n,0),&x.crd(n)(0,0),&x.dcrd(n,0)(0,0));
		}
		else {
			for(n=0;n<tri_mesh::ND;++n) {
				basis::tri(x.log2p)->proj1d(x.pnts(v0)(n),x.pnts(v1)(n),&x.crd(n)(0,0));

				for(k=0;k<basis::tri(x.log2p)->gpx();++k)
					x.dcrd(n,0)(0,k) = 0.5*(x.pnts(v1)(n)-x.pnts(v0)(n));
			}
		}

		if (basis::tri(x.log2p)->sm()) {
			for(n=0;n<x.NV;++n)
				basis::tri(x.log2p)->proj1d(x.ug.v(v0,n),x.ug.v(v1,n),&x.res(n)(0,0));

			for(k=0;k<basis::tri(x.log2p)->gpx(); ++k) {
				pt(0) = x.crd(0)(0,k);
				pt(1) = x.crd(1)(0,k);
				for(n=0;n<x.NV;++n)
					x.res(n)(0,k) -= ibc->f(n,pt,x.gbl->time);
			}
			for(n=0;n<x.NV;++n)
				basis::tri(x.log2p)->intgrt1d(&x.lf(n)(0),&x.res(n)(0,0));

			indx = sind*x.sm0;
			for(n=0;n<x.NV;++n) {
				PBTRS(uplo,basis::tri(x.log2p)->sm(),basis::tri(x.log2p)->sbwth(),1,(double *) &basis::tri(x.log2p)->sdiag1d(0,0),basis::tri(x.log2p)->sbwth()+1,&x.lf(n)(2),basis::tri(x.log2p)->sm(),info);
				for(m=0;m<basis::tri(x.log2p)->sm();++m) 
					x.ug.s(sind,m,n) = -x.lf(n)(2+m);
			}
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
	
	// x.ugtouht1d(sind);
	
	for(n=0;n<x.NV;++n)
		basis::tri(x.log2p)->proj1d(&x.uht(n)(0),&x.u(n)(0,0));
	
	for(k=0;k<basis::tri(x.log2p)->gpx();++k) {
		nrm(0) = x.dcrd(1,0)(0,k);
		nrm(1) = -x.dcrd(0,0)(0,k);                
		for(n=0;n<tri_mesh::ND;++n) {
			pt(n) = x.crd(n)(0,k);
			mvel(n) = x.gbl->bd(0)*(x.crd(n)(0,k) -dxdt(x.log2p,eind)(n,k));
		}
		x.res(0)(0,k) = RAD(x.crd(0)(0,k))*flux(x.u(0)(0,k),pt,mvel,nrm);
		
	}
	for(n=0;n<x.NV;++n)
		basis::tri(x.log2p)->intgrt1d(&x.lf(n)(0),&x.res(n)(0,0));
	
	//	for(n=0;n<x.NV;++n)
	//		x.gbl->res.v(v0,n) += x.lf(n)(0);
	//
	//	for(n=0;n<x.NV;++n)
	//		x.gbl->res.v(v1,n) += x.lf(n)(1);
	//
	//	for(k=0;k<basis::tri(x.log2p)->sm();++k) {
	//		for(n=0;n<x.NV;++n)
	//			x.gbl->res.s(sind,k,n) += x.lf(n)(k+2);
	//	}
	
	return;
}


void melt_end_pt::rsdl(int stage) {
	// THIS IS SUPPOSED TO BE CALLED FROM tri_hp::update, but these are usually commented out.
	int ebdry = base.ebdry(1);
	int sind = x.ebdry(ebdry)->seg(0);
	x.crdtocht1d(sind);
	x.ugtouht1d(sind);
	element_rsdl();
	
#ifdef petsc
	/* Store residual in r_mesh residual vector */
	r_tri_mesh::global *r_gbl = dynamic_cast<r_tri_mesh::global *>(x.gbl);
	int v0 = x.seg(sind).pnt(0);
	/* Rotate residual for better diagonal dominance */
	r_gbl->res(v0)(0) = res;
#endif	

	return;
}

void melt_end_pt::element_rsdl() {
	TinyVector<FLT,2> xp,dxpdpsi;
	FLT T;
	FLT dTdpsi;
	basis::tri(x.log2p)->ptprobe1d(2,xp.data(),dxpdpsi.data(),-1.0,&x.cht(0,0),MXTM);
	basis::tri(x.log2p)->ptprobe1d(1,&T,&dTdpsi,-1.0,&x.uht(0)(0),MXTM);
	res = dTdpsi/dxpdpsi(0);
}


void melt_end_pt::setup_preconditioner() {
	int ebdry = base.ebdry(1);
	int sind = x.ebdry(ebdry)->seg(0);
	dx = x.distance(x.seg(sind).pnt(0),x.seg(sind).pnt(1));
}

void melt_end_pt::update(int stage) {
	if (x.coarse_level)
		return;
	
	if (stage == -1) {
		pnt0 = x.pnts(base.pnt)(0);
		return;
	}
	
	FLT cflalpha = x.gbl->alpha(stage)*x.gbl->cfl(x.log2p);
	x.pnts(base.pnt) = pnt0 -cflalpha*(res/dt > 0.0 ? 1. : -1.)*min(dx,fabs(res/dt));
}

#ifdef petsc
void melt_end_pt::petsc_jacobian() {
#ifdef BZ_DEBUG
	const FLT eps_r = 0.0e-6, eps_a = 1.0e-6;  /*<< constants for debugging jacobians */
#else
	const FLT eps_r = 1.0e-6, eps_a = 1.0e-10;  /*<< constants for accurate numerical determination of jacobians */
#endif
	const int sm = basis::tri(x.log2p)->sm();
	int nvars = sm +2 +2;  // Temperature values along side plus two x-endpoint positions */
	Array<int,1> cols(nvars);
	Array<FLT,1> vals(nvars);
	
	/* Jacobian for row determine x position */
	/* Variables effecting dT/dx are strictly those on the edge */
	rsdl(0);
	FLT res0 = res;
	
	FLT dw = 0.0;
	for(int i=0;i<2;++i)
		dw = dw + fabs(x.uht(0)(i));
		
	int ebdry = base.ebdry(1);
	int sind = x.ebdry(ebdry)->seg(0);
	
	dw = dw*eps_r;
	dw += eps_a;
	FLT dx = eps_r*x.distance(x.seg(sind).pnt(0),x.seg(sind).pnt(1)) +eps_a;
	
	/* Numerically create Jacobian */
	int kcol = 0;
	for(int mode = 0; mode < sm+2; ++mode) {
		x.uht(0)(mode) += dw;
		element_rsdl();
		vals(kcol++) = (res-res0)/dw;
		x.uht(0)(mode) -= dw;
	}
	
	for(int mode = 0; mode < 2; ++mode) {
		x.cht(0,mode) += dx;
		element_rsdl();
		vals(kcol++) = (res-res0)/dx;
		x.cht(0,mode) -= dx;
	}
	
	int vdofs;
	if (x.mmovement != x.coupled_deformable)
		vdofs = x.NV;
	else
		vdofs = x.NV+x.ND;
	
	int cind = 0;
	for(int k=0;k<2;++k) {
		cols(cind++) = vdofs*x.seg(sind).pnt(k);
	}
	
	/* EDGE MODES */
	if (sm > 0) {
		int gbl_eind = x.npnt*vdofs + sind*sm*x.NV;
		for (int m = 0; m < sm; ++m) {
			cols(cind++) = gbl_eind++;
		}
	}			
	
	/* X-vertex positions */
	for(int k=0;k<2;++k) {
		cols(cind++) = vdofs*x.seg(sind).pnt(k)+1;  // For x position
	}
	

	
#ifdef MY_SPARSE
	Array<int,1> rows(1);
	
	/* Add diagonal term ? */	
	dx = x.distance(x.seg(sind).pnt(0),x.seg(sind).pnt(1));
	FLT diag = vals(cind-2);
	FLT limit = MAX(fabs(res0/diag),dx)/dx -1.0;
	*x.gbl->log << "diag " << diag << " res0 " << res0 << " limit " << limit << " dx "  << dx << " diag_addition " << diag_addition << std::endl;
	vals(cind-2) += diag_addition*limit;

	rows(0) = base.pnt*vdofs+1;
	x.J.zero_rows(1,rows);
	x.J.add_values(rows(0),nvars,cols,vals);
#else
	*x.gbl->log << "petsc not working for melt_end_pt\n";
	MatSetValuesLocal(x.petsc_J,x.NV*(sm+2),rows.data(),vdofs*2 +x.NV*sm,cols.data(),K.data(),ADD_VALUES);
#endif
}
#endif
