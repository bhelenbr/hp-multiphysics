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
			*x.gbl->log << base.idprefix << " circumference: " << circumference << std::endl;
			*x.gbl->log << base.idprefix << " total diffusive flux: " << precision(10) << diff_total << std::endl;
			*x.gbl->log << base.idprefix << " total convective flux: " << conv_total << std::endl;
			
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

void neumann::rsdl(int stage) {
	int j,k,n,v0,v1,sind;
	TinyVector<FLT,2> pt,mvel,nrm;

	for(j=0;j<base.nseg;++j) {
		sind = base.seg(j);
		v0 = x.seg(sind).pnt(0);
		v1 = x.seg(sind).pnt(1);

		x.crdtocht1d(sind);
		for(n=0;n<tri_mesh::ND;++n)
			basis::tri(x.log2p)->proj1d(&x.cht(n,0),&x.crd(n)(0,0),&x.dcrd(n,0)(0,0));

		x.crdtocht1d(sind,1);
		for(n=0;n<tri_mesh::ND;++n)
			basis::tri(x.log2p)->proj1d(&x.cht(n,0),&x.crd(n)(1,0));

		x.ugtouht1d(sind);
		for(n=0;n<x.NV;++n)
			basis::tri(x.log2p)->proj1d(&x.uht(n)(0),&x.u(n)(0,0));

		for(k=0;k<basis::tri(x.log2p)->gpx();++k) {
			nrm(0) = x.dcrd(1,0)(0,k);
			nrm(1) = -x.dcrd(0,0)(0,k);                
			for(n=0;n<tri_mesh::ND;++n) {
				pt(n) = x.crd(n)(0,k);
				mvel(n) = x.gbl->bd(0)*(x.crd(n)(0,k) -x.crd(n)(1,k));
			}

			x.res(0)(0,k) = RAD(x.crd(0)(0,k))*flux(x.u(0)(0,k),pt,mvel,nrm);
		}

		for(n=0;n<x.NV;++n)
			basis::tri(x.log2p)->intgrt1d(&x.lf(n)(0),&x.res(n)(0,0));

		for(n=0;n<x.NV;++n)
			x.gbl->res.v(v0,n) += x.lf(n)(0);

		for(n=0;n<x.NV;++n)
			x.gbl->res.v(v1,n) += x.lf(n)(1);

		for(k=0;k<basis::tri(x.log2p)->sm();++k) {
			for(n=0;n<x.NV;++n)
				x.gbl->res.s(sind,k,n) += x.lf(n)(k+2);
		}
	}

	return;
}
