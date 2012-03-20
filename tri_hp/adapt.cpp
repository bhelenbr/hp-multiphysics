/*
 *  adapt.cpp
 *  planar++
 *
 *  Created by helenbrk on Tue Oct 23 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp.h"
#include "hp_boundary.h"
#include <myblas.h>

void tri_hp::adapt() {

		
	treeinit();  // FIXME??
	gbl->pstr->copy(*this);
	tri_mesh::adapt();
	setinfo();
	
#ifdef petsc
	petsc_finalize();
	petsc_initialize();
#endif

	if (gbl->adapt_output) {
		std::ostringstream fname;
		fname << "adapted_solution" << gbl->tstep << '_' << gbl->idprefix;
		tri_mesh::output(fname.str().c_str(),tri_mesh::grid);
		tri_hp::output(fname.str().c_str(),tri_hp::tecplot);
	}
	
	*gbl->log << "# Adaptation Complete with DOF: " << npnt +nseg*sm0 +ntri*im0 << std::endl;
	
	// FOR TESTING ADAPTATION TO A SPECIFIED FUNCTION
	// tobasis(gbl->ibc);
	
	return;
}


void tri_hp::length() {
	
	const FLT alpha = 2.0*(basis::tri(log2p)->p()-1.0+ND)/static_cast<FLT>(ND);

	sim::blks.allreduce(gbl->eanda.data(),gbl->eanda_recv.data(),3,blocks::flt_msg,blocks::sum);
	FLT energy2 = gbl->eanda_recv(0);
	FLT e2to_pow = gbl->eanda_recv(1);
	FLT totalerror2 = gbl->eanda_recv(2);

	if (gbl->error_estimator == global::energy_norm) {
		*gbl->log << "# DOF: " << npnt +nseg*sm0 +ntri*im0 << " Normalized Error " << sqrt(totalerror2/energy2) << " Target " << gbl->error_target << '\n';
		
		/* Determine error target (SEE AEA Paper) */
		FLT etarget2 = gbl->error_target*gbl->error_target*energy2;
		FLT K = pow(etarget2/e2to_pow,1./(ND*alpha));
		gbl->res.v(Range(0,npnt-1),0) = 1.0;
		gbl->res_r.v(Range(0,npnt-1),0) = 0.0;
		for(int tind=0;tind<ntri;++tind) {
			FLT error2 = gbl->fltwk(tind);
			FLT ri = K*pow(error2, -1./(ND*(1.+alpha)));
			for (int j=0;j<3;++j) {
				int p0 = tri(tind).pnt(j);
				/* Calculate average at vertices */
				gbl->res.v(p0,0) *= ri;
				gbl->res_r.v(p0,0) += 1.0;
			}
		}
	}
	else if (gbl->error_estimator == global::scale_independent) {
		/* This is to maintain a constant local truncation error (independent of scale) */
		gbl->res.v(Range(0,npnt-1),0) = 1.0;
		gbl->res_r.v(Range(0,npnt-1),0) = 0.0;
		for(int tind=0;tind<ntri;++tind) {
			FLT jcb = 0.25*area(tind);
			gbl->fltwk(tind) = sqrt(gbl->fltwk(tind)/jcb)/gbl->error_target;
			FLT error = gbl->fltwk(tind);  // Magnitude of local truncation error
			FLT ri = pow(error, -1./(basis::tri(log2p)->p()));
			for (int j=0;j<3;++j) {
				int p0 = tri(tind).pnt(j);
				/* Calculate average at vertices */
				gbl->res.v(p0,0) *= ri;
				gbl->res_r.v(p0,0) += 1.0;
			}
		}	
	}
	else {
		*gbl->log << "Unknown error estimator??" << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
	
	for (int pind=0;pind<npnt;++pind) {
		FLT ri = pow(gbl->res.v(pind,0),1.0/gbl->res_r.v(pind,0));
		if (ri < 2.0 && ri > 0.5)
			ri = 1.0;
		gbl->res.v(pind,0) = ri;
	}
	
	/* This is to smooth the change to the length function */
	int iter,sind,i,j,p0,p1;
	int niter = 1;
	
	for(i=0;i<npnt;++i)
	pnt(i).info = 0;
	
	for(i=0;i<nebd;++i) {
		for(j=0;j<ebdry(i)->nseg;++j) {
			sind = ebdry(i)->seg(j);
			pnt(seg(sind).pnt(0)).info = -1;
			pnt(seg(sind).pnt(1)).info = -1;
		}
	}
	
	for(iter=0; iter< niter; ++iter) {
		/* SMOOTH POINT DISTRIBUTION IN INTERIOR*/
		for(i=0;i<npnt;++i)
			gbl->res_r.v(i,0) = 0.0;
		
		for(i=0;i<nseg;++i) {
			p0 = seg(i).pnt(0);
			p1 = seg(i).pnt(1);
			gbl->res_r.v(p0,0) += 1./gbl->res.v(p1,0);
			gbl->res_r.v(p1,0) += 1./gbl->res.v(p0,0);
		}
		
		for(i=0;i<npnt;++i) {
			if (pnt(i).info == 0) {
				gbl->res.v(i,0) = 1./(gbl->res_r.v(i,0)/pnt(i).nnbor);
			}
		}
	}
	
	/* NOW RESCALE AT VERTICES */
	FLT maxlngth = 50.0;
	FLT minlngth = 0.0;
	for (int pind=0;pind<npnt;++pind) {
		lngth(pind) *= gbl->res.v(pind,0);
		lngth(pind) = MIN(lngth(pind),maxlngth);
		lngth(pind) = MAX(lngth(pind),minlngth);
	}
	
	/* LIMIT BOUNDARY CURVATURE */
	for(int i=0;i<nebd;++i) {
		if (!(hp_ebdry(i)->is_curved())) continue;
		
		for(int j=0;j<ebdry(i)->nseg;++j) {
			int sind = ebdry(i)->seg(j);
			int v1 = seg(sind).pnt(0);
			int v2 = seg(sind).pnt(1);
			
			crdtocht1d(sind);
			
			/* FIND ANGLE BETWEEN LINEAR SIDES */
			int tind = seg(sind).tri(0);
			int k;
			for(k=0;k<3;++k)
				if (tri(tind).seg(k) == sind) break;
			
			int v0 = tri(tind).pnt(k);
			
			TinyVector<FLT,ND> dx0;
			dx0(0) = pnts(v2)(0)-pnts(v1)(0);
			dx0(1) = pnts(v2)(1)-pnts(v1)(1);
			FLT length0 = dx0(0)*dx0(0) +dx0(1)*dx0(1);
			
			TinyVector<FLT,ND> dx1;
			dx1(0) = pnts(v0)(0)-pnts(v2)(0);
			dx1(1) = pnts(v0)(1)-pnts(v2)(1);
			FLT length1 = dx1(0)*dx1(0) +dx1(1)*dx1(1);
			
			TinyVector<FLT,ND> dx2;
			dx2(0) = pnts(v1)(0)-pnts(v0)(0);
			dx2(1) = pnts(v1)(1)-pnts(v0)(1);
			FLT length2 = dx2(0)*dx2(0) +dx2(1)*dx2(1);
			
			TinyVector<FLT,2> ep, dedpsi;
			basis::tri(log2p)->ptprobe1d(2,&ep(0),&dedpsi(0),-1.0,&cht(0,0),MXTM);
			FLT lengthept = dedpsi(0)*dedpsi(0) +dedpsi(1)*dedpsi(1);
			
			FLT ang1 = acos(-(dx0(0)*dx2(0) +dx0(1)*dx2(1))/sqrt(length0*length2));
			FLT curved1 = acos((dx0(0)*dedpsi(0) +dx0(1)*dedpsi(1))/sqrt(length0*lengthept));
			
			basis::tri(log2p)->ptprobe1d(2,&ep(0),&dedpsi(0),1.0,&cht(0,0),MXTM);
			lengthept = dedpsi(0)*dedpsi(0) +dedpsi(1)*dedpsi(1);
			
			FLT ang2 = acos(-(dx0(0)*dx1(0) +dx0(1)*dx1(1))/sqrt(length0*length1));
			FLT curved2 = acos((dx0(0)*dedpsi(0) +dx0(1)*dedpsi(1))/sqrt(length0*lengthept));                            
			
			// FIXME: end points are wrong for periodic boundary or communication boundary
			FLT sum = gbl->curvature_sensitivity*(fabs(curved1/ang1) +fabs(curved2/ang2));
			lngth(v1) /= 1. +sum;
			lngth(v2) /= 1. +sum;
		}
	}
	
	//	/* AVOID HIGH ASPECT RATIOS */
	//	int nsweep = 0;
	//	int count;
	//	do {
	//		count = 0;
	//		for(int i=0;i<nseg;++i) {
	//			int v0 = seg(i).pnt(0);
	//			int v1 = seg(i).pnt(1);
	//			FLT ratio = lngth(v1)/lngth(v0);
	//
	//			if (ratio > 3.0) {
	//				lngth(v1) = 2.5*lngth(v0);
	//				++count;
	//			}
	//			else if (ratio < 0.333) {
	//				lngth(v0) = 2.5*lngth(v1);
	//				++count;
	//			}
	//		}
	//		++nsweep;
	//		*gbl->log << "#aspect ratio fixes " << nsweep << ' ' << count << std::endl;
	//	} while(count > 0 && nsweep < 5);
	
	
	if (gbl->adapt_output) {
		ostringstream fname;
		fname << "adapt_diagnostic" << gbl->tstep << '_' << gbl->idprefix;
		output(fname.str(),tri_hp::adapt_diagnostic);
	}
	
	return;
}




static int error_count = 0;

void tri_hp::updatepdata(int v0) {
	int n,tind=-1,step; 
	FLT r,s;      

	bool found = gbl->pstr->findinteriorpt(pnts(v0),tind,r,s);
	if (!found) {
		*gbl->log << "Warning #" << error_count << ": didn't find interior point in updatepdata for " << v0 << ' ' << pnts(v0) << std::endl;
		//			std::ostringstream fname;
		//			fname << "current_solution" << error_count++ << '_' << gbl->idprefix;
		//			tri_mesh::output(fname.str().c_str(),tri_mesh::grid);
		//			tri_hp::output(fname.str().c_str(),tri_hp::tecplot);
	}
	
	
	// gbl->pstr->findandmvptincurved(pnts(v0),tind,r,s);  // old bad way

	for(step=0;step<gbl->nadapt;++step) {
		gbl->pstr->ugtouht(tind,step);
		basis::tri(log2p)->ptprobe(NV,&ugbd(step).v(v0,0),&gbl->pstr->uht(0)(0),MXTM);
	}

	if (gbl->pstr->tri(tind).info > -1) {
		for(step=1;step<gbl->nadapt;++step) {
			gbl->pstr->crdtocht(tind,step);
			basis::tri(log2p)->ptprobe_bdry(ND,&vrtxbd(step)(v0)(0),&gbl->pstr->cht(0,0),MXTM);
		}
	}
	else {
		for(step=1;step<gbl->nadapt;++step) {
			for(n=0;n<ND;++n) 
				vrtxbd(step)(v0)(n) = gbl->pstr->vrtxbd(step)(gbl->pstr->tri(tind).pnt(0))(n)*(s +1.)/2.
					+gbl->pstr->vrtxbd(step)(gbl->pstr->tri(tind).pnt(1))(n)*(-r -s)/2.
					+gbl->pstr->vrtxbd(step)(gbl->pstr->tri(tind).pnt(2))(n)*(r +1.)/2.;
		}
	}

	return;
}

void tri_hp::updatepdata_bdry(int bnum, int bel, int endpt) {
	int n,sind,sidloc,v0,step;
	FLT psi;

	v0 = seg(ebdry(bnum)->seg(bel)).pnt(endpt);
	gbl->pstr->hp_ebdry(bnum)->findandmovebdrypt(pnts(v0),sidloc,psi);
	sind = gbl->pstr->ebdry(bnum)->seg(sidloc);

	for(step=0;step<gbl->nadapt;++step) {
		gbl->pstr->ugtouht1d(sind,step);
		basis::tri(log2p)->ptprobe1d(NV,&ugbd(step).v(v0,0),&gbl->pstr->uht(0)(0),MXTM);
	}

	if (hp_ebdry(bnum)->is_curved()) {
		for(step=1;step<gbl->nadapt;++step) {
			gbl->pstr->crdtocht1d(sind,step);
			basis::tri(log2p)->ptprobe1d(ND,&vrtxbd(step)(v0)(0),&gbl->pstr->cht(0,0),MXTM);
		}
	}
	else {
		for(step=1;step<gbl->nadapt;++step) {
			for(n=0;n<ND;++n) 
				vrtxbd(step)(v0)(n) = gbl->pstr->vrtxbd(step)(gbl->pstr->seg(sind).pnt(0))(n)*(1. -psi)/2.
					+gbl->pstr->vrtxbd(step)(gbl->pstr->seg(sind).pnt(1))(n)*(1. +psi)/2.;
		}
	}

	/* FOR INTERNALLY STORED DATA */
	hp_ebdry(bnum)->updatepdata_bdry(bel,endpt,gbl->pstr->hp_ebdry(bnum));

	return;
}

void tri_hp::movepdata(int from, int to) {
	int n,step;

	for(step=0;step<gbl->nadapt;++step) {
		for(n=0;n<NV;++n)
			ugbd(step).v(to,n) = ugbd(step).v(from,n);

		for(n=0;n<ND;++n)
			vrtxbd(step)(to)(n) = vrtxbd(step)(from)(n);
	}

	return;
}

void tri_hp::movepdata_bdry(int bnum,int bel,int endpt) {
	/* This is just for internal data (if any) */
	hp_ebdry(bnum)->movepdata_bdry(bel,endpt,gbl->pstr->hp_ebdry(bnum));
}


void tri_hp::updatesdata(int sind) {
	int i,m,n,v0,v1,step,info;
	FLT r,s,upt[NV];
	char uplo[] = "U";
	TinyVector<FLT,2> pt;
	bool found;
	int tind = -1;

	if (!sm0) return;

	v0 = seg(sind).pnt(0);
	v1 = seg(sind).pnt(1);

	for(n=0;n<ND;++n)
		basis::tri(log2p)->proj1d(pnts(v0)(n),pnts(v1)(n),&crd(n)(0,0));

	for(step=0;step<gbl->nadapt;++step)
		for(n=0;n<NV;++n)
			basis::tri(log2p)->proj1d(ugbd(step).v(v0,n),ugbd(step).v(v1,n),&bdwk(step,n)(0,0));

	for(i=0;i<basis::tri(log2p)->gpx();++i) {
		pt(0) = crd(0)(0,i);
		pt(1) = crd(1)(0,i);
		found = gbl->pstr->findinteriorpt(pt,tind,r,s);
		if (!found) {
			*gbl->log << "Warning #" << error_count << ": didn't find interior point in updatesdata for " << sind << ' ' << pt << std::endl;
//			std::ostringstream fname;
//			fname << "current_solution" << error_count++ << '_' << gbl->idprefix;
//			tri_mesh::output(fname.str().c_str(),tri_mesh::grid);
//			tri_hp::output(fname.str().c_str(),tri_hp::tecplot);
		}

		for(step=0;step<gbl->nadapt;++step) {
			gbl->pstr->ugtouht(tind,step);
			basis::tri(log2p)->ptprobe(NV,upt,&gbl->pstr->uht(0)(0),MXTM);
			for(n=0;n<NV;++n)    {
				bdwk(step,n)(0,i) -= upt[n];
			}
		}
	}

	for(step=0;step<gbl->nadapt;++step) {
		for(n=0;n<NV;++n)
			basis::tri(log2p)->intgrt1d(&lf(n)(0),&bdwk(step,n)(0,0));

		for(n=0;n<NV;++n) {
			PBTRS(uplo,basis::tri(log2p)->sm(),basis::tri(log2p)->sbwth(),1,(double *) &basis::tri(log2p)->sdiag1d(0,0),basis::tri(log2p)->sbwth()+1,&lf(n)(2),basis::tri(log2p)->sm(),info);
			for(m=0;m<basis::tri(log2p)->sm();++m) 
				ugbd(step).s(sind,m,n) = -lf(n)(2+m);
		}
	}
	return;
}

void tri_hp::updatesdata_bdry(int bnum,int bel) {
	int m,n,sind,v0,v1,step,stgt,info;
	TinyVector<FLT,2> pt;
	FLT psi;
	FLT upt[NV];
	char uplo[] = "U";

	if (!sm0) return;

	sind = ebdry(bnum)->seg(bel);
	v0 = seg(sind).pnt(0);
	v1 = seg(sind).pnt(1);

	for(step=0;step<gbl->nadapt;++step)
		for(n=0;n<NV;++n)
			basis::tri(log2p)->proj1d(ugbd(step).v(v0,n),ugbd(step).v(v1,n),&bdwk(step,n)(0,0));

	for(step=0;step<gbl->nadapt;++step)
		for(n=0;n<ND;++n)
			basis::tri(log2p)->proj1d(vrtxbd(step)(v0)(n),vrtxbd(step)(v1)(n),&bdwk(step,n)(1,0));


	if (hp_ebdry(bnum)->is_curved()) {

		for(m=0;m<basis::tri(log2p)->gpx();++m) {
			pt(0) = bdwk(0,0)(1,m);
			pt(1) = bdwk(0,1)(1,m);
			gbl->pstr->hp_ebdry(bnum)->findandmovebdrypt(pt,stgt,psi);
			stgt = gbl->pstr->ebdry(bnum)->seg(stgt);

			for(step=0;step<gbl->nadapt;++step) {
				gbl->pstr->ugtouht1d(stgt,step);
				basis::tri(log2p)->ptprobe1d(NV,upt,&gbl->pstr->uht(0)(0),MXTM);
				for(n=0;n<NV;++n)    
					bdwk(step,n)(0,m) -= upt[n];

				gbl->pstr->crdtocht1d(stgt,step);
				basis::tri(log2p)->ptprobe1d(ND,upt,&gbl->pstr->cht(0,0),MXTM);
				for(n=0;n<ND;++n)    
					bdwk(step,n)(1,m) -= upt[n];
			}                          
		}      

		for(step=0;step<gbl->nadapt;++step) {
			for(n=0;n<ND;++n) {
				basis::tri(log2p)->intgrt1d(&lf(n)(0),&bdwk(step,n)(1,0));
				PBTRS(uplo,basis::tri(log2p)->sm(),basis::tri(log2p)->sbwth(),1,(double *) &basis::tri(log2p)->sdiag1d(0,0),basis::tri(log2p)->sbwth()+1,&lf(n)(2),basis::tri(log2p)->sm(),info);

				for(m=0;m<basis::tri(log2p)->sm();++m)
					hp_ebdry(bnum)->crdsbd(step,bel,m,n) = -lf(n)(m+2);
			}
		}
	}
	else {
		for(m=0;m<basis::tri(log2p)->gpx();++m) {
			pt(0) = bdwk(0,0)(1,m);
			pt(1) = bdwk(0,1)(1,m);

			/* FIND PSI */                
			gbl->pstr->hp_ebdry(bnum)->findandmovebdrypt(pt,stgt,psi);
			stgt = gbl->pstr->ebdry(bnum)->seg(stgt);

			/* CALCULATE VALUE OF SOLUTION AT POINT */
			for(step=0;step<gbl->nadapt;++step) {
				gbl->pstr->ugtouht1d(stgt,step);
				basis::tri(log2p)->ptprobe1d(NV,upt,&gbl->pstr->uht(0)(0),MXTM);
				for(n=0;n<NV;++n)    
					bdwk(step,n)(0,m) -= upt[n];
			}
		}
	}

	for(step=0;step<gbl->nadapt;++step) {
		for(n=0;n<NV;++n)
			basis::tri(log2p)->intgrt1d(&lf(n)(0),&bdwk(step,n)(0,0));

		for(n=0;n<NV;++n) {
			PBTRS(uplo,basis::tri(log2p)->sm(),basis::tri(log2p)->sbwth(),1,(double *) &basis::tri(log2p)->sdiag1d(0,0),basis::tri(log2p)->sbwth()+1,&lf(n)(2),basis::tri(log2p)->sm(),info);
			for(m=0;m<basis::tri(log2p)->sm();++m) {
				ugbd(step).s(sind,m,n) = -lf(n)(2+m);
			}
		}
	}

	/* UPDATE INTERNAL INFORMATION */
	hp_ebdry(bnum)->updatesdata_bdry(bel,gbl->pstr->hp_ebdry(bnum));

	return;
}

void tri_hp::movesdata(int from, int to) {
	int step;

	if (!sm0) return;

	for(step=0;step<gbl->nadapt;++step)
		ugbd(step).s(to,Range::all(),Range::all()) = ugbd(step).s(from,Range::all(),Range::all());

	return;
}

void tri_hp::movesdata_bdry(int bnum,int bel) {

	hp_ebdry(bnum)->movesdata_bdry(bel,gbl->pstr->hp_ebdry(bnum));

	return;
}

void tri_hp::updatetdata(int tind) {
	int i,j,n,step,info;
	FLT r,s;
	FLT upt[NV];
	char uplo[] = "U";
	TinyVector<FLT,2> pt;
	bool found;
	int ttgt = -1;

	if (!im0) return;  /* FIXME NEED TO FIX THIS IN MESH SO CAN TURN OFF ENTIRE LOOP */

	for(step=0;step<gbl->nadapt;++step) {
		ugtouht_bdry(tind,step);
		for(n=0;n<NV;++n)
			basis::tri(log2p)->proj_bdry(&uht(n)(0),&bdwk(step,n)(0,0),MXGP);
	}

	crdtocht(tind);
	for(n=0;n<ND;++n)
		basis::tri(log2p)->proj_bdry(&cht(n,0),&crd(n)(0,0),MXGP);

	for (i=0; i < basis::tri(log2p)->gpx(); ++i ) {
		for (j=0; j < basis::tri(log2p)->gpn(); ++j ) {
			pt(0) = crd(0)(i,j);
			pt(1) = crd(1)(i,j);
			found = gbl->pstr->findinteriorpt(pt,ttgt,r,s);
			if (!found) {
				*gbl->log << "Warning #" << error_count << ": didn't find interior point in updatetdata for " << tind << ' ' << pt << std::endl;
				*gbl->log << "Using triangle " << ttgt << " with (r,s) = (" << r << ',' << s << ')' << std::endl;
//				std::ostringstream fname;
//				fname << "current_solution" << error_count++ << '_' << gbl->idprefix;
//				tri_mesh::output(fname.str().c_str(),tri_mesh::grid);
//				tri_hp::output(fname.str().c_str(),tri_hp::tecplot);
			}            
			for(step=0;step<gbl->nadapt;++step) {
				gbl->pstr->ugtouht(ttgt,step);
				basis::tri(log2p)->ptprobe(NV,upt,&gbl->pstr->uht(0)(0),MXTM);
				for(n=0;n<NV;++n)
					bdwk(step,n)(i,j) -= upt[n];
			}
		}
	}

	for(step=0;step<gbl->nadapt;++step) {
		for(n=0;n<NV;++n) {
			basis::tri(log2p)->intgrt(lf(n).data(),&bdwk(step,n)(0,0),MXGP);
			PBTRS(uplo,basis::tri(log2p)->im(),basis::tri(log2p)->ibwth(),1,(double *) &basis::tri(log2p)->idiag(0,0),basis::tri(log2p)->ibwth()+1,&lf(n)(basis::tri(log2p)->bm()),basis::tri(log2p)->im(),info);
			for(i=0;i<basis::tri(log2p)->im();++i)
				ugbd(step).i(tind,i,n) = -lf(n)(basis::tri(log2p)->bm()+i);
		}
	}

	return;
}

void tri_hp::movetdata(int from, int to) {
	int step;

	if (!im0) return;  /* FIXME NEED TO FIX THIS IN MESH SO CAN TURN OFF ENTIRE LOOP */

	for(step=0;step<gbl->nadapt;++step) {
		ugbd(step).i(to,Range::all(),Range::all()) = ugbd(step).i(from,Range::all(),Range::all());
	}

	return;
}




