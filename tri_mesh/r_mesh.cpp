#include "r_tri_mesh.h"
#include "r_tri_boundary.h"
#include "block.h"
#include <utilities.h>
#include <iostream>
#include <cmath>
#include <input_map.h>
#include <fstream>


void r_tri_mesh::init(input_map& input, void *gin) {
	std::string keyword;
	std::istringstream data;
	std::string filename;
	int ival;

	tri_mesh::init(input,gin);
	gbl = static_cast<global *>(gin);

	keyword = gbl->idprefix + "_r_fadd";
	if (!input.get(keyword,fadd)) {
		input.getwdefault("r_fadd",fadd,1.0);
	}

	keyword = gbl->idprefix + "_r_cfl";
	if (!input.get(keyword,r_cfl)) {
		input.getwdefault("r_cfl",r_cfl,0.5);
	}

	keyword = gbl->idprefix + "_r_output_type";
	if (input.get(keyword,ival)) {
		output_type = static_cast<tri_mesh::filetype>(ival);
	}
	else {
		if (input.get("r_output_type",ival)) {
			output_type = static_cast<tri_mesh::filetype>(ival);
		}
		else {
			output_type = tri_mesh::grid;
		}
	}

	/* local storage */
	ksprg.resize(maxpst);
	kvol.resize(maxpst);
	src.resize(maxpst);
	isfrst = false;

	/* BLOCK SHARED INFORMATION */
	gbl->diag.resize(maxpst);
	gbl->res.resize(maxpst);
	gbl->res1.resize(maxpst);

	r_sbdry.resize(nebd);
	for(int i=0;i<nebd;++i)
		r_sbdry(i) = getnewedgeobject(i,input);

	r_vbdry.resize(nvbd);
	for(int i=0;i<nvbd;++i)
		r_vbdry(i) = getnewvrtxobject(i,input);

	return;
}

void r_tri_mesh::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
	std::string keyword;
	std::istringstream data;
	std::string filename;

	tri_mesh::init(in,why,sizereduce1d);
	const r_tri_mesh& inmesh = dynamic_cast<const r_tri_mesh &>(in);
	gbl = inmesh.gbl;
	fadd = inmesh.fadd;
	r_cfl = inmesh.r_cfl;

	/* local storage */
	ksprg.resize(maxpst);
	kvol.resize(maxpst);
	src.resize(maxpst);
	pnts_frst.resize(maxpst);
	isfrst = false;

	r_sbdry.resize(nebd);
	for(int i=0;i<nebd;++i)
		r_sbdry(i) = inmesh.r_sbdry(i)->create(*this,*ebdry(i));

	r_vbdry.resize(nebd);
	for(int i=0;i<nvbd;++i)
		r_vbdry(i) = inmesh.r_vbdry(i)->create(*this,*vbdry(i));

	return;
}


r_tri_mesh::~r_tri_mesh() {
	for(int i=0;i<nebd;++i)
		delete r_sbdry(i);

	for(int i=0;i<nvbd;++i)
		delete r_vbdry(i);
}


void r_tri_mesh::rklaplace() {
	int sind,tind,p0,p1,k;
	FLT dx,dy,l;

	for(sind=0;sind<nseg;++sind)
		ksprg(sind) = 0.0;

	/* COEFFICIENTS FOR LAPLACE EQUATION */
	/* THIS REQUIRES 2 EVALUATIONS OF SIDE LENGTH FOR EACH SIDE */
	/* BUT IS LOGISTICALLY SIMPLE          */
	for(tind=0;tind<ntri;++tind) {
		for(k=0;k<3;++k) {
			sind = tri(tind).seg(k);
			p0 = seg(sind).pnt(0);
			p1 = seg(sind).pnt(1);
			dx = pnts(p1)(0) -pnts(p0)(0);
			dy = pnts(p1)(1) -pnts(p0)(1);
			l  = (dx*dx +dy*dy)/area(tind);

			ksprg(sind) -= l;
			sind = tri(tind).seg((k+1)%3);
			ksprg(sind) += l;
			sind = tri(tind).seg((k+2)%3);
			ksprg(sind) += l;
		}
	}

	return;
}

#ifdef FOURTH
void r_tri_mesh::calc_kvol() {
	int last_phase, mp_phase;

	for(int i=0;i<npnt;++i)
		kvol(i) = 0.0;

	for(int tind=0;tind<ntri;++tind)
		for(int i=0;i<3;++i)
			kvol(tri(tind).pnt(i)) += area(tind);

	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		pmsgload(boundary::all_phased,mp_phase,boundary::symmetric,kvol.data(),0,0,1);
		pmsgpass(boundary::all_phased,mp_phase, boundary::symmetric);
		last_phase = true;
		last_phase &= pmsgwait_rcv(boundary::all_phased,mp_phase, boundary::symmetric, boundary::average, kvol.data(),0,0,1);
	}

	for(int i=0;i<npnt;++i)
		kvol(i) = 1./kvol(i);
}
#endif

void r_tri_mesh::rksprg() {
	int sind,p0,p1;
	double dx,dy;

	/* 2D SPRING CONSTANTS FINE MESH*/
	for(sind=0;sind<nseg;++sind) {
		p0 = seg(sind).pnt(0);
		p1 = seg(sind).pnt(1);
		dx = pnts(p1)(0) -pnts(p0)(0);
		dy = pnts(p1)(1) -pnts(p0)(1);
		ksprg(sind) = 1.0/(dx*dx +dy*dy);
	}

	return;
}

void r_tri_mesh::rkmgrid() {
	int i,j,sind,tind,tind0,tind1,p0,p1;

	r_tri_mesh* fmesh = dynamic_cast<r_tri_mesh *>(fine);

	/* Load Diag Locally */
	Array<FLT,1> diag;
	diag.reference(gbl->diag);
	Array<TinyVector<FLT,ND>,1> res1;
	res1.reference(gbl->res1);

	/* TEMPORARILY USE DIAG TO STORE DIAGONAL SUM */
	for(i=0;i<fmesh->npnt;++i)
		diag(i) = 0.0;

	/* FORM KIJ SUM AT POINTS */
	for(sind=0;sind<fmesh->nseg;++sind) {
		p0 = fmesh->seg(sind).pnt(0);
		p1 = fmesh->seg(sind).pnt(1);
		diag(p0) += fmesh->ksprg(sind);
		diag(p1) += fmesh->ksprg(sind);
	}

	for(i=0;i<nseg;++i)
		ksprg(i) = 0.0;

	/* LOOP THROUGH FINE POINTS    */
	/* TO CALCULATE KSPRG ON COARSE MESH */
	for(i=0;i<fmesh->npnt;++i) {
		tind = fmesh->ccnnct(i).tri;
		for(j=0;j<3;++j) {
			sind = tri(tind).seg(j);
			ksprg(sind) -= fmesh->ccnnct(i).wt(j)*fmesh->ccnnct(i).wt((j+1)%3)*diag(i);
		}
	}

	/* LOOP THROUGH FINE SIDES */
	for(i=0;i<fmesh->nseg;++i) {
		p0 = fmesh->seg(i).pnt(0);
		p1 = fmesh->seg(i).pnt(1);
		tind0 = fmesh->ccnnct(p0).tri;
		tind1 = fmesh->ccnnct(p1).tri;

		/* TEMPORARILY STORE WEIGHTS FOR FINE POINTS (0,1) */
		/* FOR EACH COARSE VERTEX */
		for(j=0;j<3;++j)  {
			res1(tri(tind1).pnt(j))(0) = 0.0;
			res1(tri(tind0).pnt(j))(1) = 0.0;
		}

		for(j=0;j<3;++j)  {
			res1(tri(tind0).pnt(j))(0) = fmesh->ccnnct(p0).wt(j);
			res1(tri(tind1).pnt(j))(1) = fmesh->ccnnct(p1).wt(j);
		}

		/* LOOP THROUGH COARSE TRIANGLE 0 SIDES */
		for(j=0;j<3;++j) {
			sind = tri(tind0).seg(j);
			ksprg(sind) += fmesh->ksprg(i)*
				(res1(seg(sind).pnt(0))(0)*res1(seg(sind).pnt(1))(1)
				+res1(seg(sind).pnt(1))(0)*res1(seg(sind).pnt(0))(1));
		}

		if (tind0 != tind1) {
			for(j=0;j<3;++j) {
				sind = tri(tind1).seg(j);
				if (seg(sind).tri(0) +seg(sind).tri(1) != tind0 +tind1) {
					ksprg(sind) += fmesh->ksprg(i)*
						(res1(seg(sind).pnt(0))(0)*res1(seg(sind).pnt(1))(1)
						+res1(seg(sind).pnt(1))(0)*res1(seg(sind).pnt(0))(1));

				}
			}
		}
	}

	return;
}


void r_tri_mesh::update() {
	int i,n;

	r_tri_mesh::rsdl();

	Array<FLT,1> diag;
	diag.reference(gbl->diag);
	Array<TinyVector<FLT,ND>,1> res;
	res.reference(gbl->res);

	for(i=0;i<npnt;++i)
		for(n=0;n<ND;++n)
			pnts(i)(n) -= diag(i)*res(i)(n);
}

void r_tri_mesh::zero_source() {
	int i,n;

	for(i=0;i<npnt;++i)
		for(n=0;n<ND;++n)
			src(i)(n) = 0.0;

	return;
}

void r_tri_mesh::sumsrc() {
	int i,n;

	Array<TinyVector<FLT,ND>,1> res;
	res.reference(gbl->res);

	for(i=0;i<npnt;++i)
		for(n=0;n<ND;++n)
			src(i)(n) = -1.0*res(i)(n);

	return;
}


void r_tri_mesh::mg_restrict() {
	int i,j,n,tind,p0;
	r_tri_mesh *fmesh = dynamic_cast<r_tri_mesh *>(fine);

	Array<TinyVector<FLT,ND>,1> fres(gbl->res);
	Array<transfer,1> fccnnct(fmesh->ccnnct);
	int fnvrtx = fmesh->npnt;

	for(i=0;i<npnt;++i)
		for(n=0;n<ND;++n)
			src(i)(n) = 0.0;

	/* LOOP THROUGH FINE POINTS TO CALCULATE RESIDUAL  */
	for(i=0;i<fnvrtx;++i) {
		tind = fccnnct(i).tri;
		for(j=0;j<3;++j) {
			p0 = tri(tind).pnt(j);
			for(n=0;n<ND;++n)
				src(p0)(n) += fadd*fccnnct(i).wt(j)*fres(i)(n);
		}
	}
	
	/* Need to communicate for partition boundaries */
	for(int last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		for(i=0;i<nvbd;++i)
			vbdry(i)->vloadbuff(boundary::partitions,(FLT *) src.data(),0,ND-1,ND);
		for(i=0;i<nebd;++i)
			ebdry(i)->vloadbuff(boundary::partitions,(FLT *) src.data(),0,ND-1,ND);
		
		for(i=0;i<nebd;++i) 
			ebdry(i)->comm_prepare(boundary::partitions,mp_phase,boundary::symmetric);
		for(i=0;i<nvbd;++i)
			vbdry(i)->comm_prepare(boundary::partitions,0,boundary::symmetric);
		
		for(i=0;i<nebd;++i) 
			ebdry(i)->comm_exchange(boundary::partitions,mp_phase,boundary::symmetric);
		for(i=0;i<nvbd;++i)
			vbdry(i)->comm_exchange(boundary::partitions,0,boundary::symmetric);
		
		last_phase = true;
		for(i=0;i<nebd;++i) {
			last_phase &= ebdry(i)->comm_wait(boundary::partitions,mp_phase,boundary::symmetric);
			ebdry(i)->vfinalrcv(boundary::partitions,mp_phase,boundary::symmetric,boundary::sum,(FLT *) src.data(),0,ND-1,ND);
		}
		for(i=0;i<nvbd;++i) {
			vbdry(i)->comm_wait(boundary::partitions,0,boundary::symmetric);
			vbdry(i)->vfinalrcv(boundary::partitions,0,boundary::symmetric,boundary::average,(FLT *) src.data(),0,ND-1,ND);
		}
	}	

	/* LOOP THROUGH fv_to_ct POINTS    */
	/* TO CALCULATE POINT ON fv_to_ct MESH */
	Array<TinyVector<FLT,tri_mesh::ND>,1> fpnts(fmesh->pnts);
	Array<tristruct,1> ftri(fmesh->tri);
	for(i=0;i<npnt;++i) {
		tind = fcnnct(i).tri;

		for(n=0;n<ND;++n)
			pnts(i)(n) = 0.0;

		for(j=0;j<3;++j) {
			for(n=0;n<ND;++n)
				pnts(i)(n) += fcnnct(i).wt(j)*fpnts(ftri(tind).pnt(j))(n);
		}
	}

	for(i=0;i<npnt;++i) {
		for(n=0;n<ND;++n)
			pnts_frst(i)(n) = pnts(i)(n);
	}
	isfrst = true;

	return;
}

void r_tri_mesh::mg_prolongate() {
	int i,j,n,ind,tind;
	
	/* DETERMINE CORRECTIONS    */
	for(i=0;i<npnt;++i)
		for(n=0;n<ND;++n)
			pnts_frst(i)(n) -= pnts(i)(n);

	/* LOOP THROUGH FINE POINTS    */
	/* TO DETERMINE CHANGE IN SOLUTION */
	r_tri_mesh *fmesh = dynamic_cast<r_tri_mesh *>(fine);
	int fnpnt = fmesh->npnt;
	Array<TinyVector<FLT,ND>,1> res(gbl->res);

	for(i=0;i<fnpnt;++i) {

		for(n=0;n<ND;++n)
			res(i)(n) = 0.0;

		tind = fmesh->ccnnct(i).tri;

		for(j=0;j<3;++j) {
			ind = tri(tind).pnt(j);
			for(n=0;n<ND;++n)
				res(i)(n) -= fmesh->ccnnct(i).wt(j)*pnts_frst(ind)(n);
		}
	}
	
	/* COMMUNICATION FOR PARTITION BOUNDARIES */
	for(int last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		for(i=0;i<nvbd;++i)
			fmesh->vbdry(i)->vloadbuff(boundary::partitions,(FLT *) res.data(),0,ND-1,ND);
		for(i=0;i<nebd;++i)
			fmesh->ebdry(i)->vloadbuff(boundary::partitions,(FLT *) res.data(),0,ND-1,ND);
		
		for(i=0;i<nebd;++i) 
			fmesh->ebdry(i)->comm_prepare(boundary::partitions,mp_phase,boundary::symmetric);
		for(i=0;i<nvbd;++i)
			fmesh->vbdry(i)->comm_prepare(boundary::partitions,0,boundary::symmetric);
		
		for(i=0;i<nebd;++i) 
			fmesh->ebdry(i)->comm_exchange(boundary::partitions,mp_phase,boundary::symmetric);
		for(i=0;i<nvbd;++i)
			fmesh->vbdry(i)->comm_exchange(boundary::partitions,0,boundary::symmetric);
		
		last_phase = true;
		for(i=0;i<nebd;++i) {
			last_phase &= fmesh->ebdry(i)->comm_wait(boundary::partitions,mp_phase,boundary::symmetric);
			fmesh->ebdry(i)->vfinalrcv(boundary::partitions,mp_phase,boundary::symmetric,boundary::sum,(FLT *) res.data(),0,ND-1,ND);
		}
		for(i=0;i<nvbd;++i) {
			fmesh->vbdry(i)->comm_wait(boundary::partitions,0,boundary::symmetric);
			fmesh->vbdry(i)->vfinalrcv(boundary::partitions,0,boundary::symmetric,boundary::average,(FLT *) res.data(),0,ND-1,ND);
		}
	}

	Array<TinyVector<FLT,tri_mesh::ND>,1> fpnts(fmesh->pnts);
	for(i=0;i<fnpnt;++i)
		for(n=0;n<ND;++n)
			fpnts(i)(n) += res(i)(n);

	return;
}

FLT r_tri_mesh::maxres() {
	int i,n;
	FLT mxr[ND];
	FLT sum;

	Array<TinyVector<FLT,ND>,1> res;
	res.reference(gbl->res);

	for(n=0;n<ND;++n)
		mxr[n] = 0.0;

	for(i=0;i<npnt;++i)
		for(n=0;n<ND;++n)
			mxr[n] = MAX(mxr[n],fabs(res(i)(n)));

	sum = 0.0;
	for(n=0;n<ND;++n) {
		*gbl->log << ' ' << mxr[n] << ' ';
		sum += mxr[n];
	}


	return(sum);
}


void r_tri_mesh::tadvance() {
	if (!coarse_level) {
		rklaplace();
#ifdef FOURTH
		calc_kvol();
#endif
		zero_source();
		r_tri_mesh::setup_preconditioner();
		r_tri_mesh::rsdl();
		sumsrc();
		moveboundaries();
	}
	else {
#ifdef GEOMETRIC
		rklaplace();
#else
		/* USE MULTIGRID INTERPOLATION (ALGEBRAIC) */
		/* MUST BE DONE THIS WAY FOR SPRING METHOD */
		/* SETUP FIRST MESH */
		rkmgrid(fv_to_ct,fmesh);
#endif

#ifdef  FOURTH
		calc_kvol();
#endif
	}

	return;
}

void r_tri_mesh::moveboundaries() {

	/* MOVE BOUNDARY POSITIONS */
	for(int i=0;i<nebd;++i)
		r_sbdry(i)->tadvance();

	for(int i=0;i<nvbd;++i)
		r_vbdry(i)->tadvance();

	return;
}

void r_tri_mesh::rsdl() {
	int last_phase, mp_phase;
	int i,n,p0,p1;
	FLT dx,dy;
	Array<TinyVector<FLT,ND>,1> res;
	res.reference(gbl->res);

#ifdef FOURTH
	/*************************************/
	/* FOURTH ORDER MESH MOVEMENT SCHEME */
	/*************************************/
	Array<TinyVector<FLT,ND>,1> res1;
	res1.reference(gbl->res1);

	for(i=0;i<npnt;++i)
		for(n=0;n<ND;++n)
			res1(i)(n) = 0.0;

	for (int sind = 0; sind < nseg; ++sind) {
		p0 = seg(sind).pnt(0);
		p1 = seg(sind).pnt(1);

		dx = ksprg(sind)*(pnts(p1)(0)-pnts(p0)(0));
		dy = ksprg(sind)*(pnts(p1)(1)-pnts(p0)(1));

		res1(p0)(0) -= dx;
		res1(p0)(1) -= dy;

		res1(p1)(0) += dx;
		res1(p1)(1) += dy;
	}

	/* APPLY DIRICHLET BOUNDARY CONDITIONS */
	for(i=0;i<nebd;++i)
		r_sbdry(i)->fixdx2();

	for(i=0;i<nvbd;++i)
		r_vbdry(i)->fixdx2();

	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		pmsgload(boundary::all_phased, mp_phase, boundary::symmetric,(FLT *) gbl->res1.data(),0,1,2);
		pmsgpass(boundary::all_phased, mp_phase, boundary::symmetric);
		last_phase = true;
		last_phase &= pmsgwait_rcv(boundary::all_phased, mp_phase/3, boundary::symmetric,  boundary::average, (FLT *) gbl->res1.data(),0,1,2);
	}

	/* DIVIDE BY VOLUME FOR AN APPROXIMATION TO D^2/DX^2 */
	for(i=0;i<npnt;++i)
		for(n=0;n<ND;++n)
			res1(i)(n) *= kvol(i);

	for(i=0;i<npnt;++i)
		for(n=0;n<ND;++n)
			res(i)(n) = 0.0;

	for (int sind = 0; sind < nseg; ++sind) {
		p0 = seg(sind).pnt(0);
		p1 = seg(sind).pnt(1);

		dx = ksprg(sind)*(res1(p1)(0)-res1(p0)(0));
		dy = ksprg(sind)*(res1(p1)(1)-res1(p0)(1));

		res(p0)(0) -= dx;
		res(p0)(1) -= dy;

		res(p1)(0) += dx;
		res(p1)(1) += dy;
	}
#else
	res(Range(0,npnt-1)) = 0.0;

	int lnside = nseg;
	FLT lksprg;
	for(int sind=0;sind<lnside;++sind) {
		p0 = seg(sind).pnt(0);
		p1 = seg(sind).pnt(1);

		lksprg = ksprg(sind);
		dx = lksprg*(pnts(p1)(0)-pnts(p0)(0));
		dy = lksprg*(pnts(p1)(1)-pnts(p0)(1));

		res(p0)(0) -= dx;
		res(p0)(1) -= dy;

		res(p1)(0) += dx;
		res(p1)(1) += dy;
	}
#endif

	/* CALCULATE DRIVING TERM ON FIRST ENTRY TO COARSE MESH */
	if (isfrst) {
		for(i=0;i<npnt;++i)
			for(n=0;n<ND;++n)
				src(i)(n) -= res(i)(n);

		isfrst = false;
	}

	/* ADD IN MULTIGRID SOURCE OR FMESH SOURCE */
	for(i=0;i<npnt;++i)
		for(n=0;n<ND;++n)
			res(i)(n) += src(i)(n);

	/* APPLY DIRICHLET BOUNDARY CONDITIONS */
	for(i=0;i<nebd;++i)
		r_sbdry(i)->dirichlet();

	/* APPLY DIRICHLET BOUNDARY CONDITIONS */
	for(i=0;i<nvbd;++i)
		r_vbdry(i)->dirichlet();

	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		pmsgload(boundary::all_phased,mp_phase, boundary::symmetric,(FLT *) gbl->res.data(),0,1,2);
		pmsgpass(boundary::all_phased,mp_phase, boundary::symmetric);
		last_phase = true;
		last_phase &= pmsgwait_rcv(boundary::all_phased,mp_phase, boundary::symmetric, boundary::average, (FLT *) gbl->res.data(),0,1,2);
	}

	return;
}

void r_tri_mesh::element_jacobian(int tind, Array<FLT,2> K) {

#ifdef FOURTH
	*gbl->log << "DIDN'T DO THIS FOR ELEMENT JACOBIAN\n";
	sim::abort(__LINE__,__FILE__,gbl->log);
#else
	TinyMatrix<int,3,2> seg_pnts;
	seg_pnts =	1,2,
							2,0,
							0,1;
	
	K = 0;
	for(int s=0;s<3;++s) {
		int sind = tri(tind).seg(s);
		FLT lksprg = ksprg(sind);
		int tadj = seg(sind).tri(1);
		if (tadj >= 0) 
			lksprg *= 0.5;  // Hack for the fact that element jacobian hits each side twice
		
		K(seg_pnts(s,0)*ND,seg_pnts(s,0)*ND) += lksprg;
		K(seg_pnts(s,1)*ND,seg_pnts(s,1)*ND) += lksprg;
		K(seg_pnts(s,0)*ND,seg_pnts(s,1)*ND) -= lksprg;
		K(seg_pnts(s,1)*ND,seg_pnts(s,0)*ND) -= lksprg;	
		
		K(seg_pnts(s,0)*ND+1,seg_pnts(s,0)*ND+1) += lksprg;
		K(seg_pnts(s,1)*ND+1,seg_pnts(s,1)*ND+1) += lksprg;
		K(seg_pnts(s,0)*ND+1,seg_pnts(s,1)*ND+1) -= lksprg;
		K(seg_pnts(s,1)*ND+1,seg_pnts(s,0)*ND+1) -= lksprg;	
		
		
	}
#endif
	
	return;
}



void r_tri_mesh::setup_preconditioner() {
	int last_phase, mp_phase;
	int i,p0,p1,sind;
	Array<FLT,1> diag;
	diag.reference(gbl->diag);

#ifdef FOURTH
	/**************************************************/
	/* DETERMINE MESH MOVEMENT TIME STEP              */
	/**************************************************/
	for(i=0;i<npnt;++i)
		diag(i) = 0.0;

	for(sind=0;sind<nseg;++sind) {
		p0 = seg(sind).pnt(0);
		p1 = seg(sind).pnt(1);
		diag(p0) += fabs(ksprg(sind));
		diag(p1) += fabs(ksprg(sind));
	}


	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		pmsgload(boundary::all_phased, mp_phase, boundary::symmetric,gbl->diag.data(),0,0,1);
		pmsgpass(boundary::all_phased, mp_phase, boundary::symmetric);
		last_phase = true;
		last_phase &= pmsgwait_rcv(boundary::all_phased, mp_phase, boundary::symmetric,  boundary::average, gbl->diag.data(),0,0,1);
	}

	Array<TinyVector<FLT,ND>,1> res1;
	res1.reference(gbl->res1);

	for(i=0;i<npnt;++i)
		res1(i)(0) = diag(i)*kvol(i);

	for(i=0;i<npnt;++i)
		diag(i) = 0.0;

	for(sind=0;sind<nseg;++sind) {
		p0 = seg(sind).pnt(0);
		p1 = seg(sind).pnt(1);
		diag(p0) += fabs(ksprg(sind))*(res1(p0)(0) +fabs(ksprg(sind))*kvol(p1));
		diag(p1) += fabs(ksprg(sind))*(res1(p1)(0) +fabs(ksprg(sind))*kvol(p0));
	}
#else
	/**************************************************/
	/* DETERMINE MESH MOVEMENT TIME STEP              */
	/**************************************************/
	for(i=0;i<npnt;++i)
		diag(i) = 0.0;

	/* FORM TIME STEP FOR MV_UPDATE */
	for(sind=0;sind<nseg;++sind) {
		p0 = seg(sind).pnt(0);
		p1 = seg(sind).pnt(1);
		diag(p0) += fabs(ksprg(sind));
		diag(p1) += fabs(ksprg(sind));
	}
#endif

	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		pmsgload(boundary::all_phased,mp_phase, boundary::symmetric,gbl->diag.data(),0,0,1);
		pmsgpass(boundary::all_phased,mp_phase, boundary::symmetric);
		last_phase = true;
		last_phase &= pmsgwait_rcv(boundary::all_phased,mp_phase, boundary::symmetric, boundary::average, gbl->diag.data(),0,0,1);
	}

	for(i=0;i<npnt;++i)
		diag(i) = r_cfl/diag(i);

	return;
}
