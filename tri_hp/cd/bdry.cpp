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

//#define MPDEBUG

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
#ifdef MESH_REF_VEL
						mvel(n) += x.gbl->mesh_ref_vel(n);
#endif
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


void melt::init(input_map& inmap,void* gbl_in) {
	std::string keyword;
	
	keyword = base.idprefix + "_curved";
	inmap[keyword] = "1";
	
	keyword = base.idprefix + "_coupled";
	inmap[keyword] = "1";
	
	dirichlet::init(inmap,gbl_in);

	return;
}


void melt::tadvance() {
	dirichlet::tadvance();
	
	/* THIS IS TO RECEIVE MESSAGES SENT FROM SURFACE DURING TADVANCE RSDL CALL */
	if (x.gbl->substep == 0) {
		if (!x.coarse_level) {
			rsdl(x.gbl->nstage);
		}
	}
	
	return;
}

void melt::rsdl(int stage) {
	
#ifdef petsc
	/* Store 0.0 in vertex residual in r_mesh residual vector */
	r_tri_mesh::global *r_gbl = dynamic_cast<r_tri_mesh::global *>(x.gbl);
	int v0,sind,i = 0;
	do {
		sind = base.seg(i);
		v0 = x.seg(sind).pnt(0);
		/* Rotate residual for better diagonal dominance */
		r_gbl->res(v0)(0) = 0.0;
		r_gbl->res(v0)(1) = 0.0;
	} while (++i < base.nseg);
	v0 = x.seg(sind).pnt(1);
	r_gbl->res(v0)(0) = 0.0;
	r_gbl->res(v0)(1) = 0.0;
#endif	
		
	return;
}

void melt::update(int stage) {
	int i,m,n,msgn,count,sind,v0;
	
	if (stage < 0) return;
	
	base.comm_prepare(boundary::all,0,boundary::master_slave);
	base.comm_exchange(boundary::all,0,boundary::master_slave);
	base.comm_wait(boundary::all,0,boundary::master_slave);
	
	count = 0;
	i = base.nseg-1;
	do {
		sind = base.seg(i);
		v0 = x.seg(sind).pnt(1);
		for(n=0;n<tri_mesh::ND;++n) 
			x.pnts(v0)(n) = base.frcvbuf(0,count++);
	} while (--i >= 0);
	v0 = x.seg(sind).pnt(0);
	for(n=0;n<tri_mesh::ND;++n)
		x.pnts(v0)(n) = base.frcvbuf(0,count++);
	
	if (basis::tri(x.log2p)->sm() > 0) {
		for(i=base.nseg-1;i>=0;--i) {
			sind = base.seg(i);
			msgn = 1;
			for(m=0;m<basis::tri(x.log2p)->sm();++m) {
				for(n=0;n<tri_mesh::ND;++n)
					crds(i,m,n) = msgn*base.frcvbuf(0,count++);
				msgn *= -1;
			}
		}
	}
	return;
}

#ifdef petsc
void melt::petsc_jacobian() {
	int sm = basis::tri(x.log2p)->sm();	
	int vdofs;
	if (x.mmovement != x.coupled_deformable)
		vdofs = x.NV;
	else
		vdofs = x.NV+x.ND;
	Array<FLT,2> K(vdofs*(sm+2),vdofs*(sm+2));
	Array<FLT,1> row_store(vdofs*(sm+2));
	Array<int,1> loc_to_glo(vdofs*(sm+2));
	
	/* ZERO ROWS CREATED BY R_MESH */
	Array<int,1> indices((base.nseg+1)*tri_mesh::ND);
	int j = 0;
	int cnt = 0;
	int sind;
	do {
		sind = base.seg(j);
		indices(cnt++) = x.seg(sind).pnt(0)*vdofs+x.NV;
		indices(cnt++) = x.seg(sind).pnt(0)*vdofs+x.NV+1;
	} while (++j < base.nseg);
	indices(cnt++) = x.seg(sind).pnt(1)*vdofs+x.NV;
	indices(cnt++) = x.seg(sind).pnt(1)*vdofs+x.NV+1;
	
#ifdef MY_SPARSE
	x.J.zero_rows(cnt,indices);
	x.J_mpi.zero_rows(cnt,indices);
#else
	/* Must zero rows of jacobian created by r_mesh */
	MatAssemblyBegin(x.petsc_J,MAT_FINAL_ASSEMBLY); 
	MatAssemblyEnd(x.petsc_J,MAT_FINAL_ASSEMBLY); 
	MatZeroRows(x.petsc_J,cnt,indices.data(),PETSC_NULL);
#endif
	
	if (sm) {
		for (int j=0;j<base.nseg;++j) {
			int sind = base.seg(j);
			
			/* Now fill in effect of curvature on element resdiuals */
			Array<TinyVector<FLT,MXTM>,1> R(x.NV),Rbar(x.NV),lf_re(x.NV),lf_im(x.NV);
#ifdef BZ_DEBUG
			const FLT eps_r = 0.0e-6, eps_a = 1.0e-6;  /*<< constants for debugging jacobians */
#else
			const FLT eps_r = 1.0e-6, eps_a = 1.0e-10;  /*<< constants for accurate numerical determination of jacobians */
#endif
			
			int tind = x.seg(sind).tri(0);		
			x.ugtouht(tind);
			FLT dx = eps_r*x.distance(x.seg(sind).pnt(0),x.seg(sind).pnt(1)) +eps_a;
			const int tm = basis::tri(x.log2p)->tm();
			const int im = basis::tri(x.log2p)->im();
			Array<FLT,2> Ke(x.NV*tm,x.ND*sm);
			Array<int,1> loc_to_glo_e(x.NV*tm);
			Array<int,1> loc_to_glo_crv(sm*tri_mesh::ND);
			
			x.element_rsdl(tind,0,x.uht,lf_re,lf_im);
			for(int i=0;i<tm;++i) 
				for(int n=0;n<x.NV;++n) 
					Rbar(n)(i)=lf_re(n)(i)+lf_im(n)(i);
			
			int kcol = 0;
			for(int mode = 2; mode < sm+2; ++mode) {
				for(int var = 0; var < tri_mesh::ND; ++var) {
					crds(j,mode-2,var) += dx;
					
					x.element_rsdl(tind,0,x.uht,lf_re,lf_im);
					
					int krow = 0;
					for(int i=0;i<tm;++i)
						for(int n=0;n<x.NV;++n) 
							Ke(krow++,kcol) = (lf_re(n)(i) +lf_im(n)(i) -Rbar(n)(i))/dx;	
					
					
					++kcol;
					crds(j,mode-2,var) -= dx;
				}
			}
			
			int ind = 0;
			for (int m = 0; m < 3; ++m) {
				int gindx = vdofs*x.tri(tind).pnt(m);
				for (int n = 0; n < x.NV; ++n)
					loc_to_glo_e(ind++) = gindx++;
			}		
			
			/* EDGE MODES */
			if (sm) {
				for(int i = 0; i < 3; ++i) {
					int gindx = x.npnt*vdofs +x.tri(tind).seg(i)*sm*x.NV;
					int sgn = x.tri(tind).sgn(i);
					int msgn = 1;
					for (int m = 0; m < sm; ++m) {
						for(int n = 0; n < x.NV; ++n) {
							for(int j = 0; j < sm*x.ND; ++j) {
								Ke(ind,j) *= msgn;
							}
							loc_to_glo_e(ind++) = gindx++;
						}
						msgn *= sgn;
					}
				}
			}
			
			/* INTERIOR	MODES */
			if (tm) {
				int gindx = x.npnt*vdofs +x.nseg*sm*x.NV +tind*im*x.NV;
				for(int m = 0; m < im; ++m) {
					for(int n = 0; n < x.NV; ++n){
						loc_to_glo_e(ind++) = gindx++;
					}
				}
			}
			
			int gindxND = jacobian_start +j*tri_mesh::ND*sm;
			for(int m=0;m<sm*tri_mesh::ND;++m)
				loc_to_glo_crv(m) = gindxND++;
			
#ifdef MY_SPARSE
			x.J.add_values(tm*x.NV,loc_to_glo_e,sm*tri_mesh::ND,loc_to_glo_crv,Ke);
#else
			MatSetValuesLocal(x.petsc_J,tm*x.NV,loc_to_glo_e.data(),sm*tri_mesh::ND,loc_to_glo_crv.data(),Ke.data(),ADD_VALUES);
#endif		
		}
	}
}

void melt::petsc_matchjacobian_snd() {	
	int vdofs;
	if (x.mmovement != x.coupled_deformable)
		vdofs = x.NV;
	else
		vdofs = x.NV+x.ND;
	
	if (!base.is_comm()) return;
	
	/* Now do stuff for communication boundaries */
	int row,sind=-2;
	
	std::vector<int> c0vars;
	c0vars.push_back(0);  // Only temperature row
	 
	/* I am cheating here and sending floats and int's together */
#ifdef MY_SPARSE
	/* Send Jacobian entries for u,v but not p */
	base.sndsize() = 0;
	base.sndtype() = boundary::flt_msg;
	base.fsndbuf(base.sndsize()++) = x.jacobian_start +FLT_EPSILON;
	
	for(int i=0;i<base.nseg;++i) {
		sind = base.seg(i);
		int rowbase = x.seg(sind).pnt(0)*vdofs; 
		
		/* attach diagonal column # to allow continuity enforcement */
		base.fsndbuf(base.sndsize()++) = rowbase +FLT_EPSILON;;
		
		for(std::vector<int>::iterator it=c0vars.begin();it!=c0vars.end();++it) {
			row = rowbase + *it;
			base.fsndbuf(base.sndsize()++) = x.J._cpt(row+1) -x.J._cpt(row) +FLT_EPSILON;
#ifdef MPDEBUG
			*x.gbl->log << "sending " << x.J._cpt(row+1) -x.J._cpt(row) << " jacobian entries for vertex " << row/vdofs << " and variable " << *it << std::endl;
#endif
			for (int col=x.J._cpt(row);col<x.J._cpt(row+1);++col) {
#ifdef MPDEBUG
				*x.gbl->log << x.J._col(col) << ' ';
#endif
				base.fsndbuf(base.sndsize()++) = x.J._col(col) +FLT_EPSILON;
				base.fsndbuf(base.sndsize()++) = x.J._val(col);
			}
#ifdef MPDEBUG
			*x.gbl->log << std::endl;
#endif
		}
		
		/* Send Side Information */
		row = x.npnt*vdofs +sind*x.NV*x.sm0;
		
		/* attach diagonal column # to allow continuity enforcement */
		base.fsndbuf(base.sndsize()++) = row +FLT_EPSILON;
		
		for(int mode=0;mode<x.sm0;++mode) {
			base.fsndbuf(base.sndsize()++) = x.J._cpt(row+1) -x.J._cpt(row) +FLT_EPSILON;
#ifdef MPDEBUG
			*x.gbl->log << "sending " << x.J._cpt(row+1) -x.J._cpt(row) << " jacobian entries for side " << sind << " and variable " << 0 << std::endl;
#endif
			for (int col=x.J._cpt(row);col<x.J._cpt(row+1);++col) {
#ifdef MPDEBUG
				*x.gbl->log << x.J._col(col) << ' ';
#endif
				base.fsndbuf(base.sndsize()++) = x.J._col(col) +FLT_EPSILON;
				base.fsndbuf(base.sndsize()++) = x.J._val(col);
			}
			
#ifdef MPDEBUG
			*x.gbl->log << std::endl;
#endif
			++row;
		}
	}
	
	/* LAST POINT */
	int rowbase = x.seg(sind).pnt(1)*vdofs; 
	
	/* attach diagonal # to allow continuity enforcement */
	base.fsndbuf(base.sndsize()++) = rowbase +FLT_EPSILON;
	
	for(std::vector<int>::iterator it=c0vars.begin();it!=c0vars.end();++it) {
		row = rowbase + *it;
		base.fsndbuf(base.sndsize()++) = x.J._cpt(row+1) -x.J._cpt(row) +FLT_EPSILON;
#ifdef MPDEBUG
		*x.gbl->log << "sending " << x.J._cpt(row+1) -x.J._cpt(row) << " jacobian entries for vertex " << row/vdofs << " and variable " << *it << std::endl;
#endif
		for (int col=x.J._cpt(row);col<x.J._cpt(row+1);++col) {
#ifdef MPDEBUG
			*x.gbl->log << x.J._col(col) << ' ';
#endif
			base.fsndbuf(base.sndsize()++) = x.J._col(col) +FLT_EPSILON;
			base.fsndbuf(base.sndsize()++) = x.J._val(col);
		}
		
#ifdef MPDEBUG
		*x.gbl->log << std::endl;
#endif
	}
	/* Send index of start of curved modes */
	base.fsndbuf(base.sndsize()++) = jacobian_start;
}

void melt::petsc_matchjacobian_rcv(int phase) {
	
	if (!base.is_comm() || base.matchphase(boundary::all_phased,0) != phase) return;
	
	Array<int,1> indices(tri_mesh::ND+1 +x.sm0);
	int count = 0;
	int Jstart_mpi = static_cast<int>(base.frcvbuf(0, count++));
	
	sparse_row_major *pJ_mpi;
	if (base.is_local(0)) {
		pJ_mpi = &x.J;
		Jstart_mpi = 0;
	}
	else {
		pJ_mpi = &x.J_mpi;
	}
	
	int vdofs;
	if (x.mmovement != x.coupled_deformable)
		vdofs = x.NV;
	else
		vdofs = x.NV+x.ND;
	
	int sm = basis::tri(x.log2p)->sm();	
	
	
	/* Now do stuff for communication boundaries */
	/* Now Receive Information */		
	for (int i=base.nseg-1;i>=0;--i) {
		int sind = base.seg(i);
		int rowbase = x.seg(sind).pnt(1)*vdofs; 
		int row_mpi = static_cast<int>(base.frcvbuf(0,count++)) +Jstart_mpi;
#ifdef MPDEBUG
		*x.gbl->log << "received vertex rowbase " << row_mpi << std::endl;
#endif	
		
#ifdef MY_SPARSE
		/* ZERO ROWS AND SET EQUATIONS TO BE EQUALITY CONSTRAINTS (FOR T,X,Y FOR EACH VERTEX) */
		indices(0) = rowbase;
		indices(1) = rowbase+1;
		indices(2) = rowbase+2;
		x.J.zero_rows(3,indices);
		(*pJ_mpi).zero_rows(3,indices);
		x.J.set_diag(3,indices,1.0);
		(*pJ_mpi).set_values(rowbase,row_mpi+2,-1.0);
		(*pJ_mpi).set_values(rowbase+1,row_mpi+4,-1.0);
		(*pJ_mpi).set_values(rowbase+2,row_mpi+5,-1.0);
#else
		This is not working
#endif

		/* Now receive side Jacobian information */
		int row = x.npnt*vdofs +sind*x.NV*x.sm0;
		row_mpi = static_cast<int>(base.frcvbuf(0,count++)) +Jstart_mpi;
#ifdef MPDEBUG
		*x.gbl->log << "received side rowbase " << row_mpi << std::endl;
#endif	
		if (x.sm0) {
			for(int mode=0;mode<x.sm0;++mode)
				indices(mode) = row +mode;
				
			x.J.zero_rows(x.sm0,indices);
			(*pJ_mpi).zero_rows(x.sm0,indices);
			x.J.set_diag(x.sm0,indices,1.0);
			FLT sgn = -1;
			for(int mode=0;mode<x.sm0;++mode) {
				/* In buoyancy class, there are 4 degrees of freedom and T is 2nd */
				(*pJ_mpi).add_values(row+mode,row_mpi +4*mode +2,sgn);
				sgn *= -1.;
			}
		}
	}

	int sind = base.seg(0);
	int rowbase = x.seg(sind).pnt(0)*vdofs; 
	int row_mpi = static_cast<int>(base.frcvbuf(0,count++)) +Jstart_mpi;
#ifdef MPDEBUG
	*x.gbl->log << "received vertex rowbase " << row_mpi << std::endl;
#endif	
#ifdef MY_SPARSE
	indices(0) = rowbase;
	indices(1) = rowbase+1;
	indices(2) = rowbase+2;
	x.J.zero_rows(3,indices);
	(*pJ_mpi).zero_rows(3,indices);
	x.J.set_diag(3,indices,1.0);
	(*pJ_mpi).add_values(rowbase,row_mpi+2,-1.0);
	(*pJ_mpi).add_values(rowbase+1,row_mpi+4,-1.0);
	(*pJ_mpi).add_values(rowbase+2,row_mpi+5,-1.0);
#else
	This is not working
#endif
#endif
	
	if (sm && !base.is_frst()) {
		/* Equality of curved mode constraint */
		int ind = jacobian_start;
		int ind_mpi = static_cast<int>(base.frcvbuf(0,count++)) +(base.nseg-1)*sm*x.ND +Jstart_mpi;
#ifdef MPDEBUG
		*x.gbl->log << "received jacobian start " << ind_mpi << std::endl;
#endif	
		for(int i=0;i<base.nseg;++i) {
			int sgn = 1;
			for(int mode=0;mode<x.sm0;++mode) {
				for (int n=0;n<x.ND;++n) {
					x.J(ind,ind) = 1.0;
					(*pJ_mpi)(ind,ind_mpi) = -1.0*sgn;
					++ind;
					++ind_mpi;
				}
				sgn *= -1;
			}
			ind_mpi -= 2*sm*x.ND;		
		}
	}
}


int melt::petsc_make_1D_rsdl_vector(Array<double,1> res) {
	int sm = basis::tri(x.log2p)->sm();
	int ind = 0;
	
	/* T,x,y continuity constraint */
	
	
	/* Curvature equality constraint */
	for(int j=0;j<base.nseg;++j) {
		for(int m=0;m<sm;++m) {
			res(ind++) = 0.0;
			res(ind++) = 0.0;
		}
	}
	return(ind);
}


void melt::non_sparse(Array<int,1> &nnzero) {
	const int sm=basis::tri(x.log2p)->sm();
	const int im=basis::tri(x.log2p)->im();
	const int NV = x.NV;
	const int ND = tri_mesh::ND;
	
	int vdofs;
	if (x.mmovement != tri_hp::coupled_deformable) 
		vdofs = NV;
	else
		vdofs = ND+NV;
	
	int begin_seg = x.npnt*vdofs;
	int begin_tri = begin_seg+x.nseg*sm*NV;
	
	if(x.sm0 > 0) {
		if (base.is_frst()) {
			*x.gbl->log << "Solid melt boundary should not be first boundary condition" << std::endl;
			sim::abort(__LINE__,__FILE__,x.gbl->log);
			//nnzero(Range(jacobian_start,jacobian_start+base.nseg*sm*tri_mesh::ND-1)) = vdofs*(sm+2);
		}
		else {
			/* Just an equality constraint */
			nnzero(Range(jacobian_start,jacobian_start+base.nseg*sm*tri_mesh::ND-1)) = 1;
		}
		
		for (int i=0;i<base.nseg;++i) {
			int sind = base.seg(i);
			int tind = x.seg(sind).tri(0);
			if (im) {
				nnzero(Range(begin_tri +tind*im*NV,begin_tri +(tind+1)*im*NV-1)) += ND*sm;
			}
			
			for(int j=0;j<3;++j) {
				int sind1 = x.tri(tind).seg(j);
				nnzero(Range(begin_seg+sind1*NV*sm,begin_seg+(sind1+1)*NV*sm-1)) += ND*sm;
				if (sind1 == sind) {
					/* For opposing vertex should only be ND*sm for flow */
					/* mesh deformation equation doesn't depend on curvatures */
					int pind = x.tri(tind).pnt(j);
					nnzero(Range(pind*vdofs,pind*vdofs +x.NV -1)) += ND*sm;
				}
			}
			
			int pind = x.seg(sind).pnt(0);
			nnzero(Range(pind*vdofs,pind*vdofs +x.NV-1)) += ND*sm;
			
			pind = x.seg(sind).pnt(1);
			nnzero(Range(pind*vdofs,pind*vdofs +x.NV-1)) += ND*sm;
		}
	}
}


void melt::non_sparse_snd(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi) {
	
	if (!base.is_comm()) return;
	
	const int sm=basis::tri(x.log2p)->sm();
	const int NV = x.NV;
	const int ND = tri_mesh::ND;
	
	int vdofs;
	if (x.mmovement != tri_hp::coupled_deformable) 
		vdofs = NV;
	else
		vdofs = ND+NV;
	
	int begin_seg = x.npnt*vdofs;	
	
	std::vector<int> c0vars;
	c0vars.push_back(0);
	for(int n=x.NV;n<vdofs;++n) {
		c0vars.push_back(n);
	}		
	
	/* Going to send all jacobian entries,  Diagonal entries for matching DOF's will be merged together not individual */
	/* Send number of non-zeros to matches */
	base.sndsize() = 0;
	base.sndtype() = boundary::int_msg;
	for (int i=0;i<base.nseg;++i) {
		int sind = base.seg(i);
		int pind = x.seg(sind).pnt(0)*vdofs;
		for(std::vector<int>::iterator it=c0vars.begin();it!=c0vars.end();++it)
			base.isndbuf(base.sndsize()++) = nnzero(pind +*it);
	}
	int sind = base.seg(base.nseg-1);
	int pind = x.seg(sind).pnt(1)*vdofs;
	for(std::vector<int>::iterator it=c0vars.begin();it!=c0vars.end();++it)
		base.isndbuf(base.sndsize()++) = nnzero(pind +*it);
	
	/* Last thing to send is nnzero for edges (all the same) */
	if (sm)
		base.isndbuf(base.sndsize()++) = nnzero(begin_seg+sind*NV*sm);
	
	return;
}

void melt::non_sparse_rcv(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi) {
	
	if (!base.is_comm()) return;
	
	const int sm=basis::tri(x.log2p)->sm();
	const int NV = x.NV;
	const int ND = tri_mesh::ND;
	
	if (sm) nnzero_mpi(Range(jacobian_start,jacobian_start+base.nseg*sm*tri_mesh::ND-1)) = 1;
	
	int vdofs;
	if (x.mmovement != tri_hp::coupled_deformable) 
		vdofs = NV;
	else
		vdofs = ND+NV;
	
	int begin_seg = x.npnt*vdofs;
	
	std::vector<int> c0vars;
	c0vars.push_back(0);
	for(int n=x.NV;n<vdofs;++n) {
		c0vars.push_back(n);
	}		
	
	int count = 0;
	for (int i=base.nseg-1;i>=0;--i) {
		int sind = base.seg(i);
		int pind = x.seg(sind).pnt(1)*vdofs;
		for(std::vector<int>::iterator it=c0vars.begin();it!=c0vars.end();++it) {
			//nnzero_mpi(pind+*it) += base.ircvbuf(0,count++);
			nnzero_mpi(pind+*it) += 1;  // Just continuity constraint
			count++;
		}
	}
	int sind = base.seg(0);
	int pind = x.seg(sind).pnt(0)*vdofs;
	for(std::vector<int>::iterator it=c0vars.begin();it!=c0vars.end();++it) {
		//nnzero_mpi(pind+*it) += base.ircvbuf(0,count++);
		nnzero_mpi(pind+*it) += 1; // Just continuity constraint
		count++;
	}
	
	/* Now add to side degrees of freedom */
	if (sm) {
		// int toadd = base.ircvbuf(0,count++); 
		for (int i=0;i<base.nseg;++i) {
			int sind = base.seg(i);
			for (int mode=0;mode<sm;++mode) {
				for(int n=0;n<x.NV;++n) {
					//nnzero_mpi(begin_seg+sind*NV*sm +mode*NV +n) += toadd;
					nnzero_mpi(begin_seg+sind*NV*sm +mode*NV +n) += 1;  // Just continuity constraint

				}
			}
		}
	}
}
#endif

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

void kellerman::init(input_map& inmap,void* gbl_in) {
	melt::init(inmap,gbl_in);
	
	/* Set u,v fixed, T & P to be free */
	essential_indices.clear();
	type = natural;
	
	return;
}




