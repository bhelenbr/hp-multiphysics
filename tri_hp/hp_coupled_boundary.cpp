//
//  hp_coupled_bdry.cpp
//  tri_hp
//
//  Created by Brian Helenbrook on 12/27/14.
//
//

#include "hp_coupled_boundary.h"

//#define DEBUG

void hp_coupled_bdry::init(input_map& inmap,void* gbl_in) {
	std::string keyword,val;
	std::istringstream data;
	std::string filename;

	if (!inmap.get(base.idprefix + "_NV",NV)) {
		NV = 2; // Default is two variables (x,y)
	}
	
	keyword = base.idprefix + "_coupled";
	inmap[keyword] = "1";
	
	keyword = base.idprefix + "_curved";
	inmap[keyword] = "1";
	
	hp_edge_bdry::init(inmap,gbl_in);
	
	/* Stuff for a surface variable */
//	ugbd.resize(x.gbl->nhist+1);
//	for(int i=0;i<x.gbl->nhist+1;++i) {
//		ugbd(i).v.resize(base.maxseg+1,NV);
//		ugbd(i).s.resize(base.maxseg,x.sm0,NV);
//	}
//	ug.v.reference(ugbd(0).v);
//	ug.s.reference(ugbd(1).s);
//	
//	
//	/* UNSTEADY SOURCE TERMS */
//	dugdt.resize(x.log2pmax+1);
//#ifdef petsc
//	int start = x.log2pmax;
//#else
//	int start = 0;
//#endif
//	for(int i=start;i<=x.log2pmax;++i) {
//		dugdt(i).resize(base.maxseg,NV,basis::tri(i)->gpx());
//	}
	
	gbl = static_cast<global *>(gbl_in);
	
	std::string side_id, matching_block, matching_boundary;
	find_matching_boundary_name(inmap, matching_block, side_id);
	matching_boundary = matching_block +"_" +side_id;

	// These only matter for petsc
	if (!inmap.get(base.idprefix + "_one_sided",gbl->one_sided)) {
		if (!inmap.get(matching_boundary + "_one_sided",gbl->one_sided)) {
			gbl->one_sided = false;
		}
	}
	else if (inmap.get(matching_boundary + "_one_sided",gbl->one_sided)) {
		*x.gbl->log << "only specify one_sided for one or the other side.  doesn't matter which" << std::endl;
		sim::abort(__LINE__,__FILE__,x.gbl->log);
	}
	
	if (!inmap.get(base.idprefix + "_symmetric",gbl->symmetric)) {
		if (!inmap.get(matching_boundary + "_symmetric",gbl->symmetric)) {
			gbl->symmetric = false;
		}
	}
	else if (inmap.get(matching_boundary + "_symmetric",gbl->symmetric)) {
		*x.gbl->log << "only specify symmetric for one or the other side.  doesn't matter which" << std::endl;
		sim::abort(__LINE__,__FILE__,x.gbl->log);
	}
	
	if (!inmap.get(base.idprefix + "_precondition",gbl->precondition)) {
		if (!inmap.get(matching_boundary + "_precondition",gbl->precondition)) {
			gbl->precondition = false;
		}
	}
	else if (inmap.get(matching_boundary + "_precondition",gbl->precondition)) {
		*x.gbl->log << "only specify precondition for one or the other side.  It doesn't matter which" << std::endl;
		sim::abort(__LINE__,__FILE__,x.gbl->log);
	}

	if (gbl->precondition && (!gbl->one_sided && !gbl->symmetric)) {
		*x.gbl->log << "precondition can only be used with symmetric or one_sided" << std::endl;
		sim::abort(__LINE__,__FILE__,x.gbl->log);
	}
	
	int NVcoupled;
	inmap.getwdefault(base.idprefix + "_field_is_coupled",gbl->field_is_coupled,false);
	if (gbl->field_is_coupled) {
		NVcoupled = NV +c0_indices.size();
	}
	else {
		NVcoupled = NV; // Default is two variables (x,y)
	}

	gbl->vdt.resize(base.maxseg+1,NVcoupled,NVcoupled);
	gbl->vpiv.resize(base.maxseg+1,NVcoupled);
	gbl->sdt.resize(base.maxseg,NVcoupled,NVcoupled);
	gbl->spiv.resize(base.maxseg,NVcoupled);
	gbl->sdt2.resize(base.maxseg,x.sm0,NVcoupled,NVcoupled);
	gbl->spiv2.resize(base.maxseg,x.sm0,NVcoupled);
	gbl->meshc.resize(base.maxseg,NV);
	
	gbl->cfl.resize(x.log2p+1,NV);
	for (int n=0;n<NV;++n) {
		stringstream nstr;
		nstr << n;
		inmap.getlinewdefault(base.idprefix+"_cfl"+nstr.str(),val,std::string("2.5 1.5 1.0"));
		data.str(val);
		for (int m=0;m<x.log2p+1;++m) {
				data >> gbl->cfl(m,n);
		}
		data.clear();
	}
	
	if (x.seg(base.seg(0)).pnt(0) == x.seg(base.seg(base.nseg-1)).pnt(1)) gbl->is_loop = true;
	else gbl->is_loop = false;
	
	if (!gbl->symmetric && !is_master) return;   /* This is all that is necessary for a slave boundary */
	
	ksprg.resize(base.maxseg);

	gbl->vug0.resize(base.maxseg+1,NV);
	gbl->sug0.resize(base.maxseg,x.sm0,NV);
	
	gbl->vres0.resize(base.maxseg+1,NV);
	gbl->sres0.resize(base.maxseg,x.sm0,NV);
		
#ifdef DETAILED_MINV
	const int sm0 = x.sm0;
	gbl->ms.resize(base.maxseg,NV*sm0,NV*sm0);
	gbl->vms.resize(base.maxseg,NV,2,sm0,2);
	gbl->ipiv.resize(base.maxseg,NV*sm0);
#endif
    
	/* These are used by the slave to store transmitted preconditioner and residuals */
	/* Mostly assume that this is one way transmission for now */
	gbl->vres.resize(base.maxseg+1,NV);
	gbl->sres.resize(base.maxseg,x.sm0,NV);
		
	/* Multigrid Storage all except highest order (log2p+1)*/
	vdres.resize(x.log2p+1,base.maxseg+1,NV);
	sdres.resize(x.log2p+1,base.maxseg,x.sm0,NV);
	
	gbl->fadd.resize(NV);
	gbl->fadd = 1.0; // default multiplier is 1.0
	keyword = base.idprefix + "_fadd";
	if (inmap.getline(keyword,val)) {
		data.str(val);
		for (int n=0;n<NV;++n) {
			data >> gbl->fadd(n);
		}
		data.clear();
	}
	
	inmap.getwdefault(base.idprefix +"_adis",gbl->adis,1.0);
	
	return;
}

void hp_coupled_bdry::tadvance() {
	int i,j,m,sind;
	
	hp_edge_bdry::tadvance();
	
	/* Fixme: Need to fill in stuff for extra variables here */
	
	if (x.gbl->substep == 0) {
		if (gbl->symmetric || is_master) {
			/* SET SPRING CONSTANTS */
			for(j=0;j<base.nseg;++j) {
				sind = base.seg(j);
				ksprg(j) = 1.0/x.distance(x.seg(sind).pnt(0),x.seg(sind).pnt(1));
			}
		}
		
		if (!gbl->symmetric) {
			if (is_master) {
				/* CALCULATE TANGENT SOURCE TERM FOR FINE MESH */
				/* ZERO TANGENTIAL MESH MOVEMENT SOURCE */
				if (!x.coarse_level) {
					for(i=0;i<base.nseg+1;++i)
						vdres(x.log2p,i,0) = 0.0;
					
					for(i=0;i<base.nseg;++i)
						for(m=0;m<basis::tri(x.log2p)->sm();++m)
							sdres(x.log2p,i,m,0) = 0.0;
					
					rsdl(x.gbl->nstage);
					
					for(i=0;i<base.nseg+1;++i)
						vdres(x.log2p,i,0) = -gbl->vres(i,0);
					
					for(i=0;i<base.nseg;++i)
						for(m=0;m<basis::tri(x.log2p)->sm();++m)
							sdres(x.log2p,i,m,0) = -gbl->sres(i,m,0)*0;  /* TO KEEP SIDE MODES EQUALLY SPACED */
				}
			}
			else {
				if (!x.coarse_level) {
					rsdl(x.gbl->nstage);  // Matching call by slave for communication
				}
			}
		}
		else {
			/* CALCULATE TANGENT SOURCE TERM FOR FINE MESH */
			/* ZERO TANGENTIAL MESH MOVEMENT SOURCE */
			if (!x.coarse_level) {
				for(i=0;i<base.nseg+1;++i)
					vdres(x.log2p,i,0) = 0.0;
				
				for(i=0;i<base.nseg;++i)
					for(m=0;m<basis::tri(x.log2p)->sm();++m)
						sdres(x.log2p,i,m,0) = 0.0;
				
				rsdl(x.gbl->nstage);
				
				for(i=0;i<base.nseg+1;++i)
					vdres(x.log2p,i,0) = -gbl->vres(i,0);
				
				for(i=0;i<base.nseg;++i)
					for(m=0;m<basis::tri(x.log2p)->sm();++m)
						sdres(x.log2p,i,m,0) = -gbl->sres(i,m,0)*0;  /* TO KEEP SIDE MODES EQUALLY SPACED */
			}
		}
	}
	return;
}

void hp_coupled_bdry::rsdl(int stage) {
	
	if (!gbl->symmetric && !is_master) return;
	
	int m,n,sind,indx,v0,v1;
	TinyVector<FLT,tri_mesh::ND> norm, rp;
	Array<FLT,1> ubar(x.NV);
	Array<TinyVector<FLT,MXGP>,1> u(x.NV);
	TinyMatrix<FLT,tri_mesh::ND,MXGP> crd, dcrd, mvel;
	TinyMatrix<FLT,8,MXGP> res;
	Array<TinyVector<FLT,MXTM>,1> lf(x.NV+NV);
	
	/**************************************************/
	/* DETERMINE MESH RESIDUALS & SURFACE TENSION      */
	/**************************************************/
	for(n=0;n<NV;++n)
		gbl->vres(0,n) = 0.0;
	
	for(indx=0;indx<base.nseg;++indx) {
		sind = base.seg(indx);
		v0 = x.seg(sind).pnt(0);
		v1 = x.seg(sind).pnt(1);
		
		x.ugtouht1d(sind);
		element_rsdl(indx,lf);
		
		/* ADD FLUXES TO RESIDUAL */
		for(n=0;n<x.NV;++n)
			x.gbl->res.v(v0,n) += lf(n)(0);
		
		for(n=0;n<x.NV;++n)
			x.gbl->res.v(v1,n) += lf(n)(1);
		
		for(m=0;m<basis::tri(x.log2p)->sm();++m) {
			for(n=0;n<x.NV;++n)
				x.gbl->res.s(sind,m,n) += lf(n)(m+2);
		}
		
		/* STORE MESH-MOVEMENT RESIDUAL IN VRES/SRES */
		for(n=0;n<NV;++n) {
			gbl->vres(indx,n) += lf(x.NV+n)(0);
			gbl->vres(indx+1,n) = lf(x.NV+n)(1);
			for(m=0;m<basis::tri(x.log2p)->sm();++m)
				gbl->sres(indx,m,n) = lf(x.NV+n)(m+2);
		}
	}
	
	if (!x.coarse_flag) {
		/* ADD TANGENTIAL MESH MOVEMENT SOURCE */
		for(int i=0;i<base.nseg+1;++i)
			gbl->vres(i,0) += vdres(x.log2p,i,0);
		
		for(int i=0;i<base.nseg;++i)
			for(int m=0;m<basis::tri(x.log2p)->sm();++m)
				gbl->sres(i,m,0) += sdres(x.log2p,i,m,0);
	}
	
	
#ifdef petsc
	/* Store vertex mesh residual in r_mesh residual vector? */
	r_tri_mesh::global *r_gbl = dynamic_cast<r_tri_mesh::global *>(x.gbl);
	int i = 0;
	do {
		sind = base.seg(i);
		int v0 = x.seg(sind).pnt(0);
		r_gbl->res(v0)(0) = gbl->vres(i,0);
		r_gbl->res(v0)(1) = gbl->vres(i,1);
	} while (++i < base.nseg);
	v0 = x.seg(sind).pnt(1);
	r_gbl->res(v0)(0) = gbl->vres(i,0);
	r_gbl->res(v0)(1) = gbl->vres(i,1);
#endif
	
}

void hp_coupled_bdry::maxres() {
	if (!gbl->symmetric && !is_master) return;
	
	int i,n;
	TinyVector<FLT,tri_mesh::ND> mxr;
	
	mxr = 0.0;
	
	for(i=0;i<base.nseg+1;++i)
		for(n=0;n<NV;++n)
			mxr(n) = MAX(fabs(gbl->vres(i,n)),mxr(n));
	
	for(n=0;n<tri_mesh::ND;++n)
		*x.gbl->log << ' ' << mxr(n) << ' ';
	
	return;
}

void hp_coupled_bdry::update(int stage) {
	
	if (gbl->symmetric || is_master) {
		int i,sind,v0;
		
		if (stage < 0) {
			i = 0;
			do {
				sind = base.seg(i);
				v0 = x.seg(sind).pnt(0);
				for (int n=0;n<tri_mesh::ND;++n)
					gbl->vug0(i,n) = x.pnts(v0)(n);
			} while (++i < base.nseg);
			v0 = x.seg(sind).pnt(1);
			for (int n=0;n<tri_mesh::ND;++n)
				gbl->vug0(base.nseg,n) = x.pnts(v0)(n);
			
			if (basis::tri(x.log2p)->sm() > 0) {
				for(int i=0;i<base.nseg;++i)
					for(int m=0;m<basis::tri(x.log2p)->sm();++m)
						for(int n=0;n<tri_mesh::ND;++n)
							gbl->sug0(i,m,n) = crv(i,m)(n);
			}
			
			return;
		}
		
		minvrt();
		
#ifdef DEBUG
//		if (x.coarse_flag) {
		const int ncoupled = NV +gbl->field_is_coupled*c0_indices.size();

		for(int i=0;i<base.nseg+1;++i) {
			for(int row=0;row<ncoupled;++row) {
				*x.gbl->log << "vdt: " << i << ' ' << row << ' ';
				for(int col=0;col<ncoupled;++col) {
					*x.gbl->log << gbl->vdt(i,row,col) << ' ';
				}
				*x.gbl->log << std::endl;
			}
		}
		if (basis::tri(x.log2p)->sm() > 0) {
			for(int i=0;i<base.nseg;++i) {
				for(int row=0;row<ncoupled;++row) {
					*x.gbl->log << "sdt: " << i << ' ' << row << ' ';
					for(int col=0;col<ncoupled;++col) {
						*x.gbl->log << gbl->sdt(i,row,col) << ' ';
					}
					*x.gbl->log << std::endl;
				}
			}
		}
		
		for(int i=0;i<base.nseg+1;++i) {
			*x.gbl->log << "vres: " << i << ' ';
			for(int n=0;n<tri_mesh::ND;++n) {
				if (fabs(gbl->vres(i,n)) > 1.0e-9) *x.gbl->log << gbl->vres(i,n) << ' ';
				else *x.gbl->log << "0.0 ";
			}
			*x.gbl->log << '\n';
		}
		
		for(int i=0;i<base.nseg;++i) {
			for(int m=0;m<basis::tri(x.log2p)->sm();++m) {
				*x.gbl->log << "sres: " << i << ' ';
				for(int n=0;n<tri_mesh::ND;++n) {
					if (fabs(gbl->sres(i,m,n)) > 1.0e-9) *x.gbl->log << gbl->sres(i,m,n) << ' ';
					else *x.gbl->log << "0.0 ";
				}
				*x.gbl->log << '\n';
			}
		}
		
		i = 0;
		do {
			sind = base.seg(i);
			v0 = x.seg(sind).pnt(0);
			*x.gbl->log << "vertex positions " << v0 << ' ' << x.pnts(v0)(0) << ' ' << x.pnts(v0)(1) << '\n';
		} while(++i < base.nseg);
		v0 = x.seg(sind).pnt(1);
		*x.gbl->log << "vertex positions " << v0 << ' ' << x.pnts(v0)(0) << ' ' << x.pnts(v0)(1) << '\n';
		
		for(int i=0;i<base.nseg;++i)
			for(int m=0;m<basis::tri(x.log2p)->sm();++m)
				*x.gbl->log << "spos: " << i << ' ' << m << ' ' << crv(i,m)(0) << ' ' << crv(i,m)(1) << '\n';
//		}
#endif
		
		FLT alpha = x.gbl->alpha(stage);
        if (gbl->field_is_coupled)
            alpha *= x.gbl->cfl(x.log2p);
        
		i = 0;
		do {
			sind = base.seg(i);
			v0 = x.seg(sind).pnt(0);
			for (int n=0;n<NV;++n)
				x.pnts(v0)(n) = gbl->vug0(i,n) -alpha*gbl->vres(i,n);
		} while (++i < base.nseg);
		v0 = x.seg(sind).pnt(1);
		for (int n=0;n<NV;++n)
			x.pnts(v0)(n) = gbl->vug0(base.nseg,n) -alpha*gbl->vres(base.nseg,n);
		
		if (basis::tri(x.log2p)->sm() > 0) {
			for(int i=0;i<base.nseg;++i)
				for(int m=0;m<basis::tri(x.log2p)->sm();++m)
					for(int n=0;n<tri_mesh::ND;++n)
						crv(i,m)(n) = gbl->sug0(i,m,n) -alpha*gbl->sres(i,m,n);
		}
		
		if (gbl->field_is_coupled) {
			i = 0;
			do {
				sind = base.seg(i);
				v0 = x.seg(sind).pnt(0);
				for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
					x.ug.v(v0,*n) = x.gbl->ug0.v(v0,*n) -alpha*x.gbl->res.v(v0,*n);
				}
			} while (++i < base.nseg);
			v0 = x.seg(sind).pnt(1);
			for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
				x.ug.v(v0,*n) = x.gbl->ug0.v(v0,*n) -alpha*x.gbl->res.v(v0,*n);
			}
			
			if (basis::tri(x.log2p)->sm() > 0) {
				for(int i=0;i<base.nseg;++i) {
					sind = base.seg(i);
					for(int m=0;m<basis::tri(x.log2p)->sm();++m) {
						for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
							x.ug.s(sind,m,*n) = x.gbl->ug0.s(sind,m,*n) -alpha*x.gbl->res.s(sind,m,*n);
						}
					}
				}
			}
		}
		
		/* FIX POINTS THAT SLIDE ON CURVE */
		x.hp_vbdry(base.vbdry(1))->mvpttobdry(x.pnts(v0));
		sind = base.seg(0);
		v0 = x.seg(sind).pnt(0);
		x.hp_vbdry(base.vbdry(0))->mvpttobdry(x.pnts(v0));
		
#ifdef DEBUG
//		if (x.coarse_flag) {
		i = 0;
		do {
			sind = base.seg(i);
			v0 = x.seg(sind).pnt(0);
			*x.gbl->log << "vertex positions " << v0 << ' ' << x.pnts(v0)(0) << ' ' << x.pnts(v0)(1) << '\n';
		} while(++i < base.nseg);
		v0 = x.seg(sind).pnt(1);
		*x.gbl->log << "vertex positions " << v0 << ' ' << x.pnts(v0)(0) << ' ' << x.pnts(v0)(1) << '\n';
		
		for(int i=0;i<base.nseg;++i)
			for(int m=0;m<basis::tri(x.log2p)->sm();++m)
				*x.gbl->log << "spos: " << i << ' ' << m << ' ' << crv(i,m)(0) << ' ' << crv(i,m)(1) << '\n';
//		}
#endif
	}

	if (!gbl->symmetric) {
		if (is_master) {
			if (base.is_comm()) {
				int count = 0;
				int i = 0;
				int sind;
				do {
					sind = base.seg(i);
					int v0 = x.seg(sind).pnt(0);
					for(int n=0;n<tri_mesh::ND;++n)
						base.fsndbuf(count++) = x.pnts(v0)(n);
				} while(++i < base.nseg);
				int v0 = x.seg(sind).pnt(1);
				for(int n=0;n<tri_mesh::ND;++n)
					base.fsndbuf(count++) = x.pnts(v0)(n);
				
				for(int i=0;i<base.nseg;++i) {
					for(int m=0;m<basis::tri(x.log2p)->sm();++m) {
						for(int n=0;n<tri_mesh::ND;++n)
							base.fsndbuf(count++) = crv(i,m)(n);
					}
				}
				if (gbl->field_is_coupled) {
					i = 0;
					do {
						sind = base.seg(i);
						v0 = x.seg(sind).pnt(0);
						for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
							base.fsndbuf(count++) = x.ug.v(v0,*n);
						}
					} while (++i < base.nseg);
					v0 = x.seg(sind).pnt(1);
					for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
						base.fsndbuf(count++) = x.ug.v(v0,*n);
					}
					
					if (basis::tri(x.log2p)->sm() > 0) {
						for(int i=0;i<base.nseg;++i) {
							sind = base.seg(i);
							for(int m=0;m<basis::tri(x.log2p)->sm();++m) {
								for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
									base.fsndbuf(count++) = x.ug.s(sind,m,*n);
								}
							}
						}
					}
				}

				base.sndsize() = count;
				base.sndtype() = boundary::flt_msg;
				base.comm_prepare(boundary::all,0,boundary::master_slave);
				base.comm_exchange(boundary::all,0,boundary::master_slave);
				base.comm_wait(boundary::all,0,boundary::master_slave);
			}
		}
		else {
			int i,msgn,sind,v0;
			
			if (stage < 0) return;
			
			base.comm_prepare(boundary::all,0,boundary::master_slave);
			base.comm_exchange(boundary::all,0,boundary::master_slave);
			base.comm_wait(boundary::all,0,boundary::master_slave);
			
			int count = 0;
			i = base.nseg-1;
			do {
				sind = base.seg(i);
				v0 = x.seg(sind).pnt(1);
				for(int n=0;n<tri_mesh::ND;++n)
					x.pnts(v0)(n) = base.frcvbuf(0,count++);
			} while (--i >= 0);
			v0 = x.seg(sind).pnt(0);
			for(int n=0;n<tri_mesh::ND;++n)
				x.pnts(v0)(n) = base.frcvbuf(0,count++);
			
			if (basis::tri(x.log2p)->sm() > 0) {
				for(i=base.nseg-1;i>=0;--i) {
					msgn = 1;
					for(int m=0;m<basis::tri(x.log2p)->sm();++m) {
						for(int n=0;n<tri_mesh::ND;++n)
							crds(i,m,n) = msgn*base.frcvbuf(0,count++);
						msgn *= -1;
					}
				}
			}
			
			if (gbl->field_is_coupled) {
				i = base.nseg-1;
				do {
					sind = base.seg(i);
					v0 = x.seg(sind).pnt(1);
					for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n)
						x.ug.v(v0,*n) = base.frcvbuf(0,count++);
				} while (--i >= 0);
				v0 = x.seg(sind).pnt(0);
				for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n)
					x.ug.v(v0,*n) = base.frcvbuf(0,count++);
				
				if (basis::tri(x.log2p)->sm() > 0) {
					for(i=base.nseg-1;i>=0;--i) {
						sind = base.seg(i);
						msgn = 1;
						for(int m=0;m<basis::tri(x.log2p)->sm();++m) {
							for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n)
								x.ug.s(sind,m,*n) = msgn*base.frcvbuf(0,count++);
							msgn *= -1;
						}
					}
				}
			}
		}
	}
	
	return;
}

void hp_coupled_bdry::minvrt() {
	FLT temp;
	const int sm = basis::tri(x.log2p)->sm();
	
	/* INVERT MASS MATRIX */
	/* LOOP THROUGH SIDES */
	if (sm > 0) {
		for(int indx = 0; indx<base.nseg; ++indx) {
			/* SUBTRACT SIDE CONTRIBUTIONS TO VERTICES */
			for (int m=0; m < sm; ++m) {
				for(int n=0;n<NV;++n) {
					gbl->vres(indx,n) -= basis::tri(x.log2p)->sfmv1d(0,m)*gbl->sres(indx,m,n);
					gbl->vres(indx+1,n) -= basis::tri(x.log2p)->sfmv1d(1,m)*gbl->sres(indx,m,n);
				}
			}
		}
		// This is not necessary for flow vars because flow update will do this (sort of)
	}
	
	for(int last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		x.vbdry(base.vbdry(0))->vloadbuff(boundary::manifolds,&gbl->vres(0,0),0,1,0);
		x.vbdry(base.vbdry(1))->vloadbuff(boundary::manifolds,&gbl->vres(base.nseg,0),0,1,0);
		x.vbdry(base.vbdry(0))->comm_prepare(boundary::manifolds,mp_phase,boundary::symmetric);
		x.vbdry(base.vbdry(1))->comm_prepare(boundary::manifolds,mp_phase,boundary::symmetric);
		
		x.vbdry(base.vbdry(0))->comm_exchange(boundary::manifolds,mp_phase,boundary::symmetric);
		x.vbdry(base.vbdry(1))->comm_exchange(boundary::manifolds,mp_phase,boundary::symmetric);
		
		last_phase = true;
		
		last_phase &= x.vbdry(base.vbdry(0))->comm_wait(boundary::manifolds,mp_phase,boundary::symmetric);
		last_phase &= x.vbdry(base.vbdry(1))->comm_wait(boundary::manifolds,mp_phase,boundary::symmetric);
		
		x.vbdry(base.vbdry(0))->vfinalrcv(boundary::manifolds,mp_phase,boundary::symmetric,boundary::average,&gbl->vres(0,0),0,1,0);
		x.vbdry(base.vbdry(1))->vfinalrcv(boundary::manifolds,mp_phase,boundary::symmetric,boundary::average,&gbl->vres(base.nseg,0),0,1,0);
	}
	
	if (gbl->is_loop) {
		for(int n=0;n<NV;++n) {
			gbl->vres(0,n) = 0.5*(gbl->vres(0,n) +gbl->vres(base.nseg,n));
			gbl->vres(base.nseg,n) = gbl->vres(0,n);
		}
		gbl->vres(0,0) = 0.0;
		gbl->vres(base.nseg,0) = 0.0;
	}
	x.hp_vbdry(base.vbdry(0))->vdirichlet();
	x.hp_vbdry(base.vbdry(1))->vdirichlet();
	
	/* SOLVE FOR VERTEX MODES */
	if (gbl->field_is_coupled) {
		int info;
		char trans[] = "T";
		const int ncoupled = NV +gbl->field_is_coupled*c0_indices.size();
		Array<FLT,1> res_vec(ncoupled);

		for(int i=0;i<base.nseg;++i) {
			int sind = base.seg(i);
			int v0 = x.seg(sind).pnt(0);
			int indx = 0;
			for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
				res_vec(indx++) = x.gbl->res.v(v0,*n);
			}
			for(int n=0;n<NV;++n) {
				res_vec(indx++) = gbl->vres(i,n);
			}

#ifdef F2CFortran
			GETRS(trans,ncoupled,1,&gbl->vdt(i,0,0),gbl->vdt.length(secondDim),&gbl->vpiv(i,0),res_vec.data(),ncoupled,info);
#else
            const int one = 1;
            const int length = gbl->vdt.length(secondDim);
            dgetrs_(trans,&ncoupled,&one,&gbl->vdt(i,0,0),&length,&gbl->vpiv(i,0),res_vec.data(),&ncoupled,&info);
#endif
			if (info != 0) {
				*x.gbl->log << "DGETRS FAILED IN VERTEX PRECONDITIONER " << info << std::endl;
				sim::abort(__LINE__,__FILE__,x.gbl->log);
			}
			indx = 0;
			for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
				x.gbl->res.v(v0,*n) = res_vec(indx++);
			}
			for(int n=0;n<NV;++n) {
			 gbl->vres(i,n) = res_vec(indx++);
			}
		}
		int sind = base.seg(base.nseg-1);
		int v0 = x.seg(sind).pnt(1);
		int indx = 0;
		for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
			res_vec(indx++) = x.gbl->res.v(v0,*n);
		}
		for(int n=0;n<NV;++n) {
			res_vec(indx++) = gbl->vres(base.nseg,n);
		}
		
#ifdef F2CFortran
        GETRS(trans,ncoupled,1,&gbl->vdt(base.nseg,0,0),gbl->vdt.length(secondDim),&gbl->vpiv(base.nseg,0),res_vec.data(),ncoupled,info);
#else
        const int one = 1;
        const int length = gbl->vdt.length(secondDim);
        dgetrs_(trans,&ncoupled,&one,&gbl->vdt(base.nseg,0,0),&length,&gbl->vpiv(base.nseg,0),res_vec.data(),&ncoupled,&info);
#endif
		if (info != 0) {
			*x.gbl->log << "DGETRS FAILED IN VERTEX PRECONDITIONER " << info << std::endl;
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}
		
		indx = 0;
		for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
			x.gbl->res.v(v0,*n) = res_vec(indx++);
		}
		for(int n=0;n<NV;++n) {
			 gbl->vres(base.nseg,n) = res_vec(indx++);
		}
		
		if (sm > 0) {
			for(int j=0;j<base.nseg;++j) {
				const int sind = base.seg(j);
				
				for(int m=0;m<sm;++m) {
					int indx = 0;
					for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
						res_vec(indx++) = x.gbl->res.s(sind,m,*n);
					}
					for(int n=0;n<NV;++n) {
						res_vec(indx++) = gbl->sres(j,m,n);
					}
					
					int info;
					char trans[] = "T";
#ifdef F2CFortran
					GETRS(trans,ncoupled,1,&gbl->sdt(j,0,0),gbl->sdt.length(secondDim),&gbl->spiv(j,0),res_vec.data(),ncoupled,info);
#else
                    const int one = 1;
                    const int length = gbl->sdt.length(secondDim);
                    dgetrs_(trans,&ncoupled,&one,&gbl->sdt(j,0,0),&length,&gbl->spiv(j,0),res_vec.data(),&ncoupled,&info);
#endif
					if (info != 0) {
						*x.gbl->log << "DGETRS FAILED IN SIDE MODE PRECONDITIONER " << info << std::endl;
						sim::abort(__LINE__,__FILE__,x.gbl->log);
					}
					
					indx = 0;
					for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
						x.gbl->res.s(sind,m,*n) = res_vec(indx++);
					}
					for(int n=0;n<NV;++n) {
						gbl->sres(j,m,n) = res_vec(indx++);
					}
				}
			}
			
			// FIXME: THIS IS MESSED UP.  MINVRT FOR THE FLOW SCREWS THIS ALL UP
			for(int m=0;m<basis::tri(x.log2p)->sm();++m) {
				for(int n=0;n<NV;++n) {
					gbl->sres(indx,m,n) -= basis::tri(x.log2p)->vfms1d(0,m)*gbl->vres(indx,n);
					gbl->sres(indx,m,n) -= basis::tri(x.log2p)->vfms1d(1,m)*gbl->vres(indx+1,n);
				}
			}
		}
	}
	else {
		// FIXME: THIS ONLY WORKS FOR 2 DECOUPLED SIDE VARIABLES (X & Y)
		for(int i=0;i<base.nseg+1;++i) {
			temp                     = gbl->vres(i,0)*gbl->vdt(i,0,0) +gbl->vres(i,1)*gbl->vdt(i,0,1);
			gbl->vres(i,1) = gbl->vres(i,0)*gbl->vdt(i,1,0) +gbl->vres(i,1)*gbl->vdt(i,1,1);
			gbl->vres(i,0) = temp;
		}

		/* SOLVE FOR SIDE MODES */
		for(int indx = 0; indx<base.nseg; ++indx) {
			
#ifdef DETAILED_MINV
			for(m=0;m<basis::tri(x.log2p)->sm();++m) {
				for(n=0;n<NV;++n) {
					gbl->sres(indx,m,n) -= gbl->vms(indx,n,0,m,0)*gbl->vres(indx,0);
					gbl->sres(indx,m,n) -= gbl->vms(indx,n,0,m,1)*gbl->vres(indx+1,0);
					gbl->sres(indx,m,n) -= gbl->vms(indx,n,1,m,0)*gbl->vres(indx,1);
					gbl->sres(indx,m,n) -= gbl->vms(indx,n,1,m,1)*gbl->vres(indx+1,1);
				}
			}
			int info;
			char trans[] = "T";
			GETRS(trans,2*basis::tri(x.log2p)->sm(),1,&gbl->ms(indx,0,0),2*MAXP,&gbl->ipiv(indx,0),&gbl->sres(indx,0,0),2*MAXP,info);
			if (info != 0) {
				*x.gbl->log << "DGETRS FAILED FOR SIDE MODE UPDATE" << std::endl;
				sim::abort(__LINE__,__FILE__,x.gbl->log);
			}
#else
			/* INVERT SIDE MODES */
			DPBTRSNU2((double *) &basis::tri(x.log2p)->sdiag1d(0,0),basis::tri(x.log2p)->sbwth()+1,basis::tri(x.log2p)->sm(),basis::tri(x.log2p)->sbwth(),&(gbl->sres(indx,0,0)),NV);
			for(int m=0;m<basis::tri(x.log2p)->sm();++m) {
				temp                             = gbl->sres(indx,m,0)*gbl->sdt(indx,0,0) +gbl->sres(indx,m,1)*gbl->sdt(indx,0,1);
				gbl->sres(indx,m,1) = gbl->sres(indx,m,0)*gbl->sdt(indx,1,0) +gbl->sres(indx,m,1)*gbl->sdt(indx,1,1);
				gbl->sres(indx,m,0) = temp;
			}
			
			for(int m=0;m<basis::tri(x.log2p)->sm();++m) {
				for(int n=0;n<NV;++n) {
					gbl->sres(indx,m,n) -= basis::tri(x.log2p)->vfms1d(0,m)*gbl->vres(indx,n);
					gbl->sres(indx,m,n) -= basis::tri(x.log2p)->vfms1d(1,m)*gbl->vres(indx+1,n);
				}
			}
#endif
		}
	}
	
	return;
}

void hp_coupled_bdry::mg_restrict() {
	
	if (!gbl->symmetric && !is_master) return;

	int i,bnum,indx,tind,v0,snum,sind;
	
	if(x.p0 > 1) {
		/* TRANSFER IS ON FINEST MESH */
		gbl->vres0(Range(0,base.nseg),Range::all()) = gbl->vres(Range(0,base.nseg),Range::all());
		if (basis::tri(x.log2p)->sm() > 0) gbl->sres0(Range(0,base.nseg-1),Range(0,basis::tri(x.log2p)->sm()-1),Range::all()) = gbl->sres(Range(0,base.nseg-1),Range(0,basis::tri(x.log2p)->sm()-1),Range::all());
		return;
	}
	else {
		/* TRANSFER IS BETWEEN DIFFERENT MESHES */
		gbl->vres0(Range(0,base.nseg),Range::all()) = 0.0;
		
		/* CALCULATE COARSE RESIDUALS */
		/* DO ENDPOINTS FIRST */
		gbl->vres0(0,Range::all()) = gbl->vres(0,Range::all());
		gbl->vres0(base.nseg,Range::all()) = gbl->vres(fine->base.nseg,Range::all());
		
		tri_mesh *fmesh = dynamic_cast<tri_mesh *>(x.fine);
		for (bnum = 0; x.hp_ebdry(bnum) != this; ++bnum);
		for(i=1;i<fine->base.nseg;++i) {
			sind = fine->base.seg(i);
			v0 = fmesh->seg(sind).pnt(0);
			tind = fmesh->ccnnct(v0).tri;
			for(snum=0;snum<3;++snum)
				if (x.getbdrynum(x.tri(tind).tri(snum))  == bnum) break;
			assert(snum != 3);
			indx = x.getbdryseg(x.tri(tind).tri(snum));
			gbl->vres0(indx,Range::all()) += fmesh->ccnnct(v0).wt((snum+1)%3)*gbl->vres(i,Range::all());
			gbl->vres0(indx+1,Range::all()) += fmesh->ccnnct(v0).wt((snum+2)%3)*gbl->vres(i,Range::all());
		}
	}
	
	return;
}

void hp_coupled_bdry::mg_source() {
	/************************************************/
	/* MODIFY SURFACE RESIDUALS ON COARSER MESHES    */
	/************************************************/
	
	if(!gbl->symmetric && !is_master) return;

	if(x.coarse_flag) {
		if (x.isfrst) {
			for(int i=0;i<base.nseg+1;++i)
				for(int n=0;n<NV;++n)
					vdres(x.log2p,i,n) = gbl->fadd(n)*gbl->vres0(i,n) -gbl->vres(i,n);
			
			for(int i=0;i<base.nseg;++i)
				for(int m=0;m<basis::tri(x.log2p)->sm();++m)
					for(int n=0;n<NV;++n)
						sdres(x.log2p,i,m,n) = gbl->fadd(n)*gbl->sres0(i,m,n) -gbl->sres(i,m,n);
			
		}
		for(int i=0;i<base.nseg+1;++i)
			for(int n=0;n<NV;++n)
				gbl->vres(i,n) += vdres(x.log2p,i,n);
		
		for(int i=0;i<base.nseg;++i)
			for(int m=0;m<basis::tri(x.log2p)->sm();++m)
				for(int n=0;n<NV;++n)
					gbl->sres(i,m,n) += sdres(x.log2p,i,m,n);
	}
}


// This matches residual only
// Need to rewrite naming to make more general
void hp_coupled_bdry::smatchsolution_snd(FLT *sdata, int bgn, int end, int stride) {
	
	hp_edge_bdry::smatchsolution_snd(sdata, bgn, end, stride);
	if (!gbl->symmetric) return;
	
	const int sm = basis::tri(x.log2p)->sm();
	for(int i=0;i<base.nseg;++i) {
		for(int m=0;m<sm;++m) {
			for(int n=0;n<NV;++n) {
				base.fsndbuf(base.sndsize()++) = gbl->sres(i,m,n);
			}
		}
	}
	
	return;
}

int hp_coupled_bdry::smatchsolution_rcv(FLT *sdata, int bgn, int end, int stride) {
	int count = hp_edge_bdry::smatchsolution_rcv(sdata, bgn, end, stride);
	if (!gbl->symmetric) return(count);
	
	const int sm = basis::tri(x.log2p)->sm();

	for(int i=base.nseg-1;i>=0;--i) {
		int msgn = 1;
		for(int m=0;m<sm;++m) {
			for(int n=0;n<NV;++n) {
				gbl->sres(i,m,n) = 0.5*(gbl->sres(i,m,n) +msgn*base.frcvbuf(0,count++));
			}
			msgn *= -1;
		}
	}
	
	return(count);
}

void hp_coupled_bdry::element_jacobian(int indx, Array<FLT,2>& K) {
	int sm = basis::tri(x.log2p)->sm();

	Array<TinyVector<FLT,MXTM>,1> Rbar(x.NV+tri_mesh::ND),lf(x.NV+tri_mesh::ND);
	
	/* Calculate and store initial residual */
	int sind = base.seg(indx);
	x.ugtouht1d(sind);
	element_rsdl(indx,Rbar);
	
	Array<FLT,1> dw(x.NV);
	dw = 0.0;
	for(int i=0;i<2;++i)
		for(int n=0;n<x.NV;++n)
			dw(n) = dw(n) + fabs(x.uht(n)(i));
	
	dw *= eps_r;
	dw += eps_a;
	FLT dx = eps_r*x.distance(x.seg(sind).pnt(0),x.seg(sind).pnt(1)) +eps_a;
	
	/* Numerically create Jacobian */
	int kcol = 0;
	for(int mode = 0; mode < 2; ++mode) {
		for(int var = 0; var < x.NV; ++var) {
			x.uht(var)(mode) += dw(var);
			
			element_rsdl(indx,lf);
			
			int krow = 0;
			for(int k=0;k<sm+2;++k)
				for(int n=0;n<x.NV+NV;++n)
					K(krow++,kcol) = (lf(n)(k)-Rbar(n)(k))/dw(var);
			++kcol;
			x.uht(var)(mode) -= dw(var);
		}
		
		for(int var = 0; var < tri_mesh::ND; ++var) {
			x.pnts(x.seg(sind).pnt(mode))(var) += dx;
			
			element_rsdl(indx,lf);
			
			int krow = 0;
			for(int k=0;k<sm+2;++k)
				for(int n=0;n<x.NV+NV;++n)
					K(krow++,kcol) = (lf(n)(k)-Rbar(n)(k))/dx;
			
			++kcol;
			x.pnts(x.seg(sind).pnt(mode))(var) -= dx;
		}
	}
	
	
	for(int mode = 2; mode < sm+2; ++mode) {
		for(int var = 0; var < x.NV; ++var) {
			x.uht(var)(mode) += dw(var);
			
			element_rsdl(indx,lf);
			
			int krow = 0;
			for(int k=0;k<sm+2;++k)
				for(int n=0;n<x.NV+NV;++n)
					K(krow++,kcol) = (lf(n)(k)-Rbar(n)(k))/dw(var);

			++kcol;
			x.uht(var)(mode) -= dw(var);
		}
		
		for(int var = 0; var < tri_mesh::ND; ++var) {
			crds(indx,mode-2,var) += dx;
			
			element_rsdl(indx,lf);
			
			int krow = 0;
			for(int k=0;k<sm+2;++k)
				for(int n=0;n<x.NV+NV;++n)
					K(krow++,kcol) = (lf(n)(k)-Rbar(n)(k))/dx;
			
			++kcol;
			crds(indx,mode-2,var) -= dx;
		}
	}
		
	return;
}


#ifdef petsc
void hp_coupled_bdry::petsc_jacobian() {
    
    const int sm = basis::tri(x.log2p)->sm();
		const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;

    Array<FLT,2> K(vdofs*(sm+2),vdofs*(sm+2));
    Array<FLT,1> row_store(vdofs*(sm+2));
    Array<int,1> loc_to_glo(vdofs*(sm+2));
    
    
    /* ZERO ROWS CREATED BY R_MESH (REMOVES DIAGONAL ENTRY FOR FIXED B.C.'s) */
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
    MatZeroRows(x.petsc_J,cnt,indices.data(),PETSC_NULL,PETSC_NULL,PETSC_NULL);
#endif

    if (gbl->symmetric || is_master) {
        /* This is effect of variables u,v,p,x,y on */
        /* source terms added to flow residuals */
        /* and x,y mesh movement equations */
        for (int j=0;j<base.nseg;++j) {
            int sind = base.seg(j);
            
            /* CREATE GLOBAL NUMBERING LIST */
            int ind = 0;
            for(int mode = 0; mode < 2; ++mode) {
                int gindx = vdofs*x.seg(sind).pnt(mode);
                for(int var = 0; var < vdofs; ++var)
                    loc_to_glo(ind++) = gindx++;
            }
            
            int gindxNV = x.npnt*vdofs +x.NV*sind*sm;
            int gindxND = jacobian_start +j*tri_mesh::ND*sm;
            for(int mode = 0; mode < sm; ++mode) {
                for(int var = 0; var < x.NV; ++var)
                    loc_to_glo(ind++) = gindxNV++;
                
                for(int var = 0; var < tri_mesh::ND; ++var)
                    loc_to_glo(ind++) = gindxND++;
            }
            element_jacobian(j,K);
#ifdef MY_SPARSE
            x.J.add_values(vdofs*(sm+2),loc_to_glo,vdofs*(sm+2),loc_to_glo,K);
#else
            MatSetValuesLocal(x.petsc_J,vdofs*(sm+2),loc_to_glo.data(),vdofs*(sm+2),loc_to_glo.data(),K.data(),ADD_VALUES);
#endif
        }
    }

    if (sm) {
        for (int j=0;j<base.nseg;++j) {
            int sind = base.seg(j);
            
            /* Now fill in effect of curvature on element resdiuals */
            Array<TinyVector<FLT,MXTM>,1> R(x.NV),Rbar(x.NV),lf_re(x.NV),lf_im(x.NV);        
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
            
            /* INTERIOR	MODES */
            int gindx = x.npnt*vdofs +x.nseg*sm*x.NV +tind*im*x.NV;
            for(int m = 0; m < im; ++m) {
                for(int n = 0; n < x.NV; ++n){
                    loc_to_glo_e(ind++) = gindx++;
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

void hp_coupled_bdry::non_sparse(Array<int,1> &nnzero) {
	const int sm=basis::tri(x.log2p)->sm();
	const int im=basis::tri(x.log2p)->im();
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	const int begin_seg = x.npnt*vdofs;
	const int begin_tri = begin_seg+x.nseg*sm*x.NV;
	
	if(x.sm0 > 0) {
		if (gbl->field_is_coupled) {
			nnzero(Range(jacobian_start,jacobian_start+base.nseg*sm*NV-1)) = 3*vdofs +3*x.NV*sm +x.NV*im +NV*sm;
		}
		else {
			nnzero(Range(jacobian_start,jacobian_start+base.nseg*sm*NV-1)) = vdofs*(sm+2);
		}
		
		/* effect of curvature on flow equations */
		for (int i=0;i<base.nseg;++i) {
			int sind = base.seg(i);
			int tind = x.seg(sind).tri(0);
			if (im) {
				nnzero(Range(begin_tri +tind*im*x.NV,begin_tri +(tind+1)*im*x.NV-1)) += NV*sm;
			}
			
			
			for(int j=0;j<3;++j) {
				int sind1 = x.tri(tind).seg(j);
				nnzero(Range(begin_seg+sind1*x.NV*sm,begin_seg+(sind1+1)*x.NV*sm-1)) += NV*sm;
				if (sind1 == sind) {
					/* For opposing vertex should only be ND*sm for flow */
					/* mesh deformation equation doesn't depend on curvatures */
					int pind = x.tri(tind).pnt(j);
					nnzero(Range(pind*vdofs,pind*vdofs +x.NV -1)) += NV*sm;
				}
			}
			
			for(int endpt=0;endpt<2;++endpt) {
				int pind = x.seg(sind).pnt(endpt);
				if (!gbl->symmetric) {
					if (is_master) {
						 // if I am the master all equations affected by side modes including x,y
						nnzero(Range(pind*vdofs,(pind+1)*vdofs-1)) += NV*sm;
					}
					else {
						// Only flow equations have dependence on local side modes
						// Vertex equations are just zero then copied from other side
						nnzero(Range(pind*vdofs,pind*vdofs +x.NV-1)) += NV*sm;
					}
				}
				else {
					nnzero(Range(pind*vdofs,(pind+1)*vdofs-1)) += NV*sm;
				}
			}
		}
	}
}

void hp_coupled_bdry::non_sparse_snd(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi) {
	
	hp_edge_bdry::non_sparse_snd(nnzero,nnzero_mpi);
	
	if (!base.is_comm()) return;
	
	const int sm=basis::tri(x.log2p)->sm();
	
	/* Last thing to send is nnzero for edge equations */
	if (sm) {
		int count = jacobian_start;
		for (int i=0;i<base.nseg;++i) {
			for (int m=0;m<sm;++m) {
				for(int n=0;n<NV;++n) {
					base.isndbuf(base.sndsize()++) = nnzero(count++);
					
				}
			}
		}
	}
	
	return;
}

int hp_coupled_bdry::non_sparse_rcv(int phase, Array<int,1> &nnzero, Array<int,1> &nnzero_mpi) {
	
	if (!base.is_comm() || base.matchphase(boundary::all_phased,0) != phase) return(0);

	const int sm=basis::tri(x.log2p)->sm();
	
	int count = hp_edge_bdry::non_sparse_rcv(phase, nnzero, nnzero_mpi);
	
	/* Last thing to receive is nnzero for edge equations */
	if (sm) {
		for (int i=base.nseg-1;i>=0;--i) {
			for (int m=0;m<sm;++m) {
				for(int n=0;n<NV;++n) {
					nnzero_mpi(jacobian_start+i*sm*NV +m*NV +n) = base.ircvbuf(0, count++);
				}
			}
		}
	}
	
	if (!gbl->symmetric) {
		if (is_master) {
				/* Add to mpi (below) to allow preconditioning otherwise zero is fine */
			if (sm && !gbl->field_is_coupled) {
				nnzero_mpi(Range(jacobian_start,jacobian_start+base.nseg*sm*NV-1)) = 0;
			}
		}
		else {
			/* For anything other than symmetric, slave side modes are just followers */
			if (sm) {
				nnzero(Range(jacobian_start,jacobian_start+base.nseg*sm*NV-1)) = 1;
				nnzero_mpi(Range(jacobian_start,jacobian_start+base.nseg*sm*NV-1)) = 1;
			}
			
			if (gbl->one_sided) {
				const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;

				/* If one_sided all c0 variables are followers */
				for (int i=base.nseg-1;i>=0;--i) {
					int sind = base.seg(i);
					int pind = x.seg(sind).pnt(1)*vdofs;
					for(std::vector<int>::iterator it=c0_indices_xy.begin();it!=c0_indices_xy.end();++it) {
						nnzero_mpi(pind+*it) = 1;  // Just continuity constraint
						count++;
					}
				}
				int sind = base.seg(0);
				int pind = x.seg(sind).pnt(0)*vdofs;
				for(std::vector<int>::iterator it=c0_indices_xy.begin();it!=c0_indices_xy.end();++it) {
					nnzero_mpi(pind+*it) = 1; // Just continuity constraint
					count++;
				}
				
				if (sm) {
					const int begin_seg = x.npnt*vdofs;
					for (int i=0;i<base.nseg;++i) {
						int sind = base.seg(i);
						for (int mode=0;mode<sm;++mode) {
							for(std::vector<int>::iterator it=c0_indices.begin();it!=c0_indices.end();++it) {
								nnzero_mpi(begin_seg+sind*x.NV*sm +mode*x.NV +*it) = 1;  // Just continuity constraint
							}
						}
					}
				}
			}
		}
	}
	
	// make sure flow and coupled side variables have same number of entries
	if (is_master && gbl->field_is_coupled && gbl->precondition) {
		const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
		const int c0var = c0_indices[0];
		
		for (int i=base.nseg-1;i>=0;--i) {
			int sind = base.seg(i);
			int pind = x.seg(sind).pnt(1)*vdofs;
			int flow_vars = nnzero_mpi(pind+c0var);
			for(int n=x.NV;n<vdofs;++n) {
				nnzero_mpi(pind+n) = flow_vars;
			}
		}
		int sind = base.seg(0);
		int pind = x.seg(sind).pnt(0)*vdofs;
		int flow_vars = nnzero_mpi(pind+c0var);
		for(int n=x.NV;n<vdofs;++n) {
			nnzero_mpi(pind+n) = flow_vars;
		}

		if (sm) {
			const int begin_seg = x.npnt*vdofs;

			for (int i=0;i<base.nseg;++i) {
				int sind = base.seg(i);
				for (int mode=0;mode<sm;++mode) {
					int flow_vars = nnzero_mpi(begin_seg+sind*x.NV*sm +mode*x.NV +c0var);
					for(int n=0;n<NV;++n) {
						nnzero_mpi(jacobian_start+i*sm*NV +mode*NV +n) = flow_vars;
					}
				}
			}
		}
	}
	
	return(count);
}

void hp_coupled_bdry::petsc_matchjacobian_snd() {
	
	hp_edge_bdry::petsc_matchjacobian_snd();
	if (!gbl->symmetric) return;
	

	/* This is extra stuff for symmetric matching not used otherwise */
	const int sm = basis::tri(x.log2p)->sm();
	if (sm) {
		const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;

		// Send extra stuff for edge modes */
		int row = jacobian_start;
		for(int i=0;i<base.nseg;++i) {
			int sind = base.seg(i);
			int rowbase = x.npnt*vdofs +sind*x.NV*x.sm0;
			for(int m=0;m<sm;++m) {
				for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
					/* attach diagonal column # to allow continuity enforcement */
					base.fsndbuf(base.sndsize()++) = rowbase +*n +0.1;
				}
				rowbase += x.NV;

				/* Send side variables */
				for(int n=0;n<NV;++n) {
					/* attach diagonal column # to allow continuity enforcement */
					base.fsndbuf(base.sndsize()++) = row +0.1;
					base.fsndbuf(base.sndsize()++) = x.J._cpt(row+1) -x.J._cpt(row) +0.1;
#ifdef MPDEBUG
					*x.gbl->log << "sending " << x.J._cpt(row+1) -x.J._cpt(row) << " jacobian entries for edge " << i << " and variable " << n << std::endl;
#endif
					for (int col=x.J._cpt(row);col<x.J._cpt(row+1);++col) {
#ifdef MPDEBUG
						*x.gbl->log << x.J._col(col) << ' ';
#endif
						base.fsndbuf(base.sndsize()++) = x.J._col(col) +0.1;
						base.fsndbuf(base.sndsize()++) = x.J._val(col);
					}
#ifdef MPDEBUG
					*x.gbl->log << std::endl;
#endif
					++row;
				}
			}
		}
	}
}

int hp_coupled_bdry::petsc_matchjacobian_rcv(int phase) {
	// slave only receives indices of continuous variables
	// master receives only information for flow not mesh variables
	
	if (!base.is_comm() || base.matchphase(boundary::all_phased,0) != phase || c0_indices_xy.empty()) return(0);
    
	const int sm = x.sm0;
	int count = 0;

	if (is_master) {
		/* Master receives jacobian entries for c0 variables */
		count = hp_edge_bdry::petsc_matchjacobian_rcv(phase);
	}
	else {
		if (!gbl->one_sided) {
			count = hp_edge_bdry::petsc_matchjacobian_rcv(phase);
		}
		else {
			/* Apply matching constraint for c0 vars */
			/* This just skips the information sent and puts a -1 in J_mpi */
			
			int Jstart_mpi = static_cast<int>(base.frcvbuf(0, count++)); // Start of jacobian on matching block
			count++; // Skip index of boundary unknowns on mathcing block.   Not used in this routine
			
			sparse_row_major *pJ_mpi;
			if (base.is_local(0)) {
				pJ_mpi = &x.J;
				Jstart_mpi = 0;
			}
			else {
				pJ_mpi = &x.J_mpi;
			}
			
			const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
			
			/* Now do stuff for communication boundaries */
			int row;
			/* Now Receive Information */
			for (int i=base.nseg-1;i>=0;--i) {
				int sind = base.seg(i);
				int rowbase = x.seg(sind).pnt(1)*vdofs;
				vector<int> row_mpi_storage;
				for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
					row = rowbase + *n;
					int row_mpi = static_cast<int>(base.frcvbuf(0,count++)) +Jstart_mpi;
					row_mpi_storage.push_back(row_mpi);
					int ncol = static_cast<int>(base.frcvbuf(0,count++));
#ifdef MPDEBUG
					*x.gbl->log << "receiving " << ncol << " jacobian entries for vertex " << row/vdofs << " and variable " << *n << std::endl;
#endif
					for (int k = 0;k<ncol;++k) {
						count++; // skip column number
						count++; // skip value;
					}
				}
				/* Set equality constraint */
				std::vector<int>::iterator row_mpi=row_mpi_storage.begin();
				for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
					row = rowbase + *n;
					pJ_mpi->zero_row(row);
					x.J.zero_row(row);
					x.J.set_values(row, row, 1.0);
					pJ_mpi->set_values(row,*row_mpi,-1.0);
					++row_mpi;
				}
				row_mpi_storage.clear();
				
				/* Now receive side Jacobian information */
				rowbase = x.npnt*vdofs +sind*x.NV*x.sm0;
				for(int mode=0;mode<x.sm0;++mode) {
					for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
						row = rowbase +*n;
						int row_mpi = static_cast<int>(base.frcvbuf(0,count++)) +Jstart_mpi;
						row_mpi_storage.push_back(row_mpi);
						int ncol = static_cast<int>(base.frcvbuf(0,count++));
#ifdef MPDEBUG
						*x.gbl->log << "receiving " << ncol << " jacobian entries for side " << sind << " and variable " << *n << std::endl;
#endif
						for (int k = 0;k<ncol;++k) {
							count++; // Skip column number
							count++; // Skip value
						}
					}
					rowbase += x.NV;
				}
				
				/* Equality of C0_vars */
				rowbase = x.npnt*vdofs +sind*x.NV*x.sm0;
				row_mpi=row_mpi_storage.begin();
				int sgn_mpi = 1;
				for(int mode=0;mode<x.sm0;++mode) {
					for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
						row = rowbase +mode*x.NV +*n;
						pJ_mpi->zero_row(row);
						x.J.zero_row(row);
						x.J.set_values(row, row, 1.0);
						pJ_mpi->set_values(row,*row_mpi,-sgn_mpi);
						++row_mpi;
					}
					sgn_mpi *= -1;
				}
				row_mpi_storage.clear();
			}
			int sind = base.seg(0);
			int rowbase = x.seg(sind).pnt(0)*vdofs;
			vector<int> row_mpi_storage;
			for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
				row = rowbase + *n;
				int row_mpi = static_cast<int>(base.frcvbuf(0,count++)) +Jstart_mpi;
				row_mpi_storage.push_back(row_mpi);
				int ncol = static_cast<int>(base.frcvbuf(0,count++));
#ifdef MPDEBUG
				*x.gbl->log << "receiving " << ncol << " jacobian entries for vertex " << row/vdofs << " and variable " << *n << std::endl;
#endif
				for (int k = 0;k<ncol;++k) {
					count++; // skip column number
					count++; // skip value;
				}
			}
			/* Set equality constraint */
			std::vector<int>::iterator row_mpi=row_mpi_storage.begin();
			for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
				row = rowbase + *n;
				pJ_mpi->zero_row(row);
				x.J.zero_row(row);
				x.J.set_values(row, row, 1.0);
				pJ_mpi->set_values(row,*row_mpi,-1.0);
				++row_mpi;
			}
			row_mpi_storage.clear();
		}
	}

	if (!gbl->symmetric) {
		/* Apply matching curvature constraint */
		if (sm && !is_master) {
			int Jstart_mpi = static_cast<int>(base.frcvbuf(0,0));
			sparse_row_major *pJ_mpi;
			if (base.is_local(0)) {
				pJ_mpi = &x.J;
				Jstart_mpi = 0;
			}
			else {
				pJ_mpi = &x.J_mpi;
			}
			
			int ind = jacobian_start;
			int ind_mpi = static_cast<int>(base.frcvbuf(0,1)) +(base.nseg-1)*sm*x.ND +Jstart_mpi;
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
	else {
		/* symmetric jacobian transmission of side mode equations */
		if (sm) {
			int Jstart_mpi = static_cast<int>(base.frcvbuf(0, 0)); // Start of jacobian on matching block
			sparse_row_major *pJ_mpi;
			if (base.is_local(0)) {
				pJ_mpi = &x.J;
				Jstart_mpi = 0;
			}
			else {
				pJ_mpi = &x.J_mpi;
			}
			vector<int> row_mpi_storage;
			
			const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
			
			// Receive extra stuff for edge modes */
			for (int i=base.nseg-1;i>=0;--i) {
				int row = jacobian_start +i*x.sm0*NV;
				int sgn = 1;
				for(int m=0;m<x.sm0;++m) {
					for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
						/* diagonal column # to allow continuity enforcement */
						int row_mpi = static_cast<int>(base.frcvbuf(0,count++)) +Jstart_mpi;
						row_mpi_storage.push_back(row_mpi);
					}
					
					/* Receive side variables */
					for(int n=0;n<NV;++n) {
						int row_mpi = static_cast<int>(base.frcvbuf(0,count++)) +Jstart_mpi;
						row_mpi_storage.push_back(row_mpi);
						int ncol = static_cast<int>(base.frcvbuf(0,count++));
#ifdef MPDEBUG
						*x.gbl->log << "receiving " << ncol << " jacobian entries for side " << i << " and variable " << n << std::endl;
#endif
						for (int k = 0;k<ncol;++k) {
							int col = static_cast<int>(base.frcvbuf(0,count++));
							FLT val = sgn*base.frcvbuf(0,count++);
							if (col < INT_MAX-10 && col > -1) {
								col += Jstart_mpi;
#ifdef MPDEBUG
								*x.gbl->log  << col << ' ';
#endif
								(*pJ_mpi).add_values(row,col,val);
							}
						}
#ifdef MPDEBUG
						*x.gbl->log << std::endl;
#endif
						++row;
					}
					sgn *= -1;
				}
				
				int rowbase = jacobian_start +i*x.sm0*NV;
				int sind = base.seg(i);
				int flowbase = x.npnt*vdofs +sind*x.NV*x.sm0;

				for(int mode=0;mode<x.sm0;++mode) {
					for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
						row = flowbase +mode*x.NV+*n;
						
						int sgn_mpi = 1;
						std::vector<int>::iterator row_mpi=row_mpi_storage.begin();
						for(int mode_mpi=0;mode_mpi<x.sm0;++mode_mpi) {
							for(std::vector<int>::iterator n_mpi=c0_indices.begin();n_mpi != c0_indices.end();++n_mpi) {
#ifdef MPDEBUG
								*x.gbl->log << "side swapping " << row << ',' << *row_mpi << " for " << row << ',' << flowbase +mode_mpi*x.NV +*n_mpi << std::endl;
#endif
#ifdef DEBUG_JAC
								if (!x.gbl->jac_debug)
#endif
								{
									FLT dval = (*pJ_mpi)(row,*row_mpi);
									(*pJ_mpi)(row,*row_mpi) = 0.0;
									x.J(row,flowbase +mode_mpi*x.NV +*n_mpi) += sgn_mpi*dval;
								}
								++row_mpi;
							}
							
							for(int n_mpi=0;n_mpi<NV;++n_mpi) {
#ifdef MPDEBUG
								*x.gbl->log << "side swapping " << row << ',' << *row_mpi << " for " << row << ',' << rowbase +mode_mpi*NV +n_mpi << std::endl;
#endif
#ifdef DEBUG_JAC
								if (!x.gbl->jac_debug)
#endif
								{
									FLT dval = (*pJ_mpi)(row,*row_mpi);
									(*pJ_mpi)(row,*row_mpi) = 0.0;
									x.J(row,rowbase +mode_mpi*NV +n_mpi) += sgn_mpi*dval;
								}
								++row_mpi;
							}
							sgn_mpi *= -1;
						}
					}

					for(int n=0;n<NV;++n) {
						row = rowbase +mode*NV +n;
						int sgn_mpi = 1;
						std::vector<int>::iterator row_mpi=row_mpi_storage.begin();
						for(int mode_mpi=0;mode_mpi<x.sm0;++mode_mpi) {
							for(std::vector<int>::iterator n_mpi=c0_indices.begin();n_mpi != c0_indices.end();++n_mpi) {
#ifdef MPDEBUG
								*x.gbl->log << "side swapping " << row << ',' << *row_mpi << " for " << row << ',' << flowbase +mode_mpi*x.NV +*n_mpi << std::endl;
#endif
#ifdef DEBUG_JAC
								if (!x.gbl->jac_debug)
#endif
								{
									FLT dval = (*pJ_mpi)(row,*row_mpi);
									(*pJ_mpi)(row,*row_mpi) = 0.0;
									x.J(row,flowbase +mode_mpi*x.NV +*n_mpi) += sgn_mpi*dval;
								}
								++row_mpi;
							}

							for(int n_mpi=0;n_mpi<NV;++n_mpi) {
#ifdef MPDEBUG
								*x.gbl->log << "side swapping " << row << ',' << *row_mpi << " for " << row << ',' << rowbase +mode_mpi*NV +n_mpi << std::endl;
#endif
#ifdef DEBUG_JAC
								if (!x.gbl->jac_debug)
#endif
								{
									FLT dval = (*pJ_mpi)(row,*row_mpi);
									(*pJ_mpi)(row,*row_mpi) = 0.0;
									x.J(row,rowbase +mode_mpi*NV +n_mpi) += sgn_mpi*dval;
								}
								++row_mpi;
							}
							sgn_mpi *= -1;
						}
						x.J.multiply_row(row,0.5);
						x.J_mpi.multiply_row(row,0.5);
					}
				}
				row_mpi_storage.clear();
			}
		}
	}
	return(count);
}

void hp_coupled_bdry::petsc_make_1D_rsdl_vector(Array<double,1> res) {
	const int sm = basis::tri(x.log2p)->sm();
	int ind = jacobian_start;
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	const int ncoupled = NV +gbl->field_is_coupled*c0_indices.size();
	Array<FLT,1> res_vec(ncoupled);
	

	if ((gbl->symmetric || is_master)) {
		if (gbl->precondition)  {
			for(int j=0;j<base.nseg;++j) {
				int ind1 = x.npnt*vdofs +base.seg(j)*sm*x.NV;

				for(int m=0;m<sm;++m) {
					int indx = 0;
					if (gbl->field_is_coupled) {
						for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
							res_vec(indx++) = res(ind1+*n);
						}
					}
					for(int n=0;n<NV;++n) {
						res_vec(indx++) = gbl->sres(j,m,n);
					}
					
					int info;
					char trans[] = "T";
#ifdef F2CFortran
					GETRS(trans,ncoupled,1,&gbl->sdt2(j,m,0,0),gbl->sdt2.length(thirdDim),&gbl->spiv2(j,m,0),res_vec.data(),ncoupled,info);
#else
                    const int one = 1, length = gbl->sdt2.length(thirdDim);
                    dgetrs_(trans,&ncoupled,&one,&gbl->sdt2(j,m,0,0),&length,&gbl->spiv2(j,m,0),res_vec.data(),&ncoupled,&info);
#endif
					if (info != 0) {
						*x.gbl->log << "DGETRS FAILED IN SIDE MODE PRECONDITIONER " << info << std::endl;
						sim::abort(__LINE__,__FILE__,x.gbl->log);
					}

					indx = 0;
					if (gbl->field_is_coupled) {
						for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
							res(ind1+*n) = res_vec(indx++);
						}
					}
					for(int n=0;n<NV;++n) {
						res(ind++) = res_vec(indx++);
					}
					ind1 += x.NV;
				}
			}
		}
		else {
			for(int j=0;j<base.nseg;++j) {
				for(int m=0;m<sm;++m) {
					for(int n=0;n<NV;++n) {
						res(ind++) = gbl->sres(j,m,n);
					}
				}
			}
		}
	}
	else {
		/* Side mode equality constraint */
		for(int j=0;j<base.nseg;++j) {
			int ind1 = x.npnt*vdofs +base.seg(j)*sm*x.NV;
			for(int m=0;m<sm;++m) {
				for(int n=0;n<NV;++n)
					res(ind++) = 0.0; // all local variables assumed to be equal on both sides
				
				if (gbl->one_sided) {  // if one_sided then continuous flow variables are forced to be continuous
					for(std::vector<int>::iterator it = c0_indices.begin();it != c0_indices.end();++it)
						res(ind1+*it) = 0.0;
					ind1 += x.NV;
				}
			}
		}
	}
	
	/* Swap kinetic and energy vertex residuals */
	if ((!gbl->one_sided || is_master)) {
		if (gbl->precondition) {
			int i = 0;
			int sind;
			do {
				sind = base.seg(i);
				int row = x.seg(sind).pnt(0)*vdofs;
				int indx = 0;
				if (gbl->field_is_coupled) {
					for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
						res_vec(indx++) = res(row +*n);
					}
				}
				for(int n=0;n<NV;++n) {
					res_vec(indx++) = res(row +n +x.NV);
				}

				int info;
				char trans[] = "T";
#ifdef F2CFortran
				GETRS(trans,ncoupled,1,&gbl->vdt(i,0,0),gbl->vdt.length(secondDim),&gbl->vpiv(i,0),res_vec.data(),ncoupled,info);
#else
                const int one = 1, length = gbl->vdt.length(secondDim);
                dgetrs_(trans,&ncoupled,&one,&gbl->vdt(i,0,0),&length,&gbl->vpiv(i,0),res_vec.data(),&ncoupled,&info);
#endif
				if (info != 0) {
					*x.gbl->log << "DGETRS FAILED IN VERTEX PRECONDITIONER " << info << std::endl;
					sim::abort(__LINE__,__FILE__,x.gbl->log);
				}
				
				indx = 0;
				if (gbl->field_is_coupled) {
					for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
						res(row +*n) = res_vec(indx++);
					}
				}
				for(int n=0;n<NV;++n) {
					res(row +n +x.NV) = res_vec(indx++);
				}
			} while (++i < base.nseg);
			int row = x.seg(sind).pnt(1)*vdofs;
			int indx = 0;
			if (gbl->field_is_coupled) {
				for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
					res_vec(indx++) = res(row +*n);
				}
			}
			for(int n=0;n<NV;++n) {
				res_vec(indx++) = res(row +n +x.NV);
			}

			int info;
			char trans[] = "T";
#ifdef F2CFortran
			GETRS(trans,ncoupled,1,&gbl->vdt(i,0,0),gbl->vdt.length(secondDim),&gbl->vpiv(i,0),res_vec.data(),ncoupled,info);
#else
            const int one = 1, length = gbl->vdt.length(secondDim);
            dgetrs_(trans,&ncoupled,&one,&gbl->vdt(i,0,0),&length,&gbl->vpiv(i,0),res_vec.data(),&ncoupled,&info);
#endif
			if (info != 0) {
				*x.gbl->log << "DGETRS FAILED IN VERTEX PRECONDITIONER " << info << std::endl;
				sim::abort(__LINE__,__FILE__,x.gbl->log);
			}
			
			indx = 0;
			if (gbl->field_is_coupled) {
				for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
					res(row +*n) = res_vec(indx++);
				}
			}
			for(int n=0;n<NV;++n) {
				res(row +n +x.NV) = res_vec(indx++);
			}
		}
	}
	else {
		/* SET ALL C0 VARS TO FOLLOW */
		const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
		int i = 0;
		int sind;
		do {
			sind = base.seg(i);
			int row = x.seg(sind).pnt(0)*vdofs;
			for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
				res(row+*n) = 0.0;
			}
		} while (++i < base.nseg);
		int row = x.seg(sind).pnt(1)*vdofs;
		for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
			res(row+*n) = 0.0;
		}
	}
}

void hp_coupled_bdry::petsc_premultiply_jacobian() {
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	const int ncoupled = NV +gbl->field_is_coupled*c0_indices.size();
	Array<int,1> row_ind(ncoupled);
	
	/* Swap side equations */
	if ((gbl->symmetric || is_master) && gbl->precondition) {
		int row0 = jacobian_start; // Index of motion equations
		
		const int sm = basis::tri(x.log2p)->sm();
		for(int j=0;j<base.nseg;++j) {
			int row1 = x.npnt*vdofs +base.seg(j)*sm*x.NV;
			for(int m=0;m<sm;++m) {
				
				int ind = 0;
				if (gbl->field_is_coupled) {
					for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
						row_ind(ind++) = row1+*n;
					}
				}
				for(int n=0;n<NV;++n) {
					row_ind(ind++) = row0+n;
				}

				// Store Jacobian
				gbl->sdt2(j,m,Range::all(),Range::all()) = 0.0;
				if (gbl->precondition) {
					for(int row=0;row<ncoupled;++row) {
						for(int col=0;col<ncoupled;++col) {
							gbl->sdt2(j,m,row,col) = x.J(row_ind(row),row_ind(col));
						}
					}
				}
				else {
					// To shut off preconditioning
					for(int row=0;row<ncoupled;++row) {
						gbl->sdt2(j,m,row,row) = 1.0;
					}
				}

				int info;
#ifdef F2CFortran
				GETRF(ncoupled,ncoupled,&gbl->sdt2(j,m,0,0),gbl->sdt2.length(thirdDim),&gbl->spiv2(j,m,0),info);
#else
                const int length = gbl->sdt2.length(thirdDim);
                dgetrf_(&ncoupled,&ncoupled,&gbl->sdt2(j,m,0,0),&length,&gbl->spiv2(j,m,0),&info);
#endif
				if (info != 0) {
					*x.gbl->log << "DGETRF FAILED IN SIDE MODE PRECONDITIONER " << info << std::endl;
					sim::abort(__LINE__,__FILE__,x.gbl->log);
				}
				Array<FLT,2> A = gbl->sdt2(j,m,Range::all(),Range::all());
				Array<int,1> ipiv = gbl->spiv2(j,m,Range::all());
				x.J.match_patterns(ncoupled,row_ind);
				x.J.combine_rows(ncoupled, row_ind, A, gbl->sdt2.length(thirdDim), ipiv);
				x.J_mpi.match_patterns(ncoupled, row_ind);
				x.J_mpi.combine_rows(ncoupled, row_ind, A, gbl->sdt2.length(thirdDim), ipiv);
				row0 += NV;
				row1 += x.NV;
			}
		}
	}
	
	if ((!gbl->one_sided || is_master) && gbl->precondition) {
		/* Swap rows */
		int i = 0;
		int sind;
		do {
			sind = base.seg(i);
			int row = x.seg(sind).pnt(0)*vdofs;
			int ind = 0;
			if (gbl->field_is_coupled) {
				for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
					row_ind(ind++) = row +*n;
				}
			}
			for(int n=0;n<NV;++n) {
				row_ind(ind++) = row +n +x.NV;
			}

			// Store Jacobian
			gbl->vdt(i,Range::all(),Range::all()) = 0.0;
			if (gbl->precondition) {
				for(int row=0;row<ncoupled;++row) {
					for(int col=0;col<ncoupled;++col) {
						gbl->vdt(i,row,col) = x.J(row_ind(row),row_ind(col));
					}
				}
			}
			else {
				// To shut off preconditioning
				for(int row=0;row<ncoupled;++row) {
					gbl->vdt(i,row,row) = 1.0;
				}
			}
			int info;
#ifdef F2CFortran
			GETRF(ncoupled,ncoupled,&gbl->vdt(i,0,0),gbl->vdt.length(secondDim),&gbl->vpiv(i,0),info);
#else
            const int length = gbl->vdt.length(secondDim);
            dgetrf_(&ncoupled,&ncoupled,&gbl->vdt(i,0,0),&length,&gbl->vpiv(i,0),&info);

#endif
			if (info != 0) {
				*x.gbl->log << "DGETRF FAILED IN VERTEX MODE PRECONDITIONER\n";
				sim::abort(__LINE__,__FILE__,x.gbl->log);
			}
			Array<FLT,2> A = gbl->vdt(i,Range::all(),Range::all());
			Array<int,1> vpiv = gbl->vpiv(i,Range::all());
			x.J.match_patterns(ncoupled, row_ind);
			x.J.combine_rows(ncoupled, row_ind, A, gbl->vdt.length(secondDim), vpiv);
			x.J_mpi.match_patterns(ncoupled, row_ind);
			x.J_mpi.combine_rows(ncoupled, row_ind, A, gbl->vdt.length(secondDim), vpiv);
			
		} while (++i < base.nseg);
		int row = x.seg(sind).pnt(1)*vdofs;
		int ind = 0;
		if (gbl->field_is_coupled) {
			for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
				row_ind(ind++) = row +*n;
			}
		}
		for(int n=0;n<NV;++n) {
			row_ind(ind++) = row +n +x.NV;
		}
		
		// Store Jacobian
		gbl->vdt(i,Range::all(),Range::all()) = 0.0;
		if (gbl->precondition) {
			for(int row=0;row<ncoupled;++row) {
				for(int col=0;col<ncoupled;++col) {
					gbl->vdt(i,row,col) = x.J(row_ind(row),row_ind(col));
				}
			}
		}
		else {
			// To shut off preconditioning
			for(int row=0;row<ncoupled;++row) {
				gbl->vdt(i,row,row) = 1.0;
			}
		}
		int info;
#ifdef F2CFortran
		GETRF(ncoupled,ncoupled,&gbl->vdt(i,0,0),gbl->vdt.length(secondDim),&gbl->vpiv(i,0),info);
#else
        const int length = gbl->vdt.length(secondDim);
        dgetrf_(&ncoupled,&ncoupled,&gbl->vdt(i,0,0),&length,&gbl->vpiv(i,0),&info);
#endif
		if (info != 0) {
			*x.gbl->log << "DGETRF FAILED IN VERTEX MODE PRECONDITIONER\n";
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}
		Array<FLT,2> A = gbl->vdt(i,Range::all(),Range::all());
		Array<int,1> vpiv = gbl->vpiv(i,Range::all());
		x.J.match_patterns(ncoupled, row_ind);
		x.J.combine_rows(ncoupled, row_ind, A, gbl->vdt.length(secondDim), vpiv);
		x.J_mpi.match_patterns(ncoupled, row_ind);
		x.J_mpi.combine_rows(ncoupled, row_ind, A, gbl->vdt.length(secondDim), vpiv);
	
	}
}
#endif

void hp_deformable_fixed_pnt::init(input_map& inmap,void* gbl_in) {
	std::string keyword,val;
	std::istringstream data;
	std::string filename;
	
	multi_physics_pnt::init(inmap,gbl_in);

	
	if ((surf = dynamic_cast<hp_coupled_bdry *>(x.hp_ebdry(base.ebdry(0))))) {
		surfbdry = 0;
	}
	else if ((surf = dynamic_cast<hp_coupled_bdry *>(x.hp_ebdry(base.ebdry(1))))) {
		surfbdry = 1;
	}
	else {
		surfbdry = -1;  // Just a fixed point (used for mesh touching complicated point in 3-phase flow with free-surface)
		surf = 0;
	}

	/* Check if set manually already, otherwise use other boundary to get defaults */
    if (surfbdry > -1) {
        Array<int,1> atemp(x.NV);
        if (!inmap.get(base.idprefix+"_hp_typelist", atemp.data(), x.NV)) {
            for (int n=0;n<x.NV;++n) {
                
                // Check for essential B.C.'s on either edge-boundary
                type[n] = static_cast<bctypes>(x.hp_ebdry(base.ebdry(1-surfbdry))->type[n]);
                if (type[n] == essential) {
                    essential_indices.push_back(n);
                    continue;
                }
                
                type[n] = static_cast<bctypes>(x.hp_ebdry(base.ebdry(surfbdry))->type[n]);
                if (type[n] == essential) {
                    essential_indices.push_back(n);
                    continue;
                }
            }
        }
    }
}

void hp_deformable_fixed_pnt::vdirichlet() {
	const int nfix = tri_mesh::ND;
	
	// APPLY FLOW B.C.'S
	multi_physics_pnt::vdirichlet();
	
	if (surf) {
		if (surf->is_master) {
			if (surfbdry == 0) {
				for(int n=0;n<nfix;++n) {
					surf->gbl->vres(x.ebdry(base.ebdry(0))->nseg,n) = 0.0;
				}
			}
			else if (surfbdry == 1) {
				for(int n=0;n<nfix;++n) {
					surf->gbl->vres(0,n) = 0.0;
				}
			}
		}
	}
#ifdef petsc
	r_tri_mesh::global *r_gbl = dynamic_cast<r_tri_mesh::global *>(x.gbl);
	for(int n=0;n<nfix;++n) {
		r_gbl->res(base.pnt)(n) = 0.0;
	}
#endif
}

#ifdef petsc
void hp_deformable_fixed_pnt::petsc_jacobian_dirichlet() {
	
	multi_physics_pnt::petsc_jacobian_dirichlet();
	
	const int nfix = tri_mesh::ND;
	/* BOTH X & Y ARE FIXED */
	Array<int,1> rows(tri_mesh::ND);
	for(int n=0;n<nfix;++n)
		rows(n) = (x.NV+tri_mesh::ND)*base.pnt +x.NV +n;
	
#ifdef MY_SPARSE
	x.J.zero_rows(nfix,rows);
	x.J_mpi.zero_rows(nfix,rows);
	x.J.set_diag(nfix,rows,1.0);
#else
	MatZeroRows(x.petsc_J,nfix,rows.data(),1.0,PETSC_NULL,PETSC_NULL);
#endif
}

int hp_deformable_fixed_pnt::setup_preconditioner() {
	const int nfix = tri_mesh::ND;

	/* Turn off preconditioner for r_tri_mesh */
	x.r_tri_mesh::gbl->diag(base.pnt) = 1.0;
	
	/* Turn off side preconditioner as well */
	if (surf) {
		if (surf->is_master) {
			if (surfbdry == 0) {
				for(int n=0;n<nfix;++n) {
					surf->gbl->vdt(x.ebdry(base.ebdry(0))->nseg,n,Range::all()) = 0.0;
					surf->gbl->vdt(x.ebdry(base.ebdry(0))->nseg,n,n) = 1.0;
				}
			}
			else if (surfbdry == 1) {
				for(int n=0;n<nfix;++n) {
					surf->gbl->vdt(0,n,Range::all()) = 0.0;
					surf->gbl->vdt(0,n,n) = 1.0;
				}
			}
		}
	}
    return(0);
}

#endif
	
void hp_deformable_free_pnt::init(input_map& inmap,void* gbl_in) {
	hp_deformable_fixed_pnt::init(inmap,gbl_in);
	
	std::string input;
	if (inmap.get(base.idprefix + "_wall_type",input)) {
		if (input == "vertical") {
			wall_type = vertical;
			position = x.pnts(base.pnt)(0);
		}
		else if (input == "horizontal") {
			wall_type = horizontal;
			position = x.pnts(base.pnt)(1);
		}
        else if (input == "curved") {
			wall_type = curved;
            int bnumwall = base.ebdry(1-surfbdry);
            if (!x.hp_ebdry(bnumwall)->is_curved()) {
                *x.gbl->log << "Can't use curved hp_deformable_free_pnt on non_curved boundary" << std::endl;
                sim::abort(__LINE__,__FILE__,&std::cerr);
            }
        }
		else {
			*x.gbl->log << "Unrecognized wall type" << std::endl;
		}
	}
	else {
		wall_type = vertical;
		position = x.pnts(base.pnt)(0);
	}
}

void hp_deformable_free_pnt::mvpttobdry(TinyVector<FLT,tri_mesh::ND> &pt) {
	if (surf->is_master || surf->gbl->symmetric) {
		switch(wall_type) {
			case vertical: {
				x.pnts(base.pnt)(0) = position;
				break;
			}
			case horizontal: {
				x.pnts(base.pnt)(1)  = position;
				break;
			}
			case curved: {
				int bnumwall = base.ebdry(1-surfbdry);
				TinyVector<FLT,tri_mesh::ND> tgt = x.pnts(base.pnt);
				x.ebdry(bnumwall)->mvpttobdry(x.ebdry(bnumwall)->nseg-1, 1.0, tgt);
				break;
			}
		}
	}
	
	return;
}

void hp_deformable_free_pnt::element_rsdl(Array<FLT,1> lf) {
	lf = 0.0;
	if (surf->is_master || surf->gbl->symmetric) {
		switch(wall_type) {
			case vertical: {
				lf(x.NV) = x.pnts(base.pnt)(0) -position;
				break;
			}
			case horizontal: {
				lf(x.NV) = x.pnts(base.pnt)(1) -position;
				break;
			}
			case curved: {
                /* Calculate normal from geometry definition */
                TinyVector<FLT,tri_mesh::ND> wall_normal;
                if (surfbdry == 0) {
                    x.ebdry(base.ebdry(1))->bdry_normal(0,-1.0,wall_normal);
                }
                else {
                    x.ebdry(base.ebdry(0))->bdry_normal(x.ebdry(base.ebdry(0))->nseg-1,1.0,wall_normal);
                }

				int bnumwall = base.ebdry(1-surfbdry);
				TinyVector<FLT,tri_mesh::ND> tgt = x.pnts(base.pnt);
                if(surfbdry == 0)
                    x.ebdry(bnumwall)->mvpttobdry(0, 1.0, tgt);
                else
                    x.ebdry(bnumwall)->mvpttobdry(x.ebdry(bnumwall)->nseg-1, -1.0, tgt);
				TinyVector<FLT,tri_mesh::ND> diff = x.pnts(base.pnt)-tgt;
                lf(x.NV) = dot(diff,wall_normal);
				break;
			}
		}
	}
	
	return;
}

/* Routine to add surface tension stress or endpoint movement residual */
void hp_deformable_free_pnt::rsdl(int stage) {
	
	if (!surf->is_master && !surf->gbl->symmetric) return;
	
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	Array<FLT,1> lf(vdofs);
	element_rsdl(lf);
	
	x.gbl->res.v(base.pnt,Range::all()) += lf(Range(0,x.NV-1));
	
	int endpt;
	if (surfbdry == 0)
		endpt = x.ebdry(base.ebdry(0))->nseg;
	else
		endpt = 0;
	
#ifdef petsc
	r_tri_mesh::global *r_gbl = dynamic_cast<r_tri_mesh::global *>(x.gbl);
	/* equation for tangential poistion */
	r_gbl->res(base.pnt)(0) = lf(x.NV);
#endif
	
	// FIXME: This should be deleted eventually
	surf->gbl->vres(endpt,0) = lf(x.NV);
}

#ifdef petsc
void hp_deformable_free_pnt::petsc_jacobian() {
	
	/* ZERO TANGENT MESH MOVEMENT ROW */
	int row = (x.NV+tri_mesh::ND)*base.pnt +x.NV;
#ifdef MY_SPARSE
	x.J.zero_row(row);
	x.J_mpi.zero_row(row);
#else
	/* Must zero rows of jacobian created by r_mesh */
	MatAssemblyBegin(x.petsc_J,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(x.petsc_J,MAT_FINAL_ASSEMBLY);
	MatZeroRow(x.petsc_J,row,PETSC_NULL);
#endif

	/* Use Jacobian routine to fill in Jacobian terms */
	multi_physics_pnt::petsc_jacobian();
}
#endif


void translating_surface::init(input_map& inmap, void *gin) {
	hp_coupled_bdry::init(inmap,gin);
	if (!inmap.get(base.idprefix + "_velx",vel(0))) inmap.getwdefault("velx",vel(0),1.0);
	if (!inmap.get(base.idprefix + "_vely",vel(0))) inmap.getwdefault("vely",vel(1),0.0);
}
	
void translating_surface::element_rsdl(int indx, Array<TinyVector<FLT,MXTM>,1> lf) {
	
	if (!is_master) return;
	
	lf = 0.0;

	/* Calculate any specified fluxes to be added */
	hp_edge_bdry::element_rsdl(indx, lf);
	
	int i,n,sind;
	TinyVector<FLT,tri_mesh::ND> norm, rp;
	FLT jcb;
	TinyVector<FLT,tri_mesh::ND> u;
	TinyMatrix<FLT,tri_mesh::ND,MXGP> crd, dcrd, mvel;
	TinyMatrix<FLT,8,MXGP> res;
	
	sind = base.seg(indx);	
	x.crdtocht1d(sind);
	for(n=0;n<tri_mesh::ND;++n)
		basis::tri(x.log2p)->proj1d(&x.cht(n,0),&crd(n,0),&dcrd(n,0));
	
	for(i=0;i<basis::tri(x.log2p)->gpx();++i) {
		norm(0) =  dcrd(1,i);
		norm(1) = -dcrd(0,i);
		jcb = sqrt(norm(0)*norm(0) +norm(1)*norm(1));
		
		/* RELATIVE VELOCITY STORED IN MVEL(N)*/
		for(n=0;n<tri_mesh::ND;++n) {
			mvel(n,i) =  vel(n) -(x.gbl->bd(0)*(crd(n,i) -dxdt(x.log2p,indx)(n,i)));
#ifdef MESH_REF_VEL
			mvel(n,i) -= x.gbl->mesh_ref_vel(n);
#endif
		}
		
		/* TANGENTIAL SPACING */
		res(0,i) = -ksprg(indx)*jcb;
		/* NORMAL FLUX */
		res(1,i) = -RAD(crd(0,i))*(mvel(0,i)*norm(0) +mvel(1,i)*norm(1));
		/* UPWINDING BASED ON TANGENTIAL VELOCITY */
		res(2,i) = -res(1,i)*(-norm(1)*mvel(0,i) +norm(0)*mvel(1,i))/jcb*gbl->meshc(indx);
	}
	
	
	/* INTEGRATE & STORE MESH MOVEMENT RESIDUALS */
	basis::tri(x.log2p)->intgrtx1d(&lf(x.NV)(0),&res(0,0));
	basis::tri(x.log2p)->intgrt1d(&lf(x.NV+1)(0),&res(1,0));
	basis::tri(x.log2p)->intgrtx1d(&lf(x.NV+1)(0),&res(2,0));
	
	return;
}

int translating_surface::setup_preconditioner() {
	int indx,m,n,sind,v0,v1;
	TinyVector<FLT,tri_mesh::ND> nrm;
	FLT h, hsm;
	FLT dttang, dtnorm;
	FLT vslp;
	TinyMatrix<FLT,tri_mesh::ND,MXGP> crd, dcrd;
	TinyMatrix<FLT,4,MXGP> res;
	TinyMatrix<FLT,4,MXGP> lf;
	TinyVector<FLT,2> mvel;
	int last_phase, mp_phase;
    
	int err = hp_coupled_bdry::setup_preconditioner();
    if (!gbl->symmetric && !is_master) return(err);
	
	/**************************************************/
	/* DETERMINE MOVEMENT TIME STEP              */
	/**************************************************/
	gbl->vdt(0,Range::all(),Range::all()) = 0.0;
	
	for(indx=0; indx < base.nseg; ++indx) {
		sind = base.seg(indx);
		v0 = x.seg(sind).pnt(0);
		v1 = x.seg(sind).pnt(1);
		
		
#ifdef DETAILED_DT
		x.crdtocht1d(sind);
		for(n=0;n<tri_mesh::ND;++n)
			basis::tri(x.log2p)->proj1d(&x.cht(n,0),&crd(n,0),&dcrd(n,0));
		
		dtnorm = 1.0e99;
		dttang = 1.0e99;
		gbl->meshc(indx) = 1.0e99;
		for(i=0;i<basis::tri(x.log2p)->gpx();++i) {
			nrm(0) =  dcrd(1,i)*2;
			nrm(1) = -dcrd(0,i)*2;
			h = sqrt(nrm(0)*nrm(0) +nrm(1)*nrm(1));
			
			/* RELATIVE VELOCITY STORED IN MVEL(N)*/
			for(n=0;n<tri_mesh::ND;++n) {
				mvel(n) = vel(n) -(x.gbl->bd(0)*(crd(n,i) -dxdt(x.log2p,indx)(n,i)));
#ifdef MESH_REF_VEL
				mvel(n) -= x.gbl->mesh_ref_vel(n);
#endif
			}
			vslp = fabs(-mvel(0)*nrm(1)/h +mvel(1)*nrm(0)/h);
			hsm = h/(.25*(basis::tri(x.log2p)->p()+1)*(basis::tri(x.log2p)->p()+1));
			
			dttang = MIN(dttang,2.*ksprg(indx)*(.25*(basis::tri(x.log2p)->p()+1)*(basis::tri(x.log2p)->p()+1))/hsm);
			dtnorm = MIN(dtnorm,2.*vslp/hsm +x.gbl->bd(0));
			
			/* SET UP DISSIPATIVE COEFFICIENT */
			/* FOR UPWINDING LINEAR CONVECTIVE CASE SHOULD BE 1/|a| */
			/* RESIDUAL HAS DX/2 WEIGHTING */
			/* |a| dx/2 dv/dx  dx/2 dpsi */
			/* |a| dx/2 2/dx dv/dpsi  dpsi */
			/* |a| dv/dpsi  dpsi */
			// gbl->meshc(indx) = gbl->adis/(h*dtnorm*0.5);/* FAILED IN NATES UPSTREAM surface WAVE CASE */
			gbl->meshc(indx) = MIN(gbl->meshc(indx),gbl->adis/(h*(vslp/hsm +x.gbl->bd(0)))); /* FAILED IN MOVING UP TESTS */
		}
		nrm(0) =  (x.pnts(v1)(1) -x.pnts(v0)(1));
		nrm(1) = -(x.pnts(v1)(0) -x.pnts(v0)(0));
#else
		nrm(0) =  (x.pnts(v1)(1) -x.pnts(v0)(1));
		nrm(1) = -(x.pnts(v1)(0) -x.pnts(v0)(0));
		h = sqrt(nrm(0)*nrm(0) +nrm(1)*nrm(1));
		
		mvel(0) = vel(0)-(x.gbl->bd(0)*(x.pnts(v0)(0) -x.vrtxbd(1)(v0)(0)));
		mvel(1) = vel(1)-(x.gbl->bd(0)*(x.pnts(v0)(1) -x.vrtxbd(1)(v0)(1)));
#ifdef MESH_REF_VEL
		mvel(0) -= x.gbl->mesh_ref_vel(0);
		mvel(1) -= x.gbl->mesh_ref_vel(1);
#endif
		vslp = fabs(-mvel(0)*nrm(1)/h +mvel(1)*nrm(0)/h);
		
		mvel(0) = vel(0)-(x.gbl->bd(0)*(x.pnts(v1)(0) -x.vrtxbd(1)(v1)(0)));
		mvel(1) = vel(1)-(x.gbl->bd(0)*(x.pnts(v1)(1) -x.vrtxbd(1)(v1)(1)));
#ifdef MESH_REF_VEL
		mvel(0) -= x.gbl->mesh_ref_vel(0);
		mvel(1) -= x.gbl->mesh_ref_vel(1);
#endif
		vslp = MAX(vslp,fabs(-mvel(0)*nrm(1)/h +mvel(1)*nrm(0)/h));
		
		hsm = h/(.25*(basis::tri(x.log2p)->p()+1)*(basis::tri(x.log2p)->p()+1));
		
		dttang = 2.*ksprg(indx)*(.25*(basis::tri(x.log2p)->p()+1)*(basis::tri(x.log2p)->p()+1))/hsm;
		dtnorm = 2.*vslp/hsm +x.gbl->bd(0);
		
		/* SET UP DISSIPATIVE COEFFICIENT */
		/* FOR UPWINDING LINEAR CONVECTIVE CASE SHOULD BE 1/|a| */
		/* RESIDUAL HAS DX/2 WEIGHTING */
		/* |a| dx/2 dv/dx  dx/2 dpsi */
		/* |a| dx/2 2/dx dv/dpsi  dpsi */
		/* |a| dv/dpsi  dpsi */
		// gbl->meshc(indx) = gbl->adis/(h*dtnorm*0.5); /* FAILED IN NATES UPSTREAM surface WAVE CASE */
		gbl->meshc(indx) = gbl->adis/(h*(vslp/hsm +x.gbl->bd(0))); /* FAILED IN MOVING UP TESTS */
#endif
		
		dtnorm *= RAD(0.5*(x.pnts(v0)(0) +x.pnts(v1)(0)));
		
		nrm *= 0.5;
		
		gbl->vdt(indx,0,0) += -dttang*nrm(1)*basis::tri(x.log2p)->vdiag1d();
		gbl->vdt(indx,0,1) +=  dttang*nrm(0)*basis::tri(x.log2p)->vdiag1d();
		gbl->vdt(indx,1,0) +=  dtnorm*nrm(0)*basis::tri(x.log2p)->vdiag1d();
		gbl->vdt(indx,1,1) +=  dtnorm*nrm(1)*basis::tri(x.log2p)->vdiag1d();
		gbl->vdt(indx+1,0,0) = -dttang*nrm(1)*basis::tri(x.log2p)->vdiag1d();
		gbl->vdt(indx+1,0,1) =  dttang*nrm(0)*basis::tri(x.log2p)->vdiag1d();
		gbl->vdt(indx+1,1,0) =  dtnorm*nrm(0)*basis::tri(x.log2p)->vdiag1d();
		gbl->vdt(indx+1,1,1) =  dtnorm*nrm(1)*basis::tri(x.log2p)->vdiag1d();
		
		if (basis::tri(x.log2p)->sm()) {
			gbl->sdt(indx,0,0) = -dttang*nrm(1);
			gbl->sdt(indx,0,1) =  dttang*nrm(0);
			gbl->sdt(indx,1,0) =  dtnorm*nrm(0);
			gbl->sdt(indx,1,1) =  dtnorm*nrm(1);
			
#ifdef DETAILED_MINV
			int lsm = basis::tri(x.log2p)->sm();
			x.crdtocht1d(sind);
			for(n=0;n<tri_mesh::ND;++n)
				basis::tri(x.log2p)->proj1d(&x.cht(n,0),&crd(n,0),&dcrd(n,0));
			
			for(int m = 0; m<lsm; ++m) {
				for(int i=0;i<basis::tri(x.log2p)->gpx();++i) {
					nrm(0) =  dcrd(1,i);
					nrm(1) = -dcrd(0,i);
					res(0,i) = -dttang*nrm(1)*basis::tri(x.log2p)->gx(i,m+3);
					res(1,i) =  dttang*nrm(0)*basis::tri(x.log2p)->gx(i,m+3);
					res(2,i) =  dtnorm*nrm(0)*basis::tri(x.log2p)->gx(i,m+3);
					res(3,i) =  dtnorm*nrm(1)*basis::tri(x.log2p)->gx(i,m+3);
				}
				lf = 0;
				basis::tri(x.log2p)->intgrt1d(&lf(0,0),&res(0,0));
				basis::tri(x.log2p)->intgrt1d(&lf(1,0),&res(1,0));
				basis::tri(x.log2p)->intgrt1d(&lf(2,0),&res(2,0));
				basis::tri(x.log2p)->intgrt1d(&lf(3,0),&res(3,0));
				
				/* CFL = 0 WON'T WORK THIS WAY */
				lf(0) /= gbl->cfl(x.log2p,0);
				lf(1) /= gbl->cfl(x.log2p,0);
				lf(2) /= gbl->cfl(x.log2p,1);
				lf(3) /= gbl->cfl(x.log2p,1);
				
				for (n=0;n<lsm;++n) {
					gbl->ms(indx,2*m,2*n) = lf(0,n+2);
					gbl->ms(indx,2*m,2*n+1) = lf(1,n+2);
					gbl->ms(indx,2*m+1,2*n) = lf(2,n+2);
					gbl->ms(indx,2*m+1,2*n+1) = lf(3,n+2);
				}
				
				/* tang/norm, x/y,  mode,  vert */
				gbl->vms(indx,0,0,m,0) = lf(0,0);
				gbl->vms(indx,0,1,m,0) = lf(1,0);
				gbl->vms(indx,0,0,m,1) = lf(0,1);
				gbl->vms(indx,0,1,m,1) = lf(1,1);
				gbl->vms(indx,1,0,m,0) = lf(2,0);
				gbl->vms(indx,1,1,m,0) = lf(3,0);
				gbl->vms(indx,1,0,m,1) = lf(2,1);
				gbl->vms(indx,1,1,m,1) = lf(3,1);
			}
			
			int info;
			GETRF(2*lsm,2*lsm,&gbl->ms(indx,0,0),2*MAXP,&gbl->ipiv(indx,0),info);
			if (info != 0) {
				*x.gbl->log << "DGETRF FAILED IN SIDE MODE PRECONDITIONER\n";
                err = 1;
			}
			/*
			 \phi_n dx,dy*t = \phi_n Vt
			 \phi_t dx,dy*n = \phi_t Vn
			 */
#endif
		}
	}
	
	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		x.vbdry(base.vbdry(0))->vloadbuff(boundary::manifolds,&gbl->vdt(0,0,0),0,gbl->vdt.length(secondDim)*gbl->vdt.length(secondDim)-1,0);
		x.vbdry(base.vbdry(1))->vloadbuff(boundary::manifolds,&gbl->vdt(base.nseg,0,0),0,gbl->vdt.length(secondDim)*gbl->vdt.length(secondDim)-1,0);
		x.vbdry(base.vbdry(0))->comm_prepare(boundary::manifolds,mp_phase,boundary::symmetric);
		x.vbdry(base.vbdry(1))->comm_prepare(boundary::manifolds,mp_phase,boundary::symmetric);
		
		x.vbdry(base.vbdry(0))->comm_exchange(boundary::manifolds,mp_phase,boundary::symmetric);
		x.vbdry(base.vbdry(1))->comm_exchange(boundary::manifolds,mp_phase,boundary::symmetric);
		
		last_phase = true;
		last_phase &= x.vbdry(base.vbdry(0))->comm_wait(boundary::manifolds,mp_phase,boundary::symmetric);
		last_phase &= x.vbdry(base.vbdry(1))->comm_wait(boundary::manifolds,mp_phase,boundary::symmetric);
		x.vbdry(base.vbdry(0))->vfinalrcv(boundary::manifolds,mp_phase,boundary::symmetric,boundary::average,&gbl->vdt(0,0,0),0,gbl->vdt.length(secondDim)*gbl->vdt.length(secondDim)-1,0);
		x.vbdry(base.vbdry(1))->vfinalrcv(boundary::manifolds,mp_phase,boundary::symmetric,boundary::average,&gbl->vdt(base.nseg,0,0),0,gbl->vdt.length(secondDim)*gbl->vdt.length(secondDim)-1,0);
	}
	
	if (gbl->is_loop) {
		for(m=0;m<tri_mesh::ND;++m)
			for(n=0;n<tri_mesh::ND;++n)
				gbl->vdt(0,m,n) = 0.5*(gbl->vdt(0,m,n) +gbl->vdt(base.nseg+1,m,n));
		gbl->vdt(base.nseg+1,Range::all(),Range::all()) = gbl->vdt(0,Range::all(),Range::all());
	}
	
	FLT jcbi,temp;
	for(indx=0;indx<base.nseg+1;++indx) {
		/* INVERT VERTEX MATRIX */
		jcbi = 1.0/(gbl->vdt(indx,0,0)*gbl->vdt(indx,1,1) -gbl->vdt(indx,0,1)*gbl->vdt(indx,1,0));
		
		temp = gbl->vdt(indx,0,0)*jcbi*gbl->cfl(x.log2p,1);
		gbl->vdt(indx,0,0) = gbl->vdt(indx,1,1)*jcbi*gbl->cfl(x.log2p,0);
		gbl->vdt(indx,1,1) = temp;
		gbl->vdt(indx,0,1) *= -jcbi*gbl->cfl(x.log2p,1);
		gbl->vdt(indx,1,0) *= -jcbi*gbl->cfl(x.log2p,0);

#ifdef TEMPO
		gbl->vdt(indx,Range::all(),Range::all()) = 0.0;
		gbl->vdt(indx,0,0) = 1.0;
		gbl->vdt(indx,1,1) = 1.0;
		gbl->vdt(indx+1,Range::all(),Range::all()) = 0.0;
		gbl->vdt(indx+1,0,0) = 1.0;
		gbl->vdt(indx+1,1,1) = 1.0;
#endif
		/* DIRECT FORMATION OF vdt^{-1} theta is angle of normal from horizontal */
		//		FLT theta =  100.0*M_PI/180.0;
		//		gbl->vdt(indx,0,0) = -sin(theta);
		//		gbl->vdt(indx,1,1) =  sin(theta);
		//		gbl->vdt(indx,0,1) = cos(theta);
		//		gbl->vdt(indx,1,0) = cos(theta);
	}
	
	/* INVERT SIDE MATRIX */
	if (basis::tri(x.log2p)->sm() > 0) {
		for(indx=0;indx<base.nseg;++indx) {
			/* INVERT SIDE MVDT MATRIX */
			jcbi = 1.0/(gbl->sdt(indx,0,0)*gbl->sdt(indx,1,1) -gbl->sdt(indx,0,1)*gbl->sdt(indx,1,0));
			
			temp = gbl->sdt(indx,0,0)*jcbi*gbl->cfl(x.log2p,1);
			gbl->sdt(indx,0,0) = gbl->sdt(indx,1,1)*jcbi*gbl->cfl(x.log2p,0);
			gbl->sdt(indx,1,1) = temp;
			gbl->sdt(indx,0,1) *= -jcbi*gbl->cfl(x.log2p,1);
			gbl->sdt(indx,1,0) *= -jcbi*gbl->cfl(x.log2p,0);
		}
	}
    return(err);
}

/* Routine to make sure r_gbl residual doesn't get screwed up at triple junction */
void hp_deformable_follower_pnt::rsdl(int stage) {
#ifdef petsc
	r_tri_mesh::global *r_gbl = dynamic_cast<r_tri_mesh::global *>(x.gbl);
	/* equation for tangential poistion */
	r_gbl->res(base.pnt)(0) = 0.0;
	r_gbl->res(base.pnt)(1) = 0.0;
#endif
	
	if (!surf->is_master) return;
		
	int endpt;
	if (surfbdry == 0)
		endpt = x.ebdry(base.ebdry(0))->nseg;
	else
		endpt = 0;
	
	// FIXME: This should be deleted eventually
	surf->gbl->vres(endpt,0) = 0.0;
	surf->gbl->vres(endpt,1) = 0.0;
}

#ifdef petsc
void hp_deformable_follower_pnt::petsc_jacobian() {
	/* ZERO TANGENT AND NORMAL MESH MOVEMENT ROW BEFORE COMMUNICATION */
	int row = (x.NV+tri_mesh::ND)*base.pnt +x.NV;
#ifdef MY_SPARSE
	x.J.zero_row(row);
	x.J_mpi.zero_row(row++);
	x.J.zero_row(row);
	x.J_mpi.zero_row(row++);
#else
	/* Must zero rows of jacobian created by r_mesh */
	MatAssemblyBegin(x.petsc_J,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(x.petsc_J,MAT_FINAL_ASSEMBLY);
	MatZeroRow(x.petsc_J,row,PETSC_NULL);
	MatZeroRow(x.petsc_J,row+1,PETSC_NULL);
#endif
	
	multi_physics_pnt::petsc_jacobian();
}
#endif



