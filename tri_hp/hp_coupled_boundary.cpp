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
	
	if (!inmap.get(base.idprefix + "_nvariable",NV)) {
		*x.gbl->log << "Couldn't find number of variables for coupled boundry" << std::endl;
		sim::abort(__LINE__,__FILE__,x.gbl->log);
	}
	
	keyword = base.idprefix + "_coupled";
	inmap[keyword] = "1";
	
	hp_edge_bdry::init(inmap,gbl_in);
	
	
	if (!is_master) return;   /* This is all that is necessary for a slave boundary */
	
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
	vdres.resize(x.log2pmax,base.maxseg,NV);
	sdres.resize(x.log2pmax,base.maxseg,x.sm0,NV);
	
	if (x.seg(base.seg(0)).pnt(0) == x.seg(base.seg(base.nseg-1)).pnt(1)) gbl->is_loop = true;
	else gbl->is_loop = false;
	
	gbl->vug0.resize(base.maxseg+1,NV);
	gbl->sug0.resize(base.maxseg,x.sm0,NV);
	
	gbl->vres.resize(base.maxseg+1,NV);
	gbl->sres.resize(base.maxseg,x.sm0,NV);
	gbl->vres0.resize(base.maxseg+1,NV);
	gbl->sres0.resize(base.maxseg,x.sm0,NV);
	
	gbl->vdt.resize(base.maxseg+1,NV,NV);
	gbl->sdt.resize(base.maxseg,NV,NV);
	gbl->meshc.resize(base.maxseg,NV);
	
#ifdef DETAILED_MINV
	const int sm0 = x.sm0;
	gbl->ms.resize(base.maxseg,NV*sm0,NV*sm0);
	gbl->vms.resize(base.maxseg,NV,2,sm0,2);
	gbl->ipiv.resize(base.maxseg,NV*sm0);
#endif
		
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
	
	gbl->cfl.resize(x.log2p+1,NV);
	double CFLdflt[3] = {2.5, 1.5, 1.0};
	for (int n=0;n<NV;++n) {
		stringstream nstr;
		nstr << n;
		if (inmap.getline(base.idprefix+"_cfl"+nstr.str(),val)) {
			data.str(val);
			for (int m=0;m<x.log2p+1;++m) {
				data >> gbl->cfl(m,n);
			}
			data.clear();
		}
		else {
			for (int m=0;m<x.log2p+1;++m) {
				gbl->cfl(m,n) = CFLdflt[m];
			}
		}
	}
	
	return;
}

void hp_coupled_bdry::rsdl(int stage) {
    
    if (!is_master) return;
    
    int i,m,n,sind,indx,v0,v1;
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
    
    /* CALL VERTEX RESIDUAL HERE */
    for(i=0;i<2;++i)
        x.hp_vbdry(base.vbdry(i))->rsdl(stage);
    
    /************************************************/
    /* MODIFY SURFACE RESIDUALS ON COARSER MESHES    */
    /************************************************/
    if(x.coarse_flag) {
        if (x.isfrst) {
            for(i=0;i<base.nseg+1;++i)
                for(n=0;n<NV;++n)
                    vdres(x.log2p,i,n) = gbl->fadd(n)*gbl->vres0(i,n) -gbl->vres(i,n);
            
            for(i=0;i<base.nseg;++i)
                for(m=0;m<basis::tri(x.log2p)->sm();++m)
                    for(n=0;n<NV;++n)
                        sdres(x.log2p,i,m,n) = gbl->fadd(n)*gbl->sres0(i,m,n) -gbl->sres(i,m,n);
            
        }
        for(i=0;i<base.nseg+1;++i)
            for(n=0;n<NV;++n)
                gbl->vres(i,n) += vdres(x.log2p,i,n);
        
        for(i=0;i<base.nseg;++i)
            for(m=0;m<basis::tri(x.log2p)->sm();++m)
                for(n=0;n<NV;++n)
                    gbl->sres(i,m,n) += sdres(x.log2p,i,m,n);
    }
}



void hp_coupled_bdry::minvrt() {
    int i,m,n,indx;
    int last_phase, mp_phase;
    FLT temp;
    
    /* INVERT MASS MATRIX */
    /* LOOP THROUGH SIDES */
    if (basis::tri(x.log2p)->sm() > 0) {
        for(indx = 0; indx<base.nseg; ++indx) {
            /* SUBTRACT SIDE CONTRIBUTIONS TO VERTICES */
            for (m=0; m <basis::tri(x.log2p)->sm(); ++m) {
                for(n=0;n<NV;++n)
                    gbl->vres(indx,n) -= basis::tri(x.log2p)->sfmv1d(0,m)*gbl->sres(indx,m,n);
                for(n=0;n<NV;++n)
                    gbl->vres(indx+1,n) -= basis::tri(x.log2p)->sfmv1d(1,m)*gbl->sres(indx,m,n);
            }
        }
    }
    for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
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
            gbl->vres(0,n) = 0.5*(gbl->vres(0,n) +gbl->vres(base.nseg+1,n));
            gbl->vres(base.nseg+1,n) = gbl->vres(0,n);
        }
        //		TEMPO: Going to need to make a loop boundary condition vertex to make this general
        //		gbl->vres(0)(0) = 0.0;
        //		gbl->vres(base.nseg+1)(0) = 0.0;
    }
    x.hp_vbdry(base.vbdry(0))->vdirichlet();
    x.hp_vbdry(base.vbdry(1))->vdirichlet();
    
    
    /* SOLVE FOR VERTEX MODES */
    /* TEMPO: THIS IS NOT GENERAL ONLY FOR DEFORMABLE RIGHT NOW */
    for(i=0;i<base.nseg+1;++i) {
        temp                     = gbl->vres(i,0)*gbl->vdt(i,0,0) +gbl->vres(i,1)*gbl->vdt(i,0,1);
        gbl->vres(i,1) = gbl->vres(i,0)*gbl->vdt(i,1,0) +gbl->vres(i,1)*gbl->vdt(i,1,1);
        gbl->vres(i,0) = temp;
    }
    
    /* SOLVE FOR SIDE MODES */
    if (basis::tri(x.log2p)->sm() > 0) {
        for(indx = 0; indx<base.nseg; ++indx) {
            
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
            for(m=0;m<basis::tri(x.log2p)->sm();++m) {
                temp                             = gbl->sres(indx,m,0)*gbl->sdt(indx,0,0) +gbl->sres(indx,m,1)*gbl->sdt(indx,0,1);
                gbl->sres(indx,m,1) = gbl->sres(indx,m,0)*gbl->sdt(indx,1,0) +gbl->sres(indx,m,1)*gbl->sdt(indx,1,1);
                gbl->sres(indx,m,0) = temp;
            }
            
            for(m=0;m<basis::tri(x.log2p)->sm();++m) {
                for(n=0;n<NV;++n) {
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
	
	if (!is_master) return;
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


void hp_coupled_bdry::maxres() {
	if (!is_master) return;
	
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

#ifdef petsc
void hp_deformable_bdry::non_sparse(Array<int,1> &nnzero) {
	const int sm=basis::tri(x.log2p)->sm();
	const int im=basis::tri(x.log2p)->im();
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	const int begin_seg = x.npnt*vdofs;
	const int begin_tri = begin_seg+x.nseg*sm*x.NV;
	
	if(x.sm0 > 0) {
		if (is_master) {
			nnzero(Range(jacobian_start,jacobian_start+base.nseg*sm*tri_mesh::ND-1)) = vdofs*(sm+2);
		}
		else {
			/* Just an equality constraint */
			nnzero(Range(jacobian_start,jacobian_start+base.nseg*sm*tri_mesh::ND-1)) = 1;
		}
		
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
			
			int pind = x.seg(sind).pnt(0);
			if (is_master) {
				nnzero(Range(pind*vdofs,(pind+1)*vdofs-1)) += NV*sm;
			}
			else {
				nnzero(Range(pind*vdofs,pind*vdofs +x.NV-1)) += NV*sm;
			}
			
			pind = x.seg(sind).pnt(1);
			if (is_master) {
				nnzero(Range(pind*vdofs,(pind+1)*vdofs-1)) += NV*sm;
			}
			else {
				nnzero(Range(pind*vdofs,pind*vdofs +x.NV-1)) += NV*sm;
			}
		}
	}
}
#endif

void hp_deformable_bdry::init(input_map& inmap,void* gbl_in) {
    std::string keyword;
    
    if (!inmap.get(base.idprefix + "_nvariable",NV)) {
        NV = 2; // Default for deformable boundary (x,y positions)
        inmap[base.idprefix +"_nvariable"] = "2";
    }
    
    keyword = base.idprefix + "_curved";
    inmap[keyword] = "1";
	
		hp_coupled_bdry::init(inmap,gbl_in);
  
    if (is_master) ksprg.resize(base.maxseg);
}

void hp_deformable_bdry::tadvance() {
	int i,j,m,sind;
	
	hp_coupled_bdry::tadvance();
	
    if (x.gbl->substep == 0) {
        if (is_master) {
            /* SET SPRING CONSTANTS */
            for(j=0;j<base.nseg;++j) {
                sind = base.seg(j);
                ksprg(j) = 1.0/x.distance(x.seg(sind).pnt(0),x.seg(sind).pnt(1));
            }
            
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
	return;
}

void hp_deformable_bdry::rsdl(int stage) {
	
	hp_coupled_bdry::rsdl(stage);

    if (is_master) {
        if (!x.coarse_flag) {
            /* ADD TANGENTIAL MESH MOVEMENT SOURCE */
            for(int i=0;i<base.nseg+1;++i)
                gbl->vres(i,0) += vdres(x.log2p,i,0);
            
            for(int i=0;i<base.nseg;++i)
                for(int m=0;m<basis::tri(x.log2p)->sm();++m)
                    gbl->sres(i,m,0) += sdres(x.log2p,i,m,0);
        }
            
#ifdef petsc
        r_tri_mesh::global *r_gbl = dynamic_cast<r_tri_mesh::global *>(x.gbl);
        int i = 0;
        int sind, v0;
        /* Store rotated vertex residual in r_mesh residual vector */
        do {
            sind = base.seg(i);
            v0 = x.seg(sind).pnt(0);
            /* Rotate residual for better diagonal dominance */
            r_gbl->res(v0)(0) = gbl->vres(i,0)*gbl->vdt(i,0,0) +gbl->vres(i,1)*gbl->vdt(i,0,1);
            r_gbl->res(v0)(1) = gbl->vres(i,0)*gbl->vdt(i,1,0) +gbl->vres(i,1)*gbl->vdt(i,1,1);
        } while (++i < base.nseg);
        v0 = x.seg(sind).pnt(1);
        r_gbl->res(v0)(0) = gbl->vres(i,0)*gbl->vdt(i,0,0) +gbl->vres(i,1)*gbl->vdt(i,0,1);
        r_gbl->res(v0)(1) = gbl->vres(i,0)*gbl->vdt(i,1,0) +gbl->vres(i,1)*gbl->vdt(i,1,1);
#endif
    }

#ifndef petsc
    r_tri_mesh::global *r_gbl = dynamic_cast<r_tri_mesh::global *>(x.gbl);
    int i = 0;
    int sind, v0;
    /* Store 0.0 in vertex residual in r_mesh residual vector */
    do {
        sind = base.seg(i);
        v0 = x.seg(sind).pnt(0);
        r_gbl->res(v0)(0) = 0.0;
        r_gbl->res(v0)(1) = 0.0;
    } while (++i < base.nseg);
    v0 = x.seg(sind).pnt(1);
    r_gbl->res(v0)(0) = 0.0;
    r_gbl->res(v0)(1) = 0.0;
#endif
}

#ifdef petsc
void hp_deformable_bdry::petsc_make_1D_rsdl_vector(Array<double,1> res) {
	int sm = basis::tri(x.log2p)->sm();
	int ind = jacobian_start;
	
    if (is_master) {
        /* Side residuals rotated for better diagonal dominance */
        for(int j=0;j<base.nseg;++j) {
            for(int m=0;m<sm;++m) {
                res(ind++) = gbl->sres(j,m,0)*gbl->sdt(j,0,0) +gbl->sres(j,m,1)*gbl->sdt(j,0,1);
                res(ind++) = gbl->sres(j,m,0)*gbl->sdt(j,1,0) +gbl->sres(j,m,1)*gbl->sdt(j,1,1);
            }
        }
    }
    else {
        /* Curvature equality constraint */
        for(int j=0;j<base.nseg;++j) {
            for(int m=0;m<sm;++m) {
                res(ind++) = 0.0;
                res(ind++) = 0.0;
            }
        }
    }
}
#endif

void hp_deformable_bdry::update(int stage) {
    
    if (is_master) {
        int i,m,n,count,sind,indx,v0;
        
        if (stage < 0) {
            indx = 0;
            int i = 0;
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
                    for(m=0;m<basis::tri(x.log2p)->sm();++m)
                        for(int n=0;n<tri_mesh::ND;++n)
                            gbl->sug0(i,m,n) = crv(i,m)(n);
            }
            
            return;
        }
        
        minvrt();
        
#ifdef DEBUG
        //   if (x.coarse_flag) {
        for(i=0;i<base.nseg+1;++i)
					*x.gbl->log << "vdt: " << i << ' ' << gbl->vdt(i,0,0) << ' ' << gbl->vdt(i,0,1) << ' ' << gbl->vdt(i,1,0) << ' ' << gbl->vdt(i,1,1) << '\n';
        
        for(i=0;i<base.nseg;++i)
            *x.gbl->log << "sdt: " << i << ' ' << gbl->sdt(i,0,0) << ' ' << gbl->sdt(i,0,1) << ' ' << gbl->sdt(i,1,0) << ' ' << gbl->sdt(i,1,1) << '\n';
        
        for(i=0;i<base.nseg+1;++i) {
            *x.gbl->log << "vres: " << i << ' ';
            for(n=0;n<tri_mesh::ND;++n) {
                if (fabs(gbl->vres(i,n)) > 1.0e-9) *x.gbl->log << gbl->vres(i,n) << ' ';
                else *x.gbl->log << "0.0 ";
            }
            *x.gbl->log << '\n';
        }
        
        for(i=0;i<base.nseg;++i) {
            for(m=0;m<basis::tri(x.log2p)->sm();++m) {
                *x.gbl->log << "sres: " << i << ' ';
                for(n=0;n<tri_mesh::ND;++n) {
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
        
        for(i=0;i<base.nseg;++i)
            for(m=0;m<basis::tri(x.log2p)->sm();++m)
                *x.gbl->log << "spos: " << i << ' ' << m << ' ' << crv(i,m)(0) << ' ' << crv(i,m)(1) << '\n';
        // }
#endif
        
        i = 0;
        do {
            sind = base.seg(i);
            v0 = x.seg(sind).pnt(0);
            for (int n=0;n<NV;++n)
                x.pnts(v0)(n) = gbl->vug0(i,n) -x.gbl->alpha(stage)*gbl->vres(i,n);
        } while (++i < base.nseg);
        v0 = x.seg(sind).pnt(1);
        for (int n=0;n<NV;++n)
            x.pnts(v0)(n) = gbl->vug0(base.nseg,n) -x.gbl->alpha(stage)*gbl->vres(base.nseg,n);
        
        if (basis::tri(x.log2p)->sm() > 0) {
            for(int i=0;i<base.nseg;++i)
                for(int m=0;m<basis::tri(x.log2p)->sm();++m)
                    for(int n=0;n<tri_mesh::ND;++n)
                        crv(i,m)(n) = gbl->sug0(i,m,n) -x.gbl->alpha(stage)*gbl->sres(i,m,n);
        }
        
        /* FIX POINTS THAT SLIDE ON CURVE */
        x.hp_vbdry(base.vbdry(1))->mvpttobdry(x.pnts(v0));
        sind = base.seg(0);
        v0 = x.seg(sind).pnt(0);
        x.hp_vbdry(base.vbdry(0))->mvpttobdry(x.pnts(v0));
        
#ifdef DEBUG
        //   if (x.coarse_flag) {
        i = 0;
        do {
            sind = base.seg(i);
            v0 = x.seg(sind).pnt(0);
            *x.gbl->log << "vertex positions " << v0 << ' ' << x.pnts(v0)(0) << ' ' << x.pnts(v0)(1) << '\n';
        } while(++i < base.nseg);
        v0 = x.seg(sind).pnt(1);
        *x.gbl->log << "vertex positions " << v0 << ' ' << x.pnts(v0)(0) << ' ' << x.pnts(v0)(1) << '\n';
        
        for(i=0;i<base.nseg;++i)
            for(m=0;m<basis::tri(x.log2p)->sm();++m)
                *x.gbl->log << "spos: " << i << ' ' << m << ' ' << crv(i,m)(0) << ' ' << crv(i,m)(1) << '\n';
        // }
#endif
        
        if (base.is_comm()) {
            count = 0;
            i = 0;
            do {
                sind = base.seg(i);
                v0 = x.seg(sind).pnt(0);
                for(n=0;n<tri_mesh::ND;++n)
                    base.fsndbuf(count++) = x.pnts(v0)(n);
            } while(++i < base.nseg);
            v0 = x.seg(sind).pnt(1);
            for(n=0;n<tri_mesh::ND;++n)
                base.fsndbuf(count++) = x.pnts(v0)(n);
            
            for(i=0;i<base.nseg;++i) {
                for(m=0;m<basis::tri(x.log2p)->sm();++m) {
                    for(n=0;n<tri_mesh::ND;++n)
                        base.fsndbuf(count++) = crv(i,m)(n);
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
    }
    
	return;
}

#ifdef petsc
void hp_deformable_bdry::element_jacobian(int indx, Array<FLT,2>& K) {
	int sm = basis::tri(x.log2p)->sm();

	Array<TinyVector<FLT,MXTM>,1> Rbar(x.NV+tri_mesh::ND),lf(x.NV+tri_mesh::ND);
	
#ifdef DEBUG_JAC
	const FLT eps_r = 0.0e-6, eps_a = 1.0e-6;  /*<< constants for debugging jacobians */
#else
	const FLT eps_r = 1.0e-6, eps_a = 1.0e-10;  /*<< constants for accurate numerical determination of jacobians */
#endif
	
	/* Calculate and store initial residual */
	int sind = base.seg(indx);
	x.ugtouht1d(sind);
	element_rsdl(indx,Rbar);
	
	Array<FLT,1> dw(x.NV);
	dw = 0.0;
	for(int i=0;i<2;++i)
		for(int n=0;n<x.NV;++n)
			dw = dw + fabs(x.uht(n)(i));
	
	dw = dw*eps_r;
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



void hp_deformable_bdry::petsc_jacobian() {
    
    const int sm = basis::tri(x.log2p)->sm();
		const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;

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
    
    if (is_master) {
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
            
            /* Rotate for diagonal dominance */
            /* Rotate vertices */
            ind = x.NV;
            row_store             = K(ind,Range::all())*gbl->vdt(j,0,0) +K(ind+1,Range::all())*gbl->vdt(j,0,1);
            K(ind+1,Range::all()) = K(ind,Range::all())*gbl->vdt(j,1,0) +K(ind+1,Range::all())*gbl->vdt(j,1,1);
            K(ind,Range::all()) = row_store;
            ind += vdofs;
            
            row_store             = K(ind,Range::all())*gbl->vdt(j+1,0,0) +K(ind+1,Range::all())*gbl->vdt(j+1,0,1);
            K(ind+1,Range::all()) = K(ind,Range::all())*gbl->vdt(j+1,1,0) +K(ind+1,Range::all())*gbl->vdt(j+1,1,1);
            K(ind,Range::all()) = row_store;
            ind += vdofs;
            
            /* SAME FOR SIDE MODES */
            for(int m=0;m<basis::tri(x.log2p)->sm();++m) {
                row_store             = K(ind,Range::all())*gbl->sdt(j,0,0) +K(ind+1,Range::all())*gbl->sdt(j,0,1);
                K(ind+1,Range::all()) = K(ind,Range::all())*gbl->sdt(j,1,0) +K(ind+1,Range::all())*gbl->sdt(j,1,1);
                K(ind,Range::all()) = row_store;
                ind += vdofs;
            }
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
#ifdef DEBUG_JAC
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
    
    if (is_master) {
        /* FIXME: NOT SURE ABOUT THIS TEMPORARY */
        x.hp_vbdry(base.vbdry(0))->petsc_jacobian();
        x.hp_vbdry(base.vbdry(1))->petsc_jacobian();
        x.hp_vbdry(base.vbdry(0))->petsc_jacobian_dirichlet();
        x.hp_vbdry(base.vbdry(1))->petsc_jacobian_dirichlet();
    }
}

void hp_deformable_bdry::non_sparse_rcv(Array<int, 1> &nnzero, Array<int, 1> &nnzero_mpi) {
	hp_coupled_bdry::non_sparse_rcv(nnzero,nnzero_mpi);
	if (!is_master) nnzero_mpi(Range(jacobian_start,jacobian_start+base.nseg*x.sm0*tri_mesh::ND-1)) = 1;  // equality constraint
}

void hp_deformable_bdry::petsc_matchjacobian_rcv(int phase) {
	
	hp_coupled_bdry::petsc_matchjacobian_rcv(phase);
	
	const int sm = x.sm0;
	
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

		/* Equality of curved mode constraint */
		int ind = jacobian_start;
		int ind_mpi = static_cast<int>(base.frcvbuf(0,1)) +(base.nseg-1)*sm*x.ND +Jstart_mpi;
		*x.gbl->log << ind_mpi << ' ' << base.sndsize() << std::endl;
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

#endif

