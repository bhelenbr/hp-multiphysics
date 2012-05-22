#include "bdry_buoyancy.h"
#include <myblas.h>
#include <blitz/tinyvec-et.h>

//#define MPDEBUG

//#define DEBUG

using namespace bdry_buoyancy;

void melt::init(input_map& inmap,void* gbl_in) {
	std::string keyword,val;
	std::istringstream data;
	std::string filename;
	
	keyword = base.idprefix + "_curved";
	inmap[keyword] = "1";
	
	keyword = base.idprefix + "_coupled";
	inmap[keyword] = "1";
	
	/* Load in the heat flux */
	keyword = base.idprefix + "_hp_typelist";
	inmap[keyword] = "0 0 1 0";
	symbolic::init(inmap,gbl_in);
	
	/* Now that we have loaded flux, make temperature a dirichlet bc */
	/* And make pressure neumann */
	essential_indices.clear();
	for(int n=0;n<x.NV-1;++n)
		essential_indices.push_back(n);
	type(Range(0,x.NV-2)) = essential;
	type(x.NV-1) = natural;
	
	gbl = static_cast<global *>(gbl_in);
	ksprg.resize(base.maxseg);
	vdres.resize(x.log2pmax,base.maxseg);
	sdres.resize(x.log2pmax,base.maxseg,x.sm0);
	
	keyword = base.idprefix + "_Lf";
	if (!inmap.get(keyword,gbl->Lf)) {
		*x.gbl->log << "Missing latent heat of fusion " << keyword << std::endl;
		sim::abort(__LINE__,__FILE__,x.gbl->log);
	}
	
	if (!base.is_comm()) {
		keyword = base.idprefix + "_cp_s";
		if (!inmap.get(keyword,gbl->cp_s)) {
			*x.gbl->log << "Missing c_p of solid " << keyword << std::endl;
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}
		
		keyword = base.idprefix + "_rho_s";
		if (!inmap.get(keyword,gbl->rho_s)) {
			*x.gbl->log << "Missing rho of solid " << keyword << std::endl;
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}
	}
	else {
		gbl->cp_s = 0.0;
		gbl->rho_s = 0.0;
	}
	
	if (x.seg(base.seg(0)).pnt(0) == x.seg(base.seg(base.nseg-1)).pnt(1)) gbl->is_loop = true;
	else gbl->is_loop = false;
	
	gbl->vug0.resize(base.maxseg+1);
	gbl->sug0.resize(base.maxseg,x.sm0);
	
	gbl->vres.resize(base.maxseg+1);
	gbl->sres.resize(base.maxseg,x.sm0); 
	gbl->vres0.resize(base.maxseg+1);
	gbl->sres0.resize(base.maxseg,x.sm0); 
	
	gbl->vdt.resize(base.maxseg+1);
	gbl->sdt.resize(base.maxseg);      
	gbl->meshc.resize(base.maxseg);
	
#ifdef DETAILED_MINV
	gbl->ms.resize(base.maxseg); 
	gbl->vms.resize(base.maxseg,2,2,MAXP,2);
	gbl->ipiv.resize(base.maxseg); 
#endif
	
	/* Multigrid Storage all except highest order (log2p+1)*/
	vdres.resize(x.log2p+1,base.maxseg+1);
	sdres.resize(x.log2p+1,base.maxseg,x.sm0);
	
	keyword = base.idprefix + "_fadd";
	inmap.getlinewdefault(keyword,val,"1.0 1.0");
	data.str(val);
	data >> gbl->fadd(0) >> gbl->fadd(1);  
	data.clear(); 
	
	double CFLtdflt[3] = {2.5, 1.5, 1.0};
	inmap.getwdefault(base.idprefix + "_cfltangent",&gbl->cfl(0,0),3,CFLtdflt); 
	
	double CFLndflt[3] = {2.0, 1.25, 0.75};
	inmap.getwdefault(base.idprefix + "_cflnormal",&gbl->cfl(1,0),3,CFLndflt); 
	
	inmap.getwdefault(base.idprefix +"_adis",gbl->adis,1.0);
	
	return;
}

void melt::tadvance() {
	int i,j,m,sind;    

	symbolic::tadvance();

	if (x.gbl->substep == 0) {
		/* SET SPRING CONSTANTS */
		for(j=0;j<base.nseg;++j) {
			sind = base.seg(j);
			ksprg(j) = 1.0/x.distance(x.seg(sind).pnt(0),x.seg(sind).pnt(1));
		}

		/* CALCULATE TANGENT SOURCE TERM FOR FINE MESH */
		/* ZERO TANGENTIAL MESH MOVEMENT SOURCE */    
		if (!x.coarse_level) {
			for(i=0;i<base.nseg+1;++i)
				vdres(x.log2p,i)(0) = 0.0;

			for(i=0;i<base.nseg;++i) 
				for(m=0;m<basis::tri(x.log2p)->sm();++m)
					sdres(x.log2p,i,m)(0) = 0.0;

			rsdl(x.gbl->nstage);

			for(i=0;i<base.nseg+1;++i)
				vdres(x.log2p,i)(0) = -gbl->vres(i)(0);

			for(i=0;i<base.nseg;++i) 
				for(m=0;m<basis::tri(x.log2p)->sm();++m)
					sdres(x.log2p,i,m)(0) = -gbl->sres(i,m)(0)*0;  /* TO KEEP SIDE MODES EQUALLY SPACED */
		}
	}
	return;
}

void melt::element_rsdl(int indx, Array<TinyVector<FLT,MXTM>,1> lf) {
	int i,n,sind,v0,v1;
	TinyVector<FLT,tri_mesh::ND> norm, rp;
	Array<FLT,1> ubar(x.NV);
	FLT jcb;
	Array<TinyVector<FLT,MXGP>,1> u(x.NV);
	TinyMatrix<FLT,tri_mesh::ND,MXGP> crd, dcrd, mvel;
	TinyMatrix<FLT,8,MXGP> res;
	
	sind = base.seg(indx);
	v0 = x.seg(sind).pnt(0);
	v1 = x.seg(sind).pnt(1);
	
	x.crdtocht1d(sind);
	for(n=0;n<tri_mesh::ND;++n)
		basis::tri(x.log2p)->proj1d(&x.cht(n,0),&crd(n,0),&dcrd(n,0));
	
	for(n=0;n<x.NV;++n)
		basis::tri(x.log2p)->proj1d(&x.uht(n)(0),&u(n)(0));    
	
	for(i=0;i<basis::tri(x.log2p)->gpx();++i) {
		norm(0) =  dcrd(1,i);
		norm(1) = -dcrd(0,i);
		jcb = sqrt(norm(0)*norm(0) +norm(1)*norm(1));
		
		/* RELATIVE VELOCITY STORED IN MVEL(N)*/
		for(n=0;n<tri_mesh::ND;++n) {
			mvel(n,i) = u(n)(i) -(x.gbl->bd(0)*(crd(n,i) -dxdt(x.log2p,indx)(n,i)));  
#ifdef MESH_REF_VEL
			mvel(n,i) -= x.gbl->mesh_ref_vel(n);
#endif
		}
		
		/* TANGENTIAL SPACING */                
		res(0,i) = -ksprg(indx)*jcb;
		/* NORMAL FLUX */
		res(1,i) = RAD(crd(0,i))*x.gbl->rho*(mvel(0,i)*norm(0) +mvel(1,i)*norm(1));     
		/* UPWINDING BASED ON TANGENTIAL VELOCITY NO!!! */
		// res(2,i) = -res(1,i)*(-norm(1)*mvel(0,i) +norm(0)*mvel(1,i))/jcb*gbl->meshc(indx);
		/* Heat Flux */
		Array<FLT,1> au(x.NV), axpt(tri_mesh::ND), amv(tri_mesh::ND), anorm(tri_mesh::ND);
		for(n=0;n<x.NV;++n)
			au(n) = u(n)(i);
		axpt(0) = crd(0,i); axpt(1) = crd(1,i);
		amv(0) = (x.gbl->bd(0)*(crd(0,i) -dxdt(x.log2p,indx)(0,i))); amv(1) = (x.gbl->bd(0)*(crd(1,i) -dxdt(x.log2p,indx)(1,i)));
#ifdef MESH_REF_VEL
		amv(0) += x.gbl->mesh_ref_vel(0);
		amv(1) += x.gbl->mesh_ref_vel(1);
#endif
		anorm(0)= norm(0)/jcb; anorm(1) = norm(1)/jcb;
		res(3,i) = RAD(crd(0,i))*fluxes(2).Eval(au,axpt,amv,anorm,x.gbl->time)*jcb -gbl->Lf*res(1,i) +gbl->rho_s*gbl->cp_s*u(2)(i)*res(1,i)/x.gbl->rho;
		
		/* Heat Flux Upwinded NO!!! */
		// res(4,i) = -res(3,i)*(-norm(1)*mvel(0,i) +norm(0)*mvel(1,i))/jcb*gbl->meshc(indx);
	}
	lf = 0.0;
	
	/* INTEGRATE & STORE MESH MOVEMENT RESIDUALS */                    
	basis::tri(x.log2p)->intgrtx1d(&lf(x.NV)(0),&res(0,0)); // tangent
	basis::tri(x.log2p)->intgrt1d(&lf(x.NV-1)(0),&res(1,0)); // mass flux
	//basis::tri(x.log2p)->intgrtx1d(&lf(x.NV-1)(0),&res(2,0)); // mass flux
	basis::tri(x.log2p)->intgrt1d(&lf(x.NV-2)(0),&res(3,0)); // heat flux
	//basis::tri(x.log2p)->intgrtx1d(&lf(x.NV-2)(0),&res(4,0)); // heat flux
	
	return;
}


void melt::rsdl(int stage) {
	int i,m,n,sind,indx,v0,v1;
	TinyVector<FLT,tri_mesh::ND> norm, rp;
	Array<FLT,1> ubar(x.NV);
	Array<TinyVector<FLT,MXGP>,1> u(x.NV);
	TinyMatrix<FLT,tri_mesh::ND,MXGP> crd, dcrd, mvel;
	TinyMatrix<FLT,8,MXGP> res;
	
	Array<TinyVector<FLT,MXTM>,1> lf(x.NV+tri_mesh::ND);
	
	/**************************************************/
	/* DETERMINE MESH RESIDUALS & SURFACE TENSION      */
	/**************************************************/
	for(n=0;n<tri_mesh::ND;++n)
		gbl->vres(0)(n) = 0.0;
	
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
		for(n=0;n<x.ND;++n) {  // sets n = 1 to zero
			gbl->vres(indx)(n) += lf(x.NV+n)(0);
			gbl->vres(indx+1)(n) = lf(x.NV+n)(1);
			for(m=0;m<basis::tri(x.log2p)->sm();++m)
				gbl->sres(indx,m)(n) = lf(x.NV+n)(m+2);
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
				for(n=0;n<1;++n)
					vdres(x.log2p,i)(n) = gbl->fadd(n)*gbl->vres0(i)(n) -gbl->vres(i)(n);
			
			for(i=0;i<base.nseg;++i) 
				for(m=0;m<basis::tri(x.log2p)->sm();++m) 
					for(n=0;n<1;++n)
						sdres(x.log2p,i,m)(n) = gbl->fadd(n)*gbl->sres0(i,m)(n) -gbl->sres(i,m)(n);
			
		}
		for(i=0;i<base.nseg+1;++i) 
			for(n=0;n<1;++n)
				gbl->vres(i)(n) += vdres(x.log2p,i)(n);
		
		for(i=0;i<base.nseg;++i) 
			for(m=0;m<basis::tri(x.log2p)->sm();++m) 
				for(n=0;n<1;++n)
					gbl->sres(i,m)(n) += sdres(x.log2p,i,m)(n);
	}
	else {
		/* ADD TANGENTIAL MESH MOVEMENT SOURCE */    
		for(i=0;i<base.nseg+1;++i)
			gbl->vres(i)(0) += vdres(x.log2p,i)(0);
		
		for(i=0;i<base.nseg;++i)
			for(m=0;m<basis::tri(x.log2p)->sm();++m) 
				gbl->sres(i,m)(0) += sdres(x.log2p,i,m)(0);
	}
}

void melt::rsdl_after(int stage) {
	int i,sind,v0;
	
	
	//	i = 0;
	//	do {
	//		sind = base.seg(i);
	//		v0 = x.seg(sind).pnt(0);
	//		*x.gbl->log << x.gbl->res.v(v0,2) << std::endl;
	//	} while (++i < base.nseg);
	//	v0 = x.seg(sind).pnt(1);
	//	*x.gbl->log << x.gbl->res.v(v0,2) << std::endl;
	
	
	/* Communicate here???? */
	pmatchsolution_snd(0, x.gbl->res.v.data(), 1);
	base.comm_prepare(boundary::all,0,boundary::symmetric);
	base.comm_exchange(boundary::all,0,boundary::symmetric);
	base.comm_wait(boundary::all,0,boundary::symmetric);
	pmatchsolution_rcv(0,x.gbl->res.v.data(),1);
	
	smatchsolution_snd(x.gbl->res.s.data(),0,basis::tri(x.log2p)->sm()-1, basis::tri(x.log2p)->sm());
	base.comm_prepare(boundary::all,0,boundary::symmetric);
	base.comm_exchange(boundary::all,0,boundary::symmetric);
	base.comm_wait(boundary::all,0,boundary::symmetric);
	smatchsolution_rcv(x.gbl->res.s.data(),0,basis::tri(x.log2p)->sm()-1, basis::tri(x.log2p)->sm());
	
	/* Move Heat equation Residual */
	i = 0;
	do {
		sind = base.seg(i);
		v0 = x.seg(sind).pnt(0);
		gbl->vres(i)(1) = x.gbl->res.v(v0,2);
		for(int m=0;m<basis::tri(x.log2p)->sm();++m) 
			gbl->sres(i,m)(1) = x.gbl->res.s(sind,m,2);
	} while (++i < base.nseg);
	v0 = x.seg(sind).pnt(1);
	gbl->vres(i)(1) = x.gbl->res.v(v0,2);
	
}


void melt::vdirichlet() {
	/* Now apply dirichlet B.C.s */
	symbolic::vdirichlet();
	
#ifdef petsc
	int sind,v0;
	
	/* Store rotated vertex residual in r_mesh residual vector */
	r_tri_mesh::global *r_gbl = dynamic_cast<r_tri_mesh::global *>(x.gbl);
	int i = 0;
	do {
		sind = base.seg(i);
		v0 = x.seg(sind).pnt(0);
		/* Rotate residual for better diagonal dominance */
		r_gbl->res(v0)(0) = gbl->vres(i)(0)*gbl->vdt(i)(0,0) +gbl->vres(i)(1)*gbl->vdt(i)(0,1);
		r_gbl->res(v0)(1) = gbl->vres(i)(0)*gbl->vdt(i)(1,0) +gbl->vres(i)(1)*gbl->vdt(i)(1,1);
	} while (++i < base.nseg);
	v0 = x.seg(sind).pnt(1);
	r_gbl->res(v0)(0) = gbl->vres(i)(0)*gbl->vdt(i)(0,0) +gbl->vres(i)(1)*gbl->vdt(i)(0,1);
	r_gbl->res(v0)(1) = gbl->vres(i)(0)*gbl->vdt(i)(1,0) +gbl->vres(i)(1)*gbl->vdt(i)(1,1);
#endif		
}


#ifdef petsc
int melt::petsc_rsdl(Array<double,1> res) {
	int sm = basis::tri(x.log2p)->sm();
	int ind = 0;
	
	/* Rotated for better diagonal dominance */
	for(int j=0;j<base.nseg;++j) {
		for(int m=0;m<sm;++m) {
			res(ind++) = gbl->sres(j,m)(0)*gbl->sdt(j)(0,0) +gbl->sres(j,m)(1)*gbl->sdt(j)(0,1);
			res(ind++) = gbl->sres(j,m)(0)*gbl->sdt(j)(1,0) +gbl->sres(j,m)(1)*gbl->sdt(j)(1,1);
		}
	}
	return(ind);
}
#endif



void melt::minvrt() {
	int i,m,n,indx;
	int last_phase, mp_phase;
	FLT temp;
	
	/* INVERT MASS MATRIX */
	/* LOOP THROUGH SIDES */
	if (basis::tri(x.log2p)->sm() > 0) {
		for(indx = 0; indx<base.nseg; ++indx) {
			/* SUBTRACT SIDE CONTRIBUTIONS TO VERTICES */            
			for (m=0; m <basis::tri(x.log2p)->sm(); ++m) {
				for(n=0;n<tri_mesh::ND;++n)
					gbl->vres(indx)(n) -= basis::tri(x.log2p)->sfmv1d(0,m)*gbl->sres(indx,m)(n);
				for(n=0;n<tri_mesh::ND;++n)
					gbl->vres(indx+1)(n) -= basis::tri(x.log2p)->sfmv1d(1,m)*gbl->sres(indx,m)(n);
			}
		}
	}
	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		x.vbdry(base.vbdry(0))->vloadbuff(boundary::manifolds,&gbl->vres(0)(0),0,1,0);
		x.vbdry(base.vbdry(1))->vloadbuff(boundary::manifolds,&gbl->vres(base.nseg)(0),0,1,0);
		x.vbdry(base.vbdry(0))->comm_prepare(boundary::manifolds,mp_phase,boundary::symmetric);
		x.vbdry(base.vbdry(1))->comm_prepare(boundary::manifolds,mp_phase,boundary::symmetric);
		
		x.vbdry(base.vbdry(0))->comm_exchange(boundary::manifolds,mp_phase,boundary::symmetric);
		x.vbdry(base.vbdry(1))->comm_exchange(boundary::manifolds,mp_phase,boundary::symmetric);
		
		last_phase = true;
		
		last_phase &= x.vbdry(base.vbdry(0))->comm_wait(boundary::manifolds,mp_phase,boundary::symmetric);
		last_phase &= x.vbdry(base.vbdry(1))->comm_wait(boundary::manifolds,mp_phase,boundary::symmetric);
		
		x.vbdry(base.vbdry(0))->vfinalrcv(boundary::manifolds,mp_phase,boundary::symmetric,boundary::average,&gbl->vres(0)(0),0,1,0);
		x.vbdry(base.vbdry(1))->vfinalrcv(boundary::manifolds,mp_phase,boundary::symmetric,boundary::average,&gbl->vres(base.nseg)(0),0,1,0);
	}
	
	if (gbl->is_loop) {
		gbl->vres(0)(1) = 0.5*(gbl->vres(0)(1) +gbl->vres(base.nseg+1)(1));
		gbl->vres(base.nseg+1)(1) = gbl->vres(0)(1);
		gbl->vres(0)(0) = 0.0;
		gbl->vres(base.nseg+1)(0) = 0.0;
	}
	x.hp_vbdry(base.vbdry(0))->vdirichlet();
	x.hp_vbdry(base.vbdry(1))->vdirichlet();
	
	
	/* SOLVE FOR VERTEX MODES */
	for(i=0;i<base.nseg+1;++i) {
		temp                     = gbl->vres(i)(0)*gbl->vdt(i)(0,0) +gbl->vres(i)(1)*gbl->vdt(i)(0,1);
		gbl->vres(i)(1) = gbl->vres(i)(0)*gbl->vdt(i)(1,0) +gbl->vres(i)(1)*gbl->vdt(i)(1,1);
		gbl->vres(i)(0) = temp;
	}
	
	/* SOLVE FOR SIDE MODES */
	if (basis::tri(x.log2p)->sm() > 0) {
		for(indx = 0; indx<base.nseg; ++indx) {                   
			
#ifdef DETAILED_MINV
			for(m=0;m<basis::tri(x.log2p)->sm();++m) {
				for(n=0;n<tri_mesh::ND;++n) {
					gbl->sres(indx,m)(n) -= gbl->vms(indx,n,0,m,0)*gbl->vres(indx)(0);
					gbl->sres(indx,m)(n) -= gbl->vms(indx,n,0,m,1)*gbl->vres(indx+1)(0);
					gbl->sres(indx,m)(n) -= gbl->vms(indx,n,1,m,0)*gbl->vres(indx)(1);
					gbl->sres(indx,m)(n) -= gbl->vms(indx,n,1,m,1)*gbl->vres(indx+1)(1);                            
				}
			}
			int info;
			char trans[] = "T";
			GETRS(trans,2*basis::tri(x.log2p)->sm(),1,&gbl->ms(indx)(0,0),2*MAXP,&gbl->ipiv(indx)(0),&gbl->sres(indx,0)(0),2*MAXP,info);
			if (info != 0) {
				*x.gbl->log << "DGETRS FAILED FOR SIDE MODE UPDATE" << std::endl;
				sim::abort(__LINE__,__FILE__,x.gbl->log);
			}
#else
			
			/* INVERT SIDE MODES */
			DPBTRSNU2((double *) &basis::tri(x.log2p)->sdiag1d(0,0),basis::tri(x.log2p)->sbwth()+1,basis::tri(x.log2p)->sm(),basis::tri(x.log2p)->sbwth(),&(gbl->sres(indx,0)(0)),tri_mesh::ND);
			for(m=0;m<basis::tri(x.log2p)->sm();++m) {
				temp                             = gbl->sres(indx,m)(0)*gbl->sdt(indx)(0,0) +gbl->sres(indx,m)(1)*gbl->sdt(indx)(0,1);
				gbl->sres(indx,m)(1) = gbl->sres(indx,m)(0)*gbl->sdt(indx)(1,0) +gbl->sres(indx,m)(1)*gbl->sdt(indx)(1,1);         
				gbl->sres(indx,m)(0) = temp;
			}
			
			for(m=0;m<basis::tri(x.log2p)->sm();++m) {
				for(n=0;n<tri_mesh::ND;++n) {
					gbl->sres(indx,m)(n) -= basis::tri(x.log2p)->vfms1d(0,m)*gbl->vres(indx)(n);
					gbl->sres(indx,m)(n) -= basis::tri(x.log2p)->vfms1d(1,m)*gbl->vres(indx+1)(n);
				}
			}
#endif
		}
	}
	
	return;
}

void melt::update(int stage) {
	int i,m,n,count,sind,indx,v0;
	
	if (stage < 0) {
		indx = 0;
		int i = 0;
		do {
			sind = base.seg(i);
			v0 = x.seg(sind).pnt(0);
			gbl->vug0(i) = x.pnts(v0);
		} while (++i < base.nseg);
		v0 = x.seg(sind).pnt(1);
		gbl->vug0(base.nseg) = x.pnts(v0);
		
		if (basis::tri(x.log2p)->sm() > 0) gbl->sug0(Range(0,base.nseg-1),Range(0,basis::tri(x.log2p)->sm()-1)) = crv(Range(0,base.nseg-1),Range(0,basis::tri(x.log2p)->sm()-1));
		return;
	}
	
	minvrt();
	
#ifdef DEBUG
	//   if (x.coarse_flag) {
	for(i=0;i<base.nseg+1;++i)
		*x.gbl->log << "vdt: " << i << ' ' << gbl->vdt(i)(0,0) << ' ' << gbl->vdt(i)(0,1) << ' ' << gbl->vdt(i)(1,0) << ' ' << gbl->vdt(i)(1,1) << '\n';
	
	for(i=0;i<base.nseg;++i)
		*x.gbl->log << "sdt: " << i << ' ' << gbl->sdt(i)(0,0) << ' ' << gbl->sdt(i)(0,1) << ' ' << gbl->sdt(i)(1,0) << ' ' << gbl->sdt(i)(1,1) << '\n';
	
	for(i=0;i<base.nseg+1;++i) {
		*x.gbl->log << "vres: " << i << ' ';
		for(n=0;n<tri_mesh::ND;++n) {
			if (fabs(gbl->vres(i)(n)) > 1.0e-9) *x.gbl->log << gbl->vres(i)(n) << ' ';
			else *x.gbl->log << "0.0 ";
		}
		*x.gbl->log << '\n';
	}
	
	for(i=0;i<base.nseg;++i) {
		for(m=0;m<basis::tri(x.log2p)->sm();++m) {
			*x.gbl->log << "sres: " << i << ' ';
			for(n=0;n<tri_mesh::ND;++n) {
				if (fabs(gbl->sres(i,m)(n)) > 1.0e-9) *x.gbl->log << gbl->sres(i,m)(n) << ' ';
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
		x.pnts(v0) = gbl->vug0(i) -x.gbl->alpha(stage)*gbl->vres(i);
	} while (++i < base.nseg);
	v0 = x.seg(sind).pnt(1);
	x.pnts(v0) = gbl->vug0(base.nseg) -x.gbl->alpha(stage)*gbl->vres(base.nseg);
	
	if (basis::tri(x.log2p)->sm() > 0) crv(Range(0,base.nseg-1),Range(0,basis::tri(x.log2p)->sm()-1)) = gbl->sug0(Range(0,base.nseg-1),Range(0,basis::tri(x.log2p)->sm()-1)) -x.gbl->alpha(stage)*gbl->sres(Range(0,base.nseg-1),Range(0,basis::tri(x.log2p)->sm()-1));
	
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
	
	return;
}

void melt::mg_restrict() {
	int i,bnum,indx,tind,v0,snum,sind;
	
	if(x.p0 > 1) {
		/* TRANSFER IS ON FINEST MESH */
		gbl->vres0(Range(0,base.nseg)) = gbl->vres(Range(0,base.nseg));
		if (basis::tri(x.log2p)->sm() > 0) gbl->sres0(Range(0,base.nseg-1),Range(0,basis::tri(x.log2p)->sm()-1)) = gbl->sres(Range(0,base.nseg-1),Range(0,basis::tri(x.log2p)->sm()-1));
		return;
	}
	else {
		/* TRANSFER IS BETWEEN DIFFERENT MESHES */
		gbl->vres0(Range(0,base.nseg)) = 0.0;
		
		/* CALCULATE COARSE RESIDUALS */
		/* DO ENDPOINTS FIRST */
		gbl->vres0(0) = gbl->vres(0);
		gbl->vres0(base.nseg) = gbl->vres(fine->base.nseg);
		
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
			gbl->vres0(indx) += fmesh->ccnnct(v0).wt((snum+1)%3)*gbl->vres(i);
			gbl->vres0(indx+1) += fmesh->ccnnct(v0).wt((snum+2)%3)*gbl->vres(i);
		}
	}
	
	return;
}


void melt::maxres() {
	int i,n;
	TinyVector<FLT,tri_mesh::ND> mxr;
	
	mxr = 0.0;
	
	for(i=0;i<base.nseg+1;++i)
		for(n=0;n<tri_mesh::ND;++n)
			mxr(n) = MAX(fabs(gbl->vres(i)(n)),mxr(n));
	
	for(n=0;n<tri_mesh::ND;++n)
		*x.gbl->log << ' ' << mxr(n) << ' ';
	
	return;
}


void melt::element_jacobian(int indx, Array<FLT,2>& K) {
	int sm = basis::tri(x.log2p)->sm();	
	Array<TinyVector<FLT,MXTM>,1> Rbar(x.NV+tri_mesh::ND),lf(x.NV+tri_mesh::ND);
#ifdef BZ_DEBUG
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
				for(int n=0;n<x.NV+tri_mesh::ND;++n)
					K(krow++,kcol) = (lf(n)(k)-Rbar(n)(k))/dw(var);
			
			++kcol;
			x.uht(var)(mode) -= dw(var);
		}
		
		for(int var = 0; var < tri_mesh::ND; ++var) {
			x.pnts(x.seg(sind).pnt(mode))(var) += dx;
			
			element_rsdl(indx,lf);
			
			int krow = 0;
			for(int k=0;k<sm+2;++k)
				for(int n=0;n<x.NV+tri_mesh::ND;++n)
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
				for(int n=0;n<x.NV+tri_mesh::ND;++n)
					K(krow++,kcol) = (lf(n)(k)-Rbar(n)(k))/dw(var);
			
			++kcol;
			x.uht(var)(mode) -= dw(var);
		}
		
		for(int var = 0; var < tri_mesh::ND; ++var) {
			crds(indx,mode-2,var) += dx;
			
			element_rsdl(indx,lf);
			
			int krow = 0;
			for(int k=0;k<sm+2;++k)
				for(int n=0;n<x.NV+tri_mesh::ND;++n)
					K(krow++,kcol) = (lf(n)(k)-Rbar(n)(k))/dx;
			
			++kcol;
			
			crds(indx,mode-2,var) -= dx;
		}
	}
	
	return;
}

void melt::smatchsolution_snd(FLT *sdata, int bgn, int end, int stride) {
	int j,k,n,countup,offset;
	
	if (!base.is_comm()) return;
	
#ifdef MPDEBUG
	*x.gbl->log << base.idprefix << " In surface_snd"  << base.idnum << " " << base.is_frst() << std::endl;
#endif
	
	/* This boundary should always be first so send going up */
	countup = 0;
	for(j=0;j<base.nseg;++j) {
		offset = base.seg(j)*stride*x.NV;
		for(k=bgn;k<=end;++k) {
			for(n=x.ND;n<x.ND+1;++n) {
				base.fsndbuf(countup++) = sdata[offset +k*x.NV +n];
#ifdef MPDEBUG
				*x.gbl->log << "\t" << sdata[offset +k*x.NV +n] << std::endl;
#endif
			}
		}
	}
	base.sndsize() = countup;
	base.sndtype() = boundary::flt_msg;
	return;
}

void melt::smatchsolution_rcv(FLT *sdata, int bgn, int end, int stride) {
	
	if (!base.is_comm()) return;
	
	/* CAN'T USE sfinalrcv BECAUSE OF CHANGING SIGNS */
	int j,k,m,n,count,countup,offset,sind;
	FLT mtchinv;
	
	/* OPPOSING BOUNDARY SENDS BACKWARDS SO THIS CAN GO UP */
	
	int matches = 1;
	
	/* RELOAD FROM BUFFER */
	/* ELIMINATES V/S/F COUPLING IN ONE PHASE */
	/* FINALRCV SHOULD BE CALLED F,S,V ORDER (V HAS FINAL AUTHORITY) */    
	for(m=0;m<base.nmatches();++m) {            
		++matches;
		
		int ebp1 = end-bgn+1;
		
		countup = 0;
		for(j=0;j<base.nseg;++j) {
			for(k=0;k<ebp1;++k) {
				for(n=0;n<1;++n) {
					base.fsndbuf(countup) += base.frcvbuf(m,countup);
					++countup;
				}
			}
		}
	}
	
	if (matches > 1) {
		mtchinv = 1./matches;
		
#ifdef MPDEBUG
		*x.gbl->log << base.idprefix << " In surface_rcv"  << base.idnum << " " << base.is_frst() << std::endl;
#endif
		count = 0;
		for(j=0;j<base.nseg;++j) {
			sind = base.seg(j);
			offset = sind*stride*x.NV;
			for (k=bgn;k<=end;++k) {
				for(n=x.ND;n<x.ND+1;++n) {
					sdata[offset +k*x.NV +n] = base.fsndbuf(count++)*mtchinv;
#ifdef MPDEBUG
					*x.gbl->log << "\t" << sdata[offset +k*x.NV +n] << std::endl;
#endif
				}
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
		
		if (!sm) continue;
		
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
		
		ind = 0;
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
		
		gindxND = jacobian_start +j*tri_mesh::ND*sm;
		for(int m=0;m<sm*tri_mesh::ND;++m)
			loc_to_glo_crv(m) = gindxND++;
		
#ifdef MY_SPARSE
		x.J.add_values(tm*x.NV,loc_to_glo_e,sm*tri_mesh::ND,loc_to_glo_crv,Ke);
#else
		MatSetValuesLocal(x.petsc_J,tm*x.NV,loc_to_glo_e.data(),sm*tri_mesh::ND,loc_to_glo_crv.data(),Ke.data(),ADD_VALUES);
#endif		
	}
	
	/* FIXME: NOT SURE ABOUT THIS TEMPORARY */
	//	x.hp_vbdry(base.vbdry(0))->petsc_jacobian();
	//	x.hp_vbdry(base.vbdry(1))->petsc_jacobian();
	//	x.hp_vbdry(base.vbdry(0))->petsc_jacobian_dirichlet();
	//	x.hp_vbdry(base.vbdry(1))->petsc_jacobian_dirichlet();
}

void melt::petsc_matchjacobian_snd() {	
	
	/* Now do stuff for communication boundaries */
	/* Just sending indices of x & y vertices for equality constraint */
	int vdofs;
	if (x.mmovement != x.coupled_deformable)
		vdofs = x.NV;
	else
		vdofs = x.NV+x.ND;
	
	if (!base.is_comm()) return;
	
	int sind=-2;
	
#ifdef MY_SPARSE
	/* Send Jacobian indices for T,x,y posiions for continuity constraint */
	/* I am cheating here and sending floats and int's together */
	base.sndsize() = 0;
	base.sndtype() = boundary::flt_msg;
	base.fsndbuf(base.sndsize()++) = x.jacobian_start +FLT_EPSILON;
	
	for(int i=0;i<base.nseg;++i) {
		sind = base.seg(i);
		
		int rowbase = x.seg(sind).pnt(0)*vdofs;  // Base index		
		/* attach diagonal column # to allow vertex continuity enforcement */
		base.fsndbuf(base.sndsize()++) = rowbase +FLT_EPSILON;;
#ifdef MPDEBUG
		*x.gbl->log << "sending vertex rowbase " << rowbase << std::endl;
#endif
		
		/* Send Side Information */
		rowbase = x.npnt*vdofs +sind*x.NV*x.sm0;  // Base index
		/* attach diagonal column # to allow side continuity enforcement */
		base.fsndbuf(base.sndsize()++) = rowbase +FLT_EPSILON;
#ifdef MPDEBUG
		*x.gbl->log << "sending side rowbase " << rowbase << std::endl;
#endif
	}
	
	/* LAST POINT */
	int rowbase = x.seg(sind).pnt(1)*vdofs; 
	/* attach diagonal # to allow continuity enforcement */
	base.fsndbuf(base.sndsize()++) = rowbase +FLT_EPSILON;
#ifdef MPDEBUG
	*x.gbl->log << "sending vertex rowbase " << rowbase << std::endl;
#endif	
	
	/* Send index of start of curved modes */
	base.fsndbuf(base.sndsize()++) = jacobian_start;
#ifdef MPDEBUG
	*x.gbl->log << "sending jacobian_start " << jacobian_start << std::endl;
#endif	
	
}

void melt::petsc_matchjacobian_rcv(int phase) {
	const int ND = x.ND;
	
	if (!base.is_comm() || base.matchphase(boundary::all_phased,0) != phase) return;
	
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
	int row;
	
	std::vector<int> c0vars;
	c0vars.push_back(ND);
	
	/* Now Receive Information */		
	for (int i=base.nseg-1;i>=0;--i) {
		int sind = base.seg(i);
		int rowbase = x.seg(sind).pnt(1)*vdofs; 
		int row_mpi = static_cast<int>(base.frcvbuf(0,count++)) +Jstart_mpi;
		
		for(std::vector<int>::iterator it=c0vars.begin();it!=c0vars.end();++it) {
			row = rowbase + *it;
			int ncol = static_cast<int>(base.frcvbuf(0,count++));
#ifdef MPDEBUG
			*x.gbl->log << "receiving " << ncol << " jacobian entries for vertex " << row/vdofs << " and variable " << *it << std::endl;
#endif
			for (int k = 0;k<ncol;++k) {
				int col = static_cast<int>(base.frcvbuf(0,count++)) +Jstart_mpi;
#ifdef MPDEBUG
				*x.gbl->log  << col << ' ';
#endif
				FLT val = base.frcvbuf(0,count++);
				(*pJ_mpi).add_values(row,col,val);
			}
#ifdef MPDEBUG
			*x.gbl->log << std::endl;
#endif
			/* Shift all entries for this vertex.  Remote variables on solid are T */
			for(int n=0;n<1;++n) {
				FLT dval = x.J_mpi(row,row_mpi +n);
				(*pJ_mpi)(row,row_mpi +n) = 0.0;				
				x.J(row,rowbase +c0vars[n]) += dval;
			}
			x.J.multiply_row(row,0.5);
			x.J_mpi.multiply_row(row,0.5);
		}
		
		/* Now receive side Jacobian information */
		row = x.npnt*vdofs +sind*x.NV*x.sm0;
		row_mpi = static_cast<int>(base.frcvbuf(0,count++)) +Jstart_mpi;
		
		int mcnt = x.ND;
		int sgn = 1;
		for(int mode=0;mode<x.sm0;++mode) {
			int ncol = static_cast<int>(base.frcvbuf(0,count++));
#ifdef MPDEBUG
			*x.gbl->log << "receiving " << ncol << " jacobian entries for side " << sind << " and variable " << 2 << std::endl;
#endif
			for (int k = 0;k<ncol;++k) {
				int col = static_cast<int>(base.frcvbuf(0,count++)) +Jstart_mpi;
#ifdef MPDEBUG
				*x.gbl->log  << col << ' ';
#endif
				FLT val = sgn*base.frcvbuf(0,count++);
				(*pJ_mpi).add_values(row+mcnt,col,val);
			}
#ifdef MPDEBUG
			*x.gbl->log << std::endl;
#endif
			
			/* Shift all modes in equation */
			int mcnt_mpi = 0;
			int sgn_mpi = 1;
			for(int mode_mpi=0;mode_mpi<x.sm0;++mode_mpi) {
				FLT dval = x.J_mpi(row+mcnt,row_mpi+mcnt_mpi);
				(*pJ_mpi)(row+mcnt,row_mpi+mcnt_mpi) = 0.0;				
				x.J(row+mcnt,row+x.ND+mode_mpi*x.NV) += sgn_mpi*dval;
				sgn_mpi *= -1;
				mcnt_mpi += 1; // Only temperature on solid block
			}
			x.J.multiply_row(row+mcnt,0.5);
			x.J_mpi.multiply_row(row+mcnt,0.5);
			mcnt += x.NV;
			sgn *= -1;
		}
	}
	
	int sind = base.seg(0);
	int rowbase = x.seg(sind).pnt(0)*vdofs; 
	int row_mpi = static_cast<int>(base.frcvbuf(0,count++)) +Jstart_mpi;
	
	for(std::vector<int>::iterator it=c0vars.begin();it!=c0vars.end();++it) {
		row = rowbase + *it;
		int ncol = static_cast<int>(base.frcvbuf(0,count++));
#ifdef MPDEBUG
		*x.gbl->log << "receiving " << ncol << " jacobian entries for vertex " << row/vdofs << " and variable " << *it << std::endl;
#endif
		for (int k = 0;k<ncol;++k) {
			int col = static_cast<int>(base.frcvbuf(0,count++)) +Jstart_mpi;
#ifdef MPDEBUG
			*x.gbl->log << col << ' ';
#endif
			FLT val = base.frcvbuf(0,count++);
			(*pJ_mpi).add_values(row,col,val);
		}
#ifdef MPDEBUG
		*x.gbl->log << std::endl;
#endif
		/* Shift all entries for this vertex.  Remote variables on solid are T, x, y */
		for(int n=0;n<1;++n) {
			FLT dval = x.J_mpi(row,row_mpi +n);
			(*pJ_mpi)(row,row_mpi +n) = 0.0;				
			x.J(row,rowbase +c0vars[n]) += dval;
		}
		x.J.multiply_row(row,0.5);
		x.J_mpi.multiply_row(row,0.5);
	} 
#endif
	
	if (sm && !base.is_frst()) {
		/* Equality of curved mode constraint */
		int ind = jacobian_start;
		int ind_mpi = static_cast<int>(base.frcvbuf(0,count++)) +(base.nseg-1)*sm*x.ND +Jstart_mpi;
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


void melt::petsc_jacobian_dirichlet() {
	int sind, v0, row;
	int sm=basis::tri(x.log2p)->sm();
	Array<int,1> indices((base.nseg+1)*(x.NV-1) +base.nseg*sm*(x.NV-1));
	
	int vdofs;
	if (x.mmovement == x.coupled_deformable)
		vdofs = x.NV +tri_mesh::ND;
	else
		vdofs = x.NV;
	
	int begin_seg = x.npnt*vdofs;
	
	/****************************************/
	/* MOVE HEAT EQUATION IN J TO LAST ROW  */
	/****************************************/
	int j = 0;
	do {
		sind = base.seg(j);
		v0 = x.seg(sind).pnt(0);
		row = v0*vdofs +x.NV-2;
		int nnz1 = x.J._cpt(row+1) -x.J._cpt(row);
		int nnz2 = x.J._cpt(row+4) -x.J._cpt(row+3);
		/* SOME ERROR CHECKING TO MAKE SURE ROW SPARSENESS PATTERN IS THE SAME */
		if (nnz1 != nnz2) {
			*x.gbl->log << "zeros problem in moving heat equation to mesh movement location\n";
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}
		int cpt1 = x.J._cpt(row);
		int cpt2 = x.J._cpt(row+3);
		for(int col=0;col<nnz1;++col) {
			/* Overwrite row */
			x.J._col(cpt2) = x.J._col(cpt1);
			x.J._val(cpt2++) = x.J._val(cpt1++);
		}
		
		for (int m=0;m<sm;++m) {
			int row1 = begin_seg +sind*sm*x.NV +m*x.NV +x.NV-2;
			int row2 = jacobian_start +j*sm*x.ND +m*x.ND +1;
			int nnz1 = x.J._cpt(row1+1) -x.J._cpt(row1);
			int nnz2 = x.J._cpt(row2+1) -x.J._cpt(row2);
			/* SOME ERROR CHECKING TO MAKE SURE ROW SPARSENESS PATTERN IS THE SAME */
			if (nnz1 != nnz2) {
				*x.gbl->log << "zeros problem in moving heat equation to mesh movement location " << nnz1 << ' ' << nnz2 << std::endl;
				sim::abort(__LINE__,__FILE__,x.gbl->log);
			}
			cpt1 = x.J._cpt(row1);
			cpt2 = x.J._cpt(row2);
			for(int col=0;col<nnz1;++col) {
				/* Overwrite row */
				x.J._col(cpt2) = x.J._col(cpt1);
				x.J._val(cpt2++) = x.J._val(cpt1++);
			}
		}
	} while (++j < base.nseg);
	v0 = x.seg(sind).pnt(1);
	row = v0*vdofs +x.NV-2;
	int nnz1 = x.J._cpt(row+1) -x.J._cpt(row);
	int nnz2 = x.J._cpt(row+4) -x.J._cpt(row+3);
	/* SOME ERROR CHECKING TO MAKE SURE ROW SPARSENESS PATTERN IS THE SAME */
	if (nnz1 != nnz2) {
		*x.gbl->log << "zeros problem in moving heat equation to mesh movement location\n";
		sim::abort(__LINE__,__FILE__,x.gbl->log);
	}
	int cpt1 = x.J._cpt(row);
	int cpt2 = x.J._cpt(row+3);
	for(int col=0;col<nnz1;++col) {
		/* Overwrite row */
		x.J._col(cpt2) = x.J._col(cpt1);
		x.J._val(cpt2++) = x.J._val(cpt1++);
	}
	
	/********************************************/
	/* MOVE HEAT EQUATION IN J_mpi TO LAST ROW  */
	/********************************************/
	j = 0;
	do {
		sind = base.seg(j);
		v0 = x.seg(sind).pnt(0);
		row = v0*vdofs +x.NV-2;
		int nnz1 = x.J_mpi._cpt(row+1) -x.J_mpi._cpt(row);
		int nnz2 = x.J_mpi._cpt(row+4) -x.J_mpi._cpt(row+3);
		/* SOME ERROR CHECKING TO MAKE SURE ROW SPARSENESS PATTERN IS THE SAME */
		if (nnz1 != nnz2) {
			*x.gbl->log << "zeros problem in moving heat equation to mesh movement location " << v0 << ' ' << nnz1 << ' ' << nnz2 << std::endl;;
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}
		int cpt1 = x.J_mpi._cpt(row);
		int cpt2 = x.J_mpi._cpt(row+3);
		
		for(int col=0;col<nnz1;++col) {
			/* Overwrite row */
			x.J_mpi._col(cpt2) = x.J_mpi._col(cpt1);
			x.J_mpi._val(cpt2++) = x.J_mpi._val(cpt1++);
		}
		
		for (int m=0;m<sm;++m) {
			int row1 = begin_seg +sind*sm*x.NV +m*x.NV +x.NV-2;
			int row2 = jacobian_start +j*sm*x.ND +m*x.ND +1;
			int nnz1 = x.J_mpi._cpt(row1+1) -x.J_mpi._cpt(row1);
			int nnz2 = x.J_mpi._cpt(row2+1) -x.J_mpi._cpt(row2);
			/* SOME ERROR CHECKING TO MAKE SURE ROW SPARSENESS PATTERN IS THE SAME */
			if (nnz1 != nnz2) {
				*x.gbl->log << "zeros problem in moving heat equation to mesh movement location " << nnz1 << ' ' << nnz2 << std::endl;
				sim::abort(__LINE__,__FILE__,x.gbl->log);
			}
			int cpt1 = x.J_mpi._cpt(row1);
			int cpt2 = x.J_mpi._cpt(row2);
			for(int col=0;col<nnz1;++col) {
				/* Overwrite row */
				x.J_mpi._col(cpt2) = x.J_mpi._col(cpt1);
				x.J_mpi._val(cpt2++) = x.J_mpi._val(cpt1++);
			}
		}
	} while (++j < base.nseg);
	v0 = x.seg(sind).pnt(1);
	row = v0*vdofs +x.NV-2;
	nnz1 = x.J_mpi._cpt(row+1) -x.J_mpi._cpt(row);
	nnz2 = x.J_mpi._cpt(row+4) -x.J_mpi._cpt(row+3);
	/* SOME ERROR CHECKING TO MAKE SURE ROW SPARSENESS PATTERN IS THE SAME */
	if (nnz1 != nnz2) {
		*x.gbl->log << "zeros problem in moving heat equation to mesh movement location\n";
		sim::abort(__LINE__,__FILE__,x.gbl->log);
	}
	cpt1 = x.J_mpi._cpt(row);
	cpt2 = x.J_mpi._cpt(row+3);
	for(int col=0;col<nnz1;++col) {
		/* Overwrite row */
		x.J_mpi._col(cpt2) = x.J_mpi._col(cpt1);
		x.J_mpi._val(cpt2++) = x.J_mpi._val(cpt1++);
	}
	
	
	/***********************************/
	/* Rotate J for diagonal dominance */
	/***********************************/
	j = 0;
	do {
		sind = base.seg(j);
		v0 = x.seg(sind).pnt(0);
		row = v0*vdofs +x.NV;
		int nnz1 = x.J._cpt(row+1) -x.J._cpt(row);
		int nnz2 = x.J._cpt(row+2) -x.J._cpt(row+1);
		/* SOME ERROR CHECKING TO MAKE SURE ROW SPARSENESS PATTERN IS THE SAME */
		if (nnz1 != nnz2) {
			*x.gbl->log << "zeros problem in moving heat equation to mesh movement location\n";
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}
		int cpt1 = x.J._cpt(row);
		int cpt2 = x.J._cpt(row+1);
		for(int col=0;col<nnz2;++col) {
			if (x.J._col(cpt1) != x.J._col(cpt2)) {						
				*x.gbl->log << "zeros indexing problem in moving heat equation to mesh movement location " << col << ' ' << x.J._col(cpt1) << ' ' << x.J._col(cpt2) << std::endl;
				sim::abort(__LINE__,__FILE__,x.gbl->log);
			}	
			double row_store = x.J._val(cpt1)*gbl->vdt(j)(0,0) +x.J._val(cpt2)*gbl->vdt(j)(0,1);
			x.J._val(cpt2) = x.J._val(cpt1)*gbl->vdt(j)(1,0) +x.J._val(cpt2)*gbl->vdt(j)(1,1);
			x.J._val(cpt1) = row_store;	
			++cpt1;
			++cpt2;
		}
		
		for (int m=0;m<sm;++m) {
			int row1 = jacobian_start +j*sm*x.ND +m*x.ND;
			int row2 = jacobian_start +j*sm*x.ND +m*x.ND +1;
			int nnz1 = x.J._cpt(row1+1) -x.J._cpt(row1);
			int nnz2 = x.J._cpt(row2+1) -x.J._cpt(row2);
			/* SOME ERROR CHECKING TO MAKE SURE ROW SPARSENESS PATTERN IS THE SAME */
			if (nnz1 != nnz2) {
				*x.gbl->log << "zeros problem in moving heat equation to mesh movement location\n";
				sim::abort(__LINE__,__FILE__,x.gbl->log);
			}
			
			/* Fill in zeros */
			int cpt1 = x.J._cpt(row1);
			int cpt2 = x.J._cpt(row2);
			for(int col=0;col<nnz2;++col) {
				if (x.J._col(cpt2) < x.J._col(cpt1)) {
					x.J.set_values(row1, x.J._col(cpt2), 0.0); 
				}
				++cpt2;
				++cpt1;
			}
			
			cpt1 = x.J._cpt(row1);
			cpt2 = x.J._cpt(row2);
			for(int col=0;col<nnz1;++col) {
				if (x.J._col(cpt1) != x.J._col(cpt2)) {
					*x.gbl->log << "zeros indexing problem in moving heat equation to mesh movement location\n";
					sim::abort(__LINE__,__FILE__,x.gbl->log);
				}	
				double row_store = x.J._val(cpt1)*gbl->sdt(j)(0,0) +x.J._val(cpt2)*gbl->sdt(j)(0,1);
				x.J._val(cpt2) = x.J._val(cpt1)*gbl->sdt(j)(1,0) +x.J._val(cpt2)*gbl->sdt(j)(1,1);
				x.J._val(cpt1) = row_store;	
				++cpt1;
				++cpt2;
			}
		}
	} while (++j < base.nseg);
	v0 = x.seg(sind).pnt(1);
	row = v0*vdofs +x.NV;
	nnz1 = x.J._cpt(row+1) -x.J._cpt(row);
	nnz2 = x.J._cpt(row+2) -x.J._cpt(row+1);
	/* SOME ERROR CHECKING TO MAKE SURE ROW SPARSENESS PATTERN IS THE SAME */
	if (nnz1 != nnz2) {
		*x.gbl->log << "zeros problem in moving heat equation to mesh movement location\n";
		sim::abort(__LINE__,__FILE__,x.gbl->log);
	}
	cpt1 = x.J._cpt(row);
	cpt2 = x.J._cpt(row+1);
	for(int col=0;col<nnz1;++col) {
		if (x.J._col(cpt1) != x.J._col(cpt2)) {
			*x.gbl->log << "zeros indexing problem in deforming mesh on angled boundary\n";
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}	
		double row_store = x.J._val(cpt1)*gbl->vdt(j)(0,0) +x.J._val(cpt2)*gbl->vdt(j)(0,1);
		x.J._val(cpt2) = x.J._val(cpt1)*gbl->vdt(j)(1,0) +x.J._val(cpt2)*gbl->vdt(j)(1,1);
		x.J._val(cpt1) = row_store;	
		++cpt1;
		++cpt2;
	}
	
	
	/***************************************/
	/* Rotate J_mpi for diagonal dominance */
	/***************************************/
	/* The x-equation has no entries it is only the tangential equation */
	j = 0;
	do {
		sind = base.seg(j);
		v0 = x.seg(sind).pnt(0);
		row = v0*vdofs +x.NV;
		int nnz1 = x.J_mpi._cpt(row+1) -x.J_mpi._cpt(row);
		int nnz2 = x.J_mpi._cpt(row+2) -x.J_mpi._cpt(row+1);
		/* SOME ERROR CHECKING TO MAKE SURE ROW SPARSENESS PATTERN IS THE SAME */
		if (nnz1 != nnz2) {
			*x.gbl->log << "zeros problem in moving heat equation to mesh movement location\n";
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}
		int cpt1 = x.J_mpi._cpt(row);
		int cpt2 = x.J_mpi._cpt(row+1);
		for(int col=0;col<nnz2;++col) {
			x.J_mpi._col(cpt1) = x.J_mpi._col(cpt2);
			double row_store = x.J_mpi._val(cpt2)*gbl->vdt(j)(0,1);
			x.J_mpi._val(cpt2) = x.J_mpi._val(cpt2)*gbl->vdt(j)(1,1);
			x.J_mpi._val(cpt1) = row_store;	
			++cpt1;
			++cpt2;
		}
		
		for (int m=0;m<sm;++m) {
			int row1 = jacobian_start +j*sm*x.ND +m*x.ND;
			int row2 = jacobian_start +j*sm*x.ND +m*x.ND +1;
			int nnz1 = x.J_mpi._cpt(row1+1) -x.J_mpi._cpt(row1);
			int nnz2 = x.J_mpi._cpt(row2+1) -x.J_mpi._cpt(row2);
			/* SOME ERROR CHECKING TO MAKE SURE ROW SPARSENESS PATTERN IS THE SAME */
			if (nnz1 != nnz2) {
				*x.gbl->log << "zeros problem in moving heat equation to mesh movement location\n";
				sim::abort(__LINE__,__FILE__,x.gbl->log);
			}
			
			cpt1 = x.J_mpi._cpt(row1);
			cpt2 = x.J_mpi._cpt(row2);
			/* Fill in zeros */
			for(int col=0;col<nnz2;++col) {
				if (x.J_mpi._col(cpt2) < x.J_mpi._col(cpt1)) {
					x.J_mpi.set_values(row1, x.J_mpi._col(cpt2), 0.0); 
				}
				++cpt2;
				++cpt1;
			}
			
			cpt1 = x.J_mpi._cpt(row1);
			cpt2 = x.J_mpi._cpt(row2);
			for(int col=0;col<nnz2;++col) {
				if (x.J_mpi._col(cpt1) != x.J_mpi._col(cpt2)) {
					*x.gbl->log << "zeros indexing problem in moving heat equation to mesh movement location "  << 
					row1 << ' ' <<  x.J_mpi._col(cpt1) << ' ' << row2 << ' ' << x.J_mpi._col(cpt2) << std::endl;
					sim::abort(__LINE__,__FILE__,x.gbl->log);
				}	
				double row_store = x.J_mpi._val(cpt1)*gbl->sdt(j)(0,0) +x.J_mpi._val(cpt2)*gbl->sdt(j)(0,1);
				x.J_mpi._val(cpt2) = x.J_mpi._val(cpt1)*gbl->sdt(j)(1,0) +x.J_mpi._val(cpt2)*gbl->sdt(j)(1,1);
				x.J_mpi._val(cpt1) = row_store;	
				++cpt1;
				++cpt2;
			}
		}
	} while (++j < base.nseg);
	v0 = x.seg(sind).pnt(1);
	row = v0*vdofs +x.NV;
	nnz1 = x.J_mpi._cpt(row+1) -x.J_mpi._cpt(row);
	nnz2 = x.J_mpi._cpt(row+2) -x.J_mpi._cpt(row+1);
	/* SOME ERROR CHECKING TO MAKE SURE ROW SPARSENESS PATTERN IS THE SAME */
	if (nnz1 != nnz2) {
		*x.gbl->log << "zeros problem in moving heat equation to mesh movement location\n";
		sim::abort(__LINE__,__FILE__,x.gbl->log);
	}
	cpt1 = x.J_mpi._cpt(row);
	cpt2 = x.J_mpi._cpt(row+1);
	for(int col=0;col<nnz1;++col) {
		x.J_mpi._col(cpt1) = x.J_mpi._col(cpt2);
		double row_store = x.J_mpi._val(cpt2)*gbl->vdt(j)(0,1);
		x.J_mpi._val(cpt2) = x.J_mpi._val(cpt2)*gbl->vdt(j)(1,1);
		x.J_mpi._val(cpt1) = row_store;	
		++cpt1;
		++cpt2;
	}
	
	
	symbolic::petsc_jacobian_dirichlet();  // Sets heat equation row to be just 1 on diagonal
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
			if (base.is_frst()) {
				nnzero(Range(pind*vdofs,(pind+1)*vdofs-1)) += ND*sm;
			}
			else {
				nnzero(Range(pind*vdofs,pind*vdofs +x.NV-1)) += ND*sm;
			}
			
			pind = x.seg(sind).pnt(1);
			if (base.is_frst()) {
				nnzero(Range(pind*vdofs,(pind+1)*vdofs-1)) += ND*sm;
			}
			else {
				nnzero(Range(pind*vdofs,pind*vdofs +x.NV-1)) += ND*sm;
			}
		}
		
		/* Side mode heat equation will be copied into this location */
		nnzero(Range(jacobian_start,jacobian_start+base.nseg*sm*tri_mesh::ND-1)) = 3*vdofs +3*x.NV*sm +x.NV*im +ND*sm;
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
	
	
	/* INTERFACE TEMPERATURE EQUATION JACOBIAN WILL BE EXCHANGED */
	/* IT WILL THEN BE REPLACED BY A DIRCHLET CONDITION AND MOVED TO X/Y EQUATION */
	std::vector<int> c0vars;
	c0vars.push_back(ND);
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
	    
	int vdofs;
	if (x.mmovement != tri_hp::coupled_deformable) 
		vdofs = NV;
	else
		vdofs = ND+NV;
	
	int begin_seg = x.npnt*vdofs;
	
	std::vector<int> c0vars;
	c0vars.push_back(ND);
	for(int n=x.NV;n<vdofs;++n) {
		c0vars.push_back(n);
	}		
	
	int count = 0;
	for (int i=base.nseg-1;i>=0;--i) {
		int sind = base.seg(i);
		int pind = x.seg(sind).pnt(1)*vdofs;
		int nentry = base.ircvbuf(0,count);
		for(std::vector<int>::iterator it=c0vars.begin();it!=c0vars.end();++it) {
			// T received from conv-diffusive side;
			nnzero_mpi(pind+*it) += nentry;
			++count; // Skip sizes sent for x & y rows (must match T size)
		}
		
	}
	int sind = base.seg(0);
	int pind = x.seg(sind).pnt(0)*vdofs;
	int nentry = base.ircvbuf(0,count);
	for(std::vector<int>::iterator it=c0vars.begin();it!=c0vars.end();++it) {
		// T received from conv-diffusive side;
		nnzero_mpi(pind+*it) += nentry;
		++count; // Skip sizes sent for x & y rows (must match T size)
	}
	
	
	/* Now add to side degrees of freedom */
	if (sm) {
		int toadd = base.ircvbuf(0,count++); 
		for (int i=0;i<base.nseg;++i) {
			int sind = base.seg(i);
			for (int mode=0;mode<sm;++mode) {
				for(int n=x.ND;n<x.ND+1;++n) {
					nnzero_mpi(begin_seg+sind*NV*sm +mode*NV +n) += toadd;
				}
			}
		}
		/* Make these big enough to hold heat equation */
		nnzero_mpi(Range(jacobian_start,jacobian_start+base.nseg*sm*tri_mesh::ND-1)) += toadd;
	}
}

#endif

void melt::setup_preconditioner() {

	int indx,sind,v0,v1;
	TinyVector<FLT,tri_mesh::ND> nrm;
	FLT h, hsm;
	FLT dttang, dtnorm;
	FLT vslp;
	FLT qmax;
	TinyMatrix<FLT,tri_mesh::ND,MXGP> crd, dcrd;
	TinyMatrix<FLT,4,MXGP> res;
	TinyMatrix<FLT,4,MXGP> lf;
	TinyVector<FLT,2> mvel;
	Array<TinyVector<FLT,MXGP>,1> u(x.NV);
	int last_phase, mp_phase;

	/**************************************************/
	/* DETERMINE SURFACE MOVEMENT TIME STEP              */
	/**************************************************/
	gbl->vdt(0) = 0.0;

	for(indx=0; indx < base.nseg; ++indx) {
		sind = base.seg(indx);
		v0 = x.seg(sind).pnt(0);
		v1 = x.seg(sind).pnt(1);


#ifdef DETAILED_DT
		x.crdtocht1d(sind);
		for(n=0;n<tri_mesh::ND;++n)
			basis::tri(x.log2p)->proj1d(&x.cht(n,0),&crd(n,0),&dcrd(n,0));

		x.ugtouht1d(sind);
		for(n=0;n<tri_mesh::ND;++n)
			basis::tri(x.log2p)->proj1d(&x.uht(n)(0),&u(n)(0));    

		dtnorm = 1.0e99;
		dttang = 1.0e99;
		gbl->meshc(indx) = 1.0e99;
		for(i=0;i<basis::tri(x.log2p)->gpx();++i) {
			nrm(0) =  dcrd(1,i)*2;
			nrm(1) = -dcrd(0,i)*2;
			h = sqrt(nrm(0)*nrm(0) +nrm(1)*nrm(1));

			/* RELATIVE VELOCITY STORED IN MVEL(N)*/
			for(n=0;n<tri_mesh::ND;++n) {
				mvel(n) = u(n)(i) -(x.gbl->bd(0)*(crd(n,i) -dxdt(x.log2p,indx)(n,i))); 
#ifdef MESH_REF_VEL
				mvel(n) -= x.gbl->mesh_ref_vel(n);
#endif
			}

			vslp = fabs(-u(0)(i)*nrm(1)/h +u(1)(i)*nrm(0)/h);
			qmax = mvel(0)*mvel(0)+mvel(1)*mvel(1);
			hsm = h/(.25*(basis::tri(x.log2p)->p()+1)*(basis::tri(x.log2p)->p()+1));

			dttang = MIN(dttang,2.*ksprg(indx)*(.25*(basis::tri(x.log2p)->p()+1)*(basis::tri(x.log2p)->p()+1))/hsm);
			dtnorm = MIN(dtnorm,2.*vslp/hsm +x.gbl->bd(0))*gbl->Lf;                

			/* SET UP DISSIPATIVE COEFFICIENT */
			/* FOR UPWINDING LINEAR CONVECTIVE CASE SHOULD BE 1/|a| */
			/* RESIDUAL HAS DX/2 WEIGHTING */
			/* |a| dx/2 dv/dx  dx/2 dpsi */
			/* |a| dx/2 2/dx dv/dpsi  dpsi */
			/* |a| dv/dpsi  dpsi */
			// gbl->meshc(indx) = gbl->adis/(h*dtnorm*0.5);/* FAILED IN NATES UPSTREAM SURFACE WAVE CASE */
			//gbl->meshc(indx) = MIN(gbl->meshc(indx),gbl->adis/(h*(vslp/hsm +x.gbl->bd(0)))); /* FAILED IN MOVING UP TESTS */
			gbl->meshc(indx) = MIN(gbl->meshc(indx),gbl->adis/(h*(sqrt(qmax)/hsm +x.gbl->bd(0)))); /* SEEMS THE BEST I'VE GOT */
		}
		nrm(0) =  (x.pnts(v1)(1) -x.pnts(v0)(1));
		nrm(1) = -(x.pnts(v1)(0) -x.pnts(v0)(0));
#else
		nrm(0) =  (x.pnts(v1)(1) -x.pnts(v0)(1));
		nrm(1) = -(x.pnts(v1)(0) -x.pnts(v0)(0));
		h = sqrt(nrm(0)*nrm(0) +nrm(1)*nrm(1));

		mvel(0) = x.ug.v(v0,0)-(x.gbl->bd(0)*(x.pnts(v0)(0) -x.vrtxbd(1)(v0)(0)));
		mvel(1) = x.ug.v(v0,1)-(x.gbl->bd(0)*(x.pnts(v0)(1) -x.vrtxbd(1)(v0)(1)));
		
#ifdef MESH_REF_VEL
		mvel -= x.gbl->mesh_ref_vel;
#endif


		qmax = mvel(0)*mvel(0)+mvel(1)*mvel(1);
		vslp = fabs(-mvel(0)*nrm(1)/h +mvel(1)*nrm(0)/h);

		mvel(0) = x.ug.v(v1,0)-(x.gbl->bd(0)*(x.pnts(v1)(0) -x.vrtxbd(1)(v1)(0)));
		mvel(1) = x.ug.v(v1,1)-(x.gbl->bd(0)*(x.pnts(v1)(1) -x.vrtxbd(1)(v1)(1)));
		
#ifdef MESH_REF_VEL
		mvel -= x.gbl->mesh_ref_vel;
#endif

		qmax = MAX(qmax,mvel(0)*mvel(0)+mvel(1)*mvel(1));
		vslp = MAX(vslp,fabs(-mvel(0)*nrm(1)/h +mvel(1)*nrm(0)/h));

		hsm = h/(.25*(basis::tri(x.log2p)->p()+1)*(basis::tri(x.log2p)->p()+1));

		dttang = 2.*ksprg(indx)*(.25*(basis::tri(x.log2p)->p()+1)*(basis::tri(x.log2p)->p()+1))/hsm;
		dtnorm = (2.*vslp/hsm +x.gbl->bd(0))*gbl->Lf;  

		/* SET UP DISSIPATIVE COEFFICIENT */
		/* FOR UPWINDING LINEAR CONVECTIVE CASE SHOULD BE 1/|a| */
		/* RESIDUAL HAS DX/2 WEIGHTING */
		/* |a| dx/2 dv/dx  dx/2 dpsi */
		/* |a| dx/2 2/dx dv/dpsi  dpsi */
		/* |a| dv/dpsi  dpsi */
		// gbl->meshc(indx) = gbl->adis/(h*dtnorm*0.5); /* FAILED IN NATES UPSTREAM SURFACE WAVE CASE */
		// gbl->meshc(indx) = gbl->adis/(h*(vslp/hsm +x.gbl->bd(0))); /* FAILED IN MOVING UP TESTS */
		gbl->meshc(indx) = gbl->adis/(h*(sqrt(qmax)/hsm +x.gbl->bd(0))); /* SEEMS THE BEST I'VE GOT */
#endif
		dtnorm *= RAD(0.5*(x.pnts(v0)(0) +x.pnts(v1)(0)));
		nrm *= 0.5;

		gbl->vdt(indx)(0,0) += -dttang*nrm(1)*basis::tri(x.log2p)->vdiag1d();
		gbl->vdt(indx)(0,1) +=  dttang*nrm(0)*basis::tri(x.log2p)->vdiag1d();
		gbl->vdt(indx)(1,0) +=  dtnorm*nrm(0)*basis::tri(x.log2p)->vdiag1d();
		gbl->vdt(indx)(1,1) +=  dtnorm*nrm(1)*basis::tri(x.log2p)->vdiag1d();
		gbl->vdt(indx+1)(0,0) = -dttang*nrm(1)*basis::tri(x.log2p)->vdiag1d();
		gbl->vdt(indx+1)(0,1) =  dttang*nrm(0)*basis::tri(x.log2p)->vdiag1d();
		gbl->vdt(indx+1)(1,0) =  dtnorm*nrm(0)*basis::tri(x.log2p)->vdiag1d();
		gbl->vdt(indx+1)(1,1) =  dtnorm*nrm(1)*basis::tri(x.log2p)->vdiag1d();

		if (basis::tri(x.log2p)->sm()) {
			gbl->sdt(indx)(0,0) = -dttang*nrm(1);
			gbl->sdt(indx)(0,1) =  dttang*nrm(0);
			gbl->sdt(indx)(1,0) =  dtnorm*nrm(0);
			gbl->sdt(indx)(1,1) =  dtnorm*nrm(1);  

#ifdef DETAILED_MINV
			int lsm = basis::tri(x.log2p)->sm();
			x.crdtocht1d(sind);
			for(n=0;n<tri_mesh::ND;++n)
				basis::tri(x.log2p)->proj1d(&x.cht(n,0),&crd(n,0),&dcrd(n,0));

			for(int m = 0; m<lsm; ++m) {
				for(i=0;i<basis::tri(x.log2p)->gpx();++i) {
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
				lf(0) /= gbl->cfl(0,x.log2p);
				lf(1) /= gbl->cfl(0,x.log2p);
				lf(2) /= gbl->cfl(1,x.log2p);
				lf(3) /= gbl->cfl(1,x.log2p);                        

				for (n=0;n<lsm;++n) {
					gbl->ms(indx)(2*m,2*n) = lf(0,n+2);
					gbl->ms(indx)(2*m,2*n+1) = lf(1,n+2);
					gbl->ms(indx)(2*m+1,2*n) = lf(2,n+2);
					gbl->ms(indx)(2*m+1,2*n+1) = lf(3,n+2);
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
			GETRF(2*lsm,2*lsm,&gbl->ms(indx)(0,0),2*MAXP,&gbl->ipiv(indx)(0),info);
			if (info != 0) {
				*x.gbl->log << "DGETRF FAILED IN SIDE MODE PRECONDITIONER\n");
				sim::abort(__LINE__,__FILE__,x.gbl->log);
			}
			/*
			\phi_n dx,dy*t = \phi_n Vt
			\phi_t dx,dy*n = \phi_t Vn
			*/
#endif
		}
	}

	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		x.vbdry(base.vbdry(0))->vloadbuff(boundary::manifolds,&gbl->vdt(0)(0,0),0,3,0);
		x.vbdry(base.vbdry(1))->vloadbuff(boundary::manifolds,&gbl->vdt(base.nseg)(0,0),0,3,0);
		x.vbdry(base.vbdry(0))->comm_prepare(boundary::manifolds,mp_phase,boundary::symmetric);
		x.vbdry(base.vbdry(1))->comm_prepare(boundary::manifolds,mp_phase,boundary::symmetric);

		x.vbdry(base.vbdry(0))->comm_exchange(boundary::manifolds,mp_phase,boundary::symmetric);
		x.vbdry(base.vbdry(1))->comm_exchange(boundary::manifolds,mp_phase,boundary::symmetric);        

		last_phase = true;
		last_phase &= x.vbdry(base.vbdry(0))->comm_wait(boundary::manifolds,mp_phase,boundary::symmetric);
		last_phase &= x.vbdry(base.vbdry(1))->comm_wait(boundary::manifolds,mp_phase,boundary::symmetric);
		x.vbdry(base.vbdry(0))->vfinalrcv(boundary::manifolds,mp_phase,boundary::symmetric,boundary::average,&gbl->vdt(0)(0,0),0,3,0);
		x.vbdry(base.vbdry(1))->vfinalrcv(boundary::manifolds,mp_phase,boundary::symmetric,boundary::average,&gbl->vdt(base.nseg)(0,0),0,3,0);
	}

	if (gbl->is_loop) {
		for(int m=0;m<tri_mesh::ND;++m)
			for(int n=0;n<tri_mesh::ND;++n)
				gbl->vdt(0)(m,n) = 0.5*(gbl->vdt(0)(m,n) +gbl->vdt(base.nseg+1)(m,n));
		gbl->vdt(base.nseg+1) = gbl->vdt(0);
	}

	FLT jcbi,temp;
	for(indx=0;indx<base.nseg+1;++indx) {    
		/* INVERT VERTEX MATRIX */
		jcbi = 1.0/(gbl->vdt(indx)(0,0)*gbl->vdt(indx)(1,1) -gbl->vdt(indx)(0,1)*gbl->vdt(indx)(1,0));

		temp = gbl->vdt(indx)(0,0)*jcbi*gbl->cfl(1,x.log2p);
		gbl->vdt(indx)(0,0) = gbl->vdt(indx)(1,1)*jcbi*gbl->cfl(0,x.log2p);
		gbl->vdt(indx)(1,1) = temp;
		gbl->vdt(indx)(0,1) *= -jcbi*gbl->cfl(1,x.log2p);
		gbl->vdt(indx)(1,0) *= -jcbi*gbl->cfl(0,x.log2p);
	}

	/* INVERT SIDE MATRIX */    
	if (basis::tri(x.log2p)->sm() > 0) {
		for(indx=0;indx<base.nseg;++indx) {
			/* INVERT SIDE MVDT MATRIX */
			jcbi = 1.0/(gbl->sdt(indx)(0,0)*gbl->sdt(indx)(1,1) -gbl->sdt(indx)(0,1)*gbl->sdt(indx)(1,0));

			temp = gbl->sdt(indx)(0,0)*jcbi*gbl->cfl(1,x.log2p);
			gbl->sdt(indx)(0,0) = gbl->sdt(indx)(1,1)*jcbi*gbl->cfl(0,x.log2p);
			gbl->sdt(indx)(1,1) = temp;
			gbl->sdt(indx)(0,1) *= -jcbi*gbl->cfl(1,x.log2p);
			gbl->sdt(indx)(1,0) *= -jcbi*gbl->cfl(0,x.log2p);
		}
	}
	return;
}

#ifdef petsc
void melt_end_pt::petsc_jacobian_dirichlet() {
	 	 
	int indx;
	if (surfbdry == 0) {
	 indx = x.ebdry(base.ebdry(0))->nseg;
	}
	else {
	 indx = 0;
	}

	/* GET X & Y MESH MOVEMENT ROW */
	/* CONSTRAIN MOTION NORMAL TO BOUNDARY */
	int row = (x.NV+tri_mesh::ND)*base.pnt +x.NV; 
	int nnz1 = x.J._cpt(row+1) -x.J._cpt(row);
	int nnz2 = x.J._cpt(row+2) -x.J._cpt(row+1);
	Array<int,1> cols(nnz1);
	Array<FLT,2> vals(2,nnz1);

	/* SOME ERROR CHECKING TO MAKE SURE ROW SPARSENESS PATTERN IS THE SAME */
	if (nnz1 != nnz2) {
	 *x.gbl->log << "zeros problem in deforming mesh on angled boundary\n";
	 sim::abort(__LINE__,__FILE__,x.gbl->log);
	}
	int row1 = x.J._cpt(row);
	int row2 = x.J._cpt(row+1);
	for(int col=0;col<nnz1;++col) {
	 if (x.J._col(row1++) != x.J._col(row2++)) {
		 *x.gbl->log << "zeros indexing problem in deforming mesh on angled boundary\n";
		 sim::abort(__LINE__,__FILE__,x.gbl->log);
	 }	
	}

	/* Get back unrotated equation for normal direction */
	FLT J = surf->gbl->vdt(indx)(0,0)*surf->gbl->vdt(indx)(1,1) -surf->gbl->vdt(indx)(1,0)*surf->gbl->vdt(indx)(0,1);
	vals(0,Range(0,nnz1-1)) = -x.J._val(Range(x.J._cpt(row),x.J._cpt(row+1)-1))*surf->gbl->vdt(indx)(1,0)/J;
	vals(0,Range(0,nnz1-1)) += x.J._val(Range(x.J._cpt(row+1),x.J._cpt(row+2)-1))*surf->gbl->vdt(indx)(0,0)/J;

	/* Replace x equation with tangential position equation */
	/* Replacy y equation with normal displacement equation */
	/* Normal Equation */
	vals(1,Range::all()) = 0.0;
	for(int col=0;col<nnz1;++col) {
	 if (x.J._col(row1+col) == row) {
		 vals(1,col) = wall_normal(0);
		 break;
	 }
	}
	for(int col=0;col<nnz1;++col) {
	 if (x.J._col(row1+col) == row+1) {
		 vals(1,col) = wall_normal(1);
		 break;
	 }
	}
	 
	/* tangent = -sin(theta) i +cos(theta) j */
	/* normal = cos(theta) i + sin(theta) j */
	/* Rotate equations for diagonal dominance to match what is done to residual */
	Array<FLT,2> temp(2,nnz1);
	temp(0,Range::all()) =  vals(0,Range::all())*surf->gbl->vdt(indx)(0,1) +vals(1,Range::all())*wall_normal(0);
	temp(1,Range::all()) =  vals(0,Range::all())*surf->gbl->vdt(indx)(1,1) +vals(1,Range::all())*wall_normal(1);

	TinyVector<int,2> rows(row,row+1);
	x.J._val(Range(x.J._cpt(row),x.J._cpt(row+1)-1)) = temp(0,Range::all());
	x.J._val(Range(x.J._cpt(row+1),x.J._cpt(row+2)-1)) = temp(1,Range::all());	 
}
#endif

					 
void kellerman::init(input_map& inmap,void* gbl_in) {
	melt::init(inmap,gbl_in);
	
	/* Set u,v fixed, T & P to be free */
	essential_indices.clear();
	for(int n=0;n<x.ND;++n)
		essential_indices.push_back(n);
	type = natural;
	type(Range(0,x.ND-1)) = essential;

	return;
}

#ifdef petsc			
void kellerman::petsc_jacobian_dirichlet() {
	melt::petsc_jacobian_dirichlet();
	
	int sm=basis::tri(x.log2p)->sm();
	Array<int,1> indices((base.nseg+1)*x.ND +base.nseg*sm*x.ND);

	int vdofs = x.NV +tri_mesh::ND;
	
	int gind,v0,sind;
	int counter = 0;
	
	int j = 0;
	do {
		sind = base.seg(j);
		v0 = x.seg(sind).pnt(0);
		gind = v0*vdofs;
		for(int n=x.NV;n < vdofs;++n) {
			indices(counter++)= gind +n;
		}
	} while (++j < base.nseg);
	v0 = x.seg(sind).pnt(1);
	gind = v0*vdofs;
	for(int n=x.NV;n < vdofs;++n) {
		indices(counter++)= gind +n;
	}
	
	
	gind = jacobian_start;
	for(int i=0;i<base.nseg;++i) {
		for(int m=0; m<sm; ++m) {
			for(int n=0;n < x.ND;++n) {
				indices(counter++) = gind++;
			}
		}
	}	
	
#ifdef MY_SPARSE
	x.J.zero_rows(counter,indices);
	x.J_mpi.zero_rows(counter,indices);
	x.J.set_diag(counter,indices,1.0);
#else
	MatZeroRows(x.petsc_J,counter,indices.data(),1.0);
#endif
}
#endif
