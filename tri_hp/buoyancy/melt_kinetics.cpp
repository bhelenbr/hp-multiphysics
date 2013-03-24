#include "bdry_buoyancy.h"
#include <myblas.h>

//#define MPDEBUG

//#define DEBUG

using namespace bdry_buoyancy;

void melt_kinetics::init(input_map& inmap,void* gbl_in) {
	melt::init(inmap,gbl_in);
	
	
	/* Don't make Temperature a dirichlet B.C. */
	essential_indices.clear();
	for(int n=0;n<tri_mesh::ND;++n)
		essential_indices.push_back(n);
	type = natural;
	type(Range(0,tri_mesh::ND-1)) = essential;
	
	gbl = static_cast<global *>(gbl_in);

	inmap.getwdefault(base.idprefix + "_Krough",gbl->Krough,0.0);
	inmap.getwdefault(base.idprefix + "_Kgt",gbl->Kgt,0.0);
	inmap.getwdefault(base.idprefix + "_Ksn",gbl->Ksn,1.0);
	inmap.getwdefault(base.idprefix + "_K2Dn",gbl->K2Dn,3.0);
	inmap.getwdefault(base.idprefix + "_K2Dn_max",gbl->K2Dn_max,3.0);
	inmap.getwdefault(base.idprefix + "_A2Dn",gbl->A2Dn,1.0);
	FLT angle;
	inmap.getwdefault(base.idprefix + "_facet_angle",angle,0.0);
	gbl->facetdir(0) = cos(M_PI*angle/180.0);
	gbl->facetdir(1) = sin(M_PI*angle/180.0);
	
	gbl->vdt_kinetic.resize(base.maxseg+1);
	gbl->sdt_kinetic.resize(base.maxseg);
	
#ifdef MELT1
	*x.gbl->log << "Can't use kinetic stuff with MELT1 defined" << std::endl;
	sim::abort(__LINE__,__FILE__,x.gbl->log);
#endif
	
	return;
}

void melt_kinetics::element_rsdl(int indx, Array<TinyVector<FLT,MXTM>,1> lf) {
	int i,n,sind,v0,v1;
	TinyVector<FLT,tri_mesh::ND> norm, rp, aloc;
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
		
		Array<FLT,1> au(x.NV), axpt(tri_mesh::ND), amv(tri_mesh::ND), anorm(tri_mesh::ND);
		for(n=0;n<x.NV;++n)
			au(n) = u(n)(i);
		axpt(0) = crd(0,i); axpt(1) = crd(1,i);
		aloc(0) = crd(0,i) ; aloc(1) = crd(1,i);
		
		amv(0) = (x.gbl->bd(0)*(crd(0,i) -dxdt(x.log2p,indx)(0,i))); amv(1) = (x.gbl->bd(0)*(crd(1,i) -dxdt(x.log2p,indx)(1,i)));
#ifdef MESH_REF_VEL
		amv(0) += x.gbl->mesh_ref_vel(0);
		amv(1) += x.gbl->mesh_ref_vel(1);
#endif
		anorm(0)= norm(0)/jcb; anorm(1) = norm(1)/jcb;
		
		/* This is from Weinstein's expression */
		FLT sint = -gbl->facetdir(0)*anorm(1) +gbl->facetdir(1)*anorm(0);
		FLT DT = ibc->f(2, aloc, x.gbl->time) -u(2)(i);
		FLT K = calculate_kinetic_coefficients(DT,sint);

		/* TANGENTIAL SPACING */
		res(0,i) = -ksprg(indx)*jcb;
		/* NORMAL FLUX */
		res(1,i) = RAD(crd(0,i))*x.gbl->rho*(mvel(0,i)*norm(0) +mvel(1,i)*norm(1));     
		/* Kinetic equation for surface temperature */
		res(2,i) = RAD(crd(0,i))*x.gbl->rho*(-DT)*jcb +K*res(1,i);  // -gbl->K_gt*kappa?;
		/* Latent Heat source term and additional heat flux */
		res(3,i) = RAD(crd(0,i))*fluxes(2).Eval(au,axpt,amv,anorm,x.gbl->time)*jcb -gbl->Lf*res(1,i) +gbl->rho_s*gbl->cp_s*u(2)(i)*res(1,i)/x.gbl->rho;
	}
	lf = 0.0;
	
	/* INTEGRATE & STORE MESH MOVEMENT RESIDUALS */
	basis::tri(x.log2p)->intgrt1d(&lf(x.NV-2)(0),&res(3,0)); // heat flux
	basis::tri(x.log2p)->intgrt1d(&lf(x.NV-1)(0),&res(1,0)); // mass flux
	basis::tri(x.log2p)->intgrtx1d(&lf(x.NV)(0),&res(0,0)); // tangent
#ifdef petsc
	basis::tri(x.log2p)->intgrt1d(&lf(x.NV+1)(0),&res(2,0)); // kinetic equation
#else
	basis::tri(x.log2p)->intgrt1d(&lf(x.NV+2)(0),&res(2,0)); // kinetic equation
#endif
								
	return;
}
	
void melt_kinetics::vdirichlet() {
	int i,sind,v0;
	
	/* Put kinetic equation into slot 2 because heat equation is going into slot 1 */
	/* Do it this way to make it easier to calculate jacobian for petsc*/
#ifdef petsc
	for(int i=0;i<base.nseg;++i) {
		gbl->vres(i)(2) = gbl->vres(i)(1);
		for(int m=0;m<basis::tri(x.log2p)->sm();++m) {
					gbl->sres(i,m)(2) = gbl->sres(i,m)(1);
		}
	}
	gbl->vres(base.nseg)(2) = gbl->vres(base.nseg)(1);
#endif
	
	melt::vdirichlet();
	
#ifdef petsc
	/* Put kinetic equation residual in temperature slot */
	i = 0;
	do {
		sind = base.seg(i);
		v0 = x.seg(sind).pnt(0);
		x.gbl->res.v(v0,2) = gbl->vres(i)(2);
		for(int m=0;m<basis::tri(x.log2p)->sm();++m)
			x.gbl->res.s(sind,m,2) = gbl->sres(i,m)(2);
	} while (++i < base.nseg);
	v0 = x.seg(sind).pnt(1);
	x.gbl->res.v(v0,2) = gbl->vres(base.nseg)(2);
#else
	/* Put zero in temperature slot */
	i = 0;
	do {
		sind = base.seg(i);
		v0 = x.seg(sind).pnt(0);
		x.gbl->res.v(v0,2) = 0.0;
		for(int m=0;m<basis::tri(x.log2p)->sm();++m)
			x.gbl->res.s(sind,m,2) = 0.0;
	} while (++i < base.nseg);
	v0 = x.seg(sind).pnt(1);
	x.gbl->res.v(v0,2) = 0.0;
#endif
	
}

void melt_kinetics::minvrt() {
	int i,m,n,indx;
	int last_phase, mp_phase;
	FLT temp;
	
	/* INVERT MASS MATRIX */
	/* LOOP THROUGH SIDES */
	if (basis::tri(x.log2p)->sm() > 0) {
		for(indx = 0; indx<base.nseg; ++indx) {
			/* SUBTRACT SIDE CONTRIBUTIONS TO VERTICES */
			for (m=0; m <basis::tri(x.log2p)->sm(); ++m) {
				for(n=0;n<neq;++n)
					gbl->vres(indx)(n) -= basis::tri(x.log2p)->sfmv1d(0,m)*gbl->sres(indx,m)(n);
				for(n=0;n<neq;++n)
					gbl->vres(indx+1)(n) -= basis::tri(x.log2p)->sfmv1d(1,m)*gbl->sres(indx,m)(n);
			}
		}
	}
	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		x.vbdry(base.vbdry(0))->vloadbuff(boundary::manifolds,&gbl->vres(0)(0),0,neq-1,0);
		x.vbdry(base.vbdry(1))->vloadbuff(boundary::manifolds,&gbl->vres(base.nseg)(0),0,neq-1,0);
		x.vbdry(base.vbdry(0))->comm_prepare(boundary::manifolds,mp_phase,boundary::symmetric);
		x.vbdry(base.vbdry(1))->comm_prepare(boundary::manifolds,mp_phase,boundary::symmetric);
		
		x.vbdry(base.vbdry(0))->comm_exchange(boundary::manifolds,mp_phase,boundary::symmetric);
		x.vbdry(base.vbdry(1))->comm_exchange(boundary::manifolds,mp_phase,boundary::symmetric);
		
		last_phase = true;
		
		last_phase &= x.vbdry(base.vbdry(0))->comm_wait(boundary::manifolds,mp_phase,boundary::symmetric);
		last_phase &= x.vbdry(base.vbdry(1))->comm_wait(boundary::manifolds,mp_phase,boundary::symmetric);
		
		x.vbdry(base.vbdry(0))->vfinalrcv(boundary::manifolds,mp_phase,boundary::symmetric,boundary::average,&gbl->vres(0)(0),0,neq-1,0);
		x.vbdry(base.vbdry(1))->vfinalrcv(boundary::manifolds,mp_phase,boundary::symmetric,boundary::average,&gbl->vres(base.nseg)(0),0,neq-1,0);
	}
	
	if (gbl->is_loop) {
		for(int n=1;n<neq;++n) {
			gbl->vres(0)(n) = 0.5*(gbl->vres(0)(n) +gbl->vres(base.nseg+1)(n));
			gbl->vres(base.nseg+1)(n) = gbl->vres(0)(n);
		}
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
		gbl->vres(i)(2) *= gbl->vdt_kinetic(i);  // change here from melt
	}
	
	/* SOLVE FOR SIDE MODES */
	if (basis::tri(x.log2p)->sm() > 0) {
		for(indx = 0; indx<base.nseg; ++indx) {
			/* INVERT SIDE MODES */
			DPBTRSNU2((double *) &basis::tri(x.log2p)->sdiag1d(0,0),basis::tri(x.log2p)->sbwth()+1,basis::tri(x.log2p)->sm(),basis::tri(x.log2p)->sbwth(),&(gbl->sres(indx,0)(0)),neq);
			for(m=0;m<basis::tri(x.log2p)->sm();++m) {
				temp                             = gbl->sres(indx,m)(0)*gbl->sdt(indx)(0,0) +gbl->sres(indx,m)(1)*gbl->sdt(indx)(0,1);
				gbl->sres(indx,m)(1) = gbl->sres(indx,m)(0)*gbl->sdt(indx)(1,0) +gbl->sres(indx,m)(1)*gbl->sdt(indx)(1,1);
				gbl->sres(indx,m)(0) = temp;
				
				gbl->sres(indx,m)(2) *= gbl->sdt_kinetic(indx); // change here from melt
			}
			
			for(m=0;m<basis::tri(x.log2p)->sm();++m) {
				for(n=0;n<neq;++n) {
					gbl->sres(indx,m)(n) -= basis::tri(x.log2p)->vfms1d(0,m)*gbl->vres(indx)(n);
					gbl->sres(indx,m)(n) -= basis::tri(x.log2p)->vfms1d(1,m)*gbl->vres(indx+1)(n);
				}
			}
		}
	}
	
	return;
}

void melt_kinetics::update(int stage) {
	int i,m,n,count,sind,indx,v0;
	
	if (stage < 0) {
		indx = 0;
		int i = 0;
		do {
			sind = base.seg(i);
			v0 = x.seg(sind).pnt(0);
			for(n=0;n<tri_mesh::ND;++n)
				gbl->vug0(i)(n) = x.pnts(v0)(n);
			gbl->vug0(i)(2) = x.ug.v(v0,2);
		} while (++i < base.nseg);
		v0 = x.seg(sind).pnt(1);
		for(n=0;n<tri_mesh::ND;++n)
			gbl->vug0(base.nseg)(n) = x.pnts(v0)(n);
		gbl->vug0(base.nseg)(2) = x.ug.v(v0,2);
		
		if (basis::tri(x.log2p)->sm() > 0) {
			for (i=0;i<base.nseg;++i) {
				for(m=0;m<basis::tri(x.log2p)->sm();++m) {
					for(n=0;n<tri_mesh::ND;++n) {
						gbl->sug0(i,m)(n) = crv(i,m)(n);
					}
					sind = base.seg(i);
					gbl->sug0(i,m)(2) = x.ug.s(sind,m,2);
				}
			}
		}
		
		return;
	}
	
	minvrt();
	
#ifdef DEBUG
	//   if (x.coarse_flag) {
	for(i=0;i<base.nseg+1;++i)
		*x.gbl->log << "vdt: " << i << ' ' << gbl->vdt(i)(0,0) << ' ' << gbl->vdt(i)(0,1) << ' ' << gbl->vdt(i)(1,0) << ' ' << gbl->vdt(i)(1,1) << ' ' << gbl->vdt_kinetic(i) << '\n';
	
	for(i=0;i<base.nseg;++i)
		*x.gbl->log << "sdt: " << i << ' ' << gbl->sdt(i)(0,0) << ' ' << gbl->sdt(i)(0,1) << ' ' << gbl->sdt(i)(1,0) << ' ' << gbl->sdt(i)(1,1) << ' ' << gbl->sdt_kinetic(i) << '\n';
	
	for(i=0;i<base.nseg+1;++i) {
		*x.gbl->log << "vres: " << i << ' ';
		for(n=0;n<neq;++n) {
			if (fabs(gbl->vres(i)(n)) > 1.0e-9) *x.gbl->log << gbl->vres(i)(n) << ' ';
			else *x.gbl->log << "0.0 ";
		}
		*x.gbl->log << '\n';
	}
	
	for(i=0;i<base.nseg;++i) {
		for(m=0;m<basis::tri(x.log2p)->sm();++m) {
			*x.gbl->log << "sres: " << i << ' ';
			for(n=0;n<neq;++n) {
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
		*x.gbl->log << "vertex positions and T " << v0 << ' ' << x.pnts(v0)(0) << ' ' << x.pnts(v0)(1) << ' ' << x.ug.v(v0,2) << '\n';
	} while(++i < base.nseg);
	v0 = x.seg(sind).pnt(1);
	*x.gbl->log << "vertex positions and T " << v0 << ' ' << x.pnts(v0)(0) << ' ' << x.pnts(v0)(1) << ' ' << x.ug.v(v0,2) << '\n';
	
	for(i=0;i<base.nseg;++i)
		for(m=0;m<basis::tri(x.log2p)->sm();++m)
			*x.gbl->log << "spos: " << i << ' ' << m << ' ' << crv(i,m)(0) << ' ' << crv(i,m)(1) << ' ' << x.ug.s(base.seg(i),m,2) << '\n';
	// }
#endif
	
	i = 0;
	do {
		sind = base.seg(i);
		v0 = x.seg(sind).pnt(0);
		for(n=0;n<tri_mesh::ND;++n)
			x.pnts(v0)(n) = gbl->vug0(i)(n) -x.gbl->alpha(stage)*gbl->vres(i)(n);
		x.ug.v(v0,2) = gbl->vug0(i)(2) -x.gbl->alpha(stage)*gbl->vres(i)(2);
	} while (++i < base.nseg);
	v0 = x.seg(sind).pnt(1);
	for(n=0;n<tri_mesh::ND;++n)
		x.pnts(v0)(n) = gbl->vug0(base.nseg)(n) -x.gbl->alpha(stage)*gbl->vres(base.nseg)(n);
	x.ug.v(v0,2) = gbl->vug0(base.nseg)(2) -x.gbl->alpha(stage)*gbl->vres(base.nseg)(2);
		
	if (basis::tri(x.log2p)->sm() > 0) 
		for(i=0;i<base.nseg;++i)
			for(m=0;m<basis::tri(x.log2p)->sm();++m)
 				for(n=0;n<tri_mesh::ND;++n)
					crv(i,m)(n) = gbl->sug0(i,m)(n) -x.gbl->alpha(stage)*gbl->sres(i,m)(n);
	
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
		*x.gbl->log << "vertex positions & T " << v0 << ' ' << x.pnts(v0)(0) << ' ' << x.pnts(v0)(1) << ' ' << x.ug.v(v0,2) << '\n';
	} while(++i < base.nseg);
	v0 = x.seg(sind).pnt(1);
	*x.gbl->log << "vertex positions & T " << v0 << ' ' << x.pnts(v0)(0) << ' ' << x.pnts(v0)(1) << ' ' << x.ug.v(v0,2) << '\n';
	
	for(i=0;i<base.nseg;++i)
		for(m=0;m<basis::tri(x.log2p)->sm();++m)
			*x.gbl->log << "spos & T: " << i << ' ' << m << ' ' << crv(i,m)(0) << ' ' << crv(i,m)(1) << ' ' << x.ug.s(base.seg(i),m,2) << '\n';
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
			base.fsndbuf(count++) = x.ug.v(v0,2);
		} while(++i < base.nseg);
		v0 = x.seg(sind).pnt(1);
		for(n=0;n<tri_mesh::ND;++n)
			base.fsndbuf(count++) = x.pnts(v0)(n);
		base.fsndbuf(count++) = x.ug.v(v0,2);
		
		for(i=0;i<base.nseg;++i) {
			for(m=0;m<basis::tri(x.log2p)->sm();++m) {
				for(n=0;n<tri_mesh::ND;++n)
					base.fsndbuf(count++) = crv(i,m)(n);
				base.fsndbuf(count++) = x.ug.s(base.seg(i),m,2);
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


#ifdef petsc
void melt_kinetics::petsc_jacobian_dirichlet() {
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
	/* SWAP HEAT EQUATION WITH KINETIC EQUATION  */
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
			/* swap rows */
			if (x.J._col(cpt1) == x.J._col(cpt2)) {
				int tempInd = x.J._col(cpt2);
				FLT tempVal = x.J._val(cpt2);
				x.J._col(cpt2) = x.J._col(cpt1);
				x.J._val(cpt2) = x.J._val(cpt1);
				x.J._col(cpt1) = tempInd;
				x.J._val(cpt1) = tempVal;
			}
			else {
				x.J.set_values(row+3,x.J._col(cpt1),x.J._val(cpt1));
				x.J.set_values(row,x.J._col(cpt1),0.0);
			}
			++cpt1;
			++cpt2;
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
				/* swap rows */
				if (x.J._col(cpt1) == x.J._col(cpt2)) {
					int tempInd = x.J._col(cpt2);
					FLT tempVal = x.J._val(cpt2);
					x.J._col(cpt2) = x.J._col(cpt1);
					x.J._val(cpt2) = x.J._val(cpt1);
					x.J._col(cpt1) = tempInd;
					x.J._val(cpt1) = tempVal;
					
				}
				else {
					x.J.set_values(row2,x.J._col(cpt1),x.J._val(cpt1));
					x.J.set_values(row1,x.J._col(cpt1),0.0);
				}
				++cpt1;
				++cpt2;
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
		/* swap rows */
		if (x.J._col(cpt1) == x.J._col(cpt2)) {
			int tempInd = x.J._col(cpt2);
			FLT tempVal = x.J._val(cpt2);
			x.J._col(cpt2) = x.J._col(cpt1);
			x.J._val(cpt2) = x.J._val(cpt1);
			x.J._col(cpt1) = tempInd;
			x.J._val(cpt1) = tempVal;
		}
		else {
			x.J.set_values(row+3,x.J._col(cpt1),x.J._val(cpt1));
			x.J.set_values(row,x.J._col(cpt1),0.0);
		}
		++cpt1;
		++cpt2;
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
			/* swap rows */
			if (x.J_mpi._col(cpt1) == x.J_mpi._col(cpt2)) {
				int tempInd = x.J_mpi._col(cpt2);
				FLT tempVal = x.J_mpi._val(cpt2);
				x.J_mpi._col(cpt2) = x.J_mpi._col(cpt1);
				x.J_mpi._val(cpt2) = x.J_mpi._val(cpt1);
				x.J_mpi._col(cpt1) = tempInd;
				x.J_mpi._val(cpt1) = tempVal;
			}
			else {
				x.J_mpi.set_values(row+3,x.J_mpi._col(cpt1),x.J_mpi._val(cpt1));
				x.J_mpi.set_values(row,x.J_mpi._col(cpt1),0.0);
			}
			++cpt1;
			++cpt2;
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
				/* swap rows */
				if (x.J_mpi._col(cpt1) == x.J_mpi._col(cpt2)) {
					int tempInd = x.J_mpi._col(cpt2);
					FLT tempVal = x.J_mpi._val(cpt2);
					x.J_mpi._col(cpt2) = x.J_mpi._col(cpt1);
					x.J_mpi._val(cpt2) = x.J_mpi._val(cpt1);
					x.J_mpi._col(cpt1) = tempInd;
					x.J_mpi._val(cpt1) = tempVal;
				}
				else {
					x.J_mpi.set_values(row2,x.J_mpi._col(cpt1),x.J_mpi._val(cpt1));
					x.J_mpi.set_values(row1,x.J_mpi._col(cpt1),0.0);
				}
				++cpt1;
				++cpt2;
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
		/* swap rows */
		if (x.J_mpi._col(cpt1) == x.J_mpi._col(cpt2)) {
			int tempInd = x.J_mpi._col(cpt2);
			FLT tempVal = x.J_mpi._val(cpt2);
			x.J_mpi._col(cpt2) = x.J_mpi._col(cpt1);
			x.J_mpi._val(cpt2) = x.J_mpi._val(cpt1);
			x.J_mpi._col(cpt1) = tempInd;
			x.J_mpi._val(cpt1) = tempVal;
		}
		else {
			x.J_mpi.set_values(row+3,x.J_mpi._col(cpt1),x.J_mpi._val(cpt1));
			x.J_mpi.set_values(row,x.J_mpi._col(cpt1),0.0);
		}
		++cpt1;
		++cpt2;
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
	
	symbolic::petsc_jacobian_dirichlet();  // Sets u,v dirichlet condition

}
#endif


void melt_kinetics::setup_preconditioner() {
	
	melt::setup_preconditioner();
	
	
	int indx,sind,v0,v1;
	TinyVector<FLT,tri_mesh::ND> nrm;
	FLT h, dtnorm;
	TinyMatrix<FLT,tri_mesh::ND,MXGP> crd, dcrd;
	TinyMatrix<FLT,4,MXGP> res;
	TinyMatrix<FLT,4,MXGP> lf;
	TinyVector<FLT,2> mvel;
	Array<TinyVector<FLT,MXGP>,1> u(x.NV);
	int last_phase, mp_phase;
	TinyVector<FLT,tri_mesh::ND> aPoint;
	aPoint = 0.0;
	
	/**************************************************/
	/* DETERMINE SURFACE MOVEMENT TIME STEP              */
	/**************************************************/
	gbl->vdt_kinetic(0) = 0.0;
	
	for(indx=0; indx < base.nseg; ++indx) {
		sind = base.seg(indx);
		v0 = x.seg(sind).pnt(0);
		v1 = x.seg(sind).pnt(1);
		
		nrm(0) =  (x.pnts(v1)(1) -x.pnts(v0)(1));
		nrm(1) = -(x.pnts(v1)(0) -x.pnts(v0)(0));
		nrm *= 0.5;
		h = sqrt(nrm(0)*nrm(0) +nrm(1)*nrm(1));
		dtnorm = 1.0;
		dtnorm *= RAD(0.5*(x.pnts(v0)(0) +x.pnts(v1)(0)));

		gbl->vdt_kinetic(indx) += dtnorm*h*basis::tri(x.log2p)->vdiag1d();		
		if (basis::tri(x.log2p)->sm()) {
			gbl->sdt_kinetic(indx) = dtnorm*h;
		}
	}
	
	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		x.vbdry(base.vbdry(0))->vloadbuff(boundary::manifolds,&gbl->vdt_kinetic(0),0,0,0);
		x.vbdry(base.vbdry(1))->vloadbuff(boundary::manifolds,&gbl->vdt_kinetic(base.nseg),0,0,0);
		x.vbdry(base.vbdry(0))->comm_prepare(boundary::manifolds,mp_phase,boundary::symmetric);
		x.vbdry(base.vbdry(1))->comm_prepare(boundary::manifolds,mp_phase,boundary::symmetric);
		
		x.vbdry(base.vbdry(0))->comm_exchange(boundary::manifolds,mp_phase,boundary::symmetric);
		x.vbdry(base.vbdry(1))->comm_exchange(boundary::manifolds,mp_phase,boundary::symmetric);        
		
		last_phase = true;
		last_phase &= x.vbdry(base.vbdry(0))->comm_wait(boundary::manifolds,mp_phase,boundary::symmetric);
		last_phase &= x.vbdry(base.vbdry(1))->comm_wait(boundary::manifolds,mp_phase,boundary::symmetric);
		x.vbdry(base.vbdry(0))->vfinalrcv(boundary::manifolds,mp_phase,boundary::symmetric,boundary::average,&gbl->vdt_kinetic(0),0,0,0);
		x.vbdry(base.vbdry(1))->vfinalrcv(boundary::manifolds,mp_phase,boundary::symmetric,boundary::average,&gbl->vdt_kinetic(base.nseg),0,0,0);
	}
	
	if (gbl->is_loop) {
		gbl->vdt_kinetic(0) = 0.5*(gbl->vdt_kinetic(0) +gbl->vdt_kinetic(base.nseg+1));
		gbl->vdt_kinetic(base.nseg+1) = gbl->vdt_kinetic(0);
	}
	
	for(indx=0;indx<base.nseg+1;++indx) {    
		/* INVERT VERTEX MATRIX */
		gbl->vdt_kinetic(indx) = 1./gbl->vdt_kinetic(indx);
	}
	
	/* INVERT SIDE MATRIX */    
	if (basis::tri(x.log2p)->sm() > 0) {
		for(indx=0;indx<base.nseg;++indx) {
			/* INVERT SIDE MVDT MATRIX */
			gbl->sdt_kinetic(indx) = 1./gbl->sdt_kinetic(indx);
		}
	}
	return;
}

void melt_kinetics::output(std::ostream& fout, tri_hp::filetype typ,int tlvl) {
	melt::output(fout,typ,tlvl);
	
	TinyVector<FLT,tri_mesh::ND> norm, rp, aloc;
	Array<FLT,1> ubar(x.NV);
	FLT jcb;
	Array<TinyVector<FLT,MXGP>,1> u(x.NV);
	TinyMatrix<FLT,tri_mesh::ND,MXGP> crd, dcrd, mvel;
	
	switch (typ) {
		case(tri_hp::tecplot): {
			if (!report_flag) break;
			
			std::ostringstream fname;
			fname << "data" << x.gbl->tstep << "kinetics_" << base.idprefix << ".dat";
			std::ofstream fout;
			fout.open(fname.str().c_str());
			
			int indx = 0;
			FLT s = 0.0;
			do {
				int sind = base.seg(indx);

				x.crdtocht1d(sind);
				for(int n=0;n<tri_mesh::ND;++n)
					basis::tri(x.log2p)->proj1d(&x.cht(n,0),&crd(n,0),&dcrd(n,0));
				
				x.ugtouht1d(sind);
				for(int n=0;n<x.NV;++n)
					basis::tri(x.log2p)->proj1d(&x.uht(n)(0),&u(n)(0));
				
				for(int i=0;i<basis::tri(x.log2p)->gpx();++i) {
					norm(0) =  dcrd(1,i);
					norm(1) = -dcrd(0,i);
					jcb = sqrt(norm(0)*norm(0) +norm(1)*norm(1));
					s += basis::tri(x.log2p)->wtx(i)*jcb;
					
					/* RELATIVE VELOCITY STORED IN MVEL(N)*/
					for(int n=0;n<tri_mesh::ND;++n) {
						mvel(n,i) = u(n)(i) -(x.gbl->bd(0)*(crd(n,i) -dxdt(x.log2p,indx)(n,i)));
#ifdef MESH_REF_VEL
						mvel(n,i) -= x.gbl->mesh_ref_vel(n);
#endif
					}
					
					Array<FLT,1> au(x.NV), axpt(tri_mesh::ND), amv(tri_mesh::ND), anorm(tri_mesh::ND);
					for(int n=0;n<x.NV;++n)
						au(n) = u(n)(i);
					aloc(0) = crd(0,i); aloc(1) = crd(1,i);
					
					amv(0) = (x.gbl->bd(0)*(crd(0,i) -dxdt(x.log2p,indx)(0,i))); amv(1) = (x.gbl->bd(0)*(crd(1,i) -dxdt(x.log2p,indx)(1,i)));
#ifdef MESH_REF_VEL
					amv(0) += x.gbl->mesh_ref_vel(0);
					amv(1) += x.gbl->mesh_ref_vel(1);
#endif
					anorm(0)= norm(0)/jcb; anorm(1) = norm(1)/jcb;
					
					FLT sint = -gbl->facetdir(0)*anorm(1) +gbl->facetdir(1)*anorm(0);
					FLT DT = ibc->f(2, aloc, x.gbl->time) -u(2)(i);
					FLT K = calculate_kinetic_coefficients(DT,sint);

					fout << s << ' ' << crd(0,i) << ' ' << crd(1,i) << ' ' << K << ' ' << DT << ' ' << sint << std::endl;

				}
			} while (++indx < base.nseg);
			fout.close();
		
			break;
		}
		default:
			break;
	}
}

FLT melt_kinetics::calculate_kinetic_coefficients(FLT DT,FLT sint) {

	// K2Dn is the inverse of B and is a ratio relative to Krough
	// FLT K2Dn_exp = gbl->K2Dn_max*gbl->K2Dn/(exp(-gbl->A2Dn/(abs(DT) +100.*EPSILON))*gbl->K2Dn_max +gbl->K2Dn);
	FLT K2Dn_exp = 1./(exp(-gbl->A2Dn/(abs(DT) +100.*EPSILON))/gbl->K2Dn +1./gbl->K2Dn_max);
	
	
	//  K2Dn and Ksn are ratios relative to Krough
	// FLT K = gbl->Krough*(1. +gbl->Ksn*K2Dn_exp/sqrt(pow(K2Dn_exp*sint,2) +pow(gbl->Ksn,2)));
	FLT K = gbl->Krough*(1. + 1./sqrt(pow(sint/gbl->Ksn,2) +pow(1./K2Dn_exp,2)));

	return(K);
}

void melt_facet_pt::rsdl(int stage) {

	if (surfbdry == 0) {
		const int bnum = base.ebdry(0);
		const int seg = x.ebdry(bnum)->nseg-1;
		const int sind = x.ebdry(bnum)->seg(seg);
		
		x.crdtocht1d(sind);
		x.ugtouht1d(sind);
		
		element_rsdl();
		
		/* SET TANGENT RESIDUAL TO ZERO */
		surf->gbl->vres(seg+1)(0) = 0.0;
#ifdef petsc
		surf->gbl->vres(seg+1)(1) = res;  // Kinetic Equation
#else
		surf->gbl->vres(seg+1)(2) = res;  // Kinetic Equation
#endif
	}
	else {
		*x.gbl->log << "melt facet_pt only implemented for end point of surface" << std::endl;
		sim::abort(__LINE__,__FILE__,x.gbl->log);
	}
	return;
}

void melt_facet_pt::element_rsdl() {	
	Array<FLT,1> u(x.NV);
	TinyVector<FLT, 2> xp, dxpdpsi, mvel;
	
	/* Calculate Temperature from velocity of point & 2D nucleation coefficient only */
	/* Assuming growth at facet angle */
	int v0 = base.pnt;
	basis::tri(x.log2p)->ptprobe1d(2,xp.data(),dxpdpsi.data(),1.0,&x.cht(0,0),MXTM);
	basis::tri(x.log2p)->ptprobe1d(4,u.data(),1.0,&x.uht(0)(0),MXTM);
	/* RELATIVE VELOCITY STORED IN MVEL(N)*/
	for(int n=0;n<tri_mesh::ND;++n) {
		mvel(n) = u(n) -x.gbl->bd(0)*(xp(n) -x.vrtxbd(1)(v0)(n));
#ifdef MESH_REF_VEL
		mvel(n) -= x.gbl->mesh_ref_vel(n);
#endif
	}
	
	FLT jcb = sqrt(dxpdpsi(0)*dxpdpsi(0)+dxpdpsi(1)*dxpdpsi(1));
	FLT DT = surf->ibc->f(2, xp, x.gbl->time) -u(2);
	melt_kinetics *surf1 = dynamic_cast<melt_kinetics *>(surf);
	FLT K = surf1->calculate_kinetic_coefficients(DT,0.0);
	FLT res1 = jcb*RAD(xp(0))*x.gbl->rho*(mvel(0)*surf1->gbl->facetdir(0) +mvel(1)*surf1->gbl->facetdir(1));
	/* Kinetic equation for surface temperature */
	res = RAD(xp(0))*x.gbl->rho*(-DT)*jcb +res1*K;  // -gbl->K_gt*kappa?;
}


#ifdef petsc
void melt_facet_pt::petsc_jacobian() {
#ifdef BZ_DEBUG
	const FLT eps_r = 0.0e-6, eps_a = 1.0e-6;  /*<< constants for debugging jacobians */
#else
	const FLT eps_r = 1.0e-6, eps_a = 1.0e-10;  /*<< constants for accurate numerical determination of jacobians */
#endif
	const int sm = basis::tri(x.log2p)->sm();
	const int vdofs = x.NV+x.ND;
	int nvars = vdofs*(sm+2); 
	Array<FLT,1> vals(nvars);
	Array<int,1> col(nvars);
	
	/* Jacobian for row determine x position */
	/* Variables effecting dT/dx are strictly those on the edge */
	
	const int bnum = base.ebdry(0);
	const int seg = x.ebdry(bnum)->nseg-1;
	const int sind = x.ebdry(bnum)->seg(seg);
	
	x.crdtocht1d(sind);
	x.ugtouht1d(sind);
	element_rsdl();
	FLT res0 = res;
	
	FLT dw = 0.0;
	for(int i=0;i<2;++i)
		for(int n=0;n<x.NV;++n)
			dw = dw + fabs(x.uht(n)(i));
		
	dw = dw*eps_r;
	dw += eps_a;
	FLT dx = eps_r*x.distance(x.seg(sind).pnt(0),x.seg(sind).pnt(1)) +eps_a;
		
	/* Numerically create Jacobian */
	int kcol = 0;
	for(int mode = 0; mode < sm+2; ++mode) {
		for(int var = 0; var < x.NV; ++var) {
			x.uht(var)(mode) += dw;
			element_rsdl();
			vals(kcol++) = (res-res0)/dw;
			x.uht(var)(mode) -= dw;
		}
		
		for(int var = 0; var < tri_mesh::ND; ++var){
			x.cht(var,mode) += dx;
			element_rsdl();
			vals(kcol++) = (res-res0)/dx;
			x.cht(var,mode) -= dx;
		}		
	}
	
	/* CREATE GLOBAL NUMBERING LIST */
	int ind = 0;
	for(int mode = 0; mode < 2; ++mode) {
		int gindx = vdofs*x.seg(sind).pnt(mode);
		for(int var = 0; var < vdofs; ++var)
			col(ind++) = gindx++;
	}
	
	int gindxNV = x.npnt*vdofs +x.NV*sind*sm;
	int gindxND = surf->jacobian_start +seg*tri_mesh::ND*sm;
	for(int mode = 0; mode < sm; ++mode) {
		for(int var = 0; var < x.NV; ++var)
			col(ind++) = gindxNV++;
		
		for(int var = 0; var < tri_mesh::ND; ++var)
			col(ind++) = gindxND++;
	}
	
#ifdef MY_SPARSE
	Array<int,1> rows(1);
	rows(0) = base.pnt*vdofs+vdofs-1;  // y-row is kinetic equation.  gets moved by melt_kinetic::petsc_jacobian_dirichlet
	x.J.zero_rows(1,rows);
	x.J_mpi.zero_rows(1,rows);
	x.J.add_values(rows(0),nvars,col,vals);
#else
	*x.gbl->log << "petsc not working for melt_end_pt\n";
	MatSetValuesLocal(x.petsc_J,x.NV*(sm+2),rows.data(),vdofs*2 +x.NV*sm,cols.data(),K.data(),ADD_VALUES);
#endif
}
#endif


