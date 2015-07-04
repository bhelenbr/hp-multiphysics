#include "bdry_ins.h"
#include <myblas.h>
//#include <blitz/tinyvec-et.h>

//#define MPDEBUG

//#define DEBUG

//#define BODYFORCE


using namespace bdry_ins;

// extern FLT body[ND];

void surface_slave::init(input_map& inmap,void* gbl_in) {
	std::string keyword;
	
	keyword = base.idprefix + "_curved";
	inmap[keyword] = "1";
	
	keyword = base.idprefix + "_coupled";
	inmap[keyword] = "1";

	neumann::init(inmap,gbl_in);
	
	return;
}

void surface::init(input_map& inmap,void* gin) {
	std::string keyword,val;
	std::istringstream data;
	std::string filename;

	surface_slave::init(inmap,gin);
	
	gbl = static_cast<global *>(gin);
	ksprg.resize(base.maxseg);
	vdres.resize(x.log2pmax,base.maxseg);
	sdres.resize(x.log2pmax,base.maxseg,x.sm0);

	keyword = base.idprefix + "_sigma";
	inmap.getwdefault(keyword,gbl->sigma,0.0);
		
	keyword = base.idprefix + "_matching_block";
	if (!inmap.get(keyword,val)) {
		gbl->mu2 = 0.0;
		gbl->rho2 = 0.0;
	}
	else {
		keyword = val +"_mu";
		if (!inmap.get(keyword,gbl->mu2)) {
			*x.gbl->log << "couldn't find matching blocks viscosity" << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
		
		keyword = val +"_rho";
		if (!inmap.get(keyword,gbl->rho2)) {
			*x.gbl->log << "couldn't find matching blocks density" << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
	}
	
		
	if (x.seg(base.seg(0)).pnt(0) == x.seg(base.seg(base.nseg-1)).pnt(1)) gbl->is_loop = true;
	else gbl->is_loop = false;

	gbl->vug0.resize(base.maxseg+1);
	gbl->sug0.resize(base.maxseg,x.sm0);
	
	gbl->vres.resize(base.maxseg+1);
	gbl->sres.resize(base.maxseg,x.sm0); 
	gbl->vres0.resize(base.maxseg+1);
	gbl->sres0.resize(base.maxseg,x.sm0); 
	
#ifdef DROP        
	gbl->vvolumeflux.resize(base.maxseg+1);
	gbl->svolumeflux.resize(base.maxseg,x.sm0);
#endif

			
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

void surface_slave::tadvance() {
	/* THIS IS TO RECEIVE MESSAGES SENT FROM SURFACE DURING TADVANCE RSDL CALL */
	if (x.gbl->substep == 0) {
		if (!x.coarse_level) {
			rsdl(gbl->nstage);
		}
	}
	return;
}

void surface::tadvance() {
	int i,j,m,sind;    
	
	hp_edge_bdry::tadvance();

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
				for(m=0;m<basis::tri(x.log2p).sm;++m)
						sdres(x.log2p,i,m)(0) = 0.0;

			rsdl(gbl->nstage);

			for(i=0;i<base.nseg+1;++i)
				vdres(x.log2p,i)(0) = -gbl->vres(i)(0);

			for(i=0;i<base.nseg;++i) 
				for(m=0;m<basis::tri(x.log2p).sm;++m)
						sdres(x.log2p,i,m)(0) = -gbl->sres(i,m)(0)*0;  /* TEMPO TO KEEP SIDE MODES EQUALLY SPACED */
		}
	}
	return;
}

void surface::rsdl(int stage) {
	int i,j,m,n,sind,indx,count,v0,v1;
	TinyVector<FLT,tri_mesh::ND> norm, rp;
	Array<FLT,1> ubar(x.NV);
	FLT jcb;
	Array<TinyVector<FLT,MXGP>,1> u(x.NV);
	TinyMatrix<FLT,tri_mesh::ND,MXGP> crd, dcrd, mvel;
	TinyMatrix<FLT,8,MXGP> res;
	TinyMatrix<FLT,4,MXGP> lf,lf1;  
	
	/**************************************************/
	/* DETERMINE MESH RESIDUALS & SURFACE TENSION      */
	/**************************************************/
	for(n=0;n<tri_mesh::ND;++n)
		gbl->vres(0)(n) = 0.0;
		
#ifdef DROP
	gbl->vvolumeflux(0) = 0.0;
#endif

		for(indx=0;indx<base.nseg;++indx) {
		sind = base.seg(indx);
		v0 = x.seg(sind).pnt(0);
		v1 = x.seg(sind).pnt(1);
		
		x.crdtocht1d(sind);
		for(n=0;n<tri_mesh::ND;++n)
			basis::tri(x.log2p).proj1d(&x.cht(n,0),&crd(n,0),&dcrd(n,0));

		x.ugtouht1d(sind);
		for(n=0;n<tri_mesh::ND;++n)
			basis::tri(x.log2p).proj1d(&x.uht(n)(0),&u(n)(0));    
		
		for(i=0;i<basis::tri(x.log2p).gpx;++i) {
			norm(0) =  dcrd(1,i);
			norm(1) = -dcrd(0,i);
			jcb = sqrt(norm(0)*norm(0) +norm(1)*norm(1));

			/* RELATIVE VELOCITY STORED IN MVEL(N)*/
			for(n=0;n<tri_mesh::ND;++n) {
				mvel(n,i) = u(n)(i) -(x.gbl->bd(0)*(crd(n,i) -dxdt(x.log2p,indx)(n,i)));
#ifdef DROP
				mvel(n,i) -= tri_hp_ins::mesh_ref_vel(n);
#endif    
			}
			/* TANGENTIAL SPACING */                
			res(0,i) = -ksprg(indx)*jcb;
			/* NORMAL FLUX */
			res(1,i) = -RAD(crd(0,i))*(mvel(0,i)*norm(0) +mvel(1,i)*norm(1));     
			/* UPWINDING BASED ON TANGENTIAL VELOCITY */
			res(2,i) = -res(1,i)*(-norm(1)*mvel(0,i) +norm(0)*mvel(1,i))/jcb*gbl->meshc(indx);
#ifdef DROP
			res(3,i) = +RAD(crd(0,i))*gbl->vflux*jcb;
#endif 

			/* SURFACE TENSION SOURCE TERM X-DIRECTION */ 
			res(4,i) = +RAD(crd(0,i))*(x.gbl->rho -gbl->rho2)*x.gbl->g*crd(1,i)*norm(0);
#ifdef AXISYMMETRIC
			res(4,i) += gbl->sigma*jcb;
#endif
			/* AND INTEGRATION BY PARTS TERM */
			res(5,i) = +RAD(crd(0,i))*gbl->sigma*norm(1)/jcb;


			/* SURFACE TENSION SOURCE TERM Y-DIRECTION */
			res(6,i) = +RAD(crd(0,i))*(x.gbl->rho -gbl->rho2)*x.gbl->g*crd(1,i)*norm(1);                
			/* AND INTEGRATION BY PARTS TERM */
			res(7,i) = -RAD(crd(0,i))*gbl->sigma*norm(0)/jcb;
		}
		
		for(m=0;m<basis::tri(x.log2p).sm+2;++m)
			lf(0,m) = 0.0;

		/* INTEGRATE & STORE MESH MOVEMENT RESIDUALS */                    
		basis::tri(x.log2p).intgrtx1d(&lf(0,0),&res(0,0));
		basis::tri(x.log2p).intgrt1d(&lf(1,0),&res(1,0));
		basis::tri(x.log2p).intgrtx1d(&lf(1,0),&res(2,0));
#ifdef DROP
		basis::tri(x.log2p).intgrt1d(&lf(2,0),&res(3,0));
#endif
		
		/* STORE IN RES */
		for(n=0;n<tri_mesh::ND;++n) {
			gbl->vres(indx)(n) += lf(n,0);
			gbl->vres(indx+1)(n) = lf(n,1);
			for(m=0;m<basis::tri(x.log2p).sm;++m)
				gbl->sres(indx,m)(n) = lf(n,m+2);
		}

#ifdef DROP
		gbl->vres(indx)(1) += lf(2,0);
		gbl->vres(indx+1)(1) += lf(2,1);
		for(m=0;m<basis::tri(x.log2p).sm;++m)
			gbl->sres(indx,m)(1) += lf(2,m+2);        

		gbl->vvolumeflux(indx) += lf(2,0);
		gbl->vvolumeflux(indx+1) = lf(2,1);
		for(m=0;m<basis::tri(x.log2p).sm;++m)
			gbl->svolumeflux(indx,m) = lf(2,m+2);                    
#endif
		
		/* INTEGRATE & STORE SURFACE TENSION SOURCE TERM */
		basis::tri(x.log2p).intgrt1d(&lf1(0,0),&res(4,0));
		basis::tri(x.log2p).intgrtx1d(&lf1(0,0),&res(5,0));
		basis::tri(x.log2p).intgrt1d(&lf1(1,0),&res(6,0));
		basis::tri(x.log2p).intgrtx1d(&lf1(1,0),&res(7,0));
		
		/* MASS FLUX PRECONDITIONER */				
		for(n=2;n<x.NV-1;++n) /* For swirling case */
			for(m=0;m<basis::tri(x.log2p).sm+2;++m)
				lf1(n,m) = 0.0;

		for(m=0;m<basis::tri(x.log2p).sm+2;++m)
			lf1(x.NV-1,m) = -x.gbl->rho*lf(1,m); 

#ifndef INERTIALESS
		for (n=0;n<x.NV-1;++n) 
			ubar(n) = 0.5*(x.uht(n)(0) +x.uht(n)(1));

		for (n=0;n<x.NV-1;++n) {
			lf1(n,0) -= x.uht(n)(0)*(x.gbl->rho -gbl->rho2)*lf(1,0);
			lf1(n,1) -= x.uht(n)(1)*(x.gbl->rho -gbl->rho2)*lf(1,1);
			for(m=0;m<basis::tri(x.log2p).sm;++m)
				lf1(n,m+2) -= ubar(n)*(x.gbl->rho -gbl->rho2)*lf(1,m+2);
		}
#endif
		/* ADD TO RESIDUAL */
		for(n=0;n<x.NV;++n)
			x.gbl->res.v(v0,n) += lf1(n,0);

		for(n=0;n<x.NV;++n)
			x.gbl->res.v(v1,n) += lf1(n,1);
		
		for(m=0;m<basis::tri(x.log2p).sm;++m) {
			for(n=0;n<x.NV;++n)
				x.gbl->res.s(sind,m,n) += lf1(n,m+2);
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
				for(n=0;n<tri_mesh::ND;++n)
						vdres(x.log2p,i)(n) = gbl->fadd(n)*gbl->vres0(i)(n) -gbl->vres(i)(n);

			for(i=0;i<base.nseg;++i) 
				for(m=0;m<basis::tri(x.log2p).sm;++m) 
					for(n=0;n<tri_mesh::ND;++n)
						sdres(x.log2p,i,m)(n) = gbl->fadd(n)*gbl->sres0(i,m)(n) -gbl->sres(i,m)(n);

		}
		for(i=0;i<base.nseg+1;++i) 
			for(n=0;n<tri_mesh::ND;++n)
				gbl->vres(i)(n) += vdres(x.log2p,i)(n);
		
		for(i=0;i<base.nseg;++i) 
			for(m=0;m<basis::tri(x.log2p).sm;++m) 
				for(n=0;n<tri_mesh::ND;++n)
						gbl->sres(i,m)(n) += sdres(x.log2p,i,m)(n);
	}
	else {
		/* ADD TANGENTIAL MESH MOVEMENT SOURCE */    
		for(i=0;i<base.nseg+1;++i)
			gbl->vres(i)(0) += vdres(x.log2p,i)(0);

		for(i=0;i<base.nseg;++i)
			for(m=0;m<basis::tri(x.log2p).sm;++m) 
				gbl->sres(i,m)(0) += sdres(x.log2p,i,m)(0);
	}

	if (base.is_comm()) { 
		count = 0;
		for(j=0;j<base.nseg+1;++j) {
			base.fsndbuf(count++) = gbl->vres(j)(1)*gbl->rho2;
#ifdef MPDEBUG 
			*x.gbl->log << gbl->vres(j)(1)*gbl->rho2 << '\n';
#endif
#ifdef DROP
			base.fsndbuf(count-1) -= gbl->vvolumeflux(j)*gbl->rho2;
#endif
		}
		for(j=0;j<base.nseg;++j) {
			for(m=0;m<basis::tri(x.log2p).sm;++m) {
				base.fsndbuf(count++) = gbl->sres(j,m)(1)*gbl->rho2;
#ifdef DROP
				base.fsndbuf(count-1) -= gbl->svolumeflux(j,m)*gbl->rho2;
#endif
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

void surface_slave::rsdl(int stage) {
	int i,m,msgn,count,sind,v0;
	
	base.comm_prepare(boundary::all,0,boundary::master_slave);
	base.comm_exchange(boundary::all,0,boundary::master_slave);
	base.comm_wait(boundary::all,0,boundary::master_slave);          
	count = 0;
	for(i=base.nseg-1;i>=0;--i) {
		sind = base.seg(i);
		v0 = x.seg(sind).pnt(1);
#ifdef MPDEBUG
		*x.gbl->log << x.gbl->res.v(v0,x.NV-1) << ' ' << base.frcvbuf(0,count) << '\n';
#endif
		x.gbl->res.v(v0,x.NV-1) += base.frcvbuf(0,count++);
	}
	v0 = x.seg(sind).pnt(0);
#ifdef MPDEBUG
	*x.gbl->log << x.gbl->res.v(v0,x.NV-1) << ' ' << base.frcvbuf(0,count) << '\n';
#endif     
	x.gbl->res.v(v0,x.NV-1) += base.frcvbuf(0,count++);

	for(i=base.nseg-1;i>=0;--i) {
		sind = base.seg(i);
		msgn = 1;
		for(m=0;m<basis::tri(x.log2p).sm;++m) {
			x.gbl->res.s(sind,m,x.NV-1) += msgn*base.frcvbuf(0,count++);
			msgn *= -1;
		}
	}

	return;
}

void surface_outflow_endpt::rsdl(int stage) {
	int bnum,sind;
	TinyVector<FLT,tri_mesh::ND> ubar, tangent, rp;
	FLT jcb;
	
	/* ADD SURFACE TENSION BOUNDARY TERMS IF NECESSARY */
	/* THIS SHOULD REALLY BE PRECALCULATED AND STORED */
	bnum = base.ebdry(surfbdry);
	if (surfbdry == 0) {
		sind = x.ebdry(bnum)->seg(x.ebdry(bnum)->nseg-1);
		x.crdtocht1d(sind);
		basis::tri(x.log2p).ptprobe1d(2,&rp(0),&tangent(0),1.0,&x.cht(0,0),MXTM);
		jcb = sqrt(tangent(0)*tangent(0) +tangent(1)*tangent(1));
		x.gbl->res.v(base.pnt,0) += -RAD(rp(0))*surf->gbl->sigma*tangent(0)/jcb;
		x.gbl->res.v(base.pnt,1) += -RAD(rp(0))*surf->gbl->sigma*tangent(1)/jcb;
	}
	else {
		sind = x.ebdry(bnum)->seg(0);
		x.crdtocht1d(sind);
		basis::tri(x.log2p).ptprobe1d(2,&rp(0),&tangent(0),-1.0,&x.cht(0,0),MXTM);
		jcb = sqrt(tangent(0)*tangent(0) +tangent(1)*tangent(1));
		x.gbl->res.v(base.pnt,0) -= -RAD(rp(0))*surf->gbl->sigma*tangent(0)/jcb;
		x.gbl->res.v(base.pnt,1) -= -RAD(rp(0))*surf->gbl->sigma*tangent(1)/jcb;
	}
	return;
}

void surface::minvrt() {
	int i,m,n,indx;
	int last_phase, mp_phase;
	FLT temp;

	/* INVERT MASS MATRIX */
	/* LOOP THROUGH SIDES */
	if (basis::tri(x.log2p).sm > 0) {
		for(indx = 0; indx<base.nseg; ++indx) {
			/* SUBTRACT SIDE CONTRIBUTIONS TO VERTICES */            
			for (m=0; m <basis::tri(x.log2p).sm; ++m) {
				for(n=0;n<tri_mesh::ND;++n)
						gbl->vres(indx)(n) -= basis::tri(x.log2p).sfmv1d(0,m)*gbl->sres(indx,m)(n);
				for(n=0;n<tri_mesh::ND;++n)
						gbl->vres(indx+1)(n) -= basis::tri(x.log2p).sfmv1d(1,m)*gbl->sres(indx,m)(n);
			}
		}
	}
	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		x.vbdry(base.vbdry(0))->ploadbuff(boundary::manifolds,&gbl->vres(0)(0),0,1,0);
		x.vbdry(base.vbdry(1))->ploadbuff(boundary::manifolds,&gbl->vres(base.nseg)(0),0,1,0);
		x.vbdry(base.vbdry(0))->comm_prepare(boundary::manifolds,mp_phase,boundary::symmetric);
		x.vbdry(base.vbdry(1))->comm_prepare(boundary::manifolds,mp_phase,boundary::symmetric);

		x.vbdry(base.vbdry(0))->comm_exchange(boundary::manifolds,mp_phase,boundary::symmetric);
		x.vbdry(base.vbdry(1))->comm_exchange(boundary::manifolds,mp_phase,boundary::symmetric);
		
		last_phase = true;

		last_phase &= x.vbdry(base.vbdry(0))->comm_wait(boundary::manifolds,mp_phase,boundary::symmetric);
		last_phase &= x.vbdry(base.vbdry(1))->comm_wait(boundary::manifolds,mp_phase,boundary::symmetric);
		
		x.vbdry(base.vbdry(0))->pfinalrcv(boundary::manifolds,mp_phase,boundary::symmetric,boundary::average,&gbl->vres(0)(0),0,1,0);
		x.vbdry(base.vbdry(1))->pfinalrcv(boundary::manifolds,mp_phase,boundary::symmetric,boundary::average,&gbl->vres(base.nseg)(0),0,1,0);
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
	if (basis::tri(x.log2p).sm > 0) {
		for(indx = 0; indx<base.nseg; ++indx) {                   

#ifdef DETAILED_MINV
			for(m=0;m<basis::tri(x.log2p).sm;++m) {
				for(n=0;n<tri_mesh::ND;++n) {
						gbl->sres(indx,m)(n) -= gbl->vms(indx,n,0,m,0)*gbl->vres(indx)(0);
						gbl->sres(indx,m)(n) -= gbl->vms(indx,n,0,m,1)*gbl->vres(indx+1)(0);
						gbl->sres(indx,m)(n) -= gbl->vms(indx,n,1,m,0)*gbl->vres(indx)(1);
						gbl->sres(indx,m)(n) -= gbl->vms(indx,n,1,m,1)*gbl->vres(indx+1)(1);                            
				}
			}
			int info;
			char trans[] = "T";
			GETRS(trans,2*basis::tri(x.log2p).sm,1,&gbl->ms(indx)(0,0),2*MAXP,&gbl->ipiv(indx)(0),&gbl->sres(indx,0)(0),2*MAXP,info);
			if (info != 0) {
				printf("DGETRS FAILED FOR SIDE MODE UPDATE\n");
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
#else
	
			/* INVERT SIDE MODES */
			DPBTRSNU2(&basis::tri(x.log2p).sdiag1d(0,0),basis::tri(x.log2p).sbwth+1,basis::tri(x.log2p).sm,basis::tri(x.log2p).sbwth,&(gbl->sres(indx,0)(0)),tri_mesh::ND);
			for(m=0;m<basis::tri(x.log2p).sm;++m) {
				temp                             = gbl->sres(indx,m)(0)*gbl->sdt(indx)(0,0) +gbl->sres(indx,m)(1)*gbl->sdt(indx)(0,1);
				gbl->sres(indx,m)(1) = gbl->sres(indx,m)(0)*gbl->sdt(indx)(1,0) +gbl->sres(indx,m)(1)*gbl->sdt(indx)(1,1);         
				gbl->sres(indx,m)(0) = temp;
			}

			for(m=0;m<basis::tri(x.log2p).sm;++m) {
				for(n=0;n<tri_mesh::ND;++n) {
						gbl->sres(indx,m)(n) -= basis::tri(x.log2p).vfms1d(0,m)*gbl->vres(indx)(n);
						gbl->sres(indx,m)(n) -= basis::tri(x.log2p).vfms1d(1,m)*gbl->vres(indx+1)(n);
				}
			}
#endif
		}
	}

	return;
}

void surface::setup_preconditioner() {
	int indx,m,n,sind,v0,v1;
	TinyVector<FLT,tri_mesh::ND> nrm;
	FLT h, hsm;
	FLT dttang, dtnorm;
	FLT vslp, strss;
	FLT drho, srho, smu;
	FLT nu1, nu2;
	FLT qmax, gam1, gam2;
	TinyMatrix<FLT,tri_mesh::ND,MXGP> crd, dcrd;
	TinyMatrix<FLT,4,MXGP> res;
	TinyMatrix<FLT,4,MXGP> lf;
	TinyVector<FLT,2> mvel;
	Array<TinyVector<FLT,MXGP>,1> u(x.NV);
	int last_phase, mp_phase;

	drho = x.gbl->rho -gbl->rho2;
	srho = x.gbl->rho +gbl->rho2;
	smu = x.gbl->mu +gbl->mu2;
	nu1 = x.gbl->mu/x.gbl->rho;
	if (gbl->rho2 > 0.0) 
		nu2 = gbl->mu2/gbl->rho2;
	else
		nu2 = 0.0;

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
			basis::tri(x.log2p).proj1d(&x.cht(n,0),&crd(n,0),&dcrd(n,0));

		x.ugtouht1d(sind);
		for(n=0;n<tri_mesh::ND;++n)
			basis::tri(x.log2p).proj1d(&x.uht(n)(0),&u(n)(0));    
		
		dtnorm = 1.0e99;
		dttang = 1.0e99;
		gbl->meshc(indx) = 1.0e99;
		for(i=0;i<basis::tri(x.log2p).gpx;++i) {
			nrm(0) =  dcrd(1,i)*2;
			nrm(1) = -dcrd(0,i)*2;
			h = sqrt(nrm(0)*nrm(0) +nrm(1)*nrm(1));

			/* RELATIVE VELOCITY STORED IN MVEL(N)*/
			for(n=0;n<tri_mesh::ND;++n) {
				mvel(n) = u(n)(i) -(x.gbl->bd(0)*(crd(n,i) -dxdt(x.log2p,indx)(n,i)));
#ifdef DROP
				mvel(n) -= tri_hp_ins::mesh_ref_vel(n);
#endif    
			}



			qmax = u(0)(i)*u(0)(i) +u(1)(i)*u(1)(i);
			vslp = fabs(-u(0)(i)*nrm(1)/h +u(1)(i)*nrm(0)/h);
			hsm = h/(.25*(basis::tri(x.log2p).p+1)*(basis::tri(x.log2p).p+1));

			dttang = MIN(dttang,2.*ksprg(indx)*(.25*(basis::tri(x.log2p).p+1)*(basis::tri(x.log2p).p+1))/hsm);
#ifndef BODYFORCE
			strss =  4.*gbl->sigma/(hsm*hsm) +fabs(drho*sim::g*nrm(1)/h);
#else
			strss =  4.*gbl->sigma/(hsm*hsm) +fabs(drho*(-sim::body(0)*nrm(0) +(sim::g-sim::body(1))*nrm(1))/h);
#endif

			gam1 = 3.0*qmax +(0.5*hsm*x.gbl->bd(0) + 2.*nu1/hsm)*(0.5*hsm*x.gbl->bd(0) + 2.*nu1/hsm);
			gam2 = 3.0*qmax +(0.5*hsm*x.gbl->bd(0) + 2.*nu2/hsm)*(0.5*hsm*x.gbl->bd(0) + 2.*nu2/hsm);

			if (x.gbl->bd(0) + x.gbl->mu == 0.0) gam1 = MAX(gam1,0.1);

#ifdef INERTIALESS
			gam1 = (2.*nu1/hsm)*(2.*nu1/hsm);
			gam2 = (2.*nu2/hsm)*(2.*nu2/hsm);
#endif
			dtnorm = MIN(dtnorm,2.*vslp/hsm +x.gbl->bd(0) +1.*strss/(x.gbl->rho*sqrt(qmax +gam1) +gbl->rho2*sqrt(qmax +gam2)));                
	
			/* SET UP DISSIPATIVE COEFFICIENT */
			/* FOR UPWINDING LINEAR CONVECTIVE CASE SHOULD BE 1/|a| */
			/* RESIDUAL HAS DX/2 WEIGHTING */
			/* |a| dx/2 dv/dx  dx/2 dpsi */
			/* |a| dx/2 2/dx dv/dpsi  dpsi */
			/* |a| dv/dpsi  dpsi */
			// gbl->meshc(indx) = gbl->adis/(h*dtnorm*0.5);/* FAILED IN NATES UPSTREAM SURFACE WAVE CASE */
			// gbl->meshc(indx) = MIN(gbl->meshc(indx),gbl->adis/(h*(vslp/hsm +x.gbl->bd(0)))); /* FAILED IN MOVING UP TESTS */
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
#ifdef DROP
		mvel(0) -= tri_hp_ins::mesh_ref_vel(0);
		mvel(1) -= tri_hp_ins::mesh_ref_vel(1);
#endif

		qmax = mvel(0)*mvel(0)+mvel(1)*mvel(1);
		vslp = fabs(-mvel(0)*nrm(1)/h +mvel(1)*nrm(0)/h);

		mvel(0) = x.ug.v(v1,0)-(x.gbl->bd(0)*(x.pnts(v1)(0) -x.vrtxbd(1)(v1)(0)));
		mvel(1) = x.ug.v(v1,1)-(x.gbl->bd(0)*(x.pnts(v1)(1) -x.vrtxbd(1)(v1)(1)));
#ifdef DROP
		mvel(0) -= tri_hp_ins::mesh_ref_vel(0);
		mvel(1) -= tri_hp_ins::mesh_ref_vel(1);
#endif
		qmax = MAX(qmax,mvel(0)*mvel(0)+mvel(1)*mvel(1));
		vslp = MAX(vslp,fabs(-mvel(0)*nrm(1)/h +mvel(1)*nrm(0)/h));
		
		hsm = h/(.25*(basis::tri(x.log2p).p+1)*(basis::tri(x.log2p).p+1));

		dttang = 2.*ksprg(indx)*(.25*(basis::tri(x.log2p).p+1)*(basis::tri(x.log2p).p+1))/hsm;
#ifndef BODYFORCE
		strss =  4.*gbl->sigma/(hsm*hsm) +fabs(drho*x.gbl->g*nrm(1)/h);
#else
		strss =  4.*gbl->sigma/(hsm*hsm) +fabs(drho*(-sim::body(0)*nrm(0) +(sim::g-sim::body(1))*nrm(1))/h);
#endif

		gam1 = 3.0*qmax +(0.5*hsm*x.gbl->bd(0) + 2.*nu1/hsm)*(0.5*hsm*x.gbl->bd(0) + 2.*nu1/hsm);
		gam2 = 3.0*qmax +(0.5*hsm*x.gbl->bd(0) + 2.*nu2/hsm)*(0.5*hsm*x.gbl->bd(0) + 2.*nu2/hsm);

		if (x.gbl->bd(0) + x.gbl->mu == 0.0) gam1 = MAX(gam1,0.1);

#ifdef INERTIALESS
		gam1 = (2.*nu1/hsm)*(2.*nu1/hsm);
		gam2 = (2.*nu2/hsm)*(2.*nu2/hsm);
#endif
		dtnorm = 2.*vslp/hsm +x.gbl->bd(0) +1.*strss/(x.gbl->rho*sqrt(qmax +gam1) +gbl->rho2*sqrt(qmax +gam2));                
		
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

		gbl->vdt(indx)(0,0) += -dttang*nrm(1)*basis::tri(x.log2p).vdiag1d;
		gbl->vdt(indx)(0,1) +=  dttang*nrm(0)*basis::tri(x.log2p).vdiag1d;
		gbl->vdt(indx)(1,0) +=  dtnorm*nrm(0)*basis::tri(x.log2p).vdiag1d;
		gbl->vdt(indx)(1,1) +=  dtnorm*nrm(1)*basis::tri(x.log2p).vdiag1d;
		gbl->vdt(indx+1)(0,0) = -dttang*nrm(1)*basis::tri(x.log2p).vdiag1d;
		gbl->vdt(indx+1)(0,1) =  dttang*nrm(0)*basis::tri(x.log2p).vdiag1d;
		gbl->vdt(indx+1)(1,0) =  dtnorm*nrm(0)*basis::tri(x.log2p).vdiag1d;
		gbl->vdt(indx+1)(1,1) =  dtnorm*nrm(1)*basis::tri(x.log2p).vdiag1d;

		if (basis::tri(x.log2p).sm) {
			gbl->sdt(indx)(0,0) = -dttang*nrm(1);
			gbl->sdt(indx)(0,1) =  dttang*nrm(0);
			gbl->sdt(indx)(1,0) =  dtnorm*nrm(0);
			gbl->sdt(indx)(1,1) =  dtnorm*nrm(1);  
		
#ifdef DETAILED_MINV
			int lsm = basis::tri(x.log2p).sm;
			x.crdtocht1d(sind);
			for(n=0;n<tri_mesh::ND;++n)
				basis::tri(x.log2p).proj1d(&x.cht(n,0),&crd(n,0),&dcrd(n,0));

			for(int m = 0; m<lsm; ++m) {
				for(i=0;i<basis::tri(x.log2p).gpx;++i) {
						nrm(0) =  dcrd(1,i);
						nrm(1) = -dcrd(0,i);
						res(0,i) = -dttang*nrm(1)*basis::tri(x.log2p).gx(i,m+3);
						res(1,i) =  dttang*nrm(0)*basis::tri(x.log2p).gx(i,m+3);
						res(2,i) =  dtnorm*nrm(0)*basis::tri(x.log2p).gx(i,m+3);
						res(3,i) =  dtnorm*nrm(1)*basis::tri(x.log2p).gx(i,m+3); 
				}
				lf = 0;
				basis::tri(x.log2p).intgrt1d(&lf(0,0),&res(0,0));
				basis::tri(x.log2p).intgrt1d(&lf(1,0),&res(1,0));
				basis::tri(x.log2p).intgrt1d(&lf(2,0),&res(2,0));
				basis::tri(x.log2p).intgrt1d(&lf(3,0),&res(3,0));

				/* CFL = 0 WON'T WORK THIS WAY */
				lf(0,Range(0,lsm+1) /= gbl->cfl(0,x.log2p);
				lf(1,Range(0,lsm+1) /= gbl->cfl(0,x.log2p);
				lf(2,Range(0,lsm+1) /= gbl->cfl(1,x.log2p);
				lf(3,Range(0,lsm+1) /= gbl->cfl(1,x.log2p);                        

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
				printf("DGETRF FAILED IN SIDE MODE PRECONDITIONER\n");
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			/*
			\phi_n dx,dy*t = \phi_n Vt
			\phi_t dx,dy*n = \phi_t Vn
			*/
#endif
		}
	}

	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		x.vbdry(base.vbdry(0))->ploadbuff(boundary::manifolds,&gbl->vdt(0)(0,0),0,3,0);
		x.vbdry(base.vbdry(1))->ploadbuff(boundary::manifolds,&gbl->vdt(base.nseg)(0,0),0,3,0);
		x.vbdry(base.vbdry(0))->comm_prepare(boundary::manifolds,mp_phase,boundary::symmetric);
		x.vbdry(base.vbdry(1))->comm_prepare(boundary::manifolds,mp_phase,boundary::symmetric);

		x.vbdry(base.vbdry(0))->comm_exchange(boundary::manifolds,mp_phase,boundary::symmetric);
		x.vbdry(base.vbdry(1))->comm_exchange(boundary::manifolds,mp_phase,boundary::symmetric);        

		last_phase = true;
		last_phase &= x.vbdry(base.vbdry(0))->comm_wait(boundary::manifolds,mp_phase,boundary::symmetric);
		last_phase &= x.vbdry(base.vbdry(1))->comm_wait(boundary::manifolds,mp_phase,boundary::symmetric);
		x.vbdry(base.vbdry(0))->pfinalrcv(boundary::manifolds,mp_phase,boundary::symmetric,boundary::average,&gbl->vdt(0)(0,0),0,3,0);
		x.vbdry(base.vbdry(1))->pfinalrcv(boundary::manifolds,mp_phase,boundary::symmetric,boundary::average,&gbl->vdt(base.nseg)(0,0),0,3,0);
	}
	
	if (gbl->is_loop) {
		for(m=0;m<tri_mesh::ND;++m)
			for(n=0;n<tri_mesh::ND;++n)
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
	if (basis::tri(x.log2p).sm > 0) {
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

void surface::update(int stage) {
	int i,m,n,count,sind,indx,v0;

	if (stage < 0) {
		indx = 0;
		for(i=0;i<base.nseg;++i) {
			sind = base.seg(i);
			v0 = x.seg(sind).pnt(0);
			gbl->vug0(i) = x.pnts(v0);
		}
		v0 = x.seg(sind).pnt(1);
		gbl->vug0(base.nseg) = x.pnts(v0);

		if (basis::tri(x.log2p).sm > 0) gbl->sug0(Range(0,base.nseg-1),Range(0,basis::tri(x.log2p).sm-1)) = crv(Range(0,base.nseg-1),Range(0,basis::tri(x.log2p).sm-1));
		return;
	}

	minvrt();

#ifdef DEBUG
//   if (x.coarse_flag) {
	for(i=0;i<base.nseg+1;++i)
		printf("vdt: %d %8.4e %8.4e %8.4e %8.4e\n",i,gbl->vdt(i)(0,0),gbl->vdt(i)(0,1),gbl->vdt(i)(1,0),gbl->vdt(i)(1,1));

	for(i=0;i<base.nseg;++i)
		printf("sdt: %d %8.4e %8.4e %8.4e %8.4e\n",i,gbl->sdt(i)(0,0),gbl->sdt(i)(0,1),gbl->sdt(i)(1,0),gbl->sdt(i)(1,1));

	for(i=0;i<base.nseg+1;++i) {
		printf("vres: %d ",i);
		for(n=0;n<tri_mesh::ND;++n) {
			if (fabs(gbl->vres(i)(n)) > 1.0e-9) printf("%8.4e ",gbl->vres(i)(n));
			else printf("%8.4e ",0.0);
		}
		printf("\n");
	}
		
	for(i=0;i<base.nseg;++i) {
		for(m=0;m<basis::tri(x.log2p).sm;++m) {
			printf("sres: %d ",i);
			for(n=0;n<tri_mesh::ND;++n) {
				if (fabs(gbl->sres(i,m)(n)) > 1.0e-9) printf("%8.4e ",gbl->sres(i,m)(n));
				else printf("%8.4e ",0.0);
			}
			printf("\n");
		}
	}

	for(i=0;i<base.nseg;++i) {
		sind = base.seg(i);
		v0 = x.seg(sind).pnt(0);
		printf("vertex positions %d %8.4e %8.4e\n",v0,x.pnts(v0)(0),x.pnts(v0)(1));
	}
	v0 = x.seg(sind).pnt(1);
	printf("vertex positions %d %8.4e %8.4e\n",v0,x.pnts(v0)(0),x.pnts(v0)(1));
	
	for(i=0;i<base.nseg;++i)
		for(m=0;m<basis::tri(x.log2p).sm;++m)
			printf("spos: %d %d %8.4e %8.4e\n",i,m,crv(i,m)(0),crv(i,m)(1));
	// }
#endif

	for(i=0;i<base.nseg;++i) {
		sind = base.seg(i);
		v0 = x.seg(sind).pnt(0);
		x.pnts(v0) = gbl->vug0(i) -sim::alpha[stage]*gbl->vres(i);
	}
	v0 = x.seg(sind).pnt(1);
	x.pnts(v0) = gbl->vug0(base.nseg) -sim::alpha[stage]*gbl->vres(base.nseg);
	
	if (basis::tri(x.log2p).sm > 0) crv(Range(0,base.nseg-1),Range(0,basis::tri(x.log2p).sm-1)) = gbl->sug0(Range(0,base.nseg-1),Range(0,basis::tri(x.log2p).sm-1)) -sim::alpha[stage]*gbl->sres(Range(0,base.nseg-1),Range(0,basis::tri(x.log2p).sm-1));

	/* FIX POINTS THAT SLIDE ON CURVE */
	x.hp_vbdry(base.vbdry(1))->mvpttobdry(x.pnts(v0));
	sind = base.seg(0);
	v0 = x.seg(sind).pnt(0);
	x.hp_vbdry(base.vbdry(0))->mvpttobdry(x.pnts(v0));
	
#ifdef DEBUG
//   if (x.coarse_flag) {
	for(i=0;i<base.nseg;++i) {
		sind = base.seg(i);
		v0 = x.seg(sind).pnt(0);
		printf("vertex positions %d %e %e\n",v0,x.pnts(v0)(0),x.pnts(v0)(1));
	}
	v0 = x.seg(sind).pnt(1);
	printf("vertex positions %d %e %e\n",v0,x.pnts(v0)(0),x.pnts(v0)(1));
	
	for(i=0;i<base.nseg;++i)
		for(m=0;m<basis::tri(x.log2p).sm;++m)
			printf("spos: %d %d %e %e\n",i,m,crv(i,m)(0),crv(i,m)(1));
//   }
#endif
	
	if (base.is_comm()) {                
		count = 0;
		for(i=0;i<base.nseg;++i) {
			sind = base.seg(i);
			v0 = x.seg(sind).pnt(0);
			for(n=0;n<tri_mesh::ND;++n)
				base.fsndbuf(count++) = x.pnts(v0)(n);
		}
		v0 = x.seg(sind).pnt(1);
		for(n=0;n<tri_mesh::ND;++n)
			base.fsndbuf(count++) = x.pnts(v0)(n);                
		
		for(i=0;i<base.nseg;++i) {
			for(m=0;m<basis::tri(x.log2p).sm;++m) {
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

void surface_slave::update(int stage) {
	int i,m,n,msgn,count,sind,v0;

	if (stage < 0) return;
	
	base.comm_prepare(boundary::all,0,boundary::master_slave);
	base.comm_exchange(boundary::all,0,boundary::master_slave);
	base.comm_wait(boundary::all,0,boundary::master_slave);
	
	count = 0;
	for(i=base.nseg-1;i>=0;--i) {
		sind = base.seg(i);
		v0 = x.seg(sind).pnt(1);
		for(n=0;n<tri_mesh::ND;++n) 
			x.pnts(v0)(n) = base.frcvbuf(0,count++);
	}
	v0 = x.seg(sind).pnt(0);
	for(n=0;n<tri_mesh::ND;++n)
		x.pnts(v0)(n) = base.frcvbuf(0,count++);
	
	if (basis::tri(x.log2p).sm > 0) {
		for(i=base.nseg-1;i>=0;--i) {
			sind = base.seg(i);
			msgn = 1;
			for(m=0;m<basis::tri(x.log2p).sm;++m) {
				for(n=0;n<tri_mesh::ND;++n)
						crds(i,m,n) = msgn*base.frcvbuf(0,count++);
				msgn *= -1;
			}
		}
	}
	return;
}

void surface::mg_restrict() {
	int i,bnum,indx,tind,v0,snum,sind;

		if(x.p0 > 1) {
		/* TRANSFER IS ON FINEST MESH */
		gbl->vres0(Range(0,base.nseg)) = gbl->vres(Range(0,base.nseg));
		if (basis::tri(x.log2p).sm > 0) gbl->sres0(Range(0,base.nseg-1),Range(0,basis::tri(x.log2p).sm-1)) = gbl->sres(Range(0,base.nseg-1),Range(0,basis::tri(x.log2p).sm-1));
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

void surface_slave::smatchsolution_snd(FLT *sdata, int bgn, int end, int stride) {
	int j,k,n,countup,offset;
	
	if (!base.is_comm()) return;
	
#ifdef MPDEBUG
		*x.gbl->log << base.idprefix << " In surface_snd"  << base.idnum << " " << base.is_frst() << std::endl;
#endif
	
	countup = 0;
	for(j=0;j<base.nseg;++j) {
		offset = base.seg(j)*stride*x.NV;
		for(k=bgn;k<=end;++k) {
			for(n=0;n<x.NV-1;++n) {
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


void surface_slave::smatchsolution_rcv(FLT *sdata, int bgn, int end, int stride) {
	
	if (!base.is_comm()) return;
		
	/* CAN'T USE sfinalrcv BECAUSE OF CHANGING SIGNS */
	int j,k,m,n,count,countdn,countup,offset,sind,sign;
	FLT mtchinv;
	
	/* ASSUMES REVERSE ORDERING OF SIDES */
	/* WON'T WORK IN 3D */
	
	int matches = 1;
	
	int bgnsign = (bgn % 2 ? -1 : 1);
	
	/* RELOAD FROM BUFFER */
	/* ELIMINATES V/S/F COUPLING IN ONE PHASE */
	/* FINALRCV SHOULD BE CALLED F,S,V ORDER (V HAS FINAL AUTHORITY) */    
	for(m=0;m<base.matches();++m) {            
		++matches;
		
		int ebp1 = end-bgn+1;

		countdn = (base.nseg-1)*ebp1*(x.NV-1);
		countup = 0;
		for(j=0;j<base.nseg;++j) {
			sign = bgnsign;
			for(k=0;k<ebp1;++k) {
				for(n=0;n<x.NV-1;++n) {
						base.fsndbuf(countup++) += sign*base.frcvbuf(m,countdn++);
				}
				sign *= -1;
			}
			countdn -= 2*ebp1*(x.NV-1);
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
				for(n=0;n<x.NV-1;++n) {
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
	

void surface::maxres() {
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
