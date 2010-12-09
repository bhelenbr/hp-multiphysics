#include "bdry_ins.h"
#include <myblas.h>

/*************************************************/
/* SET DIRICHLET BOUNDARY VALUES & FLUXES ********/
/* (THINGS THAT ARE INDEPENDENT OF THE SOLUTION) */
/* BUT NOT INDEPENDENT OF TIME/MESH POSITION     */
/*************************************************/
using namespace bdry_ins;

// #define MPDEBUG

void generic::output(std::ostream& fout, tri_hp::filetype typ,int tlvl) {
	int i,m,n,ind,sind,tind,seg;
	FLT visc[tri_mesh::ND+1][tri_mesh::ND+1][tri_mesh::ND][tri_mesh::ND];
	TinyVector<FLT,tri_mesh::ND> norm, mvel;
	FLT convect,jcb;

	switch(typ) {
		case(tri_hp::text): case(tri_hp::binary): {
			hp_edge_bdry::output(fout,typ,tlvl);
			break;
		}
		case(tri_hp::tecplot): {
			if (!report_flag) return;

			conv_flux = 0.0;
			diff_flux = 0.0;
			moment = 0.0;
			circumference = 0.0;
			circulation = 0.0;
#ifdef L2_ERROR
			FLT l2error = 0.0;
			TinyVector<FLT,2> xpt;
#endif
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

					x.cjcb(0,i) = x.gbl->mu*RAD(x.crd(0)(0,i))/(x.dcrd(0,0)(0,i)*x.dcrd(1,1)(0,i) -x.dcrd(1,0)(0,i)*x.dcrd(0,1)(0,i));

					/* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
					/* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
					visc[0][0][0][0] =  x.cjcb(0,i)*(2.*x.dcrd(1,1)(0,i)*x.dcrd(1,1)(0,i) +x.dcrd(0,1)(0,i)*x.dcrd(0,1)(0,i));
					visc[0][0][1][1] =  x.cjcb(0,i)*(2.*x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i));
					visc[0][0][0][1] = -x.cjcb(0,i)*(2.*x.dcrd(1,1)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,1)(0,i)*x.dcrd(0,0)(0,i));
#define             viscI0II0II1II0I visc[0][0][0][1]

					visc[1][1][0][0] =  x.cjcb(0,i)*(x.dcrd(1,1)(0,i)*x.dcrd(1,1)(0,i) +2.*x.dcrd(0,1)(0,i)*x.dcrd(0,1)(0,i));
					visc[1][1][1][1] =  x.cjcb(0,i)*(x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i) +2.*x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i));
					visc[1][1][0][1] = -x.cjcb(0,i)*(x.dcrd(1,1)(0,i)*x.dcrd(1,0)(0,i) +2.*x.dcrd(0,1)(0,i)*x.dcrd(0,0)(0,i));
#define             viscI1II1II1II0I visc[1][1][0][1]

					visc[0][1][0][0] = -x.cjcb(0,i)*x.dcrd(0,1)(0,i)*x.dcrd(1,1)(0,i);
					visc[0][1][1][1] = -x.cjcb(0,i)*x.dcrd(0,0)(0,i)*x.dcrd(1,0)(0,i);
					visc[0][1][0][1] =  x.cjcb(0,i)*x.dcrd(0,1)(0,i)*x.dcrd(1,0)(0,i);
					visc[0][1][1][0] =  x.cjcb(0,i)*x.dcrd(0,0)(0,i)*x.dcrd(1,1)(0,i);

					/* OTHER SYMMETRIES     */                
#define             viscI1II0II0II0I visc[0][1][0][0]
#define             viscI1II0II1II1I visc[0][1][1][1]
#define             viscI1II0II0II1I visc[0][1][1][0]
#define             viscI1II0II1II0I visc[0][1][0][1]

					/* DIFFUSIVE FLUXES ( FOR EXTRA VARIABLES) */
					visc[2][2][1][0] =  x.cjcb(0,i)*(x.dcrd(1,1)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,1)(0,i)*x.dcrd(0,0)(0,i));
					visc[2][2][1][1] = -x.cjcb(0,i)*(x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i));

					for (n=tri_mesh::ND;n<x.NV-1;++n) 
						diff_flux(n) -= x.gbl->D(n)/x.gbl->mu*basis::tri(x.log2p)->wtx(i)*(-visc[2][2][1][0]*x.du(n,0)(0,i) -visc[2][2][1][1]*x.du(n,1)(0,i));

					diff_flux(0) -=    basis::tri(x.log2p)->wtx(i)*(-x.u(2)(0,i)*RAD(x.crd(0)(0,i))*x.dcrd(1,0)(0,i) 
									-viscI0II0II1II0I*x.du(0,0)(0,i) -visc[0][1][1][0]*x.du(1,0)(0,i)
									-visc[0][0][1][1]*x.du(0,1)(0,i) -visc[0][1][1][1]*x.du(1,1)(0,i));															
					diff_flux(1) -=    basis::tri(x.log2p)->wtx(i)*( x.u(2)(0,i)*RAD(x.crd(0)(0,i))*x.dcrd(0,0)(0,i)
									-viscI1II0II1II0I*x.du(0,0)(0,i) -viscI1II1II1II0I*x.du(1,0)(0,i)
									-viscI1II0II1II1I*x.du(0,1)(0,i) -visc[1][1][1][1]*x.du(1,1)(0,i));


					norm(0) = x.dcrd(1,0)(0,i);
					norm(1) = -x.dcrd(0,0)(0,i);                
					for(n=0;n<tri_mesh::ND;++n) {
						mvel(n) = x.gbl->bd(0)*(x.crd(n)(0,i) -dxdt(x.log2p,ind)(n,i));
#ifdef DROP
						mvel(n) += tri_hp_ins::mesh_ref_vel(n);
#endif
					}

					circulation += basis::tri(x.log2p)->wtx(i)*(-norm(1)*(x.u(0)(0,i)-mvel(0)) +norm(0)*(x.u(1)(0,i)-mvel(1)));

					convect = basis::tri(x.log2p)->wtx(i)*RAD(x.crd(0)(0,i))*((x.u(0)(0,i)-mvel(0))*norm(0) +(x.u(1)(0,i)-mvel(1))*norm(1));
					conv_flux(2) -= convect;
					conv_flux(0) -= x.u(0)(0,i)*convect;
					conv_flux(1) -= x.u(1)(0,i)*convect;

#ifdef L2_ERROR
					xpt(0) = x.crd(0)(0,i);
					xpt(1) = x.crd(1)(0,i);
					l2error += jcb*l2norm.Eval(xpt,x.gbl->time);
#endif
				}	
			} while (++ind < base.nseg);
			streamsize oldprecision = (*x.gbl->log).precision(10);
			*x.gbl->log << base.idprefix << " circumference: " << circumference << std::endl;
			*x.gbl->log << base.idprefix << " viscous/pressure flux: " << diff_flux << std::endl;
			*x.gbl->log << base.idprefix << " convective flux: " << conv_flux << std::endl;
			*x.gbl->log << base.idprefix << " circulation: " << circulation << std::endl;

			/* OUTPUT AUXILIARY FLUXES */
			*x.gbl->log << base.idprefix << "total fluxes: " << total_flux << std::endl;
#ifdef L2_ERROR
			*x.gbl->log << base.idprefix << "l2error: " << sqrt(l2error) << std::endl;
#endif
			(*x.gbl->log).precision(oldprecision);
			break;
		}

		case(tri_hp::auxiliary): {
			if (!report_flag) return;                

			/* AUXILIARY FLUX METHOD */
			int v0;
			total_flux = 0.0;
			ind = 0;
			do {
				sind = base.seg(ind);
				v0 = x.seg(sind).pnt(0);
				total_flux += x.gbl->res.v(v0,Range::all());
			} while (++ind < base.nseg);
			v0 = x.seg(sind).pnt(1);
			total_flux += x.gbl->res.v(v0,Range::all());
		}
		default: {
			break;
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
#ifdef DROP
			mvel(n) += tri_hp_ins::mesh_ref_vel(n);
#endif
		}

		for(n=0;n<x.NV;++n)
			u(n) = x.u(n)(0,k);

		flux(u,pt,mvel,nrm,flx);

		for(n=0;n<x.NV;++n)
			x.res(n)(0,k) = RAD(x.crd(0)(0,k))*flx(n);

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

void applied_stress::init(input_map& inmap,void* gbl_in) {
	std::string keyword;
	std::ostringstream nstr;

	neumann::init(inmap,gbl_in);

	stress.resize(tri_mesh::ND);

	for(int n=0;n<tri_mesh::ND;++n) {
		nstr.str("");
		nstr << base.idprefix << "_stress" << n << std::flush;
		if (inmap.find(nstr.str()) != inmap.end()) {
			stress(n).init(inmap,nstr.str());
		}
		else {
			*x.gbl->log << "couldn't find stress function " << nstr.str() << std::endl;
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}
	}

	return;
}

void flexible::init(input_map& inmap,void* gbl_in) {
	std::string keyword;
	std::ostringstream nstr;

	inflow::init(inmap,gbl_in);

	ibc = x.getnewibc(base.idprefix +"_fcn",inmap);

	Array<int,1> atemp(x.NV-1);
	if (!inmap.get(base.idprefix+"_ins_bcs", atemp.data(), x.NV-1)) {
		*x.gbl->log << "missing flexible specifier list (0 = essential, 1 = natural, 2 = mixed)" << std::endl;
		sim::abort(__LINE__,__FILE__,x.gbl->log);
	}
	else {
		ndirichlets = 0;
		for (int n=0;n<x.NV-1;++n) {
			type(n) = static_cast<bctypes>(atemp(n));
			if (type(n) == ess) {
				dirichlets(ndirichlets++) = n;
			}
		}
	}

	return;
}


void symmetry::tadvance() {
	int j,m,v0,sind;
	TinyVector<FLT,tri_mesh::ND> pt;

	hp_edge_bdry::tadvance();

	/* UPDATE BOUNDARY CONDITION VALUES */
	j = 0;
	do {
		sind = base.seg(j);
		v0 = x.seg(sind).pnt(0);
		x.ug.v(v0,dir) = 0.0;
	}	while (++j < base.nseg);
	v0 = x.seg(sind).pnt(1);
	x.ug.v(v0,dir) = 0.0;

	/*******************/    
	/* SET SIDE VALUES */
	/*******************/
	for(j=0;j<base.nseg;++j) {
		sind = base.seg(j);
		for(m=0;m<basis::tri(x.log2p)->sm();++m) {
			x.ug.s(sind,m,dir) = 0.0;
		}
	}

	return;
}


void characteristic::flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, Array<FLT,1>& flx) {
	FLT ul,vl,ur,vr,pl,pr,cl,cr,rho,rhoi;
	FLT s,um,v,c,den,lam0,lam1,lam2,mag,hmax;
	FLT nu,gam,qmax;
	Array<FLT,1> ub(x.NV), uvp(x.NV);

	/* CHARACTERISTIC FAR-FIELD B.C. */   

	rho = x.gbl->rho;
	nu = x.gbl->mu/x.gbl->rho;
	rhoi = 1./rho;
	mag = sqrt(norm(0)*norm(0) + norm(1)*norm(1));
	hmax = mag*2.0/(0.25*(basis::tri(x.log2p)->p() +1)*(basis::tri(x.log2p)->p()+1));
	qmax = pow(u(0)-0.5*mv(0),2.0) +pow(u(1)-0.5*mv(1),2.0);
	gam = 3.0*qmax +(0.5*hmax*x.gbl->bd(0) +2.*nu/hmax)*(0.5*hmax*x.gbl->bd(0) +2.*nu/hmax);

	norm(0) /= mag;
	norm(1) /= mag;

	ul =  u(0)*norm(0) +u(1)*norm(1);
	vl = -u(0)*norm(1) +u(1)*norm(0);
	pl =  u(x.NV-1);

	/* FREESTREAM CONDITIONS */
	for(int n=0;n<x.NV;++n)
		ub(n) = ibc->f(n,xpt,x.gbl->time);

	ur =  ub(0)*norm(0) +ub(1)*norm(1);
	vr = -ub(0)*norm(1) +ub(1)*norm(0);
	pr =  ub(x.NV-1);

	um = mv(0)*norm(0) +mv(1)*norm(1);

	cl = sqrt((ul-.5*um)*(ul-.5*um) +gam);
	cr = sqrt((ur-.5*um)*(ur-.5*um) +gam);
	c = 0.5*(cl+cr);
	s = 0.5*(ul+ur);
	v = 0.5*(vl+vr);

	den = 1./(2*c);
	lam0 = s -um;
	lam1 = s-.5*um +c; /* always positive */
	lam2 = s-.5*um -c; /* always negative */

	/* PERFORM CHARACTERISTIC SWAP */
	/* BASED ON LINEARIZATION AROUND UL,VL,PL */
	uvp(0) = ((pl-pr)*rhoi +(ul*lam1 -ur*lam2))*den;
	if (lam0 > 0.0) {
		uvp(1) = v*((pr-pl)*rhoi +lam2*(ur-ul))*den/(lam0-lam2) +vl;
		for(int n=tri_mesh::ND;n<x.NV-1;++n)
			uvp(n) = u(n);
	}
	else {
		uvp(1) = v*((pr-pl)*rhoi +lam1*(ur-ul))*den/(lam0-lam1) +vr;
		for(int n=tri_mesh::ND;n<x.NV-1;++n)
			uvp(n) = ub(n);
	}
	uvp(x.NV-1) = (rho*(ul -ur)*gam - lam2*pl +lam1*pr)*den;

	/* CHANGE BACK TO X,Y COORDINATES */
	ub(0) =  uvp(0)*norm(0) -uvp(1)*norm(1);
	ub(1) =  uvp(0)*norm(1) +uvp(1)*norm(0);

	for(int n=tri_mesh::ND;n<x.NV;++n)
		ub(n) =uvp(n);

	norm *= mag;

	flx(x.NV-1) = rho*((ub(0) -mv(0))*norm(0) +(ub(1) -mv(1))*norm(1));

	for(int n=0;n<tri_mesh::ND;++n)
		flx(n) = flx(x.NV-1)*ub(n) +ub(x.NV-1)*norm(n);

	for(int n=tri_mesh::ND;n<x.NV-1;++n)
		flx(n) = flx(x.NV-1)*ub(n);

//	*x.gbl->log << x.npnt << '\t' << u << '\t' << xpt << '\t' << mv << '\t' << norm << '\t' << flx << '\n';
//	*x.gbl->log << ul << ' ' << vl << ' ' << pl << ' ' << ur << ' ' << vr << ' ' << pr << ' '<< uvp << '\n';


	return;
}

void actuator_disc::output(std::ostream& fout, tri_hp::filetype typ,int tlvl) {
	int n,ind,sind;
	TinyVector<FLT,tri_mesh::ND> nrm, mvel, pt;
	FLT power;
	
	switch(typ) {
		case(tri_hp::tecplot): {
			if (!report_flag) return;
			
			power = 0.0;

			ind = 0;
			for (ind=0;ind<base.nseg;++ind) {
				sind = base.seg(ind);
			
				x.crdtocht1d(sind);
				for(n=0;n<tri_mesh::ND;++n)
					basis::tri(x.log2p)->proj1d(&x.cht(n,0),&x.crd(n)(0,0),&x.dcrd(n,0)(0,0));
				
				x.ugtouht1d(sind);
				for(n=0;n<x.NV;++n)
					basis::tri(x.log2p)->proj1d(&x.uht(n)(0),&x.u(n)(0,0));
				
				for(int k=0;k<basis::tri(x.log2p)->gpx();++k) {
					nrm(0) = x.dcrd(1,0)(0,k);
					nrm(1) = -x.dcrd(0,0)(0,k);                
					for(n=0;n<tri_mesh::ND;++n) {
						pt(n) = x.crd(n)(0,k);
						mvel(n) = x.gbl->bd(0)*(x.crd(n)(0,k) -dxdt(x.log2p,ind)(n,k));
					}
					FLT length = sqrt(nrm(0)*nrm(0) +nrm(1)*nrm(1));
					FLT norm_vel = ((x.u(0)(0,k) -mvel(0))*nrm(0) +(x.u(1)(0,k) -mvel(1))*nrm(1))/length;
					TinyVector<FLT,3> inpt(pt(0),pt(1),norm_vel);
					FLT delta_p = dp.Eval(inpt,x.gbl->time);			
					power +=  basis::tri(x.log2p)->wtx(k)*RAD(x.crd(0)(0,k))*length*delta_p*norm_vel;
				}
			}
			streamsize oldprecision = (*x.gbl->log).precision(10);
			*x.gbl->log << base.idprefix << " power: " << power << std::endl;
			(*x.gbl->log).precision(oldprecision);
			break;
		}
		default: {
			break;
		}
	}
	
	neumann::output(fout,typ,tlvl);
	
	return;
}


void hybrid_slave_pt::update(int stage) {

	if (stage == -1) return;

	int sendsize(6);
	base.sndsize() = sendsize;
	base.sndtype() = boundary::flt_msg;
	base.fsndbuf(0) = 0.0;
	base.fsndbuf(1) = 0.0;
	base.fsndbuf(2) = 0.0;
	base.fsndbuf(3) = 0.0;
	base.fsndbuf(4) = 0.0;
	base.fsndbuf(5) = 0.0;

	base.comm_prepare(boundary::all,0,boundary::symmetric);
	base.comm_exchange(boundary::all,0,boundary::symmetric);
	base.comm_wait(boundary::all,0,boundary::symmetric);

	for(int m=0;m<base.nmatches();++m) {
		for(int i=0;i<sendsize;++i) 
			base.fsndbuf(i) += base.frcvbuf(m,i);
	}

	if (base.fsndbuf(0)*base.fsndbuf(3) > 0.0) {
		*x.gbl->log << "uh-oh opposite characteristics at hybrid point" << std::endl;
		*x.gbl->log << "local " << base.idprefix << ' ' << base.fsndbuf(0) << "remote " << base.fsndbuf(3) << std::endl;
	}

	if (base.fsndbuf(0) > 0.0) {
		if (base.fsndbuf(4) == 2.0)
			x.pnts(base.pnt)(0) = base.fsndbuf(4) - 2.0;
		else 
			x.pnts(base.pnt)(0) = base.fsndbuf(4);
		x.pnts(base.pnt)(1) = base.fsndbuf(5);
	}
	else {
		x.pnts(base.pnt)(0) = base.fsndbuf(1);
		x.pnts(base.pnt)(1) = base.fsndbuf(2);
	}
	/* to move the other point along the edge back to the boundary */
	/*
	//std::cout<<"start moving points: SLAVE"<<endl;
	int nsides(4);
	for(int q=0;q<nsides;q++)
	{
		if(base.ebdry(q) != 0)
		{
			//std::cout<<"QQQQ: "<<q<<endl;
			//std::cout<<base.ebdry(0)<<endl;
			cout<<"NUMBER "<<q<<endl;
			for(int r=0;r<x.ebdry(q)->nseg;r++)
			{
				//std::cout<<"r: "<<r<<endl;
				x.ebdry(base.ebdry(q))->mvpttobdry(x.ebdry(base.ebdry(q))->seg(q),0.0,x.pnts(x.seg(x.ebdry(base.ebdry(q))->seg(r)).pnt(0)));
			}
		}
	}
	*/
	/*
	//std::cout<<"start moving points: SLAVE"<<endl;
	int nsides(4);
	for(int q=0;q<nsides;q++)
	{
		for(int r=0;r<x.ebdry(q)->nseg;r++)
			x.ebdry(q)->mvpttobdry(x.ebdry(q)->seg(r),0.0,x.pnts(x.seg(x.ebdry(q)->seg(r)).pnt(0)));
	}
	*/
	//std::cout<<"slave error: "<<x.pnts(base.pnt)(0)<<' '<<x.pnts(base.pnt)(1)<<' '<<x.pnts(base.pnt)(0)-0.1*sin(3.1415926*x.pnts(base.pnt)(1))-1<<std::endl;
}

void hybrid_pt::rsdl(int stage) {
	int sind,v0,v1;
	TinyVector<FLT,2> tang,vel;
	FLT tangvel;

	if (surfbdry == 0) {
		sind = x.ebdry(base.ebdry(0))->seg(x.ebdry(base.ebdry(0))->nseg-1);
		v0 = x.seg(sind).pnt(1);
		v1 = x.seg(sind).pnt(0);
	}
	else {
		sind = x.ebdry(base.ebdry(1))->seg(0);
		v0 = x.seg(sind).pnt(0);
		v1 = x.seg(sind).pnt(1);
	}


	/* TANGENT POINTS INTO DOMAIN ALONG SURFACE */
	tang(0) =  (x.pnts(v1)(0) -x.pnts(v0)(0));
	tang(1) =  (x.pnts(v1)(1) -x.pnts(v0)(1));

	vel(0) = 0.5*(x.ug.v(v0,0)-(x.gbl->bd(0)*(x.pnts(v0)(0) -x.vrtxbd(1)(v0)(0))) +
				  x.ug.v(v1,0)-(x.gbl->bd(0)*(x.pnts(v1)(0) -x.vrtxbd(1)(v1)(0))));
	vel(1) = 0.5*(x.ug.v(v0,1)-(x.gbl->bd(0)*(x.pnts(v0)(1) -x.vrtxbd(1)(v0)(1))) +
				  x.ug.v(v1,1)-(x.gbl->bd(0)*(x.pnts(v1)(1) -x.vrtxbd(1)(v1)(1))));
	tangvel = vel(0)*tang(0)+vel(1)*tang(1);

	if (tangvel > 0.0)
		fix_norm = 1;
	else
		fix_norm = 0;

	surface_outflow::rsdl(stage);
}

void hybrid_pt::update(int stage) {

	if (stage == -1) return;

	int sendsize(6);
	base.sndsize() = sendsize;
	base.sndtype() = boundary::flt_msg;
	// (0) = flow direction
	// (1) = x location
	// (2) = y location
	// see /lvlset/bdry.cpp hybrid_pt::update
	base.fsndbuf(0) = 2*fix_norm-1.0;
	base.fsndbuf(1) = x.pnts(base.pnt)(0);
	base.fsndbuf(2) = x.pnts(base.pnt)(1);
	base.fsndbuf(3) = 0.0;
	base.fsndbuf(4) = 0.0;
	base.fsndbuf(5) = 0.0;

	base.comm_prepare(boundary::all,0,boundary::symmetric);
	base.comm_exchange(boundary::all,0,boundary::symmetric);
	base.comm_wait(boundary::all,0,boundary::symmetric);

	// add information from other points
	for(int m=0;m<base.nmatches();++m) {
		for(int i=0;i<sendsize;++i) 
			base.fsndbuf(i) += base.frcvbuf(m,i);
	}

	if (base.fsndbuf(0)*base.fsndbuf(3) > 0.0) {
		*x.gbl->log << "uh-oh opposite characteristics at hybrid point" << std::endl;
		*x.gbl->log << "local "  << base.idprefix << ' ' << base.fsndbuf(0) << "remote " << base.fsndbuf(3) << std::endl;
	}
	if (fix_norm) {
		x.pnts(base.pnt)(0) = base.fsndbuf(4);
		x.pnts(base.pnt)(1) = base.fsndbuf(5);
		mvpttobdry(x.pnts(base.pnt));
	}
}

void actuator_disc::smatchsolution_snd(FLT *sdata, int bgn, int end, int stride) {
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


void actuator_disc::smatchsolution_rcv(FLT *sdata, int bgn, int end, int stride) {
	
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
	for(m=0;m<base.nmatches();++m) {            
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

#ifdef petsc
void actuator_disc::non_sparse_snd(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi) {
	
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
		
	Array<int,1> c0vars(vdofs-1);
	for(int n=0;n<x.NV-1;++n) {
		c0vars(n) = n;
	}
	for(int n=x.NV;n<vdofs;++n) {
		c0vars(n-1) = n;
	}		
	
	/* Going to send all jacobian entries,  Diagonal entries for matching DOF's will be merged together not individual */
	/* Send number of non-zeros to matches */
	base.sndsize() = 0;
	base.sndtype() = boundary::int_msg;
	for (int i=0;i<base.nseg;++i) {
		int sind = base.seg(i);
		int pind = x.seg(sind).pnt(0)*vdofs;
		for(int n=0;n<c0vars.extent(firstDim);++n)
			base.isndbuf(base.sndsize()++) = nnzero(pind+c0vars(n));
	}
	int sind = base.seg(base.nseg-1);
	int pind = x.seg(sind).pnt(1)*vdofs;
	for(int n=0;n<c0vars.extent(firstDim);++n)
		base.isndbuf(base.sndsize()++) = nnzero(pind+c0vars(n));
	
	/* Last thing to send is nnzero for edges (all the same) */
	if (sm)
		base.isndbuf(base.sndsize()++) = nnzero(begin_seg+sind*NV*sm);
	
	return;
}

void actuator_disc::non_sparse_rcv(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi) {
	
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
	
	Array<int,1> c0vars(vdofs-1);
	for(int n=0;n<x.NV-1;++n) {
		c0vars(n) = n;
	}
	for(int n=x.NV;n<vdofs;++n) {
		c0vars(n-1) = n;
	}	
	
	if (base.is_local(0)) {
		int count = 0;
		for (int i=base.nseg-1;i>=0;--i) {
			int sind = base.seg(i);
			int pind = x.seg(sind).pnt(1)*vdofs;
			for(int n=0;n<c0vars.extent(firstDim);++n) {
				nnzero(pind+c0vars(n)) += base.ircvbuf(0,count++);
			}
		}
		int sind = base.seg(0);
		int pind = x.seg(sind).pnt(0)*vdofs;
		for(int n=0;n<c0vars.extent(firstDim);++n) {
			nnzero(pind +c0vars(n)) += base.ircvbuf(0,count++);
		}
		
		
		/* Now add to side degrees of freedom */
		if (sm) {
			int toadd = base.ircvbuf(0,count++); 
			for (int i=0;i<base.nseg;++i) {
				int sind = base.seg(i);
				for (int mode=0;mode<sm;++mode) {
					for(int n=0;n<x.NV-1;++n) {
						nnzero(begin_seg+sind*NV*sm +mode*NV +n) += toadd;
					}
				}
			}
		}
	}
	else {
		int count = 0;
		for (int i=base.nseg-1;i>=0;--i) {
			int sind = base.seg(i);
			int pind = x.seg(sind).pnt(1)*vdofs;
			for(int n=0;n<c0vars.extent(firstDim);++n) {
				nnzero_mpi(pind+c0vars(n)) += base.ircvbuf(0,count++);
			}
		}
		int sind = base.seg(0);
		int pind = x.seg(sind).pnt(0)*vdofs;
		for(int n=0;n<c0vars.extent(firstDim);++n) {
			nnzero_mpi(pind +c0vars(n)) += base.ircvbuf(0,count++);
		}
		
		
		/* Now add to side degrees of freedom */
		if (sm) {
			int toadd = base.ircvbuf(0,count++); 
			for (int i=0;i<base.nseg;++i) {
				int sind = base.seg(i);
				for (int mode=0;mode<sm;++mode) {
					for(int n=0;n<x.NV-1;++n) {
						nnzero_mpi(begin_seg+sind*NV*sm +mode*NV +n) += toadd;
					}
				}
			}
		}
	}
}

void actuator_disc::petsc_matchjacobian_snd() {	
	int vdofs;
	if (x.mmovement != x.coupled_deformable)
		vdofs = x.NV;
	else
		vdofs = x.NV+x.ND;
	
	if (!base.is_comm()) return;
	
	/* Now do stuff for communication boundaries */
	int row,sind;
	
	Array<int,1> c0vars(vdofs-1);
	for(int n=0;n<x.NV-1;++n) {
		c0vars(n) = n;
	}
	for(int n=x.NV;n<vdofs;++n) {
		c0vars(n-1) = n;
	}		
	
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
		
		for (int n = 0; n <c0vars.extent(firstDim);++n) {
			row = rowbase + c0vars(n);
			base.fsndbuf(base.sndsize()++) = x.J._cpt(row+1) -x.J._cpt(row) +FLT_EPSILON;
#ifdef MPDEBUG
			*x.gbl->log << "sending " << x.J._cpt(row+1) -x.J._cpt(row) << " jacobian entries for vertex " << row/vdofs << " and variable " << c0vars(n) << std::endl;
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
			for (int n=0;n<x.NV-1;++n) {
				base.fsndbuf(base.sndsize()++) = x.J._cpt(row+1) -x.J._cpt(row) +FLT_EPSILON;
#ifdef MPDEBUG
				*x.gbl->log << "sending " << x.J._cpt(row+1) -x.J._cpt(row) << " jacobian entries for side " << sind << " and variable " << n << std::endl;
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
			++row; // Skip pressure
		}
	}
	
	/* LAST POINT */
	int rowbase = x.seg(sind).pnt(1)*vdofs; 
	
	/* attach diagonal # to allow continuity enforcement */
	base.fsndbuf(base.sndsize()++) = rowbase +FLT_EPSILON;
	
	for (int n = 0; n <c0vars.extent(firstDim);++n) {
		row = rowbase + c0vars(n);
		base.fsndbuf(base.sndsize()++) = x.J._cpt(row+1) -x.J._cpt(row) +FLT_EPSILON;
#ifdef MPDEBUG
		*x.gbl->log << "sending " << x.J._cpt(row+1) -x.J._cpt(row) << " jacobian entries for vertex " << row/vdofs << " and variable " << c0vars(n) << std::endl;
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
}

void actuator_disc::petsc_matchjacobian_rcv(int phase) {
	
	if (!base.is_comm() || base.matchphase(boundary::all_phased,0) != phase) return;
	
	int count = 0;
	int Jstart_mpi = base.frcvbuf(0, count++);
	
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
	
	/* Now do stuff for communication boundaries */
	int row;
	
	Array<int,1> c0vars(vdofs-1);
	for(int n=0;n<x.NV-1;++n) {
		c0vars(n) = n;
	}
	for(int n=x.NV;n<vdofs;++n) {
		c0vars(n-1) = n;
	}		
	
	/* Now Receive Information */		
	for (int i=base.nseg-1;i>=0;--i) {
		int sind = base.seg(i);
		int rowbase = x.seg(sind).pnt(1)*vdofs; 
		int row_mpi = base.frcvbuf(0,count++) +Jstart_mpi;
		
		for (int n = 0; n <c0vars.extent(firstDim);++n) {
			row = rowbase + c0vars(n);
			int ncol = base.frcvbuf(0,count++);
#ifdef MPDEBUG
			*x.gbl->log << "receiving " << ncol << " jacobian entries for vertex " << row/vdofs << " and variable " << c0vars(n) << std::endl;
#endif
			for (int k = 0;k<ncol;++k) {
				int col = base.frcvbuf(0,count++) +Jstart_mpi;
				FLT val = base.frcvbuf(0,count++);
				if (abs(col) < INT_MAX-10) {
#ifdef MPDEBUG
					*x.gbl->log  << col << ' ';
#endif
					(*pJ_mpi).add_values(row,col,val);
				}
			}
#ifdef MPDEBUG
			*x.gbl->log << std::endl;
#endif
			/* Shift all entries for this vertex */
			for (int n_mpi = 0; n_mpi <c0vars.extent(firstDim);++n_mpi) {
				FLT dval = (*pJ_mpi)(row,row_mpi+c0vars(n_mpi));
				(*pJ_mpi)(row,row_mpi+c0vars(n_mpi)) = 0.0;				
				x.J(row,rowbase+c0vars(n_mpi)) += dval;
			}
			x.J.multiply_row(row,0.5);
			x.J_mpi.multiply_row(row,0.5);
		}  
		
		/* Now receive side Jacobian information */
		row = x.npnt*vdofs +sind*x.NV*x.sm0;
		row_mpi = base.frcvbuf(0,count++) +Jstart_mpi;
		
		int mcnt = 0;
		int sgn = 1;
		for(int mode=0;mode<x.sm0;++mode) {
			for(int n=0;n<x.NV-1;++n) {
				int ncol = base.frcvbuf(0,count++);
#ifdef MPDEBUG
				*x.gbl->log << "receiving " << ncol << " jacobian entries for side " << sind << " and variable " << n << std::endl;
#endif
				for (int k = 0;k<ncol;++k) {
					int col = base.frcvbuf(0,count++) +Jstart_mpi;
					FLT val = sgn*base.frcvbuf(0,count++);
					if (abs(col) < INT_MAX-10) {
#ifdef MPDEBUG
						*x.gbl->log  << col << ' ';
#endif
						(*pJ_mpi).add_values(row+mcnt,col,val);
					}
				}
#ifdef MPDEBUG
				*x.gbl->log << std::endl;
#endif
				
				/* Shift all modes in equation */
				int mcnt_mpi = 0;
				int sgn_mpi = 1;
				for(int mode_mpi=0;mode_mpi<x.sm0;++mode_mpi) {
					for(int n_mpi = 0;n_mpi<x.NV-1;++n_mpi) {
						FLT dval = (*pJ_mpi)(row+mcnt,row_mpi+mcnt_mpi);
						(*pJ_mpi)(row+mcnt,row_mpi+mcnt_mpi) = 0.0;				
						x.J(row+mcnt,row+mcnt_mpi) += sgn_mpi*dval;
						++mcnt_mpi;
					}
					sgn_mpi *= -1;
					++mcnt_mpi; // Skip pressure
				}
				x.J.multiply_row(row+mcnt,0.5);
				x.J_mpi.multiply_row(row+mcnt,0.5);
				++mcnt;
			}
			++mcnt; // Skip pressure
			sgn *= -1;
		}
	}
	int sind = base.seg(0);
	int rowbase = x.seg(sind).pnt(0)*vdofs; 
	int row_mpi = base.frcvbuf(0,count++) +Jstart_mpi;
	
	for (int n = 0; n <c0vars.extent(firstDim);++n) {
		row = rowbase + c0vars(n);
		int ncol = base.frcvbuf(0,count++);
#ifdef MPDEBUG
		*x.gbl->log << "receiving " << ncol << " jacobian entries for vertex " << row/vdofs << " and variable " << c0vars(n) << std::endl;
#endif
		for (int k = 0;k<ncol;++k) {
			int col = base.frcvbuf(0,count++) +Jstart_mpi;
			FLT val = base.frcvbuf(0,count++);
			if (abs(col) < INT_MAX-10) {
#ifdef MPDEBUG
				*x.gbl->log  << col << ' ';
#endif
				(*pJ_mpi).add_values(row,col,val);
			}
		}
#ifdef MPDEBUG
		*x.gbl->log << std::endl;
#endif
		/* Shift all entries for this vertex */
		for (int n_mpi = 0; n_mpi <c0vars.extent(firstDim);++n_mpi) {
			FLT dval = (*pJ_mpi)(row,row_mpi+c0vars(n_mpi));
			(*pJ_mpi)(row,row_mpi+c0vars(n_mpi)) = 0.0;				
			x.J(row,rowbase+c0vars(n_mpi)) += dval;
		}
		x.J.multiply_row(row,0.5);
		x.J_mpi.multiply_row(row,0.5);
	} 
#endif
}	

#endif


