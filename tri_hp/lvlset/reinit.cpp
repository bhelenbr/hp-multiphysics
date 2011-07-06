#include "tri_hp_lvlset.h"
#include "../tri_hp.h"	
#include "../hp_boundary.h"
#include <math.h>
#include <blitz/tinyvec-et.h>

// #define DEBUG_TOL 1.0e-9
// #define RSDL_DEBUG
#define DEBUG_OUTPUT

void tri_hp_lvlset::reinitialize() {
	normstuff norm0list[6];
	normstuff norm1list[6];
	
	// The info gathered here is more than necessary for the implemented method, but it may be useful for another
	// First go through both end points of all boundaries and see if a bc is applied
	for(int i=0;i<nebd;++i){
		FLT nx, ny, xloc, yloc, phi, basephi;
		
		minvrt_reinit_direction(i, 0, nx, ny, xloc, yloc, basephi);
		if( !minvrt_reinit_phival(i, 0, nx, ny, xloc, yloc, basephi, 1, phi) ){
			norm0list[i].use = true;
		}
		else {
			norm0list[i].use = false;
		}
		minvrt_reinit_phival(i, 0, nx, ny, xloc, yloc, basephi, 0, phi); 
		norm0list[i].nx = nx;
		norm0list[i].ny = ny;
		norm0list[i].xloc = xloc;
		norm0list[i].yloc = yloc;
		norm0list[i].basephi = basephi;

		minvrt_reinit_direction(i, 1, nx, ny, xloc, yloc, basephi);
		if( !minvrt_reinit_phival(i, ebdry(i)->nseg-1, nx, ny, xloc, yloc, basephi, 0, phi) ){
			norm1list[i].use = true;
		}
		else {
			norm1list[i].use = false;
		}
		minvrt_reinit_phival(i, ebdry(i)->nseg-1, nx, ny, xloc, yloc, basephi, 1, phi);
		norm1list[i].nx = nx;
		norm1list[i].ny = ny;
		norm1list[i].xloc = xloc;
		norm1list[i].yloc = yloc;
		norm1list[i].basephi = basephi;

	}
	
#ifdef DEBUG_OUTPUT
	std::ostringstream nstr;
	std::string fname;
	std::stringstream s;
	nstr << gbl->tstep << "." << gbl->substep << std::flush;
#endif

	for (int ii=0; ii<reinit_iterations; ii++){
		setup_preconditioner_reinit();
		reinit(norm0list, norm1list);
		
#ifdef DEBUG_OUTPUT
		/* To output during reinitialization */
		s.str("");
		s<<ii;
		fname = "reinit" + nstr.str() + "_" + s.str() +"_" +gbl->idprefix;
		output(fname,tecplot);
#endif
	}
}

void tri_hp_lvlset::reinit(normstuff norm0list[], normstuff norm1list[]) {
	int i,m,k,n,indx,indx1;
	FLT cflalpha;

	/* STORE INITIAL VALUES FOR NSTAGE EXPLICIT SCHEME */
	gbl->ug0.v(Range(0,npnt-1),Range::all()) = ug.v(Range(0,npnt-1),Range::all());
	if (basis::tri(log2p)->sm()) {
		gbl->ug0.s(Range(0,nseg-1),Range(0,sm0-1),Range::all()) = ug.s(Range(0,nseg-1),Range::all(),Range::all());
		if (basis::tri(log2p)->im()) {
			gbl->ug0.i(Range(0,ntri-1),Range(0,im0-1),Range::all()) = ug.i(Range(0,ntri-1),Range::all(),Range::all());
		}
	}
	
	for (int stage = 0; stage < gbl->nstage; ++stage) {
		rsdl_reinit(stage);

	/*  for debugging
		for(i=0;i<npnt;++i) {
			*gbl->log << "rsdl reinit v: " << i << ' ';
			if (i==315) std::cout<<"problem: "<<gbl->res.v(i,3)<<std::endl;
			if (i==315){
				std::cout<<"This may break everything: "<<i+1<<' ' << std::endl;
				for (n=0;n<NV;++n) 
					std::cout << gbl->res.v(i+1,n)<<' ';
				std::cout<<std::endl;
			}
			for (n=0;n<NV;++n) 
				*gbl->log << gbl->res.v(i,n) << ' ';
			*gbl->log << '\n';
		}
		
		std::cout<<"nseg: "<<nseg<<std::endl;
		for(i=0;i<nseg;++i) {
			for(int m=0;m<basis::tri(log2p)->sm();++m) {
				*gbl->log << "rsdl reinit s: " << i << ' ' << m << ' '; 
				for(n=0;n<NV;++n)
					*gbl->log << gbl->res.s(i,m,n) << ' ';
				*gbl->log << '\n';
			}
		}
		
		std::cout<<"ntri: "<<ntri<<std::endl;
		for(i=0;i<ntri;++i) {
			for(int m=0;m<basis::tri(log2p)->im();++m) {
				*gbl->log << "rsdl reinit i: " << i << ' ' << m << ' ';
				for(n=0;n<NV;++n) 
					*gbl->log << gbl->res.i(i,m,n) << ' ';
				*gbl->log << '\n';
			}
		}
		*/
		
		/* SET RESDIUALS TO ZERO ON INCOMING CHARACTERISTIC BOUNDARIES */
		minvrt_reinit(norm0list, norm1list);
		
		
// TO DEBUG: copy and paste from /tri_hp/nstage/ void tri_hp::update  #ifdef DEBUG
		cflalpha = gbl->alpha(stage)*gbl->cfl(log2p);
		ug.v(Range(0,npnt-1),Range::all()) = gbl->ug0.v(Range(0,npnt-1),Range::all()) -cflalpha*gbl->res.v(Range(0,npnt-1),Range::all());

		if (basis::tri(log2p)->sm() > 0) {
			ug.s(Range(0,nseg-1),Range(0,basis::tri(log2p)->sm()-1),Range::all()) = gbl->ug0.s(Range(0,nseg-1),Range(0,basis::tri(log2p)->sm()-1),Range::all()) -cflalpha*gbl->res.s(Range(0,nseg-1),Range(0,basis::tri(log2p)->sm()-1),Range::all());

			if (basis::tri(log2p)->im() > 0) {

				for(i=0;i<ntri;++i) {
					indx = 0;
					indx1 = 0;
					for(m=1;m<basis::tri(log2p)->sm();++m) {
						for(k=0;k<basis::tri(log2p)->sm()-m;++k) {
							for(n=0;n<NV;++n) {
								ug.i(i,indx1,n) =  gbl->ug0.i(i,indx1,n) -cflalpha*gbl->res.i(i,indx,n);
							}
							++indx; ++indx1;
						}
						indx1 += sm0 -basis::tri(log2p)->sm();
					}
				}
			}
		}
	}
	FLT maxerr = 0.0;
	for(int i = 0;i<npnt;++i) {
		maxerr = MAX(maxerr,fabs(gbl->res.v(i,2)));
	}
	*gbl->log << "#reinit rsdl: " << maxerr << std::endl; 
}

void tri_hp_lvlset::rsdl_reinit(int stage) {
	/* ONLY NEED TO CALL FOR MOVEMENT BETWEEN MESHES INHERIT FROM THIS FOR SPECIFIC PHYSICS */
	FLT oneminusbeta = 1.0-gbl->beta(stage);
	gbl->res.v(Range(0,npnt-1),Range::all()) = 0.0;
	gbl->res_r.v(Range(0,npnt-1),Range::all()) *= oneminusbeta;

	if (basis::tri(log2p)->sm()) {
		gbl->res.s(Range(0,nseg-1),Range(0,basis::tri(log2p)->sm()-1),Range::all()) = 0.0;
		gbl->res_r.s(Range(0,nseg-1),Range(0,basis::tri(log2p)->sm()-1),Range::all()) *= oneminusbeta;

		if (basis::tri(log2p)->im()) {
			gbl->res.i(Range(0,ntri-1),Range(0,basis::tri(log2p)->im()-1),Range::all()) = 0.0;
			gbl->res_r.i(Range(0,ntri-1),Range(0,basis::tri(log2p)->im()-1),Range::all()) *= oneminusbeta;
		}
	}

	Array<TinyVector<FLT,MXTM>,1> lf_re(NV),lf_im(NV); 
#ifdef RSDL_DEBUG
	// for reinit debug
	for(int i=0;i<npnt;++i) {
		for(int n=0;n<NV;++n) {
			if (fabs(gbl->res.v(i,n)) > DEBUG_TOL) *gbl->log << gbl->res.v(i,n) << ' ';
			else *gbl->log << "0.0 ";
		}
		*gbl->log << '\n';
	}
#endif
	
	for(int tind = 0; tind<ntri;++tind) {
		
		/* LOAD SOLUTION COEFFICIENTS FOR THIS ELEMENT */
		ugtouht(tind);
		
		/* call rsdl for element */
		element_rsdl_reinit(tind,stage,uht,lf_re,lf_im);
		
		/* load imaginary local residual into global residual */
		for (int m = 0; m < basis::tri(log2p)->tm(); ++m) 
			for (int n = 0; n < NV; ++n)
				lf(n)(m) = lf_im(n)(m);
		
		lftog(tind,gbl->res);
		
		/* load real local residual into global residual */
		for (int m = 0; m < basis::tri(log2p)->tm(); ++m) 
			for (int n = 0; n < NV; ++n)
				lf(n)(m) = lf_re(n)(m);
		
		lftog(tind,gbl->res_r);
	}

	/* ADD IN VISCOUS/DISSIPATIVE FLUX */
	gbl->res.v(Range(0,npnt-1),Range::all()) += gbl->res_r.v(Range(0,npnt-1),Range::all());
	if (basis::tri(log2p)->sm() > 0) {
		gbl->res.s(Range(0,nseg-1),Range(0,basis::tri(log2p)->sm()-1),Range::all()) += gbl->res_r.s(Range(0,nseg-1),Range(0,basis::tri(log2p)->sm()-1),Range::all());          
		if (basis::tri(log2p)->im() > 0) {
			gbl->res.i(Range(0,ntri-1),Range(0,basis::tri(log2p)->im()-1),Range::all()) += gbl->res_r.i(Range(0,ntri-1),Range(0,basis::tri(log2p)->im()-1),Range::all());      
		}
	}

#ifdef RSDL_DEBUG
// for reinit debug {
	int i, n;
//	if (coarse_flag) {
		for(i=0;i<npnt;++i) {
			*gbl->log << "rsdl reinit v: " << i << ' ';
			for (n=0;n<NV;++n) 
				*gbl->log << gbl->res.v(i,n) << ' ';
			*gbl->log << '\n';
		}
		
		for(i=0;i<nseg;++i) {
			for(int m=0;m<basis::tri(log2p)->sm();++m) {
				*gbl->log << "rsdl reinit s: " << i << ' ' << m << ' '; 
				for(n=0;n<NV;++n)
					*gbl->log << gbl->res.s(i,m,n) << ' ';
				*gbl->log << '\n';
			}
		}
		
		for(i=0;i<ntri;++i) {
			for(int m=0;m<basis::tri(log2p)->im();++m) {
				*gbl->log << "rsdl reinit i: " << i << ' ' << m << ' ';
				for(n=0;n<NV;++n) 
					*gbl->log << gbl->res.i(i,m,n) << ' ';
				*gbl->log << '\n';
			}
		}
//	}
//	}
#endif
}

void tri_hp_lvlset::element_rsdl_reinit(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uht,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im) {
	int i,j,n;
	const int NV = 4;
	TinyVector<int,3> v;
	TinyMatrix<FLT,ND,ND> ldcrd;
	TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,NV,ND> du;
	int lgpx = basis::tri(log2p)->gpx(), lgpn = basis::tri(log2p)->gpn();
	FLT  cjcb, oneminusbeta;
	TinyMatrix<TinyMatrix<FLT,ND,ND>,NV-1,NV-1> visc;
	TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,NV-1,NV-1> cv, df;
	TinyVector<FLT,NV> tres;
	TinyMatrix<FLT,MXGP,MXGP> rho, mu;
	FLT length, phidw, signphi;
	TinyVector<FLT,ND> tang,norm;
	TinyVector<TinyMatrix<FLT,MXGP,MXGP>,ND> phivel;

	oneminusbeta = 1.0-gbl->beta(stage);
	// LOAD INDICES OF VERTEX POINTS //
	v = tri(tind).pnt;

	// IF TINFO > -1 IT IS CURVED ELEMENT //
	if (tri(tind).info > -1) {
		// LOAD ISOPARAMETRIC MAPPING COEFFICIENTS //
		crdtocht(tind);

		// PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS //
		for(n=0;n<ND;++n)
			basis::tri(log2p)->proj_bdry(&cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), &dcrd(n,1)(0,0),MXGP);
	}
	else {
		// PROJECT VERTEX COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS //
		for(n=0;n<ND;++n)
			basis::tri(log2p)->proj(pnts(v(0))(n),pnts(v(1))(n),pnts(v(2))(n),&crd(n)(0,0),MXGP);

		// CALCULATE COORDINATE DERIVATIVES A SIMPLE WAY //
		for(n=0;n<ND;++n) {
			ldcrd(n,0) = 0.5*(pnts(v(2))(n) -pnts(v(1))(n));
			ldcrd(n,1) = 0.5*(pnts(v(0))(n) -pnts(v(1))(n));
		}
	}

	// LOAD SOLUTION COEFFICIENTS FOR THIS ELEMENT //
	// PROJECT SOLUTION TO GAUSS POINTS WITH DERIVATIVES IF NEEDED FOR VISCOUS TERMS //
	ugtouht(tind);
	basis::tri(log2p)->proj(&uht(NV-2)(0),&u(NV-2)(0,0),&du(NV-2,0)(0,0),&du(NV-2,1)(0,0),MXGP);

	// lf IS WHERE I WILL STORE THE ELEMENT RESIDUAL //
	for(n=0;n<NV;++n){
		for(i=0;i<basis::tri(log2p)->tm();++i){
			lf_re(n)(i) = 0.0;
			lf_im(n)(i) = 0.0;	
		}
	}
		for(i=0;i<lgpx;++i) {
			for(j=0;j<lgpn;++j) {
				// STUFF FOR LEVEL SET //
				phidw = u(2)(i,j)/(gbl->width);
				cjcb = ldcrd(0,0)*ldcrd(1,1) -ldcrd(1,0)*ldcrd(0,1);
				norm(0) = (+ldcrd(1,1)*du(2,0)(i,j) -ldcrd(1,0)*du(2,1)(i,j))/cjcb;
				norm(1) = (-ldcrd(0,1)*du(2,0)(i,j) +ldcrd(0,0)*du(2,1)(i,j))/cjcb;
				length = sqrt(norm(0)*norm(0) +norm(1)*norm(1));
				
				if (phidw < -1.0) {
					res(2)(i,j) = -RAD(crd(0)(i,j))*cjcb*(length-1.0);
					phivel(0)(i,j) = -norm(0)/length;
					phivel(1)(i,j) = -norm(1)/length;
				}
				else if (phidw > 1.0) {
					res(2)(i,j) = RAD(crd(0)(i,j))*cjcb*(length-1.0);
					phivel(0)(i,j) = norm(0)/length;
					phivel(1)(i,j) = norm(1)/length;
				}
				else {
					// signphi = phidw;
					// signphi = u(2)(i,j) / pow( pow(u(2)(i,j),2) + pow(1-length,2) * pow(gbl->width/2.0,2), 0.5 );
					// signphi = u(2)(i,j) / pow( pow(u(2)(i,j),2) + pow(1-phidw,2) * pow(gbl->width/2.0,2), 0.5 );
					// signphi = u(2)(i,j) / pow( pow(u(2)(i,j),2) + pow(length,2) * pow(gbl->width,2), 0.5 );
					signphi = sin(M_PI*phidw/2.0);
					phivel(0)(i,j) = signphi/(fabs(signphi)+EPSILON)*norm(0)/length;
					phivel(1)(i,j) = signphi/(fabs(signphi)+EPSILON)*norm(1)/length;
					//res(2)(i,j) = RAD(crd(0)(i,j))*cjcb*(length-1.0)*signphi;
					res(2)(i,j) = RAD(crd(0)(i,j))*cjcb*(length-1.0)*signphi;

				}
			}
		}
		basis::tri(log2p)->intgrt(&lf_im(NV-2)(0),&res(NV-2)(0,0),MXGP);

		// NEGATIVE REAL TERMS //
		if (gbl->beta(stage) > 0.0) {
			for(n=0;n<basis::tri(log2p)->tm();++n)
				lf_re(NV-2)(n) = 0.0;

			// THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES //
			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpn;++j) {
					tres(2) = gbl->tau(tind,2)*res(2)(i,j);

					df(2,0)(i,j) = -(ldcrd(1,1)*phivel(0)(i,j) -ldcrd(0,1)*phivel(1)(i,j))*tres(2);
					df(2,1)(i,j) = -(-ldcrd(1,0)*phivel(0)(i,j) +ldcrd(0,0)*phivel(1)(i,j))*tres(2);
				}
			}
			basis::tri(log2p)->intgrtrs(&lf_re(2)(0),&df(2,0)(0,0),&df(2,1)(0,0),MXGP);
			
			for(i=0;i<basis::tri(log2p)->tm();++i)
				lf_re(2)(i) *= gbl->beta(stage);

		}

	return;
}

void tri_hp_lvlset::setup_preconditioner_reinit() {

	if (gbl->diagonal_preconditioner) {
		// set-up diagonal preconditioner //
		int tind,i,j,side;
		FLT jcb,jcbphi,h,hmax,lam2;
		TinyVector<int,3> v;

		/************************************/
		/* determine flow pseudo-time step **/
		/************************************/
		gbl->vprcn(Range(0,npnt-1),Range::all()) = 0.0;
		if (basis::tri(log2p)->sm() > 0) {
			gbl->sprcn(Range(0,nseg-1),Range::all()) = 0.0;
		}

		for(tind = 0; tind < ntri; ++tind) {
			jcb = 0.25*area(tind);  // area is 2 x triangle area
			v = tri(tind).pnt;
			hmax = 0.0;
			for(j=0;j<3;++j) {
				h = pow(pnts(v(j))(0) -pnts(v((j+1)%3))(0),2.0) + 
				pow(pnts(v(j))(1) -pnts(v((j+1)%3))(1),2.0);
				hmax = (h > hmax ? h : hmax);
			}
			hmax = sqrt(hmax);

			if (!(jcb > 0.0)) {  // this catches nan's too
				*gbl->log << "negative triangle area caught in tstep. problem triangle is : " << tind << std::endl;
				*gbl->log << "approximate location: " << pnts(v(0))(0) << ' ' << pnts(v(0))(1) << std::endl;
				tri_mesh::output("negative",grid);
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			h = 4.*jcb/(0.25*(basis::tri(log2p)->p() +1)*(basis::tri(log2p)->p()+1)*hmax);
			hmax = hmax/(0.25*(basis::tri(log2p)->p() +1)*(basis::tri(log2p)->p()+1));

			gbl->tau(tind,2) = adis*h/jcb;
			lam2 = 2.0/h;

			// set up diagonal preconditioner //
			jcbphi = 2.0*jcb*lam2;
			jcbphi *= RAD((pnts(v(0))(0) +pnts(v(1))(0) +pnts(v(2))(0))/3.);

			jcb *= 1.0e99*RAD((pnts(v(0))(0) +pnts(v(1))(0) +pnts(v(2))(0))/3.);

			// brute force way to freeze flow //
			gbl->tprcn(tind,0) = jcb;   
			gbl->tprcn(tind,1) = jcb;   
			gbl->tprcn(tind,3) =  jcb; 

			gbl->tprcn(tind,2) = jcbphi;      
			for(i=0;i<3;++i) {
				gbl->vprcn(v(i),Range::all())  += gbl->tprcn(tind,Range::all());
				if (basis::tri(log2p)->sm() > 0) {
					side = tri(tind).seg(i);
					gbl->sprcn(side,Range::all()) += gbl->tprcn(tind,Range::all());
				}
			}
		}
	}
	/* PREINVERT PRECONDITIONER FOR VERTICES */
	gbl->vprcn(Range(0,npnt-1),Range::all()) = 1.0/(basis::tri(log2p)->vdiag()*gbl->vprcn(Range(0,npnt-1),Range::all()));

	if (log2p) {
		/* INVERT DIAGANOL PRECONDITIONER FOR SIDES */ 
		gbl->sprcn(Range(0,nseg-1),Range::all()) = 1.0/gbl->sprcn(Range(0,nseg-1),Range::all());
	}

	return;
}

void tri_hp_lvlset::minvrt_reinit(normstuff norm0list[], normstuff norm1list[]) {
	int i,j(0),k,m,tind,sind,v0,indx,indx1,indx2,sgn,msgn;
	TinyVector<int,3> sign,side;
	Array<FLT,2> tinv(NV,NV);
	Array<FLT,1> temp(NV);
	FLT basephi;

	/* LOOP THROUGH SIDES */
	if (basis::tri(log2p)->sm() > 0) {
		indx = 0;
		for(sind = 0; sind<nseg;++sind) {
			/* SUBTRACT SIDE CONTRIBUTIONS TO VERTICES */            
			for (k=0; k <basis::tri(log2p)->sm(); ++k) {
				for (i=0; i<2; ++i) {
					v0 = seg(sind).pnt(i);
					
						gbl->res.v(v0,2) -= basis::tri(log2p)->sfmv(i,k)*gbl->res.s(sind,k,2);
				}
				++indx;
			}
		}
		
		if (basis::tri(log2p)->im() > 0) {
			/* SUBTRACT INTERIORS */
			indx = 0;
			for(tind = 0; tind<ntri;++tind) {
				indx2 = 3;
				for (i=0; i<3; ++i) {
					v0 = tri(tind).pnt(i);
					for (k=0;k<basis::tri(log2p)->im();++k)
						
							gbl->res.v(v0,2) -= basis::tri(log2p)->ifmb(i,k)*gbl->res.i(tind,k,2);
					
					sind = tri(tind).seg(i);
					sgn = tri(tind).sgn(i);
					msgn = 1;
					for (j=0;j<basis::tri(log2p)->sm();++j) {
						for (k=0;k<basis::tri(log2p)->im();++k)
							
								gbl->res.s(sind,j,2) -= msgn*basis::tri(log2p)->ifmb(indx2,k)*gbl->res.i(tind,k,2);
						msgn *= sgn;
						++indx2;
					}
				}
				indx += basis::tri(log2p)->im();
			}
		}
	}
	
	gbl->res.v(Range(0,npnt-1),Range::all()) *= gbl->vprcn(Range(0,npnt-1),Range::all());
	
//	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
//		vc0load(mp_phase,gbl->res.v.data());
//		pmsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
//		last_phase = true;
//		last_phase &= vc0wait_rcv(mp_phase,gbl->res.v.data());
//	}
	

	/* APPLY VERTEX DIRICHLET B.C.'S */
	for(i=0;i<nebd;++i) {
		FLT nx(0), ny(0), xloc(0), yloc(0), phi(0);
		//FLT nx2, ny2;
		int bdryside( minvrt_reinit_phiside(i) );
		
		// get the correct info
		// these conditions
		if (bdryside==1){
			if ( (norm1list[i].use) || (!norm1list[i].use && !norm0list[(i+1)%nebd].use) ){
				nx = norm1list[i].nx;
				ny = norm1list[i].ny;
				xloc = norm1list[i].xloc;
				yloc = norm1list[i].yloc;
				basephi = norm1list[i].basephi;
			}
			else {
				nx = norm0list[(i+1)%nebd].nx;
				ny = norm0list[(i+1)%nebd].ny;
				xloc = norm0list[(i+1)%nebd].xloc;
				yloc = norm0list[(i+1)%nebd].yloc;
				basephi = norm0list[(i+1)%nebd].basephi;
			}
		}
		if (bdryside==0){
			if ( (norm0list[i].use) || (!norm0list[i].use && !norm1list[(i-1)%nebd].use) ){
				nx = norm0list[i].nx;
				ny = norm0list[i].ny;
				xloc = norm0list[i].xloc;
				yloc = norm0list[i].yloc;
				basephi = norm0list[i].basephi;
			}
			else {
				nx = norm1list[(i-1)%nebd].nx;
				ny = norm1list[(i-1)%nebd].ny;
				xloc = norm1list[(i-1)%nebd].xloc;
				yloc = norm1list[(i-1)%nebd].yloc;
				basephi = norm1list[(i-1)%nebd].basephi;
			}
		}
		// this is to set the first point along the boundary ...
		if( minvrt_reinit_phival(i, 0, nx, ny, xloc, yloc, basephi, 0, phi) ) {
			gbl->res.v(seg(ebdry(i)->seg(0)).pnt(0),2) = 0.0;
			ug.v(seg(ebdry(i)->seg(0)).pnt(0),2) = phi;
			gbl->ug0.v(seg(ebdry(i)->seg(0)).pnt(0),2) = phi;
		}
		// ... and all the rest
		for (int j=0;j<ebdry(i)->nseg;++j) {
			 /* SET VERTEX RESIDUALS TO ZERO IF INCOMING */
		  	if( minvrt_reinit_phival(i, j, nx, ny, xloc, yloc, basephi, 1, phi) ) {
				gbl->res.v(seg(ebdry(i)->seg(j)).pnt(1),2) = 0.0;
				ug.v(seg(ebdry(i)->seg(j)).pnt(1),2) = phi;
				gbl->ug0.v(seg(ebdry(i)->seg(j)).pnt(1),2) = phi;
			}
		}
	}

	//exit(1);

	if(basis::tri(log2p)->sm() == 0) return;

	/* REMOVE VERTEX CONTRIBUTION FROM SIDE MODES */
	/* SOLVE FOR SIDE MODES */
	/* PART 1 REMOVE VERTEX CONTRIBUTIONS */
	for(tind=0;tind<ntri;++tind) {
		for(i=0;i<3;++i) {
			v0 = tri(tind).pnt(i);
			
				uht(2)(i) = gbl->res.v(v0,2)*gbl->tprcn(tind,2);
		}
		
		for(i=0;i<3;++i) {
			sind = tri(tind).seg(i);
			sgn  = tri(tind).sgn(i);
			for(j=0;j<3;++j) {
				indx1 = (i+j)%3;
				msgn = 1;
				for(k=0;k<basis::tri(log2p)->sm();++k) {
					
						gbl->res.s(sind,k,2) -= msgn*basis::tri(log2p)->vfms(j,k)*uht(2)(indx1);
					msgn *= sgn;
				}
			}
		}
	}
	
	
	for (int mode = 0; mode < basis::tri(log2p)->sm()-1; ++ mode) {
		/* SOLVE FOR SIDE MODE */
		gbl->res.s(Range(0,nseg-1),mode,Range::all()) *= gbl->sprcn(Range(0,nseg-1),Range::all())*basis::tri(log2p)->sdiag(mode);
		
		/*
		sc0load(gbl->res.s.data(),mode,mode,gbl->res.s.extent(secondDim));
		smsgpass(boundary::all,0,boundary::symmetric);
		sc0wait_rcv(gbl->res.s.data(),mode,mode,gbl->res.s.extent(secondDim));
		*/
		
		/* APPLY SIDE DIRICHLET B.C.'S */
		for(i=0;i<nebd;++i) {
			FLT nx(0), ny(0), xloc(0), yloc(0), phi(0);
			//FLT nx2, ny2;
			int bdryside( minvrt_reinit_phiside(i) );
			
			// get the correct info
			// these conditions
			if (bdryside==1){
				if ( (norm1list[i].use) || (!norm1list[i].use && !norm0list[(i+1)%nebd].use) ){
					nx = norm1list[i].nx;
					ny = norm1list[i].ny;
					xloc = norm1list[i].xloc;
					yloc = norm1list[i].yloc;
					basephi = norm1list[i].basephi;
				}
				else {
					nx = norm0list[(i+1)%nebd].nx;
					ny = norm0list[(i+1)%nebd].ny;
					xloc = norm0list[(i+1)%nebd].xloc;
					yloc = norm0list[(i+1)%nebd].yloc;
					basephi = norm0list[(i+1)%nebd].basephi;
				}
			}
			if (bdryside==0){
				if ( (norm0list[i].use) || (!norm0list[i].use && !norm1list[(i-1)%nebd].use) ){
					nx = norm0list[i].nx;
					ny = norm0list[i].ny;
					xloc = norm0list[i].xloc;
					yloc = norm0list[i].yloc;
					basephi = norm0list[i].basephi;
				}
				else {
					nx = norm1list[(i-1)%nebd].nx;
					ny = norm1list[(i-1)%nebd].ny;
					xloc = norm1list[(i-1)%nebd].xloc;
					yloc = norm1list[(i-1)%nebd].yloc;
					basephi = norm1list[(i-1)%nebd].basephi;
				}
			}
			// ... and all the rest
			for (int j=0;j<ebdry(i)->nseg;++j) {
				 /* SET SIDE RESIDUALS TO ZERO IF INCOMING and MAKE SOLUTION TO LINEAR */
			  	if( minvrt_reinit_phival(i, j, nx, ny, xloc, yloc, basephi, 1, phi) ) {
					gbl->res.s(ebdry(i)->seg(j),mode,2) = 0.0;
					ug.s(ebdry(i)->seg(j),mode,2) = 0.0;
					gbl->ug0.s(ebdry(i)->seg(j),mode,2) = 0.0;
				}
			}
		}
		
		/* REMOVE MODE FROM HIGHER MODES */
		for(tind=0;tind<ntri;++tind) {
			for(i=0;i<3;++i) {
				side(i) = tri(tind).seg(i);
				sign(i) = tri(tind).sgn(i);
				sgn      = (mode % 2 ? sign(i) : 1);
				
					uht(2)(i) = sgn*gbl->res.s(side(i),mode,2)*gbl->tprcn(tind,2);
			}
			
			/* REMOVE MODES J,K FROM MODE I,M */
			for(i=0;i<3;++i) {
				msgn = ( (mode +1) % 2 ? sign(i) : 1);
				for(m=mode+1;m<basis::tri(log2p)->sm();++m) {
					for(j=0;j<3;++j) {
						indx = (i+j)%3;
						 {
							gbl->res.s(side(i),m,2) -= msgn*basis::tri(log2p)->sfms(mode,m,j)*uht(2)(indx);
						}
					}
					msgn *= sign(i);
				}
			}
		}
	}
	/* SOLVE FOR HIGHEST MODE */
	int mode = basis::tri(log2p)->sm()-1;
	gbl->res.s(Range(0,nseg-1),mode,Range::all()) *= gbl->sprcn(Range(0,nseg-1),Range::all())*basis::tri(log2p)->sdiag(mode);

	/*
	sc0load(gbl->res.s.data(),mode,mode,gbl->res.s.extent(secondDim));
	smsgpass(boundary::all,0,boundary::symmetric);
	sc0wait_rcv(gbl->res.s.data(),mode,mode,gbl->res.s.extent(secondDim));
	*/
	
	/* APPLY VERTEX DIRICHLET B.C.'S */
	for(i=0;i<nebd;++i) {
		FLT nx(0), ny(0), xloc(0), yloc(0), phi(0);
		//FLT nx2, ny2;
		int bdryside( minvrt_reinit_phiside(i) );
		
		// get the correct info
		// these conditions
		if (bdryside==1){
			if ( (norm1list[i].use) || (!norm1list[i].use && !norm0list[(i+1)%nebd].use) ){
				nx = norm1list[i].nx;
				ny = norm1list[i].ny;
				xloc = norm1list[i].xloc;
				yloc = norm1list[i].yloc;
				basephi = norm1list[i].basephi;
			}
			else {
				nx = norm0list[(i+1)%nebd].nx;
				ny = norm0list[(i+1)%nebd].ny;
				xloc = norm0list[(i+1)%nebd].xloc;
				yloc = norm0list[(i+1)%nebd].yloc;
				basephi = norm0list[(i+1)%nebd].basephi;
			}
		}
		if (bdryside==0){
			if ( (norm0list[i].use) || (!norm0list[i].use && !norm1list[(i-1)%nebd].use) ){
				nx = norm0list[i].nx;
				ny = norm0list[i].ny;
				xloc = norm0list[i].xloc;
				yloc = norm0list[i].yloc;
				basephi = norm0list[i].basephi;
			}
			else {
				nx = norm1list[(i-1)%nebd].nx;
				ny = norm1list[(i-1)%nebd].ny;
				xloc = norm1list[(i-1)%nebd].xloc;
				yloc = norm1list[(i-1)%nebd].yloc;
				basephi = norm1list[(i-1)%nebd].basephi;
			}
		}
		// ... and all the rest
		for (int j=0;j<ebdry(i)->nseg;++j) {
			/* SET SIDE RESIDUALS TO ZERO IF INCOMING and MAKE SOLUTION TO LINEAR */
			if( minvrt_reinit_phival(i, j, nx, ny, xloc, yloc, basephi, 1, phi) ) {
				gbl->res.s(ebdry(i)->seg(j),mode,2) = 0.0;
				ug.s(ebdry(i)->seg(j),mode,2) = 0.0;
				gbl->ug0.s(ebdry(i)->seg(j),mode,2) = 0.0;
			}
		}
	}
	
	if (basis::tri(log2p)->im() == 0) return;
	
	/* SOLVE FOR INTERIOR MODES */
	for(tind = 0; tind < ntri; ++tind) {
		DPBTRSNU2((double *) &basis::tri(log2p)->idiag(0,0),basis::tri(log2p)->ibwth()+1,basis::tri(log2p)->im(),basis::tri(log2p)->ibwth(),&(gbl->res.i(tind,0,0)),NV);
		restouht_bdry(tind);
		for(k=0;k<basis::tri(log2p)->im();++k) {
			gbl->res.i(tind,k,Range::all()) /= gbl->tprcn(tind,Range::all());
			
			for (i=0;i<basis::tri(log2p)->bm();++i)
				 
					gbl->res.i(tind,k,2) -= basis::tri(log2p)->bfmi(i,k)*uht(2)(i);
		}
	}
	
	return;
}

int tri_hp_lvlset::minvrt_reinit_phiside(int edgenum){
	FLT vone, vtwo;
	int nseg(ebdry(edgenum)->nseg);
	// find which side is closer to zero level
	vone = ug.v(seg(ebdry(edgenum)->seg(0)).pnt(0),2);
	vtwo = ug.v(seg(ebdry(edgenum)->seg(nseg-1)).pnt(1),2);
	/* COMPARE MAGNITUDE OF PHI AT ENDPOINTS */
	if (fabs(vone) < fabs(vtwo) )
		return 0;
	else 
		return 1;
}

void tri_hp_lvlset::minvrt_reinit_direction(int edgenum, int side, FLT&nx, FLT&ny, FLT&xloc, FLT&yloc, FLT&basephi) {
	int sind(0),tind(0),v0(0);
	FLT vone, vtwo;
	int nseg(ebdry(edgenum)->nseg);
	vone = ug.v(seg(ebdry(edgenum)->seg(0)).pnt(0),2);
	vtwo = ug.v(seg(ebdry(edgenum)->seg(nseg-1)).pnt(1),2);
	FLT vphi;
	if ( side ==0 ){
		vphi = vone;
		sind = ebdry(edgenum)->seg(0);
		tind = seg(sind).tri(0); 
		v0 = seg(sind).pnt(0);
	}
	else {
		vphi = vtwo;
		sind = ebdry(edgenum)->seg(nseg-1);
		tind = seg(sind).tri(0); 
		v0 = seg(sind).pnt(1);
	}
	basephi = vphi;
	xloc = pnts(v0)(0);
	yloc = pnts(v0)(1);
	
	/* AT WHAT POINT OF TRIANGLE DO WE NEED TO KNOW DERIVATIVE? */
	int j;
	for (j=0;j<3;++j)
		if (tri(tind).pnt(j) == v0) break;
	assert(j < 3);
	FLT r,s;
	// from point number to r,s
	r=-1.0;
	s=-1.0;
	if (j==0) s=1.0;
	if (j==2) r=1.0;	
	/* Figure out what r & s should be based on vertex # */
	ugtouht(tind);
	crdtocht(tind);
	// 2 -> ND, 4 -> NV
	TinyVector<FLT,2> xpt,dxdr,dxds,dxmax;
	TinyVector<FLT,4> upt,dudr,duds;
	basis::tri(log2p)->ptprobe_bdry(ND,xpt.data(),dxdr.data(),dxds.data(),r,s,&cht(0,0),MXTM); 
	basis::tri(log2p)->ptprobe(NV,upt.data(),dudr.data(),duds.data(),r,s,&uht(0)(0),MXTM); 
	TinyMatrix<FLT,2,2> ldcrd;
	ldcrd(0,0) = dxdr(0);
	ldcrd(0,1) = dxds(0);
	ldcrd(1,0) = dxdr(1);
	ldcrd(1,1) = dxds(1);
	FLT cjcb = ldcrd(0,0)*ldcrd(1,1) -ldcrd(1,0)*ldcrd(0,1);
	FLT normux, normuy, lengthu;
	normux = (+ldcrd(1,1)*dudr(2) -ldcrd(1,0)*duds(2))/cjcb;
	normuy = (-ldcrd(0,1)*dudr(2) +ldcrd(0,0)*duds(2))/cjcb;
	lengthu = sqrt(normux*normux +normuy*normuy); //for projection in this direction
	nx = normux / lengthu; // phi unit normal vector x direction (in the positive phi direction)
	ny = normuy / lengthu; // phi y direction
}

bool tri_hp_lvlset::minvrt_reinit_phival(int edgenum, int segnum, FLT nx, FLT ny, FLT xloc, FLT yloc, FLT basephi, int onezero, FLT& phi) {
	int sind(0),v0(0);
	sind = ebdry(edgenum)->seg(segnum);
	v0 = seg(sind).pnt(0);
	int v1 = seg(sind).pnt(1);
	FLT xdir, ydir;
	// a vector pointing along the segment
	xdir = pnts(v0)(0)-pnts(v1)(0);
	ydir = pnts(v0)(1)-pnts(v1)(1);
	// rotate the vector so it is pointing OUT of the domain
	FLT voutx(-ydir), vouty(xdir);
	// the check will be opposite if phi is + or -
	if( (ug.v(seg(ebdry(edgenum)->seg(segnum)).pnt(onezero),2)) > 0 ) {
		// checked if phi is positive
		if ( (nx*voutx + ny*vouty) < 0.0  ){ // check if flow in
			if (onezero == 0)
				phi = (pnts(v0)(0)-xloc)*nx + (pnts(v0)(1)-yloc)*ny + basephi;
			if (onezero == 1)
				phi = (pnts(v1)(0)-xloc)*nx + (pnts(v1)(1)-yloc)*ny + basephi;
			return true;
		}
	}
	else { // phi is negative
		if ( (nx*voutx + ny*vouty) > 0.0 ) {
			if (onezero == 0)
				phi = (pnts(v0)(0)-xloc)*nx + (pnts(v0)(1)-yloc)*ny + basephi;
			if (onezero == 1)
				phi = (pnts(v1)(0)-xloc)*nx + (pnts(v1)(1)-yloc)*ny + basephi;
			return true;
		}
	}
	return false;
}
