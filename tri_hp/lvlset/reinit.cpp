#include "tri_hp_lvlset.h"
#include "../tri_hp.h"	
#include "../hp_boundary.h"
#include <math.h>
#include <blitz/tinyvec-et.h>

#ifdef REINIT

#define DEBUG_TOL 1.0e-9
//#define RSDL_DEBUG
void tri_hp_lvlset::reinitialize() {
	std::cout<<"REINITIALIZING LEVEL SET"<<std::endl;
	for (int ii=0; ii<4; ii++){
		setup_preconditioner_reinit();
//		std::cout<<"finished preconditioner"<<std::endl;
		reinit();
//		std::cout<<"finished reinit"<<std::endl;
	}
	std::cout<<"finished reinitializing"<<std::endl;
}

void tri_hp_lvlset::reinit() {
	int i,m,k,n,indx,indx1;
	FLT cflalpha;

//	std::cout<<"::reinit() 1"<<std::endl;
	/* STORE INITIAL VALUES FOR NSTAGE EXPLICIT SCHEME */
	gbl->ug0.v(Range(0,npnt-1),Range::all()) = ug.v(Range(0,npnt-1),Range::all());
	if (basis::tri(log2p)->sm()) {
		gbl->ug0.s(Range(0,nseg-1),Range(0,sm0-1),Range::all()) = ug.s(Range(0,nseg-1),Range::all(),Range::all());
		if (basis::tri(log2p)->im()) {
			gbl->ug0.i(Range(0,ntri-1),Range(0,im0-1),Range::all()) = ug.i(Range(0,ntri-1),Range::all(),Range::all());
		}
	}
//	std::cout<<"::reinit () 2"<<std::endl;
//	std::cout<<"there should be "<<gbl->nstage<<" stages"<<std::endl;
	
	for (int stage = 0; stage < gbl->nstage; ++stage) {
		rsdl_reinit(stage);
//		std::cout<<"::reinit () did rsdl_reinit("<<stage<<")"<<std::endl;

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
		minvrt();
//		std::cout<<"::reinit () did minvrt()"<<std::endl;
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
	/*
	static int display = 0;
	{
		char name [10];
		sprintf(name, "%s%d", "Reinit", display++);
		output(name,tecplot);
		std::cout<<"reinitialized "<<display++<<std::endl;
	}
	*/
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
		*gbl->log << gbl->idprefix << "rsdl_reinit 5 v2: " << i << ' ';
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
	
	/* MODIFY RESIDUALS ON COARSER MESHES */
	if(coarse_flag) {
		/* CALCULATE DRIVING TERM ON FIRST ENTRY TO COARSE MESH */
		if(isfrst) {
			dres(log2p).v(Range(0,npnt-1),Range::all()) = fadd*gbl->res0.v(Range(0,npnt-1),Range::all()) -gbl->res.v(Range(0,npnt-1),Range::all());
			if (basis::tri(log2p)->sm()) dres(log2p).s(Range(0,nseg-1),Range(0,basis::tri(log2p)->sm()-1),Range::all()) = fadd*gbl->res0.s(Range(0,nseg-1),Range(0,basis::tri(log2p)->sm()-1),Range::all()) -gbl->res.s(Range(0,nseg-1),Range(0,basis::tri(log2p)->sm()-1),Range::all());      
			if (basis::tri(log2p)->im()) dres(log2p).i(Range(0,ntri-1),Range(0,basis::tri(log2p)->im()-1),Range::all()) = fadd*gbl->res0.i(Range(0,ntri-1),Range(0,basis::tri(log2p)->im()-1),Range::all()) -gbl->res.i(Range(0,ntri-1),Range(0,basis::tri(log2p)->im()-1),Range::all());
			isfrst = false;
		}
		gbl->res.v(Range(0,npnt-1),Range::all()) += dres(log2p).v(Range(0,npnt-1),Range::all()); 
		if (basis::tri(log2p)->sm()) gbl->res.s(Range(0,nseg-1),Range(0,basis::tri(log2p)->sm()-1),Range::all()) += dres(log2p).s(Range(0,nseg-1),Range(0,basis::tri(log2p)->sm()-1),Range::all());
		if (basis::tri(log2p)->im()) gbl->res.i(Range(0,ntri-1),Range(0,basis::tri(log2p)->im()-1),Range::all()) += dres(log2p).i(Range(0,ntri-1),Range(0,basis::tri(log2p)->im()-1),Range::all());  
	}
	else {
		if (stage == gbl->nstage) {
			/* HACK FOR AUXILIARY FLUXES */
			for (int i=0;i<nebd;++i)
				hp_ebdry(i)->output(*gbl->log, tri_hp::auxiliary);
		}
	}
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
	tri_hp::setup_preconditioner();

	return;
}
#endif
