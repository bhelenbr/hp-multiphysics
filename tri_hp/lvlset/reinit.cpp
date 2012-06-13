#include "tri_hp_lvlset.h"
#include "../tri_hp.h"	
#include "../hp_boundary.h"
#include <math.h>

// #define DEBUG_TOL 1.0e-9
// #define RSDL_DEBUG
#define DEBUG_OUTPUT

// #define USECOMM
// #define ADD_DIFFUSION

void tri_hp_lvlset::reinitialize() {	
	
	const int store_log2p = log2p;  
	/* Only reinitialize vertices? */
	// log2p = 0;  // Uncomment this to only reinitialize vertex description of levelset

	for(int i=0;i<nebd;++i)
		reinit_bdry(i)->init();
				
#ifdef DEBUG_OUTPUT
	std::ostringstream nstr;
	std::string fname;
	std::stringstream s;
	nstr << gbl->tstep << "." << gbl->substep << std::flush;
#endif

	for (int ii=0; ii<reinit_iterations; ii++){
		
		/* Lets set all high-order stuff away from the interface to 0 */
		for(int i=0;i<nseg;++i) {
			if (fabs(ug.v(seg(i).pnt(0),2) +ug.v(seg(i).pnt(1),2))/(gbl->width*basis::tri(log2p)->sm()) > 1.0) {
				for(int m=0;m<basis::tri(log2p)->sm();++m) {
					ug.s(i,m,2) = 0.0;
				}
			}
		}
		
		for(int i=0;i<ntri;++i) {
			if (fabs(ug.v(tri(i).pnt(0),2) +ug.v(tri(i).pnt(1),2) +ug.v(tri(i).pnt(2),2))/(gbl->width*basis::tri(log2p)->sm()) > 1.0) {
				for(int m=0;m<basis::tri(log2p)->im();++m) {
					ug.i(i,m,2) = 0.0;
				}
			}
		}
		
		reinit_setup_preconditioner();
		reinit_update();
		
#ifdef DEBUG_OUTPUT
		/* To output during reinitialization */
		s.str("");
		s<<ii;
		fname = "reinit" + nstr.str() + "_" + s.str() +"_" +gbl->idprefix;
		output(fname,tecplot);
#endif
	}
	
	log2p = store_log2p;
}

void tri_hp_lvlset::reinit_update() {
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
		reinit_rsdl(stage);

		/* SET RESDIUALS TO ZERO ON INCOMING CHARACTERISTIC BOUNDARIES */
		reinit_minvrt();
		
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

void tri_hp_lvlset::reinit_rsdl(int stage) {
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
	
	for(int tind = 0; tind<ntri;++tind) {
		
		/* LOAD SOLUTION COEFFICIENTS FOR THIS ELEMENT */
		ugtouht(tind);
		
		/* call rsdl for element */
		reinit_element_rsdl(tind,stage,uht,lf_re,lf_im);
		
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
#endif
}

void tri_hp_lvlset::reinit_element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uht,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im) {
	int i,j,n;
	const int NV = 4;
	TinyVector<int,3> v;
	TinyMatrix<FLT,ND,ND> ldcrd;
	TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,NV,ND> du;
	int lgpx = basis::tri(log2p)->gpx(), lgpn = basis::tri(log2p)->gpn();
	FLT  cjcb, oneminusbeta;
	TinyMatrix<FLT,ND,ND> visc;
	TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,NV-1,NV-1> cv, df;
	TinyVector<FLT,NV> tres;
	TinyMatrix<FLT,MXGP,MXGP> rho, mu,e00,e01;
	FLT length, phidw, signphi;
	TinyVector<FLT,ND> tang,norm;
	TinyVector<TinyMatrix<FLT,MXGP,MXGP>,ND> phivel;

	oneminusbeta = 1.0-gbl->beta(stage);
	// LOAD INDICES OF VERTEX POINTS //
	v = tri(tind).pnt;


	// PROJECT VERTEX COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS //
	for(n=0;n<ND;++n)
		basis::tri(log2p)->proj(pnts(v(0))(n),pnts(v(1))(n),pnts(v(2))(n),&crd(n)(0,0),MXGP);

	// CALCULATE COORDINATE DERIVATIVES A SIMPLE WAY //
	for(n=0;n<ND;++n) {
		ldcrd(n,0) = 0.5*(pnts(v(2))(n) -pnts(v(1))(n));
		ldcrd(n,1) = 0.5*(pnts(v(0))(n) -pnts(v(1))(n));
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
	
	cjcb = ldcrd(0,0)*ldcrd(1,1) -ldcrd(1,0)*ldcrd(0,1);

	for(i=0;i<lgpx;++i) {
		for(j=0;j<lgpn;++j) {
			// STUFF FOR LEVEL SET //
			phidw = u(2)(i,j)/(gbl->width);
			norm(0) = (+ldcrd(1,1)*du(2,0)(i,j) -ldcrd(1,0)*du(2,1)(i,j))/cjcb;
			norm(1) = (-ldcrd(0,1)*du(2,0)(i,j) +ldcrd(0,0)*du(2,1)(i,j))/cjcb;
			length = sqrt(norm(0)*norm(0) +norm(1)*norm(1));

			if (phidw < -1.0) {
				res(2)(i,j) = -RAD(crd(0)(i,j))*cjcb*(length-1.0);
				phivel(0)(i,j) = -norm(0)/length;
				phivel(1)(i,j) = -norm(1)/length;
				signphi = -1.0;
			}
			else if (phidw > 1.0) {
				res(2)(i,j) = RAD(crd(0)(i,j))*cjcb*(length-1.0);
				phivel(0)(i,j) = norm(0)/length;
				phivel(1)(i,j) = norm(1)/length;
				signphi = 1.0;
			}
			else {
				// signphi = phidw;
				// signphi = u(2)(i,j) / pow( pow(u(2)(i,j),2) + pow(1-length,2) * pow(gbl->width/2.0,2), 0.5 );
				// signphi = u(2)(i,j) / pow( pow(u(2)(i,j),2) + pow(1-phidw,2) * pow(gbl->width/2.0,2), 0.5 );
				// signphi = u(2)(i,j) / pow( pow(u(2)(i,j),2) + pow(length,2) * pow(gbl->width,2), 0.5 );
				signphi = sin(M_PI*phidw/2.0);
				phivel(0)(i,j) = signphi/(fabs(signphi)+EPSILON)*norm(0)/length;
				phivel(1)(i,j) = signphi/(fabs(signphi)+EPSILON)*norm(1)/length;
				res(2)(i,j) = RAD(crd(0)(i,j))*cjcb*(length-1.0)*signphi;
			}			
		}
	}
	basis::tri(log2p)->intgrt(&lf_im(NV-2)(0),&res(NV-2)(0,0),MXGP);

	// NEGATIVE REAL TERMS //
	if (gbl->beta(stage) > 0.0) {
	
#ifdef ADD_DIFFUSION
		/* DIFFUSION TENSOR (LOTS OF SYMMETRY THOUGH)*/
		/* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
		FLT mu = gbl->width/cjcb;
		visc(0,0) = -mu*(ldcrd(1,1)*ldcrd(1,1) +ldcrd(0,1)*ldcrd(0,1));
		visc(1,1) = -mu*(ldcrd(1,0)*ldcrd(1,0) +ldcrd(0,0)*ldcrd(0,0));
		visc(0,1) =  mu*(ldcrd(1,1)*ldcrd(1,0) +ldcrd(0,1)*ldcrd(0,0));		
		

		/* DIFFUSIVE TERMS?  */
		for(i=0;i<basis::tri(log2p)->gpx();++i) {
			for(j=0;j<basis::tri(log2p)->gpn();++j) {
				
				phidw = u(2)(i,j)/(gbl->width);

				if (phidw < -1.0) {
					signphi = -1.0;
				}
				else if (phidw > 1.0) {
					signphi = 1.0;
				}
				else {
					signphi = sin(M_PI*phidw/2.0);
				}
				
				e00(i,j) = fabs(signphi)*(visc(0,0)*du(2,0)(i,j) +visc(0,1)*du(2,1)(i,j));
				e01(i,j) = fabs(signphi)*(viscI1II0I*du(2,0)(i,j) +visc(1,1)*du(2,1)(i,j));
			}
		}
		basis::tri(log2p)->derivr(e00.data(),&res(2)(0,0),MXGP);
		basis::tri(log2p)->derivs(e01.data(),&res(2)(0,0),MXGP);
#else
		e00 = 0.0;
		e01 = 0.0;
#endif

		// THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES //
		for(i=0;i<lgpx;++i) {
			for(j=0;j<lgpn;++j) {
				tres(2) = gbl->tau(tind,2)*res(2)(i,j);

				e00(i,j) -= (ldcrd(1,1)*phivel(0)(i,j) -ldcrd(0,1)*phivel(1)(i,j))*tres(2);
				e01(i,j) -= (-ldcrd(1,0)*phivel(0)(i,j) +ldcrd(0,0)*phivel(1)(i,j))*tres(2);
			}
		}
		basis::tri(log2p)->intgrtrs(&lf_re(2)(0),e00.data(),e01.data(),MXGP);
		
		for(i=0;i<basis::tri(log2p)->tm();++i)
			lf_re(2)(i) *= gbl->beta(stage);

	}

	return;
}

void tri_hp_lvlset::reinit_setup_preconditioner() {

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
				tri_mesh::output("negative_"+gbl->idprefix+gbl->idprefix,grid);
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
	
#ifdef USECOMM
	tri_hp::setup_preconditioner();
#else
	/* PREINVERT PRECONDITIONER FOR VERTICES */
	gbl->vprcn(Range(0,npnt-1),Range::all()) = 1.0/(basis::tri(log2p)->vdiag()*gbl->vprcn(Range(0,npnt-1),Range::all()));

	if (log2p) {
		/* INVERT DIAGANOL PRECONDITIONER FOR SIDES */ 
		gbl->sprcn(Range(0,nseg-1),Range::all()) = 1.0/gbl->sprcn(Range(0,nseg-1),Range::all());
	}
#endif

	return;
}

void tri_hp_lvlset::reinit_minvrt() {
	int i,j(0),k,m,tind,sind,v0,indx,indx1,indx2,sgn,msgn;
	TinyVector<int,3> sign,side;
	Array<FLT,2> tinv(NV,NV);
	Array<FLT,1> temp(NV);
#ifdef USECOMM
	int last_phase,mp_phase;
#endif


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
	
#ifdef USECOMM
	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		vc0load(mp_phase,gbl->res.v.data());
		pmsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
		last_phase = true;
		last_phase &= vc0wait_rcv(mp_phase,gbl->res.v.data());
	}
#endif

	/* APPLY VERTEX DIRICHLET B.C.'S */
	for(i=0;i<nebd;++i) 
		reinit_bdry(i)->vdirichlet();
		
	/* APPLY VERTEX DIRICHLET B.C.'S */
	for(i=0;i<nvbd;++i) 
		gbl->res.v(vbdry(i)->pnt,2) = 0.0;
		
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
		
#ifdef USECOMM
		sc0load(gbl->res.s.data(),mode,mode,gbl->res.s.extent(secondDim));
		smsgpass(boundary::all,0,boundary::symmetric);
		sc0wait_rcv(gbl->res.s.data(),mode,mode,gbl->res.s.extent(secondDim));
#endif
		
		/* APPLY SIDE DIRICHLET B.C.'S */
		for(i=0;i<nebd;++i) 
			reinit_bdry(i)->sdirichlet(mode);
			
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

#ifdef USECOMM
	sc0load(gbl->res.s.data(),mode,mode,gbl->res.s.extent(secondDim));
	smsgpass(boundary::all,0,boundary::symmetric);
	sc0wait_rcv(gbl->res.s.data(),mode,mode,gbl->res.s.extent(secondDim));
#endif
	
	/* APPLY VERTEX DIRICHLET B.C.'S */
	for(i=0;i<nebd;++i)
		reinit_bdry(i)->sdirichlet(mode);
	
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

void tri_hp_lvlset::reinit_bc::init() {
	// Set up reinitizialization B.C.'s
	// Find grad phi at point closest to zero level
	// Based on grad phi decide if side should be a dirchlet b.c.
	
	int sind, tind, v0;
	const int nseg(base.nseg);
	
	FLT phi0 = x.ug.v(x.seg(base.seg(0)).pnt(0),2);
	FLT phi1 = x.ug.v(x.seg(base.seg(nseg-1)).pnt(1),2);
	if (fabs(phi0) < fabs(phi1)) {
		phi_at_endpt = phi0;
		closer_endpt = 0;
		sign_phi = (phi1 > 0.0 ? 1.0 : -1.0);
		sind = base.seg(0);
		tind = x.seg(sind).tri(0); 
		v0 = x.seg(sind).pnt(0);
	}
	else {
		phi_at_endpt = phi1;
		sign_phi = (phi0 > 0.0 ? 1.0 : -1.0);
		closer_endpt = 1;
		sind = base.seg(nseg-1);
		tind = x.seg(sind).tri(0); 
		v0 = x.seg(sind).pnt(1);
	}
	
	// At what point of triangle do we need to know derivative? //
	int j;
	for (j=0;j<3;++j)
	if (x.tri(tind).pnt(j) == v0) break;
	assert(j < 3);
	
	// Figure out what r & s should be based on vertex # //
	FLT r,s;
	r=-1.0;
	s=-1.0;
	if (j==0) s=1.0;
	if (j==2) r=1.0;	
	
	/* Evaluate solution derivatives at point */
	TinyVector<FLT,tri_mesh::ND> xpt,dxdr,dxds,dxmax;
	TinyVector<FLT,4> upt,dudr,duds;
	x.ugtouht(tind);
	x.crdtocht(tind);
	basis::tri(x.log2p)->ptprobe_bdry(x.ND,xpt.data(),dxdr.data(),dxds.data(),r,s,&x.cht(0,0),MXTM); 
	basis::tri(x.log2p)->ptprobe(x.NV,upt.data(),dudr.data(),duds.data(),r,s,&x.uht(0)(0),MXTM); 
	
	/* Calculate grad phi at hybrid point */
	TinyMatrix<FLT,2,2> ldcrd;
	ldcrd(0,0) = dxdr(0);
	ldcrd(0,1) = dxds(0);
	ldcrd(1,0) = dxdr(1);
	ldcrd(1,1) = dxds(1);
	FLT cjcb = ldcrd(0,0)*ldcrd(1,1) -ldcrd(1,0)*ldcrd(0,1);
	norm(0) = (+ldcrd(1,1)*dudr(2) -ldcrd(1,0)*duds(2))/cjcb;
	norm(1) = (-ldcrd(0,1)*dudr(2) +ldcrd(0,0)*duds(2))/cjcb;
	FLT length = sqrt(norm(0)*norm(0) +norm(1)*norm(1)); //for projection in this direction
	norm(0) /= length;
	norm(1) /= length;
	
	/* Determine if normal is in or out (outward normal)*/
	int v1 = x.seg(sind).pnt(0);
	int v2 = x.seg(sind).pnt(1);
	FLT nrmx =  (x.pnts(v2)(1) -x.pnts(v1)(1));
	FLT nrmy = -(x.pnts(v2)(0) -x.pnts(v1)(0));
	length = sqrt(nrmx*nrmx +nrmy*nrmy);
	nrmx /= length;
	nrmy /= length;

	//	if ((nrmx*norm(0) +nrmy*norm(1))*sign_phi < 0.35) {
	//		active = true;
	//		*x.gbl->log << base.idprefix << ' ' << v0 << ' ' << phi_at_endpt << ' ' << norm << ' ' << sign_phi << std::endl;
	//	}
	//	else {
	//		active = false;
	//	}
		
	// This is a complete hack and is totally not general TEMPORARY
	if (v0 == 0 || v0 == 3) {
		active = true;
	}
	else {
		active = false;
	}
	
#ifdef USECOMM
	// This is for periodic b.c.'s to make non active
	active = false;
#endif

	/* Now set values if appropriate */
	set_values();
}

void tri_hp_lvlset::reinit_bc::set_values() {
	// this determines  phi along the edge
	// phi is calculated as the signed distance normal (nx,ny) to the base point, with phi value of basephi)
	if (!active) return;
	
	int v0;
	if (closer_endpt == 0) {
		v0 = x.seg(base.seg(0)).pnt(0);
	}
	else {
		v0 = x.seg(base.seg(base.nseg-1)).pnt(1);
	}
	
	/* Re-evaluate phi at inflow points */
	for(int j=0;j<base.nseg;++j) {
		int sind = base.seg(j);
		int v1 = x.seg(sind).pnt(0);
		TinyVector<FLT,tri_mesh::ND> dx;
		dx = x.pnts(v1) -x.pnts(v0);
		
		x.ug.v(v1,2) = dx(0)*norm(0) +dx(1)*norm(1) +phi_at_endpt;
		for(int m=0;m<basis::tri(x.log2p)->sm();++m)
			x.ug.s(sind,m,2) = 0.0;
	}
	int sind = base.seg(base.nseg-1);
	int v1 = x.seg(sind).pnt(1);
	TinyVector<FLT,tri_mesh::ND> dx;
	dx = x.pnts(v1) -x.pnts(v0);
	x.ug.v(v1,2) = dx(0)*norm(0) +dx(1)*norm(1) +phi_at_endpt;
	
	return;
}

void tri_hp_lvlset::reinit_bc::vdirichlet() {
	if (!active) return;
	
	int sind,j,v0;
	j = 0;
	do {
		sind = base.seg(j);
		v0 = x.seg(sind).pnt(0);
		x.gbl->res.v(v0,2) = 0.0;
	} while (++j < base.nseg);
	v0 = x.seg(sind).pnt(1);
	x.gbl->res.v(v0,2) = 0.0;
}
	
void tri_hp_lvlset::reinit_bc::sdirichlet(int mode) {
	int sind;
	
	for(int j=0;j<base.nseg;++j) {
		sind = base.seg(j);
		x.gbl->res.s(sind,mode,2) = 0.0;
	}
}
	
	
	
	

