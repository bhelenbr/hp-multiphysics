#include "tri_hp.h"
#include "hp_boundary.h"

//#define RSDL_DEBUG
//#define DEBUG
#define DEBUG_TOL 1.0e-9
// #define OP_COUNT

#ifdef OP_COUNT
#include <CHUD/CHUD.h>
#endif

void tri_hp::rsdl(int stage) {    

	/* ONLY NEED TO CALL FOR MOVEMENT BETWEEN MESHES INHERIT FROM THIS FOR SPECIFIC PHYSICS */
#ifndef petsc
	if (mmovement == coupled_deformable && stage == gbl->nstage && log2p == 0) r_tri_mesh::rsdl(); 
#else
	if (mmovement == coupled_deformable) r_tri_mesh::rsdl(); 
#endif


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

//	for(int i=0;i<nvbd;++i)
//		hp_vbdry(i)->rsdl(stage);  // CONFLICT BETWEEN ebdry's calling rsdl and vbdry's having their own rsdl (melt_end_pt)

	for(int i=0;i<nebd;++i)
		hp_ebdry(i)->rsdl(stage);
	
	helper->rsdl(stage);

	Array<TinyVector<FLT,MXTM>,1> lf_re(NV),lf_im(NV); 
	
	for(int tind = 0; tind<ntri;++tind) {
		
		/* LOAD SOLUTION COEFFICIENTS FOR THIS ELEMENT */
		ugtouht(tind);
		
		/* call rsdl for element */
		element_rsdl(tind,stage,uht,lf_re,lf_im);
		
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
	int last_phase, mp_phase;
	int i,m,n;
	
	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		vc0load(mp_phase,gbl->res.v.data());
		pmsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
		last_phase = true;
		last_phase &= vc0wait_rcv(mp_phase,gbl->res.v.data());
	}
	
	/* APPLY VERTEX DIRICHLET B.C.'S */
	for(int i=0;i<nebd;++i)
		hp_ebdry(i)->vdirichlet();
	
	for(int i=0;i<nvbd;++i)
		hp_vbdry(i)->vdirichlet2d();
			
	sc0load(gbl->res.s.data(),0,basis::tri(log2p)->sm()-1,gbl->res.s.extent(secondDim));
	smsgpass(boundary::all,0,boundary::symmetric);
	sc0wait_rcv(gbl->res.s.data(),0,basis::tri(log2p)->sm()-1,gbl->res.s.extent(secondDim));
	
	/* APPLY DIRCHLET B.C.S TO MODE */
	for(int i=0;i<nebd;++i)
		for(int sm=0;sm<basis::tri(log2p)->sm();++sm)
			hp_ebdry(i)->sdirichlet(sm);
	
	// if (coarse_level || log2p != log2pmax) {	
	for(i=0;i<npnt;++i) {
		*gbl->log << gbl->idprefix << " v: " << i << ' ';
		for(n=0;n<NV;++n) {
			if (fabs(gbl->res.v(i,n)) > DEBUG_TOL) *gbl->log << gbl->res.v(i,n) << ' ';
			else *gbl->log << "0.0 ";
		}
		*gbl->log << '\n';
	}
	
	for(i=0;i<nseg;++i) {
		for(m=0;m<basis::tri(log2p)->sm();++m) {
			*gbl->log << gbl->idprefix << " s: " << i << ' ';
			for(n=0;n<NV;++n) {
				if (fabs(gbl->res.s(i,m,n)) > DEBUG_TOL) *gbl->log << gbl->res.s(i,m,n) << ' ';
				else *gbl->log << "0.0 ";
			}
			*gbl->log << '\n';
		}
	}
	
	
	for(i=0;i<ntri;++i) {
		for(m=0;m<basis::tri(log2p)->im();++m) {
			*gbl->log << gbl->idprefix << " i: " << i << ' ';
			for(n=0;n<NV;++n) {
				if (fabs(gbl->res.i(i,m,n)) > DEBUG_TOL) *gbl->log << gbl->res.i(i,m,n) << ' ';
				else *gbl->log << "0.0 ";
			}
			*gbl->log << '\n';
		}
	}
	sim::finalize(__LINE__,__FILE__,gbl->log);		

//	}
#endif
	
	/*********************************************/
	/* MODIFY RESIDUALS ON COARSER MESHES            */
	/*********************************************/    
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
	
	for(int i=0;i<nebd;++i)
		hp_ebdry(i)->rsdl_after(stage);
	
	return;
}

void tri_hp::element_jacobian(int tind, Array<FLT,2> &K) {
	Array<TinyVector<FLT,MXTM>,1> R(NV),Rbar(NV),lf_re(NV),lf_im(NV);
	Array<FLT,1> dw(NV);
#ifdef BZ_DEBUG
	const FLT eps_r = 0.0e-6, eps_a = 1.0e-6;  /*<< constants for debugging jacobians */
#else
	const FLT eps_r = 1.0e-6, eps_a = 1.0e-10;  /*<< constants for accurate numerical determination of jacobians */
#endif

	ugtouht(tind);
	
	dw = 0.0;
	for(int i=0;i<3;++i)
		for(int n=0;n<NV;++n)
			dw = dw + fabs(uht(n)(i));
	
	dw = dw*eps_r;
	dw = dw +eps_a;
	
	element_rsdl(tind,0,uht,lf_re,lf_im);
	for(int i=0;i<basis::tri(log2p)->tm();++i) 
		for(int n=0;n<NV;++n) 
			Rbar(n)(i)=lf_re(n)(i)+lf_im(n)(i);
	
	
	if (mmovement != coupled_deformable) {
		int kcol = 0;
		for(int mode = 0; mode < basis::tri(log2p)->tm(); ++mode){
			for(int var = 0; var < NV; ++var){
				uht(var)(mode) += dw(var);
				
				element_rsdl(tind,0,uht,lf_re,lf_im);

				int krow = 0;
				for(int i=0;i<basis::tri(log2p)->tm();++i)
					for(int n=0;n<NV;++n)
						K(krow++,kcol) = (lf_re(n)(i) +lf_im(n)(i) -Rbar(n)(i))/dw(var);
				
				++kcol;
				uht(var)(mode) -= dw(var);
			}
		}
	} else {
		Array<FLT,2> r_K(3*ND,3*ND);
		/* Get deformable mesh Jacobian */
		r_tri_mesh::element_jacobian(tind,r_K);
		const TinyVector<int,ND*3> rows(NV,NV+1,2*NV+ND,2*NV+ND+1,3*NV+2*ND,3*NV+2*ND+1);
		for (int i=0;i<3*ND;++i) {
			K(rows(i),Range::all()) = 0.0;
			for (int j=0;j<3*ND;++j)
				K(rows(i),rows(j)) = r_K(i,j);
		}
				
				
		FLT dx = eps_r*sqrt(area(tind)) +eps_a;
		
		int kcol = 0;
		for(int mode = 0; mode < 3; ++mode){
			for(int var = 0; var < NV; ++var){
				uht(var)(mode) += dw(var);
				
				element_rsdl(tind,0,uht,lf_re,lf_im);

				int krow = 0;
				for(int i=0;i<3;++i) {
					for(int n=0;n<NV;++n) {
						K(krow++,kcol) = (lf_re(n)(i) +lf_im(n)(i) -Rbar(n)(i))/dw(var);
					}
					krow += ND;
				}
						
				for(int i=3;i<basis::tri(log2p)->tm();++i)
					for(int n=0;n<NV;++n) 
						K(krow++,kcol) = (lf_re(n)(i) +lf_im(n)(i) -Rbar(n)(i))/dw(var);	
				
				++kcol;
				uht(var)(mode) -= dw(var);
			}
			
			for(int n=0;n<ND;++n) {
				pnts(tri(tind).pnt(mode))(n) += dx;
				
				element_rsdl(tind,0,uht,lf_re,lf_im);
				
				int krow = 0;
				for(int i=0;i<3;++i) {
					for(int n=0;n<NV;++n) {
						K(krow++,kcol) = (lf_re(n)(i) +lf_im(n)(i) -Rbar(n)(i))/dx;
					}
					krow += ND;
				}
				
				for(int i=3;i<basis::tri(log2p)->tm();++i)
					for(int n=0;n<NV;++n) 
						K(krow++,kcol) = (lf_re(n)(i) +lf_im(n)(i) -Rbar(n)(i))/dx;	
				
				++kcol;
				pnts(tri(tind).pnt(mode))(n) -= dx;
			}
		}
		
		for(int mode = 3; mode <  basis::tri(log2p)->tm(); ++mode){
			for(int var = 0; var < NV; ++var){
				uht(var)(mode) += dw(var);
				
				element_rsdl(tind,0,uht,lf_re,lf_im);
				
				int krow = 0;
				for(int i=0;i<3;++i) {
					for(int n=0;n<NV;++n) {
						K(krow++,kcol) = (lf_re(n)(i) +lf_im(n)(i) -Rbar(n)(i))/dw(var);
					}
					krow += ND;
				}
				
				for(int i=3;i<basis::tri(log2p)->tm();++i)
					for(int n=0;n<NV;++n) 
						K(krow++,kcol) = (lf_re(n)(i) +lf_im(n)(i) -Rbar(n)(i))/dw(var);	
				
				++kcol;
				uht(var)(mode) -= dw(var);
			}
		}
	}
	
	return;
}



void tri_hp::update() {
	int i,m,k,n,indx,indx1;
	FLT cflalpha;

	// temp fix need to better incorporate more solvers
#ifdef petsc
	petsc_update();
	return;
#endif
	
	/* COUPLED MESH MOVMEMENT */
	if (mmovement == coupled_deformable  && log2p == 0) {
		r_tri_mesh::update();
	}

	/* STORE INITIAL VALUES FOR NSTAGE EXPLICIT SCHEME */
	gbl->ug0.v(Range(0,npnt-1),Range::all()) = ug.v(Range(0,npnt-1),Range::all());
	if (basis::tri(log2p)->sm()) {
		gbl->ug0.s(Range(0,nseg-1),Range(0,sm0-1),Range::all()) = ug.s(Range(0,nseg-1),Range::all(),Range::all());
		if (basis::tri(log2p)->im()) {
			gbl->ug0.i(Range(0,ntri-1),Range(0,im0-1),Range::all()) = ug.i(Range(0,ntri-1),Range::all(),Range::all());
		}
	}

	for(i=0;i<nebd;++i)
		hp_ebdry(i)->update(-1);

	for(i=0;i<nvbd;++i)
		hp_vbdry(i)->update(-1);

	helper->update(-1);


	for (int stage = 0; stage < gbl->nstage; ++stage) {

#ifdef OP_COUNT
		chudInitialize(); 
		chudAcquireRemoteAccess(); 
		chudStartRemotePerfMonitor("fpucount"); 
		for (int i=0;i<100;++i)
#endif
		rsdl(stage);
		
#ifdef OP_COUNT
		chudStopRemotePerfMonitor(); 
		chudReleaseRemoteAccess();
		sim::finalize(__LINE__,__FILE__,gbl->log);
#endif

		minvrt();


#ifdef DEBUG   
		// if (coarse_level || log2p != log2pmax) {
		*gbl->log << gbl->idprefix << " nstage: " << stage << " npnt: " << npnt << " log2p: " << log2p << '\n';

		for(i=0;i<npnt;++i) {
			*gbl->log << gbl->idprefix << " nstage: " << i << ' ';
			for(n=0;n<NV;++n) {
				if (fabs(gbl->vprcn(i,n)) > DEBUG_TOL) *gbl->log << gbl->vprcn(i,n) << ' ';
				else *gbl->log << "0.0 ";
			}
			*gbl->log << '\n';
		}

		for(i=0;i<npnt;++i) {
			*gbl->log << gbl->idprefix << " v: " << i << ' ';
			for(n=0;n<NV;++n) {
				if (fabs(gbl->res.v(i,n)) > DEBUG_TOL) *gbl->log << gbl->res.v(i,n) << ' ';
				else *gbl->log << "0.0 ";
			}
			*gbl->log << '\n';
		}

		for(i=0;i<nseg;++i) {
			for(m=0;m<basis::tri(log2p)->sm();++m) {
				*gbl->log << gbl->idprefix << " s: " << i << ' ';
				for(n=0;n<NV;++n) {
					if (fabs(gbl->res.s(i,m,n)) > DEBUG_TOL) *gbl->log << gbl->res.s(i,m,n) << ' ';
					else *gbl->log << "0.0 ";
				}
				*gbl->log << '\n';
			}
		}


		for(i=0;i<ntri;++i) {
			for(m=0;m<basis::tri(log2p)->im();++m) {
				*gbl->log << gbl->idprefix << " i: " << i << ' ';
				for(n=0;n<NV;++n) {
					if (fabs(gbl->res.i(i,m,n)) > DEBUG_TOL) *gbl->log << gbl->res.i(i,m,n) << ' ';
					else *gbl->log << "0.0 ";
				}
				*gbl->log << '\n';
			}
		}

		for(i=0;i<npnt;++i) {
			*gbl->log << gbl->idprefix << " ug.v: " << i << ' ';
			for(n=0;n<NV;++n) {
				if (fabs(ug.v(i,n)) > DEBUG_TOL) *gbl->log << ug.v(i,n) << ' ';
				else *gbl->log << "0.0 ";
			}
			*gbl->log << '\n';
		}

		for(i=0;i<nseg;++i) {
			for(m=0;m<basis::tri(log2p)->sm();++m) {
				*gbl->log << gbl->idprefix << " ug.s: " << i << ' ';
				for(n=0;n<NV;++n) {
					if (fabs(ug.s(i,m,n)) > DEBUG_TOL) *gbl->log << ug.s(i,m,n) << ' ';
					else *gbl->log << "0.0 ";
				}
				*gbl->log << '\n';
			}
		}


		for(i=0;i<ntri;++i) {
			for(m=0;m<basis::tri(log2p)->im();++m) {
				*gbl->log << gbl->idprefix << " ug.i: " << i << ' ';
				for(n=0;n<NV;++n) {
					if (fabs(ug.i(i,m,n)) > DEBUG_TOL) *gbl->log << ug.i(i,m,n) << ' ';
					else *gbl->log << "0.0 ";
				}
				*gbl->log << '\n';
			}
		}
	//  }
#endif

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

		helper->update(stage);

		for(i=0;i<nebd;++i) {
			hp_ebdry(i)->update(stage);
		}

		for(i=0;i<nvbd;++i) {
			hp_vbdry(i)->update(stage);
		}

#ifdef DEBUG
//        if (coarse_level || log2p != log2pmax) {
						sim::finalize(__LINE__,__FILE__,gbl->log);		
 //       }
#endif
	}
}
