#include "tet_hp.h"
#include "hp_boundary.h"

// #define DEBUG
#define DEBUG_TOL 1.0e-9

// #define MGRID_TEST

void tet_hp::rsdl(int stage) {  
	/* ONLY NEED TO CALL FOR MOVEMENT BETWEEN MESHES INHERIT FROM THIS FOR SPECIFIC PHYSICS */
//    if (mmovement == coupled_deformable && stage == gbl->nstage && log2p == 0) r_tet_mesh::rsdl(); 
	
	FLT oneminusbeta = 1.0-gbl->beta(stage);
	hp_gbl->res.v(Range(0,npnt-1),Range::all()) = 0.0;
	hp_gbl->res_r.v(Range(0,npnt-1),Range::all()) *= oneminusbeta;

	if (basis::tet(log2p).em) {
		hp_gbl->res.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all()) = 0.0;
		hp_gbl->res_r.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all()) *= oneminusbeta;
		
		if (basis::tet(log2p).fm) {
			hp_gbl->res.f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all()) = 0.0;
			hp_gbl->res_r.f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all()) *= oneminusbeta;
			
			if (basis::tet(log2p).im) {
				hp_gbl->res.i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all()) = 0.0;
				hp_gbl->res_r.i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all()) *= oneminusbeta;
			}
		}
	}
	
	for(int i=0;i<nfbd;++i)
		hp_fbdry(i)->rsdl(stage);

	for(int i=0;i<nebd;++i)
		hp_ebdry(i)->rsdl(stage);

	for(int i=0;i<nvbd;++i)
		hp_vbdry(i)->rsdl(stage);
		
	helper->rsdl(stage);
			
	Array<TinyVector<FLT,MXTM>,1> lf_re(NV),lf_im(NV); 
	
	for(int tind = 0; tind<ntet;++tind) {
		
		/* LOAD SOLUTION COEFFICIENTS FOR THIS ELEMENT */
		ugtouht(tind);
		
		/* call rsdl for element */
		element_rsdl(tind,stage,uht,lf_re,lf_im);
		
		/* load imaginary local residual into global residual */
		for (int m = 0; m < basis::tet(log2p).tm; ++m) 
			for (int n = 0; n < NV; ++n)
				lf(n)(m) = lf_im(n)(m);
		
		lftog(tind,hp_gbl->res);
		
		/* load real local residual into global residual */
		for (int m = 0; m < basis::tet(log2p).tm; ++m) 
			for (int n = 0; n < NV; ++n)
				lf(n)(m) = lf_re(n)(m);
		
		lftog(tind,hp_gbl->res_r);
	}
	
	
	/* ADD IN VISCOUS/DISSIPATIVE FLUX */
	hp_gbl->res.v(Range(0,npnt-1),Range::all()) += hp_gbl->res_r.v(Range(0,npnt-1),Range::all());
	if (basis::tet(log2p).em > 0) {
		hp_gbl->res.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all()) += hp_gbl->res_r.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all());          
		if (basis::tet(log2p).fm > 0) {
			hp_gbl->res.f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all()) += hp_gbl->res_r.f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all());      
			if (basis::tet(log2p).im > 0) {
				hp_gbl->res.i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all()) += hp_gbl->res_r.i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all());      
			}
		}
	}
	

#ifdef DEBUG
	if (coarse_flag) {
		for(int i=0;i<npnt;++i) {
			printf("rsdl v: %d ",i);
			for (int n=0;n<NV;++n) 
				printf("%e ",hp_gbl->res.v(i,n));
			printf("\n");
		}
		
		for(int i=0;i<nseg;++i) {
			for(int m=0;m<basis::tet(log2p).em;++m) {
				printf("rsdl s: %d %d ",i,m); 
				for(int n=0;n<NV;++n)
					printf("%e ",hp_gbl->res.e(i,m,n));
				printf("\n");
			}
		}
		
		for(int i=0;i<ntri;++i) {
			for(int m=0;m<basis::tet(log2p).fm;++m) {
				printf("rsdl i: %d %d ",i,m);
				for(int n=0;n<NV;++n) 
					printf("%e \n",hp_gbl->res.f(i,m,n));
				printf("\n");
			}
		}
	}
#endif
	
	/*********************************************/
	/* MODIFY RESIDUALS ON COARSER MESHES            */
	/*********************************************/    
	if (coarse_flag) {
		/* CALCULATE DRIVING TERM ON FIRST ENTRY TO COARSE MESH */
		if(isfrst) {
			dres(log2p).v(Range(0,npnt-1),Range::all()) = fadd*hp_gbl->res0.v(Range(0,npnt-1),Range::all()) -hp_gbl->res.v(Range(0,npnt-1),Range::all());
			if (basis::tet(log2p).em) dres(log2p).e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all()) = fadd*hp_gbl->res0.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all()) -hp_gbl->res.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all());      
			if (basis::tet(log2p).fm) dres(log2p).f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all()) = fadd*hp_gbl->res0.f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all()) -hp_gbl->res.f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all());
			if (basis::tet(log2p).im) dres(log2p).i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all()) = fadd*hp_gbl->res0.i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all()) -hp_gbl->res.i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all());
			isfrst = false;
		}
		hp_gbl->res.v(Range(0,npnt-1),Range::all()) += dres(log2p).v(Range(0,npnt-1),Range::all()); 
		if (basis::tet(log2p).em) hp_gbl->res.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all()) += dres(log2p).e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all());
		if (basis::tet(log2p).fm) hp_gbl->res.f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all()) += dres(log2p).f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all());  
		if (basis::tet(log2p).im) hp_gbl->res.i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all()) += dres(log2p).i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all());  
	}

	
	return;
}


void tet_hp::element_jacobian(int tind, Array<FLT,2> &K) {
	Array<TinyVector<FLT,MXTM>,1> R(NV),Rbar(NV),lf_re(NV),lf_im(NV);
	Array<FLT,1> dw(NV);
#ifdef BZ_DEBUG
	const FLT eps_r = 0.0e-6, eps_a = 1.0e-6;  /*<< constants for debugging jacobians */
#else
	const FLT eps_r = 1.0e-6, eps_a = 1.0e-10;  /*<< constants for accurate numerical determination of jacobians */
#endif
	
	ugtouht(tind);
	
	dw = 0.0;
	for(int i=0;i<4;++i)
		for(int n=0;n<NV;++n)
			dw = dw + fabs(uht(n)(i));
	
	dw = dw*eps_r;
	dw = dw +eps_a;
	
	element_rsdl(tind,0,uht,lf_re,lf_im);
	for(int i=0;i<basis::tet(log2p).tm;++i) 
		for(int n=0;n<NV;++n) 
			Rbar(n)(i)=lf_re(n)(i)+lf_im(n)(i);
	
	
	int kcol = 0;
	for(int mode = 0; mode < basis::tet(log2p).tm; ++mode){
		for(int var = 0; var < NV; ++var){
			uht(var)(mode) += dw(var);
			
			element_rsdl(tind,0,uht,lf_re,lf_im);
			
			int krow = 0;
			for(int i=0;i<basis::tet(log2p).tm;++i)
				for(int n=0;n<NV;++n)
					K(krow++,kcol) = (lf_re(n)(i) +lf_im(n)(i) -Rbar(n)(i))/dw(var);
			
			++kcol;
			uht(var)(mode) -= dw(var);
		}
	}
	
	return;
}



void tet_hp::update() {
	int i,j,m,k,n,indx,indx1;
	FLT cflalpha;
	
	
	
	//l2error(hp_gbl->ibc);
	
//#ifdef petsc
//	petsc_solve();
//#else
//	superlu();
//	//superilu();
//	return;
//#endif
	
	//conjugate_gradient();
	//bicgstab();
	//cgs();
	//return;

	
	/* COUPLED MESH MOVMEMENT */
//    if (mmovement == coupled_deformable  && log2p == 0) {
//        r_tet_mesh::update();
//    }

	/* STORE INITIAL VALUES FOR NSTAGE EXPLICIT SCHEME */
	hp_gbl->ug0.v(Range(0,npnt-1),Range::all()) = ug.v(Range(0,npnt-1),Range::all());
	if (basis::tet(log2p).em) {
		hp_gbl->ug0.e(Range(0,nseg-1),Range(0,em0-1),Range::all()) = ug.e(Range(0,nseg-1),Range::all(),Range::all());
		if (basis::tet(log2p).fm) {
			hp_gbl->ug0.f(Range(0,ntri-1),Range(0,fm0-1),Range::all()) = ug.f(Range(0,ntri-1),Range::all(),Range::all());
			if (basis::tet(log2p).im) {
				hp_gbl->ug0.i(Range(0,ntet-1),Range(0,im0-1),Range::all()) = ug.i(Range(0,ntet-1),Range::all(),Range::all());
			}
		}
	}

	for(i=0;i<nfbd;++i)
		hp_fbdry(i)->update(-1);
		
	for(i=0;i<nebd;++i)
		hp_ebdry(i)->update(-1);
		
	for(i=0;i<nvbd;++i)
		hp_vbdry(i)->update(-1);

	helper->update(-1);

	for (int stage = 0; stage < gbl->nstage; ++stage) {

		rsdl(stage);	
			
		for(i=0;i<nfbd;++i)
			hp_fbdry(i)->vdirichlet();
		for(i=0;i<nfbd;++i)
			hp_fbdry(i)->edirichlet();
		for(i=0;i<nfbd;++i)
			hp_fbdry(i)->fdirichlet();

#ifdef JACOBI
		jacobi_relaxation();
#else
		minvrt();
#endif
		
		
#ifdef DEBUG
		// if (coarse_level || log2p != log2pmax) {
		for(i=0;i<npnt;++i) {
			*gbl->log << gbl->idprefix << " v: " << i << ' ';
			for(n=0;n<NV;++n) {
				if (fabs(hp_gbl->res.v(i,n)) > DEBUG_TOL) *gbl->log << hp_gbl->res.v(i,n) << ' ';
				else *gbl->log << "0.0 ";
			}
			*gbl->log << '\n';
		}
		
		for(i=0;i<nseg;++i) {
			for(m=0;m<basis::tet(log2p).em;++m) {
				*gbl->log << gbl->idprefix << " s: " << i << ' ';
				for(n=0;n<NV;++n) {
					if (fabs(hp_gbl->res.e(i,m,n)) > DEBUG_TOL) *gbl->log << hp_gbl->res.e(i,m,n) << ' ';
					else *gbl->log << "0.0 ";
				}
				*gbl->log << '\n';
			}
		}
		
		
		for(i=0;i<ntri;++i) {
			for(m=0;m<basis::tet(log2p).fm;++m) {
				*gbl->log << gbl->idprefix << " i: " << i << ' ';
				for(n=0;n<NV;++n) {
					if (fabs(hp_gbl->res.f(i,m,n)) > DEBUG_TOL) *gbl->log << hp_gbl->res.f(i,m,n) << ' ';
					else *gbl->log << "0.0 ";
				}
				*gbl->log << '\n';
			}
		}
		
		for(i=0;i<ntet;++i) {
			for(m=0;m<basis::tet(log2p).im;++m) {
				*gbl->log << gbl->idprefix << " i: " << i << ' ';
				for(n=0;n<NV;++n) {
					if (fabs(hp_gbl->res.i(i,m,n)) > DEBUG_TOL) *gbl->log << hp_gbl->res.i(i,m,n) << ' ';
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
			for(m=0;m<basis::tet(log2p).em;++m) {
				*gbl->log << gbl->idprefix << " ug.e: " << i << ' ';
				for(n=0;n<NV;++n) {
					if (fabs(ug.e(i,m,n)) > DEBUG_TOL) *gbl->log << ug.e(i,m,n) << ' ';
					else *gbl->log << "0.0 ";
				}
				*gbl->log << '\n';
			}
		}
		
		
		for(i=0;i<ntri;++i) {
			for(m=0;m<basis::tet(log2p).fm;++m) {
				*gbl->log << gbl->idprefix << " ug.f: " << i << ' ';
				for(n=0;n<NV;++n) {
					if (fabs(ug.f(i,m,n)) > DEBUG_TOL) *gbl->log << ug.f(i,m,n) << ' ';
					else *gbl->log << "0.0 ";
				}
				*gbl->log << '\n';
			}
		}
		
		for(i=0;i<ntet;++i) {
			for(m=0;m<basis::tet(log2p).im;++m) {
				*gbl->log << gbl->idprefix << " ug.i: " << i << ' ';
				for(n=0;n<NV;++n) {
					if (fabs(ug.i(i,m,n)) > DEBUG_TOL) *gbl->log << ug.i(i,m,n) << ' ';
					else *gbl->log << "0.0 ";
				}
				*gbl->log << '\n';
			}
		}
		sim::finalize(__LINE__, __FILE__, gbl->log);

		//  }
#endif

		
#ifdef MGRID_TEST
		if (!coarse_flag) continue; // TO TEST COARSE GRID CORRECTION ONLY
#endif

		cflalpha = gbl->alpha(stage)*hp_gbl->cfl(log2p);
		
		ug.v(Range(0,npnt-1),Range::all()) = hp_gbl->ug0.v(Range(0,npnt-1),Range::all()) -cflalpha*hp_gbl->res.v(Range(0,npnt-1),Range::all());

		if (basis::tet(log2p).em > 0) {
			
			ug.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all()) = hp_gbl->ug0.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all()) -cflalpha*hp_gbl->res.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all());

			if (basis::tet(log2p).fm > 0) {

				for(i=0;i<ntri;++i) {
					indx = 0;
					indx1 = 0;
					for(m=1;m<=basis::tet(log2p).em;++m) {
						for(k=1;k<=basis::tet(log2p).em-m;++k) {
							for(n=0;n<NV;++n) {
								ug.f(i,indx1,n) =  hp_gbl->ug0.f(i,indx1,n) -cflalpha*hp_gbl->res.f(i,indx,n);
							}
							++indx; ++indx1;
						}
						indx1 += em0 -basis::tet(log2p).em;
					}
				}
				if (basis::tet(log2p).im > 0) {

					for(i=0;i<ntet;++i) {
						indx = 0;
						indx1 = 0;
						for(m=1;m<=basis::tet(log2p).em-1;++m) {
							for(j=1;j<=basis::tet(log2p).em-m;++j) {
								for(k=1;k<=basis::tet(log2p).em-m-j;++k) {
									for(n=0;n<NV;++n) {
										ug.i(i,indx1,n) =  hp_gbl->ug0.i(i,indx1,n) -cflalpha*hp_gbl->res.i(i,indx,n);
									}
								   ++indx; ++indx1;
								}
								indx1 += em0 -basis::tet(log2p).em;
							}
						}
					}
				}
			}
		}
			
		helper->update(stage);

		for(i=0;i<nfbd;++i) {
			hp_fbdry(i)->update(stage);
		}
		
		for(i=0;i<nebd;++i) {
			hp_ebdry(i)->update(stage);
		}

		for(i=0;i<nvbd;++i) {
			hp_vbdry(i)->update(stage);
		}

#ifdef DEBUG
//   if (coarse_level) {
#ifdef PTH
//			pth_exit(NULL);
#endif
#ifdef MPI
//			MPI_Finalize();
#endif
//	   }
#endif
	}
}
