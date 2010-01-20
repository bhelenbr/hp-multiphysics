#include "tri_hp.h"
#include "hp_boundary.h"

//#define DEBUG
#define DEBUG_TOL 1.0e-9
// #define OP_COUNT

#ifdef OP_COUNT
#include <CHUD/CHUD.h>
#endif

void tri_hp::rsdl(int stage) {    
	/* ONLY NEED TO CALL FOR MOVEMENT BETWEEN MESHES INHERIT FROM THIS FOR SPECIFIC PHYSICS */
	if (mmovement == coupled_deformable && stage == gbl->nstage && log2p == 0) r_tri_mesh::rsdl(); 


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
	int i, n;
	if (coarse_flag) {
		for(i=0;i<npnt;++i) {
			printf("rsdl v: %d ",i);
			for (n=0;n<NV;++n) 
				printf("%e ",gbl->res.v(i,n));
			printf("\n");
		}
		
		for(i=0;i<nseg;++i) {
			for(int m=0;m<basis::tri(log2p)->sm();++m) {
				printf("rsdl s: %d %d ",i,m); 
				for(n=0;n<NV;++n)
					printf("%e ",gbl->res.s(i,m,n));
				printf("\n");
			}
		}
		
		for(i=0;i<ntri;++i) {
			for(int m=0;m<basis::tri(log2p)->im();++m) {
				printf("rsdl i: %d %d ",i,m);
				for(n=0;n<NV;++n) 
					printf("%e %e %e\n",gbl->res.i(i,m,n));
				printf("\n");
			}
		}
	}
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
	
	return;
}

void tri_hp::element_jacobian(int tind, Array<FLT,2> &K) {
	Array<TinyVector<FLT,MXTM>,1> R(NV),Rbar(NV),lf_re(NV),lf_im(NV);
	FLT dw = 1.0e-4;  //dw=sqrt(eps/l2_norm(q))
	
	ugtouht(tind);
	
	element_rsdl(tind,0,uht,lf_re,lf_im);
	for(int i=0;i<basis::tri(log2p)->tm();++i)
		for(int n=0;n<NV;++n)
			Rbar(n)(i)=lf_re(n)(i)+lf_im(n)(i);
	
	
	if (mmovement != coupled_deformable) {
		int kcol = 0;
		for(int mode = 0; mode < basis::tri(log2p)->tm(); ++mode){
			for(int var = 0; var < NV; ++var){
				uht(var)(mode) += dw;
				
				element_rsdl(tind,0,uht,lf_re,lf_im);

				int krow = 0;
				for(int i=0;i<basis::tri(log2p)->tm();++i)
					for(int n=0;n<NV;++n)
						K(krow++,kcol) = (lf_re(n)(i) +lf_im(n)(i) -Rbar(n)(i))/dw;
				
				++kcol;
				uht(var)(mode) -= dw;
			}
		}
	} else {
		Array<FLT,2> r_K(3*ND,3*ND);
		/* Get deformable mesh Jacobian */
		r_tri_mesh::element_jacobian(tind,r_K);
		const TinyVector<int,ND*3> rows(NV,NV+1,2*NV+ND,2*NV+ND+1,3*NV+2*ND,3*NV+2*ND+1);
		for (int i=0;i<3*ND;++i)
			for (int j=0;j<3*ND;++j)
				K(rows(i),rows(j)) = r_K(i,j);
		
		int kcol = 0;
		for(int mode = 0; mode < 3; ++mode){
			for(int var = 0; var < NV; ++var){
				uht(var)(mode) += dw;
				
				element_rsdl(tind,0,uht,lf_re,lf_im);

				int krow = 0;
				for(int i=0;i<3;++i) {
					for(int n=0;n<NV;++n) {
						K(krow++,kcol) = (lf_re(n)(i) +lf_im(n)(i) -Rbar(n)(i))/dw;
					}
					krow += ND;
				}
						
				for(int i=3;i<basis::tri(log2p)->tm();++i)
					for(int n=0;n<NV;++n) 
						K(krow++,kcol) = (lf_re(n)(i) +lf_im(n)(i) -Rbar(n)(i))/dw;	
				
				++kcol;
				uht(var)(mode) -= dw;
			}
			kcol += ND;
		}
		
		for(int mode = 3; mode <  basis::tri(log2p)->tm(); ++mode){
			for(int var = 0; var < NV; ++var){
				uht(var)(mode) += dw;
				
				element_rsdl(tind,0,uht,lf_re,lf_im);
				
				int krow = 0;
				for(int i=0;i<3;++i) {
					for(int n=0;n<NV;++n) {
						K(krow++,kcol) = (lf_re(n)(i) +lf_im(n)(i) -Rbar(n)(i))/dw;
					}
					krow += ND;
				}
				
				for(int i=3;i<basis::tri(log2p)->tm();++i)
					for(int n=0;n<NV;++n) 
						K(krow++,kcol) = (lf_re(n)(i) +lf_im(n)(i) -Rbar(n)(i))/dw;	
				
				++kcol;
				uht(var)(mode) -= dw;
			}
		}
		
		for(int p=0;p<3;++p) {
			for(int n=0;n<ND;++n) {
				pnts(tri(tind).pnt(p))(n) += dw;
				
				element_rsdl(tind,0,uht,lf_re,lf_im);
				
				int krow = 0;
				for(int i=0;i<3;++i) {
					for(int n=0;n<NV;++n) {
						K(krow++,kcol) = (lf_re(n)(i) +lf_im(n)(i) -Rbar(n)(i))/dw;
					}
					krow += ND;
				}
				
				for(int i=3;i<basis::tri(log2p)->tm();++i)
					for(int n=0;n<NV;++n) 
						K(krow++,kcol) = (lf_re(n)(i) +lf_im(n)(i) -Rbar(n)(i))/dw;	
				
				++kcol;
				pnts(tri(tind).pnt(p))(n) -= dw;
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
		exit(1);
#endif

		minvrt();


#ifdef DEBUG   
		// if (coarse_level || log2p != log2pmax) {
		printf("%s nstage: %d npnt: %d log2p: %d\n",gbl->idprefix.c_str(),stage,npnt,log2p);

		for(i=0;i<npnt;++i) {
			printf("%s nstage: %d ",gbl->idprefix.c_str(),i);
			for(n=0;n<NV;++n) {
				if (fabs(gbl->vprcn(i,n)) > DEBUG_TOL) printf("%8.5e ",gbl->vprcn(i,n));
				else printf("%8.5e ",0.0);
			}
			printf("\n");
		}

		for(i=0;i<npnt;++i) {
			printf("%s v: %d ",gbl->idprefix.c_str(),i);
			for(n=0;n<NV;++n) {
				if (fabs(gbl->res.v(i,n)) > DEBUG_TOL) printf("%8.5e ",gbl->res.v(i,n));
				else printf("%8.5e ",0.0);
			}
			printf("\n");
		}

		for(i=0;i<nseg;++i) {
			for(m=0;m<basis::tri(log2p)->sm();++m) {
				printf("%s s: %d ",gbl->idprefix.c_str(),i);
				for(n=0;n<NV;++n) {
					if (fabs(gbl->res.s(i,m,n)) > DEBUG_TOL) printf("%8.5e ",gbl->res.s(i,m,n));
					else printf("%8.5e ",0.0);
				}
				printf("\n");
			}
		}


		for(i=0;i<ntri;++i) {
			for(m=0;m<basis::tri(log2p)->im();++m) {
				printf("%s i: %d ",gbl->idprefix.c_str(),i);
				for(n=0;n<NV;++n) {
					if (fabs(gbl->res.i(i,m,n)) > DEBUG_TOL) printf("%8.5e ",gbl->res.i(i,m,n));
					else printf("%8.5e ",0.0);
				}
				printf("\n");
			}
		}

		for(i=0;i<npnt;++i) {
			printf("%s ug.v: %d ",gbl->idprefix.c_str(),i);
			for(n=0;n<NV;++n) {
				if (fabs(ug.v(i,n)) > DEBUG_TOL) printf("%8.5e ",ug.v(i,n));
				else printf("%8.5e ",0.0);
			}
			printf("\n");
		}

		for(i=0;i<nseg;++i) {
			for(m=0;m<basis::tri(log2p)->sm();++m) {
				printf("%s ug.s: %d ",gbl->idprefix.c_str(),i);
				for(n=0;n<NV;++n) {
					if (fabs(ug.s(i,m,n)) > DEBUG_TOL) printf("%8.5e ",ug.s(i,m,n));
					else printf("%8.5e ",0.0);
				}
				printf("\n");
			}
		}


		for(i=0;i<ntri;++i) {
			for(m=0;m<basis::tri(log2p)->im();++m) {
				printf("%s ug.i: %d ",gbl->idprefix.c_str(),i);
				for(n=0;n<NV;++n) {
					if (fabs(ug.i(i,m,n)) > DEBUG_TOL) printf("%8.5e ",ug.i(i,m,n));
					else printf("%8.5e ",0.0);
				}
				printf("\n");
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
#ifdef PTH
		pth_exit(NULL);
#endif
#ifdef MPI
		MPI_Finalize();
#endif

		exit(1);
		
 //       }
#endif

	}
}
