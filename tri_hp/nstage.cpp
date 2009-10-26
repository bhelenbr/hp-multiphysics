#include "tri_hp.h"
#include "hp_boundary.h"

// #define DEBUG
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

	return;
}


void tri_hp::update() {
	int i,m,k,n,indx,indx1;
	FLT cflalpha;

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
        if (coarse_level || log2p != log2pmax) {
#ifdef PTH
		pth_exit(NULL);
#endif
#ifdef MPI
//		MPI_Finalize();
#endif
        }
#endif

	}
}
