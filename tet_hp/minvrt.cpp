#include "tet_hp.h"
#include "hp_boundary.h"
#include <myblas.h>

//#define NODAL

/************************************************/
/**********        INVERT MASS MATRIX     *******/
/************************************************/
#ifndef NODAL
void tet_hp::minvrt() {
	int i,j,k,n,m,tind,msgn,sgn,side,sind,find,v0,indx,ind,indx2,info;
	Array<FLT,2> spokemass;
	char trans[] = "T";
	int last_phase, mp_phase;
	
	/* LOOP THROUGH EDGES */
	if (basis::tet(log2p).em > 0) {
		indx = 0;
		for(int eind = 0; eind<nseg;++eind) {
			/* SUBTRACT SIDE CONTRIBUTIONS TO VERTICES */         
			for (k=0; k <basis::tet(log2p).em; ++k) {
				for (i=0; i<2; ++i) {
					v0 = seg(eind).pnt(i);
					for(n=0;n<NV;++n)
						gbl->res.v(v0,n) -= basis::tet(log2p).sfmv(i,k)*gbl->res.e(eind,k,n);
				}
			}
		}
			
		if (basis::tet(log2p).fm > 0) {
			/* SUBTRACT FACES */
			indx = 0;
			for(find = 0; find<ntri;++find) {
				indx2 = 3;
				for (i=0; i<3; ++i) {
					v0 = tri(find).pnt(i);
					for (k=0;k<basis::tet(log2p).fm;++k)
						for(n=0;n<NV;++n)
							gbl->res.v(v0,n) -= basis::tet(log2p).ffmv(i,k)*gbl->res.f(find,k,n);

				}
			}

			if (basis::tet(log2p).im > 0) {
				/* SUBTRACT INTERIOR */
				for(tind = 0; tind<ntet;++tind) {
					for (i=0; i<4; ++i) {
						v0 = tet(tind).pnt(i);
						for (k=0;k<basis::tet(log2p).im;++k)
							for(n=0;n<NV;++n)
								gbl->res.v(v0,n) -= basis::tet(log2p).ifmb(i,k)*gbl->res.i(tind,k,n);

					}
				}
			}
		}
	}
		
	if (gbl->diagonal_preconditioner) {
		gbl->res.v(Range(0,npnt-1),Range::all()) *= gbl->vprcn(Range(0,npnt-1),Range::all())*basis::tet(log2p).vdiag;
	}    
	
	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		pc0load(mp_phase,gbl->res.v.data());
		pmsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
		last_phase = true;
		last_phase &= pc0wait_rcv(mp_phase,gbl->res.v.data());
	}
		
	/* APPLY VERTEX DIRICHLET B.C.'S */
	for(i=0;i<nfbd;++i)
		hp_fbdry(i)->vdirichlet();
		
	for(i=0;i<nebd;++i)
		hp_ebdry(i)->vdirichlet3d();        
		
	for(i=0;i<nvbd;++i)
		hp_vbdry(i)->vdirichlet3d();

	if(basis::tet(log2p).em == 0) return;

	/* REMOVE VERTEX CONTRIBUTION FROM SIDE MODES */
	/* SOLVE FOR SIDE MODES */
	/* PART 1 REMOVE VERTEX CONTRIBUTIONS */
	for(tind=0;tind<ntet;++tind) {         
		if (gbl->diagonal_preconditioner) { 
			for(i=0;i<4;++i) {
				v0 = tet(tind).pnt(i);
				for(n=0;n<NV;++n)
					uht(n)(i) = gbl->res.v(v0,n)*gbl->iprcn(tind,n);
			}
			/* edges */
			for(i=0;i<6;++i) {
				sind = tet(tind).seg(i);
				sgn  = tet(tind).sgn(i);
				for(j=0;j<4;++j) {
					msgn = 1;
					for(k=0;k<basis::tet(log2p).em;++k) {
						for(n=0;n<NV;++n)
							gbl->res.e(sind,k,n) -= msgn*basis::tet(log2p).vfms(j,4+k+i*basis::tet(log2p).em)*uht(n)(j);
						msgn *= sgn;
					}
				}
			}
			/* faces */
			for(i=0;i<4;++i) {
				find = tet(tind).tri(i);
				sgn  = -tet(tind).rot(i);
				for(j=0;j<4;++j) {
					msgn = 1;
					ind = 0;
					for(m = 1; m <= basis::tet(log2p).em-1; ++m) {
						for(k = 1; k <= basis::tet(log2p).em-m; ++k) {
							for(n=0;n<NV;++n)
								gbl->res.f(find,ind,n) -= msgn*basis::tet(log2p).vfms(j,4+6*basis::tet(log2p).em+i*basis::tet(log2p).fm+ind)*uht(n)(j);
							++ind;
						}
						msgn *= sgn;
					}
				}
			}
			/* interior */
			for(j=0;j<4;++j) {
				for(k=0;k<basis::tet(log2p).im;++k) {
					for(n=0;n<NV;++n)
						gbl->res.i(tind,k,n) -= basis::tet(log2p).vfms(j,basis::tet(log2p).bm+k)*uht(n)(j);
				}
			}			
		}
	}

	/* VERTEX BALL TO SOLVE FOR EDGE MODES */
	if (basis::tet(log2p).em > 1) {
		for(sind = 0; sind< nseg;++sind) {
			for(n=0;n<NV;++n) {
				wkseg(sind,0,n)(0) = gbl->res.e(sind,0,n)+basis::tet(log2p).sfms(0)*gbl->res.e(sind,1,n);
				wkseg(sind,0,n)(1) = gbl->res.e(sind,0,n)-basis::tet(log2p).sfms(0)*gbl->res.e(sind,1,n);
				gbl->res.e(sind,0,n) = 0.0;
			}
		}
			
		for(find = 0; find<ntri;++find) {
			sind = tri(find).seg(0);
			/* i = 0 */
			for (k=0;k<basis::tet(log2p).fm;++k){
				for(n=0;n<NV;++n){
					if(seg(sind).pnt(0) == tri(find).pnt(2) ){
						wkseg(sind,0,n)(0) += basis::tet(log2p).ffms(0,k,1)*gbl->res.f(find,k,n);
						wkseg(sind,0,n)(1) += basis::tet(log2p).ffms(0,k,0)*gbl->res.f(find,k,n);
					}
					else{
						wkseg(sind,0,n)(0) += basis::tet(log2p).ffms(0,k,0)*gbl->res.f(find,k,n);
						wkseg(sind,0,n)(1) += basis::tet(log2p).ffms(0,k,1)*gbl->res.f(find,k,n);
					}						
				}					
			}
					
			for (i=1; i<3; ++i) {
				sind = tri(find).seg(i);
				for (k=0;k<basis::tet(log2p).fm;++k){
					for(n=0;n<NV;++n){
						if(seg(sind).pnt(0) == tri(find).pnt(0) ){
							wkseg(sind,0,n)(0) += basis::tet(log2p).ffms(i,k,1)*gbl->res.f(find,k,n);
							wkseg(sind,0,n)(1) += basis::tet(log2p).ffms(i,k,0)*gbl->res.f(find,k,n);
						}
						else {
							wkseg(sind,0,n)(0) += basis::tet(log2p).ffms(i,k,0)*gbl->res.f(find,k,n);
							wkseg(sind,0,n)(1) += basis::tet(log2p).ffms(i,k,1)*gbl->res.f(find,k,n);
						}						
					}					
				}
			}		
		}
		
		for(int vind = 0; vind < npnt; ++vind){
			for(n=0;n<NV;++n){
				for(int spk = 0; spk < pnt(vind).nspk; ++spk){
					spkres(spk) = wkseg(spklink(vind)(spk)(0),0,n)(spklink(vind)(spk)(1));
				}
				
				GETRS(trans,pnt(vind).nspk,1,&spkmass(vind)(0,0),pnt(vind).nspk,&spkpiv(vind)(0),&spkres(0),pnt(vind).nspk,info); 
				if (info != 0) {
					printf("DGETRS FAILED - VERTEX BALL info:%d vertex:%d \n",info,vind);
					exit(1);
				}
				for(int spk = 0; spk < pnt(vind).nspk; ++spk){
					gbl->res.e(spklink(vind)(spk)(0),0,n)+=spkres(spk)/2.0;
				}
			}
		}
	}
	
	FLT	diagcoef = 1.0;
	if(basis::tet(log2p).p == 2) {
		basis::tet(log2p).ediag(0) = diagcoef*80.0;//157
		gbl->res.e(Range(0,nseg-1),0,Range::all()) *= gbl->eprcn(Range(0,nseg-1),Range::all())*basis::tet(log2p).ediag(0);
		
		for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
			sc0load(mp_phase,gbl->res.e.data(),0,0,gbl->res.e.extent(secondDim));
			smsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
			last_phase = true;
			last_phase &= sc0wait_rcv(mp_phase,gbl->res.e.data(),0,0,gbl->res.e.extent(secondDim));
		}
	}
	
	if(basis::tet(log2p).p == 3) {
		basis::tet(log2p).ediag(1) = diagcoef*1000.0;//1890 optimize later
		basis::tet(log2p).fdiag(0) = diagcoef*3000.0;//5670
	}
	
	/* APPLY EDGE DIRICHLET B.C.'S */
    for(i=0;i<nfbd;++i)
        hp_fbdry(i)->edirichlet();	
		
	for (i=0;i<nebd;++i) 
		hp_ebdry(i)->edirichlet3d();
		
	if(basis::tet(log2p).em > 1) {
		for(tind=0;tind<ntet;++tind) {         
			if (gbl->diagonal_preconditioner) { 
				for(i=0;i<6;++i) {
					side = tet(tind).seg(i);
					for(n=0;n<NV;++n)
						uht(n)(i) = gbl->res.e(side,0,n)*gbl->iprcn(tind,n);
				}
				/* edges */
				for(i=0;i<6;++i) {
					sind = tet(tind).seg(i);
					sgn  = tet(tind).sgn(i);
					for(j=0;j<6;++j) {
						msgn = sgn;
						for(k=1;k<basis::tet(log2p).em;++k) {
							for(n=0;n<NV;++n)
								gbl->res.e(sind,k,n) -= msgn*basis::tet(log2p).mm(k+4+i*basis::tet(log2p).em,4+j*basis::tet(log2p).em)*uht(n)(j);
							msgn *= sgn;
						}
					}
				}
				/* faces */
				for(i=0;i<4;++i) {
					find = tet(tind).tri(i);
					sgn  = -tet(tind).rot(i);
					for(j=0;j<6;++j) {
						msgn = 1;
						ind = 0;
						for(m = 1; m <= basis::tet(log2p).em-1; ++m) {
							for(k = 1; k <= basis::tet(log2p).em-m; ++k) {
								for(n=0;n<NV;++n)
									gbl->res.f(find,ind,n) -= msgn*basis::tet(log2p).mm(4+6*basis::tet(log2p).em+i*basis::tet(log2p).fm+ind,4+j*basis::tet(log2p).em)*uht(n)(j);
								++ind;
							}
							msgn *= sgn;
						}
					}
				}
//				/* interior */
//				for(j=0;j<6;++j) {
//					for(k=0;k<basis::tet(log2p).im;++k) {
//						for(n=0;n<NV;++n)
//							gbl->res.i(tind,k,n) -= basis::tet(log2p).vfms(j+4,4+6*basis::tet(log2p).em+4*basis::tet(log2p).fm+k)*uht(n)(j);
//					}
//				}


				
			}
		}
		if (gbl->diagonal_preconditioner) {
			gbl->res.e(Range(0,nseg-1),1,Range::all()) *= gbl->eprcn(Range(0,nseg-1),Range::all())*basis::tet(log2p).ediag(1);
		}
		
		for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
			sc0load(mp_phase,gbl->res.e.data(),0,basis::tet(log2p).em-1,gbl->res.e.extent(secondDim));
			smsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
			last_phase = true;
			last_phase &= sc0wait_rcv(mp_phase,gbl->res.e.data(),0,basis::tet(log2p).em-1,gbl->res.e.extent(secondDim));
		}
		
		/* APPLY EDGE DIRICHLET B.C.'S */
		for(i=0;i<nfbd;++i)
			hp_fbdry(i)->edirichlet();	
		
		for (i=0;i<nebd;++i) 
			hp_ebdry(i)->edirichlet3d();
		
		for(tind=0;tind<ntet;++tind) {         
			if (gbl->diagonal_preconditioner) { 
				for(i=0;i<6;++i) {
					side = tet(tind).seg(i);
					for(n=0;n<NV;++n)
						uht(n)(i) = gbl->res.e(side,1,n)*gbl->iprcn(tind,n);
				}

				//faces
				for(i=0;i<4;++i) {
					find = tet(tind).tri(i);
					sgn  = -tet(tind).rot(i);
					for(j=0;j<6;++j) {
						msgn = 1;
						ind = 0;
						for(m = 1; m <= basis::tet(log2p).em-1; ++m) {
							for(k = 1; k <= basis::tet(log2p).em-m; ++k) {
								for(n=0;n<NV;++n)
									gbl->res.f(find,ind,n) -= msgn*basis::tet(log2p).mm(4+6*basis::tet(log2p).em+i*basis::tet(log2p).fm+ind,5+j*basis::tet(log2p).em)*uht(n)(j);
								++ind;
							}
							msgn *= sgn;
						}
					}
				}
			}
		}
	}
	
	/* ALL HIGH ORDER MODES */
	/* LOOP THROUGH EDGES */	         
	if (basis::tet(log2p).fm > 0){
		for(find = 0; find<ntri;++find) {
			for (k=0;k<basis::tet(log2p).fm;++k){
				for(n=0;n<NV;++n){
					gbl->res.f(find,k,n) *=gbl->fprcn(find,n)*basis::tet(log2p).fdiag(k); 
				}         
			}   
		}
		
		tc0load(gbl->res.f.data(),0,basis::tet(log2p).fm-1,gbl->res.f.extent(secondDim));
		tmsgpass(boundary::all,0,boundary::symmetric);
		tc0wait_rcv(gbl->res.f.data(),0,basis::tet(log2p).fm-1,gbl->res.f.extent(secondDim));
		
		for(i=0;i<nfbd;++i)
			hp_fbdry(i)->fdirichlet();

//      if (basis::tet(log2p).im > 0) {
//         for(tind = 0; tind<ntet;++tind) {
//            for(k=0;k<basis::tet(log2p).im;++k){
//               for(n=0;n<NV;++n){
//                  gbl->res.i(tind,k,n) *=gbl->iprcn(tind,n)*basis::tet(log2p).idiag(k); 
//               }
//            }
//         }
//      }
	}
		
	return;
}
#endif

#ifdef NODAL
void tet_hp::minvrt() {
	int i,j,k,n,m,tind,msgn,sgn,side,sind,find,v0,indx,ind,indx2,info;
	//cout << basis::tet(log2p).vdiag << ' ' << basis::tet(log2p).ediag << ' ' << basis::tet(log2p).fdiag << endl;
	basis::tet(log2p).vdiag = 5;
	/* VERTEX */
	if (gbl->diagonal_preconditioner) {
		gbl->res.v(Range(0,npnt-1),Range::all()) *= gbl->vprcn(Range(0,npnt-1),Range::all())*basis::tet(log2p).vdiag;
	}  
		
	/* APPLY VERTEX DIRICHLET B.C.'S */
	for(i=0;i<nfbd;++i)
		hp_fbdry(i)->vdirichlet();      
			
	//cout <<  basis::tet(log2p).vdiag <<  basis::tet(log2p).ediag <<  basis::tet(log2p).fdiag << endl;            

		
	if(basis::tet(log2p).em == 0) return;
		//cout << basis::tet(log2p).vfms(4,Range::all()) << endl;
	if(basis::tet(log2p).p == 2){
		basis::tet(log2p).ediag(0) =10;
		gbl->res.e(Range(0,nseg-1),0,Range::all()) *= gbl->eprcn(Range(0,nseg-1),Range::all())*basis::tet(log2p).ediag(0);
	
	/* APPLY VERTEX DIRICHLET B.C.'S */
	for(i=0;i<nfbd;++i)
		hp_fbdry(i)->edirichlet();
		
	}

	if(basis::tet(log2p).p == 3){
		basis::tet(log2p).ediag(1) = basis::tet(log2p).ediag(0);

		if (gbl->diagonal_preconditioner) {
			gbl->res.e(Range(0,nseg-1),1,Range::all()) *= gbl->eprcn(Range(0,nseg-1),Range::all())*basis::tet(log2p).ediag(1);
		}
		
		for(i=0;i<nfbd;++i)
			hp_fbdry(i)->edirichlet();			
	}
	
	/* FACE MODES */	         
		if (basis::tet(log2p).fm > 0){
			basis::tet(log2p).fdiag(0) = 20;
			for(find = 0; find<ntri;++find) {
				for (k=0;k<basis::tet(log2p).fm;++k){
					for(n=0;n<NV;++n){
						gbl->res.f(find,k,n) *=gbl->fprcn(find,n)*basis::tet(log2p).fdiag(k); 
					}         
				}   
			}
			for(i=0;i<nfbd;++i)
				hp_fbdry(i)->fdirichlet();

		}
	

					
	return;
}
#endif

void tet_hp::setup_preconditioner() {
	int i,last_phase,mp_phase;
	
	/* SET UP TSTEP FOR MESH MOVEMENT */
//   if (mmovement == coupled_deformable && log2p == 0) {
//      r_tet_mesh::setup_preconditioner();   
//   }
	
	/* SET UP TSTEP FOR ACTIVE BOUNDARIES */
	for(i=0;i<nfbd;++i)
		hp_fbdry(i)->setup_preconditioner();
		
	for(i=0;i<nebd;++i)
		hp_ebdry(i)->setup_preconditioner();
		
	for(i=0;i<nvbd;++i)
		hp_vbdry(i)->setup_preconditioner();
	
	/* SET UP TSTEP FOR HELPER */
	helper->setup_preconditioner();   
	
	if (gbl->diagonal_preconditioner) {
		for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
			pc0load(mp_phase,gbl->vprcn.data());
			pmsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
			last_phase = true;
			last_phase &= pc0wait_rcv(mp_phase,gbl->vprcn.data());
		}
		/* PREINVERT PRECONDITIONER FOR VERTICES */
		gbl->vprcn(Range(0,npnt-1),Range::all()) = 1.0/gbl->vprcn(Range(0,npnt-1),Range::all());
			
		if (basis::tet(log2p).em) {
			for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
				sc0load(mp_phase,gbl->eprcn.data(),0,0,1);
				smsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
				last_phase = true;
				last_phase &= sc0wait_rcv(mp_phase,gbl->eprcn.data(),0,0,1);
			}
		
			/* INVERT DIAGANOL PRECONDITIONER FOR SIDES */            
			gbl->eprcn(Range(0,nseg-1),Range::all()) = 1.0/gbl->eprcn(Range(0,nseg-1),Range::all());
			
			if (basis::tet(log2p).fm) {
				for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
					tc0load(gbl->fprcn.data(),0,0,gbl->fprcn.extent(secondDim));
					tmsgpass(boundary::all,mp_phase,boundary::symmetric);
					last_phase = true;
					last_phase &= tc0wait_rcv(gbl->fprcn.data(),0,0,gbl->fprcn.extent(secondDim));
				}	
				gbl->fprcn(Range(0,ntri-1),Range::all()) = 1.0/gbl->fprcn(Range(0,ntri-1),Range::all());

				if(basis::tet(log2p).im > 0) {
					gbl->iprcn(Range(0,ntet-1),Range::all()) = 1.0/gbl->iprcn(Range(0,ntet-1),Range::all());
				}
			}
		}
	}
//   else {
//      /* NEED STUFF HERE FOR CONTINUITY OF MATRIX PRECONDITIONER */
//      for(int stage = 0; stage<NV; ++stage) {
//         for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
//            pc0load(mp_phase,gbl->vprcn_ut.data() +stage*NV,NV);
//            pmsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
//            last_phase = true;
//            last_phase &= pc0wait_rcv(mp_phase,gbl->vprcn_ut.data()+stage*NV,NV);
//         }
//         if (log2p) {
//            sc0load(gbl->eprcn_ut.data()+stage*NV,0,0,NV);
//            smsgpass(boundary::all,0,boundary::symmetric);
//            sc0wait_rcv(gbl->eprcn_ut.data()+stage*NV,0,0,NV);
//         }
//      }
//      
//      /* FACTORIZE PRECONDITIONER FOR VERTICES ASSUMES LOWER TRIANGULAR NOTHING  */
//      for(i=0;i<npnt;++i)
//         for(int n=0;n<NV;++n)
//            gbl->vprcn_ut(i,n,n) = 1.0/(basis::tri(log2p).vdiag*gbl->vprcn_ut(i,n,n));
//     
//      if (basis::tri(log2p).sm > 0) {
//         /* INVERT DIAGANOL PRECONDITIONER FOR SIDES ASSUMES LOWER TRIANGULAR */    
//         for(i=0;i<nseg;++i)
//            for(int n=0;n<NV;++n)
//               gbl->eprcn_ut(i,n,n)= 1.0/gbl->eprcn_ut(i,n,n);
//      }
//   }

	
	if(basis::tet(log2p).p == 3)
		spoke();
	
	return;
}




void tet_hp::spoke(){
	FLT jcb;
	int sind,ind,ind2,tind,info;
	TinyVector<int,2> v0;
	TinyVector<int,3> e0;
	TinyMatrix<int,6,2> vrtxseg;

	for(int vind = 0;vind < npnt; ++vind){
		tet_mesh::vertexball(vind);
		spkmass(vind) = 0.0;
		ind = 0;
		for(int i = 0; i < pnt(vind).nnbor; ++i){
			tind = gbl->i2wk(i);
			//jcb = tet(tind).vol/8;
			jcb = gbl->iprcn(tind,0);
			ind2 = 0;
			for(int j=0; j < 6; ++j){
				sind=tet(tind).seg(j);
				v0 = seg(sind).pnt;
				for(int k = 0; k < 2; ++k){
					if(vind == v0(k)){
						if(gbl->i1wk(sind) < 0) {
							gbl->i1wk(sind) = ind;
							spklink(vind)(ind)(0) = sind;
							spklink(vind)(ind)(1) = k;					
							++ind;
						}
						e0(ind2++) = gbl->i1wk(sind);						
					}
				}
			}
			for(int j=0; j < 3; ++j){
				spkmass(vind)(e0(j),e0(j)) += jcb*basis::tet(log2p).mdiag;
				for(int k = 0; k < j; ++k)
					spkmass(vind)(e0(j),e0(k)) += jcb*basis::tet(log2p).odiag;
				for(int k = j+1; k < 3; ++k)
					spkmass(vind)(e0(j),e0(k)) += jcb*basis::tet(log2p).odiag;
			}
		}

		GETRF(pnt(vind).nspk,pnt(vind).nspk,&spkmass(vind)(0,0),pnt(vind).nspk,&spkpiv(vind)(0),info); 
		if (info != 0) {
			printf("DGETRF FAILED - VERTEX BALL info:%d vertex:%d \n",info,vind);
			exit(1);
		}

		for(int spk = 0; spk < pnt(vind).nspk; ++spk)
			gbl->i1wk(spklink(vind)(spk)(0))=-1;			
	}
	return;
}

void tet_hp::minvrt_test() {
	int i,j,k,n,find,tind,p0,p1,p2,v0;
	int stridey = MXGP;
	int stridex = MXGP*MXGP;
	FLT jcb,a,h,amax,amin,hmax,hmin,havg,maxres;
	FLT dx1,dy1,dx2,dy2,dz1,dz2,cpi,cpj,cpk;
	TinyVector<int,4> v;

	TinyVector<FLT,3> pt;
	
	gbl->vprcn(Range::all(),Range::all())=0;
	gbl->eprcn(Range::all(),Range::all())=0;
	gbl->fprcn(Range::all(),Range::all())=0;
	
	for(tind=0;tind<ntet;++tind){      
		gbl->iprcn(tind,0) = tet(tind).vol/8; 	
		for(i=0;i<4;++i) 
			gbl->vprcn(tet(tind).pnt(i),0) += gbl->iprcn(tind,0);			
		if (basis::tet(log2p).em > 0) {
			for(i=0;i<6;++i)				
				gbl->eprcn(tet(tind).seg(i),0) += gbl->iprcn(tind,0);		
			if (basis::tet(log2p).fm > 0) {
				for(i=0;i<4;++i)
					gbl->fprcn(tet(tind).tri(i),0) += gbl->iprcn(tind,0);				
			}
		}
	}
	
	hmax = 0;
	hmin = 1000000;
	havg = 0.0;
	for(tind = 0; tind < ntet; ++tind) {
		jcb = tet(tind).vol/8; 
		v = tet(tind).pnt;
		amax = 0.0;
		amin = 1000000;
		for(j=0;j<4;++j) { // FIND MAX FACE AREA AND THEN DIVIDE VOLUME BY IT 
			find = tet(tind).tri(j);
			p0 = tri(find).pnt(0);
			p1 = tri(find).pnt(1);
			p2 = tri(find).pnt(2);
	
			dx1 = pnts(p0)(0)-pnts(p1)(0);
			dy1 = pnts(p0)(1)-pnts(p1)(1);
			dz1 = pnts(p0)(2)-pnts(p1)(2);
			dx2 = pnts(p0)(0)-pnts(p2)(0);
			dy2 = pnts(p0)(1)-pnts(p2)(1);
			dz2 = pnts(p0)(2)-pnts(p2)(2);
			cpi = dy1*dz2-dz1*dy2;
			cpj = -dx1*dz2+dz1*dx2;
			cpk = dx1*dy2-dy1*dx2;
			a =	.5*sqrt(cpi*cpi+cpj*cpj+cpk*cpk);
			amax = (a > amax ? a : amax);
			amin = (a < amin ? a : amin);
		}

		havg += 4.0*jcb/amax;
		h = 4.0*jcb/(0.25*(basis::tet(log2p).p+1)*(basis::tet(log2p).p+1)*amax); // 3*8/6=4
		
		if(4.0*jcb/amin > hmax)
			hmax = 4.0*jcb/amin;
			
		if(4.0*jcb/amax < hmin)
			hmin = 4.0*jcb/amax;
	}
			
		cout << "hmin = " << hmin << "  hmax = " << hmax << endl;
		cout << hmin << ' ' << hmax << ' ' << havg/ntet << endl;
		
		
	#ifndef NODAL	
	if(basis::tet(log2p).p > 2)
		spoke();
	#endif
		
	setup_preconditioner();

	ug.v(Range(0,npnt-1),Range::all()) = 0.0;
	if(basis::tet(log2p).em > 0) {
		ug.e(Range(0,nseg-1),Range::all(),Range::all()) = 0.0;
		
		if(basis::tet(log2p).fm > 0) {
			ug.f(Range(0,ntri-1),Range::all(),Range::all()) = 0.0;
			
			if(basis::tet(log2p).im > 0) {
				ug.i(Range(0,ntet-1),Range::all(),Range::all()) = 0.0;
			}
		}
	}
	

	for(int it = 0; it < 200; ++it){
	
	   gbl->res.v(Range(0,npnt-1),Range::all()) = 0.0;
		if(basis::tet(log2p).em > 0) {
			gbl->res.e(Range(0,nseg-1),Range::all(),Range::all()) = 0.0;
		
			if(basis::tet(log2p).fm > 0) {
				gbl->res.f(Range(0,ntri-1),Range::all(),Range::all()) = 0.0;
			
				if(basis::tet(log2p).im > 0) {
					gbl->res.i(Range(0,ntet-1),Range::all(),Range::all()) = 0.0; 
				}
			}
		}

		for(tind = 0; tind<ntet;++tind) {
			if (tet(tind).info < 0) {
				//cout << "straight" <<  endl;
				for(n=0;n<ND;++n)
					basis::tet(log2p).proj(pnts(tet(tind).pnt(0))(n),pnts(tet(tind).pnt(1))(n),pnts(tet(tind).pnt(2))(n),pnts(tet(tind).pnt(3))(n),&crd(n)(0)(0)(0),stridex,stridey);

				for(i=0;i<basis::tet(log2p).gpx;++i) {
					for(j=0;j<basis::tet(log2p).gpy;++j) {
						for(k=0;k<basis::tet(log2p).gpz;++k) {
							for(n=0;n<ND;++n) {
							dcrd(n)(0)(i)(j)(k) = 0.5*(pnts(tet(tind).pnt(3))(n) -pnts(tet(tind).pnt(2))(n));
							dcrd(n)(1)(i)(j)(k) = 0.5*(pnts(tet(tind).pnt(1))(n) -pnts(tet(tind).pnt(2))(n));
							dcrd(n)(2)(i)(j)(k) = 0.5*(pnts(tet(tind).pnt(0))(n) -pnts(tet(tind).pnt(2))(n));

							}
						}
					}
				}
			}
			else {
				crdtocht(tind);
				for(n=0;n<ND;++n)
					basis::tet(log2p).proj_bdry(&cht(n)(0), &crd(n)(0)(0)(0), &dcrd(n)(0)(0)(0)(0), &dcrd(n)(1)(0)(0)(0),&dcrd(n)(2)(0)(0)(0),stridex,stridey);
				//cout << "curvy" << endl;

			}
			
			for(n=0;n<NV;++n)
				for(i=0;i<basis::tet(log2p).tm;++i)
					lf(n)(i) = 0.0;
			
			ugtouht(tind);
			for(n=0;n<NV;++n)
				basis::tet(log2p).proj(&uht(n)(0),&u(n)(0)(0)(0),stridex,stridey);
			
			for(i=0;i<basis::tet(log2p).gpx;++i) {
				for(j=0;j<basis::tet(log2p).gpy;++j) {
					for(k=0;k<basis::tet(log2p).gpz;++k) {
						pt(0) = crd(0)(i)(j)(k);
						pt(1) = crd(1)(i)(j)(k);
						pt(2) = crd(2)(i)(j)(k);
						cjcb(i)(j)(k) = dcrd(0)(0)(i)(j)(k)*(dcrd(1)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(1)(i)(j)(k))-dcrd(0)(1)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(0)(i)(j)(k))+dcrd(0)(2)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)-dcrd(1)(1)(i)(j)(k)*dcrd(2)(0)(i)(j)(k));
						for(n=0;n<NV;++n)
							res(n)(i)(j)(k) = (u(n)(i)(j)(k)-gbl->ibc->f(n,pt,gbl->time))*cjcb(i)(j)(k);
					}
				}
			}
			for(n=0;n<NV;++n)
				basis::tet(log2p).intgrt(&lf(n)(0),&res(n)(0)(0)(0),stridex,stridey);
			  
			lftog(tind,gbl->res);

		}
		
		minvrt(); 
		//cout << gbl->res.e << endl;
		maxres = 0.0;
		for(int i=0; i < npnt; ++i){
			if(fabs(gbl->res.v(i,0)) > maxres)
				maxres = fabs(gbl->res.v(i,0));
		}
		cout << it+1 << ' ' <<maxres << endl;
		/* Inversion finished */
		ug.v(Range(0,npnt-1),Range::all()) -= gbl->res.v(Range(0,npnt-1),Range::all());
		if (basis::tet(log2p).em > 0) {
			ug.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all()) -= gbl->res.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all());
		}

		if (basis::tet(log2p).fm > 0) {
			ug.f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all()) -= gbl->res.f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all());
		}

		if (basis::tet(log2p).im > 0) {
			ug.i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all()) -= gbl->res.i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all());
		}
		//l2error(gbl->ibc);
//		std::ostringstream filename;
//		filename.str("");
//		filename << "test" << it << std::flush;
//		output(filename.str(),block::display);
	}


			
	return;
}


void tet_hp::minvrt_iter() {
	FLT lcl,diag;
	int dind,tind,vind,sind,find,iind;
	vefi wk;
	wk.v.resize(maxvst,NV);
	wk.e.resize(maxvst,em0,NV);
	wk.f.resize(maxvst,fm0,NV);
	wk.i.resize(maxvst,im0,NV);
	
	int i,j,k,n;
	int stridey = MXGP;
	int stridex = MXGP*MXGP;
	TinyVector<int,4> v;
	TinyVector<FLT,3> pt;
	FLT dt=10;   
	setup_preconditioner();

	ug.v(Range(0,npnt-1),Range::all()) = 0.0;
	if(basis::tet(log2p).em > 0) {
		ug.e(Range(0,nseg-1),Range::all(),Range::all()) = 0.0;
		if(basis::tet(log2p).fm > 0) {
			ug.f(Range(0,ntri-1),Range::all(),Range::all()) = 0.0;
			
			if(basis::tet(log2p).im > 0) {
				ug.i(Range(0,ntet-1),Range::all(),Range::all()) = 0.0;
			}
		}
	}
	
	
	gbl->res.v(Range(0,npnt-1),Range::all()) = 0.0;
	if(basis::tet(log2p).em > 0) {
		gbl->res.e(Range(0,nseg-1),Range::all(),Range::all()) = 0.0;
	
		if(basis::tet(log2p).fm > 0) {
			gbl->res.f(Range(0,ntri-1),Range::all(),Range::all()) = 0.0;
		
			if(basis::tet(log2p).im > 0) {
				gbl->res.i(Range(0,ntet-1),Range::all(),Range::all()) = 0.0; 
			}
		}
	}

	for(tind = 0; tind<ntet;++tind) {
		if (tet(tind).info < 0) {
			for(n=0;n<ND;++n)
				basis::tet(log2p).proj(pnts(tet(tind).pnt(0))(n),pnts(tet(tind).pnt(1))(n),pnts(tet(tind).pnt(2))(n),pnts(tet(tind).pnt(3))(n),&crd(n)(0)(0)(0),stridex,stridey);

			for(i=0;i<basis::tet(log2p).gpx;++i) {
				for(j=0;j<basis::tet(log2p).gpy;++j) {
					for(k=0;k<basis::tet(log2p).gpz;++k) {
						for(n=0;n<ND;++n) {
							dcrd(n)(0)(i)(j)(k) = 0.5*(pnts(tet(tind).pnt(3))(n) -pnts(tet(tind).pnt(2))(n));
							dcrd(n)(1)(i)(j)(k) = 0.5*(pnts(tet(tind).pnt(1))(n) -pnts(tet(tind).pnt(2))(n));
							dcrd(n)(2)(i)(j)(k) = 0.5*(pnts(tet(tind).pnt(0))(n) -pnts(tet(tind).pnt(2))(n));

						}
					}
				}
			}
		}
		else {
			crdtocht(tind);
			for(n=0;n<ND;++n)
				basis::tet(log2p).proj_bdry(&cht(n)(0), &crd(n)(0)(0)(0), &dcrd(n)(0)(0)(0)(0), &dcrd(n)(1)(0)(0)(0),&dcrd(n)(2)(0)(0)(0),stridex,stridey);

		}
		
		for(n=0;n<NV;++n)
			for(i=0;i<basis::tet(log2p).tm;++i)
				lf(n)(i) = 0.0;
		
		for(i=0;i<basis::tet(log2p).gpx;++i) {
			for(j=0;j<basis::tet(log2p).gpy;++j) {
				for(k=0;k<basis::tet(log2p).gpz;++k) {
					pt(0) = crd(0)(i)(j)(k);
					pt(1) = crd(1)(i)(j)(k);
					pt(2) = crd(2)(i)(j)(k);
					cjcb(i)(j)(k) = dcrd(0)(0)(i)(j)(k)*(dcrd(1)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(1)(i)(j)(k))-dcrd(0)(1)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(0)(i)(j)(k))+dcrd(0)(2)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)-dcrd(1)(1)(i)(j)(k)*dcrd(2)(0)(i)(j)(k));
					for(n=0;n<NV;++n)
						res(n)(i)(j)(k) = gbl->ibc->f(n,pt,gbl->time)*cjcb(i)(j)(k);
				}
			}
		}
		for(n=0;n<NV;++n)
			basis::tet(log2p).intgrt(&lf(n)(0),&res(n)(0)(0)(0),stridex,stridey);
		  
		lftog(tind,gbl->res);

	}
	
	for(int it = 0; it < 20; ++it){
		wk.v(Range(0,npnt-1),Range::all()) = ug.v(Range(0,npnt-1),Range::all());
		if (basis::tet(log2p).em > 0) {
			wk.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all()) = ug.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all());
		}	
		if (basis::tet(log2p).fm > 0) {
			wk.f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all()) = ug.f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all());
		}	
		if (basis::tet(log2p).im > 0) {
			wk.i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all()) = ug.i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all());
		}
		for(int i = 0; i < npnt; ++i){	
			for(int n = 0; n < NV; ++n){
				lcl = 0.0;
				diag = 0.0;
				tet_mesh::vertexball(i);
				for(int j = 0; j < pnt(i).nnbor; ++j){
					tind = gbl->i2wk(j);
					for(int k = 0; k < 4; ++k){
						if(i == tet(tind).pnt(k)){
							vind = k;
						}
					}
					diag += basis::tet(log2p).mm(vind,vind);
					for(int k = 0; k < vind; ++k){
						lcl += basis::tet(log2p).mm(vind,k)*ug.v(tet(tind).pnt(k),n);
					}
					for(int k = vind+1; k < 4; ++k){
						lcl += basis::tet(log2p).mm(vind,k)*ug.v(tet(tind).pnt(k),n);
					}
					for(int k = 0; k < 6; ++k){
						for(int m = 0; m < basis::tet(log2p).em; ++m){
							lcl += basis::tet(log2p).mm(vind,4+k*basis::tet(log2p).em+m)*ug.e(tet(tind).seg(k),m,n);
						}
					}
					for(int k = 0; k < 4; ++k){
						for(int m = 0; m < basis::tet(log2p).fm; ++m){
							lcl += basis::tet(log2p).mm(vind,4+6*basis::tet(log2p).em+k*basis::tet(log2p).fm+m)*ug.f(tet(tind).tri(k),m,n);
						}
					}
					for(int m = 0; m < basis::tet(log2p).im; ++m){
						lcl += basis::tet(log2p).mm(vind,basis::tet(log2p).bm+m)*ug.i(tind,m,n);
					}
					lcl *= dt*tet(tind).vol/6;//jacobian
					diag *= dt*tet(tind).vol/6;

				}
				wk.v(i,n) = (gbl->res.v(i,n)-lcl)/diag;// jacobi iteration
			}
		}
		
		for(int i = 0; i < nseg; ++i){	
			for(int n = 0; n < NV; ++n){
				for(int m = 0; m < basis::tet(log2p).em; ++m){
					lcl = 0.0;
					diag = 0.0;
					tet_mesh::ring(i);
					for(int j = 0; j < seg(i).nnbor; ++j){
						tind = gbl->i2wk(j);
						for(int k = 0; k < 6; ++k){
							if(i == tet(tind).seg(k)){
								dind = k;
							}
						}
						sind = 4+basis::tet(log2p).em*dind+m;
						diag += basis::tet(log2p).mm(sind,sind);
						for(int k = 0; k < 4; ++k){
							lcl += basis::tet(log2p).mm(sind,k)*ug.v(tet(tind).pnt(k),n);
						}
						for(int k = 0; k < dind; ++k){
							for(int lm = 0; lm < basis::tet(log2p).em; ++lm){
								lcl += basis::tet(log2p).mm(sind,4+k*basis::tet(log2p).em+lm)*ug.e(tet(tind).seg(k),lm,n);
							}
						}
						for(int k = dind+1; k < 6; ++k){
							for(int lm = 0; lm < basis::tet(log2p).em; ++lm){
								lcl += basis::tet(log2p).mm(sind,4+k*basis::tet(log2p).em+lm)*ug.e(tet(tind).seg(k),lm,n);
							}
						}
						for(int k = 0; k < 4; ++k){
							for(int lm = 0; lm < basis::tet(log2p).fm; ++lm){
								lcl += basis::tet(log2p).mm(sind,4+6*basis::tet(log2p).em+k*basis::tet(log2p).fm+lm)*ug.f(tet(tind).tri(k),lm,n);
							}
						}
						for(int lm = 0; lm < basis::tet(log2p).im; ++lm){
							lcl += basis::tet(log2p).mm(sind,basis::tet(log2p).bm+lm)*ug.i(tind,lm,n);
						}
						lcl *= dt*tet(tind).vol/6;//jacobian
						diag *= dt*tet(tind).vol/6;

					}
					wk.e(i,m,n) = (gbl->res.e(i,m,n)-lcl)/diag;// jacobi iteration
				}
			}
		}
		
		for(int i = 0; i < ntri; ++i){	
			for(int n = 0; n < NV; ++n){
				for(int m = 0; m < basis::tet(log2p).fm; ++m){
					lcl = 0.0;
					diag = 0.0;
					for(int j = 0; j < 2; ++j){
						tind = tri(i).tet(j);
						if(tind != -1){
							for(int k = 0; k < 4; ++k){
								if(i == tet(tind).tri(k)){
									dind = k;
								}
							}
							find = 4+6*basis::tet(log2p).em+dind*basis::tet(log2p).fm+m;
							diag += basis::tet(log2p).mm(find,find);
							for(int k = 0; k < 4; ++k){
								lcl += basis::tet(log2p).mm(find,k)*ug.v(tet(tind).pnt(k),n);
							}
							for(int k = 0; k < 6; ++k){
								for(int lm = 0; lm < basis::tet(log2p).em; ++lm){
									lcl += basis::tet(log2p).mm(find,4+k*basis::tet(log2p).em+lm)*ug.e(tet(tind).seg(k),lm,n);
								}
							}
							for(int k = 0; k < dind; ++k){
								for(int lm = 0; lm < basis::tet(log2p).fm; ++lm){
									lcl += basis::tet(log2p).mm(find,4+6*basis::tet(log2p).em+k*basis::tet(log2p).fm+lm)*ug.f(tet(tind).tri(k),lm,n);
								}
							}
							for(int k = dind+1; k < 4; ++k){
								for(int lm = 0; lm < basis::tet(log2p).fm; ++lm){
									lcl += basis::tet(log2p).mm(find,4+6*basis::tet(log2p).em+k*basis::tet(log2p).fm+lm)*ug.f(tet(tind).tri(k),lm,n);
								}
							}
							for(int lm = 0; lm < basis::tet(log2p).im; ++lm){
								lcl += basis::tet(log2p).mm(find,basis::tet(log2p).bm+lm)*ug.i(tind,lm,n);
							}									
							lcl *= dt*tet(tind).vol/6;//jacobian
							diag *= dt*tet(tind).vol/6;
							
						}
					}
					wk.f(i,m,n) = (gbl->res.f(i,m,n)-lcl)/diag;// jacobi iteration
				}
			}
		}
		
		for(tind = 0; tind < ntet; ++tind){	
			for(int n = 0; n < NV; ++n){
				for(int m = 0; m < basis::tet(log2p).im; ++m){
					lcl = 0.0;
					diag = 0.0;
					iind = basis::tet(log2p).bm+m;
					diag += basis::tet(log2p).mm(iind,iind);
					for(int k = 0; k < 4; ++k){
						lcl += basis::tet(log2p).mm(iind,k)*ug.v(tet(tind).pnt(k),n);
					}
					for(int k = 0; k < 6; ++k){
						for(int lm = 0; lm < basis::tet(log2p).em; ++lm){
							lcl += basis::tet(log2p).mm(iind,4+k*basis::tet(log2p).em+lm)*ug.e(tet(tind).seg(k),lm,n);
						}
					}
					for(int k = 0; k < 4; ++k){
						for(int lm = 0; lm < basis::tet(log2p).fm; ++lm){
							lcl += basis::tet(log2p).mm(iind,4+6*basis::tet(log2p).em+k*basis::tet(log2p).fm+lm)*ug.f(tet(tind).tri(k),lm,n);
						}
					}
					for(int lm = 0; lm < m; ++lm){
						lcl += basis::tet(log2p).mm(iind,basis::tet(log2p).bm+lm)*ug.i(tind,lm,n);
					}
					for(int lm = m+1; lm < basis::tet(log2p).im; ++lm){
						lcl += basis::tet(log2p).mm(iind,basis::tet(log2p).bm+lm)*ug.i(tind,lm,n);
					}					
					lcl *= dt*tet(tind).vol/6;//jacobian
					diag *= dt*tet(tind).vol/6;
					wk.i(tind,m,n) = (gbl->res.i(tind,m,n)-lcl)/diag;// jacobi iteration
				}
			}
		}
		
		ug.v(Range(0,npnt-1),Range::all()) = wk.v(Range(0,npnt-1),Range::all());
		if (basis::tet(log2p).em > 0) {
			ug.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all()) = wk.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all());
		}	
		if (basis::tet(log2p).fm > 0) {
			ug.f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all()) = wk.f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all());
		}	
		if (basis::tet(log2p).im > 0) {
			ug.i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all()) = wk.i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all());
		}
		l2error(gbl->ibc);
		std::ostringstream filename;
		filename.str("");
		filename << "test" << it << std::flush;
		output(filename.str(),block::display);
	}
	return;

}


void tet_hp::minvrt_direct() {
	FLT lcl,jcb;
	int dind,tind,vind,sind,find,iind;
	int i,j,k,n;
	int stridey = MXGP;
	int stridex = MXGP*MXGP;
	int tm = npnt+nseg*basis::tet(log2p).em+ntri*basis::tet(log2p).fm+ntet*basis::tet(log2p).im;
	TinyVector<int,4> v;
	TinyVector<FLT,3> pt;
	setup_preconditioner();
	FLT dt = 1;
	Array<FLT,2> mm(tm,tm);
	mm=0;

	ug.v(Range(0,npnt-1),Range::all()) = 0.0;
	if(basis::tet(log2p).em > 0) {
		ug.e(Range(0,nseg-1),Range::all(),Range::all()) = 0.0;
		if(basis::tet(log2p).fm > 0) {
			ug.f(Range(0,ntri-1),Range::all(),Range::all()) = 0.0;
			
			if(basis::tet(log2p).im > 0) {
				ug.i(Range(0,ntet-1),Range::all(),Range::all()) = 0.0;
			}
		}
	}
	
	
	gbl->res.v(Range(0,npnt-1),Range::all()) = 0.0;
	if(basis::tet(log2p).em > 0) {
		gbl->res.e(Range(0,nseg-1),Range::all(),Range::all()) = 0.0;
	
		if(basis::tet(log2p).fm > 0) {
			gbl->res.f(Range(0,ntri-1),Range::all(),Range::all()) = 0.0;
		
			if(basis::tet(log2p).im > 0) {
				gbl->res.i(Range(0,ntet-1),Range::all(),Range::all()) = 0.0; 
			}
		}
	}

	for(tind = 0; tind<ntet;++tind) {
		if (tet(tind).info < 0) {
			for(n=0;n<ND;++n)
				basis::tet(log2p).proj(pnts(tet(tind).pnt(0))(n),pnts(tet(tind).pnt(1))(n),pnts(tet(tind).pnt(2))(n),pnts(tet(tind).pnt(3))(n),&crd(n)(0)(0)(0),stridex,stridey);

			for(i=0;i<basis::tet(log2p).gpx;++i) {
				for(j=0;j<basis::tet(log2p).gpy;++j) {
					for(k=0;k<basis::tet(log2p).gpz;++k) {
						for(n=0;n<ND;++n) {
							dcrd(n)(0)(i)(j)(k) = 0.5*(pnts(tet(tind).pnt(3))(n) -pnts(tet(tind).pnt(2))(n));
							dcrd(n)(1)(i)(j)(k) = 0.5*(pnts(tet(tind).pnt(1))(n) -pnts(tet(tind).pnt(2))(n));
							dcrd(n)(2)(i)(j)(k) = 0.5*(pnts(tet(tind).pnt(0))(n) -pnts(tet(tind).pnt(2))(n));

						}
					}
				}
			}
		}
		else {
			crdtocht(tind);
			for(n=0;n<ND;++n)
				basis::tet(log2p).proj_bdry(&cht(n)(0), &crd(n)(0)(0)(0), &dcrd(n)(0)(0)(0)(0), &dcrd(n)(1)(0)(0)(0),&dcrd(n)(2)(0)(0)(0),stridex,stridey);

		}
		
		for(n=0;n<NV;++n)
			for(i=0;i<basis::tet(log2p).tm;++i)
				lf(n)(i) = 0.0;
		
		for(i=0;i<basis::tet(log2p).gpx;++i) {
			for(j=0;j<basis::tet(log2p).gpy;++j) {
				for(k=0;k<basis::tet(log2p).gpz;++k) {
					pt(0) = crd(0)(i)(j)(k);
					pt(1) = crd(1)(i)(j)(k);
					pt(2) = crd(2)(i)(j)(k);
					cjcb(i)(j)(k) = dcrd(0)(0)(i)(j)(k)*(dcrd(1)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(1)(i)(j)(k))-dcrd(0)(1)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(0)(i)(j)(k))+dcrd(0)(2)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)-dcrd(1)(1)(i)(j)(k)*dcrd(2)(0)(i)(j)(k));
					for(n=0;n<NV;++n)
						res(n)(i)(j)(k) = gbl->ibc->f(n,pt,gbl->time)*cjcb(i)(j)(k);
				}
			}
		}
		for(n=0;n<NV;++n)
			basis::tet(log2p).intgrt(&lf(n)(0),&res(n)(0)(0)(0),stridex,stridey);
		  
		lftog(tind,gbl->res);

	}
	
	for(int i = 0; i < npnt; ++i){	
		for(int n = 0; n < NV; ++n){
			tet_mesh::vertexball(i);
			for(int j = 0; j < pnt(i).nnbor; ++j){
				tind = gbl->i2wk(j);
				jcb=dt*tet(tind).vol/6;
				for(int k = 0; k < 4; ++k){
					if(i == tet(tind).pnt(k)){
						vind = k;
					}
				}
				for(int k = 0; k < 4; ++k){
					mm(i,tet(tind).pnt(k)) += basis::tet(log2p).mm(vind,k)*ug.v(tet(tind).pnt(k),n);
				}
				for(int k = 0; k < 6; ++k){
					for(int m = 0; m < basis::tet(log2p).em; ++m){
						mm(i,npnt+basis::tet(log2p).em*tet(tind).seg(k)+m) += basis::tet(log2p).mm(vind,4+k*basis::tet(log2p).em+m)*ug.e(tet(tind).seg(k),m,n);
					}
				}
				for(int k = 0; k < 4; ++k){
					for(int m = 0; m < basis::tet(log2p).fm; ++m){
						mm(i,npnt+nseg*basis::tet(log2p).em+tet(tind).tri(k)*basis::tet(log2p).fm+m) += basis::tet(log2p).mm(vind,4+6*basis::tet(log2p).em+k*basis::tet(log2p).fm+m)*ug.f(tet(tind).tri(k),m,n);
					}
				}
				for(int m = 0; m < basis::tet(log2p).im; ++m){
					mm(i,npnt+nseg*basis::tet(log2p).em+ntri*basis::tet(log2p).fm+tind*basis::tet(log2p).im+m) += basis::tet(log2p).mm(vind,basis::tet(log2p).bm+m)*ug.i(tind,m,n);
				}
			}
		}
	}
	
	for(int i = 0; i < nseg; ++i){	
		for(int n = 0; n < NV; ++n){
			for(int m = 0; m < basis::tet(log2p).em; ++m){
				tet_mesh::ring(i);
				for(int j = 0; j < seg(i).nnbor; ++j){
					tind = gbl->i2wk(j);
					jcb=dt*tet(tind).vol/6;
					for(int k = 0; k < 6; ++k){
						if(i == tet(tind).seg(k)){
							dind = k;
						}
					}
					sind = 4+basis::tet(log2p).em*dind+m;
					for(int k = 0; k < 4; ++k){
						lcl += basis::tet(log2p).mm(sind,k)*ug.v(tet(tind).pnt(k),n);
					}
					for(int k = 0; k < 6; ++k){
						for(int lm = 0; lm < basis::tet(log2p).em; ++lm){
							lcl += basis::tet(log2p).mm(sind,4+k*basis::tet(log2p).em+lm)*ug.e(tet(tind).seg(k),lm,n);
						}
					}
					for(int k = 0; k < 4; ++k){
						for(int lm = 0; lm < basis::tet(log2p).fm; ++lm){
							lcl += basis::tet(log2p).mm(sind,4+6*basis::tet(log2p).em+k*basis::tet(log2p).fm+lm)*ug.f(tet(tind).tri(k),lm,n);
						}
					}
					for(int lm = 0; lm < basis::tet(log2p).im; ++lm){
						lcl += basis::tet(log2p).mm(sind,basis::tet(log2p).bm+lm)*ug.i(tind,lm,n);
					}
				}
			}
		}
	}
	
	for(int i = 0; i < ntri; ++i){	
		for(int n = 0; n < NV; ++n){
			for(int m = 0; m < basis::tet(log2p).fm; ++m){
				for(int j = 0; j < 2; ++j){
					tind = tri(i).tet(j);
					jcb=dt*tet(tind).vol/6;
					if(tind != -1){
						for(int k = 0; k < 4; ++k){
							if(i == tet(tind).tri(k)){
								dind = k;
							}
						}
						find = 4+6*basis::tet(log2p).em+dind*basis::tet(log2p).fm+m;
						for(int k = 0; k < 4; ++k){
							lcl += basis::tet(log2p).mm(find,k)*ug.v(tet(tind).pnt(k),n);
						}
						for(int k = 0; k < 6; ++k){
							for(int lm = 0; lm < basis::tet(log2p).em; ++lm){
								lcl += basis::tet(log2p).mm(find,4+k*basis::tet(log2p).em+lm)*ug.e(tet(tind).seg(k),lm,n);
							}
						}
						for(int k = 0; k < 4; ++k){
							for(int lm = 0; lm < basis::tet(log2p).fm; ++lm){
								lcl += basis::tet(log2p).mm(find,4+6*basis::tet(log2p).em+k*basis::tet(log2p).fm+lm)*ug.f(tet(tind).tri(k),lm,n);
							}
						}
						for(int lm = 0; lm < basis::tet(log2p).im; ++lm){
							lcl += basis::tet(log2p).mm(find,basis::tet(log2p).bm+lm)*ug.i(tind,lm,n);
						}					
					}
				}
			}
		}
	}
	
	for(tind = 0; tind < ntet; ++tind){	
		for(int n = 0; n < NV; ++n){
			for(int m = 0; m < basis::tet(log2p).im; ++m){
				jcb=dt*tet(tind).vol/6;
				iind = basis::tet(log2p).bm+m;
				for(int k = 0; k < 4; ++k){
					lcl += basis::tet(log2p).mm(iind,k)*ug.v(tet(tind).pnt(k),n);
				}
				for(int k = 0; k < 6; ++k){
					for(int lm = 0; lm < basis::tet(log2p).em; ++lm){
						lcl += basis::tet(log2p).mm(iind,4+k*basis::tet(log2p).em+lm)*ug.e(tet(tind).seg(k),lm,n);
					}
				}
				for(int k = 0; k < 4; ++k){
					for(int lm = 0; lm < basis::tet(log2p).fm; ++lm){
						lcl += basis::tet(log2p).mm(iind,4+6*basis::tet(log2p).em+k*basis::tet(log2p).fm+lm)*ug.f(tet(tind).tri(k),lm,n);
					}
				}
				for(int lm = 0; lm < basis::tet(log2p).im; ++lm){
					lcl += basis::tet(log2p).mm(iind,basis::tet(log2p).bm+lm)*ug.i(tind,lm,n);
				}
			}
		}
	}		

	l2error(gbl->ibc);
	return;

}

