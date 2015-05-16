#include "tri_hp_cns.h"
#include "../hp_boundary.h"
#include <myblas.h>

void tri_hp_cns::project_new_variables(){
	int info,last_phase, mp_phase;
	char uplo[] = "U";
	Array<double,1> lcl(NV), lclug(NV),lclres(NV);
	Array<TinyVector<double,MXGP>,2> P(NV,NV);
	Array<TinyMatrix<double,MXGP,MXGP>,2> P2d(NV,NV);
	Array<TinyVector<double,MXGP>,1> u1d(NV),res1d(NV),temp1d(NV);
	Array<TinyMatrix<double,MXGP,MXGP>,1> u2d(NV),res2d(NV),temp2d(NV);
	Array<TinyVector<double,MXTM>,1> ucoef(NV),rcoef(NV),tcoef(NV);
	
	
	gbl->res_temp.v = gbl->res.v;
	gbl->res_temp.s = gbl->res.s;
	gbl->res_temp.i = gbl->res.i;
	
	/* LOOP THROUGH VERTICES */
	for(int i=0;i<npnt;++i){
	
		for(int n = 0; n < NV; ++n){
			lclug(n) = ug.v(i,n);
			lclres(n) = gbl->res.v(i,n);
		}

		switch_variables(lclug,lclres);

		for(int j=0;j<NV;++j){
			FLT lcl0 = lclres(j);
			for(int k=0;k<j;++k){
				lcl0 -= gbl->vpreconditioner(i,j,k)*lclres(k);
			}
			lclres(j) = lcl0/gbl->vpreconditioner(i,j,j);
		}
		
		for(int n = 0; n < NV; ++n)
			gbl->res.v(i,n) = lclres(n);
				
	}

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
		hp_vbdry(i)->vdirichlet();
	
	if (!basis::tri(log2p)->sm()) return;

	/* LOOP THROUGH SIDES */    
    for(int sind=0;sind<nseg;++sind) {
		
		/* project linears */
		for(int n=0;n<NV;++n)
			basis::tri(log2p)->proj1d(gbl->res.v(seg(sind).pnt(0),n),gbl->res.v(seg(sind).pnt(1),n),&temp1d(n)(0));
		
		for(int m=0;m<NV;++m) 
			for(int n=0;n<NV;++n) 
				basis::tri(log2p)->proj1d(gbl->vpreconditioner(seg(sind).pnt(0),m,n),gbl->vpreconditioner(seg(sind).pnt(1),m,n),&P(m,n)(0));
		
		/* take global coefficients and put into local vector */
		for(int n=0;n<NV;++n) {
			for (int m=0; m<2; ++m) {
				ucoef(n)(m) = ug.v(seg(sind).pnt(m),n);
				rcoef(n)(m) = gbl->res_temp.v(seg(sind).pnt(m),n);
			}
		}
		
		for(int n=0;n<NV;++n) {
			for(int m=0;m<basis::tri(log2p)->sm();++m) {
				ucoef(n)(m+2) = ug.s(sind,m,n);
				rcoef(n)(m+2) = gbl->res_temp.s(sind,m,n);
			}
		}

		/* project sides */
		for(int n=0;n<NV;++n) {
			basis::tri(log2p)->proj1d(&ucoef(n)(0),&u1d(n)(0));
			basis::tri(log2p)->proj1d(&rcoef(n)(0),&res1d(n)(0));
		}
		
		for(int i=0;i<basis::tri(log2p)->gpx(); ++i) {
			
			for(int m = 0; m < NV; ++m){
				lclug(m) = u1d(m)(i);
				lclres(m) = res1d(m)(i);
			}

			switch_variables(lclug,lclres);

			for(int j=0;j<NV;++j){
				FLT lcl0 = lclres(j);
				for(int k=0;k<j;++k){
					lcl0 -= P(j,k)(i)*lclres(k);
				}
				lclres(j) = lcl0/P(j,j)(i);
			}

			for(int m = 0; m < NV; ++m)
				temp1d(m)(i) -= lclres(m);	
			
		}
		
		/* integrate right hand side: phi*P*res */
		for(int n=0;n<NV;++n)
			basis::tri(log2p)->intgrt1d(&rcoef(n)(0),&temp1d(n)(0));
		
		/* invert 1d mass matrix */
		for(int n=0;n<NV;++n) {
			PBTRS(uplo,basis::tri(log2p)->sm(),basis::tri(log2p)->sbwth(),1,(double *) &basis::tri(log2p)->sdiag1d(0,0),basis::tri(log2p)->sbwth()+1,&rcoef(n)(2),basis::tri(log2p)->sm(),info);
			for(int m=0;m<basis::tri(log2p)->sm();++m) 
				gbl->res.s(sind,m,n) = -rcoef(n)(m+2);
			
		}

	}
	
	for(int mode = 0; mode < basis::tri(log2p)->sm(); ++mode) {
		
		sc0load(gbl->res.s.data(),mode,mode,gbl->res.s.extent(secondDim));
		smsgpass(boundary::all,0,boundary::symmetric);
		sc0wait_rcv(gbl->res.s.data(),mode,mode,gbl->res.s.extent(secondDim));
		
		/* APPLY DIRCHLET B.C.S TO MODE */
		for(int i=0;i<nebd;++i)
			hp_ebdry(i)->sdirichlet(mode);
	}
	
	if (!basis::tri(log2p)->im()) return;
	
	for(int tind = 0; tind < ntri; ++tind) {
		
		for (int i=0; i<3; ++i) {
			int vrtx = tri(tind).pnt(i);
			for(int n=0; n<NV; ++n){
				ucoef(n)(i) = ug.v(vrtx,n);
				rcoef(n)(i) = gbl->res_temp.v(vrtx,n);
				tcoef(n)(i) = gbl->res.v(vrtx,n);
			}
		}
		
		/* SIDES */
		int cnt = 3;
		for(int i=0;i<3;++i) {
			int sind = tri(tind).seg(i);
			int sign = tri(tind).sgn(i);
			int msgn = 1;
			for (int m = 0; m < basis::tri(log2p)->sm(); ++m) {
				for(int n=0; n<NV; ++n){
					ucoef(n)(cnt) = msgn*ug.s(sind,m,n);
					rcoef(n)(cnt) = msgn*gbl->res_temp.s(sind,m,n);
					tcoef(n)(cnt) = msgn*gbl->res.s(sind,m,n);
				}
				msgn *= sign;
				++cnt;
			}
		}
		
		/* INTERIORS */    
		if (basis::tri(log2p)->im() > 0) {    
			int indx = 0;
			for(int m = 1; m < basis::tri(log2p)->sm(); ++m) {
				for(int k = 0; k < basis::tri(log2p)->sm()-m; ++k) {
					for(int n=0; n<NV; ++n){
						ucoef(n)(cnt) = ug.i(tind,indx,n);
						rcoef(n)(cnt) = gbl->res_temp.i(tind,indx,n);
						tcoef(n)(cnt) = 0.0;
					}
					++cnt; ++indx;
				}
				indx += sm0 -basis::tri(log2p)->sm();
			}
		}
		
		
		for(int n=0;n<NV;++n){			
			basis::tri(log2p)->proj(&ucoef(n)(0),&u2d(n)(0,0),MXGP);
			basis::tri(log2p)->proj(&rcoef(n)(0),&res2d(n)(0,0),MXGP);
			basis::tri(log2p)->proj_bdry(&tcoef(n)(0),&temp2d(n)(0,0),MXGP);
		}
		
		for(int m=0;m<NV;++m)		
			for(int n=0;n<NV;++n)
				basis::tri(log2p)->proj(gbl->vpreconditioner(tri(tind).pnt(0),m,n),gbl->vpreconditioner(tri(tind).pnt(1),m,n),gbl->vpreconditioner(tri(tind).pnt(2),m,n),&P2d(m,n)(0,0),MXGP);

		
		for (int i=0; i < basis::tri(log2p)->gpx(); ++i ) {
			for (int j=0; j < basis::tri(log2p)->gpn(); ++j ) {

				for(int m = 0; m < NV; ++m){
					lclug(m) = u2d(m)(i,j);
					lclres(m) = res2d(m)(i,j);
				}
	
				switch_variables(lclug,lclres);
				
				for(int m=0;m<NV;++m){
					FLT lcl0 = lclres(m);
					for(int n=0;n<m;++n){
						lcl0 -= P2d(m,n)(i,j)*lclres(n);
					}
					lclres(m) = lcl0/P2d(m,m)(i,j);
				}
				
				for(int m = 0; m < NV; ++m)
					temp2d(m)(i,j) -= lclres(m);	
				

			}
		}

		for(int n=0;n<NV;++n) {
			basis::tri(log2p)->intgrt(&rcoef(n)(0),&temp2d(n)(0,0),MXGP);
			DPBTRS(uplo,basis::tri(log2p)->im(),basis::tri(log2p)->ibwth(),1,(double *) &basis::tri(log2p)->idiag(0,0),basis::tri(log2p)->ibwth()+1,&rcoef(n)(basis::tri(log2p)->bm()),basis::tri(log2p)->im(),info);
			for(int i=0;i<basis::tri(log2p)->im();++i)
				gbl->res.i(tind,i,n) = -rcoef(n)(basis::tri(log2p)->bm()+i);
		}
	}
	
	return;
}

void tri_hp_cns::switch_variables(Array<double,1> pvu, Array<double,1> &a){
	
	Array<double,2> dpdc(NV,NV);
	Array<double,1> temp(NV);
	double gm1 = gbl->gamma-1.0;

	double pr = pvu(0),u = pvu(1),v = pvu(2),rt = pvu(3);
	double rho = pr/rt;
	double ke = 0.5*(u*u+v*v);
	 
	/* jacobian derivative of primitive wrt conservative */
	dpdc = ke*gm1,          -u*gm1,     -v*gm1,     gm1,
		   -u/rho,          1.0/rho,    0.0,        0.0,
		   -v/rho,          0.0,        1.0/rho,    0.0,
		   (gm1*ke-rt)/rho, -u*gm1/rho, -v*gm1/rho, gm1/rho;

	temp = 0.0;
	for(int i = 0; i < NV; ++i)
		for(int j = 0; j < NV; ++j)
			temp(i) += dpdc(i,j)*a(j);

	a = temp;

	return;
}




void tri_hp_cns::update() {
	int i,m,k,n,indx,indx1;
	FLT cflalpha;
	
	// temp fix need to better incorporate more solvers
#ifdef petsc
	tri_hp::petsc_update();
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

		tri_hp::rsdl(stage);
				
		tri_hp::minvrt();
		
		for(i=0;i<nebd;++i)
			hp_ebdry(i)->modify_boundary_residual();		
		
		project_new_variables();
			
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

	}
	
	return;
}


void tri_hp_cns::project_res_vertex() {
	
	Array<double,1> lclug(NV),lclres(NV);
	
	/* LOOP THROUGH VERTICES */
	for(int i=0;i<npnt;++i){
		
		for(int n = 0; n < NV; ++n){
			lclug(n) = ug.v(i,n);
			lclres(n) = gbl->res.v(i,n);
		}
		
		switch_variables(lclug,lclres);
		
		for(int j=0;j<NV;++j){
			FLT lcl0 = lclres(j);
			for(int k=0;k<j;++k){
				lcl0 -= gbl->vpreconditioner(i,j,k)*lclres(k);
			}
			lclres(j) = lcl0/gbl->vpreconditioner(i,j,j);
		}
		
		for(int n = 0; n < NV; ++n)
			gbl->res.v(i,n) = lclres(n);
		
	}
	
	return;
}

void tri_hp_cns::project_res_side(int mode) {
	
	int info;
	char uplo[] = "U";
	Array<double,1> lclug(NV),lclres(NV);
	Array<TinyVector<double,MXGP>,2> P(NV,NV);
	Array<TinyVector<double,MXGP>,1> u1d(NV),res1d(NV),temp1d(NV);
	Array<TinyVector<double,MXTM>,1> ucoef(NV),rcoef(NV),tcoef(NV);
	
	/* LOOP THROUGH SIDES */    
    for(int sind=0;sind<nseg;++sind) {
		
		/* project linears */
		for(int n=0;n<NV;++n)
			basis::tri(log2p)->proj1d(gbl->res.v(seg(sind).pnt(0),n),gbl->res.v(seg(sind).pnt(1),n),&temp1d(n)(0));
		
		for(int m=0;m<NV;++m) 
			for(int n=0;n<NV;++n) 
				basis::tri(log2p)->proj1d(gbl->vpreconditioner(seg(sind).pnt(0),m,n),gbl->vpreconditioner(seg(sind).pnt(1),m,n),&P(m,n)(0));
		
		/* take global coefficients and put into local vector */
		for(int n=0;n<NV;++n) {
			for (int m=0; m<2; ++m) {
				ucoef(n)(m) = ug.v(seg(sind).pnt(m),n);
				rcoef(n)(m) = gbl->res_temp.v(seg(sind).pnt(m),n);
			}
		}
		
		for(int n=0;n<NV;++n) {
			for(int m=0;m<basis::tri(log2p)->sm();++m) {
				ucoef(n)(m+2) = gbl->ug0.s(sind,m,n);
				rcoef(n)(m+2) = gbl->res_temp.s(sind,m,n);
			}
		}
		
		/* project sides */
		for(int n=0;n<NV;++n) {
			basis::tri(log2p)->proj1d(&ucoef(n)(0),&u1d(n)(0));
			basis::tri(log2p)->proj1d(&rcoef(n)(0),&res1d(n)(0));
		}
		
		for(int i=0;i<basis::tri(log2p)->gpx(); ++i) {
			
			for(int m = 0; m < NV; ++m){
				lclug(m) = u1d(m)(i);
				lclres(m) = res1d(m)(i);
			}
			
			switch_variables(lclug,lclres);
			
			for(int j=0;j<NV;++j){
				FLT lcl0 = lclres(j);
				for(int k=0;k<j;++k){
					lcl0 -= P(j,k)(i)*lclres(k);
				}
				lclres(j) = lcl0/P(j,j)(i);
			}
			
			for(int m = 0; m < NV; ++m)
				temp1d(m)(i) -= lclres(m);	
			
		}
		
		/* integrate right hand side: phi*P*res */
		for(int n=0;n<NV;++n)
			basis::tri(log2p)->intgrt1d(&rcoef(n)(0),&temp1d(n)(0));
		
		/* invert 1d mass matrix */
		for(int n=0;n<NV;++n) {
			PBTRS(uplo,basis::tri(log2p)->sm(),basis::tri(log2p)->sbwth(),1,(double *) &basis::tri(log2p)->sdiag1d(0,0),basis::tri(log2p)->sbwth()+1,&rcoef(n)(2),basis::tri(log2p)->sm(),info);
			//			for(int m=0;m<basis::tri(log2p)->sm();++m) 
			//				gbl->res.s(sind,m,n) = -rcoef(n)(m+2);
			gbl->res.s(sind,mode,n) = -rcoef(n)(mode+2);
			
		}
	}
	
	return;
}

void tri_hp_cns::project_res_interior() {
	
	int info;
	char uplo[] = "U";
	Array<double,1> lclug(NV),lclres(NV);
	Array<TinyMatrix<double,MXGP,MXGP>,2> P2d(NV,NV);
	Array<TinyMatrix<double,MXGP,MXGP>,1> u2d(NV),res2d(NV),temp2d(NV);
	Array<TinyVector<double,MXTM>,1> ucoef(NV),rcoef(NV),tcoef(NV);
	
	for(int tind = 0; tind < ntri; ++tind) {
		
		for (int i=0; i<3; ++i) {
			int vrtx = tri(tind).pnt(i);
			for(int n=0; n<NV; ++n){
				ucoef(n)(i) = ug.v(vrtx,n);
				rcoef(n)(i) = gbl->res_temp.v(vrtx,n);
				tcoef(n)(i) = gbl->res.v(vrtx,n);
			}
		}
		
		/* SIDES */
		int cnt = 3;
		for(int i=0;i<3;++i) {
			int sind = tri(tind).seg(i);
			int sign = tri(tind).sgn(i);
			int msgn = 1;
			for (int m = 0; m < basis::tri(log2p)->sm(); ++m) {
				for(int n=0; n<NV; ++n){
					ucoef(n)(cnt) = msgn*ug.s(sind,m,n);
					rcoef(n)(cnt) = msgn*gbl->res_temp.s(sind,m,n);
					tcoef(n)(cnt) = msgn*gbl->res.s(sind,m,n);
				}
				msgn *= sign;
				++cnt;
			}
		}
		
		/* INTERIORS */    
		if (basis::tri(log2p)->im() > 0) {    
			int indx = 0;
			for(int m = 1; m < basis::tri(log2p)->sm(); ++m) {
				for(int k = 0; k < basis::tri(log2p)->sm()-m; ++k) {
					for(int n=0; n<NV; ++n){
						ucoef(n)(cnt) = ug.i(tind,indx,n);
						rcoef(n)(cnt) = gbl->res_temp.i(tind,indx,n);
						tcoef(n)(cnt) = 0.0;
					}
					++cnt; ++indx;
				}
				indx += sm0 -basis::tri(log2p)->sm();
			}
		}
		
		
		for(int n=0;n<NV;++n){			
			basis::tri(log2p)->proj(&ucoef(n)(0),&u2d(n)(0,0),MXGP);
			basis::tri(log2p)->proj(&rcoef(n)(0),&res2d(n)(0,0),MXGP);
			basis::tri(log2p)->proj_bdry(&tcoef(n)(0),&temp2d(n)(0,0),MXGP);
		}
		
		for(int m=0;m<NV;++m)		
			for(int n=0;n<NV;++n)
				basis::tri(log2p)->proj(gbl->vpreconditioner(tri(tind).pnt(0),m,n),gbl->vpreconditioner(tri(tind).pnt(1),m,n),gbl->vpreconditioner(tri(tind).pnt(2),m,n),&P2d(m,n)(0,0),MXGP);
		
		
		for (int i=0; i < basis::tri(log2p)->gpx(); ++i ) {
			for (int j=0; j < basis::tri(log2p)->gpn(); ++j ) {
				
				for(int m = 0; m < NV; ++m){
					lclug(m) = u2d(m)(i,j);
					lclres(m) = res2d(m)(i,j);
				}
				
				switch_variables(lclug,lclres);
				
				for(int m=0;m<NV;++m){
					FLT lcl0 = lclres(m);
					for(int n=0;n<m;++n){
						lcl0 -= P2d(m,n)(i,j)*lclres(n);
					}
					lclres(m) = lcl0/P2d(m,m)(i,j);
				}
				
				for(int m = 0; m < NV; ++m)
					temp2d(m)(i,j) -= lclres(m);	
				
				
			}
		}
		
		for(int n=0;n<NV;++n) {
			basis::tri(log2p)->intgrt(&rcoef(n)(0),&temp2d(n)(0,0),MXGP);
			DPBTRS(uplo,basis::tri(log2p)->im(),basis::tri(log2p)->ibwth(),1,(double *) &basis::tri(log2p)->idiag(0,0),basis::tri(log2p)->ibwth()+1,&rcoef(n)(basis::tri(log2p)->bm()),basis::tri(log2p)->im(),info);
			for(int i=0;i<basis::tri(log2p)->im();++i)
				gbl->res.i(tind,i,n) = -rcoef(n)(basis::tri(log2p)->bm()+i);
		}
	}
	
	return;
}



//void tri_hp_cns::minvrt() {
//	int i,j,k,m,n,tind,sind,v0,indx,indx1,indx2,sgn,msgn;
//	TinyVector<int,3> sign,side;
//	Array<FLT,2> Pinv(NV,NV);
//	Array<FLT,1> temp(NV);
//	int last_phase, mp_phase;
//	
//	/* LOOP THROUGH SIDES */
//	if (basis::tri(log2p)->sm() > 0) {
//		indx = 0;
//		for(sind = 0; sind<nseg;++sind) {
//			/* SUBTRACT SIDE CONTRIBUTIONS TO VERTICES */            
//			for (k=0; k <basis::tri(log2p)->sm(); ++k) {
//				for (i=0; i<2; ++i) {
//					v0 = seg(sind).pnt(i);
//					for(n=0;n<NV;++n)
//						gbl->res.v(v0,n) -= basis::tri(log2p)->sfmv(i,k)*gbl->res.s(sind,k,n);
//				}
//				++indx;
//			}
//		}
//		
//		if (basis::tri(log2p)->im() > 0) {
//			/* SUBTRACT INTERIORS */
//			indx = 0;
//			for(tind = 0; tind<ntri;++tind) {
//				indx2 = 3;
//				for (i=0; i<3; ++i) {
//					v0 = tri(tind).pnt(i);
//					for (k=0;k<basis::tri(log2p)->im();++k)
//						for(n=0;n<NV;++n)
//							gbl->res.v(v0,n) -= basis::tri(log2p)->ifmb(i,k)*gbl->res.i(tind,k,n);
//					
//					sind = tri(tind).seg(i);
//					sgn = tri(tind).sgn(i);
//					msgn = 1;
//					for (j=0;j<basis::tri(log2p)->sm();++j) {
//						for (k=0;k<basis::tri(log2p)->im();++k)
//							for(n=0;n<NV;++n)
//								gbl->res.s(sind,j,n) -= msgn*basis::tri(log2p)->ifmb(indx2,k)*gbl->res.i(tind,k,n);
//						msgn *= sgn;
//						++indx2;
//					}
//				}
//				indx += basis::tri(log2p)->im();
//			}
//		}
//	}
//	
//	
//	if (gbl->diagonal_preconditioner) {
//		gbl->res.v(Range(0,npnt-1),Range::all()) *= gbl->vprcn(Range(0,npnt-1),Range::all());
//    } else {
//		
//		for(i=0;i<npnt;++i) {
//
//			Pinv = gbl->vprcn_ut(i,Range::all(),Range::all());
//			temp = gbl->res.v(i,Range::all());
//
//			int info,ipiv[NV];
//			GETRF(NV, NV, Pinv.data(), NV, ipiv, info);
//			
//			if (info != 0) {
//				*gbl->log << "DGETRF FAILED FOR CNS MINVRT" << std::endl;
//				sim::abort(__LINE__,__FILE__,gbl->log);
//			}
//			
//			char trans[] = "T";
//			GETRS(trans,NV,1,Pinv.data(),NV,ipiv,temp.data(),NV,info);
//			
//			if (info != 0) {
//				*gbl->log << "DGETRS FAILED FOR CNS MINVRT" << std::endl;
//				sim::abort(__LINE__,__FILE__,gbl->log);
//			}
//
//			gbl->res.v(i,Range::all()) = temp/basis::tri(log2p)->vdiag();
//		}
//	}
//	
//	gbl->res_temp.v = gbl->res.v;
//
////	for(i=0;i<nebd;++i)
////		hp_ebdry(i)->vdirichlet();	
//	for(i=0;i<nebd;++i)
//		hp_ebdry(i)->modify_vertex_residual();	
//
//	if(basis::tri(log2p)->sm() == 0) return;
//	
//	
//	/* REMOVE VERTEX CONTRIBUTION FROM SIDE MODES */
//	/* SOLVE FOR SIDE MODES */
//	/* PART 1 REMOVE VERTEX CONTRIBUTIONS */
//	
//	if (gbl->diagonal_preconditioner) {  // IF STATEMENT IN LOOP IS BAD FIXME
//		for(tind=0;tind<ntri;++tind) {
//			for(i=0;i<3;++i) {
//				v0 = tri(tind).pnt(i);
//				for(n=0;n<NV;++n)
//					uht(n)(i) = gbl->res.v(v0,n)*gbl->tprcn(tind,n);
//			}
//			
//			for(i=0;i<3;++i) {
//				sind = tri(tind).seg(i);
//				sgn  = tri(tind).sgn(i);
//				for(j=0;j<3;++j) {
//					indx1 = (i+j)%3;
//					msgn = 1;
//					for(k=0;k<basis::tri(log2p)->sm();++k) {
//						for(n=0;n<NV;++n)
//							gbl->res.s(sind,k,n) -= msgn*basis::tri(log2p)->vfms(j,k)*uht(n)(indx1);
//						msgn *= sgn;
//					}
//				}
//			}
//		}
//	}
//	else {
//		/* THIS IS TO USE A MATRIX PRECONDITIONER */
//		for(tind=0;tind<ntri;++tind) {
//			for(i=0;i<3;++i) {
//				v0 = tri(tind).pnt(i);
//				for(n=0;n<NV;++n)
//					uht(n)(i) = gbl->res.v(v0,n);
//			}
//			
//			for(i=0;i<3;++i) {
//				indx = tri(tind).seg(i);
//				sgn  = tri(tind).sgn(i);
//				for(j=0;j<3;++j) {
//					indx1 = (i+j)%3;
//					msgn = 1;
//					for(k=0;k<basis::tri(log2p)->sm();++k) {
//						for(n=0;n<NV;++n) {
//							for(m=0;m<NV;++m) {
//								gbl->res.s(indx,k,n) -= msgn*basis::tri(log2p)->vfms(j,k)*gbl->tprcn_ut(tind,n,m)*uht(m)(indx1);
//							}
//						}
//						msgn *= sgn;
//					}
//				}
//			}
//		}
//	} 
//
//	gbl->res_temp.s = gbl->res.s;
//
//	for (int mode = 0; mode < basis::tri(log2p)->sm()-1; ++ mode) {
//		/* SOLVE FOR SIDE MODE */
//		if (gbl->diagonal_preconditioner) {
//			gbl->res.s(Range(0,nseg-1),mode,Range::all()) *= gbl->sprcn(Range(0,nseg-1),Range::all())*basis::tri(log2p)->sdiag(mode);
//		}
//		else {
//			for(sind = 0; sind < nseg; ++sind) {
//								
//				Pinv = gbl->sprcn_ut(sind,Range::all(),Range::all());
//				temp = gbl->res.s(sind,mode,Range::all());
//				
//				int info,ipiv[NV];
//				GETRF(NV, NV, Pinv.data(), NV, ipiv, info);
//				
//				if (info != 0) {
//					*gbl->log << "DGETRF FAILED FOR CNS MINVRT" << std::endl;
//					sim::abort(__LINE__,__FILE__,gbl->log);
//				}
//				
//				char trans[] = "T";
//				GETRS(trans,NV,1,Pinv.data(),NV,ipiv,temp.data(),NV,info);
//				
//				if (info != 0) {
//					*gbl->log << "DGETRS FAILED FOR CNS MINVRT" << std::endl;
//					sim::abort(__LINE__,__FILE__,gbl->log);
//				}
//				
//				gbl->res.s(sind,mode,Range::all()) = temp/basis::tri(log2p)->sdiag(mode);
//				
//			}
//		}
//		
//		sc0load(gbl->res.s.data(),mode,mode,gbl->res.s.extent(secondDim));
//		smsgpass(boundary::all,0,boundary::symmetric);
//		sc0wait_rcv(gbl->res.s.data(),mode,mode,gbl->res.s.extent(secondDim));
//
//		for(i=0;i<nebd;++i)
//			hp_ebdry(i)->sdirichlet(mode);
//		
//		/* REMOVE MODE FROM HIGHER MODES */
//		for(tind=0;tind<ntri;++tind) {
//			
//			if (gbl->diagonal_preconditioner) {
//				for(i=0;i<3;++i) {
//					side(i) = tri(tind).seg(i);
//					sign(i) = tri(tind).sgn(i);
//					sgn      = (mode % 2 ? sign(i) : 1);
//					for(n=0;n<NV;++n)
//						uht(n)(i) = sgn*gbl->res.s(side(i),mode,n)*gbl->tprcn(tind,n);
//				}
//				
//				/* REMOVE MODES J,K FROM MODE I,M */
//				for(i=0;i<3;++i) {
//					msgn = ( (mode +1) % 2 ? sign(i) : 1);
//					for(m=mode+1;m<basis::tri(log2p)->sm();++m) {
//						for(j=0;j<3;++j) {
//							indx = (i+j)%3;
//							for(n=0;n<NV;++n) {
//								gbl->res.s(side(i),m,n) -= msgn*basis::tri(log2p)->sfms(mode,m,j)*uht(n)(indx);
//							}
//						}
//						msgn *= sign(i);
//					}
//				}
//			}
//			else {
//				for(i=0;i<3;++i) {
//					side(i) = tri(tind).seg(i);
//					sign(i) = tri(tind).sgn(i);
//					sgn      = (mode % 2 ? sign(i) : 1);
//					for(n=0;n<NV;++n)
//						uht(n)(i) = sgn*gbl->res.s(side(i),mode,n);
//				}
//				
//				/* REMOVE MODES J,K FROM MODE I,M */
//				for(i=0;i<3;++i) {
//					msgn = ( (mode +1) % 2 ? sign(i) : 1);
//					for(m=mode+1;m<basis::tri(log2p)->sm();++m) {
//						for(j=0;j<3;++j) {
//							indx = (i+j)%3;
//							for(n=0;n<NV;++n) {
//								gbl->res.s(side(i),m,n) -= msgn*basis::tri(log2p)->sfms(mode,m,j)*gbl->tprcn_ut(tind,n,0)*uht(0)(indx);
//								for(k=1;k<NV;++k) {
//									gbl->res.s(side(i),m,n) -= msgn*basis::tri(log2p)->sfms(mode,m,j)*gbl->tprcn_ut(tind,n,k)*uht(k)(indx);
//								}
//							}
//						}
//						msgn *= sign(i);
//					}
//				}
//			}
//		}
//	}
//	
//	
//	
//	/* SOLVE FOR HIGHEST MODE */
//	int mode = basis::tri(log2p)->sm()-1;		
//	if (gbl->diagonal_preconditioner) {
//		gbl->res.s(Range(0,nseg-1),mode,Range::all()) *= gbl->sprcn(Range(0,nseg-1),Range::all())*basis::tri(log2p)->sdiag(mode);
//	}
//	else {
//		for(sind = 0; sind < nseg; ++sind) {
//			Pinv = gbl->sprcn_ut(sind,Range::all(),Range::all());
//			temp = gbl->res.s(sind,mode,Range::all());
//			
//			int info,ipiv[NV];
//			GETRF(NV, NV, Pinv.data(), NV, ipiv, info);
//			
//			if (info != 0) {
//				*gbl->log << "DGETRF FAILED FOR CNS MINVRT" << std::endl;
//				sim::abort(__LINE__,__FILE__,gbl->log);
//			}
//			
//			char trans[] = "T";
//			GETRS(trans,NV,1,Pinv.data(),NV,ipiv,temp.data(),NV,info);
//			
//			if (info != 0) {
//				*gbl->log << "DGETRS FAILED FOR CNS MINVRT" << std::endl;
//				sim::abort(__LINE__,__FILE__,gbl->log);
//			}
//			
//			gbl->res.s(sind,mode,Range::all()) = temp/basis::tri(log2p)->sdiag(mode);
//		}
//	}
//	
//	sc0load(gbl->res.s.data(),mode,mode,gbl->res.s.extent(secondDim));
//	smsgpass(boundary::all,0,boundary::symmetric);
//	sc0wait_rcv(gbl->res.s.data(),mode,mode,gbl->res.s.extent(secondDim));
//	
//	gbl->res_temp.s = gbl->res.s;
//
////	for(i=0;i<nebd;++i)
////		hp_ebdry(i)->sdirichlet(mode);
//
//	for(i=0;i<nebd;++i)
//		hp_ebdry(i)->modify_edge_residual(mode);
//    
//	if (basis::tri(log2p)->im() == 0) return;
//	
//    /* SOLVE FOR INTERIOR MODES */
//    for(tind = 0; tind < ntri; ++tind) {
//		DPBTRSNU2((double *) &basis::tri(log2p)->idiag(0,0),basis::tri(log2p)->ibwth()+1,basis::tri(log2p)->im(),basis::tri(log2p)->ibwth(),&(gbl->res.i(tind,0,0)),NV);
//		tri_hp::restouht_bdry(tind);
//		if (gbl->diagonal_preconditioner) {
//			for(k=0;k<basis::tri(log2p)->im();++k) {
//				gbl->res.i(tind,k,Range::all()) /= gbl->tprcn(tind,Range::all());
//				
//				for (i=0;i<basis::tri(log2p)->bm();++i)
//					for(n=0;n<NV;++n) 
//						gbl->res.i(tind,k,n) -= basis::tri(log2p)->bfmi(i,k)*uht(n)(i);
//			}
//		}
//		else {
//			for(k=0;k<basis::tri(log2p)->im();++k) {
//				/* SUBTRACT BOUNDARY MODES (bfmi is multipled by interior inverse matrix so do this after DPBSLN) */
//				for (i=0;i<basis::tri(log2p)->bm();++i) {
//					for(n=0;n<NV;++n) {
//						for(m=0;m<NV;++m) {
//							gbl->res.i(tind,k,n) -= basis::tri(log2p)->bfmi(i,k)*uht(m)(i)*gbl->tprcn_ut(tind,n,m);
//						}
//					}
//				}
//				/* INVERT PRECONDITIONER (ASSUMES LOWER TRIANGULAR) */
//				//DGETLS(&gbl->tprcn_ut(tind,0,0), NV, NV, &gbl->res.i(tind,k,0));
//				
//				Pinv = gbl->tprcn_ut(tind,Range::all(),Range::all());
//				temp = gbl->res.i(tind,k,Range::all());
//				
//				int info,ipiv[NV];
//				GETRF(NV, NV, Pinv.data(), NV, ipiv, info);
//				
//				if (info != 0) {
//					*gbl->log << "DGETRF FAILED FOR CNS MINVRT" << std::endl;
//					sim::abort(__LINE__,__FILE__,gbl->log);
//				}
//				
//				char trans[] = "T";
//				GETRS(trans,NV,1,Pinv.data(),NV,ipiv,temp.data(),NV,info);
//				
//				if (info != 0) {
//					*gbl->log << "DGETRS FAILED FOR CNS MINVRT" << std::endl;
//					sim::abort(__LINE__,__FILE__,gbl->log);
//				}
//				
//				gbl->res.i(tind,k,Range::all()) = temp;
//			}
//		}
//	}
//
////	if(gbl->diagonal_preconditioner){
////		gbl->res_temp.s = gbl->res.s; //not sure if i need this one
////		gbl->res_temp.i = gbl->res.i;
////		
////		project_res_interior();
////	}
//	
//	
//	return;
//}
