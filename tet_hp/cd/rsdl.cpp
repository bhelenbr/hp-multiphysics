/*
 *  rsdl.cpp
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp_cd.h"
#include "../hp_boundary.h"
	
void tet_hp_cd::rsdl(int stage) {
	int i,j,k,n,tind;
	FLT fluxx,fluxy,fluxz,jcb;
	FLT tres[NV];
	TinyVector<TinyVector<TinyVector<FLT,MXGP>,MXGP>,MXGP> cv00,cv01,cv02,e00,e01,e02;
	TinyVector<TinyVector<FLT,3>,3> visc, d;

	TinyVector<FLT,ND> pt;
	TinyVector<int,4> v;
	int lgpx = basis::tet(log2p).gpx, lgpy = basis::tet(log2p).gpy, lgpz = basis::tet(log2p).gpz;
	int stridey = MXGP;
	int stridex = MXGP*MXGP; 
	
	tet_hp::rsdl(stage);

	for(tind = 0; tind<ntet;++tind) {
		/* LOAD INDICES OF VERTEX POINTS */
		v = tet(tind).pnt;

		if (tet(tind).info < 0) {
			for(n=0;n<ND;++n)
				basis::tet(log2p).proj(pnts(v(0))(n),pnts(v(1))(n),pnts(v(2))(n),pnts(v(3))(n),&crd(n)(0)(0)(0),stridex,stridey);

			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpy;++j) {
	                for(k=0;k<lgpz;++k) {
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
		

		
		/* CALCULATE MESH VELOCITY */
		for(i=0;i<lgpx;++i) {
			for(j=0;j<lgpy;++j) {
	            for(k=0;k<lgpz;++k) {
					mvel(0)(i)(j)(k) = 0;//gbl->bd(0)*(crd(0)(i)(j)(k) -dxdt(log2p,tind,0)(i)(j)(k));
					mvel(1)(i)(j)(k) = 0;//gbl->bd(0)*(crd(1)(i)(j)(k) -dxdt(log2p,tind,1)(i)(j)(k));
					mvel(2)(i)(j)(k) = 0;//gbl->bd(0)*(crd(2)(i)(j)(k) -dxdt(log2p,tind,2)(i)(j)(k));

				}
			}
		}
		
		ugtouht(tind);
		basis::tet(log2p).proj(&uht(0)(0),&u(0)(0)(0)(0),&du(0,0)(0)(0)(0),&du(0,1)(0)(0)(0),&du(0,2)(0)(0)(0),stridex, stridey);

		for(n=0;n<NV;++n)
			for(i=0;i<basis::tet(log2p).tm;++i)
				lf(n)(i) = 0.0;

		/* CONVECTION */
		for(i=0;i<lgpx;++i) {
			for(j=0;j<lgpy;++j) {
	            for(k=0;k<lgpz;++k) {
				
					d(0)(0) =  dcrd(1)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(2)(1)(i)(j)(k)*dcrd(1)(2)(i)(j)(k);
					d(0)(1) = -dcrd(1)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)+dcrd(2)(0)(i)(j)(k)*dcrd(1)(2)(i)(j)(k);
					d(0)(2) =  dcrd(1)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)-dcrd(2)(0)(i)(j)(k)*dcrd(1)(1)(i)(j)(k);
					d(1)(0) = -dcrd(0)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)+dcrd(2)(1)(i)(j)(k)*dcrd(0)(2)(i)(j)(k);
					d(1)(1) =  dcrd(0)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(2)(0)(i)(j)(k)*dcrd(0)(2)(i)(j)(k);
					d(1)(2) = -dcrd(0)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)+dcrd(2)(0)(i)(j)(k)*dcrd(0)(1)(i)(j)(k);
					d(2)(0) =  dcrd(0)(1)(i)(j)(k)*dcrd(1)(2)(i)(j)(k)-dcrd(1)(1)(i)(j)(k)*dcrd(0)(2)(i)(j)(k);
					d(2)(1) = -dcrd(0)(0)(i)(j)(k)*dcrd(1)(2)(i)(j)(k)+dcrd(1)(0)(i)(j)(k)*dcrd(0)(2)(i)(j)(k);
					d(2)(2) =  dcrd(0)(0)(i)(j)(k)*dcrd(1)(1)(i)(j)(k)-dcrd(1)(0)(i)(j)(k)*dcrd(0)(1)(i)(j)(k);

					fluxx = (gbl->ax-mvel(0)(i)(j)(k))*u(0)(i)(j)(k);
					fluxy = (gbl->ay-mvel(1)(i)(j)(k))*u(0)(i)(j)(k);
					fluxz = (gbl->az-mvel(2)(i)(j)(k))*u(0)(i)(j)(k);
					cv00(i)(j)(k) = d(0)(0)*fluxx+d(1)(0)*fluxy+d(2)(0)*fluxz;
					cv01(i)(j)(k) = d(0)(1)*fluxx+d(1)(1)*fluxy+d(2)(1)*fluxz;
					cv02(i)(j)(k) = d(0)(2)*fluxx+d(1)(2)*fluxy+d(2)(2)*fluxz;
				}
			}
		}
		
		basis::tet(log2p).intgrtrst(&lf(0)(0),&cv00(0)(0)(0),&cv01(0)(0)(0),&cv02(0)(0)(0),stridex,stridey);

		/* ASSEMBLE GLOBAL FORCING (IMAGINARY TERMS) */
		lftog(tind,gbl->res);

		/* NEGATIVE REAL TERMS */
		if (gbl->beta(stage) > 0.0) {
			
			ugtouht(tind,1);
			for(n=0;n<NV;++n)
				basis::tet(log2p).proj(&uht(n)(0),&dugdt(n)(0)(0)(0),stridex,stridey);

			/* TIME DERIVATIVE TERMS */
			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpy;++j) {
					for(k=0;k<lgpz;++k) {
						pt(0) = crd(0)(i)(j)(k);
						pt(1) = crd(1)(i)(j)(k);
						pt(2) = crd(2)(i)(j)(k);
						cjcb(i)(j)(k) = dcrd(0)(0)(i)(j)(k)*(dcrd(1)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(1)(i)(j)(k))-dcrd(0)(1)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(0)(i)(j)(k))+dcrd(0)(2)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)-dcrd(1)(1)(i)(j)(k)*dcrd(2)(0)(i)(j)(k));
						res(0)(i)(j)(k) = gbl->bd(0)*cjcb(i)(j)(k)*(u(0)(i)(j)(k)-dugdt(log2p,tind,0)(i)(j)(k));
						res(0)(i)(j)(k) -= cjcb(i)(j)(k)*gbl->src->f(0,pt,gbl->time);
//						cout << cjcb(i)(j)(k) << endl;
					}
				}
			}   
//			cout << cjcb(0)(0)(0) << endl;   
//			cout << "jcb-cjcb = " << tet(tind).vol/8 - cjcb(0)(0)(0) << endl;      
			basis::tet(log2p).intgrt(&lf(0)(0),&res(0)(0)(0)(0),stridex,stridey);

			/* DIFFUSIVE TERMS  */
			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpy;++j) {
					for(k=0;k<lgpz;++k) {
						jcb = gbl->nu/cjcb(i)(j)(k);
						
//						(ys*zt-zs*yt)*dr+(-yr*zt+zr*yt)*ds+(yr*zs-zr*ys)*dt
//						(-xs*zt+zs*xt)*dr+(xr*zt-zr*xt)*ds+(-xr*zs+zr*xs)*dt
//						(xs*yt-ys*xt)*dr+(-xr*yt+yr*xt)*ds+(xr*ys-yr*xs)*dt
						
						/* DIFFUSION TENSOR (LOTS OF SYMMETRY THOUGH)*/
						/* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
						d(0)(0) =  dcrd(1)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(2)(1)(i)(j)(k)*dcrd(1)(2)(i)(j)(k);
						d(0)(1) = -dcrd(1)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)+dcrd(2)(0)(i)(j)(k)*dcrd(1)(2)(i)(j)(k);
						d(0)(2) =  dcrd(1)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)-dcrd(2)(0)(i)(j)(k)*dcrd(1)(1)(i)(j)(k);
						d(1)(0) = -dcrd(0)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)+dcrd(2)(1)(i)(j)(k)*dcrd(0)(2)(i)(j)(k);
						d(1)(1) =  dcrd(0)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(2)(0)(i)(j)(k)*dcrd(0)(2)(i)(j)(k);
						d(1)(2) = -dcrd(0)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)+dcrd(2)(0)(i)(j)(k)*dcrd(0)(1)(i)(j)(k);
						d(2)(0) =  dcrd(0)(1)(i)(j)(k)*dcrd(1)(2)(i)(j)(k)-dcrd(1)(1)(i)(j)(k)*dcrd(0)(2)(i)(j)(k);
						d(2)(1) = -dcrd(0)(0)(i)(j)(k)*dcrd(1)(2)(i)(j)(k)+dcrd(1)(0)(i)(j)(k)*dcrd(0)(2)(i)(j)(k);
						d(2)(2) =  dcrd(0)(0)(i)(j)(k)*dcrd(1)(1)(i)(j)(k)-dcrd(1)(0)(i)(j)(k)*dcrd(0)(1)(i)(j)(k);						
						
						visc(0)(0) = -jcb*(d(0)(0)*d(0)(0)+d(1)(0)*d(1)(0)+d(2)(0)*d(2)(0));
						visc(1)(1) = -jcb*(d(0)(1)*d(0)(1)+d(1)(1)*d(1)(1)+d(2)(1)*d(2)(1));
						visc(2)(2) = -jcb*(d(0)(2)*d(0)(2)+d(1)(2)*d(1)(2)+d(2)(2)*d(2)(2));
						visc(0)(1) = -jcb*(d(0)(0)*d(0)(1)+d(1)(0)*d(1)(1)+d(2)(0)*d(2)(1));
						visc(0)(2) = -jcb*(d(0)(0)*d(0)(2)+d(1)(0)*d(1)(2)+d(2)(0)*d(2)(2));
						visc(1)(2) = -jcb*(d(0)(1)*d(0)(2)+d(1)(1)*d(1)(2)+d(2)(1)*d(2)(2));
#define					viscI1II0I  visc(0)(1)
#define					viscI2II0I  visc(0)(2)
#define					viscI2II1I  visc(1)(2)

						e00(i)(j)(k) = visc(0)(0)*du(0,0)(i)(j)(k)+visc(0)(1)*du(0,1)(i)(j)(k)+visc(0)(2)*du(0,2)(i)(j)(k);
						e01(i)(j)(k) = viscI1II0I*du(0,0)(i)(j)(k)+visc(1)(1)*du(0,1)(i)(j)(k)+visc(1)(2)*du(0,2)(i)(j)(k);
						e02(i)(j)(k) = viscI2II0I*du(0,0)(i)(j)(k)+viscI2II1I*du(0,1)(i)(j)(k)+visc(2)(2)*du(0,2)(i)(j)(k);
						
						cv00(i)(j)(k) += e00(i)(j)(k);
						cv01(i)(j)(k) += e01(i)(j)(k);
						cv02(i)(j)(k) += e02(i)(j)(k);
						
						//*gbl->log << "es " <<  e00(i)(j)(k) << ' ' << e01(i)(j)(k) << ' ' <<  e02(i)(j)(k) << std::endl;
						//*gbl->log << "vs " << visc(0)(0) << ' ' << visc(1)(1) << ' ' << visc(2)(2) << std::endl;
					}
				}
			}
			basis::tet(log2p).derivr(&cv00(0)(0)(0),&res(0)(0)(0)(0),stridex,stridey);
			basis::tet(log2p).derivs(&cv01(0)(0)(0),&res(0)(0)(0)(0),stridex,stridey);
			basis::tet(log2p).derivt(&cv02(0)(0)(0),&res(0)(0)(0)(0),stridex,stridey);
			//cout << res(0)(0)(0)(1)<< endl;

			/* THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES */
			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpy;++j) {
					for(k=0;k<lgpz;++k) {
						d(0)(0) =  dcrd(1)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(2)(1)(i)(j)(k)*dcrd(1)(2)(i)(j)(k);
						d(0)(1) = -dcrd(1)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)+dcrd(2)(0)(i)(j)(k)*dcrd(1)(2)(i)(j)(k);
						d(0)(2) =  dcrd(1)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)-dcrd(2)(0)(i)(j)(k)*dcrd(1)(1)(i)(j)(k);
						d(1)(0) = -dcrd(0)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)+dcrd(2)(1)(i)(j)(k)*dcrd(0)(2)(i)(j)(k);
						d(1)(1) =  dcrd(0)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(2)(0)(i)(j)(k)*dcrd(0)(2)(i)(j)(k);
						d(1)(2) = -dcrd(0)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)+dcrd(2)(0)(i)(j)(k)*dcrd(0)(1)(i)(j)(k);
						d(2)(0) =  dcrd(0)(1)(i)(j)(k)*dcrd(1)(2)(i)(j)(k)-dcrd(1)(1)(i)(j)(k)*dcrd(0)(2)(i)(j)(k);
						d(2)(1) = -dcrd(0)(0)(i)(j)(k)*dcrd(1)(2)(i)(j)(k)+dcrd(1)(0)(i)(j)(k)*dcrd(0)(2)(i)(j)(k);
						d(2)(2) =  dcrd(0)(0)(i)(j)(k)*dcrd(1)(1)(i)(j)(k)-dcrd(1)(0)(i)(j)(k)*dcrd(0)(1)(i)(j)(k);

						tres[0] = gbl->tau(tind)*res(0)(i)(j)(k);
						e00(i)(j)(k) -= (d(0)(0)*(gbl->ax-mvel(0)(i)(j)(k))+d(1)(0)*(gbl->ay-mvel(1)(i)(j)(k))+d(2)(0)*(gbl->az-mvel(2)(i)(j)(k)))*tres[0];
						e01(i)(j)(k) -= (d(0)(1)*(gbl->ax-mvel(0)(i)(j)(k))+d(1)(1)*(gbl->ay-mvel(1)(i)(j)(k))+d(2)(1)*(gbl->az-mvel(2)(i)(j)(k)))*tres[0];
						e02(i)(j)(k) -= (d(0)(2)*(gbl->ax-mvel(0)(i)(j)(k))+d(1)(2)*(gbl->ay-mvel(1)(i)(j)(k))+d(2)(2)*(gbl->az-mvel(2)(i)(j)(k)))*tres[0];

					}
				}
			}

			basis::tet(log2p).intgrtrst(&lf(0)(0),&e00(0)(0)(0),&e01(0)(0)(0),&e02(0)(0)(0),stridex,stridey); 

			for(n=0;n<NV;++n)
				for(i=0;i<basis::tet(log2p).tm;++i)
						lf(n)(i) *= gbl->beta(stage);
			
			lftog(tind,gbl->res_r);
		}
	}

	/* ADD IN VISCOUS/DISSIPATIVE FLUX */
	gbl->res.v(Range(0,npnt-1),Range::all()) += gbl->res_r.v(Range(0,npnt-1),Range::all());
	if (basis::tet(log2p).em > 0) {
		gbl->res.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all()) += gbl->res_r.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all());          
		if (basis::tet(log2p).fm > 0) {
			gbl->res.f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all()) += gbl->res_r.f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all());      
			if (basis::tet(log2p).im > 0) {
				gbl->res.i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all()) += gbl->res_r.i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all());      
			}
		}
	}
	
	/*********************************************/
	/* MODIFY RESIDUALS ON COARSER MESHES            */
	/*********************************************/    
	if (coarse_flag) {
	/* CALCULATE DRIVING TERM ON FIRST ENTRY TO COARSE MESH */
		if(isfrst) {
			dres(log2p).v(Range(0,npnt-1),Range::all()) = fadd*gbl->res0.v(Range(0,npnt-1),Range::all()) -gbl->res.v(Range(0,npnt-1),Range::all());
			if (basis::tet(log2p).em) dres(log2p).e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all()) = fadd*gbl->res0.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all()) -gbl->res.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all());      
			if (basis::tet(log2p).fm) dres(log2p).f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all()) = fadd*gbl->res0.f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all()) -gbl->res.f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all());
			if (basis::tet(log2p).im) dres(log2p).i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all()) = fadd*gbl->res0.i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all()) -gbl->res.i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all());
			isfrst = false;
		}
		gbl->res.v(Range(0,npnt-1),Range::all()) += dres(log2p).v(Range(0,npnt-1),Range::all()); 
		if (basis::tet(log2p).em) gbl->res.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all()) += dres(log2p).e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all());
		if (basis::tet(log2p).fm) gbl->res.f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all()) += dres(log2p).f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all());  
		if (basis::tet(log2p).im) gbl->res.i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all()) += dres(log2p).i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all());  
	}
	
	return;
}
