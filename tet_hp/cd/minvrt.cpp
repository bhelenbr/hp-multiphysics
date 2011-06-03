/*
 *  update.cpp
 *  tet_hp
 *
 *  Created by michael brazell on 1/10/11.
 *  Copyright 2011 Clarkson University. All rights reserved.
 *
 */

#include "tet_hp_cd.h"
#include "../hp_boundary.h"
#include <myblas.h>


/* run code explicitly and converge to full mass matrix */

void tet_hp_cd::minvrt() {	
	
	/* uncomment to perform implicit */
	tet_hp::minvrt();return;

	int i,k,m,n,indx,eind,cnt,find;
	int sign, msgn;
	int j,tind,p0,p1,p2;
	int stridey = MXGP;
	int stridex = MXGP*MXGP;
	FLT jcb,a,h,amax,amin,hmax,hmin,havg,maxvres,maxeres;
	FLT dx1,dy1,dx2,dy2,dz1,dz2,cpi,cpj,cpk;
	TinyVector<int,4> v;
	Array<TinyVector<TinyVector<TinyVector<double,MXGP>,MXGP>,MXGP>,1> ug0(NV);
	TinyVector<FLT,3> pt;
	
	gbl->vprcn(Range::all(),Range::all()) = 0.0;
	gbl->eprcn(Range::all(),Range::all()) = 0.0;
	gbl->fprcn(Range::all(),Range::all()) = 0.0;
	
	FLT timestepinv = 12000.0;
	for(tind=0;tind<ntet;++tind){      
		gbl->iprcn(tind,0) = 0.125*timestepinv*tet(tind).vol; 	
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
	
//	hmax = 0;
//	hmin = 1000000;
//	havg = 0.0;
//	for(tind = 0; tind < ntet; ++tind) {
//		jcb = 0.125*tet(tind).vol; 
//		v = tet(tind).pnt;
//		amax = 0.0;
//		amin = 1000000;
//		for(j=0;j<4;++j) { // FIND MAX FACE AREA AND THEN DIVIDE VOLUME BY IT 
//			find = tet(tind).tri(j);
//			p0 = tri(find).pnt(0);
//			p1 = tri(find).pnt(1);
//			p2 = tri(find).pnt(2);
//			
//			dx1 = pnts(p0)(0)-pnts(p1)(0);
//			dy1 = pnts(p0)(1)-pnts(p1)(1);
//			dz1 = pnts(p0)(2)-pnts(p1)(2);
//			dx2 = pnts(p0)(0)-pnts(p2)(0);
//			dy2 = pnts(p0)(1)-pnts(p2)(1);
//			dz2 = pnts(p0)(2)-pnts(p2)(2);
//			cpi = dy1*dz2-dz1*dy2;
//			cpj = -dx1*dz2+dz1*dx2;
//			cpk = dx1*dy2-dy1*dx2;
//			a =	.5*sqrt(cpi*cpi+cpj*cpj+cpk*cpk);
//			amax = (a > amax ? a : amax);
//			amin = (a < amin ? a : amin);
//		}
//		
//		havg += 4.0*jcb/amax;
//		h = 4.0*jcb/(0.25*(basis::tet(log2p).p+1)*(basis::tet(log2p).p+1)*amax); // 3*8/6=4
//		
//		if(4.0*jcb/amin > hmax)
//			hmax = 4.0*jcb/amin;
//		
//		if(4.0*jcb/amax < hmin)
//			hmin = 4.0*jcb/amax;
//		
//		//cout << "hmin = " << hmin << "  hmax = " << hmax << endl;
//		
//	}
//	
//	cout << "hmin   hmax   havg " << endl;
//	cout << hmin << ' ' << hmax << ' ' << havg/ntet << endl;
	
	
#ifndef NODAL	
	if(basis::tet(log2p).p > 2)
		spoke();
#endif
	
	/* preinvert the jacobian and timestep */
	tet_hp::setup_preconditioner();

	/* uncomment to perform p-1 explicit */
	//tet_hp::minvrt();return;
	
	/* initialize ug with p-1 inversion */
	//tet_hp::minvrt();
	//ug.v(Range(0,npnt-1),Range::all()) = gbl->ug0.v(Range(0,npnt-1),Range::all()) - gbl->res.v(Range(0,npnt-1),Range::all());
	
	/* put Ku-f into res_r */
	gbl->res_r.v(Range(0,npnt-1),Range::all()) = gbl->res.v(Range(0,npnt-1),Range::all());
	if (basis::tet(log2p).em > 0) {
		gbl->res_r.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all()) = gbl->res.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all());          
		if (basis::tet(log2p).fm > 0) {
			gbl->res_r.f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all()) = gbl->res.f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all());          
			if (basis::tet(log2p).im > 0) {
				gbl->res_r.i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all()) = gbl->res.i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all());          
			}
		}
	}

	/* iterate to find full mass matrix */
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
				cout << "curvy" << endl;
				
			}
			
			ugtouht(tind);
			for(n=0;n<NV;++n)
				basis::tet(log2p).proj(&uht(n)(0),&u(n)(0)(0)(0),stridex,stridey);
			
			for (i = 0; i < 4; ++i) {
				int indx = tet(tind).pnt(i);
				for(n = 0; n < NV; ++n)
					uht(n)(i) = gbl->ug0.v(indx,n);
			}	
			
			/* EDGES */
			if (basis::tet(log2p).em > 0){
				cnt = 4;
				for(i = 0; i < 6; ++i) {
					eind = tet(tind).seg(i);
					sign = tet(tind).sgn(i);
					msgn = 1;
					for (m = 0; m < basis::tet(log2p).em; ++m) {
						for(n = 0; n < NV; ++n)
							uht(n)(cnt) = msgn*gbl->ug0.e(eind,m,n);      
						msgn *= sign;
						++cnt;
					}
				}
			}
			
			/* FACES */   
			if (basis::tet(log2p).fm > 0) {
				for(int i = 0; i < 4; ++i) { 
					find = tet(tind).tri(i);
					sign = -tet(tind).rot(i);
					msgn = 1;
					indx = 0;
					for(m = 1; m <= basis::tet(log2p).em-1; ++m) {
						for(k = 1; k <= basis::tet(log2p).em-m; ++k) {
							for(n = 0; n < NV; ++n)
								uht(n)(cnt) = msgn*gbl->ug0.f(find,indx,n);
							++cnt; ++indx;
						}
						msgn *= sign;
						indx += em0 -basis::tet(log2p).em;// index shift for p-mg
					}
				}
			}
			
			/* INTERIOR */
			if (basis::tet(log2p).im > 0) {
				indx = 0;
				for(m = 1; m <= basis::tet(log2p).em-1; ++m) {
					for(k = 1; k <= basis::tet(log2p).em-m; ++k) {
						for(i = 1; i <= basis::tet(log2p).em-m-k; ++i){
							for(n = 0; n < NV; ++n){
								uht(n)(cnt) = gbl->ug0.i(tind,indx,n);
							}
							++cnt; ++indx;
						}
						indx += em0 -basis::tet(log2p).em;				
					}
					
				}
			}

			for(n=0;n<NV;++n)
				basis::tet(log2p).proj(&uht(n)(0),&ug0(n)(0)(0)(0),stridex,stridey);
			
			for(i=0;i<basis::tet(log2p).gpx;++i) {
				for(j=0;j<basis::tet(log2p).gpy;++j) {
					for(k=0;k<basis::tet(log2p).gpz;++k) {
						pt(0) = crd(0)(i)(j)(k);
						pt(1) = crd(1)(i)(j)(k);
						pt(2) = crd(2)(i)(j)(k);
						cjcb(i)(j)(k) = dcrd(0)(0)(i)(j)(k)*(dcrd(1)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(1)(i)(j)(k))-dcrd(0)(1)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(0)(i)(j)(k))+dcrd(0)(2)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)-dcrd(1)(1)(i)(j)(k)*dcrd(2)(0)(i)(j)(k));
						for(n=0;n<NV;++n)
							res(n)(i)(j)(k) = (u(n)(i)(j)(k)-ug0(n)(i)(j)(k))*cjcb(i)(j)(k)*timestepinv;
					}
				}
			}
		
			for(n=0;n<NV;++n)
				for(i=0;i<basis::tet(log2p).tm;++i)
					lf(n)(i) = 0.0;
			
			for(n=0;n<NV;++n)
				basis::tet(log2p).intgrt(&lf(n)(0),&res(n)(0)(0)(0),stridex,stridey);
			
			lftog(tind,gbl->res);
			
		}
		
		/* add back in pre integrated Ku-f from residual */
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

		tet_hp::minvrt(); 

		maxvres = 0.0;
		maxeres = 0.0;
		for(int i=0; i < npnt; ++i){
			if(fabs(gbl->res.v(i,0)) > maxvres)
				maxvres = fabs(gbl->res.v(i,0));
		}
		if (basis::tet(log2p).em > 0) {
			for(int i=0; i < nseg; ++i){
				if(fabs(gbl->res.e(i,0,0)) > maxeres)
					maxeres = fabs(gbl->res.e(i,0,0));
			}
		}
		//cout << it << ' ' << maxvres << ' ' << maxeres << endl;

		FLT cfl = 1.0;
		/* approximate Inversion finished */
		ug.v(Range(0,npnt-1),Range::all()) -= cfl*gbl->res.v(Range(0,npnt-1),Range::all());
		if (basis::tet(log2p).em > 0) {
			ug.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all()) -= cfl*gbl->res.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all());
			if (basis::tet(log2p).fm > 0) {
				ug.f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all()) -= cfl*gbl->res.f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all());
				if (basis::tet(log2p).im > 0) {
					ug.i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all()) -= cfl*gbl->res.i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all());
				}
			}
		}
		
		if (maxvres < 1e-10) {
			cout << it << endl;
			break;
		}

	}
	
	//cout << maxvres << endl;
	
	/* put change of u back into res so that update adds the change to u */
	gbl->res.v(Range(0,npnt-1),Range::all()) = gbl->ug0.v(Range(0,npnt-1),Range::all()) - ug.v(Range(0,npnt-1),Range::all());
	if (basis::tet(log2p).em > 0) {
		gbl->res.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all()) = gbl->ug0.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all())-ug.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all());
		if (basis::tet(log2p).fm > 0) {
			gbl->res.f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all()) = gbl->ug0.f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all())-ug.f(Range(0,ntri-1),Range(0,basis::tet(log2p).fm-1),Range::all());
			if (basis::tet(log2p).im > 0) {
				gbl->res.i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all()) = gbl->ug0.i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all())-ug.i(Range(0,ntet-1),Range(0,basis::tet(log2p).im-1),Range::all());
			}
		}
	}
	
	
	return;
}
