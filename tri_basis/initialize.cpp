/*
 *	initialize.cpp
 *	planar++
 *
 *	Created by helenbrk on Mon Oct 01 2001.
 *	Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */
#define DEBUG

#include <math.h>
#include <utilities.h>
#include <myblas.h>
#include "tri_basis.h"

#ifdef DEBUG
#include <blitz/tinyvec-et.h>
FLT func(FLT r,FLT s) {return(s*s*s*s);}
#endif

template<int _p, int ep> void tri_basis<_p,ep>::initialize() {	

	/*****************************/
	/* SETUP VALUES OF FUNCTIONS */
	/*****************************/
	initialize_values(); // SET UP THINGS FOR PROJECT/INTEGRATE/DERIV
	sideinfoinit(); // SET UP THINGS TO CALCULATE NORMAL DERIVATIVE TO SIDE ALONG SIDE
	lumpinv(); // SET UP THINGS FOR INVERSE OF LUMPED MASS MATRIX
	legpt(); // SET UP PROJECTION TO LEGENDRE POINTS (FOR OUTPUTING)
	
#ifdef DEBUG
	/* SOME TESTING */
	Array<double,1> uht(_tm);
	Array<double,2> u(_gpx,_gpn);
	Array<double,2> u1(_gpx,_gpn);
	Array<double,2> dx(_gpx,_gpn);
	Array<double,2> dy(_gpx,_gpn);
	Array<double,2> dx1(_gpx,_gpn);
	Array<double,2> dy1(_gpx,_gpn);
	
	for(int m=0;m<_tm;++m) {
		for(int j=0;j<_tm;++j)
			uht(j) = 0.0;
		uht(m) = 1.0;
		
		proj_side(0,&uht(0),&u(0,0),&dx(0,0),&dy(0,0));
		for (int i = 0; i<_gpx;++i)
			printf("T0A: %d %e %e %e\n",i,u(0,i),dx(0,i),dy(0,i));
		
		proj_side(1,&uht(0),&u(0,0),&dx(0,0),&dy(0,0));
		for (int i = 0; i<_gpx;++i)
			printf("T0B: %d %e %e %e\n",i,u(0,i),dx(0,i),dy(0,i));
		
		proj_side(2,&uht(0),&u(0,0),&dx(0,0),&dy(0,0));
		for (int i = 0; i<_gpx;++i)
			printf("T0C: %d %e %e %e\n",i,u(0,i),dx(0,i),dy(0,i));

		proj(&uht(0),&u(0,0),&dx(0,0),&dy(0,0),_gpn);
		for(int i=0;i<_gpx;++i) {
			for(int j=0;j<_gpn;++j) {
				dx1(i,j) = 0.0;
				dy1(i,j) = 0.0;
			}
		}
		derivr(&u(0,0),&dx1(0,0),_gpn);
		derivs(&u(0,0),&dy1(0,0),_gpn);
		for(int i=0;i<_gpx;++i)
			for(int j=0;j<_gpn;++j)
				printf("T1: %d %d %e %e %e %e\n",i,j,dx(i,j),dx1(i,j),dy(i,j),dy1(i,j));
				

		if (m < 3) {
			proj(uht(0),uht(1),uht(2),&u1(0,0),_gpn);
			for(int i=0;i<_gpx;++i)
				for(int j=0;j<_gpn;++j)
					printf("T2: %d %d %e %e\n",i,j,u1(i,j),u(i,j));
		}
		
		if (m < _bm) {
			proj_bdry(&uht(0),&u1(0,0),&dx1(0,0), &dy1(0,0),_gpn);
			for(int i=0;i<_gpx;++i)
				for(int j=0;j<_gpn;++j)
					printf("T3: %d %d %e %e %e %e %e %e\n",i,j,u1(i,j),u(i,j),dx(i,j),dx1(i,j),dy(i,j),dy1(i,j));
			
			proj_bdry(&uht(0),&u1(0,0),_gpn);
			for(int i=0;i<_gpx;++i)
				for(int j=0;j<_gpn;++j)
					printf("T4: %d %d %e %e\n",i,j,u1(i,j),u(i,j));
		}

		FLT val,valx,valy,val1,val2;
		ptprobe(1, &val, 0.25, 0.25, &uht(0), _tm);
		ptprobe_bdry(1,&val1, 0.25, 0.25,&uht(0), _tm);
		ptprobe_bdry(1, &val2, &valx, &valy, 0.25, 0.25,&uht(0),_tm);
		printf("T5: %e %e %e\n",val,val1,val2);
		printf("T5: %e %e\n",valx,valy);
		
		intgrtrs(&uht(0),&dx(0,0),&dy(0,0),_gpn);
		for(int j=0;j<_tm;++j)
			printf("T6: %d %e\n",j,uht(j));
		
			
	}


	/* 1D TESTING */
	for (int i=0;i<_gpx;++i)
		printf("_dltx: %d %e\n",i,_dltx(i));
		
	for (int i=0;i<_gpx;++i) 
		for(int n=0;n<_gpx;++n) 
			printf("_dltx1: %d %d %e\n",i,n,_dltx1(i,n));
			
			
	for(int m=0;m<_sm+2;++m) {
		for(int j=0;j<_sm+2;++j)
			uht(j) = 0.0;
		uht(m) = 1.0;
		
		proj1d(&uht(0),&u(0,0),&dx(0,0));
		for(int j=0;j<_gpx;++j)
			dx1(0,j) = 0.0;
		derivx1d(&u(0,0),&dx1(0,0));
		for(int i=0;i<_gpx;++i)
			printf("T11D: %d %e %e\n",i,dx(0,i),dx1(0,i));
		
		proj1d(&uht(0),&u1(0,0));
		for(int i=0;i<_gpx;++i)
			printf("T21D: %d %e %e\n",i,u1(0,i),u(0,i));
			
		if (m < 2) {
			proj1d(uht(0),uht(1),&u1(0,0));
			for(int i=0;i<_gpx;++i)
				printf("T31D: %d %e %e\n",i,u1(0,i),u(0,i));
		}
		
		FLT val,valx,val1;
		
		ptprobe1d(1,&val, 0.25,&uht(0),_tm);
		ptprobe1d(1,&val1,&valx, 0.25, &uht(0),_tm);
		printf("T41D: %e %e %e\n",val,val1,valx);
		
		
		intgrt1d(&uht(0),&u(0,0));
		for(int j=0;j<_sm+2;++j)
			printf("T51D: %d %e\n",j,uht(j));
			
		intgrtx1d(&uht(0),&u(0,0));
		for(int j=0;j<_sm+2;++j)
			printf("T61D: %d %e\n",j,uht(j));
	}
		
#endif
	
	return;
}

template<int _p, int ep> void tri_basis<_p,ep>::initialize_values(void)
{

	FLT al,be,x,eta;
	int i,j,k,m,n,ipoly,ierr;
	Array<FLT,1> e(_tm),e1(_tm),e2(_tm);

	/* ALLOCATE STORAGE FOR RECURSION RELATION COEFFICENTS*/
	ipoly = MAX(_gpx+1,_sm+1);
	ipoly = MAX(ipoly,_gpn+1);

	/* GENERATE RECURSION RELATION FOR LEGENDRE
	POLYNOMIALS (USED TO GENERATE GAUSS POINTS)
	IPOLY = 1 SIGNIFIES LEGENDRE POLYNOMIALS
	NEED RECURSION RELATIONS FROM 1 TO N FOR
	N POINT GAUSS FORMULA

	RECURSION RELATION IS OF THE FORM:
	_p(k+1)(x)=(x-a(k))*_p(k)(x)-b(k)*_p(k-1)(x),
	k=0,1,...,n-1,
	
	_p(-1)(x)=0,	 _p(0)(x)=1
	*/

	/******************************************/
	/* CALCULATE GAUSS POINTS / COEFFICIENTS  */
	/**************************************	   */
	ipoly = 1;
	al = 0.0;
	be = 0.0;
	ierr = recur(_gpx+1,ipoly,al,be,&_a0(0,0),&_b0(0,0));
	if (ierr != 0) {
		printf("recur #1 error %d\n",ierr);
		exit(1);
	}

	ierr = gauss(_gpx,&_a0(0,0),&_b0(0,0),EPSILON,&_x0(0),&_wtx(0),&e(0));
	if (ierr != 0) {
		printf("gauss #1 error %d\n",ierr);
		exit(1);
	}

	/***********************************************/
	/* CALCULATE FACTORS FOR LAGRANGIAN DERIVATIVE */
	/***********************************************/
	for (i=0;i<_gpx;++i) {
		_dltx(i) = 1.0;
		for(k=0;k<i;++k)
			_dltx(i) *= (_x0(i)-_x0(k));
		for(k=i+1;k<_gpx;++k)
			_dltx(i) *= (_x0(i)-_x0(k));
	}
	
	for (i=0;i<_gpx;++i) {
		for(n=0;n<i;++n) 
			_dltx1(i,n) = _dltx(i)/(_x0(i)-_x0(n));
		for(n=i+1;n<_gpx;++n) 
			_dltx1(i,n) = _dltx(i)/(_x0(i)-_x0(n));
	}

	for (i=0;i<_gpx;++i) {
		for(n=0;n<i;++n) 
			_dltn1(i,n) = _dltx(i)/(_x0(i)-_x0(n))*.5*(1+_x0(i));
		for(n=i+1;n<_gpx;++n) 
			_dltn1(i,n) = _dltx(i)/(_x0(i)-_x0(n))*.5*(1+_x0(i));
	}
		
	for (i=0;i<_gpx;++i)
		_dltx(i) = 1.0/_dltx(i);	  

	/*************************************************/ 
	/* NOW CALCULATE VALUES OF G, G' AT GAUSS POINTS */
	/*************************************************/

	/* SIDE 1 IPOLY = 6 FOR JACOBI POLYNOMIALS	*/
	ipoly = 6;
#ifdef MORTHOGONAL
	al = 2.0;
	be = 2.0;
#else
	al = 1.0;
	be = 1.0;
#endif
	ierr = recur(_sm+1,ipoly,al,be,&_a0(0,0),&_b0(0,0));
	if (ierr != 0) {
		printf("recur #3 error %d\n",ierr);
		exit(1);
	}

	for(i=0;i<_tm;++i)
		_norm(i) = 1.0;
	
	for(i = 0;i < _gpx; ++i) {
		x = _x0(i);
		_xp(i) = x;
		_x0(i) = 0.5*(1+x);
		
		ptvalues_deriv(x,0.0);

		for (m = 0;m < _nmodx;++m) {
			_gx(i,m) = pgx(m);
			_dgx(i,m) = dpgx(m);
		}
	}

	/*********************************************/
	/* CALCULATE GAUSS POINTS / COEFFICIENTS (ETA=N) */
	/****************************************	 */
	ipoly = 6;
	al = 1.0;
	be = 0.0;
#ifdef OLDWAY
	ierr = recur(_gpn,ipoly,al,be,&_a0(1,0),&_b0(1,0));
	if (ierr != 0) {
		printf("recur #2 error %d\n",ierr);
		exit(1);
	}	 
	
	ierr = radau(_gpn-1,&_a0(1,0),&_b0(1,0),-1.0,&_n0(0),&_wtn(0),&e(0),&e1(0),&e2(0));		  
	if (ierr != 0) {
		printf("gauss #3 error %d\n",ierr);
		exit(1);
	}
#else
	ierr = recur(_gpn+1,ipoly,al,be,&_a0(1,0),&_b0(1,0));
	if (ierr != 0) {
		printf("recur #2 error %d\n",ierr);
		exit(1);
	}

	ierr = gauss(_gpn,&_a0(1,0),&_b0(1,0),EPSILON,&_n0(0),&_wtn(0),&e(0));
	if (ierr != 0) {
		printf("gauss #1 error %d\n",ierr);
		exit(1);
	}
#endif
	
	for (i=0;i<_gpn;++i)
		_wtn(i) *= 0.5;

	/***********************************************/
	/* CALCULATE FACTORS FOR LAGRANGIAN DERIVATIVE */
	/***********************************************/
	for (i=0;i<_gpn;++i) {
		_dltn(i) = 1.0;
		for(k=0;k<i;++k)
			_dltn(i) *= (_n0(i)-_n0(k));
		for(k=i+1;k<_gpn;++k)
			_dltn(i) *= (_n0(i)-_n0(k));
	}
	
	for (j=0;j<_gpn;++j) {
		for(n=0;n<j;++n) 
			_dltn2(j,n) = _dltn(j)/(_n0(j)-_n0(n));
		for(n=j+1;n<_gpn;++n) 
			_dltn2(j,n) = _dltn(j)/(_n0(j)-_n0(n));
	}
	
	for (j=0;j<_gpn;++j)
		_dltn(j) = 1.0/_dltn(j);
	

	/******************************************/
	/* GENERATE JACOBI POLY FOR S DIRECTION */
	/****************************************/
	/* RECURSION RELATION FOR SIDE MODES */
	/* SIDE 2 IPOLY = 6 FOR JACOBI	  */
	ipoly = 6;
#ifdef MORTHOGONAL
	al = 2.0;
	be = 2.0;
#else
	al = 1.0;
	be = 1.0;
#endif
	ierr = recur(_sm+1,ipoly,al,be,&_a0(1,0),&_b0(1,0));
	if (ierr != 0) {
		printf("recur #3 error %d\n",ierr);
		exit(1);
	}

	/*	RECURSION RELATION FOR INTERIOR MODES */
	for(m = 2; m< _sm+1;++m) {		 
	/* CALCULATE RECURSION RELATION FOR _p^(2m-1,1)(s) SIDE 1 IPOLY = 6 FOR JACOBI	   */
		ipoly = 6;
#ifdef MORTHOGONAL
		al = 2.*m+1;
		be = 2.0;
#else
		al = 2.*m-1;
		be = 1.0;
#endif
		ierr = recur(_sm+2-m,ipoly,al,be,&_a0(m,0),&_b0(m,0));
		if (ierr != 0) {
			printf("recur #4 error %d\n",ierr);
			exit(1);
		}
	}
	
	
	for(j=0;j<_gpn;++j) {
		eta = _n0(j);
		_np(j) = eta;
		_n0(j) = 2.0/(1-eta);

		ptvalues_deriv(0.0,eta);
		
		for(m=0;m<_tm;++m) {
			_gn(j,m) = pgn(m);
			_dgn(j,m) = dpgn(m);
		}
	}

	/******************************/
	/* CALCULATE _norm				  */
	/* ************************** */
	/* SIDE & VERTEX MODES */
	for(n=0;n<3;++n)
		_norm(n) = 1.0;
	for (n = 3; n < _sm+3; ++n) {
		_norm(n) = 0.0;
		for(i = 0; i < _gpx; ++i)
			_norm(n) += _wtx(i)*_gx(i,n)*_gx(i,n);
		_norm(n) = 1./sqrt(_norm(n));
		_norm(n+_sm) = _norm(n);
		_norm(n+2*_sm) = _norm(n);
	}
	
	/* INTERIOR MODES */
	for(m = _bm; m < _tm; ++m) {
		_norm(m) = 0.0;
		for(j = 0; j < _gpn; ++j)
			  _norm(m) += _wtn(j)*_gn(j,m)*_gn(j,m);
		  _norm(m) = 1./sqrt(_norm(m));
	}

	/***************/
	/* RENORMALIZE */
	/***************/
	for (n = 3; n < _nmodx; ++n) {			 
		for (i =0;i<_gpx;++i) {
			_gx(i,n) *= _norm(n);
			_dgx(i,n) *= _norm(n);
#ifdef DEBUG
			printf("IV1: %d %d %e %e\n",i,n,_gx(i,n),_dgx(i,n));
#endif
		}
		for (j=0;j<_gpn;++j) {
			/* SIDE 2 & 3 MUST BE RENORMALIZED BY SAME CONSTANT TO MATCH */
			_gn(j,n +_sm) *= _norm(n);
			_dgn(j,n+_sm) *= _norm(n);
			_gn(j,n +2*_sm) *= _norm(n);
			_dgn(j,n +2*_sm) *= _norm(n);
#ifdef DEBUG
			printf("IV2: %d %d %e %e\n",j,n,_gn(j,n),_dgx(j,n));
#endif
		}
	}

	for(m = _bm; m < _tm; ++m) {
		for(j = 0; j < _gpn; ++j) {
			_gn(j,m) = _gn(j,m)*_norm(m);
			_dgn(j,m) = _dgn(j,m)*_norm(m);
#ifdef DEBUG
			printf("IV3: %d %d %e %e\n",j,m,_gn(j,m),_dgn(j,m));
#endif
		}
	}

	/*****************************************/
	/* PRECALCULATE THINGS FOR FAST INTGRTRS */
	/*****************************************/
	for(m=0;m<_nmodx;++m) {
		for(i=0;i<_gpx;++i) {
			_gxwtx(m,i) = _gx(i,m)*_wtx(i);
			_dgxwtx(m,i) = _dgx(i,m)*_wtx(i);
#ifdef DEBUG
			printf("%d %d %e %e\n",m,i,_gxwtx(m,i),_dgxwtx(m,i));
#endif
		}
	}
	
	for(m=0;m<_tm;++m) {
		for(j=0;j<_gpn;++j) {
			_gnwtn(m,j) = _gn(j,m)*_wtn(j);
			_gnwtnn0(m,j) = _gn(j,m)*_wtn(j)*_n0(j);
			_dgnwtn(m,j) = _dgn(j,m)*_wtn(j);
#ifdef DEBUG
			printf("%d %d %e %e %e\n",m,j,_gnwtn(m,j),_gnwtnn0(m,j),_dgnwtn(m,j));
#endif
		}
	}
	
	return;
}
	
	
template<int _p, int ep> void tri_basis<_p,ep>::sideinfoinit() {
	int i,m,n,ind;
	FLT x,eta,xp1oeta,xp1,oeta;
	
	/*	THIS IS TO CALCULATE NORMAL DERIVATIVES TO SIDE */
	/* SIDES 1 & 2 ARE ROTATED TO SIDE 0 POSITION */
	_dgnorm.resize(3,_tm,_gpx);
		
	if (_p == 0) {
		for(ind=0;ind<3;++ind)
			for(i=0;i<_gpx;++i)
				_dgnorm(ind,0,i) = 0.0;
		
		return;
	}
	
	/*	CALCULATE NORMAL DERIVATIVE VALUES ALONG SIDES */
	/*	SIDE 0 */
	for(i=0;i<_gpx;++i) {
		eta = -1.0;
		x = 2.*_x0(i) -1.0;
		ptvalues_deriv(x,eta);
		xp1oeta = _x0(i)*2.0/(1-eta);
		
		/* CALCULATE POLYNOMIALS */
		/* VERTEX 0	   */
		_dgnorm(0,0,i) = dpgn(0)*pgx(0) +xp1oeta*pgn(0)*dpgx(0);

		/* VERTEX 1	 */
		_dgnorm(0,1,i) = dpgn(1)*pgx(1) +xp1oeta*pgn(1)*dpgx(1);

		/* VERTEX 2		*/	  
		_dgnorm(0,2,i) = dpgn(2)*pgx(2) +xp1oeta*pgn(2)*dpgx(2);

		for(m = 3; m < _sm+3; ++m)
			_dgnorm(0,m,i) = dpgn(m)*pgx(m) +xp1oeta*pgn(m)*dpgx(m);
			
		for(m=_sm+3;m<2*_sm+3;++m)
			_dgnorm(0,m,i) = dpgn(m)*pgx(2) +xp1oeta*pgn(m)*dpgx(2);

		for(m=2*_sm+3;m<_bm;++m)
			_dgnorm(0,m,i) = dpgn(m)*pgx(1) +xp1oeta*pgn(m)*dpgx(1);

		/*	INTERIOR MODES	  */
		ind = _bm;
		for(m = 3; m< _sm+2;++m) {		 
			for(n = 1; n < _sm+3-m;++n) {
				_dgnorm(0,ind,i) = dpgn(ind)*pgx(m) +xp1oeta*pgn(ind)*dpgx(m);
				++ind;
			}
		}
	}
	
	/* SIDE 1 */
	for(i=0;i<_gpx;++i) {
		eta = _xp(i);
		x = 1.0;
		ptvalues_deriv(x,eta);
		oeta = 2.0/(1-eta);

		/* CALCULATE POLYNOMIALS */
		/* VERTEX 0	   */
		_dgnorm(1,0,i) = -oeta*pgn(0)*dpgx(0);

		/* VERTEX 1	 */
		_dgnorm(1,1,i) = -oeta*pgn(1)*dpgx(1);

		/* VERTEX 2		*/	  
		_dgnorm(1,2,i) = -oeta*pgn(2)*dpgx(2);

		for(m = 3; m < _sm+3; ++m)
			_dgnorm(1,m,i) = -oeta*pgn(m)*dpgx(m);
			
		for(m=_sm+3;m<2*_sm+3;++m)
			_dgnorm(1,m,i) = -oeta*pgn(m)*dpgx(2);

		for(m=2*_sm+3;m<_bm;++m)
			_dgnorm(1,m,i) = -oeta*pgn(m)*dpgx(1);

		/*	INTERIOR MODES	  */
		ind = _bm;
		for(m = 3; m< _sm+2;++m) {		 
			for(n = 1; n < _sm+3-m;++n) {
				_dgnorm(1,ind,i) = -oeta*pgn(ind)*dpgx(m);
				++ind;
			}
		}		 
	}
	
	/* SIDE 2 */
	for(i=0;i<_gpx;++i) {
		x = -1.0;
		eta = 1.0 -2.*_x0(i);
		ptvalues_deriv(x,eta);
		oeta = 2.0/(1-eta);
		xp1 = (x+1)/2.0;

		/* CALCULATE POLYNOMIALS */
		/* VERTEX 0	   */
		_dgnorm(2,0,i) = -(dpgn(0)*pgx(0) +(xp1 -1)*oeta*pgn(0)*dpgx(0));

		/* VERTEX 1	 */
		_dgnorm(2,1,i) = -(dpgn(1)*pgx(1) +(xp1 -1)*oeta*pgn(1)*dpgx(1));

		/* VERTEX 2		*/	  
		_dgnorm(2,2,i) = -(dpgn(2)*pgx(2) +(xp1 -1)*oeta*pgn(2)*dpgx(2));

		for(m = 3; m < _sm+3; ++m)
			_dgnorm(2,m,i) = -(dpgn(m)*pgx(m) +(xp1 -1)*oeta*pgn(m)*dpgx(m));
			
		for(m=_sm+3;m<2*_sm+3;++m)
			_dgnorm(2,m,i) = -(dpgn(m)*pgx(2) +(xp1 -1)*oeta*pgn(m)*dpgx(2));

		for(m=2*_sm+3;m<_bm;++m)
			_dgnorm(2,m,i) = -(dpgn(m)*pgx(1) +(xp1 -1)*oeta*pgn(m)*dpgx(1));

		/*	INTERIOR MODES	  */
		ind = _bm;
		for(m = 3; m< _sm+2;++m) {		 
			for(n = 1; n < _sm+3-m;++n) {
				_dgnorm(2,ind,i) = -(dpgn(ind)*pgx(m) +(xp1 -1)*oeta*pgn(ind)*dpgx(m));
				++ind;
			}
		}
	}

	return;
}

/************************************************/
/** CALCULATE THINGS FOR LUMPED MASS INVERSION	*/
/************************************************/
template<int _p, int ep> void tri_basis<_p,ep>::lumpinv(void) {
	int i,i1,i2,j,k,m,info,ind,ind1,n;
	Array<int,1> ipiv(2*_tm);
	Array<FLT,2> mwk(_tm,_tm),mm(_tm,_tm);
	Array<FLT,1> u(_tm),l(_tm),vwk(_tm),wk1(_tm*_gpx);
	FLT rcond=1;
	char trans[] = "T", uplo[] = "U";
	
	/* ALLOCATE MASS MATRIX INVERSION VARIABLES */
	if (_sm > 0) {
		_vfms.resize(3,_sm);
		_sfmv.resize(2,_sm);
		_sdiag.resize(_sm);
	}
	if (_sm > 1) _sfms.resize(_sm-1,_sm,3);
	if (_im > 0) {
		_ifmb.resize(_bm,_im);
		_bfmi.resize(_bm,_im);
		_idiag.resize(_im,_ibwth+1);
	}
	_msi.resize(_bm,_bm);
	
	/* ALLOCATE 1D MASS MATRIX INVERSION VARIABLES */
	if (_sm > 0) {
		_vfms1d.resize(2,_sm);
		_sfmv1d.resize(2,_sm);
		_sdiag1d.resize(_sm,_sbwth+1);
	}

	/********************************************/	  
	/* GENERATE MASS MATRIX						   */
	/********************************************/
	for(m=0;m<_tm;++m) {
		for(i=0;i<_tm;++i)
			u(i) = 0.0;
		u(m) = 1.0;

		proj(&u(0),&wk1(0),_gpn);  // PROJECT USES WK0
		intgrt(&l(0),&wk1(0),_gpn); // INTGRT USES WK0

#ifdef DEBUG
		for(i=0;i<_gpx;++i)
			for(j=0;j<_gpn;++j)
				printf("MIMM: %d %d %e\n",i,j,wk1(i*_gpn +j));
		for(i=0;i<_tm;++i)
			printf("MIMM2: %d %e\n",i,l(i));
#endif
				
		for(i=0;i<_tm;++i) 
			mm(m,i) = l(i);		   
	}

	/*******************************************************/		 
	/*	EQUATIONS TO FIND VERTEX VALUES TO _sm-1 ACCURACY */
	/*******************************************************/		 
	if (_p > 1) {		
		for(i=0;i<3;++i) {
			i1 = (i+1)%3;
			i2 = (i+2)%3;
			
			/*2 CROSS VERTEX CONSTRAINTS */
			for(k=0;k<2;++k) {
				ind = (i+1+k)%3;
				vwk(k) = mm(ind,i);
				for(j=0;j<_sm;++j) {
					mwk(k,j) = mm(ind,3+j+i1*_sm);
					mwk(k,j+_sm) = mm(ind,3+j+i2*_sm);
				}
				for(j=0;j<_im;++j)
					mwk(k,2*_sm+j) = mm(ind,j+_bm);
			}

			/*3x(_sm-3) SIDE CONSTRAINTS */
			for(k=0;k<3;++k) {
				for(m=0;m<_sm-1;++m) {	 
					vwk(m+2+k*(_sm-1)) = mm(3+m+k*_sm,i);
					for(j=0;j<_sm;++j) {
						mwk(m+2+k*(_sm-1),j) = mm(3+m+k*_sm,3+j+i1*_sm);
						mwk(m+2+k*(_sm-1),j+_sm) = mm(3+m+k*_sm,3+j+i2*_sm);
					}
					for(j=0;j<_im;++j)
						mwk(m+2+k*(_sm-1),2*_sm+j) = mm(3+m+k*_sm,j+_bm);
				}
			}
			
			
			/*(_sm-3)*(_sm-4) INTERNAL MODE CONSTRAINTS  */				 
			ind = 3*_sm-1;	 ind1 = 0;
			for(m=2;m<_sm+1;++m) {		 
				for(n = 1; n < _sm+1-m;++n) {
					vwk(ind) = mm(ind1+_bm,i);
					for(j=0;j<_sm;++j) {
						mwk(ind,j) = mm(ind1+_bm,3+j+i1*_sm);
						mwk(ind,j+_sm) = mm(ind1+_bm,3+j+i2*_sm);
					}
					for(j=0;j<_im;++j)
						mwk(ind,2*_sm+j) = mm(ind1+_bm,j+_bm);
					++ind;
					++ind1;
				}
				++ind1;
			}

			GETRF((_sm+1)*(_sm+2)/2-1,(_sm+1)*(_sm+2)/2-1,&mwk(0,0),_tm,&ipiv(0),info);
			if (info != 0) {
				printf("DGETRF FAILED - VRTX info:%d _sm:%d i:%d\n",info,(_sm+2),i);
				exit(1);
			}
			GETRS(trans,(_sm+1)*(_sm+2)/2-1,1,&mwk(0,0),_tm,&ipiv(0),&vwk(0),_tm,info);
													
			/* STORE INTERIOR VALUES */
			for(k=0;k<_im;++k)
				_ifmb(i,k) = vwk(k+2*_sm);			
		}
	
		/* STORE SIDE VALUES */
		for(k=0;k<_sm;++k) {
			_sfmv(0,k) = vwk(k+_sm);
			_sfmv(1,k) = vwk(k);
			
		}

		/* FIND VERTEX DIAGANOL ELEMENT */
		_vdiag = mm(2,2);
		for(k=0;k<_sm;++k) {
			_vdiag -= _sfmv(0,k)*mm(2,3+k+1*_sm);
			_vdiag -= _sfmv(1,k)*mm(2,3+k);
		}

		for(k=0;k<_im;++k)
			_vdiag -= _ifmb(2,k)*mm(2,k+_bm);
	}
	else {
		if (_p > 0)
			_vdiag = mm(0,0) + mm(0,1) + mm(0,2);
		else
			_vdiag = mm(0,0);
	}


	/*****************************************************************/
	/* EQUATIONS TO FIND STATIC INVERSION OF INTERIORS FROM SIDES */
	/*****************************************************************/
	/** WARNING THIS ONLY WORKS UP TO _sm = 5!!!! ***/
	if (_im > 0) {
		for(i=0;i<3;++i) {
			i1 = (i+1)%3;
			i2 = (i+2)%3;

			for(k=0;k<(_sm+2)-3;++k) {
				ind = 0;
				/* CROSS SIDE MODES */
				for(j=k;j<_sm-1;++j) {
					vwk(ind) =	mm(i*_sm+3+k,i1*_sm+3+j);
					for(m=0;m<_im;++m)  
						mwk(ind,m) = mm(_bm+m,i1*_sm+3+j);
					++ind;
				}
				
				for(j=k;j<MIN(2*k,_sm-1);++j) {
					vwk(ind) =	mm(i*_sm+3+k,i2*_sm+3+j);
					for(m=0;m<_im;++m)  
						mwk(ind,m) = mm(_bm+m,i2*_sm+3+j);
					++ind;
				}							 

				/* INTERIOR MODES */	
				ind1 = 0;
				for(m = 2; m< _sm+1;++m) {		 
					for(n = 1; n < _sm+1-m;++n) {
						vwk(ind) = mm(i*_sm+3+k,_bm+ind1);
						for(j=0;j<_im;++j) 
							mwk(ind,j) = mm(_bm+j,_bm+ind1);
						++ind;
						++ind1;
					}
					++ind1;
				}

				GETRF(_im,_im,&mwk(0,0),_tm,&ipiv(0),info);
				if (info != 0) {
					printf("DGETRF FAILED SIDE info:%d (_sm+2):%d k:%d\n",info,(_sm+2),k);
					exit(1);
				}
				else
					GETRS(trans,_im,1,&mwk(0,0),_tm,&ipiv(0),&vwk(0),_tm,info);

				for(j=0;j<_im;++j)
					_ifmb(i*_sm+3+k,j) = vwk(j);
			}

			/* FOR HIGHEST ORDER MODE - STATIC INVERT ALL INTERIOR MODES */
			for(j=0;j<_im;++j) {
				vwk(j) = mm(i*_sm +_sm +2,j+_bm);
				for(k=0;k<_im;++k)
					mwk(j,k) = mm(k+_bm,j+_bm);
			}

			GETRF(_im,_im,&mwk(0,0),_tm,&ipiv(0),info);
			if (info != 0) {
				printf("GETRF FAILED info: %d cond: %f\n",info,rcond);
				exit(1);
			}
			GETRS(trans,_im,1,&mwk(0,0),_tm,&ipiv(0),&vwk(0),_tm,info);

			for(j=0;j<_im;++j)
				_ifmb(i*_sm +_sm +2,j) = vwk(j);
		}
	}

	/*********************************************/
	/* FIND VERTEX-SIDE & SIDE-SIDE MATRICES  */
	/*********************************************/
	for(k=0;k<_sm;++k) {
		for(j=0;j<_tm;++j)
			u(j) = 0.0;
			
		u(k+3) = 1.0;
		for(j=0;j<_im;++j) {
			u(j+_bm) = -_ifmb(k+3,j);
		}

		proj(&u(0),&wk1(0),_gpn);
		intgrt(&l(0),&wk1(0),_gpn);
		
		for(j=0;j<3;++j)
			_vfms(j,k) = l(j);
			
		_sdiag(k) = 1.0/l(3+k);
		
		for(j=0;j<k;++j) {
			_sfms(j,k,0) = l(3+j);
			_sfms(j,k,1) = l(3+j+_sm);
			_sfms(j,k,2) = l(3+j+2*_sm);
		}
	}


	/************************************************/
	/* FIND MATRICES TO DETERMINE INTERIOR MODES */
	/************************************************/
	if (_im > 0) {
		/* SETUP DIAGANOL FORM OF MATRIX */
		/* ONLY NECESSARY WHEN USING DPBTRF */	  
		for(j=0;j<_im;++j) {
			i1 = (0 > j-_ibwth ? 0 : j-_ibwth);
			for(i=i1;i<=j;++i) {
				k = i - (j-_ibwth);
				_idiag(j,k) = mm(i+_bm,j+_bm);
			}
		}
		
		PBTRF(uplo,_im,_ibwth,&_idiag(0,0),_ibwth+1,info);
		if (info != 0) {
			printf("1:PBTRF FAILED info: %d\n", info);
			exit(1);
		}
					
		/* MATRIX TO REMOVE BOUNDARY MODES FROM INTERIOR */
		for(i=0; i<_bm; ++i ) {
			for(j=0; j<_im; ++j ) {
				_bfmi(i,j) = mm(i,j+_bm);
			}
			PBTRS(uplo,_im,_ibwth,1,&_idiag(0,0),_ibwth+1,&_bfmi(i,0),_im,info);
		}
		
		/* TEST INVERSE 
		for(i=0;i<_im;++i) {
			for(j=0;j<_bm;++j)
				printf(" %f ",_bfmi(j,i));	  
			printf("\n");
		}
			
		printf("_ibwth %d\n",_ibwth);
		for(i=0;i<_im;++i) {
			PBTRS(uplo,_im,_ibwth,1,&_idiag(0,0),_ibwth+1,&mm(i+_bm,_bm),_im,info);
			printf("%d: ",i);
			for(j=0;j<_im;++j)
				printf("%f ",mm(i+_bm,_bm+j));
			printf("\n");
		}
		*/
	}
		
	/* FORM A - BC^{-1}B^T */  
	/* MASS MATRIX WITH STATIC CONDENSATION OF INTERIOR MODES */
	for(i=0;i<_bm;++i) {
		for(j=0;j<_bm;++j) {
			mwk(i,j) = mm(i,j);
			for(k=0;k<_im;++k)
				mwk(i,j) -= _bfmi(i,k)*mm(j,k+_bm);
		}
	}

#ifndef SKIP
	/* SETUP BAND DIAGONAL FORM OF MATRIX */
	/* ONLY NECESSARY WHEN USING DPBTRF */ 
	for(j=0;j<_bm;++j) {
		for(i=0;i<=j;++i) {
			k = i - (j-(_bm-1));
			_msi(j,k) = mwk(i,j);
		}
	}
	PBTRF(uplo,_bm,_bm-1,&_msi(0,0),_bm,info);
	if (info != 0) {
		printf("2:PBTRF FAILED info: %d\n", info);
		exit(1);
	}
#else
	for(i=0;i<_bm;++i)
		for(j=0;j<_bm;++j) 
			_msi(i,j) = mwk(i,j);
			
	DPOTRF(uplo,_bm,&_msi(0,0),_bm,info);
	if (info != 0) {
		printf("POTRF FAILED info: %d\n", info);
		exit(1);
	}
#endif


#ifdef DEBUG
	std::cout << "_sfmv" << std::endl;
	std::cout << _sfmv << std::endl;
	std::cout << "_ifmb" << std::endl;
	std::cout << _ifmb << std::endl;
	std::cout << "_vdiag" << std::endl;
	std::cout << _vdiag << std::endl;
	std::cout << "_sdiag" << std::endl;
	std::cout << _sdiag << std::endl;
	std::cout << "_vfms" << std::endl;
	std::cout << _vfms << std::endl;
	std::cout << "_sfms" << std::endl;
	std::cout << _sfms << std::endl;
	std::cout << "_bfmi" << std::endl;
	std::cout << _bfmi << std::endl;
	std::cout << "_idiag" << std::endl;
	std::cout << _idiag << std::endl;

	/* CHECK TO MAKE SURE PREVIOUS RESULTS ARE RIGHT */
	for(i=0;i<3;++i) {
		i1 = (i+1)%3;
		i2 = (i+2)%3;

		u(i) = 1.0;
		u(i1) = 0.0;
		u(i2) = 0.0;
		
		for(j=0;j<_sm;++j) {
			u(3+j+i*_sm) = 0.0;
			u(3+j+i1*_sm) = -_sfmv(1,j);
			u(3+j+i2*_sm) = -_sfmv(0,j);
		}
		
		for(j=0;j<_im;++j)
			u(_bm+j) = -_ifmb(i,j);
			
		proj(&u(0),&wk1(0),_gpn);
		intgrt(&mwk(i,0),&wk1(0),_gpn);
		printf("%2d:",i);
		for(j=0;j<_tm;++j)
			printf("%+.4le	",mwk(i,j));
		printf("\n");
	}
	
	for(i=0;i<3;++i) {
		for(k=0;k<_sm;++k) {
			for(j=0;j<_tm;++j)
				u(j) = 0.0;
			u(i*_sm+k+3) = 1.0;
			for(j=0;j<_im;++j) {
				u(j+_bm) = -_ifmb(i*_sm+k+3,j);
			}

			proj(&u(0),&wk1(0),_gpn);
			intgrt(&mwk(i*_sm+k+3,0),&wk1(0),_gpn);
			
			printf("%2d:",3+i*_sm+k);
			for(j=0;j<_tm;++j)
				printf("%+.4le	",mwk(i*_sm+k+3,j));
			printf("\n");
		}			 
	}
	
	for(i=0;i<_im;++i) {
		for(j=0;j<_tm;++j)
			u(j) = 0.0;
		u(i +_bm) = 1.0;	 
		proj(&u(0),&wk1(0),_gpn);
		intgrt(&mwk(i+_bm,0),&wk1(0),_gpn);
	}
	
	for(i=0;i<3;++i)
		for(m=0;m<_tm;++m)
			mwk(i,m) /= _vdiag;
			
	for(i=0;i<3;++i) {
		for(j=0;j<3;++j) {
			i1 = (i+j)%3;
			for(k=0;k<_sm;++k) {
				for(m=0;m<_tm;++m)
					mwk(i*_sm+k+3,m) -= _vfms(j,k)*mwk(i1,m);
			}
		}
	}
	
	/* REMOVE MODES J,K FROM MODE I,M */
	int mode;
	for(mode = 0; mode <_sm;++mode) {
		for(i=0;i<3;++i)
			for(k=0;k<_tm;++k)
				mwk(3+i*_sm+mode,k) *= _sdiag(mode);
		for(i=0;i<3;++i) {
			for(m=mode+1;m<_sm;++m) {
				for(j=0;j<3;++j) {
					i1 = (i+j)%3;
					for(k=0;k<_tm;++k)
						mwk(3+i*_sm+m,k) -= _sfms(mode,m,j)*mwk(3+i1*_sm+mode,k);
				}
			}
		}
	}
	if (_im) DPBTRSNU1(_idiag.data(),_im,_ibwth,&mwk(_bm,0),_tm);
	for(k=0;k<_im;++k)
		for(i=0;i<_bm;++i)
			for(j=0;j<_tm;++j)
				mwk(_bm+k,j) -= _bfmi(i,k)*mwk(i,j);

	for(m=0;m<_tm;++m) {
		printf("LI1: %2d:",m);
		for(j=0;j<_tm;++j) {
			if (fabs(mwk(m,j)) > 1.0e-12)
				printf("%+.2e  ",mwk(m,j));
			else
				printf("%+.2e  ",0.0);
		}
		printf("\n");
	}
#endif	  

	lumpinv1d();

}


template<int _p, int ep> void tri_basis<_p,ep>::lumpinv1d() {
	int i,i1,j,k,m,info;
	Array<int,1> ipiv(2*_tm);
	Array<FLT,2> mwk(_tm,_tm),mm(_tm,_tm);
	Array<FLT,1> u(_tm),l(_tm),vwk(_tm),wk1(_tm*_gpx);
	FLT rcond=1;
	char trans[] = "T", uplo[] = "U";

	/********************************************************************/	  
	/* NOW SETUP SIMILAR THING FOR 1-D MATRICES										*/
	/********************************************************************/	  
	/* CALCULATE 1-D MASS MATRIX */	   
	for(m=1;m<_p+2;++m) {
		for(k=1;k<_p+2;++k) {
			mm(m-1,k-1)= 0.0;
			for(i=0;i<_gpx;++i)
				mm(m-1,k-1) += _wtx(i)*_gx(i,m)*_gx(i,k);
		}
	}

#ifdef DEBUG
	printf("\n mm MATRIX (_sm+2) = %d\n",(_sm+2));
	for(i=0;i<_sm+2;++i) {
		printf("LI2: %2d:",i);
		for(j=0;j<_sm+2;++j)
			printf("%+.4lf ",mm(i,j));
		printf("\n");
	}
#endif

	if(_p > 1) {
		/* STATIC INVERSION SIDE MATRIX	 */			   
		for(j=0;j<_sm;++j) {
			i1 = (0 > j-_sbwth ? 0 : j-_sbwth);
			for(i=i1;i<=j;++i) {
				k = i - (j-_sbwth);
				_sdiag1d(j,k) = mm(i+2,j+2);
			}
		}

		PBTRF(uplo,_sm,_sbwth,&_sdiag1d(0,0),_sbwth+1,info);
		if (info != 0 || rcond < 100.*EPSILON) {
			printf("PBTRF FAILED - 1D (_sm+2) : %d info: %d cond: %f\n",(_sm+2), info,rcond);
			exit(1);
		}

		/* MATRIX TO REMOVE VERTEX MODES FROM SIDES */
		for(i=0; i<2; ++i) {
			for(j=0; j<_sm; ++j) {
				_vfms1d(i,j) = mm(i,j+2);
			}
			PBTRS(uplo,_sm,_sbwth,1,&_sdiag1d(0,0),_sbwth+1,&_vfms1d(i,0),(_sm+2)-2,info);
		}

		/* MATRIX TO REMOVE SIDE MODES FROM VERTICES */	   
		for(i=0;i<2;++i) {		  
			/* VERTEX CONSTRAINT */
			vwk(0) = mm(i,(i+1)%2);
			for(m=0;m<_sm;++m)
				mwk(0,m) = mm((i+1)%2,m+2);

			/* SIDE CONSTRAINTS */
			for(m=0;m<_sm-1;++m) {
				vwk(m+1) = mm(i,m+2);
				for(k=0;k<_sm;++k) {
					mwk(m+1,k) = mm(m+2,k+2);
				}
			}
			GETRF(_sm,_sm,&mwk(0,0),_tm,&ipiv(0),info);
			if (info != 0) {
				printf("1D DGETRF FAILED - VRTX info:%d (_sm+2):%d i:%d\n",info,(_sm+2),i);
				exit(1);
			}
			GETRS(trans,_sm,1,&mwk(0,0),_tm,&ipiv(0),&vwk(0),_tm,info);

			/* STORE MATRIX */
			for(k=0;k<_sm;++k)
				_sfmv1d(i,k) = vwk(k);
		}

		/* FIND VERTEX DIAGANOL ELEMENT */
		_vdiag1d = mm(0,0);
		for(k=0;k<_sm;++k)
			_vdiag1d -= _sfmv1d(0,k)*mm(0,2+k);
	}
	else {
		if (_p > 0) 
			_vdiag1d = mm(0,0)+mm(0,1);
		else
			_vdiag1d = mm(0,0);
	}

#ifdef DEBUG
	/* CHECK TO MAKE SURE PREVIOUS 1D RESULTS ARE RIGHT */
	Array<FLT,1> uht(_tm);
	for(k=0;k<2;++k) {
		u(k) = 1.0;
		u((k+1)%2) = 0.0;
		for(m=0;m<_sm;++m)
			u(m+2) = -_sfmv1d(k,m);

		for(i=0;i<_gpx;++i) {
			uht(i) = u(0)*_gx(i,0);
			uht(i) += u(1)*_gx(i,1);
			for(m=0;m<_sm;++m)
				uht(i) += u(m+2)*_gx(i,3+m);
		}

		for(m=1;m<_sm+3;++m) {
			l(m) = 0.0;
			for(i=0;i<_gpx;++i) {
				l(m) += _wtx(i)*_gx(i,m)*uht(i);
			}
		}
		printf("LI3 %2d:",k);
		for(j=1;j<_sm+3;++j)
			printf("%+.4le	",l(j));
		printf("%+.4le	",_vdiag1d);
		printf("\n");
	}
#endif



	return;
}			 

template<int _p, int ep> void tri_basis<_p,ep>::legpt()
{
	FLT x,eta,r,s;
	int i,j,m,n,ind;
	
	_lgrnge1d.resize(_nmodx,_sm+1);
	_lgrnge.resize(_tm,_sm+1,_sm+1);
	
//	  /* USING LOBATTO SPACING */
//	  int ipoly = 1;
//	  double al = 0.0;
//	  double be = 0.0;
//	  ierr = recur(_sm+3,ipoly,al,be,&_a0(0),&_b0(0));
//	  if (ierr != 0) {
//		  printf("recur #1 error %d\n",ierr);
//		  exit(1);
//	  }
//	  ierr = lob(_sm,&_a0(0),&_b0(0),-1.0,1.0,&pts(0),&e3(0),&e(0),&e1(0),&e2(0));		   
//	  if (ierr != 0) {
//		  printf("gauss #3 error %d\n",ierr);
//		  exit(1);
//	  }

	/* CALCULATE PROJECTION POINTS IN INTERIOR */
	for(i=1;i<_sm;++i) {
		for(j=1;j<_sm-(i-1);++j) {
			s = -1 +2.0*((FLT) j)/(FLT)(_sm+1);
			r = -1 +2.0*((FLT) i)/(FLT)(_sm+1);
			x = 2.0*(1+r)/(1-s) -1.0;
			eta = s;
					  
			/* CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
			ptvalues(x,eta);
	  
			/* CALCULATE S POLYNOMIALS */
			/* VERTEX A	   */
			_lgrnge(0,i,j) = pgn(0)*pgx(0);

			/* VERTEX B	 */
			_lgrnge(1,i,j) = pgn(1)*pgx(1);

			/* VERTEX C		*/	  
			_lgrnge(2,i,j) = pgn(2)*pgx(2);

			/*	SIDE 1 (s)		  */
			for(m = 3; m < _sm+3; ++m)
				_lgrnge(m,i,j) = pgn(m)*pgx(m);
				
			for(m=_sm+3;m<2*_sm+3;++m)
				_lgrnge(m,i,j) = pgn(m)*pgx(2);

			for(m=2*_sm+3;m<_bm;++m)
				_lgrnge(m,i,j) = pgn(m)*pgx(1);

			/*	INTERIOR MODES	  */
			ind = _bm;
			for(m = 3; m< _sm+2;++m) {		 
				for(n = 1; n < _sm+3-m;++n) {
					_lgrnge(ind,i,j) = pgn(ind)*pgx(m);
					++ind;
				}
			}
		}
	}
	
	/* NOW CALCULATE VALUES OF G FOR SIDE PROJECTION */
	for (i=1;i<_sm+1;++i) {
		x = 2.0*(FLT) i/(FLT)(_sm+1) -1.0;
		ptvalues1d(x);
		for(m=0;m<_p+1;++m)
			_lgrnge1d(m,i) = pgx(m);
	}
	
#ifdef DEBUG	
	Array<FLT,1> test(_tm), rslt(_tm), test1(_tm), uht(_tm); 
	
	test(0) = func(-1.0,1.0);
	test(1) = func(-1.0,-1.0);
	test(2) = func(1.0,-1.0);
	
	for (i=1;i<_sm+1;++i) {
		x = 2.0*(FLT) i/(FLT)(_sm+1) -1.0;
		test(i+2) = func(x,-1.0);
		test(i+2+_sm) = func(-x,x);
		test(i+2+2*_sm) = func(-1.,-x);
	}
	
	int count = _bm;
	for(i=1;i<_sm;++i) {
			for(j=1;j<_sm-(i-1);++j) {
					s = -1 +2.0*((FLT) j)/(FLT)(_sm+1);
					r = -1 +2.0*((FLT) i)/(FLT)(_sm+1);	
					test(count++) = func(r,s);
			}
	}
	std::cout << test << std::endl;

	legtobasis(test.data(),rslt.data());
	
	std::cout << rslt << std::endl;
	
	test1(Range(0,2)) = rslt(Range(0,2));
	if (_sm) {
		for(i=0;i<3;++i) {
			uht(0) = rslt((i+1)%3);
			uht(1) = rslt((i+2)%3);
			uht(Range(2,_sm+1)) = rslt(Range(3+i*_sm,3+(i+1)*_sm-1));
			proj1d_leg(uht.data(),&test1(2+i*_sm));
		}
		
		Array<FLT,2> d1_leg(_sm,_sm);
		proj_leg(rslt.data(),d1_leg.data(), _sm);
		count = _bm;
		for(i=1;i<_sm;++i) {
				for(j=1;j<_sm-(i-1);++j) {	
						test1(count++) = d1_leg(i,j);
				}
		}
	}
	test -= test1;
	std::cout << test << std::endl;
#endif
	
	
	return;
}

template<int _p, int ep> void tri_basis<_p,ep>::legtobasis(const FLT *data, FLT *coeff) const {
	int i,j,m,n;
	char trans[] = "T";

	
	/* Vertex coefficients are the same */
	for(int i=0;i<3;++i)
		coeff[i] = data[i];
		
	/* Side coefficients */
	TinyMatrix<FLT,_tm,_tm> matrix;
	TinyVector<int,2*_tm> ipiv;
	int info;

	/* REVERSE OUTPUTING PROCESS */
	for(m=0;m<_sm;++m)
		for(n=0;n<_sm;++n)
			matrix(n,m) = _lgrnge1d(m+2,n+1);

	GETRF(_sm,_sm,matrix.data(),_tm,ipiv.data(),info);
	if (info != 0) {
		printf("DGETRF FAILED FOR INPUTING TECPLOT SIDES\n");
		exit(1);
	}
	
	TinyVector<FLT,_tm> uht;
	TinyVector<FLT,_tm> u1d;
	
	int count1 = 3;
	int count2 = 3;
	for (int i=0;i<3;++i) {
		uht = 0.0;
		uht(0) = data[(i+1)%3];
		uht(1) = data[(i+2)%3];
		proj1d_leg(uht.data(),u1d.data());
		for(m=0;m<_sm;++m) {
			u1d(m+1) -= data[count1++];
		}
		GETRS(trans,_sm,1,matrix.data(),_tm,ipiv.data(),u1d.data()+1,_tm,info);
		for(m=0;m<_sm;++m)
			coeff[count2++] = -u1d(1+m);
	}
		
	/* Interior modes */
	for(int m=0;m<_im;++m) {
		n = 0;
		for(int i=1;i<_sm;++i) {
			for(int j=1;j<_sm-(i-1);++j) {
				matrix(n++,m) = _lgrnge(m+_bm,i,j);
			}
		}
	}

	GETRF(_im,_im,matrix.data(),_tm,ipiv.data(),info);
	if (info != 0) {
		printf("DGETRF FAILED FOR INPUTING TECPLOT SIDES\n");
		exit(1);
	}

	for (int i=_bm;i<_tm;++i)
		coeff[i] = 0.0;
		
	TinyMatrix<FLT,_sm+1,_sm+1> u2d;

	proj_leg(coeff,u2d.data(),_sm);

	m = 0;
	for(i=1;i<_sm;++i) {
		for(j=1;j<_sm-(i-1);++j) {
			uht(m) = u2d(i,j) -data[m+_bm];	 
			++m;
		}
	}
	GETRS(trans,_im,1,matrix.data(),_tm,ipiv.data(),uht.data(),_tm,info);
	for(m=0;m<_im;++m)
		coeff[m+_bm] = -uht(m);
		
	return;
}

