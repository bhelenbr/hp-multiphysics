/*
 *	ptvalues.cpp
 *	spectral_hp
 *
 *	Created by Brian Helenbrook on Tue Apr 16 2002.
 *	Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_basis.h"
#include<math.h>

/* CALCULATES VALUES OF _gx POLYNOMIALS & GS POLYNOMIALS AT POINT */
template<int _p, int ep> void tri_basis<_p,ep>::ptvalues(FLT x, FLT eta) const {
	 FLT pkp,pk,pkm;
	 int k,m,ind;
	
	/* CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
	pgx(0) = 1.0;
	pgx(1) = .5*(1-x);
	pgx(2) = .5*(1+x);

	/* SIDE 0	 */
	/* CALCULATE _p, _p' USING RECURSION RELATION */
	pk = 1.0;
	pkm = 0.0;
	for (m = 0;m < _sm;++m) {
		pgx(m+3) = pgx(1)*pgx(2)*pk*_norm(m+3);
		pkp = (x-_a0(0,m))*pk - _b0(0,m)*pkm;
		pkm = pk;
		pk = pkp;
	}


	/* CALCULATE S POLYNOMIALS */
	ind = 0;
	
	/* VERTEX 0,1,2	   */
	pgn(ind++) = (1+eta)*.5;
	pgn(ind++) = (1-eta)*.5;
	pgn(ind++) = (1-eta)*.5;

	/* SIDE 0 (s)		 */
	for(m = 2; m <= _sm+1; ++m)
		pgn(ind++) = pow(.5*(1-eta),m);

	/* SIDE 1	 */
	pk = 1.0;
	pkm = 0.0;
	for(m=0;m<_sm;++m) {
		pgn(ind++) = (1.-eta)*(1.+eta)*.25*pk*_norm(m+3);
		pkp = (eta-_a0(1,m))*pk - _b0(1,m)*pkm;
		pkm = pk;
		pk = pkp;
	}

	/* SIDE 2	 */
	pk = 1.0;
	pkm = 0.0;
	for(m = 0;m<_sm;++m) {
		pgn(ind++) = (m % 2 ? -1 : 1)*(1.-eta)*(1.+eta)*.25*pk*_norm(m+3);
		pkp = (eta-_a0(1,m))*pk - _b0(1,m)*pkm;
		pkm = pk;
		pk = pkp;
	}

	/* INTERIOR MODES	 */
	for(m = 2; m< _sm+1;++m) {		 
		pk = 1.0;
		pkm = 0.0;
		for(k = 0; k < _sm+1-m;++k) {
			pgn(ind) = pow(.5*(1.-eta),m)*.5*(1.+eta)*pk*_norm(ind);
			pkp = (eta-_a0(m,k))*pk - _b0(m,k)*pkm;
			pkm = pk;
			pk = pkp;
			++ind;
		}
	}
	
	return;
}

/* CALCULATE VALUE OF G(X) & DG/DX, G(eta), DG/Deta AT POINT */
template<int _p, int ep> void tri_basis<_p,ep>::ptvalues_deriv(FLT x, FLT eta) const {
	FLT pkp,pk,pkm,dpk,dpkm,dpkp;
	int k,m,n,ind;
	
	/* CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
	pgx(0) = 1.0;
	dpgx(0) = 0.0;
	
	pgx(1) = .5*(1-x);
	dpgx(1) = -0.5;
	
	pgx(2) = .5*(1+x);
	dpgx(2) = 0.5;


	/* SIDE 0	 */
	/* CALCULATE _p, _p' USING RECURSION RELATION */
	pk = 1.0;
	pkm = 0.0;
	dpk = 0.0;
	dpkm = 0.0;
	for (m = 0;m < _sm;++m) {
		pgx(m+3) = (1.+x)*(1.-x)*.25*pk*_norm(m+3);
		dpgx(m+3) = (-x*.5*pk +(1.+x)*(1.-x)*.25*dpk)*_norm(m+3);
		pkp = (x-_a0(0,m))*pk - _b0(0,m)*pkm;
		dpkp = pk + (x-_a0(0,m))*dpk - _b0(0,m)*dpkm;
		dpkm = dpk;
		dpk = dpkp;
		pkm = pk;
		pk = pkp;
	}

	/******************************************/
	/* GENERATE JACOBI POLY FOR S DIRECTION */
	/****************************************/
	ind = 0;

	/* VERTEX 0,1,2	   */
	pgn(ind) = (1+eta)*.5;
	dpgn(ind++) = 0.5;
	
	pgn(ind) = (1-eta)*.5;
	dpgn(ind++) = -0.5;

	pgn(ind) = (1-eta)*.5;
	dpgn(ind++) = -0.5;

	/* SIDE 0 (s)		 */
	for(m = 2; m <= _sm+1; ++m) {
		pgn(ind) = pow(.5*(1-eta),m);
		dpgn(ind) = -.5*m*pow(.5*(1.-eta),m-1);
		++ind;
	}

	/* SIDE 1	 */
	pk = 1.0;
	pkm = 0.0;
	dpk = 0.0;
	dpkm = 0.0;
	for(n=0;n<_sm;++n) {
		pgn(ind) = (1.-eta)*(1.+eta)*.25*pk*_norm(n+3);
		dpgn(ind) = (-.5*eta*pk +(1.-eta)*(1.+eta)*.25*dpk)*_norm(n+3);
		pkp = (eta-_a0(1,n))*pk - _b0(1,n)*pkm;
		dpkp = pk + (eta-_a0(1,n))*dpk - _b0(1,n)*dpkm;
		dpkm = dpk;
		dpk = dpkp;
		pkm = pk;
		pk = pkp;
		++ind;
	}

	/* SIDE 2	 */
	pk = 1.0;
	pkm = 0.0;
	dpk = 0.0;
	dpkm = 0.0;

	for(n=0;n<_sm;++n) {
		pgn(ind) = (n % 2 ? -1 : 1)*(1.-eta)*(1.+eta)*.25*pk*_norm(n+3);
		dpgn(ind) = (n % 2 ? -1 : 1)*(-.5*eta*pk +(1.-eta)*(1.+eta)*.25*dpk)*_norm(n+3);
		pkp = (eta-_a0(1,n))*pk - _b0(1,n)*pkm;
		dpkp = pk + (eta-_a0(1,n))*dpk - _b0(1,n)*dpkm;
		dpkm = dpk;
		dpk = dpkp;
		pkm = pk;
		pk = pkp;
		++ind;
	}

	
	/*	INTERIOR MODES	  */
	for(m = 2; m< _sm+1;++m) {		 
		pk = 1.0;
		pkm = 0.0;
		dpk = 0.0;
		dpkm = 0.0;

		for(k = 0; k < _sm+1-m;++k) {
			pgn(ind) = (pow(.5*(1.-eta),m)*.5*(1.+eta)*pk)*_norm(ind);
			dpgn(ind) = (-.25*m*pow(.5*(1.-eta),m-1)*(1.+eta)*pk +pow(.5*(1.-eta),m)*.5*(pk + (1.+eta)*dpk))*_norm(ind);
			pkp = (eta-_a0(m,k))*pk - _b0(m,k)*pkm;
			dpkp = pk + (eta-_a0(m,k))*dpk - _b0(m,k)*dpkm;
			dpkm = dpk;
			dpk = dpkp;
			pkm = pk;
			pk = pkp;
			++ind;
		}
	}

	return;
}

/* CALCULATES VALUES OF _gx POLYNOMIALS & GS POLYNOMIALS AT POINT */
/* BOUNDARY MODES ONLY */
template<int _p, int ep> void tri_basis<_p,ep>::ptvalues_bdry(FLT x, FLT eta) const {
	 FLT pkp,pk,pkm;
	 int m,ind;
	
	/* CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
	pgx(0) = 1.0;
	pgx(1) = .5*(1-x);
	pgx(2) = .5*(1+x);

	/* SIDE 0	 */
	/* CALCULATE _p, _p' USING RECURSION RELATION */
	pk = 1.0;
	pkm = 0.0;
	for (m = 0;m < _sm;++m) {
		pgx(m+3) = pgx(1)*pgx(2)*pk*_norm(m+3);
		pkp = (x-_a0(0,m))*pk - _b0(0,m)*pkm;
		pkm = pk;
		pk = pkp;
	}


	/* CALCULATE S POLYNOMIALS */
	ind = 0;
	
	/* VERTEX 0,1,2	   */
	pgn(ind++) = (1+eta)*.5;
	pgn(ind++) = (1-eta)*.5;
	pgn(ind++) = (1-eta)*.5;

	/* SIDE 0 (s)		 */
	for(m = 2; m <= _sm+1; ++m)
		pgn(ind++) = pow(.5*(1-eta),m);

	/* SIDE 1	 */
	pk = 1.0;
	pkm = 0.0;
	for(m=0;m<_sm;++m) {
		pgn(ind++) = (1.-eta)*(1.+eta)*.25*pk*_norm(m+3);
		pkp = (eta-_a0(1,m))*pk - _b0(1,m)*pkm;
		pkm = pk;
		pk = pkp;
	}

	/* SIDE 2	 */
	pk = 1.0;
	pkm = 0.0;
	for(m = 0;m<_sm;++m) {
		pgn(ind++) = (m % 2 ? -1 : 1)*(1.-eta)*(1.+eta)*.25*pk*_norm(m+3);
		pkp = (eta-_a0(1,m))*pk - _b0(1,m)*pkm;
		pkm = pk;
		pk = pkp;
	}
	
	return;
}

/* CALCULATE VALUE OF G(X) & DG/DX, G(eta), DG/Deta AT POINT FOR ONLY BOUNDARY MODES */
template<int _p, int ep> void tri_basis<_p,ep>::ptvalues_deriv_bdry(FLT x, FLT eta) const {
	FLT pkp,pk,pkm,dpk,dpkm,dpkp;
	int m,n,ind;

	/* CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
	pgx(0) = 1.0;
	dpgx(0) = 0.0;

	pgx(1) = .5*(1-x);
	dpgx(1) = -0.5;
	
	pgx(2) = .5*(1+x);
	dpgx(2) = 0.5;


	/* SIDE 0	 */
	/* CALCULATE _p, _p' USING RECURSION RELATION */
	pk = 1.0;
	pkm = 0.0;
	dpk = 0.0;
	dpkm = 0.0;
	for (m = 0;m < _sm;++m) {
		pgx(m+3) = (1.+x)*(1.-x)*.25*pk*_norm(m+3);
		dpgx(m+3) = (-x*.5*pk +(1.+x)*(1.-x)*.25*dpk)*_norm(m+3);
		pkp = (x-_a0(0,m))*pk - _b0(0,m)*pkm;
		dpkp = pk + (x-_a0(0,m))*dpk - _b0(0,m)*dpkm;
		dpkm = dpk;
		dpk = dpkp;
		pkm = pk;
		pk = pkp;
	}

	/******************************************/
	/* GENERATE JACOBI POLY FOR S DIRECTION */
	/****************************************/
	ind = 0;

	/* VERTEX 0,1,2	   */
	pgn(ind) = (1+eta)*.5;
	dpgn(ind++) = 0.5;
	
	pgn(ind) = (1-eta)*.5;
	dpgn(ind++) = -0.5;

	pgn(ind) = (1-eta)*.5;
	dpgn(ind++) = -0.5;

	/* SIDE 0 (s)		 */
	for(m = 2; m <= _sm+1; ++m) {
		pgn(ind) = pow(.5*(1-eta),m);
		dpgn(ind) = -.5*m*pow(.5*(1.-eta),m-1);
		++ind;
	}

	/* SIDE 1	 */
	pk = 1.0;
	pkm = 0.0;
	dpk = 0.0;
	dpkm = 0.0;
	for(n=0;n<_sm;++n) {
		pgn(ind) = (1.-eta)*(1.+eta)*.25*pk*_norm(n+3);
		dpgn(ind) = (-.5*eta*pk +(1.-eta)*(1.+eta)*.25*dpk)*_norm(n+3);
		pkp = (eta-_a0(1,n))*pk - _b0(1,n)*pkm;
		dpkp = pk + (eta-_a0(1,n))*dpk - _b0(1,n)*dpkm;
		dpkm = dpk;
		dpk = dpkp;
		pkm = pk;
		pk = pkp;
		++ind;
	}

	/* SIDE 2	 */
	pk = 1.0;
	pkm = 0.0;
	dpk = 0.0;
	dpkm = 0.0;

	for(n=0;n<_sm;++n) {
		pgn(ind) = (n % 2 ? -1 : 1)*(1.-eta)*(1.+eta)*.25*pk*_norm(n+3);
		dpgn(ind) = (n % 2 ? -1 : 1)*(-.5*eta*pk +(1.-eta)*(1.+eta)*.25*dpk)*_norm(n+3);
		pkp = (eta-_a0(1,n))*pk - _b0(1,n)*pkm;
		dpkp = pk + (eta-_a0(1,n))*dpk - _b0(1,n)*dpkm;
		dpkm = dpk;
		dpk = dpkp;
		pkm = pk;
		pk = pkp;
		++ind;
	}

	return;
}

template<int _p, int ep> void tri_basis<_p,ep>::ptvalues1d(FLT x) const {
	 FLT pkp,pk,pkm;
	 int m;
	
	/* CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
	pgx(0) = .5*(1-x);
	pgx(1) = .5*(1+x);

	/* SIDE 1	 */
	/* CALCULATE _p, _p' USING RECURSION RELATION */
	pk = 1.0;
	pkm = 0.0;
	for (m = 0;m < _sm;++m) {
		pgx(m+2) = (1.+x)*(1.-x)*.25*pk*_norm(m+3);
		pkp = (x-_a0(0,m))*pk - _b0(0,m)*pkm;
		pkm = pk;
		pk = pkp;
	}
	
	return;
}

template<int _p, int ep> void tri_basis<_p,ep>::ptvalues1d_deriv(FLT x) const {
	 FLT pkp,pk,pkm,dpkp,dpk,dpkm;
	 int m;

	pgx(0)	= .5*(1-x);
	dpgx(0) = -.5;
	
	pgx(1)	= .5*(1+x);
	dpgx(1) = .5;

	/* SIDE 1 CALCULATE _p, _p' USING RECURSION RELATION */
	pk = 1.0;
	pkm = 0.0;
	dpk = 0.0;
	dpkm = 0.0;

	for (m = 0;m < _sm;++m) {
		pgx(m+2) = (1.+x)*(1.-x)*.25*pk*_norm(m+3);
		dpgx(m+2) = (-x*.5*pk +(1.+x)*(1.-x)*.25*dpk)*_norm(m+3);
		pkp = (x-_a0(0,m))*pk - _b0(0,m)*pkm;
		dpkp = pk + (x-_a0(0,m))*dpk - _b0(0,m)*dpkm;
		dpkm = dpk;
		dpk = dpkp;
		pkm = pk;
		pk = pkp;
	}
	
	return;
}
