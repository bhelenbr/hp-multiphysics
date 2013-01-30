/*
 *	tri_basis.h
 *	planar++
 *
 *	Created by helenbrk on Fri Oct 12 2001.
 *	Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _tri_basis_h_
#define _tri_basis_h_

#include <blitz/array.h>
#include <float.h>

#ifdef SINGLE
#define FLT float
#define EPSILON FLT_EPSILON
#else
#ifndef FLT
#define FLT double
#define EPSILON DBL_EPSILON
#endif
#endif

#define MORTHOGONAL

using namespace blitz;

/* VIRTUAL INTERFACE FOR BASIS */
class tri_basis_interface {
	public:
		tri_basis_interface() {}
		virtual ~tri_basis_interface() {}
		
		/* NUMBER OF MODES : VERTEX/SIDE/INTERIOR/TOTAL/BOUNDARY */
		virtual int p() = 0;
		virtual int sm() = 0;
		virtual int im() = 0;
		virtual int tm() = 0;
		virtual int bm() = 0;
		
		/* NUMBER OF X/S MODES & GAUSS POINTS */		
		virtual int nmodx() = 0;
		virtual int nmodn() = 0;
		virtual int gpx() = 0;
		virtual int gpn() = 0;
		/* BANDWITH OF INTERIOR/INTERIOR MATRIX */
		virtual int ibwth() = 0;
		/* BANDWITH OF SIDE/SIDE MATRIX (1D) */
		virtual int sbwth() = 0;
		
		/* X FUNCTIONS & DERIVATIVES */
		virtual FLT gx(int,int) = 0;
		virtual FLT dgx(int,int) = 0;
		
		/* GAUSS WEIGTS & LOCATIONS */
		virtual FLT wtx(int) = 0;
		virtual FLT xp(int) = 0;
		virtual FLT x0(int) = 0;

		/* COMBINED THINGS FOR FAST INTEGRATION */
		virtual FLT gxwtx(int,int) = 0;
		virtual FLT dgxwtx(int,int) = 0;

		/* TO TAKE X,Y DERIVATIVES OF A FUNCTION WITH VALUES ON GAUSS POINTS */
		virtual FLT dltx(int) = 0;
		virtual FLT dltx1(int,int) = 0;
		
		/* ETA FUNCTIONS & DERIVATIVES */
		virtual FLT gn(int,int) = 0;
		virtual FLT dgn(int,int) = 0;
		/* GAUSS WEIGHTS & LOCATIONS */
		virtual FLT wtn(int) = 0; 
		virtual FLT np(int) = 0;
		virtual FLT n0(int) = 0;

		/* COMBINED THINGS FOR FAST INTEGRATION */
		virtual FLT gnwtn(int,int) = 0;
		virtual FLT gnwtnn0(int,int) = 0;
		virtual FLT dgnwtn(int,int) = 0;
		
		/* TO TAKE X,Y DERIVATIVES OF A FUNCTION WITH VALUES ON GAUSS POINTS */
		virtual FLT dltn(int) = 0;
		virtual FLT dltn1(int,int) = 0;
		virtual FLT dltn2(int,int) = 0;

		/* RENORMALIZATION CONSTANTS */
		virtual FLT norm(int) = 0;
		
		/* FOR OUTPUTING TO LEGENDRE POINTS */			  
		virtual FLT lgrnge1d(int,int) = 0;
		virtual FLT lgrnge(int,int,int) = 0;
		
		/* FOR CALCULATING NORMAL DERIVATIVES ALONG SIDES */
		virtual FLT dgnorm(int,int,int) = 0;

		/* LUMPED MASS MATRIX INVERSION */
		/* REMOVE SIDE COMPONENT FROM VERTICES */
		virtual const FLT& sfmv(int,int) = 0;
		/* REMOVE INTERIOR COMPONENT FROM VERTICES & SIDES */
		virtual const FLT& ifmb(int,int) = 0;
		/* DIAGONAL FOR VERTICES AFTER REMOVAL OF SIDES & INTERIORS */
		virtual const FLT& vdiag() = 0;
		/* AFTER FINDING VERTEX REMOVE FROM SIDE MODES */
		virtual const FLT& vfms(int,int) = 0;
		/* DIAGONAL FOR SIDES */		
		virtual const FLT& sdiag(int) = 0;
		/* AFTER FINDING LOW ORDER SIDES REMOVE FROM HIGH ORDER SIDE MODES */
		virtual const FLT& sfms(int,int,int) = 0;
		/* AFTER FINDING VERTICES & ALL SIDES REMOVE FROM INTERIOR */
		virtual const FLT& bfmi(int,int) = 0;
		/* INTERIOR DIAGONAL */
		virtual const FLT& idiag(int,int) = 0;
		/* MASS MATRIX WITH STATIC INVERSION OF INTERIOR MODES */
		virtual const FLT& msi(int,int) = 0;
		
		/*1D (SIDE) MASS MATRIX INVERSION */
		/* REMOVE SIDE COMPONENT FROM VERTICES */
		virtual const FLT& sfmv1d(int,int) = 0;
		/* DIAGONAL FOR VERTICES AFTER REMOVAL OF SIDES */
		virtual const FLT& vdiag1d() = 0;
		/* AFTER FINDING VERTEX REMOVE FROM SIDE MODES */
		virtual const FLT& vfms1d(int,int) = 0;
		/* DIAGONAL FOR SIDES */		
		virtual const FLT& sdiag1d(int,int) = 0;
			
		virtual void initialize() = 0;
		/* PROJECT WITH R & S DERIVATIVES */
		virtual void proj(const FLT *lin, FLT *f, FLT *dx, FLT *dy, int stride) const = 0;
		/* PROJECT ONLY VALUE */
		virtual void proj(const FLT *lin, FLT *f, int stride) const = 0;
		/* PROJECT A LINEAR FUNCTION */
		virtual void proj(FLT u1, FLT u2, FLT u3, FLT *f, int stride) const = 0;
		/* PROJECT USING SIDE/VERTEX MODES WITH R & S DERIVATIVES */
		virtual void proj_bdry(const FLT *lin, FLT *f, FLT *dr, FLT *ds, int stride) const = 0;
		/* PROJECT USING SIDE/VERTEX MODES ONLY */
		virtual void proj_bdry(const FLT *lin, FLT *f, int stride) const = 0;
		/* PROJECT TO LEGENDRE POINTS (FOR OUTPUTING) */
		virtual void proj_leg(const FLT *lin, FLT *f, int stride) const = 0;
		/* PROJECT LINEAR FUNCTION TO LEGENDRE POINTS (FOR OUTPUTING) */
		virtual void proj_leg(FLT u1, FLT u2, FLT u3, FLT *f, int stride) const = 0;
		/* PROJECT TO LEGENDRE POINTS USING SIDE/VERTEX MODES ONLY (FOR OUTPUTING) */
		virtual void proj_bdry_leg(const FLT *lin, FLT *f, int stride) const = 0;
		/* PROJECT VALUES & TANGENT & NORMAL DERIVATIVES TO 1D SIDE GAUSS POINTS */
		/* dx is tangential derivative, dn is normal derivative to side */
		virtual void proj_side(int side,const FLT *lin, FLT *f, FLT *dx, FLT *dn) const = 0;

		/* DERIVATIVE IN R */
		virtual void derivr(const FLT *f, FLT *dr, int stride) const = 0;
		/* DERIVATIVE IN S */
		virtual void derivs(const FLT *f, FLT *ds, int stride) const = 0;
		/* INTEGRATE WITH RESPECT TO BASIS */
		virtual void intgrt(FLT *rslt,const FLT *f, int stride) const = 0;
		/* INTEGRATE FX W/RSPCT TO DG/DR & FY W/RSPCT TO DG/DS */
		/* WARNING THIS ADDS INTEGRATION TO RESULT: RESULT IS NOT CLEARED FIRST */
		virtual void intgrtrs(FLT *rslt,const FLT *dx,const FLT *dy, int stride) const = 0;

		/* SAME STUFF EXCEPT 1D	   */	 
		/* PROJECT WITH X DERIVATIVES */
		virtual void proj1d(const FLT *lin, FLT *f, FLT *dx) const = 0;
		/* PROJECT ONLY VALUE */
		virtual void proj1d(const FLT *lin, FLT *f) const = 0;
		/* PROJECT A LINEAR FUNCTION */
		virtual void proj1d(FLT u1, FLT u2, FLT *f) const = 0;
		/* PROJECT TO LEGENDRE POINTS (FOR OUTPUTING) */
		virtual void proj1d_leg(const FLT *lin, FLT *f) const = 0;
		/* PROJECT TO LEGENDRE POINTS USING SIDE/VERTEX MODES ONLY (FOR OUTPUTING) */
		virtual void proj1d_leg(FLT u1, FLT u2, FLT *f) const = 0;
		/* DERIVATIVE IN X */
		virtual void derivx1d(const FLT *f, FLT *dx) const = 0;
		/* INTEGRATE WITH RESPECT TO BASIS */
		virtual void intgrt1d(FLT *rslt,const FLT *f) const = 0;
		/* INTEGRATE W/RSPCT TO DG/DX */
		virtual void intgrtx1d(FLT *rslt,const FLT *f) const = 0;
		
		/* TO CALCULATE BASIS FUNCTIONS & DERIVATIVES */
		virtual void ptvalues(FLT xi, FLT s) const = 0; // CALCULATES GX, gn VALUES AT A POINT
		virtual void ptvalues_deriv(FLT xi, FLT s) const = 0; // CALCULATES GX, DGX, GN, DGN AT A POINT
		virtual void ptvalues_bdry(FLT xi, FLT s) const = 0; // CALCULATES GX, gn VALUES AT A POINT (BDRY MODES ONLY)
		virtual void ptvalues_deriv_bdry(FLT xi, FLT s) const = 0; // CALCULATES GX, gn & DERIV VALUES AT A POINT (BDRY MODES ONLY)
		
		void ptvalues_rs(FLT r, FLT s)	const {ptvalues(2.0*(1+r)/(1-s+10.*EPSILON) -1.0,s);} // CALCULATES GX, gn VALUES AT A POINT
		void ptvalues_deriv_rs(FLT r, FLT s) const {ptvalues_deriv(2.0*(1+r)/(1-s+10.*EPSILON) -1.0,s);} // CALCULATES GX, DGX, GN, DGN AT A POINT
		void ptvalues_bdry_rs(FLT r, FLT s) const {ptvalues_bdry(2.0*(1+r)/(1-s+10.*EPSILON) -1.0,s);} // CALCULATES GX, gn VALUES AT A POINT (BDRY MODES ONLY)
		void ptvalues_deriv_bdry_rs(FLT r, FLT s) const {ptvalues_deriv_bdry(2.0*(1+r)/(1-s+10.*EPSILON) -1.0,s);} // CALCULATES GX, gn & DERIV VALUES AT A POINT (BDRY MODES ONLY)
		
		/* 1D SIDE BASIS FUNCTIONS & DERIVATIVES */
		virtual void ptvalues1d(FLT x) const = 0;
		virtual void ptvalues1d_deriv(FLT x) const = 0;
			  
		/* POINT PROBE IN STANDARD ELEMENT FOR VECTOR */
		virtual void ptprobe(int nv, FLT *f,const FLT *lin, int stride) const = 0;	// REUSES OLD R,S
		virtual void ptprobe(int nv, FLT *f, FLT r, FLT s,const FLT *lin, int stride) const  = 0;
		virtual void ptprobe(int nv, FLT *f1, FLT *dx1, FLT *dy1, FLT r, FLT s, const FLT *lin1, int stride) const = 0;

		/* POINT PROBE FUNCTIONS USING BOUNDARY MODES ONLY */
		virtual void ptprobe_bdry(int nv, FLT *f, const FLT *lin, int stride) const = 0; // REUSES OLD R,S ONLY BDRY MODES 
		virtual void ptprobe_bdry(int nv, FLT *f, FLT *dx, FLT *dy, FLT r, FLT s, const FLT *lin, int stride) const = 0; // BOUNDARY MODES ONLY CALC'S DERIVATIVES
		virtual void ptprobe_bdry(int nv, FLT *f, FLT r, FLT s, const FLT *lin, int stride) const = 0;
		
		/* 1D SIDE PROBE FUNCTIONS */
		virtual void ptprobe1d(int nv, FLT *f, const FLT *sin, int stride) const = 0;	 // REUSES OLD VALUES OF X
		virtual void ptprobe1d(int nv, FLT *f, FLT x, const FLT *sin, int stride) const = 0; // CALCULATES MODES AT X
		virtual void ptprobe1d(int nv, FLT *f, FLT *dx, const FLT *sin, int stride) const = 0; // REUSES & CALCULATES DERIVATIVE 
		virtual void ptprobe1d(int nv, FLT *f, FLT *dx, FLT x, const FLT *sin, int stride) const = 0;  // CALCULATES MODES & DERIVS AT X 
		
		/* Utility for switching between uniform legendre representation & this basis */
		virtual void legtobasis(const FLT *data, FLT *coeff) const = 0;
};

template<int _p,int ep> class tri_basis : public tri_basis_interface {	
	private:
		/* NUMBER OF MODES : VERTEX,SIDE/INTERIOR/TOTAL/BOUNDARY */
		static const int _sm = _p-1;
		static const int _im = (_p-1)*(_p-2)/2;
		static const int _tm = (_p+1)*(_p+2)/2;
		static const int _bm = 3*_p;
		/* NUMBER OF X/S MODES & GAUSS POINTS */		
		static const int _nmodx = _p+2;
		static const int _nmodn = _tm;
		static const int _gpx = _p+1+ep;
		static const int _gpn = _p+1+ep;
		/* BANDWITH OF INTERIOR/INTERIOR MATRIX */
		static const int _ibwth = 1;
		/* BANDWITH OF SIDE/SIDE MATRIX (1D) */
		static const int _sbwth = 2;
		/* RECURSION RELATION COEFFICIENTS */
		TinyMatrix<FLT,_sm+2,_gpx+1> _a0, _b0;
		
		/* X FUNCTIONS & DERIVATIVES */
		TinyMatrix<FLT,_gpx,_nmodx> _gx, _dgx;
		/* GAUSS WEIGTS & LOCATIONS */
		TinyVector<FLT,_gpx> _wtx, _xp, _x0;
		/* COMBINED THINGS FOR FAST INTEGRATION */
		TinyMatrix<FLT,_nmodx,_gpx> _gxwtx, _dgxwtx;
		/* TO TAKE X,Y DERIVATIVES OF A FUNCTION WITH VALUES ON GAUSS POINTS */
		TinyVector<FLT,_gpx> _dltx;
		TinyMatrix<FLT,_gpx,_gpx> _dltx1;
		
		/* ETA FUNCTIONS & DERIVATIVES */
		TinyMatrix<FLT,_gpn,_tm> _gn, _dgn;
		/* GAUSS WEIGTS & LOCATIONS */
		TinyVector<FLT,_gpn> _wtn, _np, _n0;
		/* COMBINED THINGS FOR FAST INTEGRATION */
		TinyMatrix<FLT,_tm,_gpn> _gnwtn, _gnwtnn0, _dgnwtn;
		/* TO TAKE X,Y DERIVATIVES OF A FUNCTION WITH VALUES ON GAUSS POINTS */
		TinyVector<FLT,_gpn> _dltn;
		TinyMatrix<FLT,_gpx,_gpx> _dltn1;
		TinyMatrix<FLT,_gpn,_gpn> _dltn2;

		/* RENORMALIZATION CONSTANTS */
		TinyVector<FLT,_tm> _norm;
		
		/* FOR OUTPUTING TO LEGENDRE POINTS */			  
		Array<FLT,2> _lgrnge1d;
		Array<FLT,3> _lgrnge;
		
		/* FOR CALCULATING NORMAL DERIVATIVES ALONG SIDES */
		Array<FLT,3> _dgnorm;

		/* LUMPED MASS MATRIX INVERSION */
		/* REMOVE SIDE COMPONENT FROM VERTICES */
		Array<FLT,2> _sfmv;
		/* REMOVE INTERIOR COMPONENT FROM VERTICES & SIDES */
		Array<FLT,2> _ifmb;
		/* DIAGONAL FOR VERTICES AFTER REMOVAL OF SIDES & INTERIORS */
		FLT _vdiag;
		/* AFTER FINDING VERTEX REMOVE FROM SIDE MODES */
		Array<FLT,2> _vfms;
		/* DIAGONAL FOR SIDES */		
		Array<FLT,1> _sdiag;
		/* AFTER FINDING LOW ORDER SIDES REMOVE FROM HIGH ORDER SIDE MODES */
		Array<FLT,3> _sfms;
		/* AFTER FINDING VERTICES & ALL SIDES REMOVE FROM INTERIOR */
		Array<FLT,2> _bfmi;
		/* INTERIOR DIAGONAL */
		Array<FLT,2> _idiag;
		/* MASS MATRIX WITH STATIC INVERSION OF INTERIOR MODES */
		Array<FLT,2> _msi;
		
		/*1D (SIDE) MASS MATRIX INVERSION */
		/* REMOVE SIDE COMPONENT FROM VERTICES */
		Array<FLT,2> _sfmv1d;
		/* DIAGONAL FOR VERTICES AFTER REMOVAL OF SIDES */
		FLT _vdiag1d;
		/* AFTER FINDING VERTEX REMOVE FROM SIDE MODES */
		Array<FLT,2> _vfms1d;
		/* DIAGONAL FOR SIDES */		
		Array<FLT,2> _sdiag1d;
			
	public:
		/* NUMBER OF MODES : VERTEX,SIDE/INTERIOR/TOTAL/BOUNDARY */
		int p() {return(_p);}
		int sm() {return(_sm);}
		int im() {return(_im);}
		int tm() {return(_tm);}
		int bm() {return(_bm);}
		
		/* NUMBER OF X/S MODES & GAUSS POINTS */		
		int nmodx() {return(_nmodx);}
		int nmodn() {return(_nmodn);}
		int gpx() {return(_gpx);}
		int gpn() {return(_gpn);}
		/* BANDWITH OF INTERIOR/INTERIOR MATRIX */
		int ibwth() {return(_ibwth);}
		/* BANDWITH OF SIDE/SIDE MATRIX (1D) */
		int sbwth() {return(_sbwth);}
		
		/* X FUNCTIONS & DERIVATIVES */
		FLT gx(int pnt,int mod) {return(_gx(pnt,mod));}
		FLT dgx(int pnt,int mod) {return(_dgx(pnt,mod));}
		
		/* GAUSS WEIGTS & LOCATIONS */
		FLT wtx(int pnt) {return(_wtx(pnt));}
		FLT xp(int pnt) {return(_xp(pnt));}
		FLT x0(int pnt) {return(_x0(pnt));}

		/* COMBINED THINGS FOR FAST INTEGRATION */
		FLT gxwtx(int mod,int pnt) {return(_gxwtx(mod,pnt));}
		FLT dgxwtx(int mod,int pnt) {return(_dgxwtx(mod,pnt));}

		/* TO TAKE X,Y DERIVATIVES OF A FUNCTION WITH VALUES ON GAUSS POINTS */
		FLT dltx(int pnt) {return(_dltx(pnt));}
		FLT dltx1(int pnt1,int pnt2) {return(_dltx1(pnt1,pnt2));}
		
		/* ETA FUNCTIONS & DERIVATIVES */
		FLT gn(int pnt,int mod) {return(_gn(pnt,mod));}
		FLT dgn(int pnt,int mod) {return(_dgn(pnt,mod));}
		/* GAUSS WEIGTS & LOCATIONS */
		FLT wtn(int pnt) {return(_wtn(pnt));} 
		FLT np(int pnt) {return(_np(pnt));} 
		FLT n0(int pnt) {return(_n0(pnt));} 

		/* COMBINED THINGS FOR FAST INTEGRATION */
		FLT gnwtn(int pnt,int mod) {return(_gnwtn(pnt,mod));}
		FLT gnwtnn0(int pnt,int mod) {return(_gnwtnn0(pnt,mod));}
		FLT dgnwtn(int pnt,int mod) {return(_dgnwtn(pnt,mod));}
		
		/* TO TAKE X,Y DERIVATIVES OF A FUNCTION WITH VALUES ON GAUSS POINTS */
		FLT dltn(int pnt) {return(_dltn(pnt));}
		FLT dltn1(int pntx1,int pntx2) {return(_dltn1(pntx1,pntx2));}
		FLT dltn2(int pntn1,int pntn2) {return(_dltn2(pntn1,pntn2));}

		/* RENORMALIZATION CONSTANTS */
		FLT norm(int mod) {return(_norm(mod));}
		
		/* FOR OUTPUTING TO LEGENDRE POINTS */			  
		FLT lgrnge1d(int mod,int pnt) {return(_lgrnge1d(mod,pnt));}
		FLT lgrnge(int mod,int pntx,int pnty) {return(_lgrnge(mod,pntx,pnty));}
		
		/* FOR CALCULATING NORMAL DERIVATIVES ALONG SIDES */
		FLT dgnorm(int sd,int mod,int pnt) {return(_dgnorm(sd,mod,pnt));}

		/* LUMPED MASS MATRIX INVERSION */
		/* REMOVE SIDE COMPONENT FROM VERTICES */
		const FLT& sfmv(int v,int s) {return(_sfmv(v,s));}
		/* REMOVE INTERIOR COMPONENT FROM VERTICES & SIDES */
		const FLT& ifmb(int b,int i) {return(_ifmb(b,i));}
		/* DIAGONAL FOR VERTICES AFTER REMOVAL OF SIDES & INTERIORS */
		const FLT& vdiag() {return(_vdiag);}
		/* AFTER FINDING VERTEX REMOVE FROM SIDE MODES */
		const FLT& vfms(int v,int s) {return(_vfms(v,s));}
		/* DIAGONAL FOR SIDES */		
		const FLT& sdiag(int s) {return(_sdiag(s));}
		/* AFTER FINDING LOW ORDER SIDES REMOVE FROM HIGH ORDER SIDE MODES */
		const FLT& sfms(int ml,int mh,int sd) {return(_sfms(ml,mh,sd));}
		/* AFTER FINDING VERTICES & ALL SIDES REMOVE FROM INTERIOR */
		const FLT& bfmi(int b,int i) {return(_bfmi(b,i));}
		/* INTERIOR DIAGONAL */
		const FLT& idiag(int i,int band) {return(_idiag(i,band));}
		/* MASS MATRIX WITH STATIC INVERSION OF INTERIOR MODES */
		const FLT& msi(int b1,int b2) {return(_msi(b1,b2));}
		
		/*1D (SIDE) MASS MATRIX INVERSION */
		/* REMOVE SIDE COMPONENT FROM VERTICES */
		const FLT& sfmv1d(int v,int s) {return(_sfmv1d(v,s));}
		/* DIAGONAL FOR VERTICES AFTER REMOVAL OF SIDES */
		const FLT& vdiag1d() {return(_vdiag1d);}
		/* AFTER FINDING VERTEX REMOVE FROM SIDE MODES */
		const FLT& vfms1d(int v,int s) {return(_vfms1d(v,s));}
		/* DIAGONAL FOR SIDES */		
		const FLT& sdiag1d(int s,int band) {return(_sdiag1d(s,band));}	

		void initialize();
		/* PROJECT WITH R & S DERIVATIVES */
		void proj(const FLT *lin, FLT *f, FLT *dx, FLT *dy, int stride) const;
		/* PROJECT ONLY VALUE */
		void proj(const FLT *lin, FLT *f, int stride) const;
		/* PROJECT A LINEAR FUNCTION */
		void proj(FLT u1, FLT u2, FLT u3, FLT *f, int stride) const;
		/* PROJECT USING SIDE/VERTEX MODES WITH R & S DERIVATIVES */
		void proj_bdry(const FLT *lin, FLT *f, FLT *dr, FLT *ds, int stride) const;
		/* PROJECT USING SIDE/VERTEX MODES ONLY */
		void proj_bdry(const FLT *lin, FLT *f, int stride) const;
		/* PROJECT TO LEGENDRE POINTS (FOR OUTPUTING) */
		void proj_leg(const FLT *lin, FLT *f, int stride) const;
		/* PROJECT LINEAR FUNCTION TO LEGENDRE POINTS (FOR OUTPUTING) */
		void proj_leg(FLT u1, FLT u2, FLT u3, FLT *f, int stride) const;
		/* PROJECT TO LEGENDRE POINTS USING SIDE/VERTEX MODES ONLY (FOR OUTPUTING) */
		void proj_bdry_leg(const FLT *lin, FLT *f, int stride) const;
		/* PROJECT VALUES & TANGENT & NORMAL DERIVATIVES TO 1D SIDE GAUSS POINTS */
		/* dx is tangential derivative, dn is normal derivative to side */
		void proj_side(int side,const FLT *lin, FLT *f, FLT *dx, FLT *dn) const;

		/* DERIVATIVE IN R */
		void derivr(const FLT *f, FLT *dr, int stride) const;
		/* DERIVATIVE IN S */
		void derivs(const FLT *f, FLT *ds, int stride) const;
		/* INTEGRATE WITH RESPECT TO BASIS */
		void intgrt(FLT *rslt,const FLT *f, int stride) const;
		/* INTEGRATE FX W/RSPCT TO DG/DR & FY W/RSPCT TO DG/DS */
		/* WARNING THIS ADDS INTEGRATION TO RESULT: RESULT IS NOT CLEARED FIRST */
		void intgrtrs(FLT *rslt,const FLT *dx,const FLT *dy, int stride) const;

		/* SAME STUFF EXCEPT 1D	   */	 
		/* PROJECT WITH X DERIVATIVES */
		void proj1d(const FLT *lin, FLT *f, FLT *dx) const;
		/* PROJECT ONLY VALUE */
		void proj1d(const FLT *lin, FLT *f) const;
		/* PROJECT A LINEAR FUNCTION */
		void proj1d(FLT u1, FLT u2, FLT *f) const;
		/* PROJECT TO LEGENDRE POINTS (FOR OUTPUTING) */
		void proj1d_leg(const FLT *lin, FLT *f) const;
		/* PROJECT TO LEGENDRE POINTS USING SIDE/VERTEX MODES ONLY (FOR OUTPUTING) */
		void proj1d_leg(FLT u1, FLT u2, FLT *f) const;
		/* DERIVATIVE IN X */
		void derivx1d(const FLT *f, FLT *dx) const;
		/* INTEGRATE WITH RESPECT TO BASIS */
		void intgrt1d(FLT *rslt,const FLT *f) const;
		/* INTEGRATE W/RSPCT TO DG/DX */
		void intgrtx1d(FLT *rslt,const FLT *f) const;
			  
		/* POINT PROBE IN STANDARD ELEMENT FOR VECTOR */
		void ptprobe(int nv, FLT *f, FLT r, FLT s,const FLT *lin, int stride) const {
			ptvalues(2.0*(1+r)/(1-s+10.*EPSILON) -1.0,s);
			ptprobe(nv, f, lin, stride);
		}
		void ptprobe(int nv, FLT *f,const FLT *lin, int stride) const;	// REUSES OLD R,S
		void ptprobe(int nv, FLT *f1, FLT *dx1, FLT *dy1, FLT r, FLT s, const FLT *lin1, int stride) const;
		
		/* POINT PROBE FUNCTIONS USING BOUNDARY MODES ONLY */
		void ptprobe_bdry(int nv, FLT *f, FLT r, FLT s, const FLT *lin, int stride) const {
			ptvalues(2.0*(1+r)/(1-s+10.*EPSILON) -1.0,s);
			ptprobe_bdry(nv, f, lin, stride);
		}
		void ptprobe_bdry(int nv, FLT *f, const FLT *lin, int stride) const; // REUSES OLD R,S ONLY BDRY MODES 
		void ptprobe_bdry(int nv, FLT *f, FLT *dx, FLT *dy, FLT r, FLT s, const FLT *lin, int stride) const; // BOUNDARY MODES ONLY CALC'S DERIVATIVES
		
		/* 1D SIDE PROBE FUNCTIONS */
		void ptprobe1d(int nv, FLT *f, FLT x, const FLT *sin, int stride) const {	   
			ptvalues1d(x);
			ptprobe1d(nv,f,sin,stride);
		}
		void ptprobe1d(int nv, FLT *f, const FLT *sin, int stride) const;	 // REUSES OLD VALUES OF X
		void ptprobe1d(int nv, FLT *f, FLT *dx, FLT x, const FLT *sin, int stride) const {		
			ptvalues1d_deriv(x);
			ptprobe1d(nv,f,dx,sin,stride);
		}
		void ptprobe1d(int nv, FLT *f, FLT *dx, const FLT *sin, int stride) const;

		/* TO CALCULATE BASIS FUNCTIONS & DERIVATIVES */
		void ptvalues(FLT xi, FLT s) const; // CALCULATES GX, gn VALUES AT A POINT
		void ptvalues_deriv(FLT xi, FLT s) const; // CALCULATES GX, DGX, GN, DGN AT A POINT
		void ptvalues_bdry(FLT xi, FLT s) const; // CALCULATES GX, gn VALUES AT A POINT (BDRY MODES ONLY)
		void ptvalues_deriv_bdry(FLT xi, FLT s) const; // CALCULATES GX, gn & DERIV VALUES AT A POINT (BDRY MODES ONLY)
		
		void ptvalues_rs(FLT r, FLT s)	const {ptvalues(2.0*(1+r)/(1-s+10.*EPSILON) -1.0,s);} // CALCULATES GX, gn VALUES AT A POINT
		void ptvalues_deriv_rs(FLT r, FLT s) const {ptvalues_deriv(2.0*(1+r)/(1-s+10.*EPSILON) -1.0,s);} // CALCULATES GX, DGX, GN, DGN AT A POINT
		void ptvalues_bdry_rs(FLT r, FLT s) const {ptvalues_bdry(2.0*(1+r)/(1-s+10.*EPSILON) -1.0,s);} // CALCULATES GX, gn VALUES AT A POINT (BDRY MODES ONLY)
		void ptvalues_deriv_bdry_rs(FLT r, FLT s) const {ptvalues_deriv_bdry(2.0*(1+r)/(1-s+10.*EPSILON) -1.0,s);} // CALCULATES GX, gn & DERIV VALUES AT A POINT (BDRY MODES ONLY)
		
		/* 1D SIDE BASIS FUNCTIONS & DERIVATIVES */
		void ptvalues1d(FLT x) const;
		void ptvalues1d_deriv(FLT x) const;
		
		/* Utility for switching between uniform legendre representation & this basis */
		void legtobasis(const FLT *data, FLT *coeff) const;
		
	/* LOCAL STORAGE/WORK */
	private:
		mutable TinyVector<FLT,_nmodx> pgx, dpgx;
		mutable TinyVector<FLT,_tm> pgn, dpgn; // FOR POINT PROBE
		
		/* SETUP FUNCTIONS */
		void initialize_values(); // SET UP THINGS FOR PROJECT/INTEGRATE/DERIV
		void sideinfoinit(); // SET UP THINGS TO EVALUATE NORMAL DERIVATIVES ALONG SIDE
		void lumpinv(); // SET UP THINGS FOR INVERSE OF LUMPED MASS MATRIX
		void lumpinv1d();
		void legpt(); // SET UP PROJECTION TO LEGENDRE POINTS (FOR OUTPUTING)
};

/** This is an array for bases of various orders for general use 
	 The polynomial degree increases by factors of 2 for multigrid */
template<int EP> class tri_basis_array {
	private:
		TinyVector<tri_basis_interface *,3> tri;
	public: 
		tri_basis_array() {
			tri(0) = new tri_basis<1,EP>;
			tri(1) = new tri_basis<2,EP>;
			tri(2) = new tri_basis<4,EP>;	
			//tri(3) = new tri_basis<8,EP>;
		
			tri(0)->initialize();
			tri(1)->initialize();
			tri(2)->initialize();
			//tri(3)->initialize();
		}
		
		tri_basis_interface* operator()(int i) {
			return(tri(i));
		}
};
#endif



