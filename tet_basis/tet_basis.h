/*
 *  tet_basis.h
 *
 *  Created by helenbrk on Fri Oct 12 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _tet_basis_h_
#define _tet_basis_h_

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

#define MAXP 4
#define MXTM (MAXP+1)*(MAXP+2)*(MAXP+3)/6
#define MXGP (MAXP+1)
#define MORTHOGONAL

using namespace blitz;

class tet_basis {
	public:
		int p;
		/* NUMBER OF MODES: VERTEX,EDGE,FACE,INTERIOR,BOUNDARY,TOTAL*/
		int vm,em,fm,im,bm,tm;
		/* NUMBER OF X/Y/Z MODES & GAUSS POINTS*/
		int nmodx,nmody,nmodz,gpx,gpy,gpz;
		/* BANDWITH OF INTERIOR/INTERIOR MATRIX */
		int ibwth;
		/* BANDWITH OF SIDE/SIDE MATRIX (1D) */
		static const int sbwth = 2;
		/* RECURSION RELATION COEFFICIENTS */
		Array<FLT,2> a0, b0, a1, b1, a2, b2;      
		/* X FUNCTIONS & DERIVATIVES */
		Array<FLT,2> gx, dgx;
		/* GAUSS WEIGTS & LOCATIONS */
		Array<FLT,1> wtx, xp, x0;
		/* COMBINED THINGS FOR FAST INTEGRATION */
		Array<FLT,2> gxwtx, dgxwtx;
		/* TO TAKE X,Y DERIVATIVES OF A FUNCTION WITH VALUES ON GAUSS POINTS */
		Array<FLT,1> dltx;
		Array<FLT,2> dltx2,dltx3;

		/* Y FUNCTIONS & DERIVATIVES */
		Array<FLT,2> gy, dgy;
		/* GAUSS WEIGTS & LOCATIONS */
		Array<FLT,1> wty, yp, y0, y1;
		/* COMBINED THINGS FOR FAST INTEGRATION */
		Array<FLT,2> gywty, gywtyy0, dgywty;
		/* TO TAKE X,Y DERIVATIVES OF A FUNCTION WITH VALUES ON GAUSS POINTS */
		Array<FLT,1> dlty;
		Array<FLT,2> dlty2, dlty3;

		/* Z FUNCTIONS & DERIVATIVES */
		Array<FLT,2> gz, dgz;
		/* GAUSS WEIGTS & LOCATIONS */
		Array<FLT,1> wtz, zp, z0;
		/* COMBINED THINGS FOR FAST INTEGRATION */
		Array<FLT,2> gzwtz, gzwtzz0, dgzwtz;
		/* TO TAKE X,Y DERIVATIVES OF A FUNCTION WITH VALUES ON GAUSS POINTS */
		Array<FLT,1> dltz;
		Array<FLT,2> dltz2;	

		/* RENORMALIZATION CONSTANTS */
		Array<FLT,1> norm;

		/* FOR OUTPUTING TO LAGRANGE POINTS */         
		Array<FLT,2> lgrnge1d;
		Array<FLT,3> lgrnge2d;
		Array<FLT,4> lgrnge3d;

		/* FOR CALCULATING NORMAL DERIVATIVES ALONG FACE 0 */
		Array<FLT,3> dgnorm;

		/* LUMPED MASS MATRIX INVERSION */
		Array<FLT,2> mm;
		/* REMOVE SIDE and FACE COMPONENT FROM VERTICES */
		Array<FLT,2> sfmv;
		Array<FLT,2> ffmv;
		/* REMOVE INTERIOR COMPONENT FROM VERTICES, SIDES, & FACES */
		Array<FLT,2> ifmb;
		/* DIAGONAL FOR VERTICES AFTER REMOVAL OF SIDES & INTERIORS */
		/* AFTER FINDING VERTEX REMOVE FROM SIDE MODES */
		Array<FLT,2> vfms;
		/* AFTER FINDING LOW ORDER SIDES REMOVE FROM HIGH ORDER SIDE MODES */
		Array<FLT,1> sfms;//temp
		Array<FLT,3> ffms;//temp
		/* AFTER FINDING VERTICES & ALL SIDES REMOVE FROM INTERIOR */
		Array<FLT,2> bfmi;
		/* EDGE DIAGONAL */
		Array<FLT,1> diag1d;
		/* FACE DIAGONAL */
		Array<FLT,1> diag2d;
		/* TET DIAGONAL */
		Array<FLT,1> diag3d;
		/* MASS MATRIX WITH STATIC INVERSION OF INTERIOR MODES */
		Array<FLT,2> msi;
		/* Diagonal mass matrix */
		FLT vdiag,mdiag,odiag;
		Array<FLT,1> ediag,fdiag,idiag;

		/*1D (SIDE) MASS MATRIX INVERSION */
		/* REMOVE SIDE COMPONENT FROM VERTICES */
		Array<FLT,2> sfmv1d;
		/* DIAGONAL FOR VERTICES AFTER REMOVAL OF SIDES */
		FLT vdiag1d;
		/* AFTER FINDING VERTEX REMOVE FROM SIDE MODES */
		Array<FLT,2> vfms1d;
		/* DIAGONAL FOR SIDES */      
		Array<FLT,2> sdiag1d;

		void initialize(int pdegree, int gpoints);
		inline void initialize(int pdegree) { initialize(pdegree, pdegree+1);}
		static inline int smode(int p1) {return(p1-1);}
		static inline int imode(int p1) {return((p1-2)*(p1-1)/2);}
		static inline int tmode(int p1) {return((p1+2)*(p1+1)/2);}
		/* PROJECT WITH R & S DERIVATIVES */
		void proj(FLT *lin1, FLT *f1, FLT *dx, FLT *dy, FLT *dz, int stridex, int stridey);
		/* PROJECT ONLY VALUE */
    void proj(FLT *lin1, FLT *f1, int stridex, int stridey);
		/* PROJECT A LINEAR FUNCTION */
		void proj(FLT u1, FLT u2, FLT u3, FLT u4, FLT *f, int stridex, int stridey);
		/* PROJECT USING VERTEX/EDGE/FACE MODES WITH R & S DERIVATIVES */
		void proj_bdry(FLT *lin, FLT *f, FLT *dx, FLT *dy, FLT *dz, int stridex, int  stridey);
		/* PROJECT USING VERTEX/EDGE/FACE MODES ONLY */
		void proj_bdry(FLT *lin, FLT *f, int stridex, int stridey);
		/* PROJECT TO GL POINTS (FOR OUTPUTING) */
		void proj_leg(FLT *lin, FLT *f, int stridex, int stridey);
		/* PROJECT LINEAR FUNCTION TO LEGENDRE POINTS (FOR OUTPUTING) */
		void proj_leg(FLT u1, FLT u2, FLT u3, FLT u4, FLT *f, int stridex, int stridey);
		/* PROJECT TO LEGENDRE POINTS USING VERTEX/EDGE/FACE MODES ONLY (FOR OUTPUTING) */
		void proj_bdry_leg(FLT *lin, FLT *f, int stridex, int stridey);
		/* PROJECT VALUES & TANGENT & NORMAL DERIVATIVES TO 1D SIDE GAUSS POINTS */
		/* dt is tangential derivative, dn is normal derivative to face */
		void proj_face(FLT *lin, FLT *f, FLT *dt1, FLT *dt2, FLT *dn);

		/* DERIVATIVE IN R */
		void derivr(FLT *f, FLT *dr, int stridex, int stridey);
		/* DERIVATIVE IN S */
		void derivs(FLT *f, FLT *ds, int stridex, int stridey);
		/* DERIVATIVE IN T */
		void derivt(FLT *f, FLT *dt, int stridex, int stridey);

		/* INTEGRATE WITH RESPECT TO BASIS */
		void intgrt1d(FLT *rslt, FLT *f);
		void intgrt2d(FLT *rslt, FLT *f, int stride);
		void intgrt(FLT *rslt1, FLT *f1, int stridex, int stridey);
		/* INTEGRATE FX W/RSPCT TO DG/DR & FY W/RSPCT TO DG/DS */
		/* WARNING THIS ADDS INTEGRATION TO RESULT: RESULT IS NOT CLEARED FIRST */
		void intgrtrst(FLT *rslt, FLT *dx, FLT *dy, FLT *dz, int stridex, int stridey);

		void proj2d(FLT *lin1, FLT *f1, FLT *dx1, FLT *dy1, int stride); 
		void proj2d(FLT *lin1, FLT *f1, int stride); 
		void proj2d_bdry(FLT *lin1, FLT *f1, int stride); 
		void proj2d(FLT u1, FLT u2, FLT u3, FLT *f1, int stride);
		void proj2d_leg(FLT *lin1, FLT *f1, int stride);
		void proj2d_leg(FLT u1, FLT u2, FLT u3, FLT *f1, int stride);
		void proj2d_bdry_leg(FLT *lin1, FLT *f1, int stride); 

		/* SAME STUFF EXCEPT 1D   */   
		/* PROJECT WITH X DERIVATIVES */
		void proj1d(FLT *lin, FLT *f, FLT *dx);
		/* PROJECT ONLY VALUE */
		void proj1d(FLT *lin, FLT *f);
		/* PROJECT A LINEAR FUNCTION */
		void proj1d(FLT u1, FLT u2, FLT *f);
		/* PROJECT TO LEGENDRE POINTS (FOR OUTPUTING) */
		void proj1d_leg(FLT *lin, FLT *f);
		/* PROJECT TO LEGENDRE POINTS USING SIDE/VERTEX MODES ONLY (FOR OUTPUTING) */
		void proj1d_leg(FLT u1, FLT u2, FLT *f);
		/* DERIVATIVE IN X */
		void derivx1d(FLT *f, FLT *dx);
		/* INTEGRATE W/RSPCT TO DG/DX */
		void intgrtx1d(FLT *rslt, FLT *f);
		   
		/* POINT PROBE IN STANDARD ELEMENT FOR VECTOR */
		inline void ptprobe(int nv, FLT *f, FLT r, FLT s, FLT t, FLT *lin, int stride) {
			ptvalues(2.0*(1+r)/(-s-t+10.*EPSILON) -1.0, 2.0*(1+s)/(1-t+10.*EPSILON)-1.0, t);
			ptprobe(nv, f, lin, stride);
		}
      
      void ptprobe(int nv, FLT *f1, FLT *dx1, FLT *dy1, FLT *dz1, FLT r, FLT s, FLT t, FLT *lin1, int stride);
      
//        {
//         ptvalues_deriv(2.0*(1+r)/(-s-t+10.*EPSILON) -1.0, 2.0*(1+s)/(1-t+10.*EPSILON)-1.0, t);
//         ptprobe(nv, f1, dx1, dy1, dz1, lin1, stride);
//      }
		
		/* 2D SIDE PROBE FUNCTIONS */
		inline void ptprobe2d(int nv, FLT *f, FLT r, FLT s, FLT *sin, int stride) {    
			ptvalues2d(2.0*(1+r)/(1-s+10.*EPSILON) -1.0,s);
			ptprobe2d(nv,f,sin,stride);
		}

		/* 1D SIDE PROBE FUNCTIONS */
		inline void ptprobe1d(int nv, FLT *f, FLT x, FLT *sin, int stride) {    
			ptvalues1d(x);
			ptprobe1d(nv,f,sin,stride);
		}	 
		 
		void ptprobe(int nv, FLT *f, FLT *lin, int stride);  // REUSES OLD R,S,T
//      void ptprobe(int nv, FLT *f1, FLT *dx1, FLT *dy1, FLT *dz1, FLT *lin, int stride); // REUSES OLD R,S,T
		void ptprobe2d(int nv, FLT *f, FLT *lin, int stride); // REUSES OLD R,S  
		void ptprobe1d(int nv, FLT *f, FLT *sin, int stride);  // REUSES OLD VALUES OF X

		/* POINT PROBE FUNCTIONS USING BOUNDARY MODES ONLY */
		inline void ptprobe_bdry(int nv, FLT *f, FLT r, FLT s, FLT t, FLT *lin, int stride) {
			ptvalues(2.0*(1+r)/(-s-t+10.*EPSILON) -1.0, 2.0*(1+s)/(1-t+10.*EPSILON)-1.0, t);
			ptprobe_bdry(nv, f, lin, stride);
		}
		void ptprobe_bdry(int nv, FLT *f, FLT *lin, int stride); // REUSES OLD R,S ONLY BDRY MODES 
//      void ptprobe_bdry(int nv, FLT *f, FLT *dx, FLT *dy, FLT dz, FLT r, FLT s, FLT t, FLT *lin, int stride); // BOUNDARY MODES ONLY CALC'S DERIVATIVES

//      inline void ptprobe1d(int nv, FLT *f, FLT *dx, FLT x, FLT *sin, int stride) {    
//         ptvalues1d_deriv(x);
//         ptprobe1d(nv,f,dx,sin,stride);
//      }
//      void ptprobe1d(int nv, FLT *f, FLT *dx, FLT *sin, int stride);

		/* TO CALCULATE BASIS FUNCTIONS & DERIVATIVES */
		void ptvalues(FLT x, FLT y, FLT z); // CALCULATES GX, gn VALUES AT A POINT
		void ptvalues_deriv(FLT x, FLT y, FLT z); // CALCULATES GX, DGX, GN, DGN AT A POINT
		void ptvalues_bdry(FLT x, FLT y, FLT z); // CALCULATES GX, gn VALUES AT A POINT (BDRY MODES ONLY)
		void ptvalues_deriv_bdry(FLT x, FLT y, FLT z); // CALCULATES GX, gn & DERIV VALUES AT A POINT (BDRY MODES ONLY)

		void ptvalues_rst(FLT r, FLT s, FLT t) {ptvalues(2.0*(1+r)/(-s-t+10.*EPSILON) -1.0, 2.0*(1+s)/(1-t+10.*EPSILON)-1.0, t);} // CALCULATES GX, gn VALUES AT A POINT
		void ptvalues_deriv_rst(FLT r, FLT s, FLT t) {ptvalues_deriv(2.0*(1+r)/(-s-t+10.*EPSILON) -1.0, 2.0*(1+s)/(1-t+10.*EPSILON)-1.0, t);} // CALCULATES GX, DGX, GN, DGN AT A POINT
		void ptvalues_bdry_rst(FLT r, FLT s, FLT t) {ptvalues_bdry(2.0*(1+r)/(-s-t+10.*EPSILON) -1.0, 2.0*(1+s)/(1-t+10.*EPSILON)-1.0, t);} // CALCULATES GX, gn VALUES AT A POINT (BDRY MODES ONLY)
		void ptvalues_deriv_bdry_rst(FLT r, FLT s, FLT t) {ptvalues_deriv_bdry(2.0*(1+r)/(-s-t+10.*EPSILON) -1.0, 2.0*(1+s)/(1-t+10.*EPSILON)-1.0, t);} // CALCULATES GX, gn & DERIV VALUES AT A POINT (BDRY MODES ONLY)


		/* 2D SIDE BASIS FUNCTIONS & DERIVATIVES */
		void ptvalues2d(FLT x, FLT y);
		void ptvalues2d_deriv(FLT x, FLT y);

		void ptvalues2d_rst(FLT r, FLT s){ptvalues2d(2.0*(1+r)/(1-s+10.*EPSILON) -1.0,s);}
		void ptvalues2d_deriv_rst(FLT r, FLT s){ptvalues2d_deriv(2.0*(1+r)/(1-s+10.*EPSILON) -1.0,s);}
		/* 1D SIDE BASIS FUNCTIONS & DERIVATIVES */
		void ptvalues1d(FLT x);
		void ptvalues1d_deriv(FLT x);

		/* outputs assortment of stuff to check code*/
		void outputstuff(int check);
      
	/* LOCAL STORAGE/WORK */
	private:
		Array<FLT,1> pgx, dpgx, pgy, dpgy, pgz, dpgz; // FOR POINT PROBE
		/* SETUP FUNCTIONS */
		void initialize_values(); // SET UP THINGS FOR PROJECT/INTEGRATE/DERIV
		void faceinfoinit(); // SET UP THINGS TO EVALUATE NORMAL DERIVATIVES ALONG FACE
		void lumpinv(); // SET UP THINGS FOR INVERSE OF LUMPED MASS MATRIX
		void legpt(); // SET UP PROJECTION TO GL POINTS (FOR OUTPUTING)
};

/** This is an array for bases of various orders for general use 
    The polynomial degree increases by factors of 2 for multigrid */
namespace basis {
   extern Array<tet_basis,1> tet;
}

#endif



