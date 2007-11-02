/*
 *  tri_basis.h
 *  planar++
 *
 *  Created by helenbrk on Fri Oct 12 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
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

#define MAXP 4
#define MXTM (MAXP+2)*(MAXP+1)/2
#define MXGP (MAXP+2)
#define MORTHOGONAL

using namespace blitz;



class tri_basis {
    public:
        int p;
        /* NUMBER OF MODES : VERTEX,SIDE/INTERIOR/TOTAL/BOUNDARY */
        int sm,im,tm,bm;
        /* NUMBER OF X/S MODES & GAUSS POINTS */        
        int nmodx,nmodn,gpx,gpn;
        /* BANDWITH OF INTERIOR/INTERIOR MATRIX */
        int ibwth;
        /* BANDWITH OF SIDE/SIDE MATRIX (1D) */
        static const int sbwth = 2;
        /* RECURSION RELATION COEFFICIENTS */
        Array<FLT,2> a0, b0;
        
        /* X FUNCTIONS & DERIVATIVES */
        Array<FLT,2> gx, dgx;
        /* GAUSS WEIGTS & LOCATIONS */
        Array<FLT,1> wtx, xp, x0;
        /* COMBINED THINGS FOR FAST INTEGRATION */
        Array<FLT,2> gxwtx, dgxwtx;
        /* TO TAKE X,Y DERIVATIVES OF A FUNCTION WITH VALUES ON GAUSS POINTS */
        Array<FLT,1> dltx;
        Array<FLT,2> dltx1;
        
        /* ETA FUNCTIONS & DERIVATIVES */
        Array<FLT,2> gn, dgn;
        /* GAUSS WEIGTS & LOCATIONS */
        Array<FLT,1> wtn, np, n0;
        /* COMBINED THINGS FOR FAST INTEGRATION */
        Array<FLT,2> gnwtn, gnwtnn0, dgnwtn;
        /* TO TAKE X,Y DERIVATIVES OF A FUNCTION WITH VALUES ON GAUSS POINTS */
        Array<FLT,1> dltn;
        Array<FLT,2> dltn1, dltn2;

        /* RENORMALIZATION CONSTANTS */
        Array<FLT,1> norm;
        
        /* FOR OUTPUTING TO LEGENDRE POINTS */            
        Array<FLT,2> lgrnge1d;
        Array<FLT,3> lgrnge;
        
        /* FOR CALCULATING NORMAL DERIVATIVES ALONG SIDES */
        Array<FLT,3> dgnorm;

        /* LUMPED MASS MATRIX INVERSION */
        /* REMOVE SIDE COMPONENT FROM VERTICES */
        Array<FLT,2> sfmv;
        /* REMOVE INTERIOR COMPONENT FROM VERTICES & SIDES */
        Array<FLT,2> ifmb;
        /* DIAGONAL FOR VERTICES AFTER REMOVAL OF SIDES & INTERIORS */
        FLT vdiag;
        /* AFTER FINDING VERTEX REMOVE FROM SIDE MODES */
        Array<FLT,2> vfms;
        /* DIAGONAL FOR SIDES */        
        Array<FLT,1> sdiag;
        /* AFTER FINDING LOW ORDER SIDES REMOVE FROM HIGH ORDER SIDE MODES */
        Array<FLT,3> sfms;
        /* AFTER FINDING VERTICES & ALL SIDES REMOVE FROM INTERIOR */
        Array<FLT,2> bfmi;
        /* INTERIOR DIAGONAL */
        Array<FLT,2> idiag;
        /* MASS MATRIX WITH STATIC INVERSION OF INTERIOR MODES */
        Array<FLT,2> msi;
        
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
        void proj(FLT *lin, FLT *f, FLT *dx, FLT *dy, int stride);
        /* PROJECT ONLY VALUE */
        void proj(FLT *lin, FLT *f, int stride);
        /* PROJECT A LINEAR FUNCTION */
        void proj(FLT u1, FLT u2, FLT u3, FLT *f, int stride);
        /* PROJECT USING SIDE/VERTEX MODES WITH R & S DERIVATIVES */
        void proj_bdry(FLT *lin, FLT *f, FLT *dr, FLT *ds, int stride);
        /* PROJECT USING SIDE/VERTEX MODES ONLY */
        void proj_bdry(FLT *lin, FLT *f, int stride);
        /* PROJECT TO LEGENDRE POINTS (FOR OUTPUTING) */
        void proj_leg(FLT *lin, FLT *f, int stride);
        /* PROJECT LINEAR FUNCTION TO LEGENDRE POINTS (FOR OUTPUTING) */
        void proj_leg(FLT u1, FLT u2, FLT u3, FLT *f, int stride);
        /* PROJECT TO LEGENDRE POINTS USING SIDE/VERTEX MODES ONLY (FOR OUTPUTING) */
        void proj_bdry_leg(FLT *lin, FLT *f, int stride);
        /* PROJECT VALUES & TANGENT & NORMAL DERIVATIVES TO 1D SIDE GAUSS POINTS */
        /* dx is tangential derivative, dn is normal derivative to side */
        void proj_side(int side, FLT *lin, FLT *f, FLT *dx, FLT *dn);

        /* DERIVATIVE IN R */
        void derivr(FLT *f, FLT *dr, int stride);
        /* DERIVATIVE IN S */
        void derivs(FLT *f, FLT *ds, int stride);
        /* INTEGRATE WITH RESPECT TO BASIS */
        void intgrt(FLT *rslt, FLT *f, int stride);
        /* INTEGRATE FX W/RSPCT TO DG/DR & FY W/RSPCT TO DG/DS */
        /* WARNING THIS ADDS INTEGRATION TO RESULT: RESULT IS NOT CLEARED FIRST */
        void intgrtrs(FLT *rslt, FLT *dx, FLT *dy, int stride);

        /* SAME STUFF EXCEPT 1D    */    
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
        /* INTEGRATE WITH RESPECT TO BASIS */
        void intgrt1d(FLT *rslt, FLT *f);
        /* INTEGRATE W/RSPCT TO DG/DX */
        void intgrtx1d(FLT *rslt, FLT *f);
              
        /* POINT PROBE IN STANDARD ELEMENT FOR VECTOR */
        inline void ptprobe(int nv, FLT *f, FLT r, FLT s, FLT *lin, int stride) {
            ptvalues(2.0*(1+r)/(1-s+10.*EPSILON) -1.0,s);
            ptprobe(nv, f, lin, stride);
        }
        void ptprobe(int nv, FLT *f, FLT *lin, int stride);  // REUSES OLD R,S
        
        /* POINT PROBE FUNCTIONS USING BOUNDARY MODES ONLY */
        inline void ptprobe_bdry(int nv, FLT *f, FLT r, FLT s, FLT *lin, int stride) {
            ptvalues(2.0*(1+r)/(1-s+10.*EPSILON) -1.0,s);
            ptprobe_bdry(nv, f, lin, stride);
        }
        void ptprobe_bdry(int nv, FLT *f, FLT *lin, int stride); // REUSES OLD R,S ONLY BDRY MODES 
        void ptprobe_bdry(int nv, FLT *f, FLT *dx, FLT *dy, FLT r, FLT s, FLT *lin, int stride); // BOUNDARY MODES ONLY CALC'S DERIVATIVES
        
        /* 1D SIDE PROBE FUNCTIONS */
        inline void ptprobe1d(int nv, FLT *f, FLT x, FLT *sin, int stride) {     
            ptvalues1d(x);
            ptprobe1d(nv,f,sin,stride);
        }
        void ptprobe1d(int nv, FLT *f, FLT *sin, int stride);  // REUSES OLD VALUES OF X
        inline void ptprobe1d(int nv, FLT *f, FLT *dx, FLT x, FLT *sin, int stride) {     
            ptvalues1d_deriv(x);
            ptprobe1d(nv,f,dx,sin,stride);
        }
        void ptprobe1d(int nv, FLT *f, FLT *dx, FLT *sin, int stride);
        
        /* TO CALCULATE BASIS FUNCTIONS & DERIVATIVES */
        void ptvalues(FLT xi, FLT s); // CALCULATES GX, gn VALUES AT A POINT
        void ptvalues_deriv(FLT xi, FLT s); // CALCULATES GX, DGX, GN, DGN AT A POINT
        void ptvalues_bdry(FLT xi, FLT s); // CALCULATES GX, gn VALUES AT A POINT (BDRY MODES ONLY)
        void ptvalues_deriv_bdry(FLT xi, FLT s); // CALCULATES GX, gn & DERIV VALUES AT A POINT (BDRY MODES ONLY)
        
        void ptvalues_rs(FLT r, FLT s) {ptvalues(2.0*(1+r)/(1-s+10.*EPSILON) -1.0,s);} // CALCULATES GX, gn VALUES AT A POINT
        void ptvalues_deriv_rs(FLT r, FLT s) {ptvalues_deriv(2.0*(1+r)/(1-s+10.*EPSILON) -1.0,s);} // CALCULATES GX, DGX, GN, DGN AT A POINT
        void ptvalues_bdry_rs(FLT r, FLT s) {ptvalues_bdry(2.0*(1+r)/(1-s+10.*EPSILON) -1.0,s);} // CALCULATES GX, gn VALUES AT A POINT (BDRY MODES ONLY)
        void ptvalues_deriv_bdry_rs(FLT r, FLT s) {ptvalues_deriv_bdry(2.0*(1+r)/(1-s+10.*EPSILON) -1.0,s);} // CALCULATES GX, gn & DERIV VALUES AT A POINT (BDRY MODES ONLY)
        
        /* 1D SIDE BASIS FUNCTIONS & DERIVATIVES */
        void ptvalues1d(FLT x);
        void ptvalues1d_deriv(FLT x);
        
    /* LOCAL STORAGE/WORK */
    private:
        Array<FLT,1> pgx, dpgx, pgn, dpgn; // FOR POINT PROBE
        
        /* SETUP FUNCTIONS */
        void initialize_values(); // SET UP THINGS FOR PROJECT/INTEGRATE/DERIV
        void sideinfoinit(); // SET UP THINGS TO EVALUATE NORMAL DERIVATIVES ALONG SIDE
        void lumpinv(); // SET UP THINGS FOR INVERSE OF LUMPED MASS MATRIX
        void legpt(); // SET UP PROJECTION TO LEGENDRE POINTS (FOR OUTPUTING)
};

/** This is an array for bases of various orders for general use 
     The polynomial degree increases by factors of 2 for multigrid */
namespace basis {
    extern Array<tri_basis,1> tri;
}
#endif



