/*
 *  hpbasis.h
 *  planar++
 *
 *  Created by helenbrk on Fri Oct 12 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#define FLT double

#define MXTM 100

class hpbasis {
   public:
      int p;
/*		NUMBER OF MODES : SIDE/INTERIOR/TOTAL/BOUNDARY */
      int sm,im,tm,bm;
/*		NUMBER OF X/S MODES & GAUSS POINTS */      
      int nmodx,nmodn,gpx,gpn;
/*		BANDWITH OF INTERIOR/INTERIOR MATRIX */
      int ibwth;
/*		BANDWITH OF SIDE/SIDE MATRIX (1D) */
      static const int sbwth = 2;
/*		RECURSION RELATION COEFFICIENTS */
      FLT **a0, **b0;
      
/*		X FUNCTIONS & DERIVATIVES */
      FLT **gx, **dgx;
/*		GAUSS WEIGTS & LOCATIONS */
      FLT *wtx, *x0;
/*		COMBINED THINGS FOR FAST INTEGRATION */
      FLT **gxwtx, **dgxwtx;
/*		TO TAKE X,Y DERIVATIVES OF A FUNCTION WITH VALUES ON GAUSS POINTS */
      FLT *dltx, **dltx1;
      
/*		Y FUNCTIONS & DERIVATIVES */
      FLT **gn, **dgn;
/*		GAUSS WEIGTS & LOCATIONS */
      FLT *wtn, *n0;
/*		COMBINED THINGS FOR FAST INTEGRATION */
      FLT **gnwtnn0,**dgnwtn;
/*		TO TAKE X,Y DERIVATIVES OF A FUNCTION WITH VALUES ON GAUSS POINTS */
      FLT *dltn, **dltn1, **dltn2;

/*		RENORMALIZATION CONSTANTS */
      FLT *norm;
      
/*		FOR OUTPUTING TO LEGENDRE POINTS */         
      FLT **lgrnge1d, ***lgrnge;

/* 	LUMPED MASS MATRIX INVERSION */
/*		REMOVE SIDE COMPONENT FROM VERTICES */
      FLT **sfmv;
/*		REMOVE INTERIOR COMPONENT FROM VERTICES & SIDES */
      FLT **ifmb;
/* 	DIAGONAL FOR VERTICES AFTER REMOVAL OF SIDES & INTERIORS */
      FLT vdiag;
/*		AFTER FINDING VERTEX REMOVE FROM SIDE MODES */
      FLT **vfms;
/*		DIAGONAL FOR SIDES */      
      FLT *sdiag;
/*		AFTER FINDING LOW ORDER SIDES REMOVE FROM HIGH ORDER SIDE MODES */
      FLT ***sfms;
/*		AFTER FINDING VERTICES & ALL SIDES REMOVE FROM INTERIOR */
      FLT **bfmi;
/*		INTERIOR DIAGONAL */
      FLT **idiag;
      
/*		1D (SIDE) MASS MATRIX INVERSION */
/*		REMOVE SIDE COMPONENT FROM VERTICES */
      FLT **sfmv1d;
/* 	DIAGONAL FOR VERTICES AFTER REMOVAL OF SIDES */
      FLT vdiag1d;
/*		AFTER FINDING VERTEX REMOVE FROM SIDE MODES */
      FLT **vfms1d;
/*		DIAGONAL FOR SIDES */      
      FLT **sdiag1d;
         
      void initialize(int pdegree);
/*		PROJECT WITH R & S DERIVATIVES */
      void proj(FLT *lin, FLT **f, FLT **dr, FLT **ds);
/*		PROJECT ONLY VALUE */
      void proj(FLT *lin, FLT **f);
/*		PROJECT A LINEAR FUNCTION */
      void proj(FLT u1, FLT u2, FLT u3, FLT **f);
/*		PROJECT USING SIDE/VERTEX MODES WITH R & S DERIVATIVES */
      void proj_bdry(FLT *lin, FLT **f, FLT **dr, FLT **ds);
/*		PROJECT USING SIDE/VERTEX MODES ONLY VALUE */
      void proj_bdry(FLT *lin, FLT **f);
/*		PROJECT TO LEGENDRE POINTS (FOR OUTPUTING) */
      void proj_leg(FLT *lin, FLT **f);
/*		PROJECT LINEAR FUNCTION TO LEGENDRE POINTS (FOR OUTPUTING) */
      void proj_leg(FLT u1, FLT u2, FLT u3, FLT **f);
/*		PROJECT TO LEGENDRE POINTS USING SIDE/VERTEX MODES ONLY (FOR OUTPUTING) */
      void proj_bdry_leg(FLT *lin, FLT **f);

/*		DERIVATIVE IN R */
      void derivr(FLT **f, FLT **dr);
/*		DERIVATIVE IN S */
      void derivs(FLT **f, FLT **ds);
/*		INTEGRATE WITH RESPECT TO BASIS */
/* 	WARNING THESE ADD INTEGRATION TO RESULT: RESULT IS NOT CLEARED FIRST */
      void intgrt(FLT **f, FLT *rslt);
/*		INTEGRATE FX W/RSPCT TO DG/DR & FY W/RSPCT TO DG/DS */
      void intgrtrs(FLT **fr, FLT **fs, FLT *rslt);
/*		INTEGRATE W/RSPCT TO DG/DR */
      void intgrtr(FLT **f, FLT *rslt1);
/*		INTEGRATE W/RSPCT TO DG/DS */
      void intgrts(FLT **f, FLT *rslt2);

/*		SAME STUFF EXCEPT 1D   */   
/*		PROJECT WITH X DERIVATIVES */
      void proj1d(FLT *lin, FLT *f, FLT *dx);
/*		PROJECT ONLY VALUE */
      void proj1d(FLT *lin, FLT *f);
/*		PROJECT A LINEAR FUNCTION */
      void proj1d(FLT u1, FLT u2, FLT *f);
/*		PROJECT TO LEGENDRE POINTS (FOR OUTPUTING) */
      void proj1d_leg(FLT *lin, FLT *f);
/*		PROJECT TO LEGENDRE POINTS USING SIDE/VERTEX MODES ONLY (FOR OUTPUTING) */
      void proj1d_leg(FLT u1, FLT u2, FLT *f);
/*		DERIVATIVE IN X */
      void derivx1d(FLT *f, FLT *dx);
/*		INTEGRATE WITH RESPECT TO BASIS */
      void intgrt1d(FLT *f, FLT *rslt);
/*		INTEGRATE W/RSPCT TO DG/DX */
      void intgrtx1d(FLT *f, FLT *rslt);
      
/*		POINT PROBE IN STANDARD ELEMENT FOR VECTOR */
      void ptprobe(int nv, FLT **lin, FLT *f, FLT r, FLT s, FLT *pgx, FLT *pgn);
      void ptprobe_bdry(int nv, FLT **lin, FLT *f, FLT *dx, FLT *dy, FLT r, FLT s, FLT **pgx, FLT **pgn);
      void ptprobe1d(int nv, FLT **lin, FLT *f, FLT x, FLT *pgx);
      
/*	LOCAL STORAGE/WORK */
   private:
      static int wkpmax;
      static FLT **wk0,**wk1,**wk2,**wk3;
      
/*		SETUP FUNCTIONS */
      void initialize_values(); // SET UP THINGS FOR PROJECT/INTEGRATE/DERIV
      void lumpinv(); // SET UP THINGS FOR INVERSE OF LUMPED MASS MATRIX
      void legpt(); // SET UP PROJECTION TO LEGENDRE POINTS (FOR OUTPUTING)
      void ptvalues(FLT *pgx, FLT *pgn, FLT r, FLT s); // CALCULATES GX, gn VALUES AT A POINT
      void ptvalues_bdry(FLT **pgx, FLT **pgn, FLT r, FLT s); // CALCULATES GX, gn & DERIV VALUES AT A POINT (BDRY MODES ONLY)
      void ptvalues1d(FLT *pgx, FLT x);
};

#if (FLT == double)
#define GETRF DGETRF
#define GETRS DGETRS
#define PBTRF DPBTRF
#define PBTRS DPBTRS
#define EPSILON DBL_EPSILON
#else
#define GETRF SGETRF
#define GETRS SGETRS
#define PBTRF SPBTRF
#define PBTRS SPBTRS
#define EPSILON FLT_EPSILON
#endif
