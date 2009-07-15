/*
 *  initialize.cpp
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 01 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */
// #define DEBUG

#include <math.h>
#include <utilities.h>
#include <myblas.h>
#include "tri_basis.h"

const int tri_basis::sbwth;

Array<tri_basis,1> basis::tri;

#ifdef DEBUG
#include <blitz/tinyvec-et.h>
FLT func(FLT r,FLT s) {return(s*s*s*s);}
#endif


void tri_basis::initialize(int pdegree, int gpoints) {    
    if (pdegree < 1) {
        printf("error can't use 0th order basis with vertex based modes\n");
        exit(1);
    }
    
    p = pdegree;
    sm = MAX(p -1,0);
    im = (p > 0 ? (p-2)*(p-1)/2 : 0);
    bm = MAX(3*p,1);
    tm = bm + im;
    ibwth = (2*(sm-1) < im -1? 2*(sm-1) : im-1);
    ibwth = (ibwth > 0 ? ibwth : 0);
#ifdef MORTHOGONAL
    ibwth = 1;
#endif
    
    (p > 0 ? nmodx = p+2 : nmodx = 1);
    nmodn = tm;
    gpx = gpoints;
    gpn = gpoints;

    pgx.resize(nmodx);
    dpgx.resize(nmodx);
    pgn.resize(tm);
    dpgn.resize(tm);
    
    /*****************************/
    /* SETUP VALUES OF FUNCTIONS */
    /*****************************/
    initialize_values(); // SET UP THINGS FOR PROJECT/INTEGRATE/DERIV
    sideinfoinit(); // SET UP THINGS TO CALCULATE NORMAL DERIVATIVE TO SIDE ALONG SIDE
    lumpinv(); // SET UP THINGS FOR INVERSE OF LUMPED MASS MATRIX
    legpt(); // SET UP PROJECTION TO LEGENDRE POINTS (FOR OUTPUTING)
    
#ifdef DEBUG
    /* SOME TESTING */
    Array<double,1> uht(MXTM);
    Array<double,2> u(gpx,gpn);
    Array<double,2> u1(gpx,gpn);
    Array<double,2> dx(gpx,gpn);
    Array<double,2> dy(gpx,gpn);
    Array<double,2> dx1(gpx,gpn);
    Array<double,2> dy1(gpx,gpn);
    
    for(int m=0;m<tm;++m) {
        for(int j=0;j<tm;++j)
            uht(j) = 0.0;
        uht(m) = 1.0;
        
        proj_side(0,&uht(0),&u(0,0),&dx(0,0),&dy(0,0));
        for (int i = 0; i<gpx;++i)
            printf("T0A: %d %e %e %e\n",i,u(0,i),dx(0,i),dy(0,i));
        
        proj_side(1,&uht(0),&u(0,0),&dx(0,0),&dy(0,0));
        for (int i = 0; i<gpx;++i)
            printf("T0B: %d %e %e %e\n",i,u(0,i),dx(0,i),dy(0,i));
        
        proj_side(2,&uht(0),&u(0,0),&dx(0,0),&dy(0,0));
        for (int i = 0; i<gpx;++i)
            printf("T0C: %d %e %e %e\n",i,u(0,i),dx(0,i),dy(0,i));

        proj(&uht(0),&u(0,0),&dx(0,0),&dy(0,0),gpn);
        for(int i=0;i<gpx;++i) {
            for(int j=0;j<gpn;++j) {
                dx1(i,j) = 0.0;
                dy1(i,j) = 0.0;
            }
        }
        derivr(&u(0,0),&dx1(0,0),gpn);
        derivs(&u(0,0),&dy1(0,0),gpn);
        for(int i=0;i<gpx;++i)
            for(int j=0;j<gpn;++j)
                printf("T1: %d %d %e %e %e %e\n",i,j,dx(i,j),dx1(i,j),dy(i,j),dy1(i,j));
                

        if (m < 3) {
            proj(uht(0),uht(1),uht(2),&u1(0,0),gpn);
            for(int i=0;i<gpx;++i)
                for(int j=0;j<gpn;++j)
                    printf("T2: %d %d %e %e\n",i,j,u1(i,j),u(i,j));
        }
        
        if (m < bm) {
            proj_bdry(&uht(0),&u1(0,0),&dx1(0,0), &dy1(0,0),gpn);
            for(int i=0;i<gpx;++i)
                for(int j=0;j<gpn;++j)
                    printf("T3: %d %d %e %e %e %e %e %e\n",i,j,u1(i,j),u(i,j),dx(i,j),dx1(i,j),dy(i,j),dy1(i,j));
            
            proj_bdry(&uht(0),&u1(0,0),gpn);
            for(int i=0;i<gpx;++i)
                for(int j=0;j<gpn;++j)
                    printf("T4: %d %d %e %e\n",i,j,u1(i,j),u(i,j));
        }

        FLT val,valx,valy,val1,val2;
        ptprobe(1, &val, 0.25, 0.25, &uht(0), MXTM);
        ptprobe_bdry(1,&val1, 0.25, 0.25,&uht(0), MXTM);
        ptprobe_bdry(1, &val2, &valx, &valy, 0.25, 0.25,&uht(0),MXTM);
        printf("T5: %e %e %e\n",val,val1,val2);
        printf("T5: %e %e\n",valx,valy);
        
        intgrtrs(&uht(0),&dx(0,0),&dy(0,0),gpn);
        for(int j=0;j<tm;++j)
            printf("T6: %d %e\n",j,uht(j));
        
            
    }


    /* 1D TESTING */
    for (int i=0;i<gpx;++i)
        printf("dltx: %d %e\n",i,dltx(i));
        
    for (int i=0;i<gpx;++i) 
        for(int n=0;n<gpx;++n) 
            printf("dltx1: %d %d %e\n",i,n,dltx1(i,n));
            
            
    for(int m=0;m<sm+2;++m) {
        for(int j=0;j<sm+2;++j)
            uht(j) = 0.0;
        uht(m) = 1.0;
        
        proj1d(&uht(0),&u(0,0),&dx(0,0));
        for(int j=0;j<gpx;++j)
            dx1(0,j) = 0.0;
        derivx1d(&u(0,0),&dx1(0,0));
        for(int i=0;i<gpx;++i)
            printf("T11D: %d %e %e\n",i,dx(0,i),dx1(0,i));
        
        proj1d(&uht(0),&u1(0,0));
        for(int i=0;i<gpx;++i)
            printf("T21D: %d %e %e\n",i,u1(0,i),u(0,i));
            
        if (m < 2) {
            proj1d(uht(0),uht(1),&u1(0,0));
            for(int i=0;i<gpx;++i)
                printf("T31D: %d %e %e\n",i,u1(0,i),u(0,i));
        }
        
        FLT val,valx,val1;
        
        ptprobe1d(1,&val, 0.25,&uht(0),MXTM);
        ptprobe1d(1,&val1,&valx, 0.25, &uht(0),MXTM);
        printf("T41D: %e %e %e\n",val,val1,valx);
        
        
        intgrt1d(&uht(0),&u(0,0));
        for(int j=0;j<sm+2;++j)
            printf("T51D: %d %e\n",j,uht(j));
            
        intgrtx1d(&uht(0),&u(0,0));
        for(int j=0;j<sm+2;++j)
            printf("T61D: %d %e\n",j,uht(j));
    }
        
#endif
    
    return;
}

void tri_basis::initialize_values(void)
{

    FLT al,be,x,eta;
    int i,j,k,m,n,ipoly,ierr;
    Array<FLT,1> e(MXTM),e1(MXTM),e2(MXTM);

    /* ALLOCATE STORAGE FOR RECURSION RELATION COEFFICENTS*/
    ipoly = MAX(gpx+1,sm+1);
    ipoly = MAX(ipoly,gpn+1);
    a0.resize(sm+2,ipoly);
    b0.resize(sm+2,ipoly);
    
    /* ALLOCATE INTEGRATION, PROJECTION, & DERIVATIVE VARIABLES */
    gx.resize(gpx,nmodx);
    dgx.resize(gpx,nmodx);
    wtx.resize(gpx);
    xp.resize(gpx);
    x0.resize(gpx);
    gxwtx.resize(nmodx,gpx);
    dgxwtx.resize(nmodx,gpx);
    dltx.resize(gpx);
    dltx1.resize(gpx,gpx);
    
    gn.resize(gpn,tm);
    dgn.resize(gpn,tm);
    wtn.resize(gpn);
    np.resize(gpn);
    n0.resize(gpn);
    gnwtn.resize(tm,gpn);
    gnwtnn0.resize(tm,gpn);
    dgnwtn.resize(tm,gpn);
    dltn.resize(gpn);
    dltn1.resize(gpx,gpx);
    dltn2.resize(gpn,gpn);
    norm.resize(tm);
        
    /* GENERATE RECURSION RELATION FOR LEGENDRE
    POLYNOMIALS (USED TO GENERATE GAUSS POINTS)
    IPOLY = 1 SIGNIFIES LEGENDRE POLYNOMIALS
    NEED RECURSION RELATIONS FROM 1 TO N FOR
    N POINT GAUSS FORMULA

    RECURSION RELATION IS OF THE FORM:
    p(k+1)(x)=(x-a(k))*p(k)(x)-b(k)*p(k-1)(x),
    k=0,1,...,n-1,
    
    p(-1)(x)=0,  p(0)(x)=1
    */

    /******************************************/
    /* CALCULATE GAUSS POINTS / COEFFICIENTS  */
    /**************************************    */
    ipoly = 1;
    al = 0.0;
    be = 0.0;
    ierr = recur(gpx+1,ipoly,al,be,&a0(0,0),&b0(0,0));
    if (ierr != 0) {
        printf("recur #1 error %d\n",ierr);
        exit(1);
    }

    ierr = gauss(gpx,&a0(0,0),&b0(0,0),EPSILON,&x0(0),&wtx(0),&e(0));
    if (ierr != 0) {
        printf("gauss #1 error %d\n",ierr);
        exit(1);
    }

    /***********************************************/
    /* CALCULATE FACTORS FOR LAGRANGIAN DERIVATIVE */
    /***********************************************/
    for (i=0;i<gpx;++i) {
        dltx(i) = 1.0;
        for(k=0;k<i;++k)
            dltx(i) *= (x0(i)-x0(k));
        for(k=i+1;k<gpx;++k)
            dltx(i) *= (x0(i)-x0(k));
    }
    
    for (i=0;i<gpx;++i) {
        for(n=0;n<i;++n) 
            dltx1(i,n) = dltx(i)/(x0(i)-x0(n));
        for(n=i+1;n<gpx;++n) 
            dltx1(i,n) = dltx(i)/(x0(i)-x0(n));
    }

    for (i=0;i<gpx;++i) {
        for(n=0;n<i;++n) 
            dltn1(i,n) = dltx(i)/(x0(i)-x0(n))*.5*(1+x0(i));
        for(n=i+1;n<gpx;++n) 
            dltn1(i,n) = dltx(i)/(x0(i)-x0(n))*.5*(1+x0(i));
    }
        
    for (i=0;i<gpx;++i)
        dltx(i) = 1.0/dltx(i);    

    /*************************************************/ 
    /* NOW CALCULATE VALUES OF G, G' AT GAUSS POINTS */
    /*************************************************/

    /* SIDE 1 IPOLY = 6 FOR JACOBI POLYNOMIALS  */
    ipoly = 6;
#ifdef MORTHOGONAL
    al = 2.0;
    be = 2.0;
#else
    al = 1.0;
    be = 1.0;
#endif
    ierr = recur(sm+1,ipoly,al,be,&a0(0,0),&b0(0,0));
    if (ierr != 0) {
        printf("recur #3 error %d\n",ierr);
        exit(1);
    }

    for(i=0;i<tm;++i)
        norm(i) = 1.0;
    
    for(i = 0;i < gpx; ++i) {
        x = x0(i);
        xp(i) = x;
        x0(i) = 0.5*(1+x);
        
        ptvalues_deriv(x,0.0);

        for (m = 0;m < nmodx;++m) {
            gx(i,m) = pgx(m);
            dgx(i,m) = dpgx(m);
        }
    }

    /*********************************************/
    /* CALCULATE GAUSS POINTS / COEFFICIENTS (ETA=N) */
    /****************************************    */
    ipoly = 6;
    al = 1.0;
    be = 0.0;
#ifdef OLDWAY
    ierr = recur(gpn,ipoly,al,be,&a0(1,0),&b0(1,0));
    if (ierr != 0) {
        printf("recur #2 error %d\n",ierr);
        exit(1);
    }    
    
    ierr = radau(gpn-1,&a0(1,0),&b0(1,0),-1.0,&n0(0),&wtn(0),&e(0),&e1(0),&e2(0));        
    if (ierr != 0) {
        printf("gauss #3 error %d\n",ierr);
        exit(1);
    }
#else
    ierr = recur(gpn+1,ipoly,al,be,&a0(1,0),&b0(1,0));
    if (ierr != 0) {
        printf("recur #2 error %d\n",ierr);
        exit(1);
    }

    ierr = gauss(gpn,&a0(1,0),&b0(1,0),EPSILON,&n0(0),&wtn(0),&e(0));
    if (ierr != 0) {
        printf("gauss #1 error %d\n",ierr);
        exit(1);
    }
#endif
    
    for (i=0;i<gpn;++i)
        wtn(i) *= 0.5;

    /***********************************************/
    /* CALCULATE FACTORS FOR LAGRANGIAN DERIVATIVE */
    /***********************************************/
    for (i=0;i<gpn;++i) {
        dltn(i) = 1.0;
        for(k=0;k<i;++k)
            dltn(i) *= (n0(i)-n0(k));
        for(k=i+1;k<gpn;++k)
            dltn(i) *= (n0(i)-n0(k));
    }
    
    for (j=0;j<gpn;++j) {
        for(n=0;n<j;++n) 
            dltn2(j,n) = dltn(j)/(n0(j)-n0(n));
        for(n=j+1;n<gpn;++n) 
            dltn2(j,n) = dltn(j)/(n0(j)-n0(n));
    }
    
    for (j=0;j<gpn;++j)
        dltn(j) = 1.0/dltn(j);
    

    /******************************************/
    /* GENERATE JACOBI POLY FOR S DIRECTION */
    /****************************************/
    /* RECURSION RELATION FOR SIDE MODES */
    /* SIDE 2 IPOLY = 6 FOR JACOBI    */
    ipoly = 6;
#ifdef MORTHOGONAL
    al = 2.0;
    be = 2.0;
#else
    al = 1.0;
    be = 1.0;
#endif
    ierr = recur(sm+1,ipoly,al,be,&a0(1,0),&b0(1,0));
    if (ierr != 0) {
        printf("recur #3 error %d\n",ierr);
        exit(1);
    }

    /*	RECURSION RELATION FOR INTERIOR MODES */
    for(m = 2; m< sm+1;++m) {        
    /* CALCULATE RECURSION RELATION FOR P^(2m-1,1)(s) SIDE 1 IPOLY = 6 FOR JACOBI      */
        ipoly = 6;
#ifdef MORTHOGONAL
        al = 2.*m+1;
        be = 2.0;
#else
        al = 2.*m-1;
        be = 1.0;
#endif
        ierr = recur(sm+2-m,ipoly,al,be,&a0(m,0),&b0(m,0));
        if (ierr != 0) {
            printf("recur #4 error %d\n",ierr);
            exit(1);
        }
    }
    
    
    for(j=0;j<gpn;++j) {
        eta = n0(j);
        np(j) = eta;
        n0(j) = 2.0/(1-eta);

        ptvalues_deriv(0.0,eta);
        
        for(m=0;m<tm;++m) {
            gn(j,m) = pgn(m);
            dgn(j,m) = dpgn(m);
        }
    }

    /******************************/
    /* CALCULATE NORM                 */
    /* ************************** */
    /* SIDE & VERTEX MODES */
    for(n=0;n<3;++n)
        norm(n) = 1.0;
    for (n = 3; n < sm+3; ++n) {
        norm(n) = 0.0;
        for(i = 0; i < gpx; ++i)
            norm(n) += wtx(i)*gx(i,n)*gx(i,n);
        norm(n) = 1./sqrt(norm(n));
        norm(n+sm) = norm(n);
        norm(n+2*sm) = norm(n);
    }
    
    /* INTERIOR MODES */
    for(m = bm; m < tm; ++m) {
        norm(m) = 0.0;
        for(j = 0; j < gpn; ++j)
              norm(m) += wtn(j)*gn(j,m)*gn(j,m);
          norm(m) = 1./sqrt(norm(m));
    }

    /***************/
    /* RENORMALIZE */
    /***************/
    for (n = 3; n < nmodx; ++n) {            
        for (i =0;i<gpx;++i) {
            gx(i,n) *= norm(n);
            dgx(i,n) *= norm(n);
#ifdef DEBUG
            printf("IV1: %d %d %e %e\n",i,n,gx(i,n),dgx(i,n));
#endif
        }
        for (j=0;j<gpn;++j) {
            /* SIDE 2 & 3 MUST BE RENORMALIZED BY SAME CONSTANT TO MATCH */
            gn(j,n +sm) *= norm(n);
            dgn(j,n+sm) *= norm(n);
            gn(j,n +2*sm) *= norm(n);
            dgn(j,n +2*sm) *= norm(n);
#ifdef DEBUG
            printf("IV2: %d %d %e %e\n",j,n,gn(j,n),dgx(j,n));
#endif
        }
    }

    for(m = bm; m < tm; ++m) {
        for(j = 0; j < gpn; ++j) {
            gn(j,m) = gn(j,m)*norm(m);
            dgn(j,m) = dgn(j,m)*norm(m);
#ifdef DEBUG
            printf("IV3: %d %d %e %e\n",j,m,gn(j,m),dgn(j,m));
#endif
        }
    }

    /*****************************************/
    /* PRECALCULATE THINGS FOR FAST INTGRTRS */
    /*****************************************/
    for(m=0;m<nmodx;++m) {
        for(i=0;i<gpx;++i) {
            gxwtx(m,i) = gx(i,m)*wtx(i);
            dgxwtx(m,i) = dgx(i,m)*wtx(i);
#ifdef DEBUG
            printf("%d %d %e %e\n",m,i,gxwtx(m,i),dgxwtx(m,i));
#endif
        }
    }
    
    for(m=0;m<tm;++m) {
        for(j=0;j<gpn;++j) {
            gnwtn(m,j) = gn(j,m)*wtn(j);
            gnwtnn0(m,j) = gn(j,m)*wtn(j)*n0(j);
            dgnwtn(m,j) = dgn(j,m)*wtn(j);
#ifdef DEBUG
            printf("%d %d %e %e %e\n",m,j,gnwtn(m,j),gnwtnn0(m,j),dgnwtn(m,j));
#endif
        }
    }
    
    return;
}
    
    
void tri_basis::sideinfoinit() {
    int i,m,n,ind;
    FLT x,eta,xp1oeta,xp1,oeta;
    
    /*	THIS IS TO CALCULATE NORMAL DERIVATIVES TO SIDE */
    /* SIDES 1 & 2 ARE ROTATED TO SIDE 0 POSITION */
    dgnorm.resize(3,tm,gpx);
        
    if (p == 0) {
        for(ind=0;ind<3;++ind)
            for(i=0;i<gpx;++i)
                dgnorm(ind,0,i) = 0.0;
        
        return;
    }
    
    /*	CALCULATE NORMAL DERIVATIVE VALUES ALONG SIDES */
    /*	SIDE 0 */
    for(i=0;i<gpx;++i) {
        eta = -1.0;
        x = 2.*x0(i) -1.0;
        ptvalues_deriv(x,eta);
        xp1oeta = x0(i)*2.0/(1-eta);
        
        /* CALCULATE POLYNOMIALS */
        /* VERTEX 0    */
        dgnorm(0,0,i) = dpgn(0)*pgx(0) +xp1oeta*pgn(0)*dpgx(0);

        /* VERTEX 1  */
        dgnorm(0,1,i) = dpgn(1)*pgx(1) +xp1oeta*pgn(1)*dpgx(1);

        /* VERTEX 2     */    
        dgnorm(0,2,i) = dpgn(2)*pgx(2) +xp1oeta*pgn(2)*dpgx(2);

        for(m = 3; m < sm+3; ++m)
            dgnorm(0,m,i) = dpgn(m)*pgx(m) +xp1oeta*pgn(m)*dpgx(m);
            
        for(m=sm+3;m<2*sm+3;++m)
            dgnorm(0,m,i) = dpgn(m)*pgx(2) +xp1oeta*pgn(m)*dpgx(2);

        for(m=2*sm+3;m<bm;++m)
            dgnorm(0,m,i) = dpgn(m)*pgx(1) +xp1oeta*pgn(m)*dpgx(1);

        /*  INTERIOR MODES    */
        ind = bm;
        for(m = 3; m< sm+2;++m) {        
            for(n = 1; n < sm+3-m;++n) {
                dgnorm(0,ind,i) = dpgn(ind)*pgx(m) +xp1oeta*pgn(ind)*dpgx(m);
                ++ind;
            }
        }
    }
    
    /* SIDE 1 */
    for(i=0;i<gpx;++i) {
        eta = xp(i);
        x = 1.0;
        ptvalues_deriv(x,eta);
        oeta = 2.0/(1-eta);

        /* CALCULATE POLYNOMIALS */
        /* VERTEX 0    */
        dgnorm(1,0,i) = -oeta*pgn(0)*dpgx(0);

        /* VERTEX 1  */
        dgnorm(1,1,i) = -oeta*pgn(1)*dpgx(1);

        /* VERTEX 2     */    
        dgnorm(1,2,i) = -oeta*pgn(2)*dpgx(2);

        for(m = 3; m < sm+3; ++m)
            dgnorm(1,m,i) = -oeta*pgn(m)*dpgx(m);
            
        for(m=sm+3;m<2*sm+3;++m)
            dgnorm(1,m,i) = -oeta*pgn(m)*dpgx(2);

        for(m=2*sm+3;m<bm;++m)
            dgnorm(1,m,i) = -oeta*pgn(m)*dpgx(1);

        /*  INTERIOR MODES    */
        ind = bm;
        for(m = 3; m< sm+2;++m) {        
            for(n = 1; n < sm+3-m;++n) {
                dgnorm(1,ind,i) = -oeta*pgn(ind)*dpgx(m);
                ++ind;
            }
        }        
    }
    
    /* SIDE 2 */
    for(i=0;i<gpx;++i) {
        x = -1.0;
        eta = 1.0 -2.*x0(i);
        ptvalues_deriv(x,eta);
        oeta = 2.0/(1-eta);
        xp1 = (x+1)/2.0;

        /* CALCULATE POLYNOMIALS */
        /* VERTEX 0    */
        dgnorm(2,0,i) = -(dpgn(0)*pgx(0) +(xp1 -1)*oeta*pgn(0)*dpgx(0));

        /* VERTEX 1  */
        dgnorm(2,1,i) = -(dpgn(1)*pgx(1) +(xp1 -1)*oeta*pgn(1)*dpgx(1));

        /* VERTEX 2     */    
        dgnorm(2,2,i) = -(dpgn(2)*pgx(2) +(xp1 -1)*oeta*pgn(2)*dpgx(2));

        for(m = 3; m < sm+3; ++m)
            dgnorm(2,m,i) = -(dpgn(m)*pgx(m) +(xp1 -1)*oeta*pgn(m)*dpgx(m));
            
        for(m=sm+3;m<2*sm+3;++m)
            dgnorm(2,m,i) = -(dpgn(m)*pgx(2) +(xp1 -1)*oeta*pgn(m)*dpgx(2));

        for(m=2*sm+3;m<bm;++m)
            dgnorm(2,m,i) = -(dpgn(m)*pgx(1) +(xp1 -1)*oeta*pgn(m)*dpgx(1));

        /*  INTERIOR MODES    */
        ind = bm;
        for(m = 3; m< sm+2;++m) {        
            for(n = 1; n < sm+3-m;++n) {
                dgnorm(2,ind,i) = -(dpgn(ind)*pgx(m) +(xp1 -1)*oeta*pgn(ind)*dpgx(m));
                ++ind;
            }
        }
    }

    return;
}

// #define NEWWAY
#ifdef NEWWAY
/************************************************/
/** CALCULATE THINGS FOR LUMPED MASS INVERSION  */
/************************************************/
void tri_basis::lumpinv(void) {
    int i,j,k,m,info;
    Array<int,1> ipiv(2*tm);
    Array<FLT,2> mm(tm,tm);
	Array<FLT,2> mmpsi(bm,tm);
	Array<FLT,3> psi(bm,gpx,gpn);
    Array<FLT,1> uht(tm),l(tm);
	FLT psinorm;
	Array<FLT,1> apsi(p), bpsi(p);
    char trans[] = "T", uplo[] = "U";
	int ipoly,ierr;
	FLT al,be;
	FLT pkp,pk,pkm,pks;
	FLT r,s,u,v,x,n;
    
    /* ALLOCATE MASS MATRIX INVERSION VARIABLES */
    if (sm > 0) {
        vfms.resize(3,sm);
        sfmv.resize(2,sm);
        sdiag.resize(sm);
    }
    if (sm > 1) sfms.resize(sm-1,sm,3);
    if (im > 0) {
        ifmb.resize(bm,im);
        bfmi.resize(bm,im);
        idiag.resize(im,ibwth+1);
    }
    msi.resize(bm,bm);
	
    /********************************************/    
    /* GENERATE MASS MATRIX                        */
    /********************************************/
    for(m=0;m<tm;++m) {
        for(i=0;i<tm;++i)
            uht(i) = 0.0;
        uht(m) = 1.0;

        proj(&uht(0),&psi(0,0,0),gpn);  // PROJECT USES WK0
        intgrt(&l(0),&psi(0,0,0),gpn); // INTGRT USES WK0

#ifdef DEBUG
        for(i=0;i<gpx;++i)
            for(j=0;j<gpn;++j)
                printf("MIMM: %d %d %e\n",i,j,psi(0,i,j));
        for(i=0;i<tm;++i)
            printf("MIMM2: %d %e\n",i,l(i));
#endif
                
        for(i=0;i<tm;++i) 
            mm(m,i) = l(i);        
    }
	
	
	/*******************************************/
	/* GENERATE PSI BASIS VERTEX & SIDE MODES  */
	/*******************************************/
	/* VERTEX */
    ipoly = 6;
    al = 2.0;
    be = 1.0;
    ierr = recur(sm+1,ipoly,al,be,&apsi(0),&bpsi(0));
    if (ierr != 0) {
        printf("recur #3 error %d\n",ierr);
        exit(1);
    }
	
	Array<FLT,2> m1d(sm,sm);
	Array<FLT,2> psin(sm,sm);
	psin = 0.0;
	
	for (m=0;m<sm;++m) {
		m1d = 0.0;
		/* FIND FUNCTION ORTHOGONAL TO OPPOSING SIDE MODES */
		for(i=0;i<sm-1-m;++i) {
			for(j=0;j<sm-m;++j) {
				for(k=0;k<gpn;++k) {
					m1d(i+1,j) += gnwtn(3+sm+i+m,k)*pow(0.5*(1.-np(k)),2.+m)*pow(0.5*(1.+np(k)),j);
				}
			}
		}
		m1d(0,0) = 1.0;  // NORMALIZATION FOR FIRST COEFFICIENT 
		psin(m,0) = 1.0;  // NORMALIZATION CONSTRAINT 

		/* COEFFICIENTS TO FORM PSI BASIS FROM PHI BASIS */
		GETRF(sm-m,sm-m,m1d.data(),sm,&ipiv(0),info);
		if (info != 0) {
			printf("DGETRF FAILED info:%d sm:%d i:%d\n",info,(sm+2),i);
			exit(1);
		}
		GETRS(trans,sm-m,1,m1d.data(),sm,&ipiv(0),&psin(m,0),sm,info);
	}

	/* CALCULATE NORMALIZATION CONSTANT FOR VERTEX FUNCTION */
	pk = 1.0;
	pkm = 0.0;
	for(m=0;m<p-1;++m) {
		pkp = (1.0-apsi(m))*pk - bpsi(m)*pkm;
		pkm = pk;
		pk = pkp;
	}
	psinorm = pk;
	
	/* NOW CALCULATE PSI BASIS FUNCTIONS */
	for (i=0;i<gpx;++i) {
		for(j=0;j<gpn;++j) {
		
			x = xp(i);
			n = np(j);
		
			/* VERTEX 0 */
			pk = 1.0;
			pkm = 0.0;
			for(m=0;m<p-1;++m) {
				pkp = (n-apsi(m))*pk - bpsi(m)*pkm;
				pkm = pk;
				pk = pkp;
			}
			psi(0,i,j) = pk/psinorm*(1.+n)*.5;
			
			/* SIDE 0 */
			pk = 1.0;
			pkm = 0.0;
			for(m=0;m<sm;++m) {
				/* CALCULATE ETA FUNCTION */
				pks = 0.0;
				for(k=0;k<sm-m;++k) 
					pks += psin(m,k)*pow(0.5*(1.-n),2.+m)*pow(0.5*(1.+n),k);

				/* CALCULATE PSI FUNCTION */
				psi(m+3,i,j) = pks*(1.-x)*(1.+x)*0.25*pk*norm(m+3);
				pkp = (x-a0(0,m))*pk - b0(0,m)*pkm;
				pkm = pk;
				pk = pkp;  
			}
			
			/* NOW ROTATED FUNCTIONS */
			r = (x+1.)*(1.-n)*0.5 -1.0;
			s = n;

			u = s;
			v = -1.-r-s;
			
			x = 2.0*(1.+u)/(1.-v) -1.0;
			n = v;
			
			/* VERTEX 1 */
			pk = 1.0;
			pkm = 0.0;
			for(m=0;m<p-1;++m) {
				pkp = (n-apsi(m))*pk - bpsi(m)*pkm;
				pkm = pk;
				pk = pkp;
			}
			psi(1,i,j) = pk/psinorm*(1.+n)*.5;
			
			/* SIDE 1 */
			pk = 1.0;
			pkm = 0.0;
			for(m=0;m<sm;++m) {
				/* CALCULATE ETA FUNCTION */
				pks = 0.0;
				for(k=0;k<sm-m;++k) 
					pks += psin(m,k)*pow(0.5*(1.-n),2.+m)*pow(0.5*(1.+n),k);

				/* CALCULATE PSI FUNCTION */
				psi(m+3+sm,i,j) = pks*(1.-x)*(1.+x)*0.25*pk*norm(m+3);
				pkp = (x-a0(0,m))*pk - b0(0,m)*pkm;
				pkm = pk;
				pk = pkp;  
			}
			
			/* NOW ROTATED FUNCTIONS AGAIN */
			r = (x+1.)*(1.-n)*0.5 -1.0;
			s = n;

			u = s;
			v = -1.-r-s;
			
			x = 2.0*(1+u)/(1-v) -1.0;
			n = v;
			
			/* VERTEX 2 */
			pk = 1.0;
			pkm = 0.0;
			for(m=0;m<p-1;++m) {
				pkp = (n-apsi(m))*pk - bpsi(m)*pkm;
				pkm = pk;
				pk = pkp;
			}
			psi(2,i,j) = pk/psinorm*(1.+n)*.5;
			
			/* SIDE 2 */
			pk = 1.0;
			pkm = 0.0;
			for(m=0;m<sm;++m) {
				/* CALCULATE ETA FUNCTION */
				pks = 0.0;
				for(k=0;k<sm-m;++k) 
					pks += psin(m,k)*pow(0.5*(1.-n),2.+m)*pow(0.5*(1.+n),k);

				/* CALCULATE PSI FUNCTION */
				psi(m+3+2*sm,i,j) = pks*(1.-x)*(1.+x)*0.25*pk*norm(m+3);
				pkp = (x-a0(0,m))*pk - b0(0,m)*pkm;
				pkm = pk;
				pk = pkp;  
			}
		}
	}
	
	for (m=0;m<bm;++m) {
        intgrt(&l(0),&psi(m,0,0),gpn);
		mmpsi(m,Range::all()) = l;
	}
			
	if (p > 1)
		vdiag = mmpsi(0,0);
	else 
		vdiag = mm(0,0) + mm(0,1) + mm(0,2);
		
	/* SIDES LOWER TRIANGLE*/
	for(k=0;k<sm;++k) {
        for(j=0;j<3;++j)
            vfms(j,k) = mmpsi(3+k,j);
            
        sdiag(k) = 1./mmpsi(3+k,3+k);
        
        for(j=0;j<k;++j) {
            sfms(j,k,0) = mmpsi(3+k,3+j);
            sfms(j,k,1) = mmpsi(3+k,3+j+sm);
            sfms(j,k,2) = mmpsi(3+k,3+j+2*sm);
        }
    }
	
	/* REMOVAL FROM INTERIORS */
	for(i=0; i<bm; ++i ) {
		for(j=0; j<im; ++j ) {
			bfmi(i,j) = mm(j+bm,i);
		}
	}
	
	/* INTERIOR - INTERIOR MATRIX (THIS IS STUPID: Identity) */
	int i1;
	if (im > 0) {
		/* SETUP DIAGANOL FORM OF MATRIX */
		/* ONLY NECESSARY WHEN USING DPBTRF */    
		for(j=0;j<im;++j) {
			i1 = (0 > j-ibwth ? 0 : j-ibwth);
			for(i=i1;i<=j;++i) {
				k = i - (j-ibwth);
				idiag(j,k) = mm(i+bm,j+bm);
			}
		}
		
		PBTRF(uplo,im,ibwth,&idiag(0,0),ibwth+1,info);
		if (info != 0) {
			printf("1:PBTRF FAILED info: %d\n", info);
			exit(1);
		}
	}
	
	/* COEFFICIENTS TO FORM PSI BASIS FROM PHI BASIS */
	GETRF(tm,tm,&mm(0,0),tm,&ipiv(0),info);
	if (info != 0) {
		printf("DGETRF FAILED info:%d sm:%d i:%d\n",info,(sm+2),i);
		exit(1);
	}
	GETRS(trans,tm,bm,mm.data(),tm,&ipiv(0),mmpsi.data(),tm,info);
	
	/* STORE SIDE VALUES */
	for(k=0;k<sm;++k) {
		sfmv(0,k) = -mmpsi(1,k+3);
		sfmv(1,k) = -mmpsi(2,k+3);
	}
	
	for(k=0;k<bm;++k)
		for(j=0;j<im;++j)
			ifmb(k,j) = -mmpsi(k,j+bm);
			
			
	/************************************************/
	/* FIND MATRICES TO DETERMINE INTERIOR MODES */
	/************************************************/
	if (im > 0) {
		/* PREMULTIPLY MATRIX TO REMOVE BOUNDARY MODES FROM INTERIOR */
		for(i=0; i<bm; ++i ) {
			PBTRS(uplo,im,ibwth,1,&idiag(0,0),ibwth+1,&bfmi(i,0),im,info);
		}
	}
        

	
#ifdef DEBUG
	std::cout << "sfmv" << std::endl;
	std::cout << sfmv << std::endl;
	std::cout << "ifmb" << std::endl;
	std::cout << ifmb << std::endl;
	std::cout << "vdiag" << std::endl;
	std::cout << vdiag << std::endl;
	std::cout << "sdiag" << std::endl;
	std::cout << sdiag << std::endl;
	std::cout << "vfms" << std::endl;
	std::cout << vfms << std::endl;
	std::cout << "sfms" << std::endl;
	std::cout << sfms << std::endl;
	std::cout << "bfmi" << std::endl;
	std::cout << bfmi << std::endl;
	std::cout << "idiag" << std::endl;
	std::cout << idiag << std::endl;

	int i2;
	Array<FLT,2> mwk(MXTM,MXTM);
    Array<FLT,1> vwk(MXTM),wk1(MXTM*MXGP);

    /* CHECK TO MAKE SURE PREVIOUS RESULTS ARE RIGHT */
    for(i=0;i<3;++i) {
        i1 = (i+1)%3;
        i2 = (i+2)%3;

        uht(i) = 1.0;
        uht(i1) = 0.0;
        uht(i2) = 0.0;
        
        for(j=0;j<sm;++j) {
            uht(3+j+i*sm) = 0.0;
            uht(3+j+i1*sm) = -sfmv(1,j);
            uht(3+j+i2*sm) = -sfmv(0,j);
        }
        
        for(j=0;j<im;++j)
            uht(bm+j) = -ifmb(i,j);
            
        proj(&uht(0),&wk1(0),gpn);
        intgrt(&mwk(i,0),&wk1(0),gpn);
        printf("%2d:",i);
        for(j=0;j<tm;++j)
            printf("%+.4le  ",mwk(i,j));
        printf("\n");
    }
    
    for(i=0;i<3;++i) {
        for(k=0;k<sm;++k) {
            for(j=0;j<tm;++j)
                uht(j) = 0.0;
            uht(i*sm+k+3) = 1.0;
            for(j=0;j<im;++j) {
                uht(j+bm) = -ifmb(i*sm+k+3,j);
            }

            proj(&uht(0),&wk1(0),gpn);
            intgrt(&mwk(i*sm+k+3,0),&wk1(0),gpn);
            
            printf("%2d:",3+i*sm+k);
            for(j=0;j<tm;++j)
                printf("%+.4le  ",mwk(i*sm+k+3,j));
            printf("\n");
        }            
    }
    
    for(i=0;i<im;++i) {
        for(j=0;j<tm;++j)
            uht(j) = 0.0;
        uht(i +bm) = 1.0;  
        proj(&uht(0),&wk1(0),gpn);
        intgrt(&mwk(i+bm,0),&wk1(0),gpn);
    }
    
    for(i=0;i<3;++i)
        for(m=0;m<tm;++m)
            mwk(i,m) /= vdiag;
            
    for(i=0;i<3;++i) {
        for(j=0;j<3;++j) {
            i1 = (i+j)%3;
            for(k=0;k<sm;++k) {
                for(m=0;m<tm;++m)
                    mwk(i*sm+k+3,m) -= vfms(j,k)*mwk(i1,m);
            }
        }
    }
    
    /* REMOVE MODES J,K FROM MODE I,M */
    int mode;
    for(mode = 0; mode <sm;++mode) {
        for(i=0;i<3;++i)
            for(k=0;k<tm;++k)
                mwk(3+i*sm+mode,k) *= sdiag(mode);
        for(i=0;i<3;++i) {
            for(m=mode+1;m<sm;++m) {
                for(j=0;j<3;++j) {
                    i1 = (i+j)%3;
                    for(k=0;k<tm;++k)
                        mwk(3+i*sm+m,k) -= sfms(mode,m,j)*mwk(3+i1*sm+mode,k);
                }
            }
        }
    }
    DPBTRSNU1(&idiag(0,0),im,ibwth,&mwk(bm,0),MXTM);
    for(k=0;k<im;++k)
        for(i=0;i<bm;++i)
            for(j=0;j<tm;++j)
                mwk(bm+k,j) -= bfmi(i,k)*mwk(i,j);
                            
    for(m=0;m<tm;++m) {
        printf("LI1: %2d:",m);
        for(j=0;j<tm;++j) {
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
#else
/************************************************/
/** CALCULATE THINGS FOR LUMPED MASS INVERSION  */
/************************************************/
void tri_basis::lumpinv(void) {
    int i,i1,i2,j,k,m,info,ind,ind1,n;
    Array<int,1> ipiv(2*MXTM);
    Array<FLT,2> mwk(MXTM,MXTM),mm(MXTM,MXTM);
    Array<FLT,1> u(MXTM),l(MXTM),vwk(MXTM),wk1(MXTM*MXGP);
    FLT rcond=1;
    char trans[] = "T", uplo[] = "U";
    
    /* ALLOCATE MASS MATRIX INVERSION VARIABLES */
    if (sm > 0) {
        vfms.resize(3,sm);
        sfmv.resize(2,sm);
        sdiag.resize(sm);
    }
    if (sm > 1) sfms.resize(sm-1,sm,3);
    if (im > 0) {
        ifmb.resize(bm,im);
        bfmi.resize(bm,im);
        idiag.resize(im,ibwth+1);
    }
    msi.resize(bm,bm);

    /********************************************/    
    /* GENERATE MASS MATRIX                        */
    /********************************************/
    for(m=0;m<tm;++m) {
        for(i=0;i<tm;++i)
            u(i) = 0.0;
        u(m) = 1.0;

        proj(&u(0),&wk1(0),gpn);  // PROJECT USES WK0
        intgrt(&l(0),&wk1(0),gpn); // INTGRT USES WK0

#ifdef DEBUG
        for(i=0;i<gpx;++i)
            for(j=0;j<gpn;++j)
                printf("MIMM: %d %d %e\n",i,j,wk1(i*gpn +j));
        for(i=0;i<tm;++i)
            printf("MIMM2: %d %e\n",i,l(i));
#endif
                
        for(i=0;i<tm;++i) 
            mm(m,i) = l(i);        
    }
	
    /*******************************************************/        
    /*  EQUATIONS TO FIND VERTEX VALUES TO SM-1 ACCURACY */
    /*******************************************************/        
    if (p > 1) {        
        for(i=0;i<3;++i) {
            i1 = (i+1)%3;
            i2 = (i+2)%3;
            
            /*2 CROSS VERTEX CONSTRAINTS */
            for(k=0;k<2;++k) {
                ind = (i+1+k)%3;
                vwk(k) = mm(ind,i);
                for(j=0;j<sm;++j) {
                    mwk(k,j) = mm(ind,3+j+i1*sm);
                    mwk(k,j+sm) = mm(ind,3+j+i2*sm);
                }
                for(j=0;j<im;++j)
                    mwk(k,2*sm+j) = mm(ind,j+bm);
            }

            /*3x(SM-3) SIDE CONSTRAINTS */
            for(k=0;k<3;++k) {
                for(m=0;m<sm-1;++m) {    
                    vwk(m+2+k*(sm-1)) = mm(3+m+k*sm,i);
                    for(j=0;j<sm;++j) {
                        mwk(m+2+k*(sm-1),j) = mm(3+m+k*sm,3+j+i1*sm);
                        mwk(m+2+k*(sm-1),j+sm) = mm(3+m+k*sm,3+j+i2*sm);
                    }
                    for(j=0;j<im;++j)
                        mwk(m+2+k*(sm-1),2*sm+j) = mm(3+m+k*sm,j+bm);
                }
            }
            
            
            /*(SM-3)*(SM-4) INTERNAL MODE CONSTRAINTS  */                
            ind = 3*sm-1;    ind1 = 0;
            for(m=2;m<sm+1;++m) {        
                for(n = 1; n < sm+1-m;++n) {
                    vwk(ind) = mm(ind1+bm,i);
                    for(j=0;j<sm;++j) {
                        mwk(ind,j) = mm(ind1+bm,3+j+i1*sm);
                        mwk(ind,j+sm) = mm(ind1+bm,3+j+i2*sm);
                    }
                    for(j=0;j<im;++j)
                        mwk(ind,2*sm+j) = mm(ind1+bm,j+bm);
                    ++ind;
                    ++ind1;
                }
                ++ind1;
            }

            GETRF((sm+1)*(sm+2)/2-1,(sm+1)*(sm+2)/2-1,&mwk(0,0),MXTM,&ipiv(0),info);
            if (info != 0) {
                printf("DGETRF FAILED - VRTX info:%d sm:%d i:%d\n",info,(sm+2),i);
                exit(1);
            }
            GETRS(trans,(sm+1)*(sm+2)/2-1,1,&mwk(0,0),MXTM,&ipiv(0),&vwk(0),MXTM,info);
                                                    
            /* STORE INTERIOR VALUES */
            for(k=0;k<im;++k)
                ifmb(i,k) = vwk(k+2*sm);            
        }
    
        /* STORE SIDE VALUES */
        for(k=0;k<sm;++k) {
            sfmv(0,k) = vwk(k+sm);
            sfmv(1,k) = vwk(k);
            
        }

        /* FIND VERTEX DIAGANOL ELEMENT */
        vdiag = mm(2,2);
        for(k=0;k<sm;++k) {
            vdiag -= sfmv(0,k)*mm(2,3+k+1*sm);
            vdiag -= sfmv(1,k)*mm(2,3+k);
        }

        for(k=0;k<im;++k)
            vdiag -= ifmb(2,k)*mm(2,k+bm);
    }
    else {
        if (p > 0)
            vdiag = mm(0,0) + mm(0,1) + mm(0,2);
        else
            vdiag = mm(0,0);
    }


    /*****************************************************************/
    /* EQUATIONS TO FIND STATIC INVERSION OF INTERIORS FROM SIDES */
    /*****************************************************************/
    /** WARNING THIS ONLY WORKS UP TO SM = 5!!!! ***/
    if (im > 0) {
        for(i=0;i<3;++i) {
            i1 = (i+1)%3;
            i2 = (i+2)%3;

            for(k=0;k<(sm+2)-3;++k) {
                ind = 0;
                /* CROSS SIDE MODES */
                for(j=k;j<sm-1;++j) {
                    vwk(ind) =  mm(i*sm+3+k,i1*sm+3+j);
                    for(m=0;m<im;++m)  
                        mwk(ind,m) = mm(bm+m,i1*sm+3+j);
                    ++ind;
                }
                
                for(j=k;j<MIN(2*k,sm-1);++j) {
                    vwk(ind) =  mm(i*sm+3+k,i2*sm+3+j);
                    for(m=0;m<im;++m)  
                        mwk(ind,m) = mm(bm+m,i2*sm+3+j);
                    ++ind;
                }                            

                /* INTERIOR MODES */    
                ind1 = 0;
                for(m = 2; m< sm+1;++m) {        
                    for(n = 1; n < sm+1-m;++n) {
                        vwk(ind) = mm(i*sm+3+k,bm+ind1);
                        for(j=0;j<im;++j) 
                            mwk(ind,j) = mm(bm+j,bm+ind1);
                        ++ind;
                        ++ind1;
                    }
                    ++ind1;
                }

                GETRF(im,im,&mwk(0,0),MXTM,&ipiv(0),info);
                if (info != 0) {
                    printf("DGETRF FAILED SIDE info:%d (sm+2):%d k:%d\n",info,(sm+2),k);
                    exit(1);
                }
                else
                    GETRS(trans,im,1,&mwk(0,0),MXTM,&ipiv(0),&vwk(0),MXTM,info);

                for(j=0;j<im;++j)
                    ifmb(i*sm+3+k,j) = vwk(j);
            }

            /* FOR HIGHEST ORDER MODE - STATIC INVERT ALL INTERIOR MODES */
            for(j=0;j<im;++j) {
                vwk(j) = mm(i*sm +sm +2,j+bm);
                for(k=0;k<im;++k)
                    mwk(j,k) = mm(k+bm,j+bm);
            }

            GETRF(im,im,&mwk(0,0),MXTM,&ipiv(0),info);
            if (info != 0) {
                printf("GETRF FAILED info: %d cond: %f\n",info,rcond);
                exit(1);
            }
            GETRS(trans,im,1,&mwk(0,0),MXTM,&ipiv(0),&vwk(0),MXTM,info);

            for(j=0;j<im;++j)
                ifmb(i*sm +sm +2,j) = vwk(j);
        }
    }

    /*********************************************/
    /* FIND VERTEX-SIDE & SIDE-SIDE MATRICES  */
    /*********************************************/
    for(k=0;k<sm;++k) {
        for(j=0;j<tm;++j)
            u(j) = 0.0;
            
        u(k+3) = 1.0;
        for(j=0;j<im;++j) {
            u(j+bm) = -ifmb(k+3,j);
        }

        proj(&u(0),&wk1(0),gpn);
        intgrt(&l(0),&wk1(0),gpn);
        
        for(j=0;j<3;++j)
            vfms(j,k) = l(j);
            
        sdiag(k) = 1.0/l(3+k);
        
        for(j=0;j<k;++j) {
            sfms(j,k,0) = l(3+j);
            sfms(j,k,1) = l(3+j+sm);
            sfms(j,k,2) = l(3+j+2*sm);
        }
    }


    /************************************************/
    /* FIND MATRICES TO DETERMINE INTERIOR MODES */
    /************************************************/
    if (im > 0) {
        /* SETUP DIAGANOL FORM OF MATRIX */
        /* ONLY NECESSARY WHEN USING DPBTRF */    
        for(j=0;j<im;++j) {
            i1 = (0 > j-ibwth ? 0 : j-ibwth);
            for(i=i1;i<=j;++i) {
                k = i - (j-ibwth);
                idiag(j,k) = mm(i+bm,j+bm);
            }
        }
        
        PBTRF(uplo,im,ibwth,&idiag(0,0),ibwth+1,info);
        if (info != 0) {
            printf("1:PBTRF FAILED info: %d\n", info);
            exit(1);
        }
                    
        /* MATRIX TO REMOVE BOUNDARY MODES FROM INTERIOR */
        for(i=0; i<bm; ++i ) {
            for(j=0; j<im; ++j ) {
                bfmi(i,j) = mm(i,j+bm);
            }
            PBTRS(uplo,im,ibwth,1,&idiag(0,0),ibwth+1,&bfmi(i,0),im,info);
        }
        
        /* TEST INVERSE 
        for(i=0;i<im;++i) {
            for(j=0;j<bm;++j)
                printf(" %f ",bfmi(j,i));     
            printf("\n");
        }
            
        printf("ibwth %d\n",ibwth);
        for(i=0;i<im;++i) {
            PBTRS(uplo,im,ibwth,1,&idiag(0,0),ibwth+1,&mm(i+bm,bm),im,info);
            printf("%d: ",i);
            for(j=0;j<im;++j)
                printf("%f ",mm(i+bm,bm+j));
            printf("\n");
        }
        */
    }
        
    /* FORM A - BC^{-1}B^T */  
    /* MASS MATRIX WITH STATIC CONDENSATION OF INTERIOR MODES */
    for(i=0;i<bm;++i) {
        for(j=0;j<bm;++j) {
            mwk(i,j) = mm(i,j);
            for(k=0;k<im;++k)
                mwk(i,j) -= bfmi(i,k)*mm(j,k+bm);
        }
    }

#ifndef SKIP
    /* SETUP BAND DIAGONAL FORM OF MATRIX */
    /* ONLY NECESSARY WHEN USING DPBTRF */ 
    for(j=0;j<bm;++j) {
        for(i=0;i<=j;++i) {
            k = i - (j-(bm-1));
            msi(j,k) = mwk(i,j);
        }
    }
    PBTRF(uplo,bm,bm-1,&msi(0,0),bm,info);
    if (info != 0) {
        printf("2:PBTRF FAILED info: %d\n", info);
        exit(1);
    }
#else
    for(i=0;i<bm;++i)
        for(j=0;j<bm;++j) 
            msi(i,j) = mwk(i,j);
            
    DPOTRF(uplo,bm,&msi(0,0),bm,info);
    if (info != 0) {
        printf("POTRF FAILED info: %d\n", info);
        exit(1);
    }
#endif


#ifdef DEBUG
	std::cout << "sfmv" << std::endl;
	std::cout << sfmv << std::endl;
	std::cout << "ifmb" << std::endl;
	std::cout << ifmb << std::endl;
	std::cout << "vdiag" << std::endl;
	std::cout << vdiag << std::endl;
	std::cout << "sdiag" << std::endl;
	std::cout << sdiag << std::endl;
	std::cout << "vfms" << std::endl;
	std::cout << vfms << std::endl;
	std::cout << "sfms" << std::endl;
	std::cout << sfms << std::endl;
	std::cout << "bfmi" << std::endl;
	std::cout << bfmi << std::endl;
	std::cout << "idiag" << std::endl;
	std::cout << idiag << std::endl;

    /* CHECK TO MAKE SURE PREVIOUS RESULTS ARE RIGHT */
    for(i=0;i<3;++i) {
        i1 = (i+1)%3;
        i2 = (i+2)%3;

        u(i) = 1.0;
        u(i1) = 0.0;
        u(i2) = 0.0;
        
        for(j=0;j<sm;++j) {
            u(3+j+i*sm) = 0.0;
            u(3+j+i1*sm) = -sfmv(1,j);
            u(3+j+i2*sm) = -sfmv(0,j);
        }
        
        for(j=0;j<im;++j)
            u(bm+j) = -ifmb(i,j);
            
        proj(&u(0),&wk1(0),gpn);
        intgrt(&mwk(i,0),&wk1(0),gpn);
        printf("%2d:",i);
        for(j=0;j<tm;++j)
            printf("%+.4le  ",mwk(i,j));
        printf("\n");
    }
    
    for(i=0;i<3;++i) {
        for(k=0;k<sm;++k) {
            for(j=0;j<tm;++j)
                u(j) = 0.0;
            u(i*sm+k+3) = 1.0;
            for(j=0;j<im;++j) {
                u(j+bm) = -ifmb(i*sm+k+3,j);
            }

            proj(&u(0),&wk1(0),gpn);
            intgrt(&mwk(i*sm+k+3,0),&wk1(0),gpn);
            
            printf("%2d:",3+i*sm+k);
            for(j=0;j<tm;++j)
                printf("%+.4le  ",mwk(i*sm+k+3,j));
            printf("\n");
        }            
    }
    
    for(i=0;i<im;++i) {
        for(j=0;j<tm;++j)
            u(j) = 0.0;
        u(i +bm) = 1.0;  
        proj(&u(0),&wk1(0),gpn);
        intgrt(&mwk(i+bm,0),&wk1(0),gpn);
    }
    
    for(i=0;i<3;++i)
        for(m=0;m<tm;++m)
            mwk(i,m) /= vdiag;
            
    for(i=0;i<3;++i) {
        for(j=0;j<3;++j) {
            i1 = (i+j)%3;
            for(k=0;k<sm;++k) {
                for(m=0;m<tm;++m)
                    mwk(i*sm+k+3,m) -= vfms(j,k)*mwk(i1,m);
            }
        }
    }
    
    /* REMOVE MODES J,K FROM MODE I,M */
    int mode;
    for(mode = 0; mode <sm;++mode) {
        for(i=0;i<3;++i)
            for(k=0;k<tm;++k)
                mwk(3+i*sm+mode,k) *= sdiag(mode);
        for(i=0;i<3;++i) {
            for(m=mode+1;m<sm;++m) {
                for(j=0;j<3;++j) {
                    i1 = (i+j)%3;
                    for(k=0;k<tm;++k)
                        mwk(3+i*sm+m,k) -= sfms(mode,m,j)*mwk(3+i1*sm+mode,k);
                }
            }
        }
    }
    DPBTRSNU1(&idiag(0,0),im,ibwth,&mwk(bm,0),MXTM);
    for(k=0;k<im;++k)
        for(i=0;i<bm;++i)
            for(j=0;j<tm;++j)
                mwk(bm+k,j) -= bfmi(i,k)*mwk(i,j);
                            
    for(m=0;m<tm;++m) {
        printf("LI1: %2d:",m);
        for(j=0;j<tm;++j) {
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
#endif


void tri_basis::lumpinv1d() {
	int i,i1,j,k,m,info;
    Array<int,1> ipiv(2*MXTM);
    Array<FLT,2> mwk(MXTM,MXTM),mm(MXTM,MXTM);
    Array<FLT,1> u(MXTM),l(MXTM),vwk(MXTM),wk1(MXTM*MXGP);
    FLT rcond=1;
    char trans[] = "T", uplo[] = "U";

    /********************************************************************/    
    /* NOW SETUP SIMILAR THING FOR 1-D MATRICES                                     */
    /********************************************************************/    
	/* ALLOCATE 1D MASS MATRIX INVERSION VARIABLES */
    if (sm > 0) {
        vfms1d.resize(2,sm);
        sfmv1d.resize(2,sm);
        sdiag1d.resize(sm,sbwth+1);
    }
	
	
    /* CALCULATE 1-D MASS MATRIX */    
    for(m=1;m<p+2;++m) {
        for(k=1;k<p+2;++k) {
            mm(m-1,k-1)= 0.0;
            for(i=0;i<gpx;++i)
                mm(m-1,k-1) += wtx(i)*gx(i,m)*gx(i,k);
        }
    }

#ifdef DEBUG
    printf("\n mm MATRIX (sm+2) = %d\n",(sm+2));
    for(i=0;i<sm+2;++i) {
        printf("LI2: %2d:",i);
        for(j=0;j<sm+2;++j)
            printf("%+.4lf ",mm(i,j));
        printf("\n");
    }
#endif

    if(p > 1) {
        /* STATIC INVERSION SIDE MATRIX  */            
        for(j=0;j<sm;++j) {
            i1 = (0 > j-sbwth ? 0 : j-sbwth);
            for(i=i1;i<=j;++i) {
                k = i - (j-sbwth);
                sdiag1d(j,k) = mm(i+2,j+2);
            }
        }

        PBTRF(uplo,sm,sbwth,&sdiag1d(0,0),sbwth+1,info);
        if (info != 0 || rcond < 100.*EPSILON) {
            printf("PBTRF FAILED - 1D (sm+2) : %d info: %d cond: %f\n",(sm+2), info,rcond);
            exit(1);
        }

        /* MATRIX TO REMOVE VERTEX MODES FROM SIDES */
        for(i=0; i<2; ++i) {
            for(j=0; j<sm; ++j) {
                vfms1d(i,j) = mm(i,j+2);
            }
            PBTRS(uplo,sm,sbwth,1,&sdiag1d(0,0),sbwth+1,&vfms1d(i,0),(sm+2)-2,info);
        }

        /* MATRIX TO REMOVE SIDE MODES FROM VERTICES */    
        for(i=0;i<2;++i) {        
            /* VERTEX CONSTRAINT */
            vwk(0) = mm(i,(i+1)%2);
            for(m=0;m<sm;++m)
                mwk(0,m) = mm((i+1)%2,m+2);

            /* SIDE CONSTRAINTS */
            for(m=0;m<sm-1;++m) {
                vwk(m+1) = mm(i,m+2);
                for(k=0;k<sm;++k) {
                    mwk(m+1,k) = mm(m+2,k+2);
                }
            }
            GETRF(sm,sm,&mwk(0,0),MXTM,&ipiv(0),info);
            if (info != 0) {
                printf("1D DGETRF FAILED - VRTX info:%d (sm+2):%d i:%d\n",info,(sm+2),i);
                exit(1);
            }
            GETRS(trans,sm,1,&mwk(0,0),MXTM,&ipiv(0),&vwk(0),MXTM,info);

            /* STORE MATRIX */
            for(k=0;k<sm;++k)
                sfmv1d(i,k) = vwk(k);
        }

        /* FIND VERTEX DIAGANOL ELEMENT */
        vdiag1d = mm(0,0);
        for(k=0;k<sm;++k)
            vdiag1d -= sfmv1d(0,k)*mm(0,2+k);
    }
    else {
        if (p > 0) 
            vdiag1d = mm(0,0)+mm(0,1);
        else
            vdiag1d = mm(0,0);
    }

#ifdef DEBUG
    /* CHECK TO MAKE SURE PREVIOUS 1D RESULTS ARE RIGHT */
    Array<FLT,1> uht(MXTM);
    for(k=0;k<2;++k) {
        u(k) = 1.0;
        u((k+1)%2) = 0.0;
        for(m=0;m<sm;++m)
            u(m+2) = -sfmv1d(k,m);

        for(i=0;i<gpx;++i) {
            uht(i) = u(0)*gx(i,0);
            uht(i) += u(1)*gx(i,1);
            for(m=0;m<sm;++m)
                uht(i) += u(m+2)*gx(i,3+m);
        }

        for(m=1;m<sm+3;++m) {
            l(m) = 0.0;
            for(i=0;i<gpx;++i) {
                l(m) += wtx(i)*gx(i,m)*uht(i);
            }
        }
        printf("LI3 %2d:",k);
        for(j=1;j<sm+3;++j)
            printf("%+.4le  ",l(j));
        printf("%+.4le  ",vdiag1d);
        printf("\n");
    }
#endif

    return;
}            

void tri_basis::legpt()
{
    FLT x,eta,r,s;
    int i,j,m,n,ind;
    
    lgrnge1d.resize(nmodx,sm+1);
    lgrnge.resize(tm,sm+1,sm+1);
    
//    /* USING LOBATTO SPACING */
//    int ipoly = 1;
//    double al = 0.0;
//    double be = 0.0;
//    ierr = recur(sm+3,ipoly,al,be,&a0(0),&b0(0));
//    if (ierr != 0) {
//        printf("recur #1 error %d\n",ierr);
//        exit(1);
//    }
//    ierr = lob(sm,&a0(0),&b0(0),-1.0,1.0,&pts(0),&e3(0),&e(0),&e1(0),&e2(0));        
//    if (ierr != 0) {
//        printf("gauss #3 error %d\n",ierr);
//        exit(1);
//    }

    /* CALCULATE PROJECTION POINTS IN INTERIOR */
    for(i=1;i<sm;++i) {
        for(j=1;j<sm-(i-1);++j) {
            s = -1 +2.0*((FLT) j)/(FLT)(sm+1);
            r = -1 +2.0*((FLT) i)/(FLT)(sm+1);
            x = 2.0*(1+r)/(1-s) -1.0;
            eta = s;
                      
            /* CALCULATE VALUES OF PSI POLYNOMIALS AT POINT */
            ptvalues(x,eta);
      
            /* CALCULATE S POLYNOMIALS */
            /* VERTEX A    */
            lgrnge(0,i,j) = pgn(0)*pgx(0);

            /* VERTEX B  */
            lgrnge(1,i,j) = pgn(1)*pgx(1);

            /* VERTEX C     */    
            lgrnge(2,i,j) = pgn(2)*pgx(2);

            /*  SIDE 1 (s)        */
            for(m = 3; m < sm+3; ++m)
                lgrnge(m,i,j) = pgn(m)*pgx(m);
                
            for(m=sm+3;m<2*sm+3;++m)
                lgrnge(m,i,j) = pgn(m)*pgx(2);

            for(m=2*sm+3;m<bm;++m)
                lgrnge(m,i,j) = pgn(m)*pgx(1);

            /*  INTERIOR MODES    */
            ind = bm;
            for(m = 3; m< sm+2;++m) {        
                for(n = 1; n < sm+3-m;++n) {
                    lgrnge(ind,i,j) = pgn(ind)*pgx(m);
                    ++ind;
                }
            }
        }
    }
    
    /* NOW CALCULATE VALUES OF G FOR SIDE PROJECTION */
    for (i=1;i<sm+1;++i) {
        x = 2.0*(FLT) i/(FLT)(sm+1) -1.0;
        ptvalues1d(x);
        for(m=0;m<p+1;++m)
            lgrnge1d(m,i) = pgx(m);
    }
		
#ifdef DEBUG	
	Array<FLT,1> test(tm), rslt(tm), test1(tm), uht(tm); 
	
	test(0) = func(-1.0,1.0);
	test(1) = func(-1.0,-1.0);
	test(2) = func(1.0,-1.0);
	
	for (i=1;i<sm+1;++i) {
		x = 2.0*(FLT) i/(FLT)(sm+1) -1.0;
		test(i+2) = func(x,-1.0);
		test(i+2+sm) = func(-x,x);
		test(i+2+2*sm) = func(-1.,-x);
	}
	
	int count = bm;
	for(i=1;i<sm;++i) {
			for(j=1;j<sm-(i-1);++j) {
					s = -1 +2.0*((FLT) j)/(FLT)(sm+1);
					r = -1 +2.0*((FLT) i)/(FLT)(sm+1);	
					test(count++) = func(r,s);
			}
	}
	std::cout << test << std::endl;

	legtobasis(test.data(),rslt.data());
	
	std::cout << rslt << std::endl;
	
	test1(Range(0,2)) = rslt(Range(0,2));
	
	for(i=0;i<3;++i) {
		uht(0) = rslt((i+1)%3);
		uht(1) = rslt((i+2)%3);
		uht(Range(2,sm+1)) = rslt(Range(3+i*sm,3+(i+1)*sm-1));
		proj1d_leg(uht.data(),&test1(2+i*sm));
	}
	
	Array<FLT,2> d1_leg(MXGP,MXGP);
	proj_leg(rslt.data(),d1_leg.data(), MXGP);
	count = bm;
	for(i=1;i<sm;++i) {
			for(j=1;j<sm-(i-1);++j) {	
					test1(count++) = d1_leg(i,j);
			}
	}
	test -= test1;
	std::cout << test << std::endl;
#endif
	
    
    return;
}

void tri_basis::legtobasis(const FLT *data, FLT *coeff) const {
	int i,j,m,n;
	char trans[] = "T";

	
	/* Vertex coefficients are the same */
	for(int i=0;i<3;++i)
		coeff[i] = data[i];
		
	/* Side coefficients */
	TinyMatrix<FLT,MXTM,MXTM> matrix;
	TinyVector<int,2*MXTM> ipiv;
	int info;

	/* REVERSE OUTPUTING PROCESS */
	for(m=0;m<sm;++m)
		for(n=0;n<sm;++n)
			matrix(n,m) = lgrnge1d(m+2,n+1);

	GETRF(sm,sm,matrix.data(),MXTM,ipiv.data(),info);
	if (info != 0) {
		printf("DGETRF FAILED FOR INPUTING TECPLOT SIDES\n");
		exit(1);
	}
	
	TinyVector<FLT,MXTM> uht;
	TinyVector<FLT,MXTM> u1d;
	
	int count1 = 3;
	int count2 = 3;
	for (int i=0;i<3;++i) {
		uht = 0.0;
		uht(0) = data[(i+1)%3];
		uht(1) = data[(i+2)%3];
		proj1d_leg(uht.data(),u1d.data());
		for(m=0;m<sm;++m) {
			u1d(m+1) -= data[count1++];
		}
		GETRS(trans,sm,1,matrix.data(),MXTM,ipiv.data(),u1d.data()+1,MXTM,info);
		for(m=0;m<sm;++m)
			coeff[count2++] = -u1d(1+m);
	}
		
	/* Interior modes */
	for(int m=0;m<im;++m) {
		n = 0;
		for(int i=1;i<sm;++i) {
			for(int j=1;j<sm-(i-1);++j) {
				matrix(n++,m) = lgrnge(m+bm,i,j);
			}
		}
	}

	GETRF(im,im,matrix.data(),MXTM,ipiv.data(),info);
	if (info != 0) {
		printf("DGETRF FAILED FOR INPUTING TECPLOT SIDES\n");
		exit(1);
	}

	for (int i=bm;i<tm;++i)
		coeff[i] = 0.0;
		
	TinyMatrix<FLT,MXGP+2,MXGP+2> u2d;

	proj_leg(coeff,u2d.data(),MXGP+2);

	m = 0;
	for(i=1;i<sm;++i) {
		for(j=1;j<sm-(i-1);++j) {
			uht(m) = u2d(i,j) -data[m+bm];   
			++m;
		}
	}
	GETRS(trans,im,1,matrix.data(),MXTM,ipiv.data(),uht.data(),MXTM,info);
	for(m=0;m<im;++m)
		coeff[m+bm] = -uht(m);
		
	return;
}

