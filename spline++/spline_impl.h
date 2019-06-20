/*
 *  spline.cpp
 *  spline++
 *
 *  Created by Brian Helenbrook on 4/16/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "spline.h"
#include <fstream>
#include <myblas.h>
//#include <blitz/tinyvec-et.h>

template<int ND> int spline<ND>::read(std::string filename) {
	std::ifstream in;
	Array<TinyVector<double,ND>,1> yd,ydd;

    in.open(filename.c_str());
    if (!in) {
        std::cerr << "couldn't open spline input file " << filename << std::endl;
        exit(1);
    }

	/* READ # OF POINTS AND ALLOCATE */
	in.ignore(160,'\n');
	in.ignore(160,':');
	in >> npts;
			
	c.resize(npts);
	x.resize(npts);
	y.resize(npts);
	yd.resize(npts);
	ydd.resize(npts);
	in.ignore(80,'\n');
	in.ignore(80,'\n');
			
	/* NOW READ DATA */
	for (int i=0;i<npts;++i) {
		in >> x(i);

		for(int n=0;n<ND;++n)
			in >> y(i)(n);
			
	}	
	
	/* BASIC SPLINE FIRST (NO CUSPS NOT SMOOTH LOOP) */
	TinyVector<double,ND> d1,d2,d3,d4;
	double h1, h2, h3, h4;
	double h02, h24, h12, h23, h13;
	TinyVector<double,ND> yddbk, yddc, yddfw;
	TinyVector<double,ND> ydbk, ydc, ydfw;
	double a1, b1, c1;
	double sum;
	
	/* POINT 0 */
	h1=x(1)-x(0);
	d1=(y(1)-y(0))/h1;
	h2=x(2)-x(1);
	d2=(y(2)-y(1))/h2;
	ydd(0)=(d2-d1)/(0.5*(h1+h2));
	yd(0)=d1-0.5*h1*ydd(0);
	
	/* POINT 1 */
	yddc=ydd(0);
	ydc=(d1*h2+h1*d2)/(h1+h2);
	yddbk=yddc;
	ydbk=ydc;
	h3 = x(2)-x(1);
	d3 = (y(2)-y(1))/h3;
	h4 = x(3)-x(2);
	d4 = (y(3)-y(2))/h4;
	yddfw = 2.*(d4-d3)/(h4+h3);
	ydfw = d3 -0.5*h3*yddfw;
	h02 = 2*h1;
	h24 = h3+h4;
	h12 = h1;
	h23 = h2;
	
	for (int n=0;n<ND;++n) {
		a1=h23*h24*fabs(yddfw(n)*yddc(n));
		b1=2.*h02*h24*fabs(yddfw(n)*yddbk(n));
		c1=h12*h02*fabs(yddc(n)*yddbk(n));
		sum = a1 + b1 + c1;
	  
		if (sum == 0.) {
			sum = 1.0;
			if (yddc(n) == 0.0) {
				b1 = 1.;
				a1 = 0.0;
				c1 = 0.0;
			}
			else {
				b1 = 0.0;
				h13 = h12+h23;
				c1 = h12/h13;
				a1 = h23/h13;
			}
		}
		yd(1)(n) = (a1*ydbk(n) + b1*ydc(n) + c1*ydfw(n))/sum;
		ydd(1)(n) = (a1*yddbk(n) + b1*yddc(n) +c1*yddfw(n))/sum;
	}
		
	/* INTERIOR POINTS */
	for (int i=2;i<npts-2;++i) {
		h1 = x(i-1)-x(i-2);
		d1 = (y(i-1)-y(i-2))/h1;
		
		h2 = x(i)-x(i-1);
		d2 = (y(i)-y(i-1))/h2;

		h3 = x(i+1)-x(i);
		d3 = (y(i+1)-y(i))/h3;

		h4 = x(i+2)-x(i);
		d4 = (y(i+2)-y(i))/h4;

		
		yddbk = 2.*(d2-d1)/(h2+h1);
		yddc = 2.*(d3-d2)/(h3+h2);
		yddfw = 2.*(d4-d3)/(h4+h3);
		
		ydbk = d2 +0.5*h2*yddbk;
		ydc = (d2*h3+d3*h2)/(h2+h3);
		ydfw = d3 -0.5*h3*yddfw;
		
		h02 = h1+h2;
		h24 = h3+h4;
		h12 = h2;
		h23 = h3;
		
		for(int n=0;n<ND;++n) {
			a1=h23*h24*fabs(yddfw(n)*yddc(n));
			b1=2.*h02*h24*fabs(yddfw(n)*yddbk(n));
			c1=h12*h02*fabs(yddc(n)*yddbk(n));
			sum = a1 + b1 + c1;
		  
			if (sum == 0.) {
				sum = 1.0;
				if (yddc(n) == 0.0) {
					b1 = 1.;
					a1 = 0.0;
					c1 = 0.0;
				}
				else {
					b1 = 0.0;
					h13 = h12+h23;
					c1 = h12/h13;
					a1 = h23/h13;
				}
			}
			yd(i)(n) = (a1*ydbk(n) + b1*ydc(n) + c1*ydfw(n))/sum;
			ydd(i)(n) = (a1*yddbk(n) + b1*yddc(n) +c1*yddfw(n))/sum;
		}
	}

	/* POINT NPTS-1 */
	h4=x(npts-1)-x(npts-2);
	d4=(y(npts-1)-y(npts-2))/h4;
	h3=x(npts-2)-x(npts-3);
	d3=(y(npts-2)-y(npts-2))/h3;
	ydd(npts-1)=(d4-d3)/(0.5*(h4+h3));
	yd(npts-1)=d4+0.5*h4*ydd(npts-1);
	
	/* POINT NPTS-2 */
	yddc=ydd(npts-1);
	ydc=(d3*h4+h3*d4)/(h4+h3);
	yddfw=yddc;
	ydfw=ydc;
	h1 = x(npts-3)-x(npts-4);
	d1 = (y(npts-3)-y(npts-4))/h1;
	h2 = x(npts-2)-x(npts-3);
	d2 = (y(npts-2)-y(npts-3))/h2;
	yddbk = 2.*(d2-d1)/(h2+h1);
	ydbk = d2 +0.5*h2*yddbk;
	
	h02 = h1+h2;
	h24 = 2.*h4;
	h12 = h2;
	h23 = h4;
	
	for (int n=0;n<ND;++n) {
		a1=h23*h24*fabs(yddfw(n)*yddc(n));
		b1=2.*h02*h24*fabs(yddfw(n)*yddbk(n));
		c1=h12*h02*fabs(yddc(n)*yddbk(n));
		sum = a1 + b1 + c1;
	  
		if (sum == 0.) {
			sum = 1.0;
			if (yddc(n) == 0.0) {
				b1 = 1.;
				a1 = 0.0;
				c1 = 0.0;
			}
			else {
				b1 = 0.0;
				h13 = h12+h23;
				c1 = h12/h13;
				a1 = h23/h13;
			}
		}
		yd(npts-2)(n) = (a1*ydbk(n) + b1*ydc(n) + c1*ydfw(n))/sum;
		ydd(npts-2)(n) = (a1*yddbk(n) + b1*yddc(n) +c1*yddfw(n))/sum;
	}
	
	double a,b,bma,bmas2;
	TinyVector<double,ND> tem, t1, t2, t3;
	
	/* NOW DEFINE QUINTIC POLYNOMIAL COEFFICIENTS */
	for(int i=0;i<npts-1;++i) {
		a=x(i);
		b=x(i+1);
		bma=b-a;
		bmas2=bma*bma*.5;
		c(i)(0)=y(i);
		c(i)(1)=bma*yd(i);
		c(i)(2)=bmas2*ydd(i);
		tem=c(i)(1)+c(i)(2);
		t1=y(i+1)-c(i)(0)-tem;
		t2=bma*yd(i+1)-tem-c(i)(2);
		t3=bmas2*ydd(i+1)-c(i)(2);
		c(i)(5)=6.*t1-3.*t2+t3;
		c(i)(3)=4.*t1+c(i)(5)-t2;
		c(i)(4)=t1-c(i)(3)-c(i)(5);
	}
	
	return 0;
}
	
	
template<int ND> int spline<ND>::interpolate(double xptin, TinyVector<double,ND>& loc) {
	double a,b,bma,z;
	
	double xpt = xptin;
	if (xpt < x(0)) xpt=x(0);
	if (xpt > x(npts-1)) xpt=x(npts-1);

	int i;
	for (i=1;i<npts;++i)
		if (x(i) >= xpt) break;
	--i;
	
	a=x(i);
	b=x(i+1);
	bma=b-a;
	z=(xpt-a)/bma;
	loc=c(i)(0)+z*(c(i)(1)+z*(c(i)(2)+z*(c(i)(3)+z*(c(i)(4)+z*c(i)(5)))));

	return 0;
}

template<int ND> int spline<ND>::tangent(double xptin, TinyVector<double,ND>& tan) {
    double a,b,bma,z;
    
    double xpt = xptin;
    if (xpt < x(0)) xpt=x(0);
    if (xpt > x(npts-1)) xpt=x(npts-1);
    
    int i;
    for (i=1;i<npts;++i)
        if (x(i) >= xpt) break;
    --i;
    
    a=x(i);
    b=x(i+1);
    bma=b-a;
    z=(xpt-a)/bma;
    tan=(c(i)(1)+z*(2*c(i)(2)+z*(3*c(i)(3)+z*(4*c(i)(4)+z*5*c(i)(5)))))/bma;
    return 0;
}

template<int ND> int spline<ND>::curvature(double xptin, TinyVector<double,ND>& curv) {
    double a,b,bma,z;
    
    double xpt = xptin;
    if (xpt < x(0)) xpt=x(0);
    if (xpt > x(npts-1)) xpt=x(npts-1);
    
    int i;
    for (i=1;i<npts;++i)
        if (x(i) >= xpt) break;
    --i;
    
    a=x(i);
    b=x(i+1);
    bma=b-a;
    z=(xpt-a)/bma;
    curv=(2*c(i)(2)+z*(6*c(i)(3)+z*(12*c(i)(4)+z*20*c(i)(5))))/(bma*bma);
    return 0;
}



template<int ND> int spline<ND>::find(double& s0, TinyVector<double,ND>& ypt) {
    int k,sidloc,sidlocprev=0;
    double ol,psi,normdist;
    double psiloc,psiprev,normdistprev=1.0e32;
    double mindist = 1.0e32;
	TinyVector<double,ND> y1, dy, dy1, tan, curv;

    psiloc = 1.0;
    sidloc = npts-2;
	psiprev = -1.0;
	
    for(k=0;k<npts-1;++k) {
        dy = y(k+1) -y(k);
        ol = 2./dot(dy,dy);
		dy1 = ypt -y(k);
        psi = ol*(dot(dy1,dy)) -1.;
        normdist = dy(0)*dy1(1)-dy(1)*dy1(0);
        normdist *= sqrt(ol/2.);
        
        if (psi <= -1.0 && psiprev >= 1.0) {
            /* PREVIOUS & THIS SIDE ARE POTENTIAL MATCHES */
            if (fabs(normdist) < mindist) {
                mindist = fabs(normdist);
                sidloc = k;
                psiloc = -1.0;
            }
            if (fabs(normdistprev) < mindist) {
                mindist = fabs(normdistprev);
                sidloc = sidlocprev;
                psiloc = 1.0;
            }
        }
        else if (psi >= -1.0 && psi <= 1.0) {
            /* POTENTIAL SIDE */
            if (fabs(normdist) < mindist) {
                mindist = fabs(normdist);
                sidloc = k;
                psiloc = psi;
            }
        }
        psiprev = psi;
        normdistprev = normdist;
        sidlocprev = k;
    }
	
	int i = sidloc;
    assert(0 <= i && i < npts-1);
	dy = y(i+1)-y(i);
	ol = 2./dot(dy,dy);

	psi = psiloc;
	s0 = 0.5*((x(i+1)+x(i)) +psi*(x(i+1)-x(i)));
    
    /* Find s location of minimum distance) */
    for (int iter=0;iter<100;++iter) {
        interpolate(s0,y1);
        dy = y1-ypt;
        mindist = dot(dy,dy);
        tangent(s0,tan);
        curvature(s0,curv);
        double dmindist = 2*dot(tan,dy);
        double dminddist2 = 2*(dot(tan,tan)+dot(curv,dy));
        double ds = -dmindist/dminddist2;
        s0 += ds;
        //std::cout << iter << ' ' << s0 << ' ' << ds << ' ' << dmindist << std::endl;
        if (fabs(ds) < 1.0e-10) {
            ypt = y1;
            return 0;
        }
    }
    ypt = y1;
    
	return(1);
}


template<int ND> int spline3<ND>::read(std::string filename) {
	std::ifstream in;
	Array<TinyVector<double,3>,1> tri_diagonal;
	TinyVector<Array<double,1>,ND> rhs;
	
	in.open(filename.c_str());
	if (!in) {
		std::cerr << "couldn't open spline input file " << filename << std::endl;
		exit(1);
	}
	
	/* READ # OF POINTS AND ALLOCATE */
	in.ignore(80,'\n');
	in.ignore(80,':');
	in >> npts;
	
	c.resize(npts);
	x.resize(npts);
	y.resize(npts);
	tri_diagonal.resize(npts);
	for (int n=0;n<ND;++n)
		rhs(n).resize(npts);
	
	in.ignore(80,'\n');
	in.ignore(80,'\n');
	
	/* NOW READ DATA */
	for (int i=0;i<npts;++i) {
		in >> x(i);

		for(int n=0;n<ND;++n)
			in >> y(i)(n);
		
	}
	
	// http://en.wikipedia.org/wiki/Spline_interpolation#Algorithm_to_find_the_interpolating_cubic_spline
	for (int row=1;row<npts-1;++row) {
		tri_diagonal(row)(0) = 1./(x(row)-x(row-1));
		tri_diagonal(row)(1) = 2.*(1./(x(row)-x(row-1)) +1./(x(row+1)-x(row)));
		tri_diagonal(row)(2) = 1./(x(row+1)-x(row));

		for(int n=0;n<ND;++n)
			rhs(n)(row) = 3*((y(row)(n) -y(row-1)(n))/pow(x(row)-x(row-1),2) +(y(row+1)(n) -y(row)(n))/pow(x(row+1)-x(row),2));
	}
	
	/* Can either constrain slope or do natural spline */  
	/* Here we do natural spline */

	tri_diagonal(0)(0) = 0.0;
	tri_diagonal(0)(1) = 2./(x(1)-x(0));
	tri_diagonal(0)(2) = 1./(x(1)-x(0));
	for(int n=0;n<ND;++n)
		rhs(n)(0) = 3*((y(1)(n) -y(0)(n))/pow(x(1)-x(0),2));

	tri_diagonal(npts-1)(0) = 1./(x(npts-1)-x(npts-2));
	tri_diagonal(npts-1)(1) = 2./(x(npts-1)-x(npts-2));
	tri_diagonal(npts-1)(2) = 0.0;
	for(int n=0;n<ND;++n)
		rhs(n)(npts-1) = 3*((y(npts-1)(n) -y(npts-2)(n))/pow(x(npts-1)-x(npts-2),2));

	/* PUT IN BANDED FORM & INVERT */
	Array<double,2> band_matrix(npts,2);
	for(int j=0;j<npts;++j) {
		int i1 = (0 > j-1 ? 0 : j-1);
		for(int i=i1;i<=j;++i) {
			int k = i - (j-1);
			band_matrix(j,k) = tri_diagonal(j)(k);
		}
	}

	int info;
	char uplo[] = "U";
#ifdef F2CFortran
	DPBTRF(uplo,npts,1,band_matrix.data(),2,info);
#else
    const int one = 1, two = 2;
    dpbtrf_(uplo,&npts,&one,band_matrix.data(),&two,&info);
#endif
	if (info != 0) {
		printf("1:PBTRF FAILED info: %d\n", info);
		exit(1);
	}
	
    for (int n=0;n<ND;++n) {
#ifdef F2CFortran
        DPBTRS(uplo,npts,1,1,band_matrix.data(),2,rhs(n).data(),npts,info);
#else
        const int one = 1, two = 2;
        dpbtrs_(uplo,&npts,&one,&one,band_matrix.data(),&two,rhs(n).data(),&npts,&info);
#endif

    }
	
	for(int seg=0;seg<npts-1;++seg) {
		for(int n=0;n<ND;++n) {
			c(seg)(0)(n) = y(seg)(n);
			c(seg)(1)(n) = y(seg+1)(n);
			c(seg)(2)(n) = rhs(n)(seg)*(x(seg+1)-x(seg)) -(y(seg+1)(n)-y(seg)(n));
			c(seg)(3)(n) = -rhs(n)(seg+1)*(x(seg+1)-x(seg)) +(y(seg+1)(n)-y(seg)(n));
		}
	}
		
	return 0;
}

template<int ND> int spline3<ND>::interpolate(double xptin, TinyVector<double,ND>& loc) {
	double a,b,bma,t;
	
	double xpt = xptin;
	if (xpt < x(0)) xpt=x(0);
	if (xpt > x(npts-1)) xpt=x(npts-1);
	
	int i;
	for (i=1;i<npts;++i)
		if (x(i) >= xpt) break;
	--i;
	
	a=x(i);
	b=x(i+1);
	bma=b-a;
	t=(xpt-a)/bma;
	loc=c(i)(0)*(1-t) + c(i)(1)*t +t*(1-t)*(c(i)(2)*(1-t) + c(i)(3)*t);
	
	return 0;
}

template<int ND> int spline3<ND>::tangent(double xptin, TinyVector<double,ND>& tan) {
    double a,b,bma,t;
    
    double xpt = xptin;
    if (xpt < x(0)) xpt=x(0);
    if (xpt > x(npts-1)) xpt=x(npts-1);
    
    int i;
    for (i=1;i<npts;++i)
        if (x(i) >= xpt) break;
    --i;
    
    a=x(i);
    b=x(i+1);
    bma=b-a;
    t=(xpt-a)/bma;
    tan=  (c(i)(1)-c(i)(0) +(1-2*t)*(c(i)(2)*(1-t) + c(i)(3)*t) +t*(1-t)*(c(i)(3)-c(i)(2)))/bma;
    
    return 0;
}

template<int ND> int spline3<ND>::curvature(double xptin, TinyVector<double,ND>& curv) {
    double a,b,bma,t;
    
    double xpt = xptin;
    if (xpt < x(0)) xpt=x(0);
    if (xpt > x(npts-1)) xpt=x(npts-1);
    
    int i;
    for (i=1;i<npts;++i)
        if (x(i) >= xpt) break;
    --i;
    
    a=x(i);
    b=x(i+1);
    bma=b-a;
    t=(xpt-a)/bma;
    curv=  (-2*(c(i)(2)*(1-t) + c(i)(3)*t) +2*(1-2*t)*(c(i)(3)-c(i)(2)))/(bma*bma);
    
    return 0;
}

template<int ND> int spline3<ND>::find(double& s0, TinyVector<double,ND>& ypt) {
    int k,sidloc,sidlocprev=0;
    double ol,psi,normdist;
    double psiloc,psiprev,normdistprev=1.0e32;
    double mindist = 1.0e32;
    TinyVector<double,ND> y1, dy, dy1, tan, curv;
    
    psiloc = 1.0;
    sidloc = npts-2;
    psiprev = -1.0;
    
    for(k=0;k<npts-1;++k) {
        dy = y(k+1) -y(k);
        ol = 2./dot(dy,dy);
        dy1 = ypt -y(k);
        psi = ol*(dot(dy1,dy)) -1.;
        normdist = dy(0)*dy1(1)-dy(1)*dy1(0);
        normdist *= sqrt(ol/2.);
        
        if (psi <= -1.0 && psiprev >= 1.0) {
            /* PREVIOUS & THIS SIDE ARE POTENTIAL MATCHES */
            if (fabs(normdist) < mindist) {
                mindist = fabs(normdist);
                sidloc = k;
                psiloc = -1.0;
            }
            if (fabs(normdistprev) < mindist) {
                mindist = fabs(normdistprev);
                sidloc = sidlocprev;
                psiloc = 1.0;
            }
        }
        else if (psi >= -1.0 && psi <= 1.0) {
            /* POTENTIAL SIDE */
            if (fabs(normdist) < mindist) {
                mindist = fabs(normdist);
                sidloc = k;
                psiloc = psi;
            }
        }
        psiprev = psi;
        normdistprev = normdist;
        sidlocprev = k;
    }
    
    int i = sidloc;
    assert(0 <= i && i < npts-1);
    dy = y(i+1)-y(i);
    ol = 2./dot(dy,dy);
    
    psi = psiloc;
    s0 = 0.5*((x(i+1)+x(i)) +psi*(x(i+1)-x(i)));
    
    /* Find s location of minimum distance) */
    for (int iter=0;iter<100;++iter) {
        interpolate(s0,y1);
        dy = y1-ypt;
        mindist = dot(dy,dy);
        tangent(s0,tan);
        curvature(s0,curv);
        double dmindist = 2*dot(tan,dy);
        double dminddist2 = 2*(dot(tan,tan)+dot(curv,dy));
        double ds = -dmindist/dminddist2;
        s0 += ds;
        // std::cout << iter << ' ' << s0 << ' ' << ds << ' ' << dmindist << std::endl;
        if (fabs(ds) < 1.0e-10) {
            ypt = y1;
            return 0;
        }
    }
    ypt = y1;
    
    return(1);
}



	
