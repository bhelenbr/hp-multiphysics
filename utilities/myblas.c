#include<math.h>
#include<float.h>

void DPBSLN(double **abd, int ordr, int ofdg, double *b, int rhs) {

/* SAME AS DPBSL, BUT FOR RHS VARIABLES AT ONCE */
	static int i,j,frst,ba,n;		

	for(i=0;i<ordr;++i) {
		frst = (ofdg -i > 0 ? ofdg -i : 0);
		ba = i - ofdg +frst;
		for(j=frst;j<ofdg;++j)
			for(n=0;n<rhs;++n)
				b[rhs*i +n] -= abd[i][j]*b[(ba+j)*rhs +n];

		for(n=0;n<3;++n)
			b[i*rhs +n] /= abd[i][ofdg];
	}

	for(i=ordr-1;i>=0;--i) {
 		for(n=0;n<3;++n)
 			b[i*rhs +n] /= abd[i][ofdg];

 		frst = (ofdg -i > 0 ? ofdg -i : 0);
 		ba = i-ofdg +frst;
 		for(j=frst;j<ofdg;++j)
 			for(n=0;n<3;++n)
				b[rhs*(ba+j) +n] -= abd[i][j]*b[i*rhs +n];
	}
	
	return;
}


void BLCKTRI(int nblck,
double (*a)[2], double (*b)[2], double (*c)[2], double (*d)[2])
{

/*	SOLVES 2x2 BLOCK OF FORM (nx ny)(-ny nx) */
/* TRIDAGANOL SYSTEM (A,B,C) = D,E FOR SIDE MESH VELOCITIES */
	static int i;
	static double jcbi, temp;


	jcbi    = 1.0/(b[0][0]*b[0][0] + b[0][1]*b[0][1]);
	temp    = jcbi*(b[0][0]*c[0][0] +b[0][1]*c[0][1]);
	c[0][1] = jcbi*(b[0][0]*c[0][1] -b[0][1]*c[0][0]);
	c[0][0] = temp;
	temp    = jcbi*(b[0][0]*d[0][0] -b[0][1]*d[0][1]);
	d[0][1] = jcbi*(b[0][1]*d[0][0] +b[0][0]*d[0][1]);
	d[0][0] = temp;

	for (i=1;i<nblck;++i) {
		b[i][0] -=  a[i][0]*c[i-1][0] -a[i][1]*c[i-1][1];
		b[i][1] -=  a[i][0]*c[i-1][1] +a[i][1]*c[i-1][0];
		temp     =  a[i][0]*d[i-1][0] +a[i][1]*d[i-1][1];
		d[i][1] -= -a[i][1]*d[i-1][0] +a[i][0]*d[i-1][1];
		d[i][0] -= temp;		
		
		jcbi = 1.0/(b[i][0]*b[i][0] + b[i][1]*b[i][1]);
		temp    = jcbi*(b[i][0]*c[i][0] +b[i][1]*c[i][1]);
		c[i][1] = jcbi*(b[i][0]*c[i][1] -b[i][1]*c[i][0]);
		c[i][0] = temp;
		temp    = jcbi*(b[i][0]*d[i][0] -b[i][1]*d[i][1]);
		d[i][1] = jcbi*(b[i][1]*d[i][0] +b[i][0]*d[i][1]);
		d[i][0] = temp;		
	}
		
/*	NOW COMPLETE BACK SUBSTITUTION */
	for(i=nblck-2;i>=0;--i) {
		d[i][0] -= c[i][0]*d[i+1][0]  +c[i][1]*d[i+1][1];
		d[i][1] -= -c[i][1]*d[i+1][0] +c[i][0]*d[i+1][1];	
	}

	return;
}

void BLCKTRI2(int nblck,
double (*a)[2], double (*b)[2], double (*c)[2], double (*d)[2], double (*e)[2])
{

/*	SOLVES 2x2 BLOCK OF FORM (nx ny)(-ny nx) */
/* TRIDAGANOL SYSTEM (A,B,C) = D,E FOR SIDE MESH VELOCITIES */
	static int i;
	static double jcbi, temp;


	jcbi    = 1.0/(b[0][0]*b[0][0] + b[0][1]*b[0][1]);
	temp    = jcbi*(b[0][0]*c[0][0] +b[0][1]*c[0][1]);
	c[0][1] = jcbi*(b[0][0]*c[0][1] -b[0][1]*c[0][0]);
	c[0][0] = temp;
	temp    = jcbi*(b[0][0]*d[0][0] -b[0][1]*d[0][1]);
	d[0][1] = jcbi*(b[0][1]*d[0][0] +b[0][0]*d[0][1]);
	d[0][0] = temp;
	temp    = jcbi*(b[0][0]*e[0][0] -b[0][1]*e[0][1]);
	e[0][1] = jcbi*(b[0][1]*e[0][0] +b[0][0]*e[0][1]);
	e[0][0] = temp;	

	for (i=1;i<nblck;++i) {
		b[i][0] -=  a[i][0]*c[i-1][0] -a[i][1]*c[i-1][1];
		b[i][1] -=  a[i][0]*c[i-1][1] +a[i][1]*c[i-1][0];
		temp     =  a[i][0]*d[i-1][0] +a[i][1]*d[i-1][1];
		d[i][1] -= -a[i][1]*d[i-1][0] +a[i][0]*d[i-1][1];
		d[i][0] -= temp;
		temp     =  a[i][0]*e[i-1][0] +a[i][1]*e[i-1][1];
		e[i][1] -= -a[i][1]*e[i-1][0] +a[i][0]*e[i-1][1];
		e[i][0] -= temp;		
		
		
		jcbi = 1.0/(b[i][0]*b[i][0] + b[i][1]*b[i][1]);
		temp    = jcbi*(b[i][0]*c[i][0] +b[i][1]*c[i][1]);
		c[i][1] = jcbi*(b[i][0]*c[i][1] -b[i][1]*c[i][0]);
		c[i][0] = temp;
		temp    = jcbi*(b[i][0]*d[i][0] -b[i][1]*d[i][1]);
		d[i][1] = jcbi*(b[i][1]*d[i][0] +b[i][0]*d[i][1]);
		d[i][0] = temp;
		temp    = jcbi*(b[i][0]*e[i][0] -b[i][1]*e[i][1]);
		e[i][1] = jcbi*(b[i][1]*e[i][0] +b[i][0]*e[i][1]);
		e[i][0] = temp;		
	}
		
/*	NOW COMPLETE BACK SUBSTITUTION */
	for(i=nblck-2;i>=0;--i) {
		d[i][0] -= c[i][0]*d[i+1][0] + c[i][1]*d[i+1][1];
		d[i][1] -= -c[i][1]*d[i+1][0] +c[i][0]*d[i+1][1];
		e[i][0] -= c[i][0]*e[i+1][0] + c[i][1]*e[i+1][1];
		e[i][1] -= -c[i][1]*e[i+1][0] +c[i][0]*e[i+1][1];		
	}

	return;
}
double alga_(double x);
double t_(double y);
double gamma_(double x,int *ierr);

int recur(int n,int ipoly, double al, double be, double *a, double *b)
{
	/* System generated locals */
	double r__1, r__2, r__3;
	int ierr;

	/* Local variables */
	int k;
	double alpbe, t;
	double almach, be2, al2, fkm1;

/* This subroutine generates the coefficients  a(k),b(k), k=0,1,...,n-1, */
/* in the recurrence relation */

/*	   p(k+1)(x)=(x-a(k))*p(k)(x)-b(k)*p(k-1)(x), */
/*							k=0,1,...,n-1, */

/*	   p(-1)(x)=0,  p(0)(x)=1, */

/* for some classical (monic) orthogonal polynomials, and sets  b(0) */
/* equal to the total mass of the weight distribution. The results are */
/* stored in the arrays  a,b,  which hold, respectively, the coefficients */
/* a(k-1),b(k-1), k=1,2,...,n. */

/*	   Input:  n - - the number of recursion coefficients desired */
/*			   ipoly-int identifying the polynomial as follows: */
/*					 1=Legendre polynomial on (-1,1) */
/*					 2=Legendre polynomial on (0,1) */
/*					 3=Chebyshev polynomial of the first kind */
/*					 4=Chebyshev polynomial of the second kind */
/*					 5=Jacobi polynomial with parameters  al=-.5,be=.5 */
/*					 6=Jacobi polynomial with parameters  al,be */
/*					 7=generalized Laguerre polynomial with */
/*					   parameter  al */
/*					 8=Hermite polynomial */
/*			   al,be-input parameters for Jacobi and generalized */
/*					 Laguerre polynomials */

/*	   Output: a,b - arrays containing, respectively, the recursion */
/*					 coefficients  a(k-1),b(k-1), k=1,2,...,n. */
/*			   ierr -an error flag, equal to  0  on normal return, */
/*					 equal to  1  if  al  or  be  are out of range */
/*					 when  ipoly=6  or  ipoly=7, equal to  2  if  b(0) */
/*					 overflows when  ipoly=6  or  ipoly=7, equal to  3 */
/*					 if  n  is out of range, and equal to  4  if  ipoly */
/*					 is not an admissible int. In the case  ierr=2, */
/*					 the coefficient  b(0)  is set equal to the largest */
/*					 machine-representable number. */

/* The subroutine calls for the function subroutines  r1mach,gamma  and */
/* alga. The routines  gamma  and  alga, which are included in this file, */
/* evaluate respectively the gamma function and its logarithm for */
/* positive arguments. They are used only in the cases  ipoly=6  and */
/* ipoly=7. */

	/* Function Body */
	if (n < 1) {
	return(3);
	}
	almach = log(DBL_MAX);
	ierr = 0;
	for (k = 0; k < n; ++k) {
	a[k] = 0.;
	}
	
	if (ipoly == 1) {
	b[0] = 2.;
	if (n == 1) {
		return(0);
	}

	for (k = 1; k < n; ++k) {
		fkm1 = k;
		b[k] = 1. / (4. - 1. /(fkm1 * fkm1));
	}	
	return 0;
	
	} else if (ipoly == 2) {
	a[0] = .5;
	b[0] = 1.;
	if (n == 1) {
		return 0 ;
	}
	
	for (k = 1; k <= n; ++k) {
		a[k] = .5;
		fkm1 =  k;
		b[k] = .25 / (4. - 1. / (fkm1 * fkm1));
	}
	return 0;
	
	} else if (ipoly == 3) {
	b[0] = atan(1.) * 4.;
	if (n == 1) {
		return 0;
	}
	b[1] = .5;
	if (n == 2) {
		return 0;
	}
	for (k = 2; k < n; ++k) {
		b[k] = .25;
	}
	return 0;
	} else if (ipoly == 4) {
	b[0] = atan(1.) * 2.;
	if (n == 1) {
		return 0;
	}
	for (k = 1; k < n; ++k) {
		b[k] = .25;
	}
	return 0;
	} else if (ipoly == 5) {
	b[0] = atan(1.) * 4.;
	a[0] =  .5;
	if (n == 1) {
		return 0;
	}
	for (k = 1; k < n; ++k) {
		b[k] =  .25;
	}
	return 0;
	} else if (ipoly == 6) {
	if (al <=  -1. || be <=  -1.) {
		return 1;
	} else {
		alpbe = al + be;
		a[0] = (be - al) / (alpbe +  2.);
		r__1 = al +  1.;
		r__2 = be +  1.;
		r__3 = alpbe +  2.;
		t = (alpbe +  1.) * log(2.) + alga_(r__1) + alga_(
			r__2) - alga_(r__3);

		if (t > almach) {
		ierr = 2;
		b[0] = DBL_MAX;
		} else {
		b[0] = exp(t);
		}
		if (n == 1) {
		return 0;
		}
		al2 = al * al;
		be2 = be * be;
		a[1] = (be2 - al2) / ((alpbe +  2.) * (alpbe +  4.));
/* Computing 2nd power */
		r__1 = alpbe +  2.;
		b[1] = (al +  1.) *  4. * (be +  1.) / ((
			alpbe +  3.) * (r__1 * r__1));
		if (n == 2) {
		return 0;
		}
		for (k = 2; k < n; ++k) {
		fkm1 =  k;
		a[k] = (be2 - al2) *  .25 / (fkm1 * fkm1 * (alpbe * (
			float).5 / fkm1 +  1.) * ((alpbe +  2.) * 
			 .5 / fkm1 +  1.));
/* Computing 2nd power */
		r__1 = alpbe *  .5 / fkm1 +  1.;
		b[k] = (al / fkm1 +  1.) *  .25 * (be / fkm1 + (
			float)1.) * (alpbe / fkm1 +   1.) / (((alpbe + (
			float)1.) *  .5 / fkm1 +  1.) * ((alpbe - 
			 1.) *  .5 / fkm1 +  1.) * (r__1 * 
			r__1));
		}
		return 0;
	}
	} else if (ipoly == 7) {
	if (al <=  -1.) {
		return 1;
	} else {
		a[0] = al +  1.;
		r__1 = al +  1.;
		b[0] = gamma_(r__1, &ierr);
		if (ierr == 2) {
		b[0] = DBL_MAX;
		}
		if (n == 1) {
		return 0;
		}
		for (k = 0; k < n; ++k) {
		fkm1 = k;
		a[k] = fkm1 *  2. + al +  1.;
		b[k] = fkm1 * (fkm1 + al);
		}
		return 0;
	}
	} else if (ipoly == 8) {
	b[0] = sqrt(atan( 1.) *  4.);
	if (n == 1) {
		return 0;
	}
	for (k = 1; k <= n; ++k) {
		b[k] = k *  .5;
	}
	return 0;
	} else {
	ierr = 4;
	}
	return(ierr);
} /* recur_ */

double alga_(double x)
{
	/* Initialized data */

	double cnum[8] = {  4.12084318584777, 85.68982062831317,
		   243.175243524421, -261.7218583856145, 
		-922.2613728801522, -517.6383498023218, 
		-77.41064071332953, -2.208843997216182 };
	double cden[8] = {  1., 45.64677187585908, 
		377.8372484823942, 951.323597679706, 
		846.0755362020782, 262.308347026946, 
		24.43519662506312, .4097792921092615 };

	/* System generated locals */
	int i__1;
	double ret_val, r__1;
	
	/* Local variables */
	double sden, snum;
	int k, m;
	double p;
	int m0;
	double xe, xi;
	int mm1;


/* This is an auxiliary function subroutine (not optimized in any */
/* sense) evaluating the logarithm of the gamma function for positive */
/* arguments  x. It is called by the subroutine  gamma. The int  m0 */
/* in the first executable statement is the smallest int  m  such */
/* that  1*3*5* ... *(2*m+1)/(2**m)  is greater than or equal to the */
/* largest machine-representable number. The routine is based on a */
/* rational approximation valid on [.5,1.5] due to W.J. Cody and */
/* K.E. Hillstrom; see Math. Comp. 21, 1967, 198-203, in particular the */
/* case  n=7  in Table II. For the computation of  m0  it calls upon the */
/* function subroutines  t  and  r1mach. The former, appended below, */
/* evaluates the inverse function  t = t(y)  of  y = t ln t. */


/* The constants in the statement below are  exp(1.)  and  .5alog(8.). */

	r__1 = (log(DBL_MAX) -  1.03972) /  2.71828;
	m0 = t_(r__1) *  2.71828;
	xi = (int)(x); /* was trunc(x) needed to fix */
	if (x - xi >  .5) {
		xi +=  1.;
	}
	m = (int) xi - 1;
	
/* Computation of log gamma on the standard interval (1/2,3/2] */

	xe = x - (double) m;
	snum = cnum[0];
	sden = cden[0];
	for (k = 2; k <= 8; ++k) {
		snum = xe * snum + cnum[k - 1];
		sden = xe * sden + cden[k - 1];
	}
	ret_val = (xe -  1.) * snum / sden;

/* Computation of log gamma on (0,1/2] */
	if (m == -1) {
		ret_val -= log(x);
		return ret_val;
	} else if (m == 0) {
		return ret_val;
	} else {

/* Computation of log gamma on (3/2,5/2] */

		p = xe;
		if (m == 1) {
			ret_val += log(p);
			return ret_val;
		} else {

/* Computation of log gamma for arguments larger than 5/2 */

			mm1 = m - 1;

/* The else-clause in the next statement is designed to avoid possible */
/* overflow in the computation of  p  in the if-clause, at the expense */
/* of computing many logarithms. */
			if (m < m0) {
				i__1 = mm1;
				for (k = 1; k <= i__1; ++k) {
					p = (xe + (double) k) * p;
				}
				ret_val += log(p);
				return ret_val;
			} else {
				ret_val += log(xe);
				i__1 = mm1;
				for (k = 1; k <= i__1; ++k) {
					ret_val += log(xe + (double) k);
				}
				return ret_val;
			}
		}
	}
} /* alga_ */
double gamma_(double x,int *ierr)
{
	/* System generated locals */
	double ret_val;

	/* Local variables */
	double t;
	double almach;


/* This evaluates the gamma function for real positive  x, using the */
/* function subroutines  alga  and  r1mach. In case of overflow, the */
/* routine returns the largest machine-representable number and the */
/* error flag  ierr=2. */

	almach = log(DBL_MAX);
	*ierr = 0;
	t = alga_(x);
	if (t >= almach) {
	*ierr = 2;
	ret_val = DBL_MAX;
	return ret_val;
	} else {
	ret_val = exp(t);
	return ret_val;
	}
} /* gamma_ */

double t_(double y)
{
	/* System generated locals */
	double ret_val;

	/* Local variables */
	double p, z__;


/* This evaluates the inverse function  t = t(y)  of y = t ln t  for */
/* nonnegative  y  to an accuracy of about one percent. For the */
/* approximation used, see pp. 51-52 in W. Gautschi,Computational */
/* aspects of three-term recurrence relations'', SIAM Rev. 9, 1967, */
/* 24-82. */

	if (y <=  10.) {
	p = y *  5.7941e-5 -  .00176148;
	p = y * p +  .0208645;
	p = y * p -  .129013;
	p = y * p +  .85777;
	ret_val = y * p +  1.0125;
	} else {
	z__ = log(y) -  .775;
	p = ( .775 - log(z__)) / (z__ +  1.);
	p =  1. / (p +  1.);
	ret_val = y * p / z__;
	}
	return ret_val;
} /* t_ */








int radau(n, alpha, beta, end, zero, weight, e, a, b)
int n;
double *alpha, *beta, end, *zero, *weight;
double *e, *a, *b;
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    int k ,ierr;
    double epsma;
    extern int gauss(int, const double *,const double *, 
                     double, double *, double *,double *);
    double p0, p1;
    double pm1;
    int np1;


/* Given  n  and a measure  dlambda, this routine generates the */
/* (n+1)-point Gauss-Radau quadrature formula */

/*   integral over supp(dlambda) of f(t)dlambda(t) */

/*     = w(0)f(x(0)) + sum from k=1 to k=n of w(k)f(x(k)) + R(n;f). */

/* The nodes are returned as  zero(k)=x(k), the weights as  weight(k) */
/* =w(k), k=0,1,2,...,n. The user has to supply the recursion */
/* coefficients  alpha(k), beta(k), k=0,1,2,...,n, for the measure */
/* dlambda. The nodes and weights are computed as eigenvalues and */
/* in terms of the first component of the respective normalized */
/* eigenvectors of a slightly modified Jacobi matrix of order  n+1. */
/* To do this, the routine calls upon the subroutine  gauss. It also */
/* uses the function subroutine  r1mach. */

/*    Input:  n - -  the number of interior points in the Gauss-Radau */
/*                   formula; type int */
/*            alpha,beta - arrays of dimension  n+1  to be supplied with */
/*                   the recursion coefficients  alpha(k-1), beta(k-1), */
/*                   k=1,2,...,n+1; the coefficient  alpha(n+1)  is not */
/*                   used by the routine */
/*            end -  the prescribed endpoint  x(0)  of the Gauss-Radau */
/*                   formula; type real */

/*    Output: zero - array of dimension  n+1  containing the nodes (in */
/*                   increasing order)  zero(k)=x(k), k=0,1,2,...,n */
/*            weight-array of dimension  n+1  containing the weights */
/*                   weight(k)=w(k), k=0,1,2,...,n */
/*            ierr - an error flag inherited from the routine  gauss */

/* The arrays  e,a,b  are needed for working space. */


/* The arrays  alpha,beta,zero,weight,e,a,b  are assumed to have */
/* dimension  n+1. */

    /* Parameter adjustments */
    --b;
    --a;
    --e;
    --weight;
    --zero;
    --beta;
    --alpha;

    /* Function Body */
    epsma = DBL_EPSILON;

/* epsma is the machine single precision. */

    np1 = n + 1;
    i__1 = np1;
    for (k = 1; k <= i__1; ++k) {
	a[k] = alpha[k];
	b[k] = beta[k];
/* L10: */
    }
    p0 = 0.;
    p1 = 1.;
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
	pm1 = p0;
	p0 = p1;
	p1 = (end - a[k]) * p0 - b[k] * pm1;
/* L20: */
    }
    a[np1] = end - b[np1] * p0 / p1;
    ierr = gauss(np1, &a[1], &b[1], epsma, &zero[1], &weight[1], &e[1]);
    return ierr;
} /* radau_ */

int gauss(int n, const double *alpha,const double *beta, double eps,
		   double *zero, double *weight,double *e)
{
	/* Local variables */
	double b, c__, f, g;
	int i__, j, k, l, m;
	double p, r__, s;
	int ii, mml;


/* Given  n  and a measure  dlambda, this routine generates the n-point */
/* Gaussian quadrature formula */

/*	 integral over supp(dlambda) of f(x)dlambda(x) */

/*		= sum from k=1 to k=n of w(k)f(x(k)) + R(n;f). */

/* The nodes are returned as  zero(k)=x(k) and the weights as */
/* weight(k)=w(k), k=1,2,...,n. The user has to supply the recursion */
/* coefficients  alpha(k), beta(k), k=0,1,2,...,n-1, for the measure */
/* dlambda. The routine computes the nodes as eigenvalues, and the */
/* weights in term of the first component of the respective normalized */
/* eigenvectors of the n-th order Jacobi matrix associated with  dlambda. */
/* It uses a translation and adaptation of the algol procedure  imtql2, */
/* Numer. Math. 12, 1968, 377-383, by Martin and Wilkinson, as modified */
/* by Dubrulle, Numer. Math. 15, 1970, 450. See also Handbook for */
/* Autom. Comput., vol. 2 - Linear Algebra, pp.241-248, and the eispack */
/* routine  imtql2. */

/*		Input:  n - - the number of points in the Gaussian quadrature */
/*					  formula; type integer */
/*				alpha,beta - - arrays of dimension  n  to be filled */
/*					  with the values of  alpha(k-1), beta(k-1), k=1,2, */
/*					  ...,n */
/*				eps - the relative accuracy desired in the nodes */
/*					  and weights */

/*		Output: zero- array of dimension  n  containing the Gaussian */
/*					  nodes (in increasing order)  zero(k)=x(k), k=1,2, */
/*					  ...,n */
/*				weight - array of dimension  n  containing the */
/*					  Gaussian weights  weight(k)=w(k), k=1,2,...,n */
/*				ierr- an error flag equal to  0  on normal return, */
/*					  equal to  i  if the QR algorithm does not */
/*					  converge within 30 iterations on evaluating the */
/*					  i-th eigenvalue, equal to  -1  if  n  is not in */
/*					  range, and equal to  -2  if one of the beta's is */
/*					  negative. */

/* The array  e  is needed for working space. */

	if (n < 1) {
		return -1;
	}
	zero[0] = alpha[0];
	if (beta[0] < 0.) {
		return -2;
	}
	weight[0] = beta[0];
	if (n == 1) {
		return 0;
	}
	weight[0] =  1.;
	e[n-1] =  0.;
	for (k = 1; k < n; ++k) {
		zero[k] = alpha[k];
		if (beta[k] <  0.) {
			return -2;
		}
		e[k - 1] = sqrt(beta[k]);
		weight[k] =  0.;
	}
	
	for (l = 0; l < n; ++l) {
		j = 0;

/* Look for a small subdiagonal element. */
L105:	
		for (m = l; m < n; ++m) {
			if (m == n-1)
				goto L120;
			if ( fabs(e[m]) <= (eps * (fabs(zero[m]) + fabs(zero[m + 1]) )) )
				goto L120;

		}
L120:
		p = zero[l];
		if (m == l)
			goto L240;
		if (j == 30)
			goto L400;
		++j;

/* Form shift. */

		g = (zero[l + 1] - p) / (e[l] *  2.);
		r__ = sqrt(g * g +  1.);
		g = zero[m] - p + e[l] / (g + copysign(r__, g));
		s =  1.;
		c__ =  1.;
		p =  0.;
		mml = m - l;
		
/* For i=m-1 step -1 until l do ... */

		for (ii = 0; ii < mml; ++ii) {
			i__ = m - ii -1;
			f = s * e[i__];
			b = c__ * e[i__];
			if (fabs(f) < fabs(g))
				goto L150;
			c__ = g / f;
			r__ = sqrt(c__ * c__ +  1.);
			e[i__ + 1] = f * r__;
			s =  1. / r__;
			c__ *= s;
			goto L160;
L150:
			s = f / g;
			r__ = sqrt(s * s +  1.);
			e[i__ + 1] = g * r__;
			c__ =  1. / r__;
			s *= c__;
L160:
			g = zero[i__ + 1] - p;
			r__ = (zero[i__] - g) * s + c__ *  2. * b;
			p = s * r__;
			zero[i__ + 1] = g + p;
			g = c__ * r__ - b;

/* Form first component of vector. */

			f = weight[i__ + 1];
			weight[i__ + 1] = s * weight[i__] + c__ * f;
			weight[i__] = c__ * weight[i__] - s * f;
/* L200: */
		}
		zero[l] -= p;
		e[l] = g;
		e[m] =  0.;
		goto L105;
L240:
	;
	}
	for (l = 0; l < n; ++l) 
/* Order eigenvalues and eigenvectors. */


	for (ii = 1; ii < n; ++ii) {
	i__ = ii - 1;
	k = i__;
	p = zero[i__];
	for (j = ii; j < n; ++j) {
		if (zero[j] >= p)
			goto L260;
		k = j;
		p = zero[j];
L260:
		;
	}
	if (k == i__) {
		goto L300;
	}
	zero[k] = zero[i__];
	zero[i__] = p;
	p = weight[i__];
	weight[i__] = weight[k];
	weight[k] = p;
L300:
	;
	}
	for (k = 0; k < n; ++k) {
	weight[k] = beta[0] * weight[k] * weight[k];
/* L310: */
	}
	return 0;

/* Set error - no convergence to an eigenvalue after 30 iterations. */

L400:
	return l;
} /* gauss_ */
