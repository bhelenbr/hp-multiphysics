/*
 *  intgrt.cpp
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 08 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#include <hpbasis.h>
#include<stdio.h>

void hpbasis::intgrt(FLT **f, FLT *rslt) {
	static int i,j,m,n,indx;
		  
/*	INTEGRATE FUNCTION FOR EACH TEST FUNCTION	*/

	for (n = 0; n < nmodx; ++n) {
		for (j = 0; j < gpn; ++j) {
			wk0[n][j] = 0.0;
			for(i = 0; i < gpx; ++i)
				wk0[n][j] += f[i][j]*gxwtx[n][i];
		}
	}

	for(n=0;n<sm+3;++n) {
		rslt[n] = 0.0;
		for (j=0; j < gpn; ++j)
				rslt[n] += gn[n][j]*wtn[j]*wk0[n][j];
	}
	
	for(n=sm+3;n<2*sm+3;++n) {
		rslt[n] = 0.0;
		for (j=0; j < gpn; ++j)
				rslt[n] += gn[n][j]*wtn[j]*wk0[1][j];
	}	

	for(n=2*sm+3;n<bm;++n) {
		rslt[n] = 0.0;
		for (j=0; j < gpn; ++j)
				rslt[n] += gn[n][j]*wtn[j]*wk0[0][j];
	}	
	
	indx = bm;
	for(m = 3; m < sm+2; ++m) {
		for(n = 0; n < sm+2-m; ++n) {
			rslt[indx] = 0.0;
			for (j=0; j < gpn; ++j )
				rslt[indx] += gn[indx][j]*wtn[j]*wk0[m][j];
			++indx;
		}
	}
	
	return;
}


/* WARNING THIS ADDS INTEGRATION TO RESULT: RESULT IS NOT CLEARED FIRST */
void hpbasis::intgrtrs(FLT **fx, FLT **fy, FLT *rslt) {
	static int i,j,m,n,indx;
	
	for(i=0;i<gpx;++i)
		for(j=0;j<gpn;++j)
			wk2[i][j] = fx[i][j] +fy[i][j]*x0[i];
			
	for (n = 0; n < nmodx; ++n) {
		for (j = 0; j < gpn; ++j) {
			wk0[n][j] = 0.0;
			wk1[n][j] = 0.0;		
			for(i = 0; i < gpx; ++i) {
				wk0[n][j] += wk2[i][j]*dgxwtx[n][i];
				wk1[n][j] += fy[i][j]*gxwtx[n][i];
			}
		}
	}

	for (m=0; m < sm+3; ++m)
		for (j=0; j < gpn; ++j )
			rslt[m] += gnwtnn0[m][j]*wk0[m][j] +dgnwtn[m][j]*wk1[m][j];

	for (m=sm+3; m < 2*sm+3; ++m) 
		for (j=0; j < gpn; ++j )
			rslt[m] += gnwtnn0[m][j]*wk0[1][j] +dgnwtn[m][j]*wk1[1][j];
	
	for (m=2*sm+3; m < bm; ++m)
		for (j=0; j < gpn; ++j )
			rslt[m] += gnwtnn0[m][j]*wk0[0][j] +dgnwtn[m][j]*wk1[0][j];
	
	indx = bm;
	for(m = 3; m < sm+2;++m) {
		for(n = 0; n < sm+2-m; ++n) {
			for (j=0; j < gpn; ++j ) {
				rslt[indx] += gnwtnn0[indx][j]*wk0[m][j]+dgnwtn[indx][j]*wk1[m][j];
			}
			++indx;
		}
	}
	
	return;
}


void hpbasis::intgrtr(FLT **f, FLT *rslt1) {
	static int i,j,m,n,indx;
	
	for (n = 0; n < nmodx; ++n) {
		for (j = 0; j < gpn; ++j) {
			wk2[n][j] = 0.0;		
			for(i = 0; i < gpx; ++i) {
				wk2[n][j] += f[i][j]*dgx[n][i]*wtx[i];
			}
		}
	}

	for (m=0; m < sm+3; ++m) {
		for (j=0; j < gpn; ++j ) {
			rslt1[m] += wtn[j]*gn[m][j]*n0[j]*wk2[m][j];
		}
	}

	for (m=sm+3; m < 2*sm+3; ++m) {
		for (j=0; j < gpn; ++j ) {
			rslt1[m] += wtn[j]*gn[m][j]*n0[j]*wk2[1][j];
		}
	}
	
	for (m=2*sm+3; m < bm; ++m) {
		for (j=0; j < gpn; ++j ) {
			rslt1[m] += wtn[j]*gn[m][j]*n0[j]*wk2[0][j];
		}
	}
	
	indx = bm;
	for(m = 3; m < sm+2;++m) {
		for(n = 0; n < sm+2-m; ++n) {
			for (j=0; j < gpn; ++j ) {
				rslt1[indx] += wtn[j]*gn[indx][j]*n0[j]*wk2[m][j];
			}
			++indx;
		}
	}	
	
	return;
}

void hpbasis::intgrts(FLT **f, FLT *rslt2) {
	static int i,j,m,n,indx;
	
	for (n = 0; n < nmodx; ++n) {
		for (j = 0; j < gpn; ++j) {
			wk0[n][j] = 0.0;
			wk1[n][j] = 0.0;
			for(i = 0; i < gpx; ++i) {
				wk0[n][j] += f[i][j]*gx[n][i]*wtx[i];
				wk1[n][j] += f[i][j]*x0[i]*dgx[n][i]*wtx[i];
			}
		}
	}

	for (m=0; m < sm+3; ++m) {
		for (j=0; j < gpn; ++j ) {
			rslt2[m] += wtn[j]*(dgn[m][j]*wk0[m][j] 
			+gn[m][j]*wk1[m][j]*n0[j]);
		}
	}

	for (m=sm+3; m < 2*sm+3; ++m) {
		for (j=0; j < gpn; ++j ) {
			rslt2[m] += wtn[j]*(dgn[m][j]*wk0[1][j] 
			+gn[m][j]*wk1[1][j]*n0[j]);
		}	
	}

	for (m=2*sm+3; m < bm; ++m) {
		for (j=0; j < gpn; ++j ) {
			rslt2[m] += wtn[j]*(dgn[m][j]*wk0[0][j] 
			+gn[m][j]*wk1[0][j]*n0[j]);
		}	
	}
	
	indx = bm;
	for(m = 3; m < sm+2;++m) {
		for(n = 0; n < sm+2-m;++n) {
			for (j=0; j < gpn; ++j ) {
				rslt2[indx] += wtn[j]*(dgn[indx][j]*wk0[m][j] 
				+gn[indx][j]*wk1[m][j]*n0[j]);
			}
			++indx;
		}
	}
	
	return;
}
