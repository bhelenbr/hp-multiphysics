/*
 *  intgrt.cpp
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 08 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#include"hpbasis.h"
#include<stdio.h>

void hpbasis::intgrt(FLT **f, FLT *rslt) {
    int i,j,m,n,indx;
        
   /* INTEGRATE FUNCTION FOR EACH TEST FUNCTION   */

   for (n = 0; n < nmodx; ++n) {
      for (j = 0; j < gpn; ++j) {
         wk0[n][j] = 0.0;
         for(i = 0; i < gpx; ++i)
            wk0[n][j] += f[i][j]*gxwtx[n][i];
      }
   }

   for(n=0;n<sm+vm;++n) {
      rslt[n] = 0.0;
      for (j=0; j < gpn; ++j)
            rslt[n] += gn[n][j]*wtn[j]*wk0[n][j];
   }
   
   for(n=sm+3;n<2*sm+3;++n) {
      rslt[n] = 0.0;
      for (j=0; j < gpn; ++j)
            rslt[n] += gn[n][j]*wtn[j]*wk0[2][j];
   }   

   for(n=2*sm+3;n<bm;++n) {
      rslt[n] = 0.0;
      for (j=0; j < gpn; ++j)
            rslt[n] += gn[n][j]*wtn[j]*wk0[1][j];
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
/* THIS IS AN OPTIMIZED VERSION NOT IN USE RIGHT NOW */
void hpbasis::intgrtrs(FLT *fx, FLT *fy, int fsz, FLT *rslt) {
   int i,j,m,n,indx,indx1;
   FLT wk0[MXTM][MXTM],wk1[MXTM][MXTM],wk2[MXTM][MXTM];
   int fcnt,s1,s2,s3,nmax;
   FLT lcl_wk0, lcl_wk1;
   
   for (j = 0; j < gpn; ++j) {
      fcnt = j;
      wk2[0][j] = fx[fcnt] +fy[fcnt]*x0[0];
      lcl_wk0 = wk2[0][j]*dgxwtx[0][0];
      lcl_wk1 = fy[fcnt]*gxwtx[0][0];  
      fcnt += fsz;    
      for(i = 1; i < gpx; ++i) {
         wk2[i][j] = fx[fcnt] +fy[fcnt]*x0[i];
         lcl_wk0 += wk2[i][j]*dgxwtx[0][i];
         lcl_wk1 += fy[fcnt]*gxwtx[0][i];
         fcnt += fsz;
      }
      wk0[0][j] = lcl_wk0;
      wk1[0][j] = lcl_wk1;
   }

   for (n = 1; n < nmodx; ++n) {
      for (j = 0; j < gpn; ++j) {
         fcnt = j;
         lcl_wk0 = wk2[0][j]*dgxwtx[n][0];
         lcl_wk1 = fy[fcnt]*gxwtx[n][0];
         fcnt += fsz;    
         for(i = 1; i < gpx; ++i) {
            lcl_wk0 += wk2[i][j]*dgxwtx[n][i];
            lcl_wk1 += fy[fcnt]*gxwtx[n][i];
            fcnt += fsz;
         }
         wk0[n][j] = lcl_wk0;
         wk1[n][j] = lcl_wk1;
      }
   }
   
   s1 = sm+vm;
   s2 = 2*sm+vm;
   s3 = bm;
   for (j=0; j < gpn; ++j ) {
      for (m=0; m < s1; ++m)
            rslt[m] -= gnwtnn0[m][j]*wk0[m][j] +dgnwtn[m][j]*wk1[m][j];
   
      lcl_wk0 = wk0[2][j];
      lcl_wk1 = wk1[2][j];
      for (; m < s2; ++m) 
            rslt[m] -= gnwtnn0[m][j]*lcl_wk0 +dgnwtn[m][j]*lcl_wk1;
      
      lcl_wk0 = wk0[1][j];
      lcl_wk1 = wk1[1][j];
      for (; m < s3; ++m)
            rslt[m] -= gnwtnn0[m][j]*lcl_wk0 +dgnwtn[m][j]*lcl_wk1;
   }
   
   indx = s3;
   for(m = 3; m < sm+2;++m) {
      nmax = sm+2-m;
      for (j=0; j < gpn; ++j ) {
         indx1 = indx;
         lcl_wk0 = wk0[m][j];
         lcl_wk1 = wk1[m][j];
         for(n = 0; n < nmax; ++n) {
            rslt[indx1] -= gnwtnn0[indx1][j]*lcl_wk0+dgnwtn[indx1][j]*lcl_wk1;
            ++indx1;
         }
      }
      indx = indx1;
   }
   
   
   return;
}

/* WARNING THIS ADDS INTEGRATION TO RESULT: RESULT IS NOT CLEARED FIRST */
void hpbasis::intgrtrs(FLT **fx, FLT **fy, FLT *rslt) {
    int i,j,m,n,indx;
   
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

   for (j=0; j < gpn; ++j ) {
      for (m=0; m < sm+vm; ++m)
            rslt[m] -= gnwtnn0[m][j]*wk0[m][j] +dgnwtn[m][j]*wk1[m][j];
   
      for (; m < 2*sm+vm; ++m) 
            rslt[m] -= gnwtnn0[m][j]*wk0[2][j] +dgnwtn[m][j]*wk1[2][j];
      
      for (; m < bm; ++m)
            rslt[m] -= gnwtnn0[m][j]*wk0[1][j] +dgnwtn[m][j]*wk1[1][j];
   }
   
   indx = bm;
   for(m = 3; m < sm+2;++m) {
      for(n = 0; n < sm+2-m; ++n) {
         for (j=0; j < gpn; ++j ) {
            rslt[indx] -= gnwtnn0[indx][j]*wk0[m][j]+dgnwtn[indx][j]*wk1[m][j];
         }
         ++indx;
      }
   }
   
   return;
}


void hpbasis::intgrtr(FLT **f, FLT *rslt1) {
    int i,j,m,n,indx;
   
   for (n = 0; n < nmodx; ++n) {
      for (j = 0; j < gpn; ++j) {
         wk2[n][j] = 0.0;      
         for(i = 0; i < gpx; ++i) {
            wk2[n][j] += f[i][j]*dgx[n][i]*wtx[i];
         }
      }
   }

   for (m=0; m < sm+vm; ++m) {
      for (j=0; j < gpn; ++j ) {
         rslt1[m] -= wtn[j]*gn[m][j]*n0[j]*wk2[m][j];
      }
   }

   for (m=sm+3; m < 2*sm+3; ++m) {
      for (j=0; j < gpn; ++j ) {
         rslt1[m] -= wtn[j]*gn[m][j]*n0[j]*wk2[2][j];
      }
   }
   
   for (m=2*sm+3; m < bm; ++m) {
      for (j=0; j < gpn; ++j ) {
         rslt1[m] -= wtn[j]*gn[m][j]*n0[j]*wk2[1][j];
      }
   }
   
   indx = bm;
   for(m = 3; m < sm+2;++m) {
      for(n = 0; n < sm+2-m; ++n) {
         for (j=0; j < gpn; ++j ) {
            rslt1[indx] -= wtn[j]*gn[indx][j]*n0[j]*wk2[m][j];
         }
         ++indx;
      }
   }   
   
   return;
}

void hpbasis::intgrts(FLT **f, FLT *rslt2) {
    int i,j,m,n,indx;
   
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

   for (m=0; m < sm+vm; ++m) {
      for (j=0; j < gpn; ++j ) {
         rslt2[m] -= wtn[j]*(dgn[m][j]*wk0[m][j] 
         +gn[m][j]*wk1[m][j]*n0[j]);
      }
   }

   for (m=sm+3; m < 2*sm+3; ++m) {
      for (j=0; j < gpn; ++j ) {
         rslt2[m] -= wtn[j]*(dgn[m][j]*wk0[2][j] 
         +gn[m][j]*wk1[2][j]*n0[j]);
      }   
   }

   for (m=2*sm+3; m < bm; ++m) {
      for (j=0; j < gpn; ++j ) {
         rslt2[m] -= wtn[j]*(dgn[m][j]*wk0[1][j] 
         +gn[m][j]*wk1[1][j]*n0[j]);
      }   
   }
   
   indx = bm;
   for(m = 3; m < sm+2;++m) {
      for(n = 0; n < sm+2-m;++n) {
         for (j=0; j < gpn; ++j ) {
            rslt2[indx] -= wtn[j]*(dgn[indx][j]*wk0[m][j] 
            +gn[indx][j]*wk1[m][j]*n0[j]);
         }
         ++indx;
      }
   }
   
   return;
}
