/*
 *  l2error.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on Tue Jun 11 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */

#include "spectral_hp.h"

void spectral_hp::l2error(FLT (*func)(int, FLT, FLT)) {
	int i,j,n,tind,loc[NV];
	FLT err,mxr[NV],l2r[NV];

	for(n=0;n<NV;++n) {
		mxr[n] = 0.0;
		l2r[n] = 0.0;
	}
	
	for(tind=0;tind<ntri;++tind) {
      
      if (tinfo[tind] > -1) {
         crdtouht(tind);
         for(n=0;n<ND;++n)
            b.proj_bdry(uht[n], crd[n], dcrd[n][0], dcrd[n][1]);
      }
      else {
         for(n=0;n<ND;++n)
            b.proj(vrtx[tvrtx[tind][0]][n],vrtx[tvrtx[tind][1]][n],vrtx[tvrtx[tind][2]][n],crd[n]);

         for(i=0;i<b.gpx;++i) {
            for(j=0;j<b.gpn;++j) {
               for(n=0;n<ND;++n) {
                  dcrd[n][0][i][j] = 0.5*(vrtx[tvrtx[tind][1]][n] -vrtx[tvrtx[tind][0]][n]);
                  dcrd[n][1][i][j] = 0.5*(vrtx[tvrtx[tind][2]][n] -vrtx[tvrtx[tind][0]][n]);
               }
            }
         }
      }

      ugtouht(tind);
		for(n=0;n<NV;++n)
			b.proj(uht[n],u[n]);
		
 		for (i=0;i<b.gpx;++i) {	
			for (j=0;j<b.gpn;++j) {
            cjcb[i][j] = (dcrd[0][0][i][j]*dcrd[1][1][i][j] -dcrd[1][0][i][j]*dcrd[0][1][i][j]);
            for(n=0;n<NV;++n) {
               err =  fabs(u[n][i][j]-func(n,crd[0][i][j],crd[1][i][j]));
               if (err > mxr[n]) {
                  mxr[n] = err;
                  loc[n] = tind;
               }
               l2r[n] += err*err*b.wtx[i]*b.wtn[j]*cjcb[i][j];
            }
         }
      }	
	}
   
   printf("# ERROR: ");
	for(n=0;n<NV;++n) {
		l2r[n] = sqrt(l2r[n]); 
		printf("%.3e %.3e %4d ",l2r[n],mxr[n],loc[n]);
	}
	printf("\n");
		
	return;
}

