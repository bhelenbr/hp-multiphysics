/*
 *  l2error.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on Tue Jun 11 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */

#include "spectral_hp.h"
#include <utilities.h>

void spectral_hp::l2error(FLT (*func)(int, FLT, FLT)) {
	int i,j,n,tind,loc[NV];
	FLT err,mxr[NV],l2r[NV];

	for(n=0;n<NV;++n) {
		mxr[n] = 0.0;
		l2r[n] = 0.0;
	}
	
	for(tind=0;tind<ntri;++tind) {
      
      if (tinfo[tind] > -1) {
         crdtocht(tind);
         for(n=0;n<ND;++n)
            b.proj_bdry(cht[n], crd[n], dcrd[n][0], dcrd[n][1]);
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
   
	for(n=0;n<NV;++n) {
		l2r[n] = sqrt(l2r[n]); 
		printf("#L_2: %.3e L_inf %.3e %4d ",l2r[n],mxr[n],loc[n]);
	}
	printf("\n");
		
	return;
}

FLT spectral_hp::findmax(int type, FLT (*fxy)(FLT x[ND])) {
   FLT ddpsi1, ddpsi2, psil, psir;
   FLT xp[ND], dx[ND];
   int bnum, v0, v1, sind;
   
   FLT max = 0.0;
   
   for(bnum=0;bnum<nsbd;++bnum) {
      if (!(sbdry[bnum].type&type)) continue;
      
      v0 = svrtx[sbdry[bnum].el[0]][0];
      v1 = svrtx[sbdry[bnum].el[sbdry[bnum].num-1]][1];
      /* CHECK IF PERIODIC */
      ddpsi2 = 0.0;
      if (v0 == v1) {
         sind = sbdry[bnum].el[sbdry[bnum].num-1];
         crdtocht(sind);
         b.ptprobe1d(ND, cht, xp, dx, 1.0);
         ddpsi2 = (*fxy)(dx);
      }
      else {
         for(int i=0;i<nvbd;++i) {
            if (vbdry[i].type&(COMX_MASK+COMY_MASK) && vbdry[i].el[0] == v0) {
               sind = sbdry[bnum].el[sbdry[bnum].num-1];
               crdtocht1d(sind);
               b.ptprobe1d(ND, cht, xp, dx, 1.0);
               ddpsi2 = (*fxy)(dx);
            }
         } 
      }

      for(int indx=0;indx<sbdry[bnum].num;++indx) {
         sind = sbdry[bnum].el[indx];
         crdtocht1d(sind);
         b.ptprobe1d(ND, cht, xp, dx, -1.0);
         ddpsi1 = (*fxy)(dx);
         if (ddpsi1 * ddpsi2 < 0.0) {
            v0 = svrtx[sbdry[bnum].el[indx]][0];
            max = MAX(max,(*fxy)(vrtx[v0]));
            printf("#LOCAL MAX: %e %e %e\n",vrtx[v0][0],vrtx[v0][1],(*fxy)(vrtx[v0]));
         }
         b.ptprobe1d(ND, cht, xp, dx, 1.0);
         ddpsi2 = (*fxy)(dx);
         if (ddpsi1 *ddpsi2 < 0.0) {
            /* INTERIOR MAXIMUM */
            psil = -1.0;
            psir = 1.0;
            while (psir-psil > 1.0e-10) {
               b.ptprobe1d(ND, cht, xp, dx, 0.5*(psil +psir));
               if ((*fxy)(dx)*ddpsi1 < 0.0) 
                  psir = 0.5*(psil+psir);
               else
                  psil = 0.5*(psil+psir);
            }
            max = MAX(max,(*fxy)(xp));
            printf("#LOCAL MAX: %e %e %e\n",xp[0],xp[1],(*fxy)(xp));
         }  
      }
   }
   return(max);
}

static FLT fx(FLT xp[ND]) { return(xp[0]);}
static FLT fy(FLT xp[ND]) { return(xp[1]);}
FLT spectral_hp::findmaxx(int type) {return(findmax(type,&fx));}
FLT spectral_hp::findmaxy(int type) {return(findmax(type,&fy));} 
   

