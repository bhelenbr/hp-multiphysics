/*
 *  tobasis.cpp
 *  planar++
 *
 *  Created by helenbrk on Fri Oct 12 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include"spectral_hp.h"
#include<myblas.h>

void spectral_hp::tobasis(struct vsi g, FLT (*func)(int, FLT, FLT)) {
	int tind,i,j,m,n,indx,v0,v1,sind,info;
   char uplo[] = "U";
   
/*	LOOP THROUGH VERTICES */
   for(i=0;i<nvrtx;++i)
      for(n=0;n<NV;++n)
         g.v[i][n] = func(n,vrtx[i][0],vrtx[i][1]);
         
   if (b.sm <= 0) return;

/*	LOOP THROUGH SIDES */	
	for(sind=0;sind<nside;++sind) {
            
      v0 = svrtx[sind][0];
      v1 = svrtx[sind][1];
      
      if (sinfo[sind] < 0) {
         for(n=0;n<ND;++n)
            b.proj1d(vrtx[v0][n],vrtx[v1][n],crd[n][0]);
      }
      else {
         crdtouht1d(sind);
         for(n=0;n<ND;++n)
            b.proj1d(uht[n],crd[n][0]);
      }
      
      for(n=0;n<NV;++n)
         b.proj1d(g.v[v0][n],g.v[v1][n],res[n][0]);

      for(i=0;i<b.gpx; ++i)
         for(n=0;n<NV;++n)
            res[n][0][i] -= (*func)(n,crd[0][0][i],crd[1][0][i]);
            
      for(n=0;n<NV;++n)
         b.intgrt1d(res[n][0],lf[n]);
   
      indx = sind*sm0;
      for(n=0;n<NV;++n) {
         PBTRS(uplo,b.sm,b.sbwth,1,b.sdiag1d[0],b.sbwth+1,&lf[n][2],b.sm,info);
         for(m=0;m<b.sm;++m) 
            g.s[indx+m][n] = -lf[n][2+m];
      }
	}
	
	if (b.im <= 0) return;
   
   for(tind = 0; tind < ntri; ++tind) {
		ugtouht_bdry(tind);
		for(n=0;n<NV;++n)
			b.proj_bdry(uht[n],u[n]);
         
      if (tinfo[tind] < 0) {
         for(n=0;n<ND;++n)
            b.proj(vrtx[tvrtx[tind][0]][n],vrtx[tvrtx[tind][1]][n],vrtx[tvrtx[tind][2]][n],crd[n]);
      }
      else {
         crdtouht(tind);
         for(n=0;n<ND;++n)
            b.proj_bdry(uht[n],crd[n]);
      }
         
      for(n=0;n<NV;++n)
         for (i=0; i < b.gpx; ++i )
            for (j=0; j < b.gpn; ++j )
               u[n][i][j] -= (*func)(n,crd[0][i][j],crd[1][i][j]);
                     
      indx = tind*im0;
      for(n=0;n<NV;++n) {
			b.intgrt(u[n],lf[n]);
			DPBTRS(uplo,b.im,b.ibwth,1,b.idiag[0],b.ibwth+1,&lf[n][b.bm],b.im,info);
			for(i=0;i<b.im;++i)
				g.i[indx+i][n] = -lf[n][b.bm+i];
		}
	}
   
	return;
}
