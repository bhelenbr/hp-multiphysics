#include"hp_mgrid.h"

void hp_mgrid::nstage1(void)
{
	static int i,n;
/******************************************************************/
/* BEGINNING OF MULTISTAGE PSEUDO-TIME-UPDATE OF FLOW & MESH POSITION ****/
/******************************************************************/
   for(i=0;i<nvrtx;++i) 
      for(n=0;n<NV;++n)
         gbl->ug0.v[i][n] = ug.v[i][n];

	if (b.p > 1) {
      for(i=0;i<nside*sm0;++i) 
			for(n=0;n<NV;++n)
				gbl->ug0.s[i][n] = ug.s[i][n];
      
      if (b.p > 2) {
         for(i=0;i<ntri*im0;++i)
            for(n=0;n<NV;++n)
               gbl->ug0.i[i][n] = ug.i[i][n];
      }   
	}

/*	COUPLED BOUNDARY MOVEMENT EQUATIONS */   
   for(i=0;i<nsbd;++i)
      if (sbdry[i].type&(FSRF_MASK +IFCE_MASK))
         surfnstage1(i);
         
   return;
}

/********************************/
/* INTERMEDIATE & FINAL STEPS ***/
/********************************/
void hp_mgrid::nstage2(int stage) {
   int i,m,k,n,indx,indx1;
   FLT cflalpha;

   cflalpha = cfl[log2p]*alpha[stage];
   
   for(i=0;i<nvrtx;++i)
      for(n=0;n<NV;++n)
         ug.v[i][n] = gbl->ug0.v[i][n] -cflalpha*gbl->res.v[i][n];

   if (b.sm > 0) {
      indx = 0;
      indx1 = 0;
      for(i=0;i<nside;++i) {
         for (m=0;m<b.sm;++m) {
            for(n=0;n<NV;++n)
               ug.s[indx1][n] = gbl->ug0.s[indx1][n] -cflalpha*gbl->res.s[indx][n];
            ++indx;
            ++indx1;
         }
         indx1 += sm0 -b.sm;
      }			

      if (b.im > 0) {
         indx = 0;
         indx1 = 0;
         for(i=0;i<ntri;++i) {
            for(m=1;m<b.sm;++m) {
               for(k=0;k<b.sm-m;++k) {
                  for(n=0;n<NV;++n) {
                     ug.i[indx1][n] =  gbl->ug0.i[indx1][n] -cflalpha*gbl->res.i[indx][n];
                  }
                  ++indx; ++indx1;
               }
               indx1 += sm0 -b.sm;
            }
         }
      }
   }
   
/*	COUPLED BOUNDARY MOVEMENT EQUATIONS */   
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type&(FSRF_MASK +IFCE_MASK))
         surfnstage2(i,stage);
   }
   
/*   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type&FSRF_MASK) {
         for(k=0;k<sbdry[i].num;++k) {
            indx = sbdry[i].el[k];
            printf("%d %f %f %f %f\n",k,vrtx[svrtx[indx][0]][0],vrtx[svrtx[indx][0]][1],vrtx[svrtx[indx][1]][0],vrtx[svrtx[indx][1]][1]);
         }
      }
   }
*/
            
   return;
}
