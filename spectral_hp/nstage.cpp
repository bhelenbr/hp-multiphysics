#include"hp_mgrid.h"

void hp_mgrid::nstage1(void)
{
	static int i,n;

/******************************************************************/
/* BEGINNING OF MULTISTAGE PSEUDO-TIME-UPDATE OF FLOW & MESH POSITION ****/
/******************************************************************/
   for(i=0;i<nvrtx;++i) 
      for(n=0;n<NV;++n)
         gbl.vug0[i][n] = vug[i][n];

	if (b.p > 1) {
      for(i=0;i<nside*sm0;++i) 
			for(n=0;n<NV;++n)
				gbl.sug0[i][n] = sug[i][n];
      
      if (b.p > 2) {
         for(i=0;i<ntri*im0;++i)
            for(n=0;n<NV;++n)
               gbl.iug0[i][n] = iug[i][n];
      }   
	}

   return;
}

/********************************/
/* INTERMEDIATE & FINAL STEPS ***/
/********************************/
void hp_mgrid::nstage2(int stage) {
   int i,m,k,n,indx,indx1;
   FLT cflalpha;

   cflalpha = gbl.flowcfl[log2p]*alpha[stage];
   
   for(i=0;i<nvrtx;++i)
      for(n=0;n<NV;++n)
         vug[i][n] = gbl.vug0[i][n] -cflalpha*gbl.vres[i][n];

   if (b.sm > 0) {
      indx = 0;
      indx1 = 0;
      for(i=0;i<nside;++i) {
         for (m=0;m<b.sm;++m) {
            for(n=0;n<NV;++n)
               sug[indx1][n] = gbl.sug0[indx1][n] -cflalpha*gbl.sres[indx][n];
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
                     iug[indx1][n] =  gbl.iug0[indx1][n] -cflalpha*gbl.ires[indx][n];
                  }
                  ++indx; ++indx1;
               }
               indx1 += sm0 -b.sm;
            }
         }
      }
   }
         
   return;
}