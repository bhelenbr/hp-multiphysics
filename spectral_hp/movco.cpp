#include"hp_mgrid.h"

void hp_mgrid::getfres() {
   int i,j,k,m,n,tind,v0,indx,indx1;
   class hp_mgrid *fmesh;
   
   isfrst = true;
   
   /* TRANSFER COUPLED SURFACE RESIDUALS */
   for(i=0;i<nsbd;++i)
      if(sbdry[i].type&(FSRF_MASK+IFCE_MASK))
         surfgetfres(i);

   if(p0 > 1) {
      /* TRANSFER IS ON FINE MESH */
      for(i=0;i<nvrtx;++i)
         for(n=0;n<NV;++n)
            gbl->res0.v[i][n] = gbl->res.v[i][n];

      if (b.p > 1) {
         indx = 0;
         indx1 = 0;
         for(i=0;i<nside;++i) {
            for (j=0;j<b.sm;++j) {
               for(n=0;n<NV;++n)
                  gbl->res0.s[indx][n] = gbl->res.s[indx1][n];
               ++indx;
               ++indx1;
            }
            indx1 += b.p;
         }
         
         if (b.p > 2) {
            indx = 0;
            indx1 = 0;
            for(tind=0;tind<ntri;++tind) {
               for(m=1;m<b.sm;++m) {
                  for(k=0;k<b.sm-m;++k) {
                     for(n=0; n<NV; ++n)
                        gbl->res0.i[indx][n] = gbl->res.i[indx1][n];
                     ++indx;
                     ++indx1;
                  }
                  indx1 += b.p;
               }
            }
         }
      }
      
      return;
   }
   
   fmesh = static_cast<class hp_mgrid *>(fmpt);
   
   /* TRANSFER IS BETWEEN DIFFERENT MESHES */
   for(i=0;i<nvrtx;++i)
      for(n=0;n<NV;++n)
         gbl->res0.v[i][n] = 0.0;
         
   /* LOOP THROUGH FINE VERTICES TO CALCULATE RESIDUAL  */
   for(i=0;i<fmesh->nvrtx;++i) {
      tind = fmesh->coarse[i].tri;
      for(j=0;j<3;++j) {
         v0 = tvrtx[tind][j];
         for(n=0;n<NV;++n)
            gbl->res0.v[v0][n] += fmesh->coarse[i].wt[j]*gbl->res.v[i][n];
      }
   }
   
   /* LOOP THROUGH COARSE VERTICES   */
   /* TO CALCULATE VUG ON COARSE MESH */
   for(i=0;i<nvrtx;++i) {
      tind = fine[i].tri;

      for(n=0;n<NV;++n)
         ug.v[i][n] = 0.0;
         
      for(j=0;j<3;++j) {
         for(n=0;n<NV;++n)
            ug.v[i][n] += fine[i].wt[j]*fmesh->ug.v[fmesh->tvrtx[tind][j]][n];
      }
      
      for(n=0;n<NV;++n)
         vug_frst[i][n] = ug.v[i][n];
   }

   return;
}

