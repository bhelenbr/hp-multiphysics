#include"hp_mgrid.h"

void hp_mgrid::getcchng(void) {
   int i,j,n,ind,tind;
   class hp_mgrid *cmesh;
   
   if(b.p > 1) {
      return;
   }
   
   cmesh = static_cast<class hp_mgrid *>(cmpt);
   
   /* TRANFER COUPLED BOUNDARY RESIDUALS */
   for(i=0;i<nsbd;++i)
      if (sbdry[i].type&(FSRF_MASK+IFCE_MASK))
         surfgetcchng(i);

   /* DETERMINE CORRECTIONS ON COARSE MESH   */   
   for(i=0;i<cmesh->nvrtx;++i)
      for(n=0;n<NV;++n) 
         cmesh->vug_frst[i][n] -= cmesh->ug.v[i][n];
   
#ifdef SKIP
   /* SMOOTH CORRECTIONS??? */
   int iter,sind,v0,v1; 
   int niter = 10;  
   for(iter=0; iter< niter; ++iter) {
      /* SMOOTH POINT DISTRIBUTION X*/
      for(i=0;i<cmesh->nvrtx;++i)
         for(n=0;n<NV;++n)
            gbl->res0.v[i][n] = 0.0;

      for(i=0;i<cmesh->nside;++i) {
         v0 = cmesh->svrtx[i][0];
         v1 = cmesh->svrtx[i][1];
         for(n=0;n<NV;++n) {
            gbl->res0.v[v0][n] += cmesh->vug_frst[v1][n];
            gbl->res0.v[v1][n] += cmesh->vug_frst[v0][n];
         }
      }
      
      /* RESET BOUNDARY VALUES SO THEY DON'T CHANGE */
      for(i=0;i<cmesh->nsbd;++i) {
         for(j=0;j<cmesh->sbdry[i].num;++j) {
            sind = cmesh->sbdry[i].el[j];
            v0 = cmesh->svrtx[sind][0];
            for(n=0;n<NV;++n) {
               gbl->res0.v[v0][n] = cmesh->vug_frst[v0][n]*cmesh->nnbor[i];
            }
         }
      }
      

      for(i=0;i<cmesh->nvrtx;++i) {
         for(n=0;n<NV;++n) {
            cmesh->vug_frst[i][n] = gbl->res0.v[i][n]/cmesh->nnbor[i];
         }
      }
   }
#endif

   /* LOOP THROUGH FINE VERTICES   */
   /* TO DETERMINE CHANGE IN SOLUTION */   
   for(i=0;i<nvrtx;++i) {
      
      for(n=0;n<NV;++n)
         gbl->res.v[i][n] = 0.0;
      
      tind = coarse[i].tri;
      
      for(j=0;j<3;++j) {
         ind = cmesh->tvrtx[tind][j];
         for(n=0;n<NV;++n) 
            gbl->res.v[i][n] -= coarse[i].wt[j]*cmesh->vug_frst[ind][n];
      }
   }
   
   
   /* ADD CORRECTION */
   for(i=0;i<nvrtx;++i)
      for(n=0;n<NV;++n) 
         ug.v[i][n] += gbl->res.v[i][n];

   return;
}


   
   
   
   
      
      
