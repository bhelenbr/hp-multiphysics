#include"hp_mgrid.h"

int hp_mgrid::setfine(class hp_mgrid& tgt) {
   fmesh = &tgt;
   if (fine == NULL)
      fine = new struct mg_trans[maxvst];
   return(mgconnect(fine,tgt));
}

void hp_mgrid::getfres() {
   int i,j,k,m,n,tind,v0,indx,indx1;
   
//   fmesh->rcv(YDIR_MP,(FLT *) rg.res,0,1,2); 
   
   
   if(p0 > 1) {
/* 	TRANSFER IS ON FINE MESH */
		for(i=0;i<nvrtx;++i)
			for(n=0;n<NV;++n)
				gbl.vres0[i][n] = gbl.vres[i][n];

      if (b.p > 1) {
         indx = 0;
         indx1 = 0;
         for(i=0;i<nside;++i) {
            for (j=0;j<b.sm;++j) {
               for(n=0;n<NV;++n)
                  gbl.sres0[indx][n] = gbl.sres[indx1][n];
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
								gbl.ires0[indx][n] = gbl.ires[indx1][n];
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
   
 /*TRANSFER IS BETWEEN DIFFERENT MESHES */
   for(i=0;i<nvrtx;++i)
      for(n=0;n<NV;++n)
         gbl.vres0[i][n] = 0.0;
         
/* LOOP THROUGH FINE VERTICES TO CALCULATE RESIDUAL  */
   for(i=0;i<fmesh->nvrtx;++i) {
      tind = fmesh->coarse[i].tri;
      for(j=0;j<3;++j) {
         v0 = tvrtx[tind][j];
         for(n=0;n<NV;++n)
            gbl.vres0[v0][n] += gbl.fadd*fmesh->coarse[i].wt[j]*gbl.vres[i][n];
      }
   }
   
/* LOOP THROUGH COARSE VERTICES   */
/* TO CALCULATE VUG ON COARSE MESH */
   for(i=0;i<nvrtx;++i) {
      tind = fine[i].tri;

      for(n=0;n<NV;++n)
         vug[i][n] = 0.0;
         
      for(j=0;j<3;++j) {
         for(n=0;n<NV;++n)
            vug[i][n] += fine[i].wt[j]*fmesh->vug[fmesh->tvrtx[tind][j]][n];
      }
      
      for(n=0;n<NV;++n)
         vug_frst[i][n] = vug[i][n];
   }

   isfrst = true;

   return;
}