#include"hp_mgrid.h"

int hp_mgrid::setcoarse(class hp_mgrid& tgt) {
   cmesh = &tgt;
   if (coarse == NULL)
      coarse = new struct mg_trans[maxvst];
   return(mgconnect(coarse,tgt));
}

void hp_mgrid::getcchng(void) {
   int i,j,n,ind,tind;
	
	if(b.p > 1) {
		return;
	}

/* DETERMINE CORRECTIONS ON COARSE MESH   */   
   for(i=0;i<cmesh->nvrtx;++i)
      for(n=0;n<ND;++n) 
         cmesh->vug_frst[i][n] -= cmesh->vug[i][n];

/* LOOP THROUGH FINE VERTICES   */
/* TO DETERMINE CHANGE IN SOLUTION */   
   for(i=0;i<nvrtx;++i) {
      
      for(n=0;n<NV;++n)
         gbl.vres[i][n] = 0.0;
      
      tind = coarse[i].tri;
      
      for(j=0;j<3;++j) {
         ind = cmesh->tvrtx[tind][j];
         for(n=0;n<NV;++n) 
            gbl.vres[i][n] -= coarse[i].wt[j]*cmesh->vug_frst[ind][n];
      }
   }
   
   for(i=0;i<nvrtx;++i)
      for(n=0;n<ND;++n) 
         vug[i][n] += gbl.vres[i][n];
                  
   return;
}


	
	
	
	
		
		
