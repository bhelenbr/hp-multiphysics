#include"mesh.h"
#include"utilities.h"
#include<assert.h>


void mesh::copy(const mesh& tgt) {
   int i,j,n;
      
   if (!initialized) {
/*		VERTEX STORAGE ALLOCATION */
      maxvst = (int) (tgt.maxvst);
      vrtx = (FLT (*)[ND]) xmalloc(maxvst*ND*sizeof(FLT));
      vlngth = new FLT[maxvst];
      vinfo = new int[maxvst+1];
      ++vinfo;  //  ALLOWS US TO ACCES VINFO[-1]
      nnbor = new int[maxvst];
      vtri = new int[maxvst];
	
/*		VERTEX BOUNDARY STORAGE INFORMATION */		
      nvbd = tgt.nvbd;
      maxvbel = tgt.maxvbel;
      for(i=0;i<nvbd;++i)
         vbdry[i].el = new int[maxvbel];

/*		SIDE STORAGE ALLOCATION */	
      svrtx  = (int (*)[2]) xmalloc(maxvst*2*sizeof(int));
      stri   = (int (*)[2]) xmalloc(maxvst*2*sizeof(int));
      sinfo = new int[maxvst+1];
      ++sinfo; // ALLOWS US TO ACCESS SINFO[-1]


/*		SIDE BOUNDARY STORAGE ALLOCATION */
      nsbd = tgt.nsbd;
      maxsbel = tgt.maxsbel;
      for(i=0;i<nsbd;++i)
         sbdry[i].el = new int[maxsbel];      

/*		TRIANGLE ALLOCATION */			
      tvrtx = (int (*)[3]) xmalloc(maxvst*3*sizeof(int));
      ttri = (int (*)[3]) xmalloc(maxvst*3*sizeof(int));
      tside = new struct tsidedata[maxvst];
      tinfo = new int[maxvst+1];
      ++tinfo; // ALLOWS US TO ACCESS TINFO(-1)
      initialized = 1;
   }
   else {
/*		CHECK IF BIG ENOUGH */
      assert(maxvst >= tgt.maxvst);
      assert(maxvbel >= tgt.maxvbel);
      assert(maxsbel >= tgt.maxsbel);
   }
   
/*	COPY VERTEX INFO OVER */
   nvrtx = tgt.nvrtx;
   for(i=0;i<nvrtx;++i)
      for(n=0;n<ND;++n)
         vrtx[i][n] = tgt.vrtx[i][n];

   for(i=0;i<nvrtx;++i)
      vlngth[i] = tgt.vlngth[i];

   for(i=0;i<nvrtx;++i)
      vinfo[i] = tgt.vinfo[i];

   for(i=0;i<nvrtx;++i)
      nnbor[i] = tgt.nnbor[i];
      
   for(i=0;i<nvrtx;++i)
      vtri[i] = tgt.vtri[i];
      
/*	tgt VERTEX BOUNDARY INFO */
   nvbd = tgt.nvbd;
   for(i=0;i<nvbd;++i) {
      vbdry[i].type = tgt.vbdry[i].type;
      vbdry[i].num = tgt.vbdry[i].num;
      for(j=0;j<vbdry[i].num;++j) 
         vbdry[i].el[j] = tgt.vbdry[i].el[j];
   }
   
/* tgt SIDE INFORMATION */
   nside = tgt.nside;
   for(i=0;i<nside;++i)
      for(n=0;n<2;++n)
         svrtx[i][n] = tgt.svrtx[i][n];

   for(i=0;i<nside;++i)
      for(n=0;n<2;++n)
         stri[i][n] = tgt.stri[i][n];
         
   for(i=0;i<nside;++i)
      sinfo[i] = tgt.sinfo[i];
      
/*	tgt SIDE BOUNDARY INFO */
   nsbd = tgt.nsbd;
   for(i=0;i<nsbd;++i) {
      sbdry[i].type = tgt.sbdry[i].type;
      sbdry[i].num = tgt.sbdry[i].num;
      for(j=0;j<sbdry[i].num;++j) 
         sbdry[i].el[j] = tgt.sbdry[i].el[j];
   }
   
/*	tgt ELEMENT DATA */
   ntri = tgt.ntri;
   for(i=0;i<ntri;++i)
      for(n=0;n<3;++n)
         tvrtx[i][n] = tgt.tvrtx[i][n];
 
   for(i=0;i<ntri;++i)
      for(n=0;n<3;++n)
         ttri[i][n] = tgt.ttri[i][n];
   
   for(i=0;i<ntri;++i) {
      for(n=0;n<3;++n) {
         tside[i].side[n] = tgt.tside[i].side[n];
         tside[i].sign[n] = tgt.tside[i].sign[n];
      }
   }
   
   for(i=0;i<ntri;++i)
      tinfo[i] = tgt.tinfo[i];
      
   qtree.copy(tgt.qtree);
   qtree.change_vptr(vrtx);
   
   return;  
}

