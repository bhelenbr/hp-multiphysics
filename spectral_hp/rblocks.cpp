/*
 *  rblocks.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on Mon Dec 31 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "blocks.h"
   
void blocks::r_jacobi(int niter, int lvl) {
   int i,iter;
      
   /*****************************************/
   /* JACOBI-ITERATION FOR MESH POSITION ****/
   /*****************************************/

   for(i=0;i<nblocks;++i)
      blk[i].grd[lvl].r_mesh::vddt();

#ifdef FOURTH
   for(i=0;i<nblocks;++i)
      blk[i].grd[lvl].r_mesh::vddt1_mp();   
   
   for(i=0;i<nblocks;++i)
      blk[i].grd[lvl].r_mesh::vddt1();
#endif

   for(i=0;i<nblocks;++i)
      blk[i].grd[lvl].r_mesh::vddt_mp();     
   
   for(i=0;i<nblocks;++i)
      blk[i].grd[lvl].r_mesh::vddti();

   for(iter=0;iter<niter;++iter) {
   
      for(i=0;i<nblocks;++i)
         blk[i].grd[lvl].r_mesh::rsdl();

#ifdef FOURTH
      for(i=0;i<nblocks;++i)
         blk[i].grd[lvl].r_mesh::rsdl1_mp();     
      
      for(i=0;i<nblocks;++i)
         blk[i].grd[lvl].r_mesh::rsdl1();
#endif

      for(i=0;i<nblocks;++i)
         blk[i].grd[lvl].r_mesh::rsdl_mp();

      for(i=0;i<nblocks;++i) 
         blk[i].grd[lvl].r_mesh::update();
   }

   return;
}

void blocks::r_cycle(int vw, int lvl) {
   int vcount,j,crscntr=0;
#ifdef SOLVECOARSE
   FLT mxr[NV], emax = 0.0, err;
#endif
   
   for (vcount=0;vcount<vw;++vcount) {
      r_jacobi(ndown,lvl);
   
      if (lvl == mgrids-1) {
#ifdef SOLVECOARSE
         /* DO A GOOD JOB ON COARSEST MESH */
         if (mgrids != 1) {

            for(int i=0;i<nblocks;++i)
               err = blk[i].grd[lvl].r_mesh::maxres(mxr);
            
            emax = MAX(emax,err);
            --vcount;
            if (err/emax < 3.0e-1 || err < 1.0e-11 || crscntr > 100) {
               // printf("# r_mesh coarsest grid iterations %d\n",crscntr);
               vcount+=2;
            }
            ++crscntr;
         }
#endif      
         continue;
      }
      
      for(j=0;j<nblocks;++j)
         blk[j].grd[lvl].r_mesh::rsdl();
         
#ifdef FOURTH
      for(j=0;j<nblocks;++j)
         blk[j].grd[lvl].r_mesh::rsdl1_mp();      

      for(j=0;j<nblocks;++j)
         blk[j].grd[lvl].r_mesh::rsdl1();
#endif
      for(j=0;j<nblocks;++j)
         blk[j].grd[lvl].r_mesh::rsdl_mp(); 

      for(j=0;j<nblocks;++j)
         blk[j].grd[lvl+1].r_mesh::mg_getfres();
      
      r_cycle(vw, lvl+1);

      for(j=0;j<nblocks;++j)
         blk[j].grd[lvl].r_mesh::mg_getcchng();
   }
   
   r_jacobi(nup,lvl);

   return;
}

void blocks::r_ksrc() {
   int i,j;

#ifdef GEOMETRIC   
   /* SETUP SPRING CONSTANTS  */
   for(i=0;i<mgrids;++i) {
      for(j=0;j<nblocks;++j)
         blk[j].grd[i].rklaplace();
      
      for(j=0;j<nblocks;++j)
         blk[j].grd[i].kvol_mp();
               
      for(j=0;j<nblocks;++j)
         blk[j].grd[i].kvoli();
   }
#else
   /* USE MULTIGRID INTERPOLATION (ALGEBRAIC) */
   /* MUST BE DONE THIS WAY FOR SPRING METHOD */
   /* SETUP FIRST MESH */
   for(j=0;j<nblocks;++j) 
      blk[j].grd[0].rklaplace();
   
   for(j=0;j<nblocks;++j)
      blk[j].grd[0].kvol_mp();
               
   for(j=0;j<nblocks;++j)
      blk[j].grd[0].kvoli();
   
   /* SETUP COARSE GRIDS */
   for(i=1;i<mgrids;++i) {
      for(j=0;j<nblocks;++j)
         blk[j].grd[i].rkmgrid();
      
      for(j=0;j<nblocks;++j)
         blk[j].grd[i].rkmgrid_mp();
               
      for(j=0;j<nblocks;++j)
         blk[j].grd[i].rkmgridi();
   }
#endif
      
/* CALCULATE SOURCE TERM ON FINEST MESH */
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].r_mesh::source();

#ifdef FOURTH
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].r_mesh::rsdl1_mp();  

   for(i=0;i<nblocks;++i)
      blk[i].grd[0].r_mesh::rsdl1();
#endif
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].r_mesh::rsdl_mp();
   
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].r_mesh::sumsrc();
   
   return;
}

