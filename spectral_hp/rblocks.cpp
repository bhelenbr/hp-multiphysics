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
   static int i,iter;
      
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

void blocks::r_cycle(int vw, int lvl = 0) {
   int i,j;
   
   for (i=0;i<vw;++i) {
      r_jacobi(1,lvl);
      if (lvl == mgrids-1) return;
      
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

   return;
}

void blocks::r_ksrc() {
   int i,j;

#define GEOMETRIC

#ifdef GEOMETRIC   
/*	SETUP SPRING CONSTANTS  */
   for(i=0;i<mgrids;++i) {
      for(j=0;j<nblocks;++j)
         blk[j].grd[i].rklaplace();
      
      for(j=0;j<nblocks;++j)
         blk[j].grd[i].kvol_mp();
               
      for(j=0;j<nblocks;++j)
         blk[j].grd[i].kvoli();
   }
#else
/*	USE MULTIGRID INTERPOLATION (ALGEBRAIC) */
/*	MUST BE DONE THIS WAY FOR SPRING METHOD */
/*	SETUP FIRST MESH */
   for(j=0;j<nblocks;++j) 
      blk[j].grd[0].rklaplace();
   
   for(j=0;j<nblocks;++j)
      blk[j].grd[0].kvol_mp();
               
	for(j=0;j<nblocks;++j)
      blk[j].grd[0].kvoli();
   
/*	SETUP COARSE GRIDS */
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

