#include"blocks.h"
#include<cstring>
#include<cstdio>
#include<utilities.h>

void blocks::init(int n, int mg, char **filename, FILETYPE filetype = easymesh, FLT grwfac = 1) {
   int i,j,k,match;

   nblocks = n;
   mglvls = mg;
   blk = new class block[nblocks];

   for (i=0;i<nblocks;++i) 
      blk[i].init(mg,filename[i],filetype,grwfac);


   /* MATCH BOUNDARIES */
   for(i=0;i<mg;++i) {
      for(j=0;j<nblocks;++j) {
         match = 0;
         for(k=0;k<nblocks;++k)
            match += blk[j].grd[i].findmatch(blk[k].grd[i]);
            
         if (match != blk[j].grd[i].alld_mp()) {
            printf("error in matching boundaries %d: %d %d\n",j,match,blk[j].grd[i].alld_mp());
            exit(1);
         }
      }
   }

   ksrc();

   return;
}

void blocks::out_mesh(char *filename, FILETYPE filetype = easymesh) {
   int i;   
   char fnmcat[80],outname[80];

   /* ASSUME FOR NOW MESHES ARE LABELED a,b,c... */
   /* I HAVEN'T FIGURED OUT HOW THIS IS GOING TO WORK IN THE TOTALLY GENERAL CASE */
   if (nblocks > 1) {
      for (i=0;i<nblocks;++i) {
         strcpy(fnmcat,filename);
         strcat(fnmcat,".");
         number_str(outname, fnmcat, i, 1);
         blk[i].grd[0].setbcinfo();
         blk[i].grd[0].out_mesh(outname,filetype);
      }
   }
   else {
      blk[0].grd[0].setbcinfo();
      blk[0].grd[0].out_mesh(filename,filetype);
   }
   
   return;
}

void blocks::out_mesh(char **filename, FILETYPE filetype = easymesh) {
   int i;   

   /* ASSUME FOR NOW MESHES ARE LABELED a,b,c... */
   /* I HAVEN'T FIGURED OUT HOW THIS IS GOING TO WORK IN THE TOTALLY GENERAL CASE */
   for (i=0;i<nblocks;++i) {
      blk[i].grd[0].setbcinfo();
      blk[i].grd[0].out_mesh(filename[i],filetype);
   }
 
   return;
}

void blocks::jacobi(int niter, int lvl) {
   int i,iter;
      
/*****************************************/
/* JACOBI-ITERATION FOR MESH POSITION ****/
/*****************************************/

   for(i=0;i<nblocks;++i)
      blk[i].grd[lvl].vddt();

#ifdef FOURTH
   for(i=0;i<nblocks;++i)
      blk[i].grd[lvl].vddt1_mp();   
   
   for(i=0;i<nblocks;++i)
      blk[i].grd[lvl].vddt1();
#endif

   for(i=0;i<nblocks;++i)
      blk[i].grd[lvl].vddt_mp();     
   
   for(i=0;i<nblocks;++i)
      blk[i].grd[lvl].vddti();

   for(iter=0;iter<niter;++iter) {
   
      for(i=0;i<nblocks;++i)
         blk[i].grd[lvl].rsdl();

#ifdef FOURTH
      for(i=0;i<nblocks;++i)
         blk[i].grd[lvl].rsdl1_mp();     
      
      for(i=0;i<nblocks;++i)
         blk[i].grd[lvl].rsdl1();
#endif

      for(i=0;i<nblocks;++i)
         blk[i].grd[lvl].rsdl_mp();

      for(i=0;i<nblocks;++i) 
         blk[i].grd[lvl].update();
   }

   return;
}

void blocks::cycle(int vw, int lvl = 0) {
   int i,j;  // DON'T MAKE THESE SCREWS UP RECURSION
   
   for (i=0;i<vw;++i) {
      jacobi(1,lvl);
      if (lvl == mglvls-1) return;
      
      for(j=0;j<nblocks;++j)
         blk[j].grd[lvl].rsdl();
         
#ifdef FOURTH
      for(j=0;j<nblocks;++j)
         blk[j].grd[lvl].rsdl1_mp();      

      for(j=0;j<nblocks;++j)
         blk[j].grd[lvl].rsdl1();
#endif
      for(j=0;j<nblocks;++j)
         blk[j].grd[lvl].rsdl_mp();      
      
      for(j=0;j<nblocks;++j)
         blk[j].grd[lvl+1].mg_getfres();
      
      cycle(vw, lvl+1);

      for(j=0;j<nblocks;++j)
         blk[j].grd[lvl].mg_getcchng();
   }

   return;
}

void blocks::ksrc() {
   int i,j;

#ifdef GEOMETRIC   
   /* SETUP SPRING CONSTANTS  */
   for(i=0;i<mglvls;++i) {
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
   for(i=1;i<mglvls;++i) {
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
      blk[i].grd[0].source();

#ifdef FOURTH
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].rsdl1_mp();  

   for(i=0;i<nblocks;++i)
      blk[i].grd[0].rsdl1();
#endif
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].rsdl_mp();
   
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].sumsrc();
   
   return;
}



   