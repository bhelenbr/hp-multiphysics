#include"blocks.h"
#include"utilities.h"
#include<assert.h>
#include<string.h>

extern FLT f1(int n, FLT x, FLT y);

class hp_mgrid block::temp_hp;

void block::initialize(char *inputfile, int grds, class hpbasis *bin, int lg2p) {
   FILE *fp;
   char grd_nm[200];
   int i,j,nsurf,surfid,bnum,sflag;
   int fmt;
   FLT grwfac;
   
   fp = fopen(inputfile,"r");
   if (fp == NULL) {
      printf("couldn't open file %s\n",inputfile);
      exit(1);
   }
   
   fscanf(fp,"%*[^\n]%s\n",grd_nm);
   printf("#GRIDNAME\n#%s\n",grd_nm);
   fscanf(fp,"%*[^\n]%d %lf\n",&fmt,&grwfac);
   printf("#FORMAT\t\tGROWTH\n#%d\t\t\t%.1f\n",fmt,grwfac);

   /* LOAD MESH & CREATE COARSER MESHES */   
   ngrid = grds;
   grd = new class hp_mgrid[ngrid];
   grd[0].in_mesh(grd_nm,static_cast<FILETYPE>(fmt),grwfac);
   reconnect();
      
   /* ALLOCATE R-DEFORMABLE MESH STORAGE */
   grd[0].r_mesh::allocate(0,&rgbl);
   for(i = 1; i< ngrid; ++i)
      grd[i].r_mesh::allocate(1,&rgbl);

   /* ALLOCATE SPECTRAL_HP VECTORS */
   hpbase = bin;
   lg2pmax = lg2p;
   
   grd[0].spectral_hp::allocate(&hpbase[lg2pmax]);
   grd[0].init_comm_buf(3*NV*(hpbase[lg2pmax].sm +2));
   for(i=1;i<ngrid;++i) {
      grd[i].spectral_hp::allocate(&hpbase[0]);
      grd[i].init_comm_buf(NV*3);
   }

   /* INITIALIZE HP_MGRID STORAGE FOR MULTIGRID SOLUTION */ 
   grd[0].allocate(0,&gbl);  // 0 DENOTES FINEST LEVEL ALLOCATES GLOBAL STORAGE
   for(i = lg2pmax -1; i >= 0; --i) {
      grd[0].loadbasis(&hpbase[i]);
      grd[0].allocate(1,&gbl); // 1 DENOTES MGRID LEVEL
   }
   grd[0].loadbasis(&hpbase[lg2pmax]);
   for(i=1;i<ngrid;++i)
      grd[i].allocate(1,&gbl);

   /* LOADS RHO, MU */
   fscanf(fp,"%*[^\n]%lf %lf\n",&gbl.rho,&gbl.mu);
   printf("#RHO\t\tMU\n#%.3f\t%.3f\n",gbl.rho,gbl.mu);
   gbl.rhoi = 1.0/gbl.rho;
   gbl.nu = gbl.mu/gbl.rho;
   
   /* LOAD FUNCTION INFO??? */
   gbl.func = &f1;
   
   sflag = 0;
   for(i=0;i<grd[0].nsbd;++i)
      if (grd[0].sbdry[i].type & (IFCE_MASK +FSRF_MASK)) sflag = 1;

   if (sflag) {
      /* LOAD SURFACE INFORMATION */
      fscanf(fp,"%*[^\n]%d\n",&nsurf);
      printf("#FIRST SURFACES\n#%d\n",nsurf);
      
      /* SET SURFACE POINTERS TO NULL FOR UNUSED SURFACES */
      for(i=0;i<grd[0].nsbd;++i) 
         if(grd[0].sbdry[i].type&(CURV_MASK+IFCE_MASK))
            for(j=0;j<ngrid;++j)
               grd[j].sbdry[i].misc = NULL;
   
      /* ALLOCATE GLOBAL STORAGE FOR EACH FIRST SURFACE */   
      sgbl = new struct surface_glbls[nsurf];
                  
      /* CREATE NSURF SURFACE OBJECTS FOR EACH GRID */
      for (j=0;j<ngrid;++j)
         grd[j].srf = new class surface[nsurf];
   
      /* ALLOCATE SURFACE OBJECTS */
      for(i=0;i<nsurf;++i) {
         fscanf(fp,"%*[^\n]%d\n",&surfid);
         printf("#SURFACE ID\n#%d\n",surfid);
   
         for(bnum=0;bnum<grd[0].nsbd;++bnum) 
            if (grd[0].sbdry[bnum].type == surfid) break;
         
         if (bnum==grd[0].nsbd) {
            printf("Side ID doesn't match input file\n");
            exit(1);
         }
   
         /* ALLOCATE FINE SURFACE STORAGE */ 
         grd[0].srf[i].alloc(grd[0].maxsbel, lg2pmax, 0, 1, &sgbl[i]);
         for(j = lg2pmax -1; j >= 0; --j) {
            grd[0].srf[i].alloc(grd[0].maxsbel, j, 1, 1, &sgbl[i]);
         }
         grd[0].sbdry[bnum].misc = static_cast<void *>(&grd[0].srf[i]);
   
         /* ALLOCATE COARSE SURFACE STORAGE */
         for(j=1;j<ngrid;++j) {
            grd[j].srf[i].alloc(grd[j].maxsbel, 0, 1, 0, &sgbl[i]);
            grd[j].sbdry[bnum].misc = static_cast<void *>(&grd[j].srf[i]);
         }
         
         /* READ SURFACE PHYSICS INFO */
         fscanf(fp,"%*[^\n]%lf %lf %lf\n",&sgbl[i].sigma,&sgbl[i].rho2,&sgbl[i].mu2);
         printf("#SIGMA\t\tRHO2\t\tMU2\n#%f\t\t%f\t\t%f\n",sgbl[i].sigma,sgbl[i].rho2,sgbl[i].mu2);
      }
   }

   grd[0].loadbasis(&hpbase[lg2pmax]);

   return;
}

#ifdef OLDRECONNECT
void block::reconnect() {
   int i;
   
#ifdef ALIGNSQUARE
   FLT xmax[ND];
   
   for (int n=0;n<ND;++n) 
      xmax[n] = grd[0].vrtx[0][n];
      
   for(int k=0;k<grd[0].nvrtx;++k) 
      for(int n=0;n<ND;++n) 
         xmax[n] = MAX(xmax[n],grd[0].vrtx[k][n]);
   
   /* RESCALE TO SQUARE */
   for(int k=0;k<grd[0].nvrtx;++k) 
      for(int n=0;n<ND;++n) 
         grd[0].vrtx[k][n] /= xmax[n];   
#endif
         

   for(i = 1; i< ngrid; ++i) {
    
      grd[i].coarsen(1.6,grd[i-1]);
      
#ifdef ALIGNSQUARE
      /* TO ALIGN EDGES ON CARTESIAN */
      for(int k=0;k<grd[i].nvrtx;++k)
         grd[i].vrtx[k][1] += 0.01*grd[i].vrtx[k][0];
      grd[i].swap();
      for(int k=0;k<grd[i].nvrtx;++k)
         grd[i].vrtx[k][1] -= 0.01*grd[i].vrtx[k][0];
#endif
      grd[i].setbcinfo();
      grd[i].setfine(grd[i-1]);
      grd[i-1].setcoarse(grd[i]);
   }
   
#ifdef ALIGNSQUARE
   /* RESCALE BACK */
   for(i=0;i<ngrid;++i)
      for(int k=0;k<grd[i].nvrtx;++k) 
         for(int n=0;n<ND;++n) 
            grd[i].vrtx[k][n] *= xmax[n];
#endif

   return;
}
#else
void block::reconnect() {
   int i;

   for(i = 1; i< ngrid; ++i) {
      grd[i].coarsen2(1.5,grd[i-1],temp_hp);
      grd[i].setbcinfo();
      grd[i].setfine(grd[i-1]);
      grd[i-1].setcoarse(grd[i]);
   }

   return;
}
#endif



void block::coarsenchk(const char *fname) {
   int i;
   char name[100];

   for(i = 0; i< ngrid; ++i) {
      number_str(name,fname,i,1);
      grd[i].mesh::setbcinfo();
      grd[i].checkintegrity();
      grd[i].out_mesh(name,grid);
      grd[i].setbcinfo();
      grd[i].testconnect(name);
   }
   return;
}




