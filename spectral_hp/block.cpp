#include"blocks.h"
#include"utilities.h"
#include<assert.h>

extern FLT f1(int n, FLT x, FLT y);

void block::initialize(char *inputfile, int grds, class hpbasis *bin, int lg2p) {
   FILE *fp;
   char grd_nm[200];
   int i,j,nsurf,surfid,bnum;
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

/*	LOAD MESH & CREATE COARSER MESHES */   
   ngrid = grds;
   grd = new class hp_mgrid[ngrid];
   grd[0].in_mesh(grd_nm,static_cast<FILETYPE>(fmt),grwfac);
   for(i = 1; i< ngrid; ++i) {
      grd[i].coarsen(grd[i-1]);
      grd[i].setfine(grd[i-1]);
      grd[i-1].setcoarse(grd[i]);
   }

/*	INITIALIZE SPECTRAL_HP VECTORS */
   hpbase = bin;
   lg2pmax = lg2p;
   
   grd[0].spectral_hp::allocate(hpbase[lg2pmax]);
   grd[0].init_comm_buf(NV*(hpbase[lg2pmax].sm +2));
   for(i=1;i<ngrid;++i) {
      grd[i].spectral_hp::allocate(hpbase[0]);
      grd[i].init_comm_buf(NV*2);
   }

/*	INITIALIZE HP_MGRID STORAGE FOR MULTIGRID SOLUTION */ 
	grd[0].allocate(0,&gbl);  // 0 DENOTES FINEST LEVEL ALLOCATES GLOBAL STORAGE
   for(i = lg2pmax -1; i >= 0; --i) {
      grd[0].loadbasis(hpbase[i]);
      grd[0].allocate(1,&gbl); // 1 DENOTES MGRID LEVEL
   }
   grd[0].loadbasis(hpbase[lg2pmax]);
   for(i=1;i<ngrid;++i)
      grd[i].allocate(1,&gbl);

/*	LOADS RHO, MU */
   fscanf(fp,"%*[^\n]%lf %lf\n",&gbl.rho,&gbl.mu);
   printf("#RHO\t\tMU\n#%.3f\t%.3f\n",gbl.rho,gbl.mu);
   gbl.rhoi = 1.0/gbl.rho;
   gbl.nu = gbl.mu/gbl.rho;
   
/* LOAD FUNCTION INFO??? */
   gbl.func = &f1;

/*	LOAD SURFACE INFORMATION */
   fscanf(fp,"%*[^\n]%d\n",&nsurf);
   printf("#FIRST SURFACES\n#%d\n",nsurf);
   
/*	SET SURFACE POINTERS TO NULL FOR UNUSED SURFACES */
   for(i=0;i<grd[0].nsbd;++i) 
      if(grd[0].sbdry[i].type&(CURV_MASK+IFCE_MASK))
         for(j=0;j<ngrid;++j)
            grd[j].sbdry[i].misc = NULL;

/*	ALLOCATE GLOBAL STORAGE FOR EACH FIRST SURFACE */   
   sgbl = new struct surface_glbls[nsurf];
               
/* CREATE NSURF SURFACE OBJECTS FOR EACH GRID */
   for (j=0;j<ngrid;++j)
      grd[j].srf = new class surface[nsurf];

/*	ALLOCATE SURFACE OBJECTS */
   for(i=0;i<nsurf;++i) {
      fscanf(fp,"%*[^\n]%d\n",&surfid);
      printf("#SURFACE ID\n#%d\n",surfid);

      for(bnum=0;bnum<grd[0].nsbd;++bnum) 
         if (grd[0].sbdry[bnum].type == surfid) break;
      assert(bnum != grd[0].nsbd);

/*		ALLOCATE FINE SURFACE STORAGE */ 
      grd[0].srf[i].alloc(grd[0].maxsbel, lg2pmax, 0, 1, &sgbl[i]);
      for(j = lg2pmax -1; j >= 0; --j) {
         grd[0].srf[i].alloc(grd[0].maxsbel, j, 1, 1, &sgbl[i]);
      }
      grd[0].sbdry[i].misc = static_cast<void *>(&grd[0].srf[i]);

/*		ALLOCATE COARSE SURFACE STORAGE */
      for(j=1;j<ngrid;++j) {
         grd[j].srf[i].alloc(grd[j].maxsbel, 0, 1, 0, &sgbl[i]);
         grd[j].sbdry[i].misc = static_cast<void *>(&grd[j].srf[i]);
      }
      
/*		READ SURFACE PHYSICS INFO */
      fscanf(fp,"%*[^\n]%lf %lf %lf\n",&sgbl[i].sigma,&sgbl[i].rho2,&sgbl[i].mu2);
      printf("#SIGMA\t\tRHO2\t\tMU2\n#%f\t\t%f\t\t%f\n",sgbl[i].sigma,sgbl[i].rho2,sgbl[i].mu2);
   }

//	scanf("%*[^\n]%d %d\n",&rf, &rfnum);
//	printf("#READ FILE?\tFILE #\n#%1d\t\t%d\n",rf,rfnum);

   grd[0].loadbasis(hpbase[lg2pmax]);
   grd[0].curvinit();
   grd[0].tobasis(&f1);
   grd[0].surfvrttoug();

   return;
}

void block::tadvance() {
   int j;
   
   grd[0].tadvance();
   
   for(j=0;j<ngrid;++j)
      grd[j].setksprg1d();
   
}

void block::reconnect() {
   int i,j;

   for(i = 1; i< ngrid; ++i) {
      grd[i].coarsen(grd[i-1]);
      grd[i].setfine(grd[i-1]);
      grd[i-1].setcoarse(grd[i]);
   }

   return;
}


