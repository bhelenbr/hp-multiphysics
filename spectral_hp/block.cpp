#include"blocks.h"
#include"utilities.h"

void block::initialize(char *inputfile, int grds, class hpbasis *bin, int lg2p) {
   FILE *fp;
   char grd_nm[200];
   int fmt;
   int nsurf, surfid[MXSB];
   
   fp = fopen(inputfile,"r");
   
   fscanf(fp,"%*[^\n]%s %lf %d\n",grd_nm,&grwfac,&fmt);
	printf("#NAME\t\tFORMAT\t\tGROWTH\n#%s\t\t%d\t\t%.1f\n",grd_nm,&fmt,grwfac);

/*	LOAD MESH & CREATE COARSER MESHES */   
   ngrid = grds;
   grd = new class hp_mgrid[ngrid];
   grd[0].in_mesh(grd_nm,static_cast<FILETYPE> fmt,grwfac);
   for(i = 1; i< ngrid; ++i) {
      grd[i].coarsen(grd[i-1]);
/*    grd[i].smooth_cofa(2); */
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

/*	INITIALIZE STORAGE FOR MULTIGRID SOLUTION */ 
	grd[0].allocate(0,&gbl);  // 0 DENOTES FINEST LEVEL ALLOCATES GLOBAL STORAGE
   for(i = lg2pmax -1; i >= 0; --i) {
      grd[0].loadbasis(hpbase[i]);
      grd[0].allocate(1,&gbl); // 1 DENOTES MGRID LEVEL
   }
   grd[0].loadbasis(hpbase[lg2pmax]);
   for(i=1;i<ngrid;++i)
      grd[i].allocate(1,&gbl);

/*	LOAD ITERATIVE CONSTANTS FOR BLOCK */
/*	ITERATION CFL NUMBERS */
   fscanf(fp,"%*[^\n]");
   printf("#FLOWCFL\n");
   for(i=0;i<=lg2pmax;++i) {
      fscanf(fp,"%lf",&gbl.flowcfl[i]);
      printf("%.2f\t",gbl.flowcfl[i]);
   }
   printf("\n");
   
   fscanf(fp,"%*[^\n]%lf %lf %lf",&gbl.fadd,&gbl.adis,&gbl.charyes);
   printf("FADD\t\tADIS\t\tCHRCTR\n%.2f\t\t%.2f\t\t%d\n,",gbl.fadd,gbl.adis,gbl.charyes);

   fscanf(fp,"%*[^\n]%lf %lf",&gbl.rho,&gbl.mu);
   printf("RHO\t\tMU\n%.3f\t\t%.3f\n,",gbl.rho,gbl.mu);

   gbl.rhoi = 1.0/rho;
   gbl.nu = mu/rho;
/*	GBLS FOR ALL BLOCKS - dt, time, g, f_init, *hpbasis, ngrids, lg2pmax */

/*	LOAD SURFACE INFORMATION */
   nsurf = 0;
   for(i=0;i<grd[0].nsbd;++i)
      if (grd[0].sbdry[i].type&(FSRF_MASK +IFCE_MASK))
         surfids[nsurf++] = grd[0].sbdry[i].type;
         
   sgbl = new struct surface_gbls[nsurf];
   
   for(i=0;i<nsurf;++i) {
      fscanf(fp,"%*[^\n],%d %d",idnum);
      for(j=0;j<nsurf;++j)
         if(surfids[j] == idnum) break;
      assert(j != nsurf);
      
      
/*			LOAD ITERATIVE INFORMATION */
         fscanf(fp,"%*[^\n]");
         printf("#NORMAL SURFCFLS\n");
         for(i=0;i<=lg2pmax;++i) {
            fscanf(fp,"%lf",&gbl.flowcfl[i]);
            printf("%.2f\t",gbl.flowcfl[i]);
         }
         printf("\n");



	printf("#FLOWCFL\tSURFTAN\t\tSURFNORM\tMESHCFL\n");
	printf("#%.1lf,%.1lf,%.1lf\t%.1lf,%.1lf,%.1lf\t%.1lf,%.1lf,%.1lf\t%.1lf\n",   

   
	scanf("%*[^\n]%d %d\n",&rf, &rfnum);
	printf("#READ FILE?\tFILE #\n#%1d\t\t%d\n",rf,rfnum);		
	scanf("%*[^\n]%lf %d %d\n",&dt,&nstep,&skip);
	printf("#1/DT\t\tNSTEP\t\tSKIP\n#%.6le\t%04d\t\t%03d\n",dt,nstep,skip);
	scanf("%*[^\n]%d %d\n",&i_cyc,&iskip);

	mg_flowcfl[0],mg_flowcfl[1],mg_flowcfl[2],
	mg_surfcfl[0][0],mg_surfcfl[1][0],mg_surfcfl[2][0],
	mg_surfcfl[0][1],mg_surfcfl[1][1],mg_surfcfl[2][1],meshcfl);
	scanf("%*[^\n] %lf %lf %lf %lf\n",&ffadd,&sfadd0,&sfadd1,&mfadd);
	printf("#FLOFADD\tSURFFADD[0]\tSURFFADD[1]\tMESHFADD\n");
	printf("#%.2lf\t\t%.2lf\t\t%.2lf\t\t%.2lf\n",ffadd,sfadd0,sfadd1,mfadd);	
	scanf("%*[^\n] %d %lf %lf\n",&mgcyc,&mvfsm,&updsm);
	printf("#MGCYC\t\tMVFSM\t\tUPDSM\n#%1d\t\t%.2lf\t\t%.2lf\n"
	,mgcyc,mvfsm,updsm);		
	scanf("%*[^\n]%lf %lf %d\n",&u2max,&adis,&charyes);
	printf("#U2MAX\t\tADIS\t\tCHAR\n#%.2le\t%.2le\t%1d\n",u2max,adis,charyes);
	scanf("%*[^\n]%lf %lf %lf %lf\n",&rho[0], &rho[1], &mu[0], &mu[1]);
	printf("#RHO1\t\tRHO2\t\tMU1\t\tMU2\n#%.6le\t%.6le\t%.6le\t%.6le\n"
	,rho[0],rho[1],mu[0],mu[1]);
	scanf("%*[^\n]%lf %lf\n",&g,&sigma);
	printf("#G\t\tSIGMA\n#%.6le\t%.6le\n",g,sigma);
	ierr = scanf("%*[^\n]%d\n",&ncnstnt);
	if (ierr == 1 && ncnstnt > 0) {
		printf("#NCONSTANTS\n#%d\n",ncnstnt);
		vect_alloc(cnstnt,ncnstnt,double);
		for(i=0;i<ncnstnt;++i) {
			ierr = scanf("%[^\n]%lf\n",title,cnstnt+i);
			printf("#%s\n#%lf\n",title,cnstnt[i]);
			if (ierr != 2) printf("trouble with constants\n");
		}
	}
   return;
}


void block::reconnect() {
   int i;
   
   for(i = 1; i< ngrid; ++i) {
      grd[i].coarsen(grd[i-1]);
/*    grd[i].smooth_cofa(2); */
      grd[i].setfine(grd[i-1]);
      grd[i-1].setcoarse(grd[i]);
   }

   return;
}


