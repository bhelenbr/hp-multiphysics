#include"blocks.h"
#include<cstring>
#include<utilities.h>

extern FLT f1(int n, FLT x, FLT y);

void blocks::init(int nb, int mg, int lg2p, char *filename, FILETYPE filetype = easymesh, FLT grwfac = 1) {
   int i,j,k,p,match;
   FLT cfl[MXLG2P];
   char fnmcat[80];
   char app[2];

   nblocks = nb;
   mglvls = mg;  // AT LEAST ONE ALWAYS
   lg2pmax = lg2p;
   mgrids = MIN(mglvls-lg2pmax,1);
 
/*	INITIALIZE BASIS FUNCTIONS */
   p = 1;
   for(i=0;i<=lg2pmax;++i) {
      base[i].initialize(p);
      p *= 2;
   }
   
   blk = new class block[nblocks];

/*	ASSUME FOR NOW MESHES ARE LABELED a,b,c... */
/*	I HAVEN'T FIGURED OUT HOW THIS IS GOING TO WORK IN THE TOTALLY GENERAL CASE */
   if (nblocks > 1) {
      for (i=0;i<nblocks;++i) {
         strcpy(fnmcat,filename);
         app[0] = 'a'+i;
         app[1] = '\0';
         strcat(fnmcat,app);
         blk[i].meshinit(mgrids,fnmcat,filetype,grwfac);
      }
   }
   else
      blk[0].meshinit(mgrids,filename,filetype,grwfac);

/*	MATCH BOUNDARIES */
   for(i=0;i<mgrids;++i) {
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

/*	WHERE IS THIS GOING TO BE STORED? */   
/*	LOAD PHYSICAL CONSTANTS & ITERATIVE THINGS FOR EACH BLOCK */
   blk[0].setphysics(1.0, 1.0, 0.0, &f1);
   cfl[0] = 2.5;
   cfl[1] = 1.25;
   cfl[2] = 0.75;
   blk[0].setiter(0.75,cfl,1.0,0);
   
   for(i=0;i<nblocks;++i)
       blk[i].hpinit(base, lg2pmax);

/*	NEED TO INITIALIZE EACH WITH INITIALIZATION FUNCTION */
   for(i=0;i<nblocks;++i) {
      blk[i].grd[0].loadbasis(base[lg2pmax]);
      blk[i].grd[0].tobasis(&f1);
   }

   return;
}

void blocks::nstage(int grdnum, int sm) {
   static int i,stage,mode;
      
/*****************************************/
/* NSTAGE UPDATE OF FLOW VARIABLES    ****/
/*****************************************/
   for(i=0;i<nblocks;++i)
      blk[i].grd[grdnum].tstep1();

   for(i=0;i<nblocks;++i)
      blk[i].grd[grdnum].tstep2();

   for(i=0;i<nblocks;++i)
      blk[i].grd[grdnum].nstage1();
      
   for(stage=0;stage<NSTAGE;++stage) {

/*		CALCULATE RESIDUAL */   
      for(i=0;i<nblocks;++i)
         blk[i].grd[grdnum].rsdl(stage,grdnum);
         
/*		INVERT MASS MATRIX (4 STEP PROCESS) */
      for(i=0;i<nblocks;++i)
         blk[i].grd[grdnum].minvrt1();
         
      for(i=0;i<nblocks;++i)
         blk[i].grd[grdnum].minvrt2();
      
      for(mode=0;mode<sm-1;++mode)
         for(i=0;i<nblocks;++i)
            blk[i].grd[grdnum].minvrt3(mode);
      
      if (sm) {
         for(i=0;i<nblocks;++i)
            blk[i].grd[grdnum].minvrt4();
      }

      
      for(i=0;i<nblocks;++i) 
         blk[i].grd[grdnum].nstage2(stage);
   }

   return;
}

void blocks::cycle(int vw, int lvl = 0) {
   static int i,j;
   int grid,bsnum;
   
   grid = lvl -lg2pmax;
   bsnum =0;

/* ASSUMES WE ENTER WITH THE CORRECT BASIS LOADED */   
   if (lvl <= lg2pmax) {
      grid = 0;
      bsnum = lg2pmax-lvl;
   }
                  
   for (i=0;i<vw;++i) {

      nstage(lvl,base[bsnum].sm);

      if (lvl == mglvls-1) return;
      
      for(j=0;j<nblocks;++j)
         blk[j].grd[grid].rsdl(NSTAGE,lvl);
      
      if (bsnum == 0) {
         for(j=0;j<nblocks;++j)
            blk[j].grd[grid+1].getfres();
      }
      else {
         for(j=0;j<nblocks;++j)
            blk[j].grd[0].loadbasis(base[bsnum -1]);
            
         for(j=0;j<nblocks;++j)
            blk[j].grd[0].getfres();
      }
      
      cycle(vw, lvl+1);
      
      if (bsnum)
         for(j=0;j<nblocks;++j)
            blk[j].grd[0].loadbasis(base[bsnum]);

      for(j=0;j<nblocks;++j)
         blk[j].grd[grid].getcchng();
   }

   return;
}

void blocks::output(char *filename, FILETYPE filetype = text) {
   int i;   
   char fnmcat[80];
   char app[2];

/*	ASSUME FOR NOW MESHES ARE LABELED a,b,c... */
/*	I HAVEN'T FIGURED OUT HOW THIS IS GOING TO WORK IN THE TOTALLY GENERAL CASE */
   if (nblocks > 1) {
      for (i=0;i<nblocks;++i) {
         strcpy(fnmcat,filename);
         app[0] = 'a'+i;
         app[1] = '\0';
         strcat(fnmcat,app);
         blk[i].grd[0].output(fnmcat,filetype);
      }
   }
   else {
      blk[0].grd[0].output(filename,filetype);
   }
   
   return;
}
   