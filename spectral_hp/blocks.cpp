#include"blocks.h"
#include<cstring>
#include<utilities.h>
#include<stdio.h>

extern FLT f1(int n, FLT x, FLT y);

void blocks::init(int nb, int mg, int lg2p, char *filename) {
   int i,j,k,p,match;
   char fnmcat[80];
   char app[2];

   nblocks = nb;
   mglvls = mg;  // AT LEAST ONE ALWAYS
   lg2pmax = lg2p;
   mgrids = MAX(mglvls-lg2pmax,1);
    
/*	INITIALIZE BASIS FUNCTIONS */
   p = 1;
   for(i=0;i<lg2pmax;++i)
      p = p<<1;
   for(i=lg2pmax;i>=0;--i) {
      base[i].initialize(p);
      p = p>>1;
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
         blk[i].initialize(fnmcat, mgrids, base, lg2pmax);
      }
   }
   else
      blk[0].initialize(filename, mgrids, base, lg2pmax);

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
   
   hp_mgrid::setstatics(0.0,0.0,1.0);
   
/*	NEED TO INITIALIZE EACH WITH INITIALIZATION FUNCTION */
   for(i=0;i<nblocks;++i) {
      blk[i].grd[0].loadbasis(base[lg2pmax]);
      blk[i].grd[0].tobasis(&f1);
      blk[i].grd[0].surfvrttoug();
   }

   return;
}

void blocks::tadvance() {
   int i,j;
   
   time = time +dt;
   
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].setinflow();
      
   for(i=0;i<nblocks;++i)
      for(j=0;j<mgrids;++j)
         blk[i].grd[j].setksprg1d();
         
   return;
}

void blocks::nstage(int grdnum, int sm, int mgrid) {
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
         blk[i].grd[grdnum].rsdl(stage,mgrid);
         
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
   static int i,j,n;
   int grid,bsnum;
   
   grid = lvl -lg2pmax;
   bsnum =0;

/* ASSUMES WE ENTER WITH THE CORRECT BASIS LOADED */   
   if (lvl <= lg2pmax) {
      grid = 0;
      bsnum = lg2pmax-lvl;
   }
                  
   for (i=0;i<vw;++i) {

      nstage(grid,base[bsnum].sm,lvl);

/*		CALCULATE MAX RESIDUAL ON FINEST MESH */
      if (lvl == 0) {
         for(n=0;n<NV;++n)
            mxr[n] = 0.0;
      
         for (i=0;i<nblocks;++i) {
            blk[i].grd[0].maxres(mxr);
            blk[i].grd[0].surfugtovrt1();
         }
      }
      
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
   
   