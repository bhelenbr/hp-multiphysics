#include"blocks.h"
#include<cstring>
#include<utilities.h>
#include<stdio.h>

extern FLT f1(int n, FLT x, FLT y);

void blocks::init(char *file) {
   int i,j,k,p,match;
   FLT temp;
   FILE *fp;
   char blockname[100];
   
   fp = fopen(file,"r");
   if (fp == NULL) {
      printf("couldn't open file %s\n",file);
      exit(1);
   }

/*	READ POLYNOMIAL DEGREE */
   fscanf(fp,"%*[^\n]%d\n",&lg2pmax);
   printf("#LOG_2 PMAX\n#%d\n",lg2pmax);

/*	READ FLOW ITERATIVE INFORMATION */
/*	ITERATION CFL NUMBERS */
   fscanf(fp,"%*[^\n]");
   printf("#FLOWCFL\n#");
   for(i=0;i<=lg2pmax;++i) {
      fscanf(fp,"%lf",&hp_mgrid::cfl[i]);
      printf("%.2f\t",hp_mgrid::cfl[i]);
   }
   fscanf(fp,"%*[^\n]\n");
   printf("\n");

/*	LOAD FADD, ADIS, CHARACTERITIC FLAG */   
   fscanf(fp,"%*[^\n]%lf %lf %d\n",&hp_mgrid::fadd,&hp_mgrid::adis,&hp_mgrid::charyes);
   printf("#FADD\t\tADIS\t\tCHRCTR\n#%.2f\t\t%.2f\t\t%d\n",hp_mgrid::fadd,hp_mgrid::adis,hp_mgrid::charyes);
   
/* LOAD ADAPTATION INFORMATION */
   fscanf(fp,"%*[^\n]%d %lf %lf\n",&adapt,&hp_mgrid::trncerr,&hp_mgrid::tol);
   printf("#ADAPT\t\tTRNCERR\t\tTOLERANCE\n#%d\t\t%.2f\t\t%.2f\n",adapt,hp_mgrid::trncerr,hp_mgrid::tol);

 
 
/*	READ SURFACE ITERATIVE INFORMATION */
   fscanf(fp,"%*[^\n]");
   printf("#TANGENT/NORMAL SURFCFLS\n#");
   for(i=0;i<=lg2pmax;++i) {
      fscanf(fp,"%lf,%lf ",&surface::cfl[i][0],&surface::cfl[i][1]);
      printf("%.2f,%.2f\t",surface::cfl[i][0],surface::cfl[i][1]);
   }
   fscanf(fp,"%*[^\n]\n");
   printf("\n");
   
   fscanf(fp,"%*[^\n]%lf %lf\n",&surface::fadd[0],&surface::fadd[1]);
   printf("#TANGENT/NORMAL FADD\n#%0.2f\t%0.2f\n",surface::fadd[0],surface::fadd[1]); 
   
/*	READ MESH MOVEMENT ITERATIVE INFORMATION */
   fscanf(fp,"%*[^\n]%lf\n",&temp); //&r_mesh::cfl);
   printf("#MESH MOVEMENT CFL\n#%f\n",temp); //r_mesh::cfl);

/*	READ MESH MOVEMENT FADD */
   fscanf(fp,"%*[^\n]%lf\n",&temp); //&r_mesh::fadd);
   printf("#MESH MOVEMENT FADD\n#%f\n",temp); //r_mesh::fadd);

/*	READ IN PHYSICAL PARAMETERS */
   fscanf(fp,"%*[^\n]%lf %lf\n",&hp_mgrid::dti,&hp_mgrid::g);
   printf("#1/DT\t\tgravity\n#%0.2f\t%0.2f\n",hp_mgrid::dti,hp_mgrid::g);

/*	SET BACKWARDS DIFFERENCE CONSTANTS */
   for(i=0;i<NSTEP+1;++i)
      hp_mgrid::bd[i] = 0.0;
   hp_mgrid::bd[0] = hp_mgrid::dti;
   hp_mgrid::bd[1] = -hp_mgrid::dti;

/*	READ IN NUMBER OF TIME STEP/OUTPUT INTERVAL */
   fscanf(fp,"%*[^\n]%d %d\n",&ntstep,&out_intrvl);
   printf("#NTSTEP\t\tOUTPUT INTERVAL\n#%d\t\t%d\n",ntstep,out_intrvl);
   
/*	READ IN ITERATION PARAMETERS */
   fscanf(fp,"%*[^\n]%d %d %d\n",&mglvls,&vwcycle,&ncycle);
   printf("#MGLEVELS\t\tVWCYCLE\t\t# OF CYCLES\n#%d\t\t%d\t\t%d\n",mglvls,vwcycle,ncycle);
   
/*	READ IN NUMBER OF BLOCKS */
   fscanf(fp,"%*[^\n]%d\n",&nblocks);  
   printf("#NBLOCKS\n#%d\n",nblocks);
   

/*	BEGIN INITIALIZATION OF EACH BLOCK */   
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
   for (i=0;i<nblocks;++i) {
      fscanf(fp,"%s\n",blockname);
      printf("#\n#opening block %d: %s\n#\n",i,blockname);
      blk[i].initialize(blockname, mgrids, base, lg2pmax);
   }

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
   
   return;
}


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

FLT outertime = 0.0;

void blocks::tadvance() {
   int i;
   
   if (hp_mgrid::dti > 0.0) 
      hp_mgrid::time = hp_mgrid::time +1./hp_mgrid::dti;
      outertime = hp_mgrid::time;  // TEMPORARY
   
   for(i=0;i<nblocks;++i)
      blk[i].tadvance();
   
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

void blocks::output(int number, FILETYPE type=text) {
   int i;
   char outname[20], bname[20];

   number_str(outname, "data", number, 3);
   strcat(outname, ".");
   for(i=0;i<nblocks;++i) {
      number_str(outname, "data", i, 1);
      strcat(outname, ".");
      number_str(bname,outname,number,3);
      blk[i].grd[0].output(bname,type);
		if (adapt) blk[i].grd[0].out_mesh(bname,gambit);            
   }
   
   return;
}

void blocks::go() {
   int i,tstep, iter;
   
   for(tstep=0;tstep<ntstep;++tstep) {
      
      tadvance();
      
      for(iter=0;iter<ncycle;++iter) {
         cycle(vwcycle);
         printf("%d ",iter);
         print_maxres();
         printf("\n");
      }
         
      if (!(tstep%out_intrvl)) {
         output(tstep,tecplot);
      }
      
      if (adapt) {
         adaptation();
         blk[0].grd[0].out_mesh("test");
         exit(1);
         for(i=0;i<nblocks;++i)
            blk[i].reconnect();
      }
      
      if (tstep == 0) {
			hp_mgrid::bd[0] =  1.5*hp_mgrid::dti;
         hp_mgrid::bd[1] = -2.0*hp_mgrid::dti;
         hp_mgrid::bd[2] =  0.5*hp_mgrid::dti;
		}
//		else {
//			hp_mgrid::bd[0] = 11./6*hp_mgrid::dti;
//			hp_mgrid::bd[1] = -3.*hp_mgrid::dti;
//			hp_mgrid::bd[2] = 1.5*hp_mgrid::dti;
//			hp_mgrid::bd[3] = -1./3.*hp_mgrid::dti;
//		}
   }
   
   return;
}
   
void blocks::adaptation() {
   int i;
   
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].density1();
      
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].density2();
      
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].adapt(temp_hp,0.66);
      
   return;
}