#include"blocks.h"
#include<cstring>
#include<utilities.h>
#include<stdio.h>

extern FLT f1(int n, FLT x, FLT y); //INITIALIZATION FUNCTIONS
extern FLT f2(int n, FLT x, FLT y);
extern int startup;  // USED IN MOVEPTTOBDRY TO SWITCH FROM INITIALIZATION TO ADAPTION

void blocks::init(char *file) {
   int i,j,p;
   FILE *fp;
   char blockname[100], outname[20], bname[20];
   
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
   fscanf(fp,"%lf%lf%lf\n",&hp_mgrid::cfl[0],&hp_mgrid::cfl[1],&hp_mgrid::cfl[2]);
   printf("#FLOWCFL\n#%.2f\t%.2f\t%.2f\n",hp_mgrid::cfl[0],hp_mgrid::cfl[1],hp_mgrid::cfl[2]);

/*	LOAD FADD, ADIS, CHARACTERITIC FLAG */   
   fscanf(fp,"%*[^\n]%lf %lf %d\n",&hp_mgrid::fadd,&hp_mgrid::adis,&hp_mgrid::charyes);
   printf("#FADD\t\tADIS\t\tCHRCTR\n#%.2f\t\t%.2f\t\t%d\n",hp_mgrid::fadd,hp_mgrid::adis,hp_mgrid::charyes);
   
/* LOAD ADAPTATION INFORMATION */
   fscanf(fp,"%*[^\n]%d %lf %lf\n",&adapt,&hp_mgrid::trncerr,&hp_mgrid::tol);
   printf("#ADAPT\t\tTRNCERR\t\tTOLERANCE\n#%d\t\t%.2f\t\t%.2f\n",adapt,hp_mgrid::trncerr,hp_mgrid::tol);
 
/*	READ SURFACE ITERATIVE INFORMATION */
   fscanf(fp,"%*[^\n]");
   fscanf(fp,"%lf,%lf%lf,%lf%lf,%lf\n",&surface::cfl[0][0],&surface::cfl[0][1],&surface::cfl[1][0],&surface::cfl[1][1],&surface::cfl[2][0],&surface::cfl[2][1]);
   printf("#TANGENT/NORMAL SURFCFLS\n#%.2f,%.2f\t%.2f,%.2f\t%.2f,%.2f\n",surface::cfl[0][0],surface::cfl[0][1],surface::cfl[1][0],surface::cfl[1][1],surface::cfl[2][0],surface::cfl[2][1]);
      
   fscanf(fp,"%*[^\n]%lf %lf\n",&surface::fadd[0],&surface::fadd[1]);
   printf("#TANGENT/NORMAL FADD\n#%0.2f\t%0.2f\n",surface::fadd[0],surface::fadd[1]); 
   
/*	READ MESH MOVEMENT ITERATIVE INFORMATION */
   fscanf(fp,"%*[^\n]%lf\n",&r_mesh::vnn); //&r_mesh::cfl);
   printf("#MESH MOVEMENT CFL\n#%f\n",r_mesh::vnn); //r_mesh::cfl);

/*	READ MESH MOVEMENT FADD */
   fscanf(fp,"%*[^\n]%lf\n",&r_mesh::fadd); //&r_mesh::fadd);
   printf("#MESH MOVEMENT FADD\n#%f\n",r_mesh::fadd); //r_mesh::fadd);

/*	READ IN PHYSICAL PARAMETERS */
   fscanf(fp,"%*[^\n]%lf %lf\n",&hp_mgrid::dti,&hp_mgrid::g);
   printf("#1/DT\t\tgravity\n#%0.2f\t%0.2f\n",hp_mgrid::dti,hp_mgrid::g);

/*	SET BACKWARDS DIFFERENCE CONSTANTS */
   for(i=0;i<MXSTEP+1;++i)
      hp_mgrid::bd[i] = 0.0;
   hp_mgrid::bd[0] = hp_mgrid::dti;
   hp_mgrid::bd[1] = -hp_mgrid::dti;

/*	READ IN NUMBER OF TIME STEP/OUTPUT INTERVAL */
   fscanf(fp,"%*[^\n]%d %d\n",&ntstep,&out_intrvl);
   printf("#NTSTEP\t\tOUTPUT INTERVAL\n#%d\t\t%d\n",ntstep,out_intrvl);
   
/*	READ IN ITERATION PARAMETERS */
   fscanf(fp,"%*[^\n]%d %d %d\n",&mglvls,&vwcycle,&ncycle);
   printf("#MGLEVELS\t\tVWCYCLE\t\t# OF CYCLES\n#%d\t\t%d\t\t%d\n",mglvls,vwcycle,ncycle);

/*	READ IN INITIALIZATION PARAMETER */
   fscanf(fp,"%*[^\n]%d\n",&readin);
   printf("#READ FILE #\n#%d\n",readin);
   
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
/*	OPEN EACH BLOCK FILE AND READ IN DATA */
   for (i=0;i<nblocks;++i) {
      fscanf(fp,"%s\n",blockname);
      printf("#\n#opening block %d: %s\n#\n",i,blockname);
      blk[i].initialize(blockname, mgrids, base, lg2pmax);
   }

/*	MATCH BOUNDARIES FOR EACH MGRID LEVEL */
   for(i=0;i<mgrids;++i) 
      findmatch(i);


/*	INITIALIZE SOLUTION FOR EACH BLOCK */
   if (readin >= 0) {
      startup = 0;
      ntstep += readin +1;
      
      for(i=0;i<nblocks;++i) {
      
/*			INPUT MESH */
         number_str(bname,"mesh",readin,3);
         strcat(bname, ".");
         number_str(outname, bname, i, 1);
         printf("#Reading mesh: %s\n",outname);
         blk[i].grd[0].in_mesh(outname,easymesh);  
         blk[i].grd[0].spectral_hp::setbcinfo();

/*			FOR ADAPTIVE MESH */         
         if (adapt) {
            number_str(bname,"vlgth",readin,3);
            strcat(bname, ".");
            number_str(outname, bname, i, 1);
            printf("#Reading vlength: %s\n",outname);
            blk[i].grd[0].inlength(outname);
         }

/*			READ SOLUTION */ 
         number_str(bname,"data",readin,3);
         strcat(bname, ".");
         number_str(outname, bname, i, 1); 
         printf("#Reading solution: %s\n",outname);
         blk[i].grd[0].spectral_hp::input(outname,text);

/*			INPUT UNSTEADY TIME HISTORY */
         number_str(outname,"rstrtdata",readin,3);
         strcat(outname, ".");
         number_str(bname, outname, i, 1);
         for(j=0;j<MXSTEP-1;++j) {
            number_str(outname,bname,j,1);
            printf("#Reading restart data: %s\n",outname);
            blk[i].grd[0].input(blk[i].gbl.ugbd[j],blk[i].gbl.vrtxbd[j],blk[i].gbl.binfobd[j],outname,text);
         }
   
         number_str(outname,"rstrtvrtx",readin,3);
         strcat(outname, ".");
         number_str(bname, outname, i, 1);
         for(j=0;j<MXSTEP-1;++j) {
            number_str(outname,bname,j,1);
            printf("#Reading restart mesh: %s\n",outname);
            blk[i].grd[0].in_mesh(blk[i].gbl.vrtxbd[j],outname,text);
         }
      }
      
/*		REFIND BOUNDARIES ON FINE MESH */
      findmatch(0);

/*		DO ADAPTATION FOR NEXT TIME STEP */      
      if (adapt) {
         adaptation();

/*			CREATE COARSE MESHES */
         for(i=0;i<nblocks;++i)
            blk[i].reconnect();

/*			REFIND BOUNDARIES FOR COARSE MESHES */            
         for(i=1;i<mgrids;++i) 
            findmatch(i);
      }
   }
   else {
      for(i=0;i<nblocks;++i) {
         blk[i].grd[0].curvinit();
         blk[i].grd[0].tobasis(&f1);
      }
      startup = 0;
   }   

   return;
}

void blocks::findmatch(int i) {
   int j,k,match;
   
   for(j=0;j<nblocks;++j) {
      match = 0;
      for(k=0;k<nblocks;++k)
         match += blk[j].grd[i].findmatch(blk[k].grd[i]);
         
      if (match != blk[j].grd[i].alld_mp()) {
         printf("error in matching boundaries %d: %d %d\n",j,match,blk[j].grd[i].alld_mp());
         exit(1);
      }
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
      
   r_ksrc();
   
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
      blk[i].grd[grdnum].tstep_mp();

   for(i=0;i<nblocks;++i)
      blk[i].grd[grdnum].tstep2();

   for(i=0;i<nblocks;++i)
      blk[i].grd[grdnum].nstage1();
      
   for(stage=0;stage<NSTAGE;++stage) {

/*		CALCULATE RESIDUAL */   
      for(i=0;i<nblocks;++i)
         blk[i].grd[grdnum].rsdl(stage,mgrid);
                           
/*		INVERT MASS MATRIX (4 STEP PROCESS) */
/*		HAVE TO BE VERY CAREFUL WITH COMMUNICATION */
/* 	USE MESSAGE BEFORE SENDING NEXT */
/*		VERTICES MUST BE 2 PART PROCESS BECAUSE OF CORNERS */
      for(i=0;i<nblocks;++i)
         blk[i].grd[grdnum].minvrt1();  //SEND Y
         
      for(i=0;i<nblocks;++i)
         blk[i].grd[grdnum].bdry_mp(); //RCV Y SEND X
         
      for(i=0;i<nblocks;++i)
         blk[i].grd[grdnum].minvrt2(); //RCV X SEND Y & X
                  
      for(mode=0;mode<sm-1;++mode) {
         for(i=0;i<nblocks;++i)
            blk[i].grd[grdnum].minvrt3(mode); // RCV Y & X
            
         for(i=0;i<nblocks;++i)
            blk[i].grd[grdnum].minvrt3_mp(mode+1); //SEND Y & X
      }
      
      if (sm) {
         for(i=0;i<nblocks;++i)
            blk[i].grd[grdnum].minvrt4(); //RCV Y & X
      }
      
      for(i=0;i<nblocks;++i) 
         blk[i].grd[grdnum].nstage2(stage);
   }

   return;
}

void blocks::cycle(int vw, int lvl = 0) {
   int i,j;  // DON'T MAKE THESE STATIC SCREWS UP RECURSION
   int grid,bsnum;

/* ASSUMES WE ENTER WITH THE CORRECT BASIS LOADED */   
   if (lvl <= lg2pmax) {
      grid = 0;
      bsnum = lg2pmax-lvl;
   }
   else {
      grid = lvl -lg2pmax;
      bsnum =0;
   }
                  
   for (i=0;i<vw;++i) {

      nstage(grid,base[bsnum].sm,lvl);
      
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

void blocks::endcycle() {
   int i;
/*	CALCULATE MAX RESIDUAL ON FINEST MESH */
/*	MOVE INTERFACE VERTICES */

   for(i=0;i<nblocks;++i)
      blk[i].grd[0].maxres();
      
   for(i=0;i<nblocks;++i) {
      blk[i].grd[0].surfmaxres();
      blk[i].grd[0].surfugtovrt1();
   }
   
   for (i=0;i<nblocks;++i)
      blk[i].grd[0].surfugtovrt2();
}

void blocks::go() {
   int i,tstep, iter;
   
   for(tstep=readin;tstep<ntstep;++tstep) {
      
      tadvance();
      
      printf("#\n#TIMESTEP NUMBER %d\n",tstep+1);
      printf("ZONE\n");
      
      for(iter=0;iter<ncycle;++iter) {
         cycle(vwcycle);
         printf("%d ",iter);
         endcycle();
         r_cycle(vwcycle);
         r_maxres();
         printf("\n");
      }

//      if (tstep == 0 && hp_mgrid::dti > 0.0) {
//         hp_mgrid::nstep = 2;
//			hp_mgrid::bd[0] =  1.5*hp_mgrid::dti;
//         hp_mgrid::bd[1] = -2.0*hp_mgrid::dti;
//         hp_mgrid::bd[2] =  0.5*hp_mgrid::dti;
//		}
//		else {
//			hp_mgrid::nstep = 3;
//			hp_mgrid::bd[0] = 11./6*hp_mgrid::dti;
//			hp_mgrid::bd[1] = -3.*hp_mgrid::dti;
//			hp_mgrid::bd[2] = 1.5*hp_mgrid::dti;
//			hp_mgrid::bd[3] = -1./3.*hp_mgrid::dti;
//		}

      if (!(tstep%out_intrvl)) {
         output(tstep+1,text);
         output(tstep+1,tecplot);
      }
      
      if (adapt && tstep != ntstep-1) {
         adaptation();
         for(i=0;i<nblocks;++i)
            blk[i].reconnect();
      }
   }
   
   return;
}
   
void blocks::adaptation() {
   int i;
   
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].length1();
      
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].length_mp();
      
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].length2();
      
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].adapt(temp_hp,0.66);


   return;
}

void blocks::output(int number, FILETYPE type=text) {
   int i,j;
   char outname[20], bname[20];

   for(i=0;i<nblocks;++i) {
      number_str(bname,"data",number,3);
      strcat(bname, ".");
      number_str(outname, bname, i, 1);

/*		OUTPUT SOLUTION */
      blk[i].grd[0].output(outname,type);
      
/*		OUTPUT MESH */
      number_str(bname,"mesh",number,3);
      strcat(bname, ".");
      number_str(outname, bname, i, 1);         
      blk[i].grd[0].mesh::setbcinfo();
      blk[i].grd[0].out_mesh(outname,easymesh);  // THIS MUST BE EASYMESH OTHERWISE SIDE ORDERING WILL CHANGE
      blk[i].grd[0].spectral_hp::setbcinfo();
      
      if (adapt) {
         number_str(bname,"vlgth",number,3);
         strcat(bname, ".");
         number_str(outname, bname, i, 1);
         blk[i].grd[0].outlength(outname,text);      
      }
      
/*		OUTPUT UNSTEADY TIME HISTORY */
      number_str(outname,"rstrtdata",number,3);
      strcat(outname, ".");
      number_str(bname, outname, i, 1);
      for(j=0;j<MXSTEP-1;++j) {
         number_str(outname,bname,j,1);
         blk[i].grd[0].output(blk[i].gbl.ugbd[j],blk[i].gbl.vrtxbd[j],blk[i].gbl.binfobd[j],outname,type);
      }

      number_str(outname,"rstrtvrtx",number,3);
      strcat(outname, ".");
      number_str(bname, outname, i, 1);
      for(j=0;j<MXSTEP-1;++j) {
         number_str(outname,bname,j,1);
         blk[i].grd[0].out_mesh(blk[i].gbl.vrtxbd[j],outname,text);
      }
   }
   
   return;
}