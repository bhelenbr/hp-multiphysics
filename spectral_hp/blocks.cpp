#include"blocks.h"
#include<cstring>
#include<utilities.h>
#include<stdio.h>
#include"osdepend.h"

extern FLT f1(int n, FLT x, FLT y); //INITIALIZATION FUNCTIONS
extern FLT f2(int n, FLT x, FLT y);
extern int startup;  // USED IN MOVEPTTOBDRY TO SWITCH FROM INITIALIZATION TO ADAPTION

static int iter;

void blocks::init(char *file) {
   int i,j,p;
   FILE *fp;
   char blockname[100], outname[20], bname[20];
   
   fp = fopen(file,"r");
   if (fp == NULL) {
      printf("couldn't open file %s\n",file);
      exit(1);
   }

   /* READ POLYNOMIAL DEGREE */
   fscanf(fp,"%*[^\n]%d\n",&lg2pmax);
   printf("#LOG_2 PMAX\n#%d\n",lg2pmax);

   /* READ FLOW ITERATIVE INFORMATION */
   /* ITERATION CFL NUMBERS */
   fscanf(fp,"%*[^\n]");
   fscanf(fp,"%lf%lf%lf\n",&hp_mgrid::cfl[0],&hp_mgrid::cfl[1],&hp_mgrid::cfl[2]);
   printf("#FLOWCFL\n#%.2f\t%.2f\t%.2f\n",hp_mgrid::cfl[0],hp_mgrid::cfl[1],hp_mgrid::cfl[2]);

   /* LOAD FADD, ADIS, CHARACTERITIC FLAG */   
   fscanf(fp,"%*[^\n]%lf %lf %d\n",&hp_mgrid::fadd,&hp_mgrid::adis,&hp_mgrid::charyes);
   printf("#FADD\t\tADIS\t\tCHRCTR\n#%.2f\t\t%.2f\t\t%d\n",hp_mgrid::fadd,hp_mgrid::adis,hp_mgrid::charyes);
   
   /* LOAD ADAPTATION INFORMATION */
   fscanf(fp,"%*[^\n]%d %lf %lf %lf %lf\n",&adapt,&hp_mgrid::trncerr,&hp_mgrid::invbdryerr,&hp_mgrid::vlngth_tol,&hp_mgrid::adapt_tol);
   printf("#ADAPT\t\tTRNCERR\t\tBDRYERR\t\tVLNGTH_TOL\t\tADAPT_TOL\n#%d\t\t%.2e\t\t%.2e\t\t%.2f\t\t%.2f\n",adapt,hp_mgrid::trncerr,hp_mgrid::invbdryerr,hp_mgrid::vlngth_tol,hp_mgrid::adapt_tol);
 
   /* READ SURFACE ITERATIVE INFORMATION */
   fscanf(fp,"%*[^\n]");
   fscanf(fp,"%lf,%lf%lf,%lf%lf,%lf\n",&surface::cfl[0][0],&surface::cfl[0][1],&surface::cfl[1][0],&surface::cfl[1][1],&surface::cfl[2][0],&surface::cfl[2][1]);
   printf("#TANGENT/NORMAL SURFCFLS\n#%.2f,%.2f\t%.2f,%.2f\t%.2f,%.2f\n",surface::cfl[0][0],surface::cfl[0][1],surface::cfl[1][0],surface::cfl[1][1],surface::cfl[2][0],surface::cfl[2][1]);
      
   fscanf(fp,"%*[^\n]%lf %lf\n",&surface::fadd[0],&surface::fadd[1]);
   printf("#TANGENT/NORMAL FADD\n#%0.2f\t%0.2f\n",surface::fadd[0],surface::fadd[1]); 
   
   /* READ MESH MOVEMENT ITERATIVE INFORMATION */
   fscanf(fp,"%*[^\n]%lf\n",&r_mesh::vnn); //&r_mesh::cfl);
   printf("#MESH MOVEMENT CFL\n#%f\n",r_mesh::vnn); //r_mesh::cfl);

   /* READ MESH MOVEMENT FADD */
   fscanf(fp,"%*[^\n]%lf\n",&r_mesh::fadd); //&r_mesh::fadd);
   printf("#MESH MOVEMENT FADD\n#%f\n",r_mesh::fadd); //r_mesh::fadd);

   /* READ IN PHYSICAL PARAMETERS */
   fscanf(fp,"%*[^\n]%lf %lf\n",&hp_mgrid::dti,&hp_mgrid::g);
   printf("#1/DT\t\tgravity\n#%0.2f\t%0.2f\n",hp_mgrid::dti,hp_mgrid::g);

   /* SET BACKWARDS DIFFERENCE CONSTANTS */
   hp_mgrid::setbd(1);
   hp_mgrid::extrap = 0; /* DON'T EXTRAPOLATE ON FIRST CALL TO TADVANCE */

   /* READ IN NUMBER OF TIME STEP/OUTPUT INTERVAL */
   fscanf(fp,"%*[^\n]%d %d\n",&ntstep,&out_intrvl);
   printf("#NTSTEP\t\tOUTPUT INTERVAL\n#%d\t\t%d\n",ntstep,out_intrvl);
   
   /* READ IN ITERATION PARAMETERS */
   fscanf(fp,"%*[^\n]%d %d %d\n",&mglvls,&vwcycle,&ncycle);
   printf("#MGLEVELS\t\tVWCYCLE\t\t# OF CYCLES\n#%d\t\t%d\t\t%d\n",mglvls,vwcycle,ncycle);

   /* READ IN INITIALIZATION PARAMETER */
   fscanf(fp,"%*[^\n]%d\n",&readin);
   printf("#READ FILE #\n#%d\n",readin);
   
   /* READ IN NUMBER OF BLOCKS */
   fscanf(fp,"%*[^\n]%d\n",&nblocks);  
   printf("#NBLOCKS\n#%d\n",nblocks);
   
   /* BEGIN INITIALIZATION OF EACH BLOCK */   
   mgrids = MAX(mglvls-lg2pmax,1);   

   /* INITIALIZE BASIS FUNCTIONS */
   p = 1;
   for(i=0;i<lg2pmax;++i)
      p = p<<1;
   for(i=lg2pmax;i>=0;--i) {
#ifdef AXISYMMETRIC
      base[i].initialize(p,p+2);
#else
      base[i].initialize(p);
#endif
      p = p>>1;
   }
   
   blk = new class block[nblocks];
   /* OPEN EACH BLOCK FILE AND READ IN DATA */
   for (i=0;i<nblocks;++i) {
      fscanf(fp,"%s\n",blockname);
      printf("#\n#opening block %d: %s\n#\n",i,blockname);
      blk[i].initialize(blockname, mgrids, base, lg2pmax);
   }

   /* MATCH BOUNDARIES FOR EACH MGRID LEVEL */
   for(i=0;i<mgrids;++i) 
      findmatch(i);


   /* INITIALIZE SOLUTION FOR EACH BLOCK */
   if (readin > 0) {
      startup = 0;
      ntstep += readin;
      
      for(i=0;i<nblocks;++i) {
      
         /* INPUT MESH */
         number_str(bname,"mesh",readin,3);
         strcat(bname, ".");
         number_str(outname, bname, i, 1);
         printf("#Reading mesh: %s\n",outname);
         blk[i].grd[0].in_mesh(outname,grid);  
         blk[i].grd[0].spectral_hp::setbcinfo();

         /* FOR ADAPTIVE MESH */ 
         if (adapt) {
            number_str(bname,"vlgth",readin,3);
            strcat(bname, ".");
            number_str(outname, bname, i, 1);
            printf("#Reading vlength: %s\n",outname);
            blk[i].grd[0].inlength(outname);
         }
         /* READ SOLUTION */ 
         number_str(bname,"data",readin,3);
         strcat(bname, ".");
         number_str(outname, bname, i, 1); 
         printf("#Reading solution: %s\n",outname);
         blk[i].grd[0].spectral_hp::input(outname,text);

         /* INPUT UNSTEADY TIME HISTORY */
         number_str(outname,"rstrtdata",readin,3);
         strcat(outname, ".");
         for(j=0;j<MXSTEPM1;++j) {
            number_str(bname,outname,j,1);
            strcat(bname, ".");
            number_str(bname, bname, i, 1);
            printf("#Reading restart data: %s\n",bname);
            blk[i].grd[0].input(blk[i].gbl.ugbd[j],blk[i].gbl.vrtxbd[j],blk[i].gbl.binfobd[j],bname,text);
         }
         
         number_str(outname,"rstrtvrtx",readin,3);
         strcat(outname, ".");
         for(j=0;j<MXSTEPM1;++j) {
            number_str(bname,outname,j,1);
            strcat(bname, ".");
            number_str(bname, bname, i, 1);
            printf("#Reading restart mesh: %s\n",bname);
            blk[i].grd[0].in_mesh(blk[i].gbl.vrtxbd[j],bname,text);
         }
      }
      hp_mgrid::setbd(MXSTEP);
            
      /* REFIND BOUNDARIES ON FINE MESH */
      findmatch(0);

      /* DO ADAPTATION FOR NEXT TIME STEP */      
      if (adapt) {
         adaptation();

         /* REFIND BOUNDARIES FOR COARSE MESHES */
         /* JUST IN CASE BDRY ORDERING CHANGED */            
         for(i=1;i<mgrids;++i) 
            findmatch(i);
      }
      else {
         /* THIS MOVES THE COARSE MESH VERTICES TO NEW POSITIONS */
         for(i = 1; i < mgrids; ++i)
            for(j=0;j<nblocks;++j)
               blk[j].grd[i].r_mesh::mg_getfres();
      }
   }
   else {
      for(i=0;i<nblocks;++i) {
         blk[i].grd[0].curvinit();
         blk[i].grd[0].tobasis(&f1);
      }
      startup = 0;
   }

   for(i=0;i<nblocks;++i) {
      number_str(bname,"multigrid.",i,1);
      strcat(bname,".");
      blk[i].coarsenchk(bname);
   }
   
#ifdef PV3
   viz_init();
#endif

   return;
}

void blocks::findmatch(int grdlvl) {
   int j,k,match;
   
   for(j=0;j<nblocks;++j)
      blk[j].grd[grdlvl].zerobdryisfrst();
   
   for(j=0;j<nblocks;++j) {
      match = 0;
      for(k=0;k<nblocks;++k)
         match += blk[j].grd[grdlvl].findmatch(blk[k].grd[grdlvl]);
         
      if (match != blk[j].grd[grdlvl].alld_mp()) {
         printf("error in matching boundaries %d: %d %d\n",j,match,blk[j].grd[grdlvl].alld_mp());
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
      outertime = hp_mgrid::time;
   
   r_ksrc();

   for(i=0;i<nblocks;++i)
      blk[i].tadvance();
   
   return;
}

void blocks::nstage(int grdnum, int sm, int mgrid) {
   int i,stage,mode;
      
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

      /* CALCULATE RESIDUAL */   
      for(i=0;i<nblocks;++i)
         blk[i].grd[grdnum].rsdl(stage,mgrid);
                           
      /* INVERT MASS MATRIX (4 STEP PROCESS) */
      /* HAVE TO BE VERY CAREFUL WITH COMMUNICATION */
      /* USE MESSAGE BEFORE SENDING NEXT */
      /* VERTICES MUST BE 2 PART PROCESS BECAUSE OF CORNERS */
      for(i=0;i<nblocks;++i)
         blk[i].grd[grdnum].minvrt1();  //SEND Y
         
      for(i=0;i<nblocks;++i)
         blk[i].grd[grdnum].bdry_mp(); //RCV Y SEND X
         
      for(i=0;i<nblocks;++i)
         blk[i].grd[grdnum].minvrt2(); //RCV X SEND SIDE MODE 0 Y & X
                  
      for(mode=0;mode<sm-1;++mode) {
         for(i=0;i<nblocks;++i)
            blk[i].grd[grdnum].minvrt3(mode); // RCV SIDE MODE Y & X
            
         for(i=0;i<nblocks;++i)
            blk[i].grd[grdnum].minvrt3_mp(mode+1); //SEND SIDE Y & X
      }
      
      if (sm) {
         for(i=0;i<nblocks;++i)
            blk[i].grd[grdnum].minvrt4(); //RCV Y & X
      }
      
      for(i=0;i<nblocks;++i) 
         blk[i].grd[grdnum].nstage2(stage);
   }
   
//   if (mgrid == 0) 
//      for(i=0;i<nblocks;++i) 
//         blk[i].grd[grdnum].maxres();

   return;
}

void blocks::cycle(int vw, int lvl = 0) {
   int i,j;  // DON'T MAKE THESE STATIC SCREWS UP RECURSION
   int grid,bsnum;
#ifdef PV3
   float pvtime = 0.0;
#endif

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
      
#ifdef PV3
      if (lvl == 0) {
         pvtime = 1.0*iter;
         pV_UPDATE(&pvtime);
      }
#endif
      
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

   /* CALCULATE MAX RESIDUAL ON FINEST MESH */
   /* MOVE INTERFACE VERTICES */
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
   int i,j,tstep;
   char outname[20], bname[20];
   
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
      
      hp_mgrid::setbd(MIN(MXSTEP,tstep+2));
      
      blk[0].grd[0].drag(66564);
      if (!(tstep%out_intrvl)) {
         output(tstep+1,text);
         output(tstep+1,tecplot);
      }
       
#define NDEBUG
#ifdef DEBUG
      /*	THIS IS TO ALLOW MACHINE EXACT RESTARTS WHEN DEBUGGING */
      for(i=0;i<nblocks;++i) {
         readin = tstep +1;
         /* INPUT MESH */
         number_str(bname,"mesh",readin,3);
         strcat(bname, ".");
         number_str(outname, bname, i, 1);
         printf("#Reading mesh: %s\n",outname);
         blk[i].grd[0].in_mesh(outname,grid);  
         blk[i].grd[0].spectral_hp::setbcinfo();

         /* FOR ADAPTIVE MESH */ 
         if (adapt) {
            number_str(bname,"vlgth",readin,3);
            strcat(bname, ".");
            number_str(outname, bname, i, 1);
            printf("#Reading vlength: %s\n",outname);
            blk[i].grd[0].inlength(outname);
         }
         /* READ SOLUTION */ 
         number_str(bname,"data",readin,3);
         strcat(bname, ".");
         number_str(outname, bname, i, 1); 
         printf("#Reading solution: %s\n",outname);
         blk[i].grd[0].spectral_hp::input(outname,text);

         /* INPUT UNSTEADY TIME HISTORY */
         number_str(outname,"rstrtdata",readin,3);
         strcat(outname, ".");
         for(j=0;j<MXSTEPM1;++j) {
            number_str(bname,outname,j,1);
            strcat(bname, ".");
            number_str(bname, bname, i, 1);
            printf("#Reading restart data: %s\n",bname);
            blk[i].grd[0].input(blk[i].gbl.ugbd[j],blk[i].gbl.vrtxbd[j],blk[i].gbl.binfobd[j],bname,text);
         }
         
         number_str(outname,"rstrtvrtx",readin,3);
         strcat(outname, ".");
         for(j=0;j<MXSTEPM1;++j) {
            number_str(bname,outname,j,1);
            strcat(bname, ".");
            number_str(bname, bname, i, 1);
            printf("#Reading restart mesh: %s\n",bname);
            blk[i].grd[0].in_mesh(blk[i].gbl.vrtxbd[j],bname,text);
         }
      }
#endif

      
      if (adapt && tstep != ntstep-1)  adaptation();
   }
   
   return;
}
   
void blocks::adaptation() {
   int i;
   char adaptfile[20];
   
   /* THIS IS TO ENSURE COMMUNICATION BOUNDARES GET */
   /* COARSENED THE SAME WAY */
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].matchboundaries1();
      
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].matchboundaries2();

   /* SET-UP LENGTH FUNCTION */
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].length1();
      
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].length_mp();
      
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].length2();

   /*	ADAPT MESH */      
   for(i=0;i<nblocks;++i) {
      number_str(adaptfile,"adapt",i,1);
      blk[i].grd[0].adapt(block::temp_hp,adaptfile);
   }

   /* RECONNECT COARSE MESHES */
   for(i=0;i<nblocks;++i)
      blk[i].reconnect();

   return;
}

void blocks::output(int number, FILETYPE type=text) {
   int i,j;
   char outname[20], bname[20];

   for(i=0;i<nblocks;++i) {
      number_str(bname,"data",number,3);
      strcat(bname, ".");
      number_str(outname, bname, i, 1);

      /* OUTPUT SOLUTION */
      blk[i].grd[0].output(outname,type);
            
      /* OUTPUT MESH */
      number_str(bname,"mesh",number,3);
      strcat(bname, ".");
      number_str(outname, bname, i, 1);         
      blk[i].grd[0].mesh::setbcinfo();
      blk[i].grd[0].out_mesh(outname,grid);
      blk[i].grd[0].spectral_hp::setbcinfo();
      
      if (adapt) {
         number_str(bname,"vlgth",number,3);
         strcat(bname, ".");
         number_str(outname, bname, i, 1);
         blk[i].grd[0].outlength(outname,type);      
      }
      
      /* OUTPUT UNSTEADY TIME HISTORY */
      /*	FIRST INDEX IS HISTORY NUMBER */
      /* SECOND IS BLOCK NUMBER */
      number_str(outname,"rstrtdata",number,3);
      strcat(outname, ".");
      for(j=0;j<MXSTEPM1;++j) {
         number_str(bname,outname,j,1);
         strcat(bname, ".");
         number_str(bname, bname, i, 1);
         blk[i].grd[0].output(blk[i].gbl.ugbd[j],blk[i].gbl.vrtxbd[j],blk[i].gbl.binfobd[j],bname,type);
      }

      /* OUTPUT INTERPOLATED MESH TIME HISTORY */
      number_str(outname,"rstrtvrtx",number,3);
      strcat(outname, ".");
      for(j=0;j<MXSTEPM1;++j) {
         number_str(bname,outname,j,1);
         strcat(bname, ".");
         number_str(bname, bname, i, 1);
         blk[i].grd[0].out_mesh(blk[i].gbl.vrtxbd[j],bname,type);
      }
   }
   
   return;
}
