#include"blocks.h"
#include<cstring>
#include<utilities.h>
#include<stdio.h>
#include<time.h>
#include"pV3.h"

extern FLT f1(int n, FLT x, FLT y); //INITIALIZATION FUNCTIONS
extern FLT f2(int n, FLT x, FLT y);
extern int startup;  // USED IN MOVEPTTOBDRY TO SWITCH FROM INITIALIZATION TO ADAPTION

#ifdef LAYER
extern FLT mux[LAYER];
extern FLT rhox[LAYER];
#endif
static int iter;
extern FLT amp,lam,theta;

#ifdef DROP
extern FLT dydt;
static FLT factor;
#endif

void blocks::init(char *file, int start_sim) {
   int i,j,p;
   FILE *fp;
   char blockname[100],outname[20], bname[20];
   
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
#ifndef DROP
   fscanf(fp,"%*[^\n]%lf %lf\n",&hp_mgrid::dti,&hp_mgrid::g);
   printf("#1/DT\t\tgravity\n#%8.4e\t%8.4e\n",hp_mgrid::dti,hp_mgrid::g);
#else
   fscanf(fp,"%*[^\n]%lf %lf %lf\n",&hp_mgrid::dti,&hp_mgrid::g,&factor);
   printf("#1/DT\t\tgravity\n#%8.4e\t%8.4e\t%8.4e\n",hp_mgrid::dti,hp_mgrid::g,factor);
   factor = pow(10.0,factor);
#endif

   /* SET BACKWARDS DIFFERENCE CONSTANTS */
   hp_mgrid::setbd(1);
   hp_mgrid::extrap = 0; /* DON'T EXTRAPOLATE ON FIRST CALL TO TADVANCE */

   /* READ IN NUMBER OF TIME STEP/OUTPUT INTERVAL */
   fscanf(fp,"%*[^\n]%d %d %d\n",&ntstep,&out_intrvl,&rstrt_intrvl);
   printf("#NTSTEP\t\tOUTPUT INTERVAL\t\tRSTRT INTERVAL\n#%d\t\t%d\t\t%d\n",ntstep,out_intrvl,rstrt_intrvl);
   
   /* READ IN ITERATION PARAMETERS */
   fscanf(fp,"%*[^\n]%d %d %d\n",&mglvls,&vwcycle,&ncycle);
   printf("#MGLEVELS\t\tVWCYCLE\t\t# OF CYCLES\n#%d\t\t%d\t\t%d\n",mglvls,vwcycle,ncycle);

   /* READ IN ITERATION PARAMETERS */
   fscanf(fp,"%*[^\n]%d %d\n",&nup,&ndown);
   printf("#NSWEEPUP\t\tNSWEEPDN\n#%d\t\t%d\n",nup,ndown);

   /* READ IN INITIALIZATION PARAMETER */
   fscanf(fp,"%*[^\n]%d\n",&readin);
   printf("#READ FILE #\n#%d\n",readin);

   /* READ IN AMPLITUDE / WAVELENGTH / THETA */
   fscanf(fp,"%*[^\n]%lf%lf%lf\n",&amp,&lam,&theta);
   printf("#AMP WAVELENGTH THETA #\n#%e\t%e\t%e\n",amp,lam,theta);
   theta *= M_PI/180.0;
#ifdef DROP
   dydt = lam;
#endif
   
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
   
#ifdef LAYER
   for(i=0;i<nblocks;++i) {
      rhox[i] = blk[i].gbl.rho;
      mux[i] = blk[i].gbl.mu;
   }
#endif

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
         for(j=0;j<TMADAPT;++j) {
            number_str(bname,outname,j,1);
            strcat(bname, ".");
            number_str(bname, bname, i, 1);
            printf("#Reading restart data: %s\n",bname);
            blk[i].grd[0].input(blk[i].gbl.ugbd[j],blk[i].gbl.vrtxbd[j],blk[i].gbl.binfobd[j],bname,text);
         }
         
         number_str(outname,"rstrtvrtx",readin,3);
         strcat(outname, ".");
         for(j=0;j<TMADAPT;++j) {
            number_str(bname,outname,j,1);
            strcat(bname, ".");
            number_str(bname, bname, i, 1);
            printf("#Reading restart mesh: %s\n",bname);
            blk[i].grd[0].in_mesh(blk[i].gbl.vrtxbd[j],bname,text);
         }
      }
      hp_mgrid::setbd(MIN(TMSCHEME,readin));
            
      /* REFIND BOUNDARIES ON FINE MESH */
      findmatch(0);
   }
   else {
      for(i=0;i<nblocks;++i) {
          blk[i].grd[0].curvinit();
         // blk[i].grd[0].smooth_cofa(2);
         blk[i].grd[0].tobasis(&f1); 
      }
#ifdef DROP
      blk[1].grd[0].tobasis(&f2); 
#endif
      startup = 0;
   }
   
   /* DO ADAPTATION FOR CONTINUATION OF NEXT TIME STEP */
   if (start_sim)  {
      if (readin) {
         if (adapt) {
            /* DO ADAPTATION FOR NEXT TIME STEP */
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
      for(i=0;i<nblocks;++i) {
         number_str(bname,"multigrid.",i,1);
         strcat(bname,".");
         blk[i].coarsenchk(bname);         
      }
   }   

   
#ifdef PV3
   if (!start_sim)
      viz_init(0);
   else if(adapt)
      viz_init(-3);
   else
      viz_init(-2);
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

void blocks::tadvance(int stage) {
   int i,grid,lvl,bsnum;
   
#ifdef BACKDIFF
   r_ksrc();

   if (hp_mgrid::dti > 0.0) 
      hp_mgrid::time = hp_mgrid::time +1.0/hp_mgrid::dti;
      outertime = hp_mgrid::time;
      
   for(i=0;i<nblocks;++i) {
      for(lvl=0;lvl<mglvls;++lvl) {
         if (lvl <= lg2pmax) {
            grid = 0;
            bsnum = lg2pmax-lvl;
         }
         else {
            grid = lvl -lg2pmax;
            bsnum =0;
         }
         blk[i].grd[grid].loadbasis(base[bsnum]);
         blk[i].grd[grid].unsteady_sources(lvl);
      }
      blk[i].grd[0].loadbasis(base[lg2pmax]);
      blk[i].grd[0].shift();
   }
#else
   if (stage == 0) r_ksrc();
   
   if (hp_mgrid::dti > 0.0) 
      hp_mgrid::time += hp_mgrid::cdirk[stage]/hp_mgrid::dti;
   outertime = hp_mgrid::time;
   
   for(i=0;i<nblocks;++i)
      blk[i].tadvance(stage);
#endif

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
         blk[i].grd[grdnum].surfrsdl(mgrid);
      
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
         
      for(i=0;i<nblocks;++i) 
         blk[i].grd[grdnum].surfugtovrt1();
      
      for (i=0;i<nblocks;++i)
         blk[i].grd[grdnum].surfugtovrt2();
   }
   
   return;
}

void blocks::cycle(int vw, int lvl) {
   int j,vcount;  // DON'T MAKE THESE STATIC SCREWS UP RECURSION
   int grid,bsnum;
#if (defined(TWO_LEVEL) || defined(SOLVECOARSE))
   FLT mxr[NV], emax = 0.0, err;
   int crscntr = 0;
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


   for (vcount=0;vcount<vw;++vcount) {
   
      /* TEMPORARY TO TEST MULTIGRID */
//      if (lvl != 0) {
//           hp_mgrid::cfl[0] = 2.5;
//           surface::cfl[bsnum][1] = 0.0;
//        }
//        else {
//           hp_mgrid::cfl[0] = 0.5;
//           surface::cfl[bsnum][1] = 2.5;
//        }      
//      char fname[100];
//      static int count[3] = {0, 0, 0};
//      number_str(fname,"iterate",lvl,1);
//      strcat(fname,".");
//      number_str(fname,fname,count[lvl]++,3);
//      blk[0].grd[grid].output(fname,tecplot);

      for (int i=ndown;i--;)
         nstage(grid,base[bsnum].sm,lvl);

#ifdef TWO_LEVEL
      if (lvl == 1) {
         err = 0.0;
         for(int i=0;i<nblocks;++i) 
            err = MAX(err,blk[i].grd[grid].maxres(mxr));
                     
         emax = MAX(emax,err);
         if (err/emax > 3.0e-5 && err > 1.0e-11) --vcount;
         printf("# second level %e %e\n",emax,err);
         
//         char fname[100],tname[100];
//         static int count;
//         number_str(fname,"coarse",count++,3);
//         strcat(fname,".");
//         for(int i=0;i<nblocks;++i) {
//            number_str(tname,fname,i,1);
//            blk[i].grd[grid].output(tname,tecplot);
//         }
      }
#endif

      if (lvl == mglvls-1) {
#ifdef SOLVECOARSE
         /* DO A GOOD JOB ON COARSEST MESH */
         if (mglvls != 1) {

            err = 0.0;
            for(int i=0;i<nblocks;++i) 
               err = MAX(err,blk[i].grd[grid].maxres(mxr));
            
            emax = MAX(emax,err);
            --vcount;
            if (err/emax < 3.0e-1 || err < 1.0e-11 || crscntr > SOLVECOARSE) {
               printf("# Coarsest grid iterations %d\n",crscntr);
               vcount+=2;
            }
            ++crscntr;
         }
#endif
         continue;
      }
      
      for(j=0;j<nblocks;++j)
         blk[j].grd[grid].surfrsdl(lvl);
        
      for(j=0;j<nblocks;++j)
         blk[j].grd[grid].rsdl(NSTAGE,lvl);

      if (bsnum == 0) {
         for(j=0;j<nblocks;++j)
            blk[j].grd[grid+1].getfres();
            
         for(j=0;j<nblocks;++j) 
            blk[j].grd[grid+1].surfugtovrt1();
         
         for (j=0;j<nblocks;++j)
            blk[j].grd[grid+1].surfugtovrt2();
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
      
//      char fname[100];
//      static int count[3] = {0, 0, 0};
//      number_str(fname,"coarsechng",lvl,1);
//      number_str(fname,fname,count[lvl]++,1);
//      number_str(fname,fname,0,1);
//      blk[0].grd[grid].output(blk[0].grd[grid].gbl->res,blk[0].grd[grid].vrtx,blk[0].grd[grid].binfo,fname,tecplot);

         
      for(j=0;j<nblocks;++j) 
         blk[j].grd[grid].surfugtovrt1();
      
      for (j=0;j<nblocks;++j)
         blk[j].grd[grid].surfugtovrt2();
   }
   
   for (int i=nup;i--;)
      nstage(grid,base[bsnum].sm,lvl);

   return;
}

FLT blocks::printerror() {
   int i,n;
   FLT mxr[NV], biggest = 0.0;
   
   /* CALCULATE MAX RESIDUAL ON FINEST MESH */
   for(i=0;i<nblocks;++i) {
      blk[i].grd[0].maxres(mxr);
      for(n=0;n<NV;++n) {
         printf("%.3e  ",mxr[n]);
         biggest = MAX(biggest,mxr[n]);
      }
   }
      
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].surfmaxres();

   return(biggest);
}

#ifdef DROP
extern FLT dydt;
#endif

void blocks::go() {
   int tstep,total = 0;
   FLT error;
#ifdef PV3
   float pvtime = 0.0;
#endif
#ifdef PARAMETERLOOP
   FILE *fout;
   fout = fopen("dropdata.dat","r");
   if (fout) {
      fclose(fout);
      fout = fopen("dropdata.dat","a");
   }
   else {
      fout = fopen("dropdata.dat","w");
      fprintf(fout,"VARIABLES=\"File\" \"Dens\" \"mu\" \"Re\" \"We\" \"Oh\" \"Cd\" \"Eo\" \"E\" \"g\" \"dydt\" \"Scap\"\n");
   }
   FLT aratioold = 1.0, aratio = 1.0;
   FLT gold = hp_mgrid::g/factor;
   FLT ratio;
#endif
   
   for(tstep=readin;tstep<ntstep;++tstep) {
#ifdef BACKDIFF
      tadvance();
      {
#else 
      for(int s=0;s<TMSCHEME;++s) {
         tadvance(s);
#endif
         printf("#\n#TIMESTEP NUMBER %d\n",tstep+1);      
         for(iter=0;iter<ncycle;++iter) {
#ifdef DEFORM
            r_cycle(vwcycle);
#endif
            cycle(vwcycle);
            printf("%d %ld ",total++,clock());
            error = printerror();
#ifdef DEFORM
            r_printerror();
#endif
            printf("\n");
#ifdef PV3
            pvtime = 1.0*iter;
            pV_UPDATE(&pvtime);
#endif
            // output(iter+1000,tecplot);

         }
         // output(tstep*TMSCHEME +s,tecplot);
      }
      
      blk[0].grd[0].drag(1028);

      if (!(tstep%out_intrvl)) {
         output(tstep+1,tecplot,0); 
         if (!(tstep%(rstrt_intrvl*out_intrvl))) {
            output(tstep+1,text,1);
         }
         
#ifdef PARAMETERLOOP
         /* FIND MAXIMUM WIDTH */
         FLT maxw;
         int scap;
         int ind;
         scap = 0;
         maxw = 0.0;
         for(int temp = 0; temp < blk[1].grd[0].sbdry[0].num;++temp) {
            ind = blk[1].grd[0].sbdry[0].el[temp];
            ind = blk[1].grd[0].svrtx[ind][0];
            maxw = MAX(maxw,blk[1].grd[0].vrtx[ind][0]);
         }
         if (blk[1].grd[0].vrtx[ind][1] > blk[1].grd[0].vrtx[blk[1].grd[0].vbdry[1].el[0]][1]) scap = 1;      
         aratio = (blk[1].grd[0].vrtx[blk[1].grd[0].vbdry[1].el[0]][1]-blk[1].grd[0].vrtx[blk[1].grd[0].vbdry[0].el[0]][1])/(2.*maxw);
         fprintf(fout,"%d %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %d\n",tstep+1,
         blk[1].gbl.rho,blk[1].gbl.mu/blk[0].gbl.mu,(1.-dydt)/blk[0].gbl.mu,
         (1.-dydt)*(1-dydt)/blk[1].sgbl[0].sigma,blk[1].gbl.mu/sqrt(blk[1].gbl.rho*blk[1].sgbl[0].sigma),hp_mgrid::g*(blk[1].gbl.rho-blk[0].gbl.rho)/(12.*(1-dydt)*(1-dydt)),hp_mgrid::g*(blk[1].gbl.rho-blk[0].gbl.rho)/blk[1].sgbl[0].sigma,
         aratio,hp_mgrid::g,dydt,scap);

         ratio = fabs((aratioold-aratio)/aratio);
         if (ratio > 0.05) {
            factor = pow(factor,0.1);
            printf("#factor change %e\n",factor);
	 }
         aratioold = aratio;
         gold = hp_mgrid::g;
         hp_mgrid::g=hp_mgrid::g*factor;
#endif
      }
       
#ifdef DEBUG
      int i,j;
      char outname[20], bname[20];
      
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
         for(j=0;j<TMADAPT;++j) {
            number_str(bname,outname,j,1);
            strcat(bname, ".");
            number_str(bname, bname, i, 1);
            printf("#Reading restart data: %s\n",bname);
            blk[i].grd[0].input(blk[i].gbl.ugbd[j],blk[i].gbl.vrtxbd[j],blk[i].gbl.binfobd[j],bname,text);
         }
         
         number_str(outname,"rstrtvrtx",readin,3);
         strcat(outname, ".");
         for(j=0;j<TMADAPT;++j) {
            number_str(bname,outname,j,1);
            strcat(bname, ".");
            number_str(bname, bname, i, 1);
            printf("#Reading restart mesh: %s\n",bname);
            blk[i].grd[0].in_mesh(blk[i].gbl.vrtxbd[j],bname,text);
         }
      }
#endif
      
      if (adapt && tstep != ntstep-1)  adaptation();
      
      /* THE FOLLOWING IS TO CHECK ADAPTATION ERROR */
      // if (!(tstep%out_intrvl)) {
      //   output((tstep+1)*100,tecplot); 
      // }

      hp_mgrid::setbd(MIN(TMSCHEME,tstep+2));
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
#ifdef ENERGY
   FLT e = 0.0, a = 0.0;
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].energy(e,a);
   e /= a;
#else
   FLT e = 1.0;
#endif

   for(i=0;i<nblocks;++i)
      blk[i].grd[0].length1(e);

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

void blocks::output(int number, FILETYPE type, int verbose) {
   int i,j;
   char outname[30], bname[30];

   if (!verbose) {
      for(i=0;i<nblocks;++i) {
         number_str(bname,"data",number,3);
         strcat(bname, ".");
         number_str(outname, bname, i, 1);
   
         /* OUTPUT SOLUTION */
         blk[i].grd[0].output(outname,type);
      }
   }
   else {
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
         for(j=0;j<TMADAPT;++j) {
            number_str(bname,outname,j,1);
            strcat(bname, ".");
            number_str(bname, bname, i, 1);
            blk[i].grd[0].output(blk[i].gbl.ugbd[j],blk[i].gbl.vrtxbd[j],blk[i].gbl.binfobd[j],bname,type);
         }
   
         /* OUTPUT INTERPOLATED MESH TIME HISTORY */
         number_str(outname,"rstrtvrtx",number,3);
         strcat(outname, ".");
         for(j=0;j<TMADAPT;++j) {
            number_str(bname,outname,j,1);
            strcat(bname, ".");
            number_str(bname, bname, i, 1);
            blk[i].grd[0].out_mesh(blk[i].gbl.vrtxbd[j],bname,type);
         }
      }
   }
      
   return;
}

void blocks::minvrt_test(int iter, FLT (*func)(int, FLT, FLT)) {
   int i,stage,mode;
      
   /*****************************************/
   /* NSTAGE UPDATE OF FLOW VARIABLES    ****/
   /*****************************************/
   
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].minvrt_test_tstep();
      
   for(i=0;i<nblocks;++i)
      blk[i].grd[0].tstep_mp();

   for(i=0;i<nblocks;++i)
      blk[i].grd[0].tstep2();
      
   for(stage=0;stage<iter;++stage) {

      /* CALCULATE RESIDUAL */   
      for(i=0;i<nblocks;++i)
         blk[i].grd[0].minvrt_test_bgn(func);
                           
      /* INVERT MASS MATRIX (4 STEP PROCESS) */
      /* HAVE TO BE VERY CAREFUL WITH COMMUNICATION */
      /* USE MESSAGE BEFORE SENDING NEXT */
      /* VERTICES MUST BE 2 PART PROCESS BECAUSE OF CORNERS */
      for(i=0;i<nblocks;++i)
         blk[i].grd[0].minvrt1();  //SEND Y
         
      for(i=0;i<nblocks;++i)
         blk[i].grd[0].bdry_mp(); //RCV Y SEND X
         
      for(i=0;i<nblocks;++i)
         blk[i].grd[0].minvrt2(); //RCV X SEND SIDE MODE 0 Y & X
                  
      for(mode=0;mode<base[lg2pmax].sm-1;++mode) {
         for(i=0;i<nblocks;++i)
            blk[i].grd[0].minvrt3(mode); // RCV SIDE MODE Y & X
            
         for(i=0;i<nblocks;++i)
            blk[i].grd[0].minvrt3_mp(mode+1); //SEND SIDE Y & X
      }
      
      if (base[lg2pmax].sm) {
         for(i=0;i<nblocks;++i)
            blk[i].grd[0].minvrt4(); //RCV Y & X
      }
      
      for(i=0;i<nblocks;++i) 
         blk[i].grd[0].minvrt_test_end();
      
      int n;
      FLT mxr[NV];
      printf("%d ",stage);
      for(i=0;i<nblocks;++i) {
         blk[i].grd[0].maxres(mxr);
         for(n=0;n<ND;++n)
            printf("%.3e  ",mxr[n]);
      }
      printf("\n");
      
   }

   return;
}

   
