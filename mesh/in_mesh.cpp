#include"mesh.h"
#include"utilities.h"
#include<cstdio>
#include<cstdlib>
#include<cstring>

FLT *mesh::fltwk;
int *mesh::intwk1, *mesh::intwk2, *mesh::intwk3;
int mesh::gblmaxvst = 0;

int mesh::in_mesh(char *filename, FILETYPE filetype = easymesh, FLT grwfac = 1) {
    int i,j,sind,count;
    int ierr;
    char grd_app[100];
    FILE *grd;
    int v[3],s[3],e[3];
    
    initialized = 1;
    
    switch (filetype) {            
        case(easymesh):
/*          LOAD SIDE INFORMATION */
            strcpy(grd_app,filename);
            strcat(grd_app,".s");
            grd = fopen(grd_app,"r");
            if (grd==NULL) {
                    printf("couldn't open file: %s\n",grd_app);
                    exit(1);
            }
            ierr = fscanf(grd,"%d\n",&nside);
            if(ierr != 1) {
                    printf("error in side file: %s\n",grd_app);
                    exit(1);
            }
            maxvst = nside + (int) (grwfac*nside);
            svrtx  = (int (*)[2]) xmalloc(maxvst*2*sizeof(int));
            stri   = (int (*)[2]) xmalloc(maxvst*2*sizeof(int));
            sinfo = new int[maxvst+1];
            ++sinfo; // ALLOWS US TO ACCESS SINFO[-1]
            
            for(i=0;i<maxvst;++i) {
                    stri[i][0] = -1;
                    stri[i][1] = -1;
            }
    
            for(i=0;i<nside;++i) {
               ierr = fscanf(grd,"%*d:%d%d%*d%*d%d\n"
               ,&svrtx[i][0],&svrtx[i][1],&sinfo[i]);
               if(ierr != 3) {
                  printf("error in side file %s\n",grd_app);
                  exit(1);
               }
            }
            fclose(grd);	
            
/*          FIGURE OUT HOW MANY BOUNDARY GROUPS THERE ARE & HOW MANY OF EACH */
            nsbd = 0;
            for(i=0;i<nside;++i) {
               if (sinfo[i]) {
                  for (j = 0; j <nsbd;++j) {
                     if (sinfo[i] == sbdry[j].type) {
                              ++sbdry[j].num;
                              goto next1;
                     }
                  }
/*						NEW SIDE */
                  sbdry[nsbd].type = sinfo[i];
                  sbdry[nsbd].num = 1;
                  ++nsbd;
               }
next1:      continue;
            }
            
/*	    		ALLOCATE STORAGE */	
            maxsbel = 0;
            for(i=0;i<nsbd;++i) {
               maxsbel = MAX(sbdry[i].num,maxsbel);
            }
            maxsbel = maxsbel + (int) (grwfac*maxsbel);

            for(i=0;i<nsbd;++i) {
               sbdry[i].el = new int[maxsbel];
               sbdry[i].num = 0;
            }
 
/*				STORE BOUNDARY INFO */
            for(i=0;i<nside;++i) {
               if (sinfo[i]) {
                  for (j = 0; j <nsbd;++j) {
                     if (sinfo[i] == sbdry[j].type) {
                        sbdry[j].el[sbdry[j].num++] = i;
                        break;
                     }
                  }
               }
            }	
               
/*	    		LOAD VERTEX INFORMATION 				  */
            strcpy(grd_app,filename);
            strcat(grd_app,".n");
            grd = fopen(grd_app,"r");
            if (!grd) {printf("trouble opening grid %s\n",grd_app); exit(1);}
    
            ierr = fscanf(grd,"%d\n",&nvrtx);
            if(ierr != 1) {
               printf("1: error in grid %s\n",grd_app);
               exit(1);
            }
            vrtx = (FLT (*)[ND]) xmalloc(maxvst*ND*sizeof(FLT));
            vlngth = new FLT[maxvst];
            vinfo = new int[maxvst+1];
            ++vinfo;  //  ALLOWS US TO ACCES VINFO[-1]
            nnbor = new int[maxvst];
            vtri = new int[maxvst];
    
/*	    		ERROR %lf SHOULD BE FLT */
            for(i=0;i<nvrtx;++i) {
               ierr = fscanf(grd,"%*d:%lf%lf%d\n",&vrtx[i][0],&vrtx[i][1],&vinfo[i]);
               if (ierr != 3) { printf("2: error in grid\n"); exit(1); }
            }
            fclose(grd);

/*	    		FIGURE OUT HOW MANY VERTEX BOUNDARY GROUPS THERE ARE & HOW MANY OF EACH */
            nvbd = 0;
            for(i=0;i<nvrtx;++i) {
               if (vinfo[i]) {
                  for (j = 0; j <nvbd;++j) {
                     if (vinfo[i] == vbdry[j].type) {
                        ++vbdry[j].num;
                        goto next2;
                     }
                  }
/*			    		NEW SIDE */
                  vbdry[nvbd].type = vinfo[i];
                  vbdry[nvbd].num = 1;
                  ++nvbd;
               }
next2:      continue;
            }

/*				ALLOCATE STORAGE */    
            maxvbel = 0;
            for(i=0;i<nvbd;++i) {
               maxvbel = MAX(vbdry[i].num,maxvbel);	
            }
            maxvbel = maxvbel + (int) (grwfac*maxvbel);
    
/*	    		ALLOCATE STORAGE */	
            for(i=0;i<nvbd;++i) {
               vbdry[i].el = new int[maxvbel];
               vbdry[i].num = 0;
            }
    
            for(i=0;i<nvrtx;++i) {
               if (vinfo[i]) {
                  for (j = 0; j <nvbd;++j) {
                     if (vinfo[i] == vbdry[j].type) {
                        vbdry[j].el[vbdry[j].num++] = i;
                        break;
                     }
                  }
               }
            }
                        
/*	    		LOAD ELEMENT INFORMATION */
            strcpy(grd_app,filename);
            strcat(grd_app,".e");
            grd = fopen(grd_app,"r");
            if (grd==NULL) {
               printf("trouble opening %s\n",grd_app);
               exit(1);
            }
            ierr = fscanf(grd,"%d\n",&ntri);
            if(ierr != 1) {
               printf("error in file %s\n",grd_app);
               exit(1);
            }
            tvrtx = (int (*)[3]) xmalloc(maxvst*3*sizeof(int));
            ttri = (int (*)[3]) xmalloc(maxvst*3*sizeof(int));
            tside = new struct tsidedata[maxvst];
    
            tinfo = new int[maxvst+1];
            ++tinfo; // ALLOWS US TO ACCESS TINFO[-1]
    
            for(i=0;i<ntri;++i) {
                ierr = fscanf(grd,"%*d:%d%d%d%d%d%d%d%d%d%*f%*f%d\n"
                ,v,v+1,v+2,e,e+1,e+2,s,s+1,s+2,&tinfo[i]);
    
                for (j=0;j<3;++j) {
                    tvrtx[i][j] = v[j];
                    if(svrtx[s[j]][0] == v[(j+1)%3]) {
/* 							SIDE DEFINED IN SAME DIRECTION */
                        stri[s[j]][0] = i;
                        tside[i].side[j] = s[j];
                        tside[i].sign[j] = 1;
                    }
                    else {
/* 							SIDE DEFINED IN OPPOSITE DIRECTION */
                        stri[s[j]][1] = i;
                        tside[i].side[j] = s[j];
                        tside[i].sign[j] = -1;
                    }
                }
            }
            break;
        
        case(gambit):
            strcpy(grd_app,filename);
            strcat(grd_app,".FDNEUT");
            grd = fopen(grd_app,"r");
            if (grd == NULL) {
                    printf("trouble opening %s\n",filename);
                    exit(1);
            }
            
            for(i=0;i<5;++i)
                fscanf(grd,"%*[^\n]\n");
                
            fscanf(grd,"%d%d%d\n",&nvrtx,&maxvst,&nsbd);
            nsbd -= 1;
 /*			MAXVST IS APPROXIMATELY THE NUMBER OF ELEMENTS  */
 /*			FOR EACH ELEMENT THERE ARE APPROXIMATELY 3/2 SIDES */
            maxvst = 3*maxvst/2; 
            maxvst = nvrtx +(int) (grwfac*nvrtx);
            vrtx = (FLT (*)[ND]) xmalloc(maxvst*ND*sizeof(FLT));
            vlngth = new FLT[maxvst];
            vinfo = new int[maxvst+1];
            ++vinfo;  //  ALLOWS US TO ACCES VINFO[-1]
            nnbor = new int[maxvst];
            vtri = new int[maxvst];
    
            for(i=0;i<8;++i)
                fscanf(grd,"%*[^\n]\n");

/*				READ VERTEX DATA */    
            for(i=0;i<nvrtx;++i) {
                fscanf(grd,"%*d %le %le\n",&vrtx[i][0],&vrtx[i][1]);
                vinfo[i] = -1;
            }
                
            for(i=0;i<2;++i)
                fscanf(grd,"%*[^\n]\n");

/*				READ ELEMENT DATA */
            fscanf(grd,"%*[^0-9]%*d%*[^0-9]%d%*[^\n]\n",&ntri);
            tvrtx = (int (*)[3]) xmalloc(maxvst*3*sizeof(int));
            ttri = (int (*)[3]) xmalloc(maxvst*3*sizeof(int));
            tside = new struct tsidedata[maxvst];
            tinfo = new int[maxvst+1];
            ++tinfo;

            fscanf(grd,"%*[^\n]");
                    
            for(i=0;i<ntri;++i) {
                fscanf(grd,"%*d %d %d %d",&tvrtx[i][0],&tvrtx[i][1],&tvrtx[i][2]);

                --tvrtx[i][0];
                --tvrtx[i][1];
                --tvrtx[i][2];
                tinfo[i] = 0;
            }

/*				READ BOUNDARY DATA STORE TEMPORARILY */            
            maxsbel = 0;
            int (*svrtxbtemp[MAXSB])[ND];
    
            for(i=0;i<nsbd;++i) {
               fscanf(grd,"%*[^0-9]%*d%*[^0-9]%d%*[^0-9]%*d%*[^0-9]%*d%*[^0-9]%d\n"
                  ,&sbdry[i].num,&sbdry[i].type);
               
               maxsbel = MAX(sbdry[i].num,maxsbel);
               fscanf(grd,"%*[^\n]\n");
               
               svrtxbtemp[i] = (int (*)[2]) xmalloc(sbdry[i].num*2*sizeof(int));
                    
               for(j=0;j<sbdry[i].num;++j) {
                  fscanf(grd,"%*d %d %d\n"
                     ,&svrtxbtemp[i][j][0],&svrtxbtemp[i][j][1]);
                  --svrtxbtemp[i][j][0];
                  --svrtxbtemp[i][j][1];
                  vinfo[svrtxbtemp[i][j][0]] = sbdry[i].type;
                  vinfo[svrtxbtemp[i][j][1]] = sbdry[i].type;
               }
            }
                                    
/*	    		ALLOCATE BOUNDARY STORAGE */
            maxsbel = maxsbel + (int) (grwfac*maxsbel);	
            for(i=0;i<nsbd;++i) 
               sbdry[i].el = new int[maxsbel];
    
            count = 0;
            for(i=0;i<nsbd;++i) 
               count += sbdry[i].num;
            
            svrtx  = (int (*)[2]) xmalloc(maxvst*2*sizeof(int));
            stri   = (int (*)[2]) xmalloc(maxvst*2*sizeof(int));
            sinfo = new int[maxvst+1];
            ++sinfo; // ALLOWS US TO ACCESS SINFO[-1]


/*				CREATE SIDE INFORMATION */
            createsideinfo();
                
/*				FIND ALL BOUNDARY SIDES */
/*				STORE LOCATION BY VERTEX NUMBER */
            for(i=0;i<nside;++i) 
               if (stri[i][1] < 0)
                  vinfo[svrtx[i][0]] = i;
                  
            for(i=0;i<nside;++i)
               sinfo[i] = 0;

/*				MATCH BOUNDARY SIDES TO GROUPS */    
            for(i=0;i<nsbd;++i) {
               for(j=0;j<sbdry[i].num;++j) {
                  sind = vinfo[svrtxbtemp[i][j][0]];
                  if (sind < 0) {
                     printf("error in boundary information %d %d\n",i,j);
                     exit(1);
                  }
                  if (svrtx[sind][1] == svrtxbtemp[i][j][1]) {
                     sinfo[sind] = sbdry[i].type;
                     sbdry[i].el[j] = sind;
                     vinfo[svrtx[sind][0]] = 0;
                  }
                  else {
                     printf("Error: boundary sides are not counterclockwise %d %d\n",
                        svrtxbtemp[i][j][0],svrtxbtemp[i][j][1]);
                     exit(1);
                  }
               }
            }
            
            for(i=0;i<nvrtx;++i)
               vinfo[i] = 0;

            break;
            
         default:
            printf("That filetype is not supported\n");
            exit(1);
   }

/*	ALLOCATE WORK ARRAYS USED BY ALL MESHES */
/*	ALWAYS MORE SIDES THAN TRI'S and VERTICES */    
   if (maxvst > gblmaxvst) {
      if (gblmaxvst > 0) {
         printf("better to allocate from largest mesh to smallest\n");
         exit(1);
      }
      fltwk = new FLT[maxvst];
      intwk1 = new int[maxvst];
      intwk2 = new int[maxvst];
      intwk3 = new int[maxvst];
      gblmaxvst = maxvst;

/*		INTWK1,2 SHOULD ALWAYS BE RESET TO NEGATIVE 1 AFTER USE */
      for(i=0;i<maxvst;++i)
         intwk1[i] = -1;
      for(i=0;i<maxvst;++i)
         intwk2[i] = -1;
      for(i=0;i<maxvst;++i)
         intwk3[i] = -1;
   }
    
/*	REORDER SIDE BOUNDARY POINTERS TO BE SEQUENTIAL */
   for(i=0;i<nsbd;++i) 
      bdrysidereorder(i);

/*	REORDER SIDE BOUNDARY GROUPS TO BE SEQUENTIAL */
   bdrygroupreorder();

   createttri();
   createvtri();
   cnt_nbor();
   qtree.allocate(vrtx,maxvst);
   treeinit();
   initvlngth();

   return(1);
}

/* MAPPING FROM OLD DEFINITIONS TO NEW */
#define NOLDBTYPE 15

#define FSRF_OLD 1
#define FCE1_OLD 2
#define FCE2_OLD 4
#define INFC_OLD 8
#define INFS_OLD 32
#define SYMC_OLD 16
#define PDX1_OLD 64
#define PDX2_OLD 128
#define FXMV_OLD 256
#define OUTF_OLD 512
#define INVC_OLD 1024
#define INVS_OLD 2048
#define PDY1_OLD 4096
#define PDY2_OLD 8192
#define OUT2_OLD 513

void mesh::convertbtypes(const int (*old)[2] = NULL, int nold = 0) {

   int oldbtype[NOLDBTYPE][2] = 
                             {{FSRF_OLD,FSRF_MASK+CURV_MASK},
                             {FCE1_OLD,IFCE_MASK+CURV_MASK},
                             {FCE2_OLD,IFCE_MASK+CURV_MASK},
                             {INFC_OLD,INFL_MASK+CURV_MASK},
                             {INVC_OLD,EULR_MASK+CURV_MASK},
                             {INFS_OLD,INFL_MASK},
                             {FXMV_OLD+INFS_OLD,INFL_MASK},
                             {SYMC_OLD,SYMM_MASK},
                             {PDX1_OLD,PRDX_MASK+COMX_MASK},
                             {PDX2_OLD,PRDX_MASK+COMX_MASK},
                             {OUTF_OLD,OUTF_MASK},
                             {INVS_OLD,EULR_MASK},
                             {PDY1_OLD,PRDY_MASK+COMY_MASK},
                             {PDY2_OLD,PRDY_MASK+COMY_MASK},
                             {OUT2_OLD,OUTF_MASK+(1<<16)}};
   int i,j;
   
   if (old == NULL) {
      old = oldbtype;
      nold = NOLDBTYPE;
   }
   
   for(i=0;i<nvbd;++i) {
      for(j=0;j<NOLDBTYPE;++j) {
         if (vbdry[i].type == old[j][0]) {
            vbdry[i].type = old[j][1];
            break;
         }
      }
   }
   
   for(i=0;i<nsbd;++i) {
      for(j=0;j<NOLDBTYPE;++j) {
         if (sbdry[i].type == old[j][0]) {
            sbdry[i].type = old[j][1];
            break;
         }
      }
   }   
   return;
}