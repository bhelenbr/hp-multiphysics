#include"mesh.h"
#include"utilities.h"
#include<cstdio>
#include<cstdlib>
#include<cstring>

FLT *mesh::fltwk;
int *mesh::intwk1, *mesh::intwk2,*mesh::intwk3;
int mesh::gblmaxvst = 0;

int mesh::in_mesh(FLT (*vin)[ND], const char *filename, FILETYPE filetype, FLT grwfac) {
    int i,j,sind,count,temp,tind,v0,v1,sign;
    int ierr;
    char grd_app[100];
    FILE *grd;
    int v[3],s[3],e[3];
        
    switch (filetype) {            
        case(easymesh):
            /* LOAD SIDE INFORMATION */
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
           
            if (!initialized) {
               allocate(nside + (int) (grwfac*nside));
               vin = vrtx;
            }
            else if (nside > maxvst) {
               printf("mesh is too large\n");
               exit(1);
            }
           
            count = 0;
            for(i=0;i<nside;++i) {
               ierr = fscanf(grd,"%*d:%d%d%*d%*d%d\n"
               ,&svrtx[i][0],&svrtx[i][1],&sinfo[i]);
               if(ierr != 3) {
                  printf("error in side file %s\n",grd_app);
                  exit(1);
               }
               if (sinfo[i]) ++count;
            }
            fclose(grd);
            
            for(i=0;i<maxvst;++i) {
               stri[i][0] = -1;
               stri[i][1] = -1;
            }
            
            /* ORGANIZE BOUNDARY GROUPS */
            nsbd = 0;
            for(i=0;i<nside;++i) {
               if (sinfo[i]) {
                  for (j = 0; j <nsbd;++j) {
                     if (sinfo[i] == sbdry[j]->idnty()) {
                        sbdry[j]->sd(sbdry[j]->nsd()++) = i;
                        goto next1;
                     }
                  }
                  /* NEW SIDE */
                  getnewsideobject(nsbd,sinfo[i]);
                  sbdry[nsbd]->alloc(static_cast<int>(count*grwfac));
                  sbdry[nsbd]->nsd() = 1;
                  sbdry[nsbd++]->sd(0)= i;
                  if (nsbd > MAXSB) {
                     printf("too many different side boundaries: increase MAXSB\n");
                     exit(1);
                  }
               }
next1:      continue;
            }
               
            /* LOAD VERTEX INFORMATION               */
            strcpy(grd_app,filename);
            strcat(grd_app,".n");
            grd = fopen(grd_app,"r");
            if (!grd) {printf("trouble opening grid %s\n",grd_app); exit(1);}
    
            ierr = fscanf(grd,"%d\n",&nvrtx);
            if(ierr != 1) {
               printf("1: error in grid %s\n",grd_app);
               exit(1);
            }
    
            /* ERROR %lf SHOULD BE FLT */
            for(i=0;i<nvrtx;++i) {
               ierr = fscanf(grd,"%*d:%lf%lf%d\n",&vin[i][0],&vin[i][1],&vinfo[i]);
               if (ierr != 3) { printf("2: error in grid\n"); exit(1); }
            }
            fclose(grd);

            /* THIS IS GOING TO HAVE TO CHANGE */
            /* COUNT VERTEX BOUNDARY GROUPS  */
            nvbd = 0;
            for(i=0;i<nvrtx;++i) {
               if (vinfo[i]) {
                  /* NEW VRTX B.C. */
                  vbdry[nvbd].type = vinfo[i];
                  vbdry[nvbd].num = 1;
                  vbdry[nvbd].el[0] = i;
                  ++nvbd;
                  if (nvbd >= MAXSB) {
                     printf("Too many vertex boundary conditions: increase MAXSB %d\n",nvbd);
                     exit(1);
                  }
               }
            }
                        
            /* LOAD ELEMENT INFORMATION */
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
    
            for(i=0;i<ntri;++i) {
                ierr = fscanf(grd,"%*d:%d%d%d%d%d%d%d%d%d%*f%*f%d\n"
                ,v,v+1,v+2,e,e+1,e+2,s,s+1,s+2,&tinfo[i]);
    
                for (j=0;j<3;++j) {
                    tvrtx[i][j] = v[j];
                    if(svrtx[s[j]][0] == v[(j+1)%3]) {
                        /* SIDE DEFINED IN SAME DIRECTION */
                        stri[s[j]][0] = i;
                        tside[i].side[j] = s[j];
                        tside[i].sign[j] = 1;
                    }
                    else {
                        /* SIDE DEFINED IN OPPOSITE DIRECTION */
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
                
            fscanf(grd,"%d%d%d\n",&nvrtx,&i,&nsbd);
            nsbd -= 1;
            /* MAXVST IS APPROXIMATELY THE NUMBER OF ELEMENTS  */
            /* FOR EACH ELEMENT THERE ARE APPROXIMATELY 3/2 SIDES */
            if (!initialized) {
               maxvst = (3*i)/2; 
               maxvst = nvrtx +(int) (grwfac*nvrtx);
               allocate(maxvst);
               vin = vrtx;
            }
            else if ((3*i)/2 > maxvst) {
               printf("mesh is too large\n");
               exit(1);
            }
    
            for(i=0;i<8;++i)
                fscanf(grd,"%*[^\n]\n");

            /* READ VERTEX DATA */    
            for(i=0;i<nvrtx;++i) {
                fscanf(grd,"%*d %le %le\n",&vin[i][0],&vin[i][1]);
                vinfo[i] = -1;
            }
                
            for(i=0;i<2;++i)
                fscanf(grd,"%*[^\n]\n");

            /* READ ELEMENT DATA */
            fscanf(grd,"%*[^0-9]%*d%*[^0-9]%d%*[^\n]\n",&ntri);
            fscanf(grd,"%*[^\n]");
                    
            for(i=0;i<ntri;++i) {
                fscanf(grd,"%*d %d %d %d",&tvrtx[i][0],&tvrtx[i][1],&tvrtx[i][2]);
                --tvrtx[i][0];
                --tvrtx[i][1];
                --tvrtx[i][2];
                tinfo[i] = 0;
            }

            /* READ BOUNDARY DATA STORE TEMPORARILY */            
            int (*svrtxbtemp[MAXSB])[ND];
    
            for(i=0;i<nsbd;++i) {
               fscanf(grd,"%*[^0-9]%*d%*[^0-9]%d%*[^0-9]%*d%*[^0-9]%*d%*[^0-9]%d\n"
                  ,&count,&temp);
               
               getnewsideobject(i,temp);
               sbdry[i]->alloc(static_cast<int>(count*grwfac));
               sbdry[i]->nsd() = count;
               
               fscanf(grd,"%*[^\n]\n");
               
               svrtxbtemp[i] = (int (*)[2]) xmalloc(sbdry[i]->nsd()*2*sizeof(int));
                    
               for(j=0;j<sbdry[i]->nsd();++j) {
                  fscanf(grd,"%*d %d %d\n"
                     ,&svrtxbtemp[i][j][0],&svrtxbtemp[i][j][1]);
                  --svrtxbtemp[i][j][0];
                  --svrtxbtemp[i][j][1];
                  vinfo[svrtxbtemp[i][j][0]] = sbdry[i]->idnty();
                  vinfo[svrtxbtemp[i][j][1]] = sbdry[i]->idnty();
               }
            }

            /* CREATE SIDE INFORMATION */
            createsideinfo();
                
            /* FIND ALL BOUNDARY SIDES */
            /* STORE LOCATION BY VERTEX NUMBER */
            for(i=0;i<nside;++i) 
               if (stri[i][1] < 0)
                  vinfo[svrtx[i][0]] = i;
                  
            for(i=0;i<nside;++i)
               sinfo[i] = 0;

            /* MATCH BOUNDARY SIDES TO GROUPS */    
            for(i=0;i<nsbd;++i) {
               for(j=0;j<sbdry[i]->nsd();++j) {
                  sind = vinfo[svrtxbtemp[i][j][0]];
                  if (sind < 0) {
                     printf("error in boundary information %d %d\n",i,j);
                     exit(1);
                  }
                  if (svrtx[sind][1] == svrtxbtemp[i][j][1]) {
                     sinfo[sind] = sbdry[i]->idnty();
                     sbdry[i]->sd(j) = sind;
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
            
         case(grid):
            strcpy(grd_app,filename);
            strcat(grd_app,".grd");
            grd = fopen(grd_app,"r");
            if (grd==NULL) {
                    printf("couldn't open file: %s\n",grd_app);
                    exit(1);
            }
            /* HEADER LINES */
            fscanf(grd,"nvrtx: %d\t nside: %d\t ntri: %d\n",&nvrtx,&nside,&ntri);
            
            if (!initialized) {
               allocate(nside + (int) (grwfac*nside));
               vin = vrtx;
            }
            else if (nside > maxvst) {
               printf("mesh is too large\n");
               exit(1);
            }
   
            /* VRTX INFO */                        
            for(i=0;i<nvrtx;++i)
               fscanf(grd,"%lf %lf\n",&vin[i][0],&vin[i][1]);
                     
            /* SIDE INFO */
            for(i=0;i<nside;++i)
               fscanf(grd,"%d %d\n",&svrtx[i][0],&svrtx[i][1]);
   
            /* THEN TRI INFO */
            for(i=0;i<ntri;++i)
               fscanf(grd,"%d %d %d\n",&tvrtx[i][0],&tvrtx[i][1],&tvrtx[i][2]);
               
            /* CREATE TSIDE & STRI */
            createtsidestri();
   
            /* SIDE BOUNDARY INFO HEADER */
            fscanf(grd,"nsbd: %d\n",&nsbd);

            count = 0;
            for(i=0;i<nsbd;++i) {
               fscanf(grd,"type: %d\n",&temp);
               getnewsideobject(i,temp);
               sbdry[i]->input(grd,grwfac);
            }
                  
            /* VERTEX BOUNDARY INFO HEADER */
            fscanf(grd,"nvbd: %d\n",&nvbd);
            for(i=0;i<nvbd;++i)
               fscanf(grd,"%d %d\n",&vbdry[i].type,&vbdry[i].num);
            
            for(i=0;i<nvbd;++i)
               for(j=0;j<vbdry[i].num;++j)
                  fscanf(grd,"%d\n",&vbdry[i].el[j]);
            
            break;
            
         case(mavriplis):
            strcpy(grd_app,filename);
            strcat(grd_app,".mvp");
            grd = fopen(grd_app,"r");
            if (grd==NULL) {
                    printf("couldn't open file: %s\n",grd_app);
                    exit(1);
            }
              
            for(i=0;i<12;++i) {
               fscanf(grd,"%*[^\n]\n");
            }          

            ierr = fscanf(grd,"%d%d%d%*d%*d%d\n",&nsbd,&nvrtx,&nside,&ntri);
            if (ierr != 4) {
               printf("trouble with mavriplis format\n");
               exit(1);
            }
            
            if (!initialized) {
               allocate(nside + (int) (grwfac*nside));
               vin = vrtx;
            }
            else if (nside > maxvst) {
               printf("mesh is too large\n");
               exit(1);
            }
            
            fscanf(grd,"%*[^\n]\n");
            
            fscanf(grd,"%*d%d%*d%d%*d%*d\n",&temp,&count);
            
            /* EXTERNAL BOUNDARY */
            getnewsideobject(0,OUTF_MASK);
            sbdry[0]->alloc(static_cast<int>(grwfac*temp));
            sbdry[0]->nsd() = temp;
            for(i=0;i<sbdry[0]->nsd();++i)
               sbdry[0]->sd(i) = i;
            
            ++nsbd;
            
            for(i=1;i<nsbd;++i) {
               getnewsideobject(i,EULR_MASK+CURV_MASK);
               fscanf(grd,"%*[^\n]\n");
               fscanf(grd,"%d%d%*[^\n]\n",&temp,&sbdry[i]->nsd());
               sbdry[i]->alloc(static_cast<int>(grwfac*sbdry[i]->nsd()));
               sbdry[i]->sd(0) = temp-1;
               sbdry[i]->nsd() -= sbdry[i]->sd(0);
               for(j=1;j<sbdry[i]->nsd();++j)
                  sbdry[i]->sd(j) = j +sbdry[i]->sd(0);
            }
            
            fscanf(grd,"%*[^\n]\n");
               
            for(i=0;i<nside;++i) {
               fscanf(grd,"%d%d%d%d%*d\n",&svrtx[i][0],&svrtx[i][1],&stri[i][0],&stri[i][1]);
               --svrtx[i][0];--svrtx[i][1];--stri[i][0];--stri[i][1];
               
               if (stri[i][1] >= ntri) stri[i][1] = -1;
            }
            
            for(i=0;i<nvrtx;++i)
               fscanf(grd,"%lf%lf%*[^\n]\n",&vrtx[i][0],&vrtx[i][1]);
               
            for(i=0;i<ntri;++i)
               for(j=0;j<3;++j)
                  tside[i].side[j] = -1;
                  
            for(i=0;i<nside;++i) {
               tind = stri[i][0];
               j = 0;
               while (tside[tind].side[j] > 0)
                  ++j;
               tside[tind].side[j] = i;
               tside[tind].sign[j] = 1;
               
               tind = stri[i][1];
               if (tind > -1) {
                  j = 0;
                  while (tside[tind].side[j] > 0)
                     ++j;
                  tside[tind].side[j] = i;
                  tside[tind].sign[j] = -1;
               }
            }
            
            /* REORDER SIDES TO BE COUNTERCLOCKWISE */
            /* FILL IN TVRTX */
            for(tind=0;tind<ntri;++tind) {
               v0 = svrtx[tside[tind].side[0]][(tside[tind].sign[0]+1)/2];
               v1 = svrtx[tside[tind].side[1]][(1-tside[tind].sign[1])/2];
               if (v0 != v1) {
                  /* SWITCH SIDES */
                  j = tside[tind].side[1];
                  tside[tind].side[1] = tside[tind].side[2];
                  tside[tind].side[2] = j;
                  j = tside[tind].sign[1];
                  tside[tind].sign[1] = tside[tind].sign[2];
                  tside[tind].sign[2] = j;
               }
               
               tvrtx[tind][2] = v0;
               sind = tside[tind].side[2];
               sign = tside[tind].sign[2];
               tvrtx[tind][1] = svrtx[sind][(sign+1)/2];
               tvrtx[tind][0] = svrtx[sind][(1-sign)/2];
            }
            nvbd = 0;
                        
            break;
            
         case(text):
            if (!initialized) {
               printf("to read in vertex positions only must first load mesh structure\n");
               exit(1);
            }

            /* LOAD VERTEX POSITIONS               */
            strcpy(grd_app,filename);
            strcat(grd_app,".txt");
            grd = fopen(grd_app,"r");
            if (!grd) {printf("trouble opening grid %s\n",grd_app); exit(1);}
    
            ierr = fscanf(grd,"%d\n",&temp);
            if(ierr != 1) {
               printf("1: error in grid %s\n",grd_app);
               exit(1);
            }
            if (temp != nvrtx) {
               printf("grid doesn't match vertex list\n");
               exit(1);
            }
    
            /* ERROR %lf SHOULD BE FLT */
            for(i=0;i<nvrtx;++i) {
               ierr = fscanf(grd,"%*d:%lf%lf\n",&vin[i][0],&vin[i][1]);
               if (ierr != 2) { printf("2: error in grid\n"); exit(1); }
            }
            fclose(grd);
            
            if (vin == vrtx) treeinit();
            
            return(1);

         default:
            printf("That filetype is not supported\n");
            exit(1);
   }

   /* ALLOCATE WORK ARRAYS USED BY ALL MESHES */
   /* ALWAYS MORE SIDES THAN TRI'S and VERTICES */    
   if (maxvst > gblmaxvst) {
      if (gblmaxvst > 0) {
         printf("#Warning: better to allocate from largest mesh to smallest\n");
      }
      fltwk = new FLT[maxvst];
      intwk1 = new int[maxvst];
      intwk2 = new int[maxvst];
      intwk3 = new int[maxvst];
      gblmaxvst = maxvst;

      /* INTWK1,2 SHOULD ALWAYS BE RESET TO NEGATIVE 1 AFTER USE */
      for(i=0;i<maxvst;++i)
         intwk1[i] = -1;
      for(i=0;i<maxvst;++i)
         intwk2[i] = -1;
      for(i=0;i<maxvst;++i)
         intwk3[i] = -1;
   }
    
   for(i=0;i<nsbd;++i) {
      /* CREATES NEW BOUNDARY FOR DISCONNECTED SEGMENTS OF SAME TYPE */
      sbdry[i]->reorder();
      sbdry[i]->getgeometryfrommesh();
   }
   
   bdrylabel();  // CHANGES STRI / TTRI ON BOUNDARIES TO POINT TO GROUP/ELEMENT

   printf("#Boundaries\n");
   for(i=0;i<nsbd;++i)
      sbdry[i]->summarize();

   createttri();
   createvtri();
   cnt_nbor();
   if (!initialized) qtree.allocate(vin,maxvst);
   treeinit();
   initvlngth();

   initialized = 1;
   
   return(1);
}

void mesh::allocate(int mxsize) {
   
   /* SIDE INFO */
   maxvst = mxsize;
   svrtx  = (int (*)[2]) xmalloc(maxvst*2*sizeof(int));
   stri   = (int (*)[2]) xmalloc(maxvst*2*sizeof(int));
   sinfo = new int[maxvst+1];
   ++sinfo; // ALLOWS US TO ACCESS SINFO[-1]

   /* VERTEX INFO */                  
   vrtx = (FLT (*)[ND]) xmalloc(maxvst*ND*sizeof(FLT));
   vlngth = new FLT[maxvst];
   vinfo = new int[maxvst+1];
   ++vinfo;  //  ALLOWS US TO ACCES VINFO[-1]
   nnbor = new int[maxvst];
   vtri = new int[maxvst];
   
   /* TRI INFO */               
   tvrtx = (int (*)[3]) xmalloc(maxvst*3*sizeof(int));
   ttri = (int (*)[3]) xmalloc(maxvst*3*sizeof(int));
   tside = new struct tsidedata[maxvst];
   tinfo = new int[maxvst+1];
   ++tinfo; // ALLOWS US TO ACCESS TINFO[-1]

   return;
}
