#include"mesh.h"
#include"utilities.h"
#include<cstdlib>
#include<cstring>

FLT *mesh<2>::fltwk;
int *mesh<2>::intwk1;
int *mesh<2>::intwk2;
int *mesh<2>::intwk3;
int mesh<2>::gblmaxvst = 0;

FLT *mesh<3>::fltwk;
int *mesh<3>::intwk1;
int *mesh<3>::intwk2;
int *mesh<3>::intwk3;
int mesh<3>::gblmaxvst = 0;

/* CREATE INSTANTIATIONS OF THESE FUNCTIONS */
template int mesh<2>::in_mesh(FLT (*vin)[2], const char *filename, FTYPE filetype, FLT grwfac);
template int mesh<3>::in_mesh(FLT (*vin)[3], const char *filename, FTYPE filetype, FLT grwfac);

template<int ND> int mesh<ND>::in_mesh(FLT (*vin)[ND], const char *filename, FTYPE filetype, FLT grwfac) {
    int i,j,n,sind,count,temp,tind,v0,v1,sign;
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
                    *log << "error: couldn't open file: " << grd_app << std::endl;
                    exit(1);
            }
            ierr = fscanf(grd,"%d\n",&nside);
            if(ierr != 1) {
                    *log << "error: in side file: " << grd_app << std::endl;
                    exit(1);
            }
           
            if (!initialized) {
               allocate(nside + (int) (grwfac*nside));
               vin = vrtx;
            }
            else if (nside > maxvst) {
               *log << "error: mesh is too large" << std::endl;
               exit(1);
            }
           
            count = 0;
            for(i=0;i<nside;++i) {
               ierr = fscanf(grd,"%*d:%d%d%*d%*d%d\n"
               ,&svrtx[i][0],&svrtx[i][1],&sinfo[i]);
               if(ierr != 3) {
                  *log << "error: in side file " << grd_app << std::endl;
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
                     *log << "error: too many different side boundaries: increase MAXSB" << std::endl;
                     exit(1);
                  }
               }
next1:      continue;
            }
               
            /* LOAD VERTEX INFORMATION               */
            strcpy(grd_app,filename);
            strcat(grd_app,".n");
            grd = fopen(grd_app,"r");
            if (!grd) { *log << "trouble opening grid" << grd_app << std::endl; exit(1);}
    
            ierr = fscanf(grd,"%d\n",&nvrtx);
            if(ierr != 1) {
               *log << "1: error in grid: " << grd_app << std::endl;
               exit(1);
            }
            
            /* ERROR %lf SHOULD BE FLT */
            for(i=0;i<nvrtx;++i) {
               fscanf(grd,"%*d:");
               if (DIM == 3) {
                  fgets(grd_app,100,grd);
                  ierr = sscanf(grd_app,"%lf%lf%lf%d",&vin[i][0],&vin[i][1],&vin[i][2],&vinfo[i]);
                  if (ierr != ND+1) {
                     ierr = sscanf(grd_app,"%lf%lf%d",&vin[i][0],&vin[i][1],&vinfo[i]);
                     if (ierr != ND) { *log << "2a: error in grid" << std::endl; exit(1); }
                  }
               }
               else {
                  ierr = 0;
                  for(n=0;n<ND;++n) {
                     ierr += fscanf(grd,"%lf",&vin[i][n]);
                  }
                  ierr += fscanf(grd,"%d",&vinfo[i]);
                  if (ierr != ND+1)  { *log << "2b: error in grid" << std::endl; exit(1); }
               }
            }
            fclose(grd);

            /* COUNT VERTEX BOUNDARY GROUPS  */
            nvbd = 0;
            for(i=0;i<nvrtx;++i) {
               if (vinfo[i]) {
                  /* NEW VRTX B.C. */
                  getnewvrtxobject(nvbd,vinfo[i]);
                  vbdry[nvbd]->alloc(1);
                  vbdry[nvbd]->v() = i;
                  ++nvbd;
                  if (nvbd >= MAXVB) {
                     *log << "Too many vertex boundary conditions: increase MAXSB: " << nvbd << std::endl;
                     exit(1);
                  }
               }
            }
                        
            /* LOAD ELEMENT INFORMATION */
            strcpy(grd_app,filename);
            strcat(grd_app,".e");
            grd = fopen(grd_app,"r");
            if (grd==NULL) {
               *log << "trouble opening " << grd_app << std::endl;
               exit(1);
            }
            ierr = fscanf(grd,"%d\n",&ntri);
            if(ierr != 1) {
               *log << "error in file " << grd_app << std::endl;
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
                    *log << "trouble opening " << filename << std::endl;
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
               *log << "mesh is too large" << std::endl;
               exit(1);
            }
    
            for(i=0;i<8;++i)
                fscanf(grd,"%*[^\n]\n");

            /* READ VERTEX DATA */    
            for(i=0;i<nvrtx;++i) {
               fscanf(grd,"%*d");
               for(n=0;n<ND;++n)
                  fscanf(grd,"%le ",&vin[i][n]);
               fscanf(grd,"\n");
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
            int (*svrtxbtemp[MAXSB])[2];
    
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
                     *log << "error in boundary information " << i << j << std::endl;
                     exit(1);
                  }
                  if (svrtx[sind][1] == svrtxbtemp[i][j][1]) {
                     sinfo[sind] = sbdry[i]->idnty();
                     sbdry[i]->sd(j) = sind;
                     vinfo[svrtx[sind][0]] = 0;
                  }
                  else {
                     *log << "Error: boundary sides are not counterclockwise " << 
                     svrtxbtemp[i][j][0] << svrtxbtemp[i][j][1] << std::endl;
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
                    *log << "couldn't open file: " << grd_app << std::endl;
                    exit(1);
            }
            /* HEADER LINES */
            fscanf(grd,"nvrtx: %d\t nside: %d\t ntri: %d\n",&nvrtx,&nside,&ntri);
            
            if (!initialized) {
               allocate(nside + (int) (grwfac*nside));
               vin = vrtx;
            }
            else if (nside > maxvst) {
               *log << "mesh is too large" << std::endl;
               exit(1);
            }
   
            /* VRTX INFO */                        
            for(i=0;i<nvrtx;++i) {
               for(n=0;n<ND;++n)
               	fscanf(grd,"%lf",&vin[i][n]);
               fscanf(grd,"\n");
            }
                     
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
               fscanf(grd,"idnty: %d\n",&temp);
               getnewsideobject(i,temp);
               sbdry[i]->input(grd,grwfac);
            }
            
            /* VERTEX BOUNDARY INFO HEADER */
            fscanf(grd,"nvbd: %d\n",&nvbd);
            for(i=0;i<nvbd;++i) {
               fscanf(grd,"idnty: %d\n",&temp);
               getnewvrtxobject(i,temp);
               vbdry[i]->alloc(1);
               vbdry[i]->input(grd);
            }
            
            break;
            
         case(mavriplis):
            strcpy(grd_app,filename);
            strcat(grd_app,".mvp");
            grd = fopen(grd_app,"r");
            if (grd==NULL) {
                    *log << "couldn't open file: " << grd_app << std::endl;
                    exit(1);
            }
              
            for(i=0;i<12;++i) {
               fscanf(grd,"%*[^\n]\n");
            }          

            ierr = fscanf(grd,"%d%d%d%*d%*d%d\n",&nsbd,&nvrtx,&nside,&ntri);
            if (ierr != 4) {
               *log << "trouble with mavriplis format" << std::endl;
               exit(1);
            }
            
            if (!initialized) {
               allocate(nside + (int) (grwfac*nside));
               vin = vrtx;
            }
            else if (nside > maxvst) {
               *log << "mesh is too large" << std::endl;
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
            
            for(i=0;i<nvrtx;++i) {
               for(n=0;n<ND;++n)
               	fscanf(grd,"%lf",&vrtx[i][n]);
               fscanf(grd,"%*[^\n]\n");
            }
               
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
               *log << "to read in vertex positions only must first load mesh structure" << std::endl;
               exit(1);
            }

            /* LOAD VERTEX POSITIONS               */
            strcpy(grd_app,filename);
            strcat(grd_app,".txt");
            grd = fopen(grd_app,"r");
            if (!grd) { *log << "trouble opening grid" << grd_app << std::endl; exit(1);}
    
            ierr = fscanf(grd,"%d\n",&temp);
            if(ierr != 1) {
               *log << "1: error in grid " << grd_app << std::endl;
               exit(1);
            }
            if (temp != nvrtx) {
               *log << "grid doesn't match vertex list" << std::endl;
               exit(1);
            }
    
            /* ERROR %lf SHOULD BE FLT */
            for(i=0;i<nvrtx;++i) {
               fscanf(grd,"%*d:");
               ierr = 0;
               for(n=0;n<ND;++n)
                  ierr = fscanf(grd,"%lf",&vin[i][n]);
               fscanf(grd,"\n");
               if (ierr != ND) { *log << "2c: error in grid" << std::endl; exit(1); }
            }
            fclose(grd);
            
            if (vin == vrtx) treeinit();
            
            return(1);

         default:
            *log << "That filetype is not supported" << std::endl;
            exit(1);
   }

   /* ALLOCATE WORK ARRAYS USED BY ALL MESHES */
   /* ALWAYS MORE SIDES THAN TRI'S and VERTICES */    
   if (maxvst > gblmaxvst) {
      if (gblmaxvst > 0) {
         *log << "#Warning: better to allocate from largest mesh to smallest" << std::endl;
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

   *log << "#Boundaries" << std::endl;
   for(i=0;i<nsbd;++i)
      sbdry[i]->summarize(*log);

   createttri();
   createvtri();
   cnt_nbor();
   if (!initialized) qtree.allocate(vin,maxvst);
   treeinit();
   initvlngth();

   initialized = 1;
   
   return(1);
}

template void mesh<2>::allocate(int mxsize);
template void mesh<3>::allocate(int mxsize);

template<int ND> void mesh<ND>::allocate(int mxsize) {
   
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
