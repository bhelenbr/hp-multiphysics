#include "mesh.h"
#include "boundary.h"
#include "utilities.h"
#include <cstdlib>
#include <cstring>
#include <input_map.h>

int *mesh::i1wk = 0, *mesh::i2wk = 0, *mesh::i3wk = 0;
int mesh::maxlst, mesh::maxsrch;

sharedmem* mesh::input(const char *filename, ftype::name filetype, const char *bdryfile, FLT grwfac, sharedmem *win) {
   int i,j,n,sind,count,temp,tind,v0,v1,sign;
   int ierr;
   char grd_app[100];
   FILE *grd;
   int v[3],s[3],e[3];
   int bcntr[MAXSB];
   std::map<std::string,std::string> bdrymap;
   std::map<std::string,std::string> *pbdrymap = &bdrymap;
   
   if (bdryfile) {
      input_map(bdrymap,bdryfile);
   }
   else {
      strcpy(grd_app,filename);
      strcat(grd_app,".bdry.inpt");
      grd = fopen(grd_app,"r");
      if (grd) {
         fclose(grd);
         input_map(bdrymap,grd_app);
      }
      else {
         pbdrymap = 0;
      }
   }
        
    switch (filetype) {            
        case(ftype::easymesh):
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
               allocate(nside + (int) (grwfac*nside),win);
            }
            else if (nside > maxvst) {
               *log << "error: mesh is too large" << std::endl;
               exit(1);
            }
           
            for(i=0;i<nside;++i) {
               ierr = fscanf(grd,"%*d:%d%d%*d%*d%d\n"
               ,&sd[i].vrtx[0],&sd[i].vrtx[1],&sd[i].info);
               if(ierr != 3) {
                  *log << "error: in side file " << grd_app << std::endl;
                  exit(1);
               }
            }
            fclose(grd);
            
            for(i=0;i<maxvst;++i) {
               sd[i].tri[0] = -1;
               sd[i].tri[1] = -1;
            }
            
            /* COUNT BOUNDARY GROUPS */
            nsbd = 0;
            for(i=0;i<nside;++i) {
               if (sd[i].info) {
                  for (j = 0; j <nsbd;++j) {
                     if (sd[i].info == sbdry[j]->idnum) {
                        ++bcntr[j];
                        goto next1;
                     }
                  }
                  /* NEW SIDE */
                  sbdry[nsbd] = getnewsideobject(sd[i].info,pbdrymap);
                  bcntr[nsbd++] = 1;
                  if (nsbd > MAXSB) {
                     *log << "error: too many different side boundaries: increase MAXSB" << std::endl;
                     exit(1);
                  }
               }
next1:      continue;
            }
            
            for(i=0;i<nsbd;++i) {
               sbdry[i]->alloc(static_cast<int>(bcntr[i]*grwfac));
               sbdry[i]->nel = 0;
            }
            
            
            for(i=0;i<nside;++i) {
               if (sd[i].info) {
                  for (j = 0; j <nsbd;++j) {
                     if (sd[i].info == sbdry[j]->idnum) {
                        sbdry[j]->el[sbdry[j]->nel++] = i;
                        goto next1a;
                     }
                  }
                  printf("Big error\n");
                  exit(1);
               }
next1a:     continue;
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
                  ierr = sscanf(grd_app,"%lf%lf%lf%d",&vrtx[i][0],&vrtx[i][1],&vrtx[i][2],&vd[i].info);
                  if (ierr != ND+1) {
                     ierr = sscanf(grd_app,"%lf%lf%d",&vrtx[i][0],&vrtx[i][1],&vd[i].info);
                     if (ierr != ND) { *log << "2a: error in grid" << std::endl; exit(1); }
                  }
               }
               else {
                  ierr = 0;
                  for(n=0;n<ND;++n) {
                     ierr += fscanf(grd,"%lf",&vrtx[i][n]);
                  }
                  ierr += fscanf(grd,"%d",&vd[i].info);
                  if (ierr != ND+1)  { *log << "2b: error in grid" << std::endl; exit(1); }
               }
            }
            fclose(grd);

            /* COUNT VERTEX BOUNDARY GROUPS  */
            nvbd = 0;
            for(i=0;i<nvrtx;++i) {
               if (vd[i].info) {
                  /* NEW VRTX B.C. */
                  vbdry[nvbd] = getnewvrtxobject(vd[i].info,pbdrymap);
                  vbdry[nvbd]->alloc(1);
                  vbdry[nvbd]->v0 = i;
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
                ,v,v+1,v+2,e,e+1,e+2,s,s+1,s+2,&td[i].info);
    
                for (j=0;j<3;++j) {
                    td[i].vrtx[j] = v[j];
                    if(sd[s[j]].vrtx[0] == v[(j+1)%3]) {
                        /* SIDE DEFINED IN SAME DIRECTION */
                        sd[s[j]].tri[0] = i;
                        td[i].side[j] = s[j];
                        td[i].sign[j] = 1;
                    }
                    else {
                        /* SIDE DEFINED IN OPPOSITE DIRECTION */
                        sd[s[j]].tri[1] = i;
                        td[i].side[j] = s[j];
                        td[i].sign[j] = -1;
                    }
                }
            }
            break;
        
        case(ftype::gambit):
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
               maxvst = static_cast<int>((grwfac*3*i)/2); 
               allocate(maxvst,win);
               vrtx = vrtx;
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
                  fscanf(grd,"%le ",&vrtx[i][n]);
               fscanf(grd,"\n");
               vd[i].info = -1;
            }
                
            for(i=0;i<2;++i)
                fscanf(grd,"%*[^\n]\n");

            /* READ ELEMENT DATA */
            fscanf(grd,"%*[^0-9]%*d%*[^0-9]%d%*[^\n]\n",&ntri);
            fscanf(grd,"%*[^\n]");
                    
            for(i=0;i<ntri;++i) {
                fscanf(grd,"%*d %d %d %d",&td[i].vrtx[0],&td[i].vrtx[1],&td[i].vrtx[2]);
                --td[i].vrtx[0];
                --td[i].vrtx[1];
                --td[i].vrtx[2];
                td[i].info = 0;
            }

            /* READ BOUNDARY DATA STORE TEMPORARILY */            
            int (*svrtxbtemp[MAXSB])[2];
    
            for(i=0;i<nsbd;++i) {
               fscanf(grd,"%*[^0-9]%*d%*[^0-9]%d%*[^0-9]%*d%*[^0-9]%*d%*[^0-9]%d\n"
                  ,&count,&temp);
               
               sbdry[i] = getnewsideobject(temp,pbdrymap);
               sbdry[i]->alloc(static_cast<int>(count*grwfac));
               sbdry[i]->nel = count;
               
               fscanf(grd,"%*[^\n]\n");
               
               svrtxbtemp[i] = (int (*)[2]) xmalloc(sbdry[i]->nel*2*sizeof(int));
                    
               for(j=0;j<sbdry[i]->nel;++j) {
                  fscanf(grd,"%*d %d %d\n"
                     ,&svrtxbtemp[i][j][0],&svrtxbtemp[i][j][1]);
                  --svrtxbtemp[i][j][0];
                  --svrtxbtemp[i][j][1];
                  vd[svrtxbtemp[i][j][0]].info = sbdry[i]->idnum;
                  vd[svrtxbtemp[i][j][1]].info = sbdry[i]->idnum;
               }
            }

            /* CREATE SIDE INFORMATION */
            createsideinfo();
                
            /* FIND ALL BOUNDARY SIDES */
            /* STORE LOCATION BY VERTEX NUMBER */
            for(i=0;i<nside;++i) 
               if (sd[i].tri[1] < 0)
                  vd[sd[i].vrtx[0]].info = i;
                  
            for(i=0;i<nside;++i)
               sd[i].info = 0;

            /* MATCH BOUNDARY SIDES TO GROUPS */    
            for(i=0;i<nsbd;++i) {
               for(j=0;j<sbdry[i]->nel;++j) {
                  sind = vd[svrtxbtemp[i][j][0]].info;
                  if (sind < 0) {
                     *log << "error in boundary information " << i << j << std::endl;
                     exit(1);
                  }
                  if (sd[sind].vrtx[1] == svrtxbtemp[i][j][1]) {
                     sd[sind].info = sbdry[i]->idnum;
                     sbdry[i]->el[j] = sind;
                     vd[sd[sind].vrtx[0]].info = 0;
                  }
                  else {
                     *log << "Error: boundary sides are not counterclockwise " << 
                     svrtxbtemp[i][j][0] << svrtxbtemp[i][j][1] << std::endl;
                     exit(1);
                  }
               }
            }
            
            for(i=0;i<nvrtx;++i)
               vd[i].info = 0;

            break;
            
         case(ftype::grid):
            strcpy(grd_app,filename);
            strcat(grd_app,".grd");
            grd = fopen(grd_app,"r");
            if (grd==NULL) {
               *log << "couldn't open grid file: " << grd_app << std::endl;
               exit(1);
            }
            /* HEADER LINES */
            fscanf(grd,"nvrtx: %d\t nside: %d\t ntri: %d\n",&nvrtx,&nside,&ntri);
            
            if (!initialized) {
               allocate(nside + (int) (grwfac*nside),win);
               vrtx = vrtx;
            }
            else if (nside > maxvst) {
               *log << "mesh is too large" << std::endl;
               exit(1);
            }
   
            /* VRTX INFO */                        
            for(i=0;i<nvrtx;++i) {
               fscanf(grd,"%*d:");
               for(n=0;n<ND;++n)
               	fscanf(grd,"%lf",&vrtx[i][n]);
               fscanf(grd,"\n");
            }
                     
            /* SIDE INFO */
            for(i=0;i<nside;++i)
               fscanf(grd,"%*d: %d %d\n",&sd[i].vrtx[0],&sd[i].vrtx[1]);
   
            /* THEN TRI INFO */
            for(i=0;i<ntri;++i)
               fscanf(grd,"%*d: %d %d %d\n",&td[i].vrtx[0],&td[i].vrtx[1],&td[i].vrtx[2]);
               
            /* CREATE TSIDE & STRI */
            createtdstri();
   
            /* SIDE BOUNDARY INFO HEADER */
            fscanf(grd,"%*[^:]:%d",&nsbd);
            count = 0;
            for(i=0;i<nsbd;++i) {
               fscanf(grd,"%*[^:]:%d",&temp);
               sbdry[i] = getnewsideobject(temp,pbdrymap);
               fscanf(grd,"%*[^:]:%d\n",&sbdry[i]->nel);
               if (!sbdry[i]->maxel) sbdry[i]->alloc(static_cast<int>(grwfac*sbdry[i]->nel));
               else assert(sbdry[i]->nel < sbdry[i]->maxel);
               for(int j=0;j<sbdry[i]->nel;++j) {
                  fscanf(grd,"%*d:%d\n",&sbdry[i]->el[j]);
               }
            }
            
            /* VERTEX BOUNDARY INFO HEADER */
            fscanf(grd,"%*[^:]:%d",&nvbd);
            for(i=0;i<nvbd;++i) {
               fscanf(grd,"%*[^:]:%d",&temp);
               vbdry[i] = getnewvrtxobject(temp,pbdrymap);
               fscanf(grd,"point: %d\n",&vbdry[i]->v0);
            }
            
            break;
            
         case(ftype::mavriplis):
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
               allocate(nside + (int) (grwfac*nside),win);
               vrtx = vrtx;
            }
            else if (nside > maxvst) {
               *log << "mesh is too large" << std::endl;
               exit(1);
            }
            
            fscanf(grd,"%*[^\n]\n");
            
            fscanf(grd,"%*d%d%*d%d%*d%*d\n",&temp,&count);
            
            /* EXTERNAL BOUNDARY */
            sbdry[0] = getnewsideobject(1,pbdrymap);
            sbdry[0]->alloc(static_cast<int>(grwfac*temp));
            sbdry[0]->nel = temp;
            for(i=0;i<sbdry[0]->nel;++i)
               sbdry[0]->el[i] = i;
            
            ++nsbd;
            
            for(i=1;i<nsbd;++i) {
               sbdry[i] = getnewsideobject(i+1,pbdrymap);
               fscanf(grd,"%*[^\n]\n");
               fscanf(grd,"%d%d%*[^\n]\n",&temp,&sbdry[i]->nel);
               sbdry[i]->alloc(static_cast<int>(grwfac*sbdry[i]->nel));
               sbdry[i]->el[0] = temp-1;
               sbdry[i]->nel -= sbdry[i]->el[0];
               for(j=1;j<sbdry[i]->nel;++j)
                  sbdry[i]->el[j] = j +sbdry[i]->el[0];
            }
            
            fscanf(grd,"%*[^\n]\n");
               
            for(i=0;i<nside;++i) {
               fscanf(grd,"%d%d%d%d%*d\n",&sd[i].vrtx[0],&sd[i].vrtx[1],&sd[i].tri[0],&sd[i].tri[1]);
               --sd[i].vrtx[0];--sd[i].vrtx[1];--sd[i].tri[0];--sd[i].tri[1];
               
               if (sd[i].tri[1] >= ntri) sd[i].tri[1] = -1;
            }
            
            for(i=0;i<nvrtx;++i) {
               for(n=0;n<ND;++n)
               	fscanf(grd,"%lf",&vrtx[i][n]);
               fscanf(grd,"%*[^\n]\n");
            }
               
            for(i=0;i<ntri;++i)
               for(j=0;j<3;++j)
                  td[i].side[j] = -1;
                  
            for(i=0;i<nside;++i) {
               tind = sd[i].tri[0];
               j = 0;
               while (td[tind].side[j] > 0)
                  ++j;
               td[tind].side[j] = i;
               td[tind].sign[j] = 1;
               
               tind = sd[i].tri[1];
               if (tind > -1) {
                  j = 0;
                  while (td[tind].side[j] > 0)
                     ++j;
                  td[tind].side[j] = i;
                  td[tind].sign[j] = -1;
               }
            }
            
            /* REORDER SIDES TO BE COUNTERCLOCKWISE */
            /* FILL IN TVRTX */
            for(tind=0;tind<ntri;++tind) {
               v0 = sd[td[tind].side[0]].vrtx[(td[tind].sign[0]+1)/2];
               v1 = sd[td[tind].side[1]].vrtx[(1-td[tind].sign[1])/2];
               if (v0 != v1) {
                  /* SWITCH SIDES */
                  j = td[tind].side[1];
                  td[tind].side[1] = td[tind].side[2];
                  td[tind].side[2] = j;
                  j = td[tind].sign[1];
                  td[tind].sign[1] = td[tind].sign[2];
                  td[tind].sign[2] = j;
               }
               
               td[tind].vrtx[2] = v0;
               sind = td[tind].side[2];
               sign = td[tind].sign[2];
               td[tind].vrtx[1] = sd[sind].vrtx[(sign+1)/2];
               td[tind].vrtx[0] = sd[sind].vrtx[(1-sign)/2];
            }
            nvbd = 0;
                        
            break;
            
         case(ftype::text):
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
                  ierr = fscanf(grd,"%lf",&vrtx[i][n]);
               fscanf(grd,"\n");
               if (ierr != ND) { *log << "2c: error in grid" << std::endl; exit(1); }
            }
            fclose(grd);
            
            treeinit();
            
            return(scratch);

#ifdef CAPRI         
         case(ftype::BRep):

            /* READ VOLUME & FACE NUMBER FROM FILENAME STRING */
            sscanf(filename,"%d%d",&cpri_vol,&cpri_face);
                        
            status = gi_dGetVolume(cpri_vol,&cpri_nnode,&cpri_nedge,&cpri_nface,&cpri_nbound, &cpri_name);
            *log << " gi_uGetVolume status =" << status << std::endl;
            if (status != CAPRI_SUCCESS) exit(1);
            *log << "  # Edges = " << cpri_nedge << " # Faces = " << cpri_nface << std::endl;

            status = gi_dTesselFace(cpri_vol, cpri_face, &ntri, &cpri_tris, &cpri_tric, &nvrtx, &cpri_points, 
               &cpri_ptype, &cpri_pindex, &cpri_uv);
            if (status != CAPRI_SUCCESS) {
               *log << "gi_dTesselFace status = " << status << std::endl;
               exit(1);
            }

            /* ALLOCATE BASIC STORAGE */
            if (!initialized) {
               maxvst = static_cast<int>((grwfac*3*ntri)/2); 
               allocate(maxvst,win);
               vrtx = vrtx;
            }
            else if ((3*ntri)/2 > maxvst) {
               *log << "mesh is too large" << std::endl;
               exit(1);
            }
            
            /* LOAD VERTEX INFO */
            for (i=0;i<nvrtx;++i)
               for(n=0;n<ND;++n)
                  vrtx[i][n] = cpri_points[i*3+n];
            
            /* LOAD TRI INFORMATION */
            for (i=0;i<ntri;++i)
               for(n=0;n<3;++n)
                  td[i].vrtx[n] = cpri_tris[i*3+n] -1;
                              
            /* CREATE SIDE INFORMATION */
            createsideinfo();
            
            nvbd = 0;
            /* CREATE VERTEX BOUNDARIES */
            for(i=0;i<nvrtx;++i) {
               if (cpri_ptype[i] == 0) {
                  vbdry[nvbd] = getnewvrtxobject(cpri_pindex[i],pbdrymap);
                  vbdry[nvbd]->alloc(1);
                  vbdry[nvbd]->v0 = i;
                  ++nvbd;
                  if (nvbd >= MAXVB) {
                     *log << "Too many vertex boundary conditions: increase MAXSB: " << nvbd << std::endl;
                     exit(1);
                  }
               }
            }
            
            for(i=0;i<nside;++i)
               sd[i].info = 0;
            
            /* COUNT BOUNDARY GROUPS */
            nsbd = 0;
            for(i=0;i<nside;++i) {
               if (sd[i].tri[1] < 0) {
                  v0 = sd[i].vrtx[0];
                  /* FIGURE OUT EDGE INDEX */
                  if (cpri_ptype[v0] > 0) {
                     sd[i].info = cpri_pindex[v0];
                  }
                  else {
                     v0 = sd[i].vrtx[1];
                     if (cpri_ptype[v0] > 0) {
                        sd[i].info = cpri_pindex[v0];
                     }
                     else {
                        *log << "Error in BRep Boundary Groups\n";
                        exit(1);
                     }
                  }
                  for (j = 0; j <nsbd;++j) {
                     if (sd[i].info == (sbdry[j]->idnum&0xFFFF)) {
                        ++bcntr[j];
                        goto next1b;
                     }
                  }
                  /* NEW SIDE */
                  sbdry[nsbd] = getnewsideobject(sd[i].info,pbdrymap);
                  bcntr[nsbd++] = 1;
                  if (nsbd > MAXSB) {
                     *log << "error: too many different side boundaries: increase MAXSB" << nsbd << std::endl;
                     exit(1);
                  }
               }
next1b:      continue;
            }
            
            for(i=0;i<nsbd;++i) {
               sbdry[i]->alloc(static_cast<int>(bcntr[i]*grwfac));
               sbdry[i]->nel = 0;
            }
            
            for(i=0;i<nside;++i) {
               if (sd[i].info) {
                  for (j = 0; j <nsbd;++j) {
                     if (sd[i].info == (sbdry[j]->idnum&0xFFFF)) {
                        sbdry[j]->el[sbdry[j]->nel++] = i;
                        goto next1c;
                     }
                  }
                  printf("Big error\n");
                  exit(1);
               }
next1c:     continue;
            }
            
            break;
#endif

         case(ftype::tecplot):
            strcpy(grd_app,filename);
            strcat(grd_app,".dat");
            grd = fopen(grd_app,"r");
            if (grd==NULL) {
               *log << "couldn't open grid file: " << grd_app << std::endl;
               exit(1);
            }
                
                
//            /* SKIP FIRST ZONE? */
//            char mychar;
//            do {
//                fscanf(grd,"%*[^\n]\n");
//                mychar = fgetc(grd);
//            } while (mychar != 'Z');
            
            
            /* HEADER LINES */
            fscanf(grd,"%*[^,],%*[^,],%*[^=]=%d%*[^=]=%d\n",&nvrtx,&ntri);
            
            /* MAXVST IS APPROXIMATELY THE NUMBER OF ELEMENTS  */
            /* FOR EACH ELEMENT THERE ARE APPROXIMATELY 3/2 SIDES */
            if (!initialized) {
               maxvst = static_cast<int>(grwfac*3*ntri);
               allocate(maxvst,win);
               vrtx = vrtx;
            }
            else if ((3*ntri) > maxvst) {
               *log << "mesh is too large" << std::endl;
               exit(1);
            }
               
            /* READ VERTEX DATA */
            for(i=0;i<nvrtx;++i) {
               for(n=0;n<ND;++n)
                  fscanf(grd,"%le ",&vrtx[i][n]);
               fscanf(grd,"%*[^\n]\n");
               vd[i].info = -1;
            }

            fscanf(grd,"%*[^0-9]");

            for(i=0;i<ntri;++i) {
               fscanf(grd,"%d %d %d",&td[i].vrtx[0],&td[i].vrtx[1],&td[i].vrtx[2]);
               --td[i].vrtx[0];
               --td[i].vrtx[1];
               --td[i].vrtx[2];
               td[i].info = 0;
            }

            /* CREATE SIDE INFORMATION */
            createsideinfo();

            /* FIND ALL BOUNDARY SIDES */
            count = 0;
            for(i=0;i<nside;++i)
               if (sd[i].tri[1] < 0) ++count;

            nsbd = 0;
            sbdry[nsbd] = getnewsideobject(nsbd+1,pbdrymap);
            sbdry[0]->alloc(static_cast<int>(grwfac*count));
            sbdry[0]->nel = count;
            count = 0;
            for(i=0;i<nside;++i)
               if (sd[i].tri[1] < 0) sbdry[0]->el[count++] = i;
            ++nsbd;
            nvbd = 0;
            
            break;
            
         case(ftype::boundary):
            /* LOAD BOUNDARY INFORMATION */
            strcpy(grd_app,filename);
            strcat(grd_app,".d");
            grd = fopen(grd_app,"r");
            if (grd == NULL) { 
               printf("couldn't open %s for reading\n",grd_app);
               exit(1);
            }

            fscanf(grd,"%d\n",&nvrtx);

            maxvst = static_cast<int>(grwfac*pow(nvrtx,2));
            allocate(maxvst,win);
            qtree.allocate(vrtx,maxvst);
            initialized = 1;

            for(i=0;i<nvrtx;++i) {
               fscanf(grd,"%*d:");
               for(n=0;n<ND;++n)
                  fscanf(grd,"%lf",&vrtx[i][n]);
               fscanf(grd,"%lf%d\n",&vlngth[i],&vd[i].info);
            }

            /* COUNT VERTEX BOUNDARY GROUPS  */
            nvbd = 0;
            for(i=0;i<nvrtx;++i) {
               if (vd[i].info) {
                  /* NEW VRTX B.C. */
                  vbdry[nvbd] = getnewvrtxobject(vd[i].info,pbdrymap);
                  vbdry[nvbd]->v0 = i;
                  ++nvbd;
                  if (nvbd >= MAXVB) {
                     *log << "Too many vertex boundary conditions: increase MAXSB: " << nvbd << std::endl;
                     exit(1);
                  }
               }
            }

            fscanf(grd,"%d\n",&nside);

            for(i=0;i<nside;++i)
               fscanf(grd,"%*d:%d%d%d\n",&sd[i].vrtx[0],&sd[i].vrtx[1],&sd[i].info);
               
            /* COUNT BOUNDARY GROUPS */
            nsbd = 0;
            for(i=0;i<nside;++i) {
               if (sd[i].info) {
                  for (j = 0; j <nsbd;++j) {
                     if (sd[i].info == sbdry[j]->idnum) {
                        ++bcntr[j];
                        goto bdnext1;
                     }
                  }
                  /* NEW SIDE */
                  sbdry[nsbd] = getnewsideobject(sd[i].info,pbdrymap);
                  bcntr[nsbd++] = 1;
                  if (nsbd > MAXSB) {
                     *log << "error: too many different side boundaries: increase MAXSB" << std::endl;
                     exit(1);
                  }
               }
            bdnext1:      continue;
            }

            for(i=0;i<nsbd;++i) {
               sbdry[i]->alloc(static_cast<int>(bcntr[i]*grwfac));
               sbdry[i]->nel = 0;
            }

            for(i=0;i<nside;++i) {
               if (sd[i].info) {
                  for (j = 0; j <nsbd;++j) {
                     if (sd[i].info == sbdry[j]->idnum) {
                        sbdry[j]->el[sbdry[j]->nel++] = i;
                        goto bdnext1a;
                     }
                  }
                  printf("Big error\n");
                  exit(1);
               }
            bdnext1a:     continue;
            }
            
            for(i=0;i<nside;++i)
               i1wk[i] = i+1;
              
            ntri = 0;
            triangulate(nside);
            break;
       
         default:
            *log << "That filetype is not supported" << std::endl;
            exit(1);
   }
    
   for(i=0;i<nsbd;++i) {
      /* CREATES NEW BOUNDARY FOR DISCONNECTED SEGMENTS OF SAME TYPE */
      sbdry[i]->reorder();
      sbdry[i]->setupcoordinates();
   }
   
   bdrylabel();  // CHANGES STRI / TTRI ON BOUNDARIES TO POINT TO GROUP/ELEMENT

   createttri();
   createvtri();
   cnt_nbor();
   if (!initialized) qtree.allocate(vrtx,maxvst);
   treeinit(); 
   if (filetype != ftype::boundary) initvlngth();
   checkintegrity();

   initialized = 1;
   
   return(scratch);
}

sharedmem* mesh::allocate(int mxsize, sharedmem* wkin) {
   int i;
   
   /* SIDE INFO */
   maxvst = mxsize;
   sd = new sstruct[maxvst+1];
   ++sd; // ALLOWS US TO ACCESS SINFO[-1]

   /* VERTEX INFO */                  
   vrtx = (FLT (*)[ND]) xmalloc(maxvst*ND*sizeof(FLT));
   vlngth = new FLT[maxvst];
   vd = new vstruct[maxvst+1];
   ++vd;
   
   /* TRI INFO */ 
   td = new tstruct[maxvst+1];
   ++td;
                 
   /* ALLOCATE WORK ARRAYS USED BY ALL MESHES */
   /* ALWAYS MORE SIDES THAN TRI'S and VERTICES */ 
   if (wkin) {   
      scratch = wkin;
      if (wkin->size < maxvst*sizeof(FLT)) {
         *log << "Warning: work variable is not large enough for mesh: " << wkin->size << ' ' << maxvst*sizeof(FLT) << std::endl;
      }
   }
   else {
      scratch = new sharedmem;
      scratch->allocate(maxvst*sizeof(FLT));
   }
   load_scratch_pointers();

   if (!i1wk) {
      i1wk = new int[maxvst];
      i2wk = new int[maxvst];
      i3wk = new int[maxvst];
      maxlst = 100;
      maxsrch = 100;
   
      for(i=0;i<maxvst;++i)
         i1wk[i] = -1;
      for(i=0;i<maxvst;++i)
         i2wk[i] = -1;
      for(i=0;i<maxvst;++i)
         i3wk[i] = -1;
      
      
   }
      
   return(scratch);
}

void mesh::load_scratch_pointers() {
   /* USE BASE FOR FWK */
   fwk = static_cast<FLT *>(scratch->p);
   return;
}
