#include "mesh.h"
#include "boundaries.h"
#include <utilities.h>
#include <cstdlib>
#include <cstring>
#include <input_map.h>
#include <iostream>

Array<int,1> mesh::i1wk, mesh::i2wk, mesh::i2wk_lst1, mesh::i2wk_lst2, mesh::i2wk_lst3;
Array<FLT,1> mesh::fscr1;
int mesh::maxsrch;

void mesh::input(const std::string &filename, mesh::filetype filetype, FLT grwfac,input_map& bdrymap) {
   int i,j,n,sind,count,temp;
   std::string grd_nm, bdry_nm, grd_app;
   TinyVector<int,3> v,s,e;
   ifstream in;
   FLT fltskip;
   int intskip;

   if (filename.substr(0,7) == "${HOME}") {
      grd_nm = getenv("HOME") +filename.substr(7,filename.length());
   }
   else 
      grd_nm = filename;

    switch (filetype) {            
        case(easymesh):
            /* LOAD SIDE INFORMATION */
            grd_app = grd_nm + ".s";
            in.open(grd_app.c_str());
            if (!in.good()) {
                    *sim::log << "error: couldn't open file: " << grd_app << std::endl;
                    exit(1);
            }
            if (!(in >> nside)) { 
               *sim::log << "error: in side file: " << grd_app << std::endl;
               exit(1);
            }
           
            if (!initialized) {
               allocate(nside + (int) (grwfac*nside));
            }
            else if (nside > maxvst) {
               *sim::log << "error: mesh is too large" << std::endl;
               exit(1);
            }
           
            for(i=0;i<nside;++i) {
               in.ignore(80,':');
               in >> sd(i).vrtx(0) >> sd(i).vrtx(1) >> sd(i).tri(0) >> sd(i).tri(1) >> sd(i).info;
               if(in.fail()) {
                  *sim::log << "error: in side file " << grd_app << std::endl;
                  exit(1);
               }
            }
            
            /* Don't Trust File So Will Calculate tri myself */
            for(i=0;i<maxvst;++i) {
               sd(i).tri(0) = -1;
               sd(i).tri(1) = -1;
            }
            
            /* COUNT BOUNDARY GROUPS */
            nsbd = 0;
            for(i=0;i<nside;++i) {
               if (sd(i).info) {
                  for (j = 0; j <nsbd;++j) {
                     if (sd(i).info == i1wk(j)) {
                        ++i2wk(j);
                        goto next1;
                     }
                  }
                  /* NEW SIDE */
                  i1wk(nsbd) = sd(i).info;
                  i2wk(nsbd++) = 1;
               }
next1:      continue;
            }
            
            sbdry.resize(nsbd);
            for(i=0;i<nsbd;++i) {
               sbdry(i) = getnewsideobject(i1wk(i),bdrymap);
               sbdry(i)->alloc(static_cast<int>(i2wk(i)*grwfac));
               sbdry(i)->nel = 0;
               i1wk(i) = -1;
               i2wk(i) = -1;
            }
            
            
            for(i=0;i<nside;++i) {
               if (sd(i).info) {
                  for (j = 0; j <nsbd;++j) {
                     if (sd(i).info == sbdry(j)->idnum) {
                        sbdry(j)->el(sbdry(j)->nel++) = i;
                        goto next1a;
                     }
                  }
                  *sim::log << "Big error\n";
                  exit(1);
               }
next1a:     continue;
            }
            in.close();
               
            /* LOAD VERTEX INFORMATION               */
            grd_app = grd_nm + ".n";
            in.open(grd_app.c_str());
            if (!in.good()) { *sim::log << "trouble opening grid" << grd_app << std::endl; exit(1);}
    
            if(!(in >> nvrtx)) {
               *sim::log << "1: error in grid: " << grd_app << std::endl;
               exit(1);
            }
            
            for(i=0;i<nvrtx;++i) {
               in.ignore(80,':');
               for(n=0;n<ND;++n) {
                  in >> vrtx(i)(n);
               }
               in >> vd(i).info;
               if (in.fail())  { *sim::log << "2b: error in grid" << std::endl; exit(1); }
            }
            in.close();

            /* COUNT VERTEX BOUNDARY GROUPS  */
            nvbd = 0;
            for(i=0;i<nvrtx;++i)
               if (vd(i).info) ++nvbd;
            vbdry.resize(nvbd);
            
            nvbd = 0;
            for(i=0;i<nvrtx;++i) {
               if (vd(i).info) {
                  /* NEW VRTX B.C. */
                  vbdry(nvbd) = getnewvrtxobject(vd(i).info,bdrymap);
                  vbdry(nvbd)->alloc(1);
                  vbdry(nvbd)->v0 = i;
                  ++nvbd;
               }
            }
                        
            /* LOAD ELEMENT INFORMATION */
            grd_app = grd_nm + ".e";
            in.open(grd_app.c_str());
            if (!in.good()) {
               *sim::log << "trouble opening " << grd_app << std::endl;
               exit(1);
            }
            if(!(in >> ntri)) {
               *sim::log << "error in file " << grd_app << std::endl;
               exit(1);
            }
    
            for(i=0;i<ntri;++i) {
               in.ignore(80,':');
                in >> v(0) >> v(1) >> v(2) >> e(0) >> e(1) >> e(2) >> s(0) >> s(1) >> s(2) >> fltskip >> fltskip >> td(i).info;
    
                for (j=0;j<3;++j) {
                    td(i).vrtx(j) = v(j);
                    if(sd(s(j)).vrtx(0) == v((j+1)%3)) {
                        /* SIDE DEFINED IN SAME DIRECTION */
                        sd(s(j)).tri(0) = i;
                        td(i).side(j) = s(j);
                        td(i).sign(j) = 1;
                    }
                    else {
                        /* SIDE DEFINED IN OPPOSITE DIRECTION */
                        sd(s(j)).tri(1) = i;
                        td(i).side(j) = s(j);
                        td(i).sign(j) = -1;
                    }
                }
            }
            in.close();
            break;
        
        case(gambit):
            grd_app = grd_nm +".FDNEUT";
            in.open(grd_app.c_str());
            if (!in.good()) {
                    *sim::log << "trouble opening " << grd_nm << std::endl;
                    exit(1);
            }
            
            for(i=0;i<5;++i)
               in.ignore(160,'\n');
                
            in >> nvrtx >> i >> nsbd;
            nsbd -= 1;
            sbdry.resize(nsbd);
            /* MAXVST IS APPROXIMATELY THE NUMBER OF ELEMENTS  */
            /* FOR EACH ELEMENT THERE ARE APPROXIMATELY 3/2 SIDES */
            if (!initialized) {
               maxvst = static_cast<int>((grwfac*3*i)/2); 
               allocate(maxvst);
            }
            else if ((3*i)/2 > maxvst) {
               *sim::log << "mesh is too large" << std::endl;
               exit(1);
            }
    
            for(i=0;i<8;++i)
               in.ignore(160,'\n');

            /* READ VERTEX DATA */    
            for(i=0;i<nvrtx;++i) {
               in >> intskip;
               for(n=0;n<ND;++n)
                  in >> vrtx(i)(n);
               vd(i).info = -1;
            }
                
            for(i=0;i<3;++i)
                in.ignore(160,'\n');

            /* READ ELEMENT DATA */
            // fscanf(grd,"%*[^0-9]%*d%*[^0-9]%d%*[^\n]\n",&ntri); TEMPORARY
            in.ignore(160,'\n'); 
                               
            for(i=0;i<ntri;++i) {
               in >> intskip >> td(i).vrtx(0) >> td(i).vrtx(1) >> td(i).vrtx(2);
                --td(i).vrtx(0);
                --td(i).vrtx(1);
                --td(i).vrtx(2);
                td(i).info = 0;
            }

            /* READ BOUNDARY DATA STORE TEMPORARILY */            
            int (*svrtxbtemp[10])[2];
    
            for(i=0;i<nsbd;++i) {
               //fscanf(grd,"%*[^0-9]%*d%*[^0-9]%d%*[^0-9]%*d%*[^0-9]%*d%*[^0-9]%d\n",&count,&temp); TEMPORARY
               
               sbdry(i) = getnewsideobject(temp,bdrymap);
               sbdry(i)->alloc(static_cast<int>(count*grwfac));
               sbdry(i)->nel = count;
               
               in.ignore(160,'\n');
                              
               svrtxbtemp[i] = (int (*)[2]) xmalloc(sbdry(i)->nel*2*sizeof(int));
                    
               for(j=0;j<sbdry(i)->nel;++j) {
                  in >> intskip >> svrtxbtemp[i][j][0] >> svrtxbtemp[i][j][1];
                  --svrtxbtemp[i][j][0];
                  --svrtxbtemp[i][j][1];
                  vd(svrtxbtemp[i][j][0]).info = sbdry(i)->idnum;
                  vd(svrtxbtemp[i][j][1]).info = sbdry(i)->idnum;
               }
            }

            /* CREATE SIDE INFORMATION */
            createsideinfo();
                
            /* FIND ALL BOUNDARY SIDES */
            /* STORE LOCATION BY VERTEX NUMBER */
            for(i=0;i<nside;++i) 
               if (sd(i).tri(1) < 0)
                  vd(sd(i).vrtx(0)).info = i;
                  
            for(i=0;i<nside;++i)
               sd(i).info = 0;

            /* MATCH BOUNDARY SIDES TO GROUPS */    
            for(i=0;i<nsbd;++i) {
               for(j=0;j<sbdry(i)->nel;++j) {
                  sind = vd(svrtxbtemp[i][j][0]).info;
                  if (sind < 0) {
                     *sim::log << "error in boundary information " << i << j << std::endl;
                     exit(1);
                  }
                  if (sd(sind).vrtx(1) == svrtxbtemp[i][j][1]) {
                     sd(sind).info = sbdry(i)->idnum;
                     sbdry(i)->el(j) = sind;
                     vd(sd(sind).vrtx(0)).info = 0;
                  }
                  else {
                     *sim::log << "Error: boundary sides are not counterclockwise " << 
                     svrtxbtemp[i][j][0] << svrtxbtemp[i][j][1] << std::endl;
                     exit(1);
                  }
               }
            }
            
            for(i=0;i<nvrtx;++i)
               vd(i).info = 0;
               
            in.close();

            break;
            
         case(grid):
            grd_app = grd_nm + ".grd";
            in.open(grd_app.c_str());
            if (!in.good()) {
               *sim::log << "couldn't open grid file: " << grd_app << std::endl;
               exit(1);
            }
            /* HEADER LINES */
            in.ignore(10,':');
            in >> nvrtx;
            in.ignore(10,':');
            in >> nside;
            in.ignore(10,':');
            in >> ntri;
            
            if (!initialized) {
               allocate(nside + (int) (grwfac*nside));
               nsbd = 0;
               nvbd = 0;
            }
            else if (nside > maxvst) {
               *sim::log << "mesh is too large" << std::endl;
               exit(1);
            }
   
            /* VRTX INFO */                        
            for(i=0;i<nvrtx;++i) {
               in.ignore(80,':');
               for(n=0;n<ND;++n)
               	in >> vrtx(i)(n);
            }
                     
            /* SIDE INFO */
            for(i=0;i<nside;++i) {
               in.ignore(80,':');
               in >> sd(i).vrtx(0) >> sd(i).vrtx(1);
            }
            
            /* THEN TRI INFO */
            for(i=0;i<ntri;++i) {
               in.ignore(80,':');
               in >> td(i).vrtx(0) >> td(i).vrtx(1) >> td(i).vrtx(2);
            }
               
            /* CREATE TSIDE & STRI */
            createtdstri();
   
            /* SIDE BOUNDARY INFO HEADER */
            in.ignore(80,':');
            int newnsbd;
            in >> newnsbd;
            if (nsbd == 0) {
               nsbd = newnsbd;
               sbdry.resize(nsbd);
               sbdry = 0;
            }
            else if (nsbd != newnsbd) {
               *sim::log << "reloading incompatible meshes" << std::endl;
               exit(1);
            }
            
            count = 0;
            for(i=0;i<nsbd;++i) {
               in.ignore(80,':');
               in >> temp;
               if (!sbdry(i)) sbdry(i) = getnewsideobject(temp,bdrymap);
               in.ignore(80,':');
               in >> sbdry(i)->nel;
               if (!sbdry(i)->maxel) sbdry(i)->alloc(static_cast<int>(grwfac*sbdry(i)->nel));
               else assert(sbdry(i)->nel < sbdry(i)->maxel);
               for(int j=0;j<sbdry(i)->nel;++j) {
                  in.ignore(80,':');
                  in >> sbdry(i)->el(j);
                  in.ignore(80,'\n');
               }
            }
            
            /* VERTEX BOUNDARY INFO HEADER */
            in.ignore(80,':');
            int newnvbd;
            in >> newnvbd;
            if (nvbd == 0) {
               nvbd = newnvbd;
               vbdry.resize(nvbd);
               vbdry = 0;
            }
            else if (nvbd != newnvbd) {
               *sim::log << "re-inputting into incompatible mesh object" << std::endl;
               exit(1);
            }
            
            for(i=0;i<nvbd;++i) {
               in.ignore(80,':');
               in >> temp;
               if (!vbdry(i)) {
                  vbdry(i) = getnewvrtxobject(temp,bdrymap);
                  vbdry(i)->alloc(1);
               }
               in.ignore(80,':');
               in >> vbdry(i)->v0;
            }
            in.close();
            
            break;
                        
         case(text):
            if (!initialized) {
               *sim::log << "to read in vertex positions only must first load mesh structure" << std::endl;
               exit(1);
            }

            /* LOAD VERTEX POSITIONS               */
            grd_app = grd_nm + ".txt";
            in.open(grd_app.c_str());
            if (!in.good()) { *sim::log << "trouble opening grid" << grd_app << std::endl; exit(1);}
    
            if(!(in >> temp)) {
               *sim::log << "1: error in grid " << grd_app << std::endl;
               exit(1);
            }
            if (temp != nvrtx) {
               *sim::log << "grid doesn't match vertex list" << std::endl;
               exit(1);
            }
    
            /* ERROR %lf SHOULD BE FLT */
            for(i=0;i<nvrtx;++i) {
               in.ignore(80,':');
               for(n=0;n<ND;++n)
                  in >> vrtx(i)(n);
               if (in.fail()) { *sim::log << "2c: error in grid" << std::endl; exit(1); }
            }
            in.close();
            
            treeinit();
            
            return;

#ifdef CAPRI         
         case(BRep):

            /* READ VOLUME & FACE NUMBER FROM FILENAME STRING */
            sscanf(filename,"%d%d",&cpri_vol,&cpri_face);
            status = gi_dGetVolume(cpri_vol,&cpri_nnode,&cpri_nedge,&cpri_nface,&cpri_nbound, &cpri_name);
            *sim::log << " gi_uGetVolume status =" << status << std::endl;
            if (status != CAPRI_SUCCESS) exit(1);
            *sim::log << "  # Edges = " << cpri_nedge << " # Faces = " << cpri_nface << std::endl;

            status = gi_dTesselFace(cpri_vol, cpri_face, &ntri, &cpri_tris, &cpri_tric, &nvrtx, &cpri_points, 
               &cpri_ptype, &cpri_pindex, &cpri_uv);
            if (status != CAPRI_SUCCESS) {
               *sim::log << "gi_dTesselFace status = " << status << std::endl;
               exit(1);
            }

            /* ALLOCATE BASIC STORAGE */
            if (!initialized) {
               maxvst = static_cast<int>((grwfac*3*ntri)/2); 
               allocate(maxvst);
            }
            else if ((3*ntri)/2 > maxvst) {
               *sim::log << "mesh is too large" << std::endl;
               exit(1);
            }
            
            /* LOAD VERTEX INFO */
            for (i=0;i<nvrtx;++i)
               for(n=0;n<ND;++n)
                  vrtx(i)(n) = cpri_points[i*3+n];
            
            /* LOAD TRI INFORMATION */
            for (i=0;i<ntri;++i)
               for(n=0;n<3;++n)
                  td(i).vrtx(n) = cpri_tris[i*3+n] -1;
                              
            /* CREATE SIDE INFORMATION */
            createsideinfo();
            
            nvbd = 0;
            for(i=0;i<nvrtx;++i) 
               if (cpri_ptype[i] == 0)
                  ++nvbd;
            vbdry.resize(nvbd);
            
            nvbd = 0;
            /* CREATE VERTEX BOUNDARIES */
            for(i=0;i<nvrtx;++i) {
               if (cpri_ptype[i] == 0) {
                  vbdry(nvbd) = getnewvrtxobject(cpri_pindex[i],bdrymap);
                  vbdry(nvbd)->alloc(1);
                  vbdry(nvbd)->v0 = i;
                  ++nvbd;
                  if (nvbd >= MAXVB) {
                     *sim::log << "Too many vertex boundary conditions: increase MAXSB: " << nvbd << std::endl;
                     exit(1);
                  }
               }
            }
            
            for(i=0;i<nside;++i)
               sd(i).info = 0;
            
            /* COUNT BOUNDARY GROUPS */
            nsbd = 0;
            for(i=0;i<nside;++i) {
               if (sd(i).tri(1) < 0) {
                  v0 = sd(i).vrtx(0);
                  /* FIGURE OUT EDGE INDEX */
                  if (cpri_ptype[v0] > 0) {
                     sd(i).info = cpri_pindex[v0];
                  }
                  else {
                     v0 = sd(i).vrtx(1);
                     if (cpri_ptype[v0] > 0) {
                        sd(i).info = cpri_pindex[v0];
                     }
                     else {
                        *sim::log << "Error in BRep Boundary Groups\n";
                        exit(1);
                     }
                  }
                  for (j = 0; j <nsbd;++j) {
                     if (sd(i).info == (i1wk(j)&0xFFFF)) {
                        ++i2wk(j);
                        goto next1b;
                     }
                  }
                  /* NEW SIDE */
                  i1wk(nsbd) = sd(i).info;
                  i2wk(nsbd++) = 1;
               }
next1b:      continue;
            }
            
            sbdry.resize(nsbd);
            for(i=0;i<nsbd;++i) {
               sbdry(i) = getnewsideobject(i1wk(i),bdrymap);
               sbdry(i)->alloc(static_cast<int>(i2wk(i)*grwfac));
               sbdry(i)->nel = 0;
               i1wk(i) = -1;
               i2wk(i) = -1;
            }
            
            for(i=0;i<nside;++i) {
               if (sd(i).info) {
                  for (j = 0; j <nsbd;++j) {
                     if (sd(i).info == (sbdry(j)->idnum&0xFFFF)) {
                        sbdry(j)->el(sbdry(j)->nel++) = i;
                        goto next1c;
                     }
                  }
                  *sim::log << "Big error\n";
                  exit(1);
               }
next1c:     continue;
            }
            
            break;
#endif

         case(tecplot):
            grd_app = grd_nm +".dat";
            in.open(grd_app.c_str());
            if (!in.good()) {
               *sim::log << "couldn't open tecplot file: " << grd_app << std::endl;
               exit(1);
            }
            
            /* HEADER LINES */
            in.ignore(80,',');
            in.ignore(80,',');
            in.ignore(80,'=');
            in >> nvrtx;
            in.ignore(80,'=');
            in >> ntri;
            
            /* MAXVST IS APPROXIMATELY THE NUMBER OF ELEMENTS  */
            /* FOR EACH ELEMENT THERE ARE APPROXIMATELY 3/2 SIDES */
            if (!initialized) {
               maxvst = static_cast<int>(grwfac*3*ntri);
               allocate(maxvst);
            }
            else if ((3*ntri) > maxvst) {
               *sim::log << "mesh is too large" << std::endl;
               exit(1);
            }
               
            /* READ VERTEX DATA */
            for(i=0;i<nvrtx;++i) {
               for(n=0;n<ND;++n)
                  in >> vrtx(i)(n);
               in.ignore(80,'\n');
               vd(i).info = -1;
            }
            in.ignore(80,'\n');
            while (in.get() == '#');
               in.ignore(80,'\n');
      
            for(i=0;i<ntri;++i) {
               in >> td(i).vrtx(0) >> td(i).vrtx(1) >> td(i).vrtx(2);
               --td(i).vrtx(0);
               --td(i).vrtx(1);
               --td(i).vrtx(2);
               td(i).info = 0;
            }
            


            /* CREATE SIDE INFORMATION */
            createsideinfo();

            /* FIND ALL BOUNDARY SIDES */
            count = 0;
            for(i=0;i<nside;++i)
               if (sd(i).tri(1) < 0) ++count;

            nsbd = 1;
            sbdry.resize(1);
            sbdry(0) = getnewsideobject(1,bdrymap);
            sbdry(0)->alloc(static_cast<int>(grwfac*count));
            sbdry(0)->nel = count;
            count = 0;
            for(i=0;i<nside;++i)
               if (sd(i).tri(1) < 0) sbdry(0)->el(count++) = i;
            
            nvbd = 0;
            in.close();
            
            break;
            
         case(boundary):
            /* LOAD BOUNDARY INFORMATION */
            grd_app = grd_nm +".d";
            in.open(grd_app.c_str());
            if (!in.good()) { 
               *sim::log << "couldn't open " << grd_app << "for reading\n";
               exit(1);
            }

            in >> nvrtx;

            maxvst = static_cast<int>(grwfac*nvrtx*nvrtx);
            allocate(maxvst);
            initialized = 1;

            for(i=0;i<nvrtx;++i) {
               in.ignore(80,':');
               for(n=0;n<ND;++n)
                  in >> vrtx(i)(n);
               in >> vlngth(i) >> vd(i).info;
            }

            /* COUNT VERTEX BOUNDARY GROUPS  */
            nvbd = 0;
            for(i=0;i<nvrtx;++i) 
               if (vd(i).info)
                  ++nvbd;
            vbdry.resize(nvbd);
               
            nvbd = 0;
            for(i=0;i<nvrtx;++i) {
               if (vd(i).info) {
                  /* NEW VRTX B.C. */
                  vbdry(nvbd) = getnewvrtxobject(vd(i).info,bdrymap);
                  vbdry(nvbd)->alloc(1);
                  vbdry(nvbd)->v0 = i;
                  ++nvbd;
               }
            }

            in >> nside;

            for(i=0;i<nside;++i) {
               in.ignore(80,':');
               in >> sd(i).vrtx(0) >> sd(i).vrtx(1) >> sd(i).info;
            }
               
            /* COUNT BOUNDARY GROUPS */
            nsbd = 0;
            for(i=0;i<nside;++i) {
               if (sd(i).info) {
                  for (j = 0; j <nsbd;++j) {
                     if (sd(i).info == i1wk(j)) {
                        ++i2wk(j);
                        goto bdnext1;
                     }
                  }
                  /* NEW SIDE */
                  i1wk(nsbd) = sd(i).info;
                  i2wk(nsbd++) = 1;
               }
            bdnext1:      continue;
            }


            sbdry.resize(nsbd);
            for(i=0;i<nsbd;++i) {
               sbdry(i) = getnewsideobject(i1wk(i),bdrymap);
               sbdry(i)->alloc(static_cast<int>(i2wk(i)*grwfac));
               sbdry(i)->nel = 0;
               i1wk(i) = -1;
               i2wk(i) = -1;
            }

            for(i=0;i<nside;++i) {
               if (sd(i).info) {
                  for (j = 0; j <nsbd;++j) {
                     if (sd(i).info == sbdry(j)->idnum) {
                        sbdry(j)->el(sbdry(j)->nel++) = i;
                        goto bdnext1a;
                     }
                  }
                  *sim::log << "Big error\n";
                  exit(1);
               }
            bdnext1a:     continue;
            }
            
            for(i=0;i<nside;++i)
               i2wk_lst1(i) = i+1;
              
            ntri = 0;
            triangulate(nside);
            
            in.close();
            
            break;
       
         default:
            *sim::log << "That filetype is not supported" << std::endl;
            exit(1);
   }
    
   for(i=0;i<nsbd;++i) {
      /* CREATES NEW BOUNDARY FOR DISCONNECTED SEGMENTS OF SAME TYPE */
      sbdry(i)->reorder();
   }
   
   /* FIND ENDPOINT MATCHES */
   for(i=0;i<nvbd;++i) {
      /* Find two connecting boundary sides */
      for(j=0;j<nsbd;++j) {
         if (sd(sbdry(j)->el(0)).vrtx(0) == vbdry(i)->v0) {
            vbdry(i)->sbdry(1) = j;
            sbdry(j)->vbdry(0) = i;
         }
         if (sd(sbdry(j)->el(sbdry(j)->nel-1)).vrtx(1) == vbdry(i)->v0) {
            vbdry(i)->sbdry(0) = j;
            sbdry(j)->vbdry(1) = i;
         }
      }
   }


   
   bdrylabel();  // CHANGES STRI / TTRI ON BOUNDARIES TO POINT TO GROUP/ELEMENT
   

   createttri();
   createvtri();
   cnt_nbor();
   treeinit(); 
   
   grd_app = grd_nm + ".vlngth";
   in.open(grd_app.c_str());
   if (in) {
      for(i=0;i<nvrtx;++i) in >> vlngth(i);
      in.close();
   }
   else if (filetype != boundary) initvlngth();
   
   checkintegrity();

   initialized = 1;
   
   return;
}

void mesh::allocate(int mxsize) {
   
   /* SIDE INFO */
   maxvst = mxsize;
   sd.resize(Range(-1,maxvst));

   /* VERTEX INFO */                  
   vrtx.resize(maxvst);
   vlngth.resize(maxvst);
   vd.resize(Range(-1,maxvst));
   
   /* TRI INFO */ 
   td.resize(Range(-1,maxvst));
   
   qtree.allocate((FLT (*)[ND]) vrtx(0).data(), maxvst);

                 
   /* ALLOCATE WORK ARRAYS USED BY ALL MESHES */
   /* ALWAYS MORE SIDES THAN TRI'S and VERTICES */ 
#ifdef USE_SIMSCRATCH
   if (sim::scratch.size() < needed_scratch_size()) {
      sim::scratch.resize(needed_scratch_size());
   }
#endif
   reload_scratch_pointers();

   if (i1wk.ubound(firstDim) < maxvst) {
      // i1wk should always be kept initialized to -1
      i1wk.resize(Range(-1,maxvst));
      i1wk = -1;
      maxsrch = 100;
   }
   initialized = 1;
   
   return;
}

void mesh::allocate_duplicate(FLT sizereduce1d,const class mesh& inmesh) {
   int i;

   if (!initialized) {
      maxvst =  MAX((int) (inmesh.maxvst/(sizereduce1d*sizereduce1d)),10);
      allocate(maxvst);
      nsbd = inmesh.nsbd;
      sbdry.resize(nsbd);
      for(i=0;i<nsbd;++i) {
         sbdry(i) = inmesh.sbdry(i)->create(*this);
         sbdry(i)->alloc(MAX(static_cast<int>(inmesh.sbdry(i)->maxel/sizereduce1d),10));
      }
      nvbd = inmesh.nvbd;
      vbdry.resize(nvbd);
      for(i=0;i<nvbd;++i) {
         vbdry(i) = inmesh.vbdry(i)->create(*this);
         vbdry(i)->alloc(1);
      }
      qtree.allocate((FLT (*)[ND]) vrtx(0).data(), maxvst);
      initialized = 1;
   }
}

void mesh::reload_scratch_pointers() {
   
#ifdef USE_SIMSCRATCH
   /* THIS PLACES IT IN GLOBALLY SHARED MEMORY (RISKY BECAUSE CAN OVERLAP) */
   if (fscr1.data() != (FLT *) sim::scratch.data() || maxvst > fscr1.extent(firstDim)) {

      /* THIS IS PLACEMENT NEW. IT STICKS THE OBJECT AT THE ADDRESS IN MEMORY */
      FLT *base = new(sim::scratch.data()) FLT;
      Array<FLT,1> temp(base, maxvst, neverDeleteData);
      fscr1.reference(temp);
      
      /* TEMPORARY NOT SURE IF THIS GUARANTEES ALIGNMENT??? */
      int *ibase = new(sim::scratch.data()+sizeof(FLT)*maxvst) int;
      Array<int,1> temp0(ibase, maxvst+1, neverDeleteData);
      i2wk.reference(temp0);
      i2wk.reindexSelf(TinyVector<int,1>(-1));
#else
   /* THIS IS A MORE STANDARD STATIC ALLOCATION */;
   if (maxvst > fscr1.extent(firstDim)) {
      fscr1.resize(maxvst);
      i2wk.resize(maxvst+1);
      i2wk.reindexSelf(TinyVector<int,1>(-1));
#endif
      // some smaller lists using i2 storage
      int mvst3 = maxvst/3;
      Array<int,1> temp1(i2wk.data(),mvst3,neverDeleteData);
      i2wk_lst1.reference(temp1);
      i2wk_lst1.reindexSelf(TinyVector<int,1>(-1));
      Array<int,1> temp2(i2wk.data()+1+mvst3,mvst3-1,neverDeleteData);
      i2wk_lst2.reference(temp2);
      i2wk_lst2.reindexSelf(TinyVector<int,1>(-1));
      Array<int,1> temp3(i2wk.data()+1+2*mvst3,mvst3-1,neverDeleteData);
      i2wk_lst3.reference(temp3);
      i2wk_lst3.reindexSelf(TinyVector<int,1>(-1));
   }

   return;
}

size_t mesh::needed_scratch_size() {
   return(maxvst*(sizeof(FLT)+sizeof(int)));
}

mesh::~mesh() {
   for(int i=0;i<nvbd;++i)
      delete vbdry(i);
   for(int i=0;i<nsbd;++i)
      delete sbdry(i);
}
   

