#include "tri_mesh.h"
#include "tri_boundary.h"
#include <utilities.h>
#include <cstdlib>
#include <cstring>
#include <input_map.h>
#include <iostream>

/* NAMING CONVENTION */
// pt, segment, triangle, tetrahedral 

/* Boundaries */
// vertex, edge, face

void tri_mesh::init(input_map &input, void *gin) {
    std::string keyword;
    std::istringstream data;
    std::string filename;
    std::string bdryfile;

    
    if (gin != 0) {
        gbl = static_cast<global *>(gin);   
    }
    else {
        /* gbl has not been set */
        /* so create internal gbl_struct */
        gbl = new global;
        gbl->idprefix = "";
        gbl->log = &std::cout;
        gbl->adapt_flag = true;
        gbl->tolerance = 1.25;
    }
    
    if (!initialized) {
        FLT grwfac;
        keyword = gbl->idprefix + "_growth factor";
        if (!input.get(keyword,grwfac)) {
            input.getwdefault("growth factor",grwfac,2.0);
        }
        
        int filetype;
        keyword = gbl->idprefix + "_filetype";
        if (!input.get(keyword,filetype)) {
            input.getwdefault("filetype",filetype,static_cast<int>(tri_mesh::grid));
        }
        
        keyword = gbl->idprefix + "_mesh";
        if (!input.get(keyword,filename)) {
            if (input.get("mesh",filename)) {
                filename = filename +"_" +gbl->idprefix;
            }
            else {
                *gbl->log << "no mesh name\n";
                exit(1);
            }
        }
        coarse_level = 0;
        tri_mesh::input(filename.c_str(),static_cast<tri_mesh::filetype>(filetype),grwfac,input);
    }
}

void tri_mesh::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
    int i;
   
     if (!initialized) {
        const tri_mesh& inmesh = dynamic_cast<const tri_mesh &>(in);
        gbl = inmesh.gbl;
        maxpst =  MAX((int) (inmesh.maxpst/(sizereduce1d*sizereduce1d)),10);
        allocate(maxpst);
        nebd = inmesh.nebd;
        ebdry.resize(nebd);
        for(i=0;i<nebd;++i) {
            ebdry(i) = inmesh.ebdry(i)->create(*this);
            ebdry(i)->alloc(MAX(static_cast<int>(inmesh.ebdry(i)->maxseg/sizereduce1d),10));
        }
        nvbd = inmesh.nvbd;
        vbdry.resize(nvbd);
        for(i=0;i<nvbd;++i) {
            vbdry(i) = inmesh.vbdry(i)->create(*this);
            vbdry(i)->alloc(4);
        }
        qtree.allocate((FLT (*)[ND]) pnts(0).data(), maxpst);
        if (why == multigrid) coarse_level = inmesh.coarse_level+1;
        initialized = 1;
    }
}

void tri_mesh::allocate(int mxsize) {
    
    /* SIDE INFO */
    maxpst = mxsize;
    seg.resize(Range(-1,maxpst));

    /* VERTEX INFO */                        
    pnts.resize(maxpst);
    lngth.resize(maxpst);
    pnt.resize(Range(-1,maxpst));
    
    /* TRI INFO */ 
    tri.resize(Range(-1,maxpst));

    qtree.allocate((FLT (*)[ND]) pnts(0).data(), maxpst);

    if (!gbl) {
        /* gbl has not been set */
        /* so create internal gbl_struct */
        gbl = new global;
        gbl->idprefix = "";
        gbl->log = &std::cout;
        gbl->adapt_flag = true;
        gbl->tolerance = 1.25;
    }
    
    if (gbl->intwk.ubound(firstDim) < maxpst) {
        // gbl->intwk should always be kept initialized to -1
        gbl->intwk.resize(Range(-1,maxpst));
        gbl->intwk = -1;
        gbl->maxsrch = 1000;
        
        gbl->fltwk.resize(maxpst);
        gbl->i2wk.resize(maxpst+1);
        gbl->i2wk.reindexSelf(TinyVector<int,1>(-1));
        // some smaller lists using i2 storage
        int mvst3 = maxpst/3;
        Array<int,1> temp1(gbl->i2wk.data(),mvst3,neverDeleteData);
        gbl->i2wk_lst1.reference(temp1);
        gbl->i2wk_lst1.reindexSelf(TinyVector<int,1>(-1));
        Array<int,1> temp2(gbl->i2wk.data()+1+mvst3,mvst3-1,neverDeleteData);
        gbl->i2wk_lst2.reference(temp2);
        gbl->i2wk_lst2.reindexSelf(TinyVector<int,1>(-1));
        Array<int,1> temp3(gbl->i2wk.data()+1+2*mvst3,mvst3-1,neverDeleteData);
        gbl->i2wk_lst3.reference(temp3);
        gbl->i2wk_lst3.reindexSelf(TinyVector<int,1>(-1));
    }
    
    initialized = 1;
    
    return;
}


void tri_mesh::input(const std::string &filename, tri_mesh::filetype filetype, FLT grwfac,input_map& bdrymap) {
    int i,j,n,sind,count,temp;
    std::string grd_nm, bdry_nm, grd_app;
    TinyVector<int,3> v,s,e;
    ifstream in;
    FLT fltskip;
    int intskip;
	Array<Array<TinyVector<int,2>,1>,1> svrtxbtemp;  // TEMPORARY FOR LOADING GAMBIT

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
                          *gbl->log << "error: couldn't open file: " << grd_app << std::endl;
                          exit(1);
                }
                if (!(in >> nseg)) { 
                    *gbl->log << "error: in segment file: " << grd_app << std::endl;
                    exit(1);
                }
              
                if (!initialized) {
                    allocate(nseg + (int) (grwfac*nseg));
                }
                else if (nseg > maxpst) {
                    *gbl->log << "error: mesh is too large" << std::endl;
                    exit(1);
                }
              
                for(i=0;i<nseg;++i) {
                    in.ignore(80,':');
                    in >> seg(i).pnt(0) >> seg(i).pnt(1) >> seg(i).tri(0) >> seg(i).tri(1) >> seg(i).info;
                    if(in.fail()) {
                        *gbl->log << "error: in segment file " << grd_app << std::endl;
                        exit(1);
                    }
                }
                
                /* Don't Trust File So Will Calculate tri myself */
                for(i=0;i<maxpst;++i) {
                    seg(i).tri(0) = -1;
                    seg(i).tri(1) = -1;
                }
                
                /* COUNT BOUNDARY GROUPS */
                nebd = 0;
                for(i=0;i<nseg;++i) {
                    if (seg(i).info) {
                        for (j = 0; j <nebd;++j) {
                            if (seg(i).info == gbl->intwk(j)) {
                                ++gbl->i2wk(j);
                                goto next1;
                            }
                        }
                        /* NEW SIDE */
                        gbl->intwk(nebd) = seg(i).info;
                        gbl->i2wk(nebd++) = 1;
                    }
next1:        continue;
                }
                
                ebdry.resize(nebd);
                for(i=0;i<nebd;++i) {
                    ebdry(i) = getnewedgeobject(gbl->intwk(i),bdrymap);
                    ebdry(i)->alloc(static_cast<int>(gbl->i2wk(i)*grwfac));
                    ebdry(i)->nseg = 0;
                    gbl->intwk(i) = -1;
                    gbl->i2wk(i) = -1;
                }
                
                
                for(i=0;i<nseg;++i) {
                    if (seg(i).info) {
                        for (j = 0; j <nebd;++j) {
                            if (seg(i).info == ebdry(j)->idnum) {
                                ebdry(j)->seg(ebdry(j)->nseg++) = i;
                                goto next1a;
                            }
                        }
                        *gbl->log << "Big error\n";
                        exit(1);
                    }
next1a:      continue;
                }
                in.close();
                    
                /* LOAD VERTEX INFORMATION                    */
                grd_app = grd_nm + ".n";
                in.open(grd_app.c_str());
                if (!in.good()) { *gbl->log << "trouble opening grid" << grd_app << std::endl; exit(1);}
     
                if(!(in >> npnt)) {
                    *gbl->log << "1: error in grid: " << grd_app << std::endl;
                    exit(1);
                }
                
                for(i=0;i<npnt;++i) {
                    in.ignore(80,':');
                    for(n=0;n<ND;++n) {
                        in >> pnts(i)(n);
                    }
                    in >> pnt(i).info;
                    if (in.fail())  { *gbl->log << "2b: error in grid" << std::endl; exit(1); }
                }
                in.close();

                /* COUNT VERTEX BOUNDARY GROUPS  */
                nvbd = 0;
                for(i=0;i<npnt;++i)
                    if (pnt(i).info) ++nvbd;
                vbdry.resize(nvbd);
                
                nvbd = 0;
                for(i=0;i<npnt;++i) {
                    if (pnt(i).info) {
                        /* NEW VRTX B.C. */
                        vbdry(nvbd) = getnewvrtxobject(pnt(i).info,bdrymap);
                        vbdry(nvbd)->alloc(4);
                        vbdry(nvbd)->pnt = i;
                        ++nvbd;
                    }
                }
                                
                /* LOAD ELEMENT INFORMATION */
                grd_app = grd_nm + ".e";
                in.open(grd_app.c_str());
                if (!in.good()) {
                    *gbl->log << "trouble opening " << grd_app << std::endl;
                    exit(1);
                }
                if(!(in >> ntri)) {
                    *gbl->log << "error in file " << grd_app << std::endl;
                    exit(1);
                }
     
                for(i=0;i<ntri;++i) {
                    in.ignore(80,':');
                     in >> v(0) >> v(1) >> v(2) >> e(0) >> e(1) >> e(2) >> s(0) >> s(1) >> s(2) >> fltskip >> fltskip >> tri(i).info;
     
                     for (j=0;j<3;++j) {
                          tri(i).pnt(j) = v(j);
                          if(seg(s(j)).pnt(0) == v((j+1)%3)) {
                                /* SIDE DEFINED IN SAME DIRECTION */
                                seg(s(j)).tri(0) = i;
                                tri(i).seg(j) = s(j);
                                tri(i).sgn(j) = 1;
                          }
                          else {
                                /* SIDE DEFINED IN OPPOSITE DIRECTION */
                                seg(s(j)).tri(1) = i;
                                tri(i).seg(j) = s(j);
                                tri(i).sgn(j) = -1;
                          }
                     }
                }
                in.close();
                break;
          
          case(gambit):
                grd_app = grd_nm +".FDNEUT";
                in.open(grd_app.c_str());
                if (!in.good()) {
                          *gbl->log << "trouble opening " << grd_nm << std::endl;
                          exit(1);
                }
                
                for(i=0;i<5;++i)
                    in.ignore(160,'\n');
                     
                in >> npnt >> i >> nebd;
                nebd -= 1;
                ebdry.resize(nebd);
                /* MAXVST IS APPROXIMATELY THE NUMBER OF ELEMENTS  */
                /* FOR EACH ELEMENT THERE ARE APPROXIMATELY 3/2 SIDES */
                if (!initialized) {
                    maxpst = static_cast<int>((grwfac*3*i)/2); 
                    allocate(maxpst);
                }
                else if ((3*i)/2 > maxpst) {
                    *gbl->log << "mesh is too large" << std::endl;
                    exit(1);
                }
     
                for(i=0;i<8;++i)
                    in.ignore(160,'\n');

                /* READ VERTEX DATA */     
                for(i=0;i<npnt;++i) {
                    in >> intskip;
                    for(n=0;n<ND;++n)
                        in >> pnts(i)(n);
                    pnt(i).info = -1;
                }
                     
                for(i=0;i<3;++i)
                     in.ignore(160,'\n');

                /* READ ELEMENT DATA */
                // fscanf(grd,"%*[^0-9]%*d%*[^0-9]%d%*[^\n]\n",&ntri); TEMPORARY
                in.ignore(160,'\n'); 
                                         
                for(i=0;i<ntri;++i) {
                    in >> intskip >> tri(i).pnt(0) >> tri(i).pnt(1) >> tri(i).pnt(2);
                     --tri(i).pnt(0);
                     --tri(i).pnt(1);
                     --tri(i).pnt(2);
                     tri(i).info = 0;
                }

                /* READ BOUNDARY DATA STORE TEMPORARILY */                
				svrtxbtemp.resize(10);
     
                for(i=0;i<nebd;++i) {
                    //fscanf(grd,"%*[^0-9]%*d%*[^0-9]%d%*[^0-9]%*d%*[^0-9]%*d%*[^0-9]%d\n",&count,&temp); TEMPORARY
                    
                    ebdry(i) = getnewedgeobject(temp,bdrymap);
                    ebdry(i)->alloc(static_cast<int>(count*grwfac));
                    ebdry(i)->nseg = count;
                    
                    in.ignore(160,'\n');
                                        
					svrtxbtemp(i).resize(ebdry(i)->nseg);
                          
                    for(j=0;j<ebdry(i)->nseg;++j) {
                        in >> intskip >> svrtxbtemp(i)(j)(0) >> svrtxbtemp(i)(j)(1);
                        --svrtxbtemp(i)(j)(0);
                        --svrtxbtemp(i)(j)(1);
                        pnt(svrtxbtemp(i)(j)(0)).info = ebdry(i)->idnum;
                        pnt(svrtxbtemp(i)(j)(1)).info = ebdry(i)->idnum;
                    }
                }

                /* CREATE SIDE INFORMATION */
                createseg();
                     
                /* FIND ALL BOUNDARY SIDES */
                /* STORE LOCATION BY VERTEX NUMBER */
                for(i=0;i<nseg;++i) 
                    if (seg(i).tri(1) < 0)
                        pnt(seg(i).pnt(0)).info = i;
                        
                for(i=0;i<nseg;++i)
                    seg(i).info = 0;

                /* MATCH BOUNDARY SIDES TO GROUPS */     
                for(i=0;i<nebd;++i) {
                    for(j=0;j<ebdry(i)->nseg;++j) {
                        sind = pnt(svrtxbtemp(i)(j)(0)).info;
                        if (sind < 0) {
                            *gbl->log << "error in boundary information " << i << j << std::endl;
                            exit(1);
                        }
                        if (seg(sind).pnt(1) == svrtxbtemp(i)(j)(1)) {
                            seg(sind).info = ebdry(i)->idnum;
                            ebdry(i)->seg(j) = sind;
                            pnt(seg(sind).pnt(0)).info = 0;
                        }
                        else {
                            *gbl->log << "Error: boundary sides are not counterclockwise " << 
                            svrtxbtemp(i)(j)(0) << svrtxbtemp(i)(j)(1) << std::endl;
                            exit(1);
                        }
                    }
                }
                
                for(i=0;i<npnt;++i)
                    pnt(i).info = 0;
                    
                in.close();
				~svrtxbtemp;

                break;
                
            case(grid):
                grd_app = grd_nm + ".grd";
                in.open(grd_app.c_str());
                if (!in.good()) {
                    *gbl->log << "couldn't open grid file: " << grd_app << std::endl;
                    exit(1);
                }
                /* HEADER LINES */
                in.ignore(10,':');
                in >> npnt;
                in.ignore(10,':');
                in >> nseg;
                in.ignore(10,':');
                in >> ntri;
                
                if (!initialized) {
                    allocate(nseg + (int) (grwfac*nseg));
                    nebd = 0;
                    nvbd = 0;
                }
                else if (nseg > maxpst) {
                    *gbl->log << "mesh is too large" << std::endl;
                    exit(1);
                }
    
                /* VRTX INFO */                                
                for(i=0;i<npnt;++i) {
                    in.ignore(80,':');
                    for(n=0;n<ND;++n)
                    	in >> pnts(i)(n);
                }
                            
                /* SIDE INFO */
                for(i=0;i<nseg;++i) {
                    in.ignore(80,':');
                    in >> seg(i).pnt(0) >> seg(i).pnt(1);
                }
                
                /* THEN TRI INFO */
                for(i=0;i<ntri;++i) {
                    in.ignore(80,':');
                    in >> tri(i).pnt(0) >> tri(i).pnt(1) >> tri(i).pnt(2);
                }
                    
                /* CREATE TSIDE & STRI */
                createsegtri();
    
                /* SIDE BOUNDARY INFO HEADER */
                in.ignore(80,':');
                int newnsbd;
                in >> newnsbd;
                if (nebd == 0) {
                    nebd = newnsbd;
                    ebdry.resize(nebd);
                    ebdry = 0;
                }
                else if (nebd != newnsbd) {
                    *gbl->log << "reloading incompatible meshes" << std::endl;
                    exit(1);
                }
                
                count = 0;
                for(i=0;i<nebd;++i) {
                    in.ignore(80,':');
                    in >> temp;
                    if (!ebdry(i)) ebdry(i) = getnewedgeobject(temp,bdrymap);
                    in.ignore(80,':');
                    in >> ebdry(i)->nseg;
                    if (!ebdry(i)->maxseg) ebdry(i)->alloc(static_cast<int>(grwfac*ebdry(i)->nseg));
                    else assert(ebdry(i)->nseg < ebdry(i)->maxseg);
                    ebdry(i)->input(in,grid);
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
                    *gbl->log << "re-inputting into incompatible mesh object" << std::endl;
                    exit(1);
                }
                
                for(i=0;i<nvbd;++i) {
                    in.ignore(80,':');
                    in >> temp;
                    if (!vbdry(i)) {
                        vbdry(i) = getnewvrtxobject(temp,bdrymap);
                        vbdry(i)->alloc(4);
                    }
                    in.ignore(80,':');
                    in >> vbdry(i)->pnt;
                }
                in.close();
                
                break;
                                
            case(text):
                if (!initialized) {
                    *gbl->log << "to read in point positions only must first load mesh structure" << std::endl;
                    exit(1);
                }

                /* LOAD VERTEX POSITIONS                    */
                grd_app = grd_nm + ".txt";
                in.open(grd_app.c_str());
                if (!in.good()) { *gbl->log << "trouble opening grid" << grd_app << std::endl; exit(1);}
     
                if(!(in >> temp)) {
                    *gbl->log << "1: error in grid " << grd_app << std::endl;
                    exit(1);
                }
                if (temp != npnt) {
                    *gbl->log << "grid doesn't match point list" << std::endl;
                    exit(1);
                }
     
                /* ERROR %lf SHOULD BE FLT */
                for(i=0;i<npnt;++i) {
                    in.ignore(80,':');
                    for(n=0;n<ND;++n)
                        in >> pnts(i)(n);
                    if (in.fail()) { *gbl->log << "2c: error in grid" << std::endl; exit(1); }
                }
                in.close();
                
                treeinit();
                
                return;

#ifdef CAPRI            
            case(BRep):

                /* READ VOLUME & FACE NUMBER FROM FILENAME STRING */
                sscanf(filename,"%d%d",&cpri_vol,&cpri_face);
                status = gi_dGetVolume(cpri_vol,&cpri_nnode,&cpri_nedge,&cpri_nface,&cpri_nbound, &cpri_name);
                *gbl->log << " gi_uGetVolume status =" << status << std::endl;
                if (status != CAPRI_SUCCESS) exit(1);
                *gbl->log << "  # Edges = " << cpri_nedge << " # Faces = " << cpri_nface << std::endl;

                status = gi_dTesselFace(cpri_vol, cpri_face, &ntri, &cpri_tris, &cpri_tric, &npnt, &cpri_points, 
                    &cpri_ptype, &cpri_pindex, &cpri_uv);
                if (status != CAPRI_SUCCESS) {
                    *gbl->log << "gi_dTesselFace status = " << status << std::endl;
                    exit(1);
                }

                /* ALLOCATE BASIC STORAGE */
                if (!initialized) {
                    maxpst = static_cast<int>((grwfac*3*ntri)/2); 
                    allocate(maxpst);
                }
                else if ((3*ntri)/2 > maxpst) {
                    *gbl->log << "mesh is too large" << std::endl;
                    exit(1);
                }
                
                /* LOAD VERTEX INFO */
                for (i=0;i<npnt;++i)
                    for(n=0;n<ND;++n)
                        pnts(i)(n) = cpri_points[i*3+n];
                
                /* LOAD TRI INFORMATION */
                for (i=0;i<ntri;++i)
                    for(n=0;n<3;++n)
                        tri(i).pnt(n) = cpri_tris[i*3+n] -1;
                                        
                /* CREATE SIDE INFORMATION */
                createseg();
                
                nvbd = 0;
                for(i=0;i<npnt;++i) 
                    if (cpri_ptype[i] == 0)
                        ++nvbd;
                vbdry.resize(nvbd);
                
                nvbd = 0;
                /* CREATE VERTEX BOUNDARIES */
                for(i=0;i<npnt;++i) {
                    if (cpri_ptype[i] == 0) {
                        vbdry(nvbd) = getnewvrtxobject(cpri_pindex[i],bdrymap);
                        vbdry(nvbd)->alloc(4);
                        vbdry(nvbd)->pnt = i;
                        ++nvbd;
                        if (nvbd >= MAXVB) {
                            *gbl->log << "Too many vertex boundary conditions: increase MAXSB: " << nvbd << std::endl;
                            exit(1);
                        }
                    }
                }
                
                for(i=0;i<nseg;++i)
                    seg(i).info = 0;
                
                /* COUNT BOUNDARY GROUPS */
                nebd = 0;
                for(i=0;i<nseg;++i) {
                    if (seg(i).tri(1) < 0) {
                        p0 = seg(i).pnt(0);
                        /* FIGURE OUT EDGE INDEX */
                        if (cpri_ptype[p0] > 0) {
                            seg(i).info = cpri_pindex[p0];
                        }
                        else {
                            p0 = seg(i).pnt(1);
                            if (cpri_ptype[p0] > 0) {
                                seg(i).info = cpri_pindex[p0];
                            }
                            else {
                                *gbl->log << "Error in BRep Boundary Groups\n";
                                exit(1);
                            }
                        }
                        for (j = 0; j <nebd;++j) {
                            if (seg(i).info == (gbl->intwk(j)&0xFFFF)) {
                                ++gbl->i2wk(j);
                                goto next1b;
                            }
                        }
                        /* NEW SIDE */
                        gbl->intwk(nebd) = seg(i).info;
                        gbl->i2wk(nebd++) = 1;
                    }
next1b:        continue;
                }
                
                ebdry.resize(nebd);
                for(i=0;i<nebd;++i) {
                    ebdry(i) = getnewedgeobject(gbl->intwk(i),bdrymap);
                    ebdry(i)->alloc(static_cast<int>(gbl->i2wk(i)*grwfac));
                    ebdry(i)->nseg = 0;
                    gbl->intwk(i) = -1;
                    gbl->i2wk(i) = -1;
                }
                
                for(i=0;i<nseg;++i) {
                    if (seg(i).info) {
                        for (j = 0; j <nebd;++j) {
                            if (seg(i).info == (ebdry(j)->idnum&0xFFFF)) {
                                ebdry(j)->seg(ebdry(j)->nseg++) = i;
                                goto next1c;
                            }
                        }
                        *gbl->log << "Big error\n";
                        exit(1);
                    }
next1c:      continue;
                }
                
                break;
#endif

            case(tecplot):
                grd_app = grd_nm +".dat";
                in.open(grd_app.c_str());
                if (!in.good()) {
                    *gbl->log << "couldn't open tecplot file: " << grd_app << std::endl;
                    exit(1);
                }
                
                /* HEADER LINES */
                in.ignore(80,',');
                in.ignore(80,',');
                in.ignore(80,'=');
                in >> npnt;
                in.ignore(80,'=');
                in >> ntri;
                
                /* MAXVST IS APPROXIMATELY THE NUMBER OF ELEMENTS  */
                /* FOR EACH ELEMENT THERE ARE APPROXIMATELY 3/2 SIDES */
                if (!initialized) {
                    maxpst = static_cast<int>(grwfac*3*ntri);
                    allocate(maxpst);
                }
                else if ((3*ntri) > maxpst) {
                    *gbl->log << "mesh is too large" << std::endl;
                    exit(1);
                }
                    
                /* READ VERTEX DATA */
                for(i=0;i<npnt;++i) {
                    for(n=0;n<ND;++n)
                        in >> pnts(i)(n);
                    in.ignore(80,'\n');
                    pnt(i).info = -1;
                }
                in.ignore(80,'\n');
                while (in.get() == '#');
                    in.ignore(80,'\n');
        
                for(i=0;i<ntri;++i) {
                    in >> tri(i).pnt(0) >> tri(i).pnt(1) >> tri(i).pnt(2);
                    --tri(i).pnt(0);
                    --tri(i).pnt(1);
                    --tri(i).pnt(2);
                    tri(i).info = 0;
                }
                


                /* CREATE SIDE INFORMATION */
                createseg();

                /* FIND ALL BOUNDARY SIDES */
                count = 0;
                for(i=0;i<nseg;++i)
                    if (seg(i).tri(1) < 0) ++count;

                nebd = 1;
                ebdry.resize(1);
                ebdry(0) = getnewedgeobject(1,bdrymap);
                ebdry(0)->alloc(static_cast<int>(grwfac*count));
                ebdry(0)->nseg = count;
                count = 0;
                for(i=0;i<nseg;++i)
                    if (seg(i).tri(1) < 0) ebdry(0)->seg(count++) = i;
                
                nvbd = 0;
                in.close();
                
                break;
                
            case(boundary):
                /* LOAD BOUNDARY INFORMATION */
                grd_app = grd_nm +".d";
                in.open(grd_app.c_str());
                if (!in.good()) { 
                    *gbl->log << "couldn't open " << grd_app << "for reading\n";
                    exit(1);
                }

                in >> npnt;

                maxpst = static_cast<int>(grwfac*npnt*npnt);
                allocate(maxpst);
                initialized = 1;

                for(i=0;i<npnt;++i) {
                    in.ignore(80,':');
                    for(n=0;n<ND;++n)
                        in >> pnts(i)(n);
                    in >> lngth(i) >> pnt(i).info;
                }

                /* COUNT VERTEX BOUNDARY GROUPS  */
                nvbd = 0;
                for(i=0;i<npnt;++i) 
                    if (pnt(i).info)
                        ++nvbd;
                vbdry.resize(nvbd);
                    
                nvbd = 0;
                for(i=0;i<npnt;++i) {
                    if (pnt(i).info) {
                        /* NEW VRTX B.C. */
                        vbdry(nvbd) = getnewvrtxobject(pnt(i).info,bdrymap);
                        vbdry(nvbd)->alloc(4);
                        vbdry(nvbd)->pnt = i;
                        ++nvbd;
                    }
                }

                in >> nseg;

                for(i=0;i<nseg;++i) {
                    in.ignore(80,':');
                    in >> seg(i).pnt(0) >> seg(i).pnt(1) >> seg(i).info;
                }
                    
                /* COUNT BOUNDARY GROUPS */
                nebd = 0;
                for(i=0;i<nseg;++i) {
                    if (seg(i).info) {
                        for (j = 0; j <nebd;++j) {
                            if (seg(i).info == gbl->intwk(j)) {
                                ++gbl->i2wk(j);
                                goto bdnext1;
                            }
                        }
                        /* NEW SIDE */
                        gbl->intwk(nebd) = seg(i).info;
                        gbl->i2wk(nebd++) = 1;
                    }
                bdnext1:        continue;
                }


                ebdry.resize(nebd);
                for(i=0;i<nebd;++i) {
                    ebdry(i) = getnewedgeobject(gbl->intwk(i),bdrymap);
                    ebdry(i)->alloc(static_cast<int>(gbl->i2wk(i)*grwfac));
                    ebdry(i)->nseg = 0;
                    gbl->intwk(i) = -1;
                    gbl->i2wk(i) = -1;
                }

                for(i=0;i<nseg;++i) {
                    if (seg(i).info) {
                        for (j = 0; j <nebd;++j) {
                            if (seg(i).info == ebdry(j)->idnum) {
                                ebdry(j)->seg(ebdry(j)->nseg++) = i;
                                goto bdnext1a;
                            }
                        }
                        *gbl->log << "Big error\n";
                        exit(1);
                    }
                bdnext1a:      continue;
                }
                
                for(i=0;i<nseg;++i)
                    gbl->i2wk_lst1(i) = i+1;
                  
                ntri = 0;
                triangulate(nseg);
                
                /* SOME TESTING FOR SPLINES */
//                for(i=0;i<nebd;++i) {
//                    ebdry(i)->input(in,boundary);
//                }
                
                in.close();
                
                break;
                
            case(vlength): {
                grd_app = grd_nm +".lngth";
                in.open(grd_app.c_str());
                if (!in) {
                    *gbl->log << "couldn't open vlength input file" << grd_app  << endl;
                    exit(1);
                }

                for(i=0;i<npnt;++i)
                    in >> lngth(i);

                return;
            }
             
            default:
                *gbl->log << "That filetype is not supported" << std::endl;
                exit(1);
    }
     
    for(i=0;i<nebd;++i) {
        /* CREATES NEW BOUNDARY FOR DISCONNECTED SEGMENTS OF SAME TYPE */
        ebdry(i)->reorder(); 
    }
    
    /* FIND ENDPOINT MATCHES */
    for(i=0;i<nvbd;++i) {
        /* Find two connecting boundary sides */
        for(j=0;j<nebd;++j) {
            if (seg(ebdry(j)->seg(0)).pnt(0) == vbdry(i)->pnt) {
                vbdry(i)->ebdry(1) = j;
                ebdry(j)->vbdry(0) = i;
            }
            if (seg(ebdry(j)->seg(ebdry(j)->nseg-1)).pnt(1) == vbdry(i)->pnt) {
                vbdry(i)->ebdry(0) = j;
                ebdry(j)->vbdry(1) = i;
            }
        }
    }


    
    bdrylabel();  
    

    createtritri();
    createpnttri();
    cnt_nbor();
    treeinit(); 
    
    grd_app = grd_nm + ".lngth";
    in.open(grd_app.c_str());
    if (in) {
        for(i=0;i<npnt;++i) in >> lngth(i);
        in.close();
    }
    else if (filetype != boundary) initlngth();
    
    
    tri_mesh::setinfo();
    checkintegrity();

    initialized = 1;
    
    return;
}

tri_mesh::~tri_mesh() {
    for(int i=0;i<nvbd;++i)
        delete vbdry(i);
    for(int i=0;i<nebd;++i)
        delete ebdry(i);
}
    

