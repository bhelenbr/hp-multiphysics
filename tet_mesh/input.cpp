#include "tet_mesh.h"
#include <utilities.h>
#include <cstdlib>
#include <cstring>
#include <input_map.h>
#include <iostream>

void tet_mesh::init(input_map &input, void *gin) {
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
            input.getwdefault("filetype",filetype,static_cast<int>(tet_mesh::grid));
        }
        
        keyword = gbl->idprefix + "_mesh";
        if (!input.get(keyword,filename)) {
            if (!input.get("mesh",filename))  {
                *gbl->log << "no mesh name\n";
                exit(1);
            }
        }
        coarse_level = 0;
        tet_mesh::input(filename.c_str(),static_cast<tet_mesh::filetype>(filetype),grwfac,input);
    }
	 
	 /* SET-UP BOUNDARY COMMUNICATIONS */
	 findmatch(gbl,coarse_level);
}

void tet_mesh::init(const multigrid_interface& mgin, init_purpose why, FLT sizereduce1d) {
    int i;
   
     if (!initialized) {
        const tet_mesh& inmesh = dynamic_cast<const tet_mesh &>(mgin);
        gbl = inmesh.gbl;
        maxvst =  MAX((int) (inmesh.maxvst/(sizereduce1d*sizereduce1d)),10);
        allocate(maxvst);
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
        otree.allocate((FLT (*)[ND]) pnts(0).data(), maxvst);
        coarse_level = inmesh.coarse_level+1;
        initialized = 1;
    }
	 
	 /* SET-UP BOUNDARY COMMUNICATIONS */
	 if (why == multigrid) findmatch(gbl,coarse_level);
}

void tet_mesh::allocate(int mxsize) {
    
    /* SIDE INFO */
    maxvst = mxsize;
    seg.resize(Range(-1,maxvst));

    /* VERTEX INFO */                        
    pnts.resize(maxvst);
    lngth.resize(maxvst);
    pnt.resize(Range(-1,maxvst));
    
    /* TRI INFO */ 
    tri.resize(Range(-1,maxvst));
    
    /* TET INFO */
    tet.resize(Range(-1,maxvst));

    otree.allocate((FLT (*)[ND]) pnts(0).data(), maxvst);

    if (!gbl) {
        /* gbl has not been set */
        /* so create internal gbl_struct */
        gbl = new global;
        gbl->idprefix = "";
        gbl->log = &std::cout;
        gbl->adapt_flag = true;
        gbl->tolerance = 1.25;
    }
    
    if (gbl->i1wk.ubound(firstDim) < maxvst) {
        // gbl->i1wk should always be kept initialized to -1
        gbl->i1wk.resize(Range(-1,maxvst));
        gbl->i1wk = -1;
        gbl->maxsrch = 100;
        
        gbl->fltwk.resize(maxvst);
        gbl->i2wk.resize(maxvst+1);
        gbl->i2wk.reindexSelf(TinyVector<int,1>(-1));
        // some smaller lists using i2 storage
        int mvst3 = maxvst/3;
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




void tet_mesh::input(const std::string &filename, tet_mesh::filetype filetype, FLT grwfac,input_map& bdrymap) {
    int temp,nbfaces;
    std::string grd_nm, bdry_nm, grd_app;
    TinyVector<int,3> v,s,e;
    ifstream in;
    int intskip;
    

    if (filename.substr(0,7) == "${HOME}") {
        grd_nm = getenv("HOME") +filename.substr(7,filename.length());
    }
    else 
        grd_nm = filename;

     switch (filetype) { 
                    
        case(baker):
            /* LOAD TETRAHEDRAL MESH DATA (BAKER FORMAT) */            
            grd_app = grd_nm +".d";
            in.open(grd_app.c_str());
            if (!in.good()) {
                *gbl->log << "error loading" << grd_nm << ".d" << std::endl;
                exit(1);
            }                    
                            
            /* LOAD NUMBER OF BDRY PTS, BDRY FACES, PTS, TETS, AND EDGES  */
            in.ignore(160,'\n');                         
            in >> intskip >> nbfaces >> npnt >> ntet >> nseg;    
        
            nfbd = 1; //TEMPORARY
            fbdry.resize(1);
                    
            maxvst = (4*ntet+nbfaces)/2;  // CHECK THIS TEMPORARY!!!!
            allocate(maxvst);
            
            /* LOAD PTS (npnt,3) */
            in.ignore(160,'\n'); 
            in.ignore(160,'\n');  
            
            fbdry(0) = getnewfaceobject(0,bdrymap);
            fbdry(0)->alloc(3*nbfaces/2);
            nvbd = 0;
            nebd = 0;
 
            for(int i = 0; i < npnt; ++i) {
                in >> pnts(i)(0) >> pnts(i)(1) >> pnts(i)(2) >> intskip;
            }
                                        
            /* LOAD SURFACE FACE DATA (nfaces,3) */
            in.ignore(160,'\n'); 
            in.ignore(160,'\n'); 
                                    
            for(int i = 0; i < nbfaces; ++i) {                
                in >> fbdry(0)->tri(i).pnt(0) >> fbdry(0)->tri(i).pnt(1) >> fbdry(0)->tri(i).pnt(2) >> intskip ;
                --fbdry(0)->tri(i).pnt(0);
                --fbdry(0)->tri(i).pnt(1);
                --fbdry(0)->tri(i).pnt(2);
            }
                
            fbdry(0)->ntri = nbfaces;
            
                 
            /* LOAD TETRAHEDRAL CONNECTION DATA (ntet,4) */            
            in.ignore(160,'\n'); 
            in.ignore(160,'\n');
            for(int i = 0; i < ntet; ++i){
                in >> tet(i).pnt(0) >> tet(i).pnt(1) >> tet(i).pnt(2) >> tet(i).pnt(3);
                --tet(i).pnt(0);
                --tet(i).pnt(1);
                --tet(i).pnt(2); 
                --tet(i).pnt(3); 
            } 
            

            /* LOAD seg CONNECTION DATA (nedg,2) */
            in.ignore(160,'\n'); 
            for(int i = 0; i < nseg; ++i){
                in >> seg(i).pnt(0) >> seg(i).pnt(1);
                --seg(i).pnt(0);
                --seg(i).pnt(1);
            }    
            
            /* LOAD TETRAHEDRAL seg DATA (ntet,6) */
            in.ignore(160,'\n'); 
            for(int i = 0; i < ntet; ++i){
                in >> tet(i).seg(0) >> tet(i).seg(1) >> tet(i).seg(2) >> tet(i).seg(3)>> tet(i).seg(4)>> tet(i).seg(5);
                --tet(i).seg(0);
                --tet(i).seg(1);
                --tet(i).seg(2);
                --tet(i).seg(3);
                --tet(i).seg(4);
                --tet(i).seg(5);
            }    
                        
            in.close(); 
            
            tet_mesh::fixvertexinfo();
            tet_mesh::createedgeinfo();
            tet_mesh::createfaceinfo();
            tet_mesh::morefaceinfo();
            tet_mesh::createtetinfo();
            tet_mesh::vertexnnbor();
            fbdry(0)->gbltolclvrtx();
            fbdry(0)->createsideinfo(); 
            fbdry(0)->createttri();
            fbdry(0)->gbltolcltri();
            fbdry(0)->gbltolclside();        
            fbdry(0)->cnt_nbor();
            fbdry(0)->createvtri();        
        
                    
            break;
            
        case (grid):
            /* LOAD TETRAHEDRAL ALL MESH DATA */            
            grd_app = grd_nm +".grd";
            in.open(grd_app.c_str());    
            if (!in.good()) {
                *gbl->log << "error loading" << grd_nm << ".grd" << std::endl;
                exit(1);
            }            

            in.ignore(80,':');
            in >> npnt;
            in.ignore(80,':');
            in >> nseg;
            in.ignore(80,':');
            in >> ntri;
            in.ignore(80,':');
            in >> ntet;
                        
            maxvst = ntri;  // CHECK THIS TEMPORARY!!!!
            allocate(maxvst);
            
            for(int i = 0; i < npnt; ++i) {
                in.ignore(80,':');
                in >> pnts(i)(0) >>  pnts(i)(1) >> pnts(i)(2);
            }
            
            for(int i = 0; i < nseg; ++i) {
                in.ignore(80,':');
                in >> seg(i).pnt(0)>> seg(i).pnt(1);
            }
    
            for(int i = 0; i < ntri; ++i){   
                in.ignore(80,':');
                in >> tri(i).pnt(0) >> tri(i).pnt(1)>> tri(i).pnt(2);                                
            }
                    
            for(int i = 0; i < ntet; ++i){
                in.ignore(80,':');
                for(int j = 0; j < 4; ++j)
                    in >> tet(i).pnt(j) ;
            }
            
                    
            in.ignore(80,':');
            in >> nvbd;
            vbdry.resize(nvbd);
            
            for(int i = 0; i < nvbd; ++i){
                in.ignore(80,':');
                in >> temp;
					 vbdry(i) = getnewvrtxobject(temp,bdrymap);
					 vbdry(i)->alloc(4);
                in.ignore(80,':');
                in >> vbdry(i)->pnt;

            }
            
            in.ignore(80,':');
            in >> nebd;
            ebdry.resize(nebd);
            for(int i = 0; i < nebd; ++i){
                in.ignore(80,':');
                in >> temp;
                ebdry(i) = getnewedgeobject(temp,bdrymap);
                in.ignore(80,':');
                in >> ebdry(i)->nseg;
                ebdry(i)->alloc(static_cast<int>(grwfac*ebdry(i)->nseg));
                for(int j=0;j<ebdry(i)->nseg;++j) {
                    in.ignore(80,':');
                    in >> ebdry(i)->seg(j).gindx;
                }
					 ebdry(i)->setup_next_prev();
					 ebdry(i)->reorder();
				}

            in.ignore(80,':');
            in >> nfbd;
            fbdry.resize(nfbd);
            
            for(int i = 0; i < nfbd; ++i){
                in.ignore(80,':');
                in >> temp;
                fbdry(i) = getnewfaceobject(temp,bdrymap);
                in.ignore(80,':');
                in >> fbdry(i)->npnt;
                in.ignore(80,':');
                in >> fbdry(i)->nseg;
                in.ignore(80,':');
                in >> fbdry(i)->ntri;                
                fbdry(i)->alloc(grwfac*3*(fbdry(i)->ntri));        //temporary fix        
                
                for(int j = 0; j < fbdry(i)->npnt; ++j){
                    in.ignore(80,':');
                    in >> fbdry(i)->pnt(j).gindx;
                }
                for(int j = 0; j < fbdry(i)->nseg; ++j){
                    in.ignore(80,':');
                    in >> fbdry(i)->seg(j).gindx;
                }
                for(int j = 0; j < fbdry(i)->ntri; ++j){
                    in.ignore(80,':');
                    in >> fbdry(i)->tri(j).gindx;
                }

            }    
            
            in.close();
            
            tet_mesh::vertexnnbor();
            tet_mesh::edgeinfo();
            tet_mesh::faceinfo();
            tet_mesh::morefaceinfo();
            tet_mesh::createtetinfo();
            
            
            for(int i = 0; i < nfbd; ++i){
                fbdry(i)->allinfo(); 
                fbdry(i)->cnt_nbor();
                fbdry(i)->createvtri();    
            }
            
            break;

            
          case(gambit):
                grd_app = grd_nm +".FDNEUT";
                in.open(grd_app.c_str());
                if (!in.good()) {
                          *gbl->log << "trouble opening " << grd_nm << std::endl;
                          exit(1);
                }
                
                for(int i=0;i<5;++i)
                    in.ignore(160,'\n');
                     
                in >> npnt >> ntet >> nfbd;
                nfbd -= 1;
                
                fbdry.resize(nfbd);
                /* MAXVST IS APPROXIMATELY THE NUMBER OF ELEMENTS  */
                /* FOR EACH ELEMENT THERE ARE APPROXIMATELY 3/2 SIDES */
                maxvst = static_cast<int>((grwfac*3*ntet)); //temporary fix 
                allocate(maxvst);
       
                for(int i=0;i< 8;++i)
                    in.ignore(160,'\n');

                /* READ VERTEX DATA */     
                for(int i=0;i<npnt;++i) {
                    in >> intskip;
                    for(int n=0;n<ND;++n)
                        in >> pnts(i)(n);
                    pnt(i).info = -1;
                    in.ignore(160,'\n');
                }
                     
                for(int i=0;i<3;++i)
                     in.ignore(160,'\n');
                     
                in.ignore(20,':');
                in.ignore(20,':');

                in >> ntet;    
                            
                for(int i=0;i<2;++i)
                    in.ignore(160,'\n');            
                                         
                for(int i=0;i<ntet;++i) {
                    in >> intskip >> tet(i).pnt(0) >> tet(i).pnt(1) >> tet(i).pnt(2)>> tet(i).pnt(3);
                     --tet(i).pnt(0);
                     --tet(i).pnt(1);
                     --tet(i).pnt(2);
                     --tet(i).pnt(3);
                     tet(i).info = -1;
                }
                                
                in.ignore(160,'\n');    
                                            
                for(int i = 0; i < nfbd; ++i){
                    fbdry(i) = getnewfaceobject(i,bdrymap);
                    in.ignore(20,':');
                    in >> intskip;
                    in.ignore(20,':');
                    in >> fbdry(i)->ntri;
                    fbdry(i)->alloc(grwfac*3*(fbdry(i)->ntri));    //temporary fix
                    in.ignore(160,'\n');
                    in.ignore(160,'\n');
                    for(int j = 0; j < fbdry(i)->ntri; ++j){
                        in >> intskip >> fbdry(i)->tri(j).pnt(0) >> fbdry(i)->tri(j).pnt(1) >> fbdry(i)->tri(j).pnt(2);
                        --fbdry(i)->tri(j).pnt(0);
                        --fbdry(i)->tri(j).pnt(1);
                        --fbdry(i)->tri(j).pnt(2);
                    }
                }
     
                in.close();
                 
                fixvertexinfo();
                createedgeinfo();
                createfaceinfo();
                morefaceinfo();
                createtetinfo();
                vertexnnbor();
                for(int i = 0; i < nfbd; ++i){
                    fbdry(i)->gbltolclvrtx();
                    fbdry(i)->createsideinfo(); 
                    fbdry(i)->createttri();
                    fbdry(i)->gbltolcltri();
                    fbdry(i)->gbltolclside();        
                    fbdry(i)->cnt_nbor();
                    fbdry(i)->createvtri();    
                }
                
                //tet_mesh::createfaceorientation();                  
                //tet_mesh::test();
                break;
                
  

            case(tecplot): {
                grd_app = grd_nm +".dat";
                in.open(grd_app.c_str());
                if (!in.good()) {
                    *gbl->log << "couldn't open tecplot file: " << grd_app<< std::endl;
                    exit(1);
                }
                
                /* HEADER LINES */
                in.ignore(80,',');
                in.ignore(80,',');
                in.ignore(80,'=');
                in >> npnt;
                in.ignore(80,'=');
                in >> ntet;
                
                /* MAXVST IS APPROXIMATELY THE NUMBER OF ELEMENTS  */
                /* FOR EACH ELEMENT THERE ARE APPROXIMATELY 3/2 SIDES */
                if (!initialized) {
                    maxvst = static_cast<int>(grwfac*3*ntet);
                    allocate(maxvst);
                }
                else if ((3*ntet) > maxvst) {
                    *gbl->log << "mesh is too large" << std::endl;
                    exit(1);
                }
                    
                /* READ VERTEX DATA */
                for(int i=0;i<npnt;++i) {
                    for(int n=0;n<ND;++n)
                        in >> pnts(i)(n);
                    in.ignore(80,'\n');
                    pnt(i).info = -1;
                }

                while (in.get() == '#');
                    in.ignore(80,'\n');
        
                for(int i=0;i<ntet;++i) {
                    in >> tet(i).pnt(0) >> tet(i).pnt(1) >> tet(i).pnt(2)>> tet(i).pnt(3);
                    --tet(i).pnt(0);
                    --tet(i).pnt(1);
                    --tet(i).pnt(2);
                    --tet(i).pnt(3);
                    tet(i).info = -1;
                    std::cout << tet(i).pnt(0) << tet(i).pnt(1) << tet(i).pnt(2) << tet(i).pnt(3) << std::endl;
                    
                }
                
                /*    CREATE ALL MESH INFO */
                tet_mesh::fixvertexinfo();
                tet_mesh::createedgeinfo();
                tet_mesh::createfaceinfo();
                tet_mesh::morefaceinfo();
                tet_mesh::createtetinfo();

                /* FIND ALL BOUNDARY SIDES */
                int count = 0;
                for(int i=0;i<ntri;++i)
                    if (tri(i).tet(1) < 0) ++count;

                nfbd = 1;
                fbdry.resize(1);
                fbdry(0) = getnewfaceobject(1,bdrymap);
                fbdry(0)->alloc(static_cast<int>(3*grwfac*count));
                fbdry(0)->ntri = count;
                count = 0;
                for(int i=0;i<ntri;++i) {
                    if (tri(i).tet(1) < 0) {
                        for (int j=0;j<3;++j)
                            fbdry(0)->tri(count++).pnt(j) = tri(i).pnt(j);
                    }
                }
                
                for(int i = 0; i < nfbd; ++i){
                    fbdry(i)->gbltolclvrtx();
                    fbdry(i)->createsideinfo(); 
                    fbdry(i)->createttri();
                    fbdry(i)->gbltolcltri();
                    fbdry(i)->gbltolclside();        
                    fbdry(i)->cnt_nbor();
                    fbdry(i)->createvtri();    
                }
                
                nvbd = 0;
                nebd = 0;
                in.close();
                
                break;
            }
                

             
            default:
                *gbl->log << "That filetype is not supported" << std::endl;
                exit(1);
    }
    
    /* FIND ENDPOINT MATCHES */
    for(int i=0;i<nvbd;++i) {
        /* Find two connecting boundary sides */
        for(int j=0;j<nebd;++j) {
            if (seg(ebdry(j)->seg(0).gindx).pnt(0) == vbdry(i)->pnt) {
                ebdry(j)->vbdry(0) = i;
            }
            if (seg(ebdry(j)->seg(ebdry(j)->nseg-1).gindx).pnt(1) == vbdry(i)->pnt) {
                ebdry(j)->vbdry(1) = i;
            }
        }
    }

    
    bdrylabel();  // MAKES BOUNDARY ELEMENTS POINT TO BOUNDARY GROUP/ELEMENT
    treeinit(); 
    
    grd_app = grd_nm + ".lngth";
    in.open(grd_app.c_str());
    if (in) {
        for(int i=0;i<npnt;++i) in >> lngth(i);
        in.close();
    }
    else if (filetype != boundary) initlngth();
    
    tet_mesh::setinfo();
    checkintegrity();

    initialized = 1;
    
    return;
}

tet_mesh::~tet_mesh() {
    for(int i=0;i<nvbd;++i)
        delete vbdry(i);
    for(int i=0;i<nebd;++i)
        delete ebdry(i);
    for(int i=0;i<nfbd;++i)
        delete fbdry(i);
}
    
