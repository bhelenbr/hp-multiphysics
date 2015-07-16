#include "tet_mesh.h"
#include <utilities.h>
#include <cstdlib>
#include <cstring>
#include <input_map.h>
#include <iostream>
#ifdef USING_MADLIB
#include "MAdLibInterface.h"
#endif

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
		gbl->adapt_interval = 1;
		gbl->tolerance = 1.25;
	}
	
	if (!initialized) {
		FLT grwfac;
		keyword = gbl->idprefix + "_growth factor";
		if (!input.get(keyword,grwfac)) {
			input.getwdefault("growth factor",grwfac,2.0);
		}
		if (grwfac < 1.0) {
			*gbl->log << "growth factor must be greater than one\n";
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
		
		keyword = gbl->idprefix + "_mesh";
		if (!input.get(keyword,filename)) {
			if (input.get("mesh",filename)) {
				filename = filename +"_" +gbl->idprefix;
			}
			else {
				*gbl->log << "no mesh name" << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
		}
		
		int filetype;
		keyword = gbl->idprefix + "_filetype";
		if (!input.get(keyword,filetype)) {
			input.getwdefault("filetype",filetype,static_cast<int>(tet_mesh::grid));
		}
		
		coarse_level = 0;
		tet_mesh::input(filename.c_str(),static_cast<tet_mesh::filetype>(filetype),grwfac,input);
	}
	 
	/* SET-UP BOUNDARY COMMUNICATIONS */
	findmatch(gbl,coarse_level);
}

void tet_mesh::init(const multigrid_interface& mgin, init_purpose why, FLT sizereduce1d) {
	int i;
	const tet_mesh& inmesh = dynamic_cast<const tet_mesh &>(mgin);

	if (!initialized) {
		gbl = inmesh.gbl;
		maxvst =  MAX((int) (inmesh.maxvst/(sizereduce1d*sizereduce1d*sizereduce1d)),10);
		maxvst =  inmesh.maxvst;

		allocate(maxvst);
		nfbd = inmesh.nfbd;
		fbdry.resize(nfbd);
		for(i=0;i<nfbd;++i) {
			fbdry(i) = inmesh.fbdry(i)->create(*this);
			//fbdry(i)->alloc(MAX(static_cast<int>(inmesh.fbdry(i)->maxpst/(sizereduce1d*sizereduce1d)),10));
			fbdry(i)->alloc(inmesh.fbdry(i)->maxpst);
		}		
		nebd = inmesh.nebd;
		ebdry.resize(nebd);
		for(i=0;i<nebd;++i) {
			ebdry(i) = inmesh.ebdry(i)->create(*this);
			//ebdry(i)->alloc(MAX(static_cast<int>(inmesh.ebdry(i)->maxseg/sizereduce1d),10));
			ebdry(i)->alloc(inmesh.ebdry(i)->maxseg);
		}
		nvbd = inmesh.nvbd;
		vbdry.resize(nvbd);
		for(i=0;i<nvbd;++i) {
			vbdry(i) = inmesh.vbdry(i)->create(*this);
			vbdry(i)->alloc(4);
		}
		otree.allocate(pnts, maxvst);
		if (why == multigrid) coarse_level = inmesh.coarse_level+1;
		else coarse_level = 0;
		initialized = 1;
	}
	 
	/* SET-UP BOUNDARY COMMUNICATIONS */
	if (why == multigrid) {
		coarsen(2.0,inmesh);
//		input_map empty;
//		input("/Users/helenbrk/Codes/ucl/Meshes/Squares/tet8.msh",gmsh,1.0,empty);
		findmatch(gbl,coarse_level);
	}
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

	otree.allocate(pnts, maxvst);

	if (!gbl) {
		/* gbl has not been set */
		/* so create internal gbl_struct */
		gbl = new global;
		gbl->idprefix = "";
		gbl->log = &std::cout;
		gbl->adapt_interval = 1;
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




void tet_mesh::input(const std::string &filename, tet_mesh::filetype filetype, FLT grwfac,input_map& inmap) {
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

	/* Override filetype based on ending? */
	size_t dotloc;
	dotloc = grd_nm.find_last_of('.');
	string ending;
	ending = grd_nm.substr(dotloc+1);
	if (ending == "d") {
		filetype = baker;
		grd_nm = grd_nm.substr(0,dotloc);
	}
	else if (ending == "grd") {
		filetype = grid;
		grd_nm = grd_nm.substr(0,dotloc);
	}
	else if (ending == "FDNEUT") {
		filetype = gambit;
		grd_nm = grd_nm.substr(0,dotloc);
	}
	else if (ending == "dat") {
		filetype = tecplot;
		grd_nm = grd_nm.substr(0,dotloc);
	}
	else if (ending == "msh") {
		filetype = gmsh;
		grd_nm = grd_nm.substr(0,dotloc);
	}

	switch (filetype) { 
					
		case(baker):
			/* LOAD TETRAHEDRAL MESH DATA (BAKER FORMAT) */            
			grd_app = grd_nm +".d";
			in.open(grd_app.c_str());
			if (!in) {
				*gbl->log << "error loading" << grd_nm << ".d" << std::endl;
				exit(1);
			}                    
							
			/* LOAD NUMBER OF BDRY PTS, BDRY FACES, PTS, TETS, AND EDGES  */
			in.ignore(160,'\n');                         
			in >> intskip >> nbfaces >> npnt >> ntet >> nseg;    
		
			nfbd = 1; //TEMPORARY
			fbdry.resize(1);
					
			maxvst = (4*ntet+nbfaces)/2;  // CHECK THIS TEMPORARY!!!!
			allocate(static_cast<int>(grwfac*maxvst));
			
			/* LOAD PTS (npnt,3) */
			in.ignore(160,'\n'); 
			in.ignore(160,'\n');  
			
			fbdry(0) = getnewfaceobject(0,inmap);
			fbdry(0)->alloc(static_cast<int>(grwfac*3*nbfaces/2));
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
			for(int i = 0; i < ntet; ++i) {
				in >> tet(i).pnt(0) >> tet(i).pnt(1) >> tet(i).pnt(2) >> tet(i).pnt(3);
				--tet(i).pnt(0);
				--tet(i).pnt(1);
				--tet(i).pnt(2); 
				--tet(i).pnt(3); 
			} 
			

			/* LOAD seg CONNECTION DATA (nedg,2) */
			in.ignore(160,'\n'); 
			for(int i = 0; i < nseg; ++i) {
				in >> seg(i).pnt(0) >> seg(i).pnt(1);
				--seg(i).pnt(0);
				--seg(i).pnt(1);
			}    
			
			/* LOAD TETRAHEDRAL seg DATA (ntet,6) */
			in.ignore(160,'\n'); 
			for(int i = 0; i < ntet; ++i) {
				in >> tet(i).seg(0) >> tet(i).seg(1) >> tet(i).seg(2) >> tet(i).seg(3)>> tet(i).seg(4)>> tet(i).seg(5);
				--tet(i).seg(0);
				--tet(i).seg(1);
				--tet(i).seg(2);
				--tet(i).seg(3);
				--tet(i).seg(4);
				--tet(i).seg(5);
			}    
						
			in.close(); 
			
			reorient_tets();
			create_from_tet_definitions();
			fbdry(0)->create_from_pnt();      
					
			break;
			
		case (grid):
			/* LOAD TETRAHEDRAL ALL MESH DATA */            
			grd_app = grd_nm +".grd";
			in.open(grd_app.c_str());    
			if (!in) {
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
			allocate(static_cast<int>(grwfac*maxvst));
			
			for(int i = 0; i < npnt; ++i) {
				in.ignore(80,':');
				in >> pnts(i)(0) >>  pnts(i)(1) >> pnts(i)(2);
			}
			
			for(int i = 0; i < nseg; ++i) {
				in.ignore(80,':');
				in >> seg(i).pnt(0)>> seg(i).pnt(1);
			}
	
			for(int i = 0; i < ntri; ++i) {   
				in.ignore(80,':');
				in >> tri(i).pnt(0) >> tri(i).pnt(1)>> tri(i).pnt(2);                                
			}
					
			for(int i = 0; i < ntet; ++i) {
				in.ignore(80,':');
				for(int j = 0; j < 4; ++j)
					in >> tet(i).pnt(j) ;
			}
			
					
			in.ignore(80,':');
			in >> nvbd;
			vbdry.resize(nvbd);
			
			for(int i = 0; i < nvbd; ++i) {
				in.ignore(80,':');
				in >> temp;
				vbdry(i) = getnewvrtxobject(temp,inmap);
				vbdry(i)->alloc(4);
				in.ignore(80,':');
				in >> vbdry(i)->pnt;
			}
			
			in.ignore(80,':');
			in >> nebd;
			ebdry.resize(nebd);
			for(int i = 0; i < nebd; ++i) {
				in.ignore(80,':');
				in >> temp;
				ebdry(i) = getnewedgeobject(temp,inmap);
				in.ignore(80,':');
				in >> ebdry(i)->nseg;
				ebdry(i)->alloc(static_cast<int>(20*grwfac*ebdry(i)->nseg));
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
				
			for(int i = 0; i < nfbd; ++i) {
				in.ignore(80,':');
				in >> temp;
				fbdry(i) = getnewfaceobject(temp,inmap);
				in.ignore(80,':');
				in >> fbdry(i)->npnt;
				in.ignore(80,':');
				in >> fbdry(i)->nseg;
				in.ignore(80,':');
				in >> fbdry(i)->ntri;                
				
				fbdry(i)->alloc(static_cast<int>(grwfac*10*(fbdry(i)->ntri)));//temporary fix        

				for(int j = 0; j < fbdry(i)->npnt; ++j) {
					in.ignore(80,':');
					in >> fbdry(i)->pnt(j).gindx;
				}
				for(int j = 0; j < fbdry(i)->nseg; ++j) {
					in.ignore(80,':');
					in >> fbdry(i)->seg(j).gindx;
				}
				for(int j = 0; j < fbdry(i)->ntri; ++j) {
					in.ignore(80,':');
					in >> fbdry(i)->tri(j).gindx;
				}
			}    
			
			in.close();
			
			create_from_pnt_definitions();
			
			for(int i = 0; i < nfbd; ++i) {
				fbdry(i)->create_from_gindx();
			}
			
			break;

			
		case(gambit):
			grd_app = grd_nm +".FDNEUT";
			in.open(grd_app.c_str());
			if (!in) {
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

			for(int i = 0; i < nfbd; ++i) {
				fbdry(i) = getnewfaceobject(i,inmap);
				in.ignore(20,':');
				in >> intskip;
				in.ignore(20,':');
				in >> fbdry(i)->ntri;
				
				fbdry(i)->alloc(static_cast<int>(grwfac*3*(fbdry(i)->ntri)));    //temporary fix
				in.ignore(160,'\n');
				in.ignore(160,'\n');
				for(int j = 0; j < fbdry(i)->ntri; ++j) {
					in >> intskip >> fbdry(i)->tri(j).pnt(0) >> fbdry(i)->tri(j).pnt(1) >> fbdry(i)->tri(j).pnt(2);
					--fbdry(i)->tri(j).pnt(0);
					--fbdry(i)->tri(j).pnt(1);
					--fbdry(i)->tri(j).pnt(2);
				}
			}

			in.close();

			reorient_tets();
			create_from_tet_definitions();
			for(int i = 0; i < nfbd; ++i) 
				fbdry(i)->create_from_pnt(); 
			
			break;
			


		case(tecplot): {
			grd_app = grd_nm +".dat";
			in.open(grd_app.c_str());
			if (!in) {
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
			}
			
			/* CREATE ALL MESH INFO */
			reorient_tets();
			create_from_tet_definitions();

			/* FIND ALL BOUNDARY SIDES */
			int count = 0;
			for(int i=0;i<ntri;++i)
				if (tri(i).tet(1) < 0) ++count;

			nfbd = 1;
			fbdry.resize(1);
			fbdry(0) = getnewfaceobject(1,inmap);
			fbdry(0)->alloc(static_cast<int>(3*grwfac*count));
			fbdry(0)->ntri = count;
			count = 0;
			for(int i=0;i<ntri;++i) {
				if (tri(i).tet(1) < 0) {
					for (int j=0;j<3;++j)
						fbdry(0)->tri(count++).pnt(j) = tri(i).pnt(j);
				}
			}
			
			for(int i = 0; i < nfbd; ++i)
				fbdry(i)->create_from_pnt();  
			
			nvbd = 0;
			nebd = 0;
			in.close();
			
			break;
		}

#ifdef USING_MADLIB
		case(gmsh): {
			MAdLib_input(grd_nm, grwfac, inmap);
			break;
		}
#endif
			


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
	for (int i=0;i<nfbd;++i)
		fbdry(i)->treeinit();
	
	grd_app = grd_nm + ".lngth";
	in.open(grd_app.c_str());
	if (in) {
		for(int i=0;i<npnt;++i)
			in >> lngth(i);
		in.close();
	}
	else if (filetype != boundary) initlngth();
	
	tet_mesh::setinfo();
	output("error",tet_mesh::grid);
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
	
	// FIXME: WHO IS RESPONSIBLE FOR THIS?
	if (gbl) delete gbl;
	gbl = 0;
}
	

