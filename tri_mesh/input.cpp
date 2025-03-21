#include "tri_mesh.h"
#include <cstdlib>
#include <cstring>
#include <input_map.h>
#include <iostream>
#ifdef libbinio
#include <libbinio/binfile.h>
#endif
#include <netcdf.h>

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
#define ERR(e) {*gbl->log << "netCDF error " <<  nc_strerror(e); sim::abort(__LINE__,__FILE__,gbl->log);}


/* NAMING CONVENTION */
// pt, segment, triangle, tetrahedral

/* Boundaries */
// vertex, edge, face

void tri_mesh::init(input_map &input, shared_ptr<block_global> gin) {
	std::string keyword;
	std::istringstream data;
	std::string filename;
	std::string bdryfile;

	if (gin != 0) {
		gbl = gin;
	}
	else {
		/* gbl has not been set */
		/* so create internal gbl_struct */
        gbl = make_shared<block_global>();
		gbl->idprefix = std::string("b0");
		gbl->log = &std::cout;
		gbl->adapt_interval	= 1;
		gbl->tolerance = sqrt(2.0);
	}
    
    tri_gbl = make_shared<tri_global>();

	if (!initialized) {
		FLT grwfac;
		keyword = gbl->idprefix + "_growth factor";
		if (!input.get(keyword,grwfac)) {
			input.getwdefault("growth factor",grwfac,1.0);
		}
		if (grwfac > 0.0 && grwfac < 1.0) {
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
			input.getwdefault("filetype",filetype,static_cast<int>(tri_mesh::grid));
		}
		
		coarse_level = 0;
		tri_mesh::input(filename.c_str(),static_cast<tri_mesh::filetype>(filetype),grwfac,input);

		if (!input.get(gbl->idprefix+"_maximum_length",gbl->max_length)) {
			input.getwdefault("maximum_length",gbl->max_length,0.5*(qtree.xmax(0)-qtree.xmin(0) +qtree.xmax(1)-qtree.xmin(1)));
		}
		
		if (!input.get(gbl->idprefix+"_minimum_length",gbl->min_length)) {
			input.getwdefault("minimum_length",gbl->min_length,-1.0);
		}
	}
    findmatch(gbl,coarse_level);
}

void tri_mesh::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
	const tri_mesh& inmesh = dynamic_cast<const tri_mesh &>(in);

    gbl = inmesh.gbl;
    tri_gbl = inmesh.tri_gbl;

	if (!initialized) {
		maxpst =  MAX((int) (inmesh.maxpst/(sizereduce1d*sizereduce1d)),20);
		allocate(maxpst);
		nebd = inmesh.nebd;
		ebdry.resize(nebd);
		for(int i=0;i<nebd;++i) {
			ebdry(i) = inmesh.ebdry(i)->create(*this);
			ebdry(i)->alloc(MAX(static_cast<int>(inmesh.ebdry(i)->maxseg/sizereduce1d),10));
		}
		nvbd = inmesh.nvbd;
		vbdry.resize(nvbd);
		for(int i=0;i<nvbd;++i) {
			vbdry(i) = inmesh.vbdry(i)->create(*this);
			vbdry(i)->alloc(4);
		}
		qtree.allocate(pnts, maxpst);
		if (why == multigrid) coarse_level = inmesh.coarse_level+1;
		else coarse_level = 0;
		initialized = 1;
	}

	if (why == multigrid) {
		coarsen(1.6,inmesh);
		findmatch(gbl,coarse_level);
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

	qtree.allocate(pnts, maxpst);

	if (!gbl) {
		/* gbl has not been set */
		/* so create internal gbl_struct */
		gbl = make_shared<block_global>();
		gbl->idprefix = "";
		gbl->log = &std::cout;
		gbl->adapt_interval = 1;
		gbl->tolerance = sqrt(2.0);
        tri_gbl = make_shared<tri_global>();
	}

	if (tri_gbl->intwk.ubound(firstDim) < maxpst) {
		// tri_gbl->intwk should always be kept initialized to -1
		tri_gbl->intwk.resize(Range(-1,maxpst));
		tri_gbl->intwk = -1;
		tri_gbl->maxsrch = 1000;

		tri_gbl->fltwk.resize(maxpst);
		tri_gbl->i2wk.resize(maxpst+1);
		tri_gbl->i2wk.reindexSelf(TinyVector<int,1>(-1));
		// some smaller lists using i2 storage
//		int mvst3 = maxpst/3;
//		Array<int,1> temp1(tri_gbl->i2wk.data(),mvst3,neverDeleteData);
//		tri_gbl->i2wk_lst1.reference(temp1);
//		tri_gbl->i2wk_lst1.reindexSelf(TinyVector<int,1>(-1));
//		Array<int,1> temp2(tri_gbl->i2wk.data()+1+mvst3,mvst3-1,neverDeleteData);
//		tri_gbl->i2wk_lst2.reference(temp2);
//		tri_gbl->i2wk_lst2.reindexSelf(TinyVector<int,1>(-1));
//		Array<int,1> temp3(tri_gbl->i2wk.data()+1+2*mvst3,mvst3-1,neverDeleteData);
//		tri_gbl->i2wk_lst3.reference(temp3);
//		tri_gbl->i2wk_lst3.reindexSelf(TinyVector<int,1>(-1));
        
        tri_gbl->i2wk_lst1.resize(maxpst+1);
        tri_gbl->i2wk_lst1.reindexSelf(TinyVector<int,1>(-1));
        tri_gbl->i2wk_lst2.resize(maxpst+1);
        tri_gbl->i2wk_lst2.reindexSelf(TinyVector<int,1>(-1));
        tri_gbl->i2wk_lst3.resize(maxpst+1);
        tri_gbl->i2wk_lst3.reindexSelf(TinyVector<int,1>(-1));
	}

	initialized = 1;
    *gbl->log << "#Allocated mesh with maxpst = " << mxsize << std::endl;

	return;
}


void tri_mesh::input(const std::string &filename, tri_mesh::filetype filetype, FLT grwfac,input_map& bdrymap) {
	int sind,count,temp,interior_pts=0;
	std::string grd_nm, bdry_nm, grd_app;
	TinyVector<int,3> v,s,e;
	ifstream in;
#ifdef libbinio
	binifstream bin;
#endif
	FLT fltskip;
	int intskip;
	Array<Array<TinyVector<int,2>,1>,1> svrtxbtemp;  // TEMPORARY FOR LOADING GAMBIT

	if (filename.substr(0,7) == "${HOME}") {
		grd_nm = getenv("HOME") +filename.substr(7,filename.length());
	}
	else
		grd_nm = filename;
    
    FLT grwfac1d = sqrt(grwfac);
		
	/* Override filetype based on ending? */
	size_t dotloc;
	dotloc = grd_nm.find_last_of('.');
	string ending;
	ending = grd_nm.substr(dotloc+1);
	if (ending == "grd") {
		filetype = grid;
		grd_nm = grd_nm.substr(0,dotloc);
	}
	else if (ending == "d") {
		filetype = boundary;
		grd_nm = grd_nm.substr(0,dotloc);
	}
	else if (ending == "bin") {
		filetype = binary;
		grd_nm = grd_nm.substr(0,dotloc);
	}
	else if (ending == "dat") {
		filetype = tecplot;
		grd_nm = grd_nm.substr(0,dotloc);
	}
	else if (ending == "nc") {
		filetype = netcdf;
		grd_nm = grd_nm.substr(0,dotloc);
	}
	

	switch (filetype) {
		case(easymesh): {
			/* LOAD SIDE INFORMATION */
			grd_app = grd_nm + ".s";
			in.open(grd_app.c_str());
			if (!in) {
				*gbl->log << "error: couldn't open file: " << grd_app << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			if (!(in >> nseg)) {
				*gbl->log << "error: in segment file: " << grd_app << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			
            if (!initialized) {
                if (grwfac > 0.0) {
                    /* Scaling factor */
                    allocate(nseg + (int) (grwfac*nseg));
                }
                else {
                    /* absolute size */
                    allocate(-(int)(grwfac));
                }
			}
			else if (nseg > maxpst) {
				*gbl->log << "error: mesh is too large" << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			
			for(int i=0;i<nseg;++i) {
				in.ignore(80,':');
				in >> seg(i).pnt(0) >> seg(i).pnt(1) >> seg(i).tri(0) >> seg(i).tri(1) >> seg(i).info;
				if(in.fail()) {
					*gbl->log << "error: in segment file " << grd_app << std::endl;
					sim::abort(__LINE__,__FILE__,gbl->log);
				}
			}
			
			/* Don't Trust File So Will Calculate tri myself */
			for(int i=0;i<maxpst;++i) {
				seg(i).tri(0) = -1;
				seg(i).tri(1) = -1;
			}
			
			/* COUNT BOUNDARY GROUPS */
			nebd = 0;
			for(int i=0;i<nseg;++i) {
				if (seg(i).info) {
					for (int j = 0; j <nebd;++j) {
						if (seg(i).info == tri_gbl->intwk(j)) {
							++tri_gbl->i2wk(j);
							goto next1;
						}
					}
					/* NEW SIDE */
					tri_gbl->intwk(nebd) = seg(i).info;
					tri_gbl->i2wk(nebd++) = 1;
				}
			next1:        continue;
			}
			
			ebdry.resize(nebd);
			for(int i=0;i<nebd;++i) {
				ebdry(i) = getnewedgeobject(tri_gbl->intwk(i),bdrymap);
				ebdry(i)->alloc(static_cast<int>(tri_gbl->i2wk(i)*4*grwfac1d));
				ebdry(i)->nseg = 0;
				tri_gbl->intwk(i) = -1;
				tri_gbl->i2wk(i) = -1;
			}
			
			
			for(int i=0;i<nseg;++i) {
				if (seg(i).info) {
					for (int j = 0; j <nebd;++j) {
						if (seg(i).info == ebdry(j)->idnum) {
							ebdry(j)->seg(ebdry(j)->nseg++) = i;
							goto next1a;
						}
					}
					*gbl->log << "Big error\n";
					sim::abort(__LINE__,__FILE__,gbl->log);
				}
			next1a:      continue;
			}
			in.close();
			
			/* LOAD VERTEX INFORMATION                    */
			grd_app = grd_nm + ".n";
			in.open(grd_app.c_str());
			if (!in) { *gbl->log << "trouble opening grid" << grd_app << std::endl; sim::abort(__LINE__,__FILE__,gbl->log);}
			
			if(!(in >> npnt)) {
				*gbl->log << "1: error in grid: " << grd_app << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			
			for(int i=0;i<npnt;++i) {
				in.ignore(80,':');
				for(int n=0;n<ND;++n) {
					in >> pnts(i)(n);
				}
				in >> pnt(i).info;
				if (in.fail())  { *gbl->log << "2b: error in grid" << std::endl; sim::abort(__LINE__,__FILE__,gbl->log); }
			}
			in.close();
			
			/* COUNT VERTEX BOUNDARY GROUPS  */
			nvbd = 0;
			for(int i=0;i<npnt;++i)
				if (pnt(i).info) ++nvbd;
			vbdry.resize(nvbd);
			
			nvbd = 0;
			for(int i=0;i<npnt;++i) {
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
			if (!in) {
				*gbl->log << "trouble opening " << grd_app << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			if(!(in >> ntri)) {
				*gbl->log << "error in file " << grd_app << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			
			for(int i=0;i<ntri;++i) {
				in.ignore(80,':');
				in >> v(0) >> v(1) >> v(2) >> e(0) >> e(1) >> e(2) >> s(0) >> s(1) >> s(2) >> fltskip >> fltskip >> tri(i).info;
				
				for (int j=0;j<3;++j) {
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
		}
		case(gambit): {
			grd_app = grd_nm +".FDNEUT";
			in.open(grd_app.c_str());
			if (!in) {
				*gbl->log << "trouble opening " << grd_nm << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			
			for(int i=0;i<5;++i)
				in.ignore(160,'\n');
			
            int i;
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
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			
			for(int i=0;i<8;++i)
				in.ignore(160,'\n');
			
			/* READ VERTEX DATA */
			for(int i=0;i<npnt;++i) {
				in >> intskip;
				for(int n=0;n<ND;++n)
					in >> pnts(i)(n);
				pnt(i).info = -1;
			}
			
			for(int i=0;i<3;++i)
				in.ignore(160,'\n');
			
			/* READ ELEMENT DATA */
			// fscanf(grd,"%*[^0-9]%*d%*[^0-9]%d%*[^\n]\n",&ntri); TEMPORARY
			in.ignore(160,'\n');
			
			for(int i=0;i<ntri;++i) {
				in >> intskip >> tri(i).pnt(0) >> tri(i).pnt(1) >> tri(i).pnt(2);
				--tri(i).pnt(0);
				--tri(i).pnt(1);
				--tri(i).pnt(2);
				tri(i).info = 0;
			}
			
			/* READ BOUNDARY DATA STORE TEMPORARILY */
			svrtxbtemp.resize(10);
			
			for(int i=0;i<nebd;++i) {
				// fscanf(grd,"%*[^0-9]%*d%*[^0-9]%d%*[^0-9]%*d%*[^0-9]%*d%*[^0-9]%d\n",&count,&temp); // FIXME
				count = 0;
				
				temp = 2;
				ebdry(i) = getnewedgeobject(temp,bdrymap);
				ebdry(i)->alloc(static_cast<int>(count*4*grwfac1d));
				ebdry(i)->nseg = count;
				
				in.ignore(160,'\n');
				
				svrtxbtemp(i).resize(ebdry(i)->nseg);
				
				for(int j=0;j<ebdry(i)->nseg;++j) {
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
			for(int i=0;i<nseg;++i)
				if (seg(i).tri(1) < 0)
					pnt(seg(i).pnt(0)).info = i;
			
			for(int i=0;i<nseg;++i)
				seg(i).info = 0;
			
			/* MATCH BOUNDARY SIDES TO GROUPS */
			for(int i=0;i<nebd;++i) {
				for(int j=0;j<ebdry(i)->nseg;++j) {
					sind = pnt(svrtxbtemp(i)(j)(0)).info;
					if (sind < 0) {
						*gbl->log << "error in boundary information " << i << j << std::endl;
						sim::abort(__LINE__,__FILE__,gbl->log);
					}
					if (seg(sind).pnt(1) == svrtxbtemp(i)(j)(1)) {
						seg(sind).info = ebdry(i)->idnum;
						ebdry(i)->seg(j) = sind;
						pnt(seg(sind).pnt(0)).info = 0;
					}
					else {
						*gbl->log << "Error: boundary sides are not counterclockwise " <<
						svrtxbtemp(i)(j)(0) << svrtxbtemp(i)(j)(1) << std::endl;
						sim::abort(__LINE__,__FILE__,gbl->log);
					}
				}
			}
			
			for(int i=0;i<npnt;++i)
				pnt(i).info = 0;
			
			in.close();
			~svrtxbtemp;
			
			break;
		}
		case(grid): {
			grd_app = grd_nm + ".grd";
			in.open(grd_app.c_str());
			if (!in) {
				*gbl->log << "couldn't open grid file: " << grd_app << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
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
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			
			/* VRTX INFO */
			for(int i=0;i<npnt;++i) {
				in.ignore(80,':');
				for(int n=0;n<ND;++n)
					in >> pnts(i)(n);
			}
			
			/* SIDE INFO */
			for(int i=0;i<nseg;++i) {
				in.ignore(80,':');
				in >> seg(i).pnt(0) >> seg(i).pnt(1);
			}
			
			/* THEN TRI INFO */
			for(int i=0;i<ntri;++i) {
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
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			
			count = 0;
			for(int i=0;i<nebd;++i) {
				in.ignore(80,':');
				in >> temp;
				if (!ebdry(i)) ebdry(i) = getnewedgeobject(temp,bdrymap);
				in.ignore(80,':');
				in >> ebdry(i)->nseg;
				if (!ebdry(i)->maxseg) ebdry(i)->alloc(static_cast<int>(4*grwfac1d*ebdry(i)->nseg));
				else assert(ebdry(i)->nseg < ebdry(i)->maxseg);
				for(int j=0;j<ebdry(i)->nseg;++j) {
					in.ignore(80,':');
					in >> ebdry(i)->seg(j);
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
				*gbl->log << "re-inputting into incompatible mesh object" << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			
			for(int i=0;i<nvbd;++i) {
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
		}

#ifdef libbinio
		case(binary): {
			grd_app = grd_nm + ".bin";
			bin.open(grd_app.c_str());
			if (bin.error()) {
				*gbl->log << "couldn't open grid file: " << grd_app << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			
			bin.setFlag(binio::BigEndian,bin.readInt(1));
			bin.setFlag(binio::FloatIEEE,bin.readInt(1));
			npnt = bin.readInt(sizeof(int));
			nseg = bin.readInt(sizeof(int));
			ntri = bin.readInt(sizeof(int));
			
			if (!initialized) {
				allocate(nseg + (int) (grwfac*nseg));
				nebd = 0;
				nvbd = 0;
			}
			else if (nseg > maxpst) {
				*gbl->log << "mesh is too large" << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			
			/* VRTX INFO */
			for(int i=0;i<npnt;++i) {
				for(n=0;n<ND;++n)
					pnts(i)(n) = bin.readFloat(binio::Double);
				lngth(i) = bin.readFloat(binio::Double);
			}
			
			/* SIDE INFO */
			for(int i=0;i<nseg;++i) {
				seg(i).pnt(0) = bin.readInt(sizeof(int));
				seg(i).pnt(1) = bin.readInt(sizeof(int));
			}
			
			/* THEN TRI INFO */
			for(int i=0;i<ntri;++i) {
				tri(i).pnt(0) = bin.readInt(sizeof(int));
				tri(i).pnt(1) = bin.readInt(sizeof(int));
				tri(i).pnt(2) = bin.readInt(sizeof(int));
			}
			
			/* CREATE TSIDE & STRI */
			createsegtri();
			
			/* SIDE BOUNDARY INFO HEADER */
			in.ignore(80,':');
			int newnsbd1;
			newnsbd1 = bin.readInt(sizeof(int));
			if (nebd == 0) {
				nebd = newnsbd1;
				ebdry.resize(nebd);
				ebdry = 0;
			}
			else if (nebd != newnsbd1) {
				*gbl->log << "reloading incompatible meshes" << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			
			count = 0;
			for(int i=0;i<nebd;++i) {
				temp = bin.readInt(sizeof(int));
				if (!ebdry(i)) ebdry(i) = getnewedgeobject(temp,bdrymap);
				ebdry(i)->nseg = bin.readInt(sizeof(int));
				if (!ebdry(i)->maxseg) ebdry(i)->alloc(static_cast<int>(4*grwfac1d*ebdry(i)->nseg));
				else assert(ebdry(i)->nseg < ebdry(i)->maxseg);
				for(int j=0;j<ebdry(i)->nseg;++j)
					ebdry(i)->seg(j) = bin.readInt(sizeof(int));
			}
			
			/* VERTEX BOUNDARY INFO HEADER */
			int newnvbd1;
			newnvbd1 = bin.readInt(sizeof(int));
			if (nvbd == 0) {
				nvbd = newnvbd1;
				vbdry.resize(nvbd);
				vbdry = 0;
			}
			else if (nvbd != newnvbd1) {
				*gbl->log << "re-inputting into incompatible mesh object" << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			
			for(int i=0;i<nvbd;++i) {
				temp = bin.readInt(sizeof(int));
				if (!vbdry(i)) {
					vbdry(i) = getnewvrtxobject(temp,bdrymap);
					vbdry(i)->alloc(4);
				}
				vbdry(i)->pnt = bin.readInt(sizeof(int));
			}
			bin.close();
			
			break;
		}
#endif
			
		case(netcdf): {
			grd_app = grd_nm +".nc";
			/* Create the file. The NC_CLOBBER parameter tells netCDF to
			 * overwrite this file, if it already exists.*/
			int retval, ncid, dim_id;
			size_t dimreturn;
			
			if ((retval = nc_open(grd_app.c_str(), NC_NOWRITE, &ncid))) ERR(retval);
			
			if ((retval = nc_inq_dimid(ncid, "npnt", &dim_id))) ERR(retval);
			if ((retval = nc_inq_dimlen(ncid, dim_id, &dimreturn))) ERR(retval);
			npnt = dimreturn;
			
			if ((retval = nc_inq_dimid(ncid, "nseg", &dim_id))) ERR(retval);
			if ((retval = nc_inq_dimlen(ncid, dim_id, &dimreturn))) ERR(retval);
			nseg = dimreturn;
			
			if ((retval = nc_inq_dimid(ncid, "ntri", &dim_id))) ERR(retval);
			if ((retval = nc_inq_dimlen(ncid, dim_id, &dimreturn))) ERR(retval);
			ntri = dimreturn;
			
			if (!initialized) {
				allocate(nseg + (int) (grwfac*nseg));
				nebd = 0;
				nvbd = 0;
			}
			else if (nseg > maxpst) {
				*gbl->log << "mesh is too large" << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			
			int var_id;
			if ((retval = nc_inq_varid (ncid, "pnts", &var_id))) ERR(retval);
			size_t index[2];
			/* POINT INFO */
			for(int i=0;i<npnt;++i) {
				index[0]= i;
				for(int n=0;n<ND;++n) {
					index[1] = n;
					nc_get_var1_double(ncid,var_id,index,&pnts(i)(n));
				}
				index[1] = 2;
				nc_get_var1_double(ncid,var_id,index,&lngth(i));
				
			}
			
			/* SEG INFO */
			if ((retval = nc_inq_varid (ncid, "segs", &var_id))) ERR(retval);
			for(int i=0;i<nseg;++i) {
				index[0]= i;
				for(int n=0;n<2;++n) {
					index[1] = n;
					nc_get_var1_int(ncid,var_id,index,&seg(i).pnt(n));
				}
			}
			
			/* TRI INFO */
			if ((retval = nc_inq_varid (ncid, "tris", &var_id))) ERR(retval);
			
			for(int i=0;i<ntri;++i) {
				index[0]= i;
				for(int n=0;n<3;++n) {
					index[1] = n;
					nc_get_var1_int(ncid,var_id,index,&tri(i).pnt(n));
				}
			}
			
			/* CREATE TSIDE & STRI */
			createsegtri();
			
			/* SIDE BOUNDARY INFO HEADER */
			if ((retval = nc_inq_dimid(ncid, "nebd", &dim_id))) ERR(retval);
			if ((retval = nc_inq_dimlen(ncid, dim_id, &dimreturn))) ERR(retval);
			int newnsbd1 = dimreturn;
			
			if (nebd == 0) {
				nebd = newnsbd1;
				ebdry.resize(nebd);
				ebdry = 0;
			}
			else if (nebd != newnsbd1) {
				*gbl->log << "reloading incompatible meshes" << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			
			/* SIDE BOUNDARY INFO HEADER */
			for(int i=0;i<nebd;++i) {
				ostringstream nstr;
				int temp;
				nstr << "edge" << i;
				/* Get the varid of the data variable, based on its name. */
				if ((retval = nc_inq_varid(ncid, nstr.str().c_str(), &var_id))) ERR(retval);
				if ((retval = nc_get_att_int(ncid, var_id, "id", &temp))) ERR(retval);
				if (!ebdry(i)) ebdry(i) = getnewedgeobject(temp,bdrymap);
				
				nstr.str("");
				nstr << "edge" << i << "_nseg";
				if ((retval = nc_inq_dimid(ncid, nstr.str().c_str(), &dim_id))) ERR(retval);
				if ((retval = nc_inq_dimlen(ncid, dim_id, &dimreturn))) ERR(retval);
				ebdry(i)->nseg = dimreturn;
				if (!ebdry(i)->maxseg) ebdry(i)->alloc(static_cast<int>(4*grwfac1d*ebdry(i)->nseg));
				
				/* Get the varid of the data variable, based on its name. */
				if ((retval = nc_get_var_int(ncid, var_id, &ebdry(i)->seg(0)))) ERR(retval);
			}
			
			/* VERTEX BOUNDARY INFO HEADER */
			if ((retval = nc_inq_dimid(ncid, "nvbd", &dim_id))) ERR(retval);
			if ((retval = nc_inq_dimlen(ncid, dim_id, &dimreturn))) ERR(retval);
			int newnvbd1 = dimreturn;
			
			if (nvbd == 0) {
				nvbd = newnvbd1;
				vbdry.resize(nvbd);
				vbdry = 0;
			}
			else if (nvbd != newnvbd1) {
				*gbl->log << "re-inputting into incompatible mesh object" << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			if ((retval = nc_inq_varid(ncid, "vrtx", &var_id))) ERR(retval);
			
			for(int i=0;i<nvbd;++i) {
				index[0] = i;
				index[1] = 0;
				int temp;
				nc_get_var1_int(ncid,var_id,index,&temp);
				if (!vbdry(i)) {
					vbdry(i) = getnewvrtxobject(temp,bdrymap);
					vbdry(i)->alloc(4);
				}
				index[1] = 1;
				nc_get_var1_int(ncid,var_id,index,&vbdry(i)->pnt);
			}
			
			break;
		}
			
			
			
			
		case(text): {
			if (!initialized) {
				*gbl->log << "to read in point positions only must first load mesh structure" << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			
			/* LOAD VERTEX POSITIONS                    */
			grd_app = grd_nm + ".txt";
			in.open(grd_app.c_str());
			if (!in) { *gbl->log << "trouble opening grid" << grd_app << std::endl; sim::abort(__LINE__,__FILE__,gbl->log);}
			
			if(!(in >> temp)) {
				*gbl->log << "1: error in grid " << grd_app << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			if (temp != npnt) {
				*gbl->log << "grid doesn't match point list" << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			
			/* ERROR %lf SHOULD BE FLT */
			for(int i=0;i<npnt;++i) {
				in.ignore(80,':');
				for(int n=0;n<ND;++n)
					in >> pnts(i)(n);
				if (in.fail()) { *gbl->log << "2c: error in grid" << std::endl; sim::abort(__LINE__,__FILE__,gbl->log); }
			}
			in.close();
			
			treeinit();
			
			return;
		}
			
#ifdef CAPRI
		case(BRep): {
			
			/* READ VOLUME & FACE NUMBER FROM FILENAME STRING */
			sscanf(filename,"%d%d",&cpri_vol,&cpri_face);
			status = gi_dGetVolume(cpri_vol,&cpri_nnode,&cpri_nedge,&cpri_nface,&cpri_nbound, &cpri_name);
			*gbl->log << " gi_uGetVolume status =" << status << std::endl;
			if (status != CAPRI_SUCCESS) sim::abort(__LINE__,__FILE__,gbl->log);
			*gbl->log << "  # Edges = " << cpri_nedge << " # Faces = " << cpri_nface << std::endl;
			
			status = gi_dTesselFace(cpri_vol, cpri_face, &ntri, &cpri_tris, &cpri_tric, &npnt, &cpri_points,
															&cpri_ptype, &cpri_pindex, &cpri_uv);
			if (status != CAPRI_SUCCESS) {
				*gbl->log << "gi_dTesselFace status = " << status << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			
			/* ALLOCATE BASIC STORAGE */
			if (!initialized) {
				maxpst = static_cast<int>((grwfac*3*ntri)/2);
				allocate(maxpst);
			}
			else if ((3*ntri)/2 > maxpst) {
				*gbl->log << "mesh is too large" << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
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
			for(int i=0;i<npnt;++i)
				if (cpri_ptype[i] == 0)
					++nvbd;
			vbdry.resize(nvbd);
			
			nvbd = 0;
			/* CREATE VERTEX BOUNDARIES */
			for(int i=0;i<npnt;++i) {
				if (cpri_ptype[i] == 0) {
					vbdry(nvbd) = getnewvrtxobject(cpri_pindex[i],bdrymap);
					vbdry(nvbd)->alloc(4);
					vbdry(nvbd)->pnt = i;
					++nvbd;
					if (nvbd >= MAXVB) {
						*gbl->log << "Too many vertex boundary conditions: increase MAXSB: " << nvbd << std::endl;
						sim::abort(__LINE__,__FILE__,gbl->log);
					}
				}
			}
			
			for(int i=0;i<nseg;++i)
				seg(i).info = 0;
			
			/* COUNT BOUNDARY GROUPS */
			nebd = 0;
			for(int i=0;i<nseg;++i) {
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
							sim::abort(__LINE__,__FILE__,gbl->log);
						}
					}
					for (j = 0; j <nebd;++j) {
						if (seg(i).info == (tri_gbl->intwk(j)&0xFFFF)) {
							++tri_gbl->i2wk(j);
							goto next1b;
						}
					}
					/* NEW SIDE */
					tri_gbl->intwk(nebd) = seg(i).info;
					tri_gbl->i2wk(nebd++) = 1;
				}
			next1b:        continue;
			}
			
			ebdry.resize(nebd);
			for(int i=0;i<nebd;++i) {
				ebdry(i) = getnewedgeobject(tri_gbl->intwk(i),bdrymap);
				ebdry(i)->alloc(static_cast<int>(tri_gbl->i2wk(i)*4*grwfac1d));
				ebdry(i)->nseg = 0;
				tri_gbl->intwk(i) = -1;
				tri_gbl->i2wk(i) = -1;
			}
			
			for(int i=0;i<nseg;++i) {
				if (seg(i).info) {
					for (j = 0; j <nebd;++j) {
						if (seg(i).info == (ebdry(j)->idnum&0xFFFF)) {
							ebdry(j)->seg(ebdry(j)->nseg++) = i;
							goto next1c;
						}
					}
					*gbl->log << "Big error\n";
					sim::abort(__LINE__,__FILE__,gbl->log);
				}
			next1c:      continue;
			}
			
			break;
		}
#endif
			
		case(tecplot): {
			grd_app = grd_nm +".dat";
			in.open(grd_app.c_str());
			if (!in) {
				*gbl->log << "couldn't open tecplot file: " << grd_app << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
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
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			
			/* READ VERTEX DATA */
			for(int i=0;i<npnt;++i) {
				for(int n=0;n<ND;++n)
					in >> pnts(i)(n);
				in.ignore(80,'\n');
				pnt(i).info = -1;
			}
			in.ignore(80,'\n');
			while (in.peek() == '#')
				in.ignore(80,'\n');
			
			for(int i=0;i<ntri;++i) {
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
			for(int i=0;i<nseg;++i)
				if (seg(i).tri(1) < 0) ++count;
			
			nebd = 1;
			ebdry.resize(1);
			ebdry(0) = getnewedgeobject(1,bdrymap);
			ebdry(0)->alloc(static_cast<int>(4*grwfac1d*count));
			ebdry(0)->nseg = count;
			count = 0;
			for(int i=0;i<nseg;++i)
				if (seg(i).tri(1) < 0) ebdry(0)->seg(count++) = i;
			
			nvbd = 0;
			in.close();
			
			break;
		}
		case(boundary): {
			/* LOAD BOUNDARY INFORMATION */
			grd_app = grd_nm +".d";
			in.open(grd_app.c_str());
			if (!in) {
				*gbl->log << "couldn't open " << grd_app << "for reading\n";
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			in >> npnt;
            for (int i = 0; i < npnt+1; ++i) {
                in.ignore(numeric_limits<streamsize>::max(), in.widen('\n'));
            }
            in >> nseg;
            maxpst = MAX(static_cast<int>(grwfac*nseg*nseg),static_cast<int>(grwfac*npnt));

            in.close();
            in.open(grd_app.c_str());
            in >> npnt;
			
            nseg = 0;

            allocate(maxpst);
			initialized = 1;
			std::string symbolic_string;
			istringstream data;
			for(int i=0;i<npnt;++i) {
				in.ignore(80,':');
				for(int n=0;n<ND;++n) {
					in >> symbolic_string;
					data.str(symbolic_string);
					if (!(data >> pnts(i)(n)) || !data.eof()) {
						bdrymap["mesh_formula"] = symbolic_string;
						if (!bdrymap.get("mesh_formula",pnts(i)(n))) {
							*gbl->log << "couldn't read formula " << symbolic_string << "from .d file\n";
							sim::abort(__LINE__,__FILE__,gbl->log);
						}
					}
					data.clear();
				}
				in >> symbolic_string;
				data.str(symbolic_string);
				if (!(data >> lngth(i)) || !data.eof()) {
					bdrymap["mesh_formula"] = symbolic_string;
					if (!bdrymap.get("mesh_formula",lngth(i))) {
						*gbl->log << "couldn't read formula " << symbolic_string << "from .d file\n";
						sim::abort(__LINE__,__FILE__,gbl->log);
					}
				}
				data.clear();
				in >> pnt(i).info;
			}
			
			
			/* COUNT VERTEX BOUNDARY GROUPS  */
			nvbd = 0;
			for(int i=0;i<npnt;++i)
				if (pnt(i).info)
					++nvbd;
			vbdry.resize(nvbd);
			
			nvbd = 0;
			for(int i=0;i<npnt;++i) {
				if (pnt(i).info) {
					/* NEW VRTX B.C. */
					vbdry(nvbd) = getnewvrtxobject(pnt(i).info,bdrymap);
					vbdry(nvbd)->alloc(4);
					vbdry(nvbd)->pnt = i;
					++nvbd;
				}
			}
			
			in >> nseg;
			
			for(int i=0;i<nseg;++i) {
				in.ignore(80,':');
				in >> seg(i).pnt(0) >> seg(i).pnt(1) >> seg(i).info;
			}
			
			/* COUNT BOUNDARY GROUPS */
			nebd = 0;
			for(int i=0;i<nseg;++i) {
				if (seg(i).info) {
					for (int j = 0; j <nebd;++j) {
						if (seg(i).info == tri_gbl->intwk(j)) {
							++tri_gbl->i2wk(j);
							goto bdnext1;
						}
					}
					/* NEW SIDE */
					tri_gbl->intwk(nebd) = seg(i).info;
					tri_gbl->i2wk(nebd++) = 1;
				}
				else {
					*gbl->log << "All sides should be boundary sides in a .d file" << std::endl;
					sim::abort(__LINE__,__FILE__,gbl->log);
				}
			bdnext1:        continue;
			}
			
			
			ebdry.resize(nebd);
			for(int i=0;i<nebd;++i) {
				ebdry(i) = getnewedgeobject(tri_gbl->intwk(i),bdrymap);
				ebdry(i)->alloc(static_cast<int>(tri_gbl->i2wk(i)*4*grwfac1d));
				ebdry(i)->nseg = 0;
				tri_gbl->intwk(i) = -1;
				tri_gbl->i2wk(i) = -1;
			}
			
			for(int i=0;i<nseg;++i) {
				if (seg(i).info) {
					for (int j = 0; j <nebd;++j) {
						if (seg(i).info == ebdry(j)->idnum) {
							ebdry(j)->seg(ebdry(j)->nseg++) = i;
							goto bdnext1a;
						}
					}
					*gbl->log << "Big error\n";
					sim::abort(__LINE__,__FILE__,gbl->log);
				}
			bdnext1a:      continue;
			}
			
			for(int i=0;i<nseg;++i)
				tri_gbl->i2wk_lst1(i) = i+1;
			
			ntri = 0;
			triangulate(nseg);
			
			int bpnt = 0;
			for (int i=0;i<nebd;++i)
				bpnt += ebdry(i)->nseg;
			
			interior_pts = npnt - bpnt;
			npnt = bpnt;
			
			in.close();
			
			break;
		}
			
		case(vlength): {
			grd_app = grd_nm +".lngth";
			in.open(grd_app.c_str());
			if (!in) {
				*gbl->log << "couldn't open vlength input file" << grd_app  << endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			
			for(int i=0;i<npnt;++i)
				in >> lngth(i);
			
			return;
		}
			
		default: {
			*gbl->log << "That filetype is not supported" << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
	}

	for(int i=0;i<nebd;++i) {
		/* CREATES NEW BOUNDARY FOR DISCONNECTED SEGMENTS OF SAME TYPE */
		ebdry(i)->reorder();
	}

	/* FIND ENDPOINT MATCHES */
	for(int i=0;i<nvbd;++i) {
		/* Find two connecting boundary sides */
		for(int j=0;j<nebd;++j) {
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
	
	if (filetype == boundary) {
		/* Check if there are additional points to insert */
		/* This only happens for .d files */
		int pnear,tind,err;
		bool found;
		for(int i=0;i<interior_pts;++i) {
			qtree.addpt(npnt);
			qtree.nearpt(npnt,pnear);
			found = findtri(pnts(npnt),pnear,tind);
			assert(found);
			err = insert(npnt,tind);
			++npnt;
		}
		cnt_nbor();
	}


	grd_app = grd_nm + ".lngth";
	in.open(grd_app.c_str());
	if (in) {
		for(int i=0;i<npnt;++i) in >> lngth(i);
		in.close();
	}
	else if (filetype != boundary && filetype != binary) initlngth();


	tri_mesh::setinfo();  // Don't call virtual because other objects must be set-up first
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


