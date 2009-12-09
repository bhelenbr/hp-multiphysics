/*
 *  getnewr_bdry.cpp
 *  mblock
 *
 *  Created by helenbrk on Wed Oct 24 2004.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "r_tri_mesh.h"
#include "r_tri_boundary.h"

class r_stype {
	public:
		static const int ntypes = 6;
		enum ids {free=1, fixed, fixed_angled, translating, oscillating, deforming};
		const static char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i)
				if (!strcmp(nin,names[i])) return(i+1);
			return(-1);
		}
};

const char r_stype::names[ntypes][40] = {"free", "fixed", "fixed_angled", "translating", "oscillating", "deforming"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
r_side_bdry* r_tri_mesh::getnewedgeobject(int bnum, input_map& in_map) {
	std::string typ_str;
	int type;
	r_side_bdry *temp;

	if (in_map.get(ebdry(bnum)->idprefix + "_r_type",typ_str)) {
		type = r_stype::getid(typ_str.c_str());
		if (type < 0)  {
			*gbl->log << ebdry(bnum)->idprefix << "_r_type: unknown type " << typ_str << std::endl;
			exit(1);
		}
	}
	else {
		/* SOME DEFAULTS FOR VARIOUS BOUNDARY TYPES */
		if (ebdry(bnum)->mytype == "comm") {
			type = r_stype::free;
		} else if (ebdry(bnum)->mytype == "partition") {
			type = r_stype::free;
		} else if (ebdry(bnum)->mytype == "prdc") {
			type = r_stype::fixed;
			int dir;
			in_map.getwdefault(ebdry(bnum)->idprefix + "_dir",dir,0);
			if (dir == 0)
				in_map[ebdry(bnum)->idprefix+"_r_dir"] = "0 0";
			else
				in_map[ebdry(bnum)->idprefix+"_r_dir"] = "1 1";
		}
		else {
			type = r_stype::fixed;
		}
		*gbl->log << '#' << ebdry(bnum)->idprefix << "_r_type: " <<  r_stype::names[type-1] << std::endl;
	}

	switch(type) {
		case r_stype::free: {
			temp = new r_side_bdry(*this,*ebdry(bnum));
			break;
		}
		case r_stype::fixed: {
			temp = new r_fixed(*this,*ebdry(bnum));
			break;
		}
		case r_stype::fixed_angled: {
			temp = new r_fixed_angled(*this,*ebdry(bnum));
			break;
		}
		case r_stype::translating: {
			temp = new r_translating(*this,*ebdry(bnum));
			break;
		}
		case r_stype::oscillating: {
			temp = new r_oscillating(*this,*ebdry(bnum));
			break;
		}
		case r_stype::deforming: {
			temp = new r_deforming(*this,*ebdry(bnum));
			break;
		}
		default: {
			temp = new r_fixed(*this,*ebdry(bnum));
			std::cout << "Don't know this r_side_bdry type\n";
		}
	}

	temp->input(in_map);

	return(temp);
}

class r_vtype {
	public:
		static const int ntypes = 3;
		enum ids {free=1, fixed, moving};
		const static char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i)
				if (!strcmp(nin,names[i])) return(i+1);
			return(-1);
		}
};

const char r_vtype::names[ntypes][40] = {"free", "fixed", "moving"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
r_vrtx_bdry* r_tri_mesh::getnewvrtxobject(int bnum, input_map& in_map) {
	std::string typ_str;
	int type;
	r_vrtx_bdry *temp;

	in_map.getwdefault(vbdry(bnum)->idprefix + "_r_type",typ_str,std::string("free"));
	type = r_vtype::getid(typ_str.c_str());
	if (type < 0)  {
		*gbl->log << vbdry(bnum)->idprefix << "_r_type: unknown type " << typ_str << std::endl;
		exit(1);
	}

	switch(type) {
		case r_vtype::free: {
			temp = new r_vrtx_bdry(*this,*vbdry(bnum));
			break;
		}
		case r_vtype::fixed: {
			temp = new r_vfixed(*this,*vbdry(bnum));
			break;
		}
		case r_vtype::moving: {
			temp = new r_vmoving(*this,*vbdry(bnum));
			break;
		}
		default: {
			temp = new r_vrtx_bdry(*this,*vbdry(bnum));
			std::cout << "Don't know this r_vrtx_bdry type\n";

		}
	}

	temp->input(in_map);

	return(temp);
}
