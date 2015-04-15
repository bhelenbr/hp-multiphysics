//
//  hp_getnewbdry.cpp
//  tri_hp
//
//  Created by Brian Helenbrook on 3/8/15.
//
//

#include <stdio.h>
#include "hp_boundary.h"
#include "hp_coupled_boundary.h"

class tri_hp_vtype {
public:
	static const int ntypes = 3;
	enum ids {unknown=-1,plain,hp_deformable_fixed_pnt,hp_deformable_free_pnt};
	const static char names[ntypes][40];
	static int getid(const char *nin) {
		for(int i=0;i<ntypes;++i)
			if (!strcmp(nin,names[i])) return(i);
		return(unknown);
	}
};

const char tri_hp_vtype::names[ntypes][40] = {"plain","hp_deformable_fixed_pnt","hp_deformable_free_pnt"};

hp_vrtx_bdry* tri_hp::getnewvrtxobject(int bnum, input_map &bdrydata) {
	std::string keyword,val;
	std::istringstream data;
	int type;
	hp_vrtx_bdry *temp;
	
	keyword =  vbdry(bnum)->idprefix + "_hp_type";
	if (!bdrydata.get(keyword,val)) {
		*gbl->log << "missing vertex type:" << keyword << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
	else {
		type = tri_hp_vtype::getid(val.c_str());
	}
	
	
	switch(type) {
		case tri_hp_vtype::hp_deformable_fixed_pnt: {
			temp = new hp_deformable_fixed_pnt(*this,*vbdry(bnum));
			break;
		}
		case tri_hp_vtype::hp_deformable_free_pnt: {
			temp = new hp_deformable_free_pnt(*this,*vbdry(bnum));
			break;
		}
		default: {
			temp = new hp_vrtx_bdry(*this,*vbdry(bnum));
		}
	}
	gbl->vbdry_gbls(bnum) = temp->create_global_structure();
	return(temp);
}

class tri_hp_stype {
public:
	static const int ntypes = 3;
	enum ids {unknown=-1,plain,symbolic,symbolic_with_integration_by_parts};
	static const char names[ntypes][40];
	static int getid(const char *nin) {
		for(int i=0;i<ntypes;++i)
			if (!strcmp(nin,names[i])) return(i);
		return(-1);
	}
};

const char tri_hp_stype::names[ntypes][40] = {"plain","symbolic","symbolic_ibp"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_edge_bdry* tri_hp::getnewsideobject(int bnum, input_map &bdrydata) {
	std::string keyword,val;
	std::istringstream data;
	int type;
	hp_edge_bdry *temp;
	
	
	keyword =  ebdry(bnum)->idprefix + "_hp_type";
	if (!bdrydata.get(keyword,val)) {
		*gbl->log << "missing side type:" << keyword << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
	else {
		type = tri_hp_stype::getid(val.c_str());
	}
	
	switch(type) {
		case tri_hp_stype::symbolic: {
			temp = new symbolic(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_stype::symbolic_with_integration_by_parts: {
			temp = new symbolic_with_integration_by_parts(*this,*ebdry(bnum));
			break;
		}
		default: {
			temp = new hp_edge_bdry(*this,*ebdry(bnum));
		}
	}
	gbl->ebdry_gbls(bnum) = temp->create_global_structure();
	
	return(temp);
}
