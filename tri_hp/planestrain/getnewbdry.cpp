/*
 *  getnewbdry.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 1/12/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
#include "tri_hp_ps.h"
#include "bdry_ps.h"
#include <input_map.h>

using namespace bdry_ps;

/** \brief Helper object for edge_bdry 
 *
 * \ingroup boundary
 * Contains list of all side_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tri_hp_ps_stype {
	public:
		static const int ntypes = 3;
		enum ids {unknown=-1,dirichlet,neumann,friction_wall};
		static const char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i)
				if (!strcmp(nin,names[i])) return(i);
			return(unknown);
		}
};

const char tri_hp_ps_stype::names[ntypes][40] = {"dirichlet","neumann","friction_wall"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_edge_bdry* tri_hp_ps::getnewsideobject(int bnum, std::string name) {
	std::string keyword,val;
	std::istringstream data;
	int type;          
	hp_edge_bdry *temp;  

	type = tri_hp_ps_stype::getid(name.c_str());
	switch(type) {
		case tri_hp_ps_stype::dirichlet: {
			temp = new dirichlet(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_ps_stype::neumann: {
			temp = new neumann(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_ps_stype::friction_wall: {
			temp = new friction_wall(*this,*ebdry(bnum));  // FIXME NOT WORKING YET
			break;
		}
		default: {
			return(tri_hp::getnewsideobject(bnum,name));
		}
	}
	gbl->ebdry_gbls(bnum) = temp->create_global_structure();

	return(temp);
}


