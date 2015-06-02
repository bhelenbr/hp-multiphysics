/*
 *  getnewswirl_bdry.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 1/12/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
#include <input_map.h>

#include "tri_hp_swirl.h"
#include "bdry_swirl.h"

using namespace bdry_swirl;

/** \brief Helper object for edge_bdry 
 *
 * \ingroup boundary
 * Contains list of all side_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tri_hp_swirl_stype {
	public:
		static const int ntypes = 1;
		enum ids {unknown=-1,symmetry};
		static const char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i)
				if (!strcmp(nin,names[i])) return(i);
			return(unknown);
		}
};

const char tri_hp_swirl_stype::names[ntypes][40] = {"symmetry"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_edge_bdry* tri_hp_swirl::getnewsideobject(int bnum, std::string name) {
	std::string keyword,val;
	std::istringstream data;
	int type;          
	hp_edge_bdry *temp;  

	type = tri_hp_swirl_stype::getid(name.c_str());
	switch(type) {
		case tri_hp_swirl_stype::symmetry: {
			temp = new symmetry(*this,*ebdry(bnum));
			break;
		}	
		default: {
			return(tri_hp_ins::getnewsideobject(bnum,name));
		}
	}
	gbl->ebdry_gbls(bnum) = temp->create_global_structure();

	return(temp);
}


