/*
 *  getnewswe_bdry.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 1/12/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
#include "tri_hp_swe.h"
#include "bdry_swe.h"
#include <input_map.h>

using namespace bdry_swe;

/** \brief Helper object for edge_bdry 
 *
 * \ingroup boundary
 * Contains list of all side_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tri_hp_swe_etype {
	public:
		static const int ntypes = 2;
		enum ids {unknown=-1,wall,characteristic};
		static const char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i)
				if (!strcmp(nin,names[i])) return(i);
			return(unknown);
		}
};

const char tri_hp_swe_etype::names[ntypes][40] = {"wall","characteristic"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_edge_bdry* tri_hp_swe::getnewedgeobject(int bnum, std::string name) {
	std::string keyword,val;
	std::istringstream data;
	int type;          
	hp_edge_bdry *temp;  

	type = tri_hp_swe_etype::getid(name.c_str());
	switch(type) {
		case tri_hp_swe_etype::wall: {
			temp = new wall(*this,*ebdry(bnum));
			break;
		}	
		case tri_hp_swe_etype::characteristic: {
			temp = new characteristic(*this,*ebdry(bnum));
			break;
		}	
		default: {
			return(tri_hp_ins::getnewedgeobject(bnum,name));
			break;
		}
	}

	return(temp);
}


