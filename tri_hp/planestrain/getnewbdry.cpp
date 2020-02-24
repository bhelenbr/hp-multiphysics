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
class tri_hp_ps_etype {
	public:
		static const int ntypes = 4;
		enum ids {unknown=-1,dirichlet,neumann,friction_wall,curve_edges};
		static const char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i)
				if (!strcmp(nin,names[i])) return(i);
			return(unknown);
		}
};

const char tri_hp_ps_etype::names[ntypes][40] = {"dirichlet","neumann","friction_wall","curve_edges"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_edge_bdry* tri_hp_ps::getnewedgeobject(int bnum, std::string name) {
	std::string keyword,val;
	std::istringstream data;
	int type;          
	hp_edge_bdry *temp;  

	type = tri_hp_ps_etype::getid(name.c_str());
	switch(type) {
		case tri_hp_ps_etype::dirichlet: {
			temp = new dirichlet(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_ps_etype::neumann: {
			temp = new neumann(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_ps_etype::friction_wall: {
			temp = new friction_wall(*this,*ebdry(bnum));  // FIXME NOT WORKING YET
			break;
		}
        case tri_hp_ps_etype::curve_edges: {
            temp = new curve_edges(*this,*ebdry(bnum));
            break;
        }
		default: {
			return(tri_hp::getnewedgeobject(bnum,name));
		}
	}
	gbl->ebdry_gbls(bnum) = temp->create_global_structure();

	return(temp);
}


