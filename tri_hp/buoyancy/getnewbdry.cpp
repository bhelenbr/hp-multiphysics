//
//  getnewbdry.cpp
//  tri_hp
//
//  Created by Brian Helenbrook on 3/22/15.
//
//

#include <stdio.h>
#include "bdry_buoyancy.h"
#include <tri_boundary.h>

using namespace bdry_buoyancy;
#include "bdry_buoyancy.h"
#include "melt_buoyancy.h"

/** \brief Helper object for vrtx_bdry
 *
 * \ingroup boundary
 * Contains list of all vrtx_bdys's by name
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tri_hp_buoyancy_vtype {
public:
	static const int ntypes = 2;
	enum ids {unknown=-1,melt_facet_pt2,triple_junction};
	const static char names[ntypes][40];
	static int getid(const char *nin) {
		for(int i=0;i<ntypes;++i)
			if (!strcmp(nin,names[i])) return(i);
		return(unknown);
	}
};

const char tri_hp_buoyancy_vtype::names[ntypes][40] = {"melt_facet_pt2","triple_junction"};

hp_vrtx_bdry* tri_hp_buoyancy::getnewvrtxobject(int bnum, std::string name) {
	std::string keyword,val;
	std::istringstream data;
	int type;
	hp_vrtx_bdry *temp;
	
	type = tri_hp_buoyancy_vtype::getid(name.c_str());
	switch(type) {
		case tri_hp_buoyancy_vtype::melt_facet_pt2: {
			temp = new melt_facet_pt2(*this,*vbdry(bnum));
			break;
		}
		case tri_hp_buoyancy_vtype::triple_junction: {
			temp = new triple_junction(*this,*vbdry(bnum));
			break;
		}
		case tri_hp_buoyancy_vtype::unknown: {
			return(tri_hp_ins::getnewvrtxobject(bnum,name));
		}
	}
	gbl->vbdry_gbls(bnum) = temp->create_global_structure();
	return(temp);
}



class tri_hp_buoyancy_etype {
public:
	static const int ntypes = 4;
	enum ids {unknown=-1,surface_marangoni,solid_fluid,characteristic,melt2};
	static const char names[ntypes][40];
	static int getid(const char *nin) {
		for(int i=0;i<ntypes;++i)
			if (!strcmp(nin,names[i])) return(i);
		return(-1);
	}
};

const char tri_hp_buoyancy_etype::names[ntypes][40] = {"surface_marangoni","solid_fluid","characteristic","melt2"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_edge_bdry* tri_hp_buoyancy::getnewedgeobject(int bnum, std::string name) {
	std::string keyword,val;
	std::istringstream data;
	int type;
	hp_edge_bdry *temp;
	
	type = tri_hp_buoyancy_etype::getid(name.c_str());
	switch(type) {
		case tri_hp_buoyancy_etype::surface_marangoni: {
			if (dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))) {
				temp = new surface_marangoni(*this,*ebdry(bnum));
				dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))->physics = temp;
			}
			else {
				std::cerr << "use coupled physics for surface boundary" << std::endl;
				sim::abort(__LINE__,__FILE__,&std::cerr);
				assert(0);
			}
			break;
		}

		case tri_hp_buoyancy_etype::solid_fluid: {
			if (dynamic_cast<eboundary_with_geometry<ecomm,symbolic_shape<tri_mesh::ND> > *>(ebdry(bnum))) {
				temp = new solid_fluid(*this,*ebdry(bnum));
			}
			else {
				std::cerr << "use symbolic_comm for solid_liquid boundary" << std::endl;
				sim::abort(__LINE__,__FILE__,&std::cerr);
				assert(0);
			}
			break;
		}
		case tri_hp_buoyancy_etype::characteristic: {
			temp = new characteristic(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_buoyancy_etype::melt2: {
			if (dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))) {
				temp = new melt_buoyancy(*this,*ebdry(bnum));
				dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))->physics = temp;
			}
			else {
				std::cerr << "use coupled physics for melt_kinetics boundary" << std::endl;
				sim::abort(__LINE__,__FILE__,&std::cerr);
				assert(0);
			}
			break;
		}
			
		default: {
			return(tri_hp_ins::getnewedgeobject(bnum,name));
			break;
		}
	}
	gbl->ebdry_gbls(bnum) = temp->create_global_structure();
	
	return(temp);
}
