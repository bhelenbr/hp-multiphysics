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
	static const int ntypes = 4;
	enum ids {unknown=-1,melt_end,melt_inflow,melt_facet_pt,melt_facet_pt2};
	const static char names[ntypes][40];
	static int getid(const char *nin) {
		for(int i=0;i<ntypes;++i)
			if (!strcmp(nin,names[i])) return(i);
		return(unknown);
	}
};

const char tri_hp_buoyancy_vtype::names[ntypes][40] = {"melt_end","melt_inflow","melt_facet_pt","melt_facet_pt2"};

hp_vrtx_bdry* tri_hp_buoyancy::getnewvrtxobject(int bnum, input_map &bdrydata) {
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
		type = tri_hp_buoyancy_vtype::getid(val.c_str());
	}
	
	
	switch(type) {
		case tri_hp_buoyancy_vtype::melt_end: {
			temp = new melt_end_pt(*this,*vbdry(bnum));
			break;
		}
		case tri_hp_buoyancy_vtype::melt_inflow: {
			temp = new melt_inflow_pt(*this,*vbdry(bnum));
			break;
		}
		case tri_hp_buoyancy_vtype::melt_facet_pt: {
			temp = new melt_facet_pt(*this,*vbdry(bnum));
			break;
		}
		case tri_hp_buoyancy_vtype::melt_facet_pt2: {
			temp = new melt_facet_pt2(*this,*vbdry(bnum));
			break;
		}
		case tri_hp_buoyancy_vtype::unknown: {
			return(tri_hp_ins::getnewvrtxobject(bnum,bdrydata));
		}
	}
	gbl->vbdry_gbls(bnum) = temp->create_global_structure();
	return(temp);
}



class tri_hp_buoyancy_stype {
public:
	static const int ntypes = 7;
	enum ids {unknown=-1,surface,surface_marangoni,melt,melt_kinetics,solid_fluid,characteristic,melt2};
	static const char names[ntypes][40];
	static int getid(const char *nin) {
		for(int i=0;i<ntypes;++i)
			if (!strcmp(nin,names[i])) return(i);
		return(-1);
	}
};

const char tri_hp_buoyancy_stype::names[ntypes][40] = {"surface","surface_marangoni","melt","melt_kinetics","solid_fluid","characteristic","melt2"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_edge_bdry* tri_hp_buoyancy::getnewsideobject(int bnum, input_map &bdrydata) {
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
		type = tri_hp_buoyancy_stype::getid(val.c_str());
	}
	
	switch(type) {
		case tri_hp_buoyancy_stype::surface: {
			if (dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))) {
				temp = new surface9(*this,*ebdry(bnum));
				dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))->physics = temp;
			}
			else {
				std::cerr << "use coupled physics for surface boundary" << std::endl;
				sim::abort(__LINE__,__FILE__,&std::cerr);
				assert(0);
			}
			break;
		}
		case tri_hp_buoyancy_stype::surface_marangoni: {
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
		case tri_hp_buoyancy_stype::melt: {
			if (dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))) {
				temp = new melt(*this,*ebdry(bnum));
				dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))->physics = temp;
			}
			else {
				std::cerr << "use coupled physics for melt boundary" << std::endl;
				sim::abort(__LINE__,__FILE__,&std::cerr);
				assert(0);
			}
			break;
		}
		case tri_hp_buoyancy_stype::melt_kinetics: {
			if (dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))) {
				temp = new melt_kinetics(*this,*ebdry(bnum));
				dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))->physics = temp;
			}
			else {
				std::cerr << "use coupled physics for melt_kinetics boundary" << std::endl;
				sim::abort(__LINE__,__FILE__,&std::cerr);
				assert(0);
			}
			break;
		}
		case tri_hp_buoyancy_stype::solid_fluid: {
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
		case tri_hp_buoyancy_stype::characteristic: {
			temp = new characteristic(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_buoyancy_stype::melt2: {
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
			return(tri_hp_ins::getnewsideobject(bnum,bdrydata));
			break;
		}
	}
	gbl->ebdry_gbls(bnum) = temp->create_global_structure();
	
	return(temp);
}
