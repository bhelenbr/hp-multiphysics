/*
 *  getnewins_bdry.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 1/12/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
#include "tri_hp_ins.h"
#include "bdry_ins.h"
#include "surface.h"
#include <input_map.h>
#include <tri_boundary.h>

using namespace bdry_ins;

/** \brief Helper object for vrtx_bdry 
 *
 * \ingroup boundary
 * Contains list of all vrtx_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tri_hp_ins_vtype {
	public:
		static const int ntypes = 2;
		enum ids {unknown=-1,surface_inflow,surface_outflow};
		const static char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i) 
				if (!strcmp(nin,names[i])) return(i);
			return(unknown);
		}
};

const char tri_hp_ins_vtype::names[ntypes][40] = {"surface_inflow","surface_outflow"};

hp_vrtx_bdry* tri_hp_ins::getnewvrtxobject(int bnum, std::string name) {
	std::string keyword,val;
	std::istringstream data;
	int type;          
	hp_vrtx_bdry *temp;  

	type = tri_hp_ins_vtype::getid(name.c_str());
	switch(type) {
		case tri_hp_ins_vtype::surface_inflow: {
			temp = new hp_deformable_fixed_pnt(*this,*vbdry(bnum));
			break;
		}
		case tri_hp_ins_vtype::surface_outflow: {
			temp = new surface_outflow(*this,*vbdry(bnum));
			break;
		}
		case tri_hp_ins_vtype::unknown: {
			return(tri_hp::getnewvrtxobject(bnum,name));
		}
	} 
	gbl->vbdry_gbls(bnum) = temp->create_global_structure();
	return(temp);
}


/** \brief Helper object for edge_bdry 
 *
 * \ingroup boundary
 * Contains list of all side_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tri_hp_ins_etype {
	public:
		static const int ntypes = 12;
		enum ids {unknown=-1,inflow,outflow,characteristic,euler,
			symmetry,applied_stress,surface,force_coupling,friction_slip,actuator_disc};
		static const char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i)
				if (!strcmp(nin,names[i])) return(i);
			return(-1);
		}
};

const char tri_hp_ins_etype::names[ntypes][40] = {"inflow","outflow","characteristic","euler",
    "symmetry","applied_stress","surface","force_coupling","friction_slip","actuator_disc"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_edge_bdry* tri_hp_ins::getnewedgeobject(int bnum, std::string name) {
	std::string keyword,val;
	std::istringstream data;
	int type;          
	hp_edge_bdry *temp;  


	type = tri_hp_ins_etype::getid(name.c_str());
	switch(type) {
		case tri_hp_ins_etype::inflow: {
			temp = new inflow(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_ins_etype::outflow: {
			temp = new generic(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_ins_etype::characteristic: {
			temp = new characteristic(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_ins_etype::euler: {
			temp = new euler(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_ins_etype::symmetry: {
			temp = new symmetry(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_ins_etype::applied_stress: {
			temp = new applied_stress(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_ins_etype::surface: {
			if (dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))) {
				temp = new surface(*this,*ebdry(bnum));
				dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))->physics = temp;
			}
			else {
				std::cerr << "use coupled physics for surface boundary" << std::endl;
				sim::abort(__LINE__,__FILE__,&std::cerr);
				assert(0);
			}
			break;
		}
		case tri_hp_ins_etype::force_coupling: {
			temp = new force_coupling(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_ins_etype::friction_slip: {
			temp = new friction_slip(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_ins_etype::actuator_disc: {
			temp = new actuator_disc(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_ins_etype::unknown: {
			return(tri_hp::getnewedgeobject(bnum,name));
		}
	}    
	gbl->ebdry_gbls(bnum) = temp->create_global_structure();

	return(temp);
}


