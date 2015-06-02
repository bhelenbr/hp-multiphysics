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
#include "surface2.h"
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
		static const int ntypes = 8;
		enum ids {unknown=-1,surface_inflow,surface_outflow,surface_inflow2,surface_outflow2,inflow,pressure,hybrid_slave_point,hybrid_point};
		const static char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i) 
				if (!strcmp(nin,names[i])) return(i);
			return(unknown);
		}
};

const char tri_hp_ins_vtype::names[ntypes][40] = {"surface_inflow","surface_outflow","surface_inflow2","surface_outflow2","inflow","pressure","hybrid_slave_point","hybrid_point"};

hp_vrtx_bdry* tri_hp_ins::getnewvrtxobject(int bnum, std::string name) {
	std::string keyword,val;
	std::istringstream data;
	int type;          
	hp_vrtx_bdry *temp;  

	type = tri_hp_ins_vtype::getid(name.c_str());
	switch(type) {
		case tri_hp_ins_vtype::surface_inflow: {
			temp = new surface_fixed_pt(*this,*vbdry(bnum));
			break;
		}
		case tri_hp_ins_vtype::surface_outflow: {
			temp = new surface_outflow(*this,*vbdry(bnum));
			break;
		}
		case tri_hp_ins_vtype::surface_inflow2: {
			temp = new hp_deformable_fixed_pnt(*this,*vbdry(bnum));
			break;
		}
		case tri_hp_ins_vtype::surface_outflow2: {
			temp = new surface_outflow2(*this,*vbdry(bnum));
			break;
		}
		case tri_hp_ins_vtype::inflow: {
			temp = new inflow_pt(*this,*vbdry(bnum));
			break;
		}
		case tri_hp_ins_vtype::pressure: {
			temp = new pressure_pt(*this,*vbdry(bnum));
			break;
		}
		case tri_hp_ins_vtype::hybrid_slave_point: {
			temp = new hybrid_slave_pt(*this,*vbdry(bnum));
			break;
		}
		case tri_hp_ins_vtype::hybrid_point: {
			temp = new hybrid_pt(*this,*vbdry(bnum));
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
class tri_hp_ins_stype {
	public:
		static const int ntypes = 12;
		enum ids {unknown=-1,inflow,outflow,characteristic,euler,
			symmetry,applied_stress,surface,surface_slave,surface2,force_coupling,friction_slip,actuator_disc};
		static const char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i)
				if (!strcmp(nin,names[i])) return(i);
			return(-1);
		}
};

const char tri_hp_ins_stype::names[ntypes][40] = {"inflow","outflow","characteristic","euler",
    "symmetry","applied_stress","surface","surface_slave","surface2","force_coupling","friction_slip","actuator_disc"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_edge_bdry* tri_hp_ins::getnewsideobject(int bnum, std::string name) {
	std::string keyword,val;
	std::istringstream data;
	int type;          
	hp_edge_bdry *temp;  


	type = tri_hp_ins_stype::getid(name.c_str());
	switch(type) {
		case tri_hp_ins_stype::inflow: {
			temp = new inflow(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_ins_stype::outflow: {
			temp = new generic(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_ins_stype::characteristic: {
			temp = new characteristic(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_ins_stype::euler: {
			temp = new euler(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_ins_stype::symmetry: {
			temp = new symmetry(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_ins_stype::applied_stress: {
			temp = new applied_stress(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_ins_stype::surface: {
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
		case tri_hp_ins_stype::surface_slave: {
			temp = new surface_slave(*this,*ebdry(bnum));
			dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))->physics = temp;
			break;
		}
		case tri_hp_ins_stype::surface2: {
			if (dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))) {
				temp = new surface2(*this,*ebdry(bnum));
				dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))->physics = temp;
			}
			else {
				std::cerr << "use coupled physics for surface boundary" << std::endl;
				sim::abort(__LINE__,__FILE__,&std::cerr);
				assert(0);
			}
			break;
		}
		case tri_hp_ins_stype::force_coupling: {
			temp = new force_coupling(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_ins_stype::friction_slip: {
			temp = new friction_slip(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_ins_stype::actuator_disc: {
			temp = new actuator_disc(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_ins_stype::unknown: {
			return(tri_hp::getnewsideobject(bnum,name));
		}
	}    
	gbl->ebdry_gbls(bnum) = temp->create_global_structure();

	return(temp);
}


