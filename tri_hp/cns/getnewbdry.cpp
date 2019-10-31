/*
 *  getnewins_bdry.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 1/12/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
#include "tri_hp_cns.h"
#include "bdry_cns.h"
#include "shock.h"
#include <input_map.h>
#include <tri_boundary.h>

using namespace bdry_cns;

/** \brief Helper object for vrtx_bdry 
 *
 * \ingroup boundary
 * Contains list of all vrtx_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tri_hp_cns_vtype {
	public:
		static const int ntypes = 7;
		enum ids {unknown=-1,surface_inflow,surface_periodic,surface_outflow,surface_outflow_planar,
		inflow,hybrid_slave_point,hybrid_point};
		const static char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i) 
				if (!strcmp(nin,names[i])) return(i);
			return(unknown);
		}
};

const char tri_hp_cns_vtype::names[ntypes][40] = {"surface_inflow","surface_periodic","surface_outflow","surface_outflow_planar",
	"inflow","hybrid_slave_point","hybrid_point"};

hp_vrtx_bdry* tri_hp_cns::getnewvrtxobject(int bnum, std::string name) {
	std::string keyword,val;
	std::istringstream data;
	int type;          
	hp_vrtx_bdry *temp;  

	type = tri_hp_cns_vtype::getid(name.c_str());
	switch(type) {
//		case tri_hp_cns_vtype::surface_inflow: {
//			temp = new surface_fixed_pt(*this,*vbdry(bnum));
//			break;
//		}
//		case tri_hp_cns_vtype::surface_periodic: {
//			temp = new surface_periodic_pt(*this,*vbdry(bnum));
//			break;
//		}
//		case tri_hp_cns_vtype::surface_outflow: {
//			temp = new surface_outflow_endpt(*this,*vbdry(bnum));
//			break;
//		}
//		case tri_hp_cns_vtype::surface_outflow_planar: {
//			temp = new surface_outflow_planar(*this,*vbdry(bnum));
//			break;
//		}
		case tri_hp_cns_vtype::inflow: {
			temp = new inflow_pt(*this,*vbdry(bnum));
			break;
		}
//		case tri_hp_cns_vtype::hybrid_slave_point: {
//			temp = new hybrid_slave_pt(*this,*vbdry(bnum));
//			break;
//		}
//		case tri_hp_cns_vtype::hybrid_point: {
//			temp = new hybrid_pt(*this,*vbdry(bnum));
//			break;
//		}
		case tri_hp_cns_vtype::unknown: {
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
class tri_hp_cns_etype {
	public:
		static const int ntypes = 10;
		enum ids {unknown=-1,inflow,outflow,characteristic,euler,
			symmetry,applied_stress,force_coupling,adiabatic,shock,outflow_supersonic};
		static const char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i)
				if (!strcmp(nin,names[i])) return(i);
			return(-1);
		}
};

const char tri_hp_cns_etype::names[ntypes][40] = {"inflow","outflow","characteristic","euler",
    "symmetry","applied_stress","force_coupling","adiabatic","shock","outflow_supersonic"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_edge_bdry* tri_hp_cns::getnewedgeobject(int bnum, std::string name) {
	std::string keyword,val;
	std::istringstream data;
	int type;          
	hp_edge_bdry *temp;  

	type = tri_hp_cns_etype::getid(name.c_str());
	switch(type) {
		case tri_hp_cns_etype::inflow: {
			temp = new inflow(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_cns_etype::outflow: {
			temp = new generic(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_cns_etype::characteristic: {
			temp = new characteristic(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_cns_etype::euler: {
			temp = new euler(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_cns_etype::symmetry: {
			temp = new symmetry(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_cns_etype::applied_stress: {
			temp = new applied_stress(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_cns_etype::shock: {
			if (dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))) {
				temp = new shock(*this,*ebdry(bnum));
				dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))->physics = temp;
			}
			else {
				std::cerr << "use coupled physics for shock boundary" << std::endl;
				sim::abort(__LINE__,__FILE__,&std::cerr);
				assert(0);
			}
			break;
		}
//		case tri_hp_cns_etype::force_coupling: {
//			temp = new force_coupling(*this,*ebdry(bnum));
//			break;
//		}
		case tri_hp_cns_etype::adiabatic: {
			temp = new adiabatic(*this,*ebdry(bnum));
			break;
		}
        case tri_hp_cns_etype::outflow_supersonic: {
            temp = new outflow_supersonic(*this,*ebdry(bnum));
            break;
        }
		default: {
			return(tri_hp::getnewedgeobject(bnum,name));
		}
	}    
	gbl->ebdry_gbls(bnum) = temp->create_global_structure();

	return(temp);
}


