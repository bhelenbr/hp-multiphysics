/*
 *  getnewins_bdry.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 1/12/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
#include "tri_hp_cns_explicit.h"
#include "bdry_cns_explicit.h"
#include <input_map.h>
#include <tri_boundary.h>

using namespace bdry_cns_explicit;

/** \brief Helper object for vrtx_bdry 
 *
 * \ingroup boundary
 * Contains list of all vrtx_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tri_hp_cns_explicit_vtype {
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

const char tri_hp_cns_explicit_vtype::names[ntypes][40] = {"surface_inflow","surface_periodic","surface_outflow","surface_outflow_planar",
	"inflow","hybrid_slave_point","hybrid_point"};

hp_vrtx_bdry* tri_hp_cns_explicit::getnewvrtxobject(int bnum, input_map &bdrydata) {
	std::string keyword,val;
	std::istringstream data;
	int type;          
	hp_vrtx_bdry *temp;  

	keyword = vbdry(bnum)->idprefix + "_cns_explicit_type";
	if (bdrydata.get(keyword,val)) {
		type = tri_hp_cns_explicit_vtype::getid(val.c_str());
		if (type == tri_hp_cns_explicit_vtype::unknown)  {
			*gbl->log << "unknown vertex type:" << val << std::endl;
			exit(1);
		}
	}
	else {
		type = tri_hp_cns_explicit_vtype::unknown;
	}


	switch(type) {
//		case tri_hp_cns_explicit_vtype::surface_inflow: {
//			temp = new surface_fixed_pt(*this,*vbdry(bnum));
//			break;
//		}
//		case tri_hp_cns_explicit_vtype::surface_periodic: {
//			temp = new surface_periodic_pt(*this,*vbdry(bnum));
//			break;
//		}
//		case tri_hp_cns_explicit_vtype::surface_outflow: {
//			temp = new surface_outflow_endpt(*this,*vbdry(bnum));
//			break;
//		}
//		case tri_hp_cns_explicit_vtype::surface_outflow_planar: {
//			temp = new surface_outflow_planar(*this,*vbdry(bnum));
//			break;
//		}
		case tri_hp_cns_explicit_vtype::inflow: {
			temp = new inflow_pt(*this,*vbdry(bnum));
			break;
		}
//		case tri_hp_cns_explicit_vtype::hybrid_slave_point: {
//			temp = new hybrid_slave_pt(*this,*vbdry(bnum));
//			break;
//		}
//		case tri_hp_cns_explicit_vtype::hybrid_point: {
//			temp = new hybrid_pt(*this,*vbdry(bnum));
//			break;
//		}
		default: {
			temp = tri_hp::getnewvrtxobject(bnum,bdrydata);
			break;
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
class tri_hp_cns_explicit_stype {
	public:
		static const int ntypes = 10;
		enum ids {unknown=-1,plain,inflow,outflow,characteristic,euler,
			symmetry,applied_stress,surface,surface_slave,force_coupling};
		static const char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i)
				if (!strcmp(nin,names[i])) return(i);
			return(-1);
		}
};

const char tri_hp_cns_explicit_stype::names[ntypes][40] = {"plain","inflow","outflow","characteristic","euler",
    "symmetry","applied_stress","surface","surface_slave","force_coupling"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_edge_bdry* tri_hp_cns_explicit::getnewsideobject(int bnum, input_map &bdrydata) {
	std::string keyword,val;
	std::istringstream data;
	int type;          
	hp_edge_bdry *temp;  


	keyword =  ebdry(bnum)->idprefix + "_cns_explicit_type";
	if (bdrydata.get(keyword,val)) {
		type = tri_hp_cns_explicit_stype::getid(val.c_str());
		if (type == tri_hp_cns_explicit_stype::unknown)  {
			*gbl->log << "unknown side type:" << val << std::endl;
			exit(1);
		}
	}
	else {
		type = tri_hp_cns_explicit_stype::unknown;
	}

	switch(type) {
		case tri_hp_cns_explicit_stype::plain: {
			temp = new generic(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_cns_explicit_stype::inflow: {
			temp = new inflow(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_cns_explicit_stype::outflow: {
			temp = new neumann(*this,*ebdry(bnum));
			break;
		}
//		case tri_hp_cns_explicit_stype::characteristic: {
//			temp = new characteristic(*this,*ebdry(bnum));
//			break;
//		}
//		case tri_hp_cns_explicit_stype::euler: {
//			temp = new euler(*this,*ebdry(bnum));
//			break;
//		}
//		case tri_hp_cns_explicit_stype::symmetry: {
//			temp = new symmetry(*this,*ebdry(bnum));
//			break;
//		}
		case tri_hp_cns_explicit_stype::applied_stress: {
			temp = new applied_stress(*this,*ebdry(bnum));
			break;
		}
//		case tri_hp_cns_explicit_stype::surface: {
//			if (dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))) {
//				temp = new surface(*this,*ebdry(bnum));
//				dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))->physics = temp;
//			}
//			else {
//				std::cerr << "use coupled physics for surface boundary" << std::endl;
//				exit(1);
//			}
//			break;
//		}
//		case tri_hp_cns_explicit_stype::surface_slave: {
//			temp = new surface_slave(*this,*ebdry(bnum));
//			dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))->physics = temp;
//			break;
//		}
//		case tri_hp_cns_explicit_stype::force_coupling: {
//			temp = new force_coupling(*this,*ebdry(bnum));
//			break;
//		}
		default: {
			temp = tri_hp::getnewsideobject(bnum,bdrydata);
			break;
		}
	}    
	gbl->ebdry_gbls(bnum) = temp->create_global_structure();

	return(temp);
}

