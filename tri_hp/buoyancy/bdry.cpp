/*
 *  bdry.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 12/13/10.
 *  Copyright 2010 Clarkson University. All rights reserved.
 *
 */

#include "bdry_buoyancy.h"
#include <tri_boundary.h>

using namespace bdry_buoyancy;
#include "bdry_buoyancy.h"
#include <myblas.h>

//#define MPDEBUG

//#define DEBUG

using namespace bdry_buoyancy;

#ifdef SKIP

void solid_fluid::init(input_map& inmap,void* gbl_in) {
	std::string keyword,val;
	std::istringstream data;
	std::string filename;
	
	/* Load in the heat flux for radiation */
	keyword = base.idprefix + "_hp_typelist";
	inmap[keyword] = "0 0 1 0";
	symbolic::init(inmap,gbl_in);
	
	/* Now that we have loaded flux, make pressure neumann */
	essential_indices.clear();
	for(int n=0;n<x.ND;++n)
		essential_indices.push_back(n);
	type(Range(0,x.ND)) = essential;
	type(Range(x.ND+1,x.NV-1)) = natural;
		
	return;
}
#endif

/** \brief Helper object for vrtx_bdry 
 *
 * \ingroup boundary
 * Contains list of all vrtx_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tri_hp_buoyancy_vtype {
public:
	static const int ntypes = 3;
	enum ids {unknown=-1,melt_end,melt_inflow,melt_facet_pt};
	const static char names[ntypes][40];
	static int getid(const char *nin) {
		for(int i=0;i<ntypes;++i) 
			if (!strcmp(nin,names[i])) return(i);
		return(unknown);
	}
};

const char tri_hp_buoyancy_vtype::names[ntypes][40] = {"melt_end","melt_inflow","melt_facet_pt"};

hp_vrtx_bdry* tri_hp_buoyancy::getnewvrtxobject(int bnum, input_map &bdrydata) {
	std::string keyword,val;
	std::istringstream data;
	int type;          
	hp_vrtx_bdry *temp;  
	
	keyword = vbdry(bnum)->idprefix + "_buoyancy_type";
	if (bdrydata.get(keyword,val)) {
		type = tri_hp_buoyancy_vtype::getid(val.c_str());
		if (type == tri_hp_buoyancy_vtype::unknown)  {
			*gbl->log << "unknown vertex type:" << val << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
	}
	else {
		type = tri_hp_buoyancy_vtype::unknown;
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
		default: {
			return(tri_hp_ins::getnewvrtxobject(bnum,bdrydata));
		}
	} 
	gbl->vbdry_gbls(bnum) = temp->create_global_structure();
	return(temp);
}



class tri_hp_buoyancy_stype {
	public:
		static const int ntypes = 4;
		enum ids {unknown=-1,surface,melt,melt_kinetics,kellerman};
		static const char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i)
				if (!strcmp(nin,names[i])) return(i);
			return(-1);
		}
};

const char tri_hp_buoyancy_stype::names[ntypes][40] = {"surface","melt","melt_kinetics","kellerman"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_edge_bdry* tri_hp_buoyancy::getnewsideobject(int bnum, input_map &bdrydata) {
	std::string keyword,val;
	std::istringstream data;
	int type;          
	hp_edge_bdry *temp;  
	
	if (bdrydata.get(ebdry(bnum)->idprefix + "_buoyancy_type",val)) {
		type = tri_hp_buoyancy_stype::getid(val.c_str());
		if (type == tri_hp_buoyancy_stype::unknown)  {
			*gbl->log << "unknown side type:" << val << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
	}
	else {
		type = tri_hp_buoyancy_stype::unknown;
	}
	
	switch(type) {
		case tri_hp_buoyancy_stype::surface: {
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
		case tri_hp_buoyancy_stype::kellerman: {
			if (dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))) {
				temp = new kellerman(*this,*ebdry(bnum));
				dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))->physics = temp;
			}
			else {
				std::cerr << "use coupled physics for kellerman boundary" << std::endl;
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

