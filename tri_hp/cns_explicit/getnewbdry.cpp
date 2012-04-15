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
		static const int ntypes = 1;
		enum ids {unknown=-1,inflow};
		const static char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i) 
				if (!strcmp(nin,names[i])) return(i);
			return(unknown);
		}
};

const char tri_hp_cns_explicit_vtype::names[ntypes][40] = {"inflow"};

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
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
	}
	else {
		type = tri_hp_cns_explicit_vtype::unknown;
	}


	switch(type) {

		case tri_hp_cns_explicit_vtype::inflow: {
			temp = new inflow_pt(*this,*vbdry(bnum));
			break;
		}
		default: {
			return(tri_hp::getnewvrtxobject(bnum,bdrydata));
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
		static const int ntypes = 6;
		enum ids {unknown=-1,plain,inflow,outflow,characteristic,applied_stress,adiabatic};
		static const char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i)
				if (!strcmp(nin,names[i])) return(i);
			return(-1);
		}
};

const char tri_hp_cns_explicit_stype::names[ntypes][40] = {"plain","inflow","outflow","characteristic","applied_stress","adiabatic"};

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
			sim::abort(__LINE__,__FILE__,gbl->log);
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
			temp = new generic(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_cns_explicit_stype::characteristic: {
			temp = new characteristic(*this,*ebdry(bnum));
			break;
		}

		case tri_hp_cns_explicit_stype::applied_stress: {
			temp = new applied_stress(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_cns_explicit_stype::adiabatic: {
			temp = new adiabatic(*this,*ebdry(bnum));
			break;
		}
		default: {
			return(tri_hp::getnewsideobject(bnum,bdrydata));
		}
	}    
	gbl->ebdry_gbls(bnum) = temp->create_global_structure();

	return(temp);
}


