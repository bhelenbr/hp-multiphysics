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

hp_vrtx_bdry* tri_hp_cns_explicit::getnewvrtxobject(int bnum, std::string name) {
	std::string keyword,val;
	std::istringstream data;
	int type;          
	hp_vrtx_bdry *temp;  

	type = tri_hp_cns_explicity_vtype::getid(name.c_str());
	switch(type) {

		case tri_hp_cns_explicit_vtype::inflow: {
			temp = new inflow_pt(*this,*vbdry(bnum));
			break;
		}
		case tri_hp_cns_explicit_vtype::unknown: {
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
class tri_hp_cns_explicit_stype {
	public:
		static const int ntypes = 5;
		enum ids {unknown=-1,inflow,outflow,characteristic,applied_stress,adiabatic};
		static const char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i)
				if (!strcmp(nin,names[i])) return(i);
			return(-1);
		}
};

const char tri_hp_cns_explicit_stype::names[ntypes][40] = {"inflow","outflow","characteristic","applied_stress","adiabatic"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_edge_bdry* tri_hp_cns_explicit::getnewsideobject(int bnum, std::string name) {
	std::string keyword,val;
	std::istringstream data;
	int type;          
	hp_edge_bdry *temp;  


	type = tri_hp_cns_explicit_stype::getid(name.c_str());
	switch(type) {
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
			return(tri_hp::getnewsideobject(bnum,name));
		}
	}    
	gbl->ebdry_gbls(bnum) = temp->create_global_structure();

	return(temp);
}


