//
//  hp_getnewbdry.cpp
//  tri_hp
//
//  Created by Brian Helenbrook on 3/8/15.
//
//

#include <stdio.h>
#include "hp_boundary.h"
#include "hp_coupled_boundary.h"
#include <tri_boundary.h>

class tri_hp_vtype {
public:
	static const int ntypes = 5;
	enum ids {unknown=-1,plain,multi_physics_pnt,hp_deformable_fixed_pnt,hp_deformable_free_pnt,hp_deformable_follower_pnt};
	const static char names[ntypes][40];
	static int getid(const char *nin) {
		for(int i=0;i<ntypes;++i)
			if (!strcmp(nin,names[i])) return(i);
		return(unknown);
	}
};

const char tri_hp_vtype::names[ntypes][40] = {"plain","multi_physics_pnt","hp_deformable_fixed_pnt","hp_deformable_free_pnt","hp_deformable_follower_pnt"};

hp_vrtx_bdry* tri_hp::getnewvrtxobject(int bnum, std::string name) {
	std::string keyword,val;
	std::istringstream data;
	int type;
	hp_vrtx_bdry *temp;
	
	type = tri_hp_vtype::getid(name.c_str());
	if (type == tri_hp_vtype::unknown) {
		*gbl->log << "unknown vrtx type:" << keyword << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
	switch(type) {
		case tri_hp_vtype::multi_physics_pnt: {
			temp = new multi_physics_pnt(*this,*vbdry(bnum));
			break;
		}
		case tri_hp_vtype::hp_deformable_fixed_pnt: {
			temp = new hp_deformable_fixed_pnt(*this,*vbdry(bnum));
			break;
		}
		case tri_hp_vtype::hp_deformable_free_pnt: {
			temp = new hp_deformable_free_pnt(*this,*vbdry(bnum));
			break;
		}
		case tri_hp_vtype::hp_deformable_follower_pnt: {
			temp = new hp_deformable_follower_pnt(*this,*vbdry(bnum));
			break;
		}
		default: {
			temp = new hp_vrtx_bdry(*this,*vbdry(bnum));
		}
	}
	gbl->vbdry_gbls(bnum) = temp->create_global_structure();
	return(temp);
}

class tri_hp_etype {
public:
	static const int ntypes = 5;
	enum ids {unknown=-1,plain,symbolic_with_integration_by_parts,translating_surface,partition};
	static const char names[ntypes][40];
	static int getid(const char *nin) {
		for(int i=0;i<ntypes;++i)
			if (!strcmp(nin,names[i])) return(i);
		return(-1);
	}
};

const char tri_hp_etype::names[ntypes][40] = {"plain","symbolic_ibp","translating_surface","partition"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_edge_bdry* tri_hp::getnewedgeobject(int bnum, std::string name) {
	std::string keyword,val;
	std::istringstream data;
	int type;
	hp_edge_bdry *temp;
	
	
	type = tri_hp_etype::getid(name.c_str());
	
	if (type == tri_hp_etype::plain && ebdry(bnum)->in_group(boundary::partitions)) {
		type = tri_hp_etype::partition;
	}
	if (type == tri_hp_etype::unknown) {
		*gbl->log << "unknown side type:" << keyword << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
	
	
	switch(type) {
		case tri_hp_etype::symbolic_with_integration_by_parts: {
			temp = new symbolic_with_integration_by_parts(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_etype::translating_surface:  {
			if (dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))) {
				temp = new  translating_surface(*this,*ebdry(bnum));
				dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))->physics = temp;
			}
			else {
				std::cerr << "use coupled physics for translating boundary" << std::endl;
				sim::abort(__LINE__,__FILE__,&std::cerr);
				assert(0);
			}
			break;
		}
		case tri_hp_etype::partition: {
			temp = new hp_partition(*this,*ebdry(bnum));
			break;
		}
		default: {
			temp = new hp_edge_bdry(*this,*ebdry(bnum));
		}
	}
	gbl->ebdry_gbls(bnum) = temp->create_global_structure();
	
	return(temp);
}
