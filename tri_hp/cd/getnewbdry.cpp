#include "tri_hp_cd.h"
#include "bdry_cd.h"
#include <input_map.h>
#include <tri_boundary.h>
#include "melt_cd.h"

using namespace bdry_cd;

class tri_hp_cd_vtype {
public:
	static const int ntypes = 1;
	enum ids {unknown=-1,melt_facet_pt,triple_junction};
	const static char names[ntypes][40];
	static int getid(const char *nin) {
		for(int i=0;i<ntypes;++i)
			if (!strcmp(nin,names[i])) return(i);
		return(unknown);
	}
};

const char tri_hp_cd_vtype::names[ntypes][40] = {"melt_facet_pt"};

hp_vrtx_bdry* tri_hp_cd::getnewvrtxobject(int bnum, std::string name) {
	std::string keyword,val;
	std::istringstream data;
	int type;
	hp_vrtx_bdry *temp;
	
	type = tri_hp_cd_vtype::getid(name.c_str());
	switch(type) {
		case tri_hp_cd_vtype::melt_facet_pt: {
			temp = new bdry_cd::melt_facet_pt(*this,*vbdry(bnum));
			break;
		}
		case tri_hp_cd_vtype::unknown: {
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
class tri_hp_cd_etype {
    public:
		static const int ntypes = 4;
		enum ids {unknown=-1,dirichlet,adiabatic,characteristic,melt};
		static const char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i)
				if (!strcmp(nin,names[i])) return(i);
			return(unknown);
		}
};

const char tri_hp_cd_etype::names[ntypes][40] = {"dirichlet","adiabatic","characteristic","melt"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_edge_bdry* tri_hp_cd::getnewedgeobject(int bnum, std::string name) {
	std::string keyword,val;
	std::istringstream data;
	int type;          
	hp_edge_bdry *temp;
	
	type = tri_hp_cd_etype::getid(name.c_str());
	switch(type) {
		case tri_hp_cd_etype::dirichlet: {
			temp = new dirichlet(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_cd_etype::adiabatic: {
			temp = new generic(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_cd_etype::characteristic: {
			temp = new characteristic(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_cd_etype::melt:  {
			if (dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))) {
				temp = new  melt_cd(*this,*ebdry(bnum));
				dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))->physics = temp;
			}
			else {
				std::cerr << "use coupled physics for melt boundary" << std::endl;
				sim::abort(__LINE__,__FILE__,&std::cerr);
				assert(0);
			}
			break;
		}
		default: {
			return(tri_hp::getnewedgeobject(bnum,name));
		}
	}
	gbl->ebdry_gbls(bnum) = temp->create_global_structure();
	
	return(temp);
}

