#include "tri_hp_cd.h"
#include "bdry_cd.h"
#include <input_map.h>

using namespace bdry_cd;

/** \brief Helper object for edge_bdry 
 *
 * \ingroup boundary
 * Contains list of all side_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tri_hp_cd_stype {
    public:
		static const int ntypes = 5;
		enum ids {unknown=-1,plain,dirichlet,adiabatic,characteristic,mixed_cd};
		static const char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i)
				if (!strcmp(nin,names[i])) return(i);
			return(unknown);
		}
};

const char tri_hp_cd_stype::names[ntypes][40] = {"plain","dirichlet","adiabatic","characteristic","mixed"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_edge_bdry* tri_hp_cd::getnewsideobject(int bnum, input_map &bdrydata) {
	std::string keyword,val;
	std::istringstream data;
	int type;          
	hp_edge_bdry *temp;  


	keyword =  ebdry(bnum)->idprefix + "_cd_type";
	if (bdrydata.get(keyword,val)) {
		type = tri_hp_cd_stype::getid(val.c_str());
		if (type == tri_hp_cd_stype::unknown)  {
			*gbl->log << "unknown side type:" << val << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
	}
	else {
		type = tri_hp_cd_stype::unknown;
	}


	switch(type) {
		case tri_hp_cd_stype::plain: {
			temp = new hp_edge_bdry(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_cd_stype::dirichlet: {
			temp = new dirichlet(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_cd_stype::adiabatic: {
			temp = new neumann(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_cd_stype::characteristic: {
			temp = new characteristic(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_cd_stype::mixed_cd: {
			temp = new mixed(*this,*ebdry(bnum));
			break;
		}
		default: {
			temp = tri_hp::getnewsideobject(bnum,bdrydata);
			break;
		}
	}

	return(temp);
}
