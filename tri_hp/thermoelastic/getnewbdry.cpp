/*
 *  getnewbdry.cpp
 *  tri_hp
 *
 */
#include "tri_hp_thermoelastic.h"
#include "bdry_thermoelastic.h"
#include "../cd/ibc_cd.h"
#include <input_map.h>

/** \brief Helper object for edge_bdry
 *
 * \ingroup boundary
 * Contains list of all side_bdys's by name
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tri_hp_thermoelastic_etype {
	public:
		static const int ntypes = 7;
		enum ids {unknown=-1,dirichlet,neumann,adiabatic,displacement,temperature,characteristic,curve_edges};
		static const char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i)
				if (!strcmp(nin,names[i])) return(i);
			return(unknown);
		}
};

const char tri_hp_thermoelastic_etype::names[ntypes][40] = {"dirichlet","neumann","adiabatic","displacement","temperature","characteristic","curve_edges"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_edge_bdry* tri_hp_thermoelastic::getnewedgeobject(int bnum, std::string name) {
	hp_edge_bdry *temp;

	int type = tri_hp_thermoelastic_etype::getid(name.c_str());
	switch(type) {
		case tri_hp_thermoelastic_etype::dirichlet: {
			temp = new bdry_thermoelastic::dirichlet(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_thermoelastic_etype::neumann: case tri_hp_thermoelastic_etype::adiabatic: {
			temp = new bdry_ps::neumann(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_thermoelastic_etype::displacement: {
			temp = new bdry_ps::dirichlet(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_thermoelastic_etype::temperature: {
			temp = new bdry_thermoelastic::temperature(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_thermoelastic_etype::characteristic: {
			temp = new bdry_thermoelastic::characteristic(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_thermoelastic_etype::curve_edges: {
			temp = new bdry_ps::curve_edges(*this,*ebdry(bnum));
			break;
		}
		default: {
			return(tri_hp::getnewedgeobject(bnum,name));
		}
	}

	return(temp);
}

/* SAME HEAT SOURCE FUNCTIONS AS cd */
init_bdry_cndtn *tri_hp_thermoelastic::getnewibc(std::string name) {
	if (name == "power")
		return(new ibc_cd::power_src);
	if (name == "soi")
		return(new ibc_cd::soi_src(*this));

	return(tri_hp::getnewibc(name));
}
