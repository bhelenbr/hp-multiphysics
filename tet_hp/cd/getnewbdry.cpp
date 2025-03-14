#include "tet_hp_cd.h"
#include "bdry_cd.h"
#include <input_map.h>

using namespace bdry_cd;


class tet_hp_cd_vtype {
public:
	static const int ntypes = 3;
	enum ids {unknown=-1,plain,dirichlet,adiabatic};
	const static char names[ntypes][40];
	static int getid(const char *nin) {
		for(int i=0;i<ntypes;++i) 
			if (!strcmp(nin,names[i])) return(i);
		return(unknown);
	}
};

const char tet_hp_cd_vtype::names[ntypes][40] = {"plain","dirichlet","adiabatic"};
	
hp_vrtx_bdry* tet_hp_cd::getnewvrtxobject(int bnum, std::string name) {
    std::string keyword,val;
    std::istringstream data;
    int type;          
    hp_vrtx_bdry *temp;  
    

		type = tet_hp_cd_vtype::getid(name.c_str());
    switch(type) {
		case tet_hp_cd_vtype::plain: {
			temp = new hp_vrtx_bdry(*this,*vbdry(bnum));
			break;
		}
		case tet_hp_cd_vtype::dirichlet: {
			temp = new dirichlet_pt(*this,*vbdry(bnum));
			break;
		}
		case tet_hp_cd_vtype::adiabatic: {
			temp = new neumann_pt(*this,*vbdry(bnum));
			break;
		}
        default: {
            temp = tet_hp::getnewvrtxobject(bnum,name);
            break;
        }
    } 
    return(temp);
}


/** \brief Helper object for edge_bdry 
 *
 * \ingroup boundary
 * Contains list of all side_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tet_hp_cd_etype {
public:
	static const int ntypes = 3;
	enum ids {unknown=-1,plain,dirichlet,adiabatic};
	static const char names[ntypes][40];
	static int getid(const char *nin) {
		for(int i=0;i<ntypes;++i)
			if (!strcmp(nin,names[i])) return(i);
		return(-1);
	}
};

const char tet_hp_cd_etype::names[ntypes][40] = {"plain","dirichlet","adiabatic"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_edge_bdry* tet_hp_cd::getnewedgeobject(int bnum, std::string name) {
    std::string keyword,val;
    std::istringstream data;
    int type;          
    hp_edge_bdry *temp;  
    
	
		type = tet_hp_cd_etype::getid(name.c_str());
    switch(type) {
		case tet_hp_cd_etype::plain: {
			temp = new hp_edge_bdry(*this,*ebdry(bnum));
			break;
		}
		case tet_hp_cd_etype::dirichlet: {
			temp = new dirichlet_edge(*this,*ebdry(bnum));
			break;
		}
		case tet_hp_cd_etype::adiabatic: {
			temp = new neumann_edge(*this,*ebdry(bnum));
			break;
		}
        default: {
            temp = tet_hp::getnewedgeobject(bnum,name);
            break;
        }
    }    
	
    return(temp);
}


class tet_hp_cd_ftype {
	public:
		static const int ntypes = 3;
		enum ids {unknown=-1,plain,dirichlet,adiabatic};
		static const char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i)
				if (!strcmp(nin,names[i])) return(i);
			return(unknown);
		}
};

const char tet_hp_cd_ftype::names[ntypes][40] = {"plain","dirichlet","adiabatic"};

hp_face_bdry* tet_hp_cd::getnewfaceobject(int bnum, std::string name) {
	std::string keyword,val;
	std::istringstream data;
	int type;          
	hp_face_bdry *temp;  
	

	type = tet_hp_cd_ftype::getid(name.c_str());
	switch(type) {
		case tet_hp_cd_ftype::plain: {
			temp = new hp_face_bdry(*this,*fbdry(bnum));
			break;
		}
		case tet_hp_cd_ftype::dirichlet: {
			temp = new dirichlet(*this,*fbdry(bnum));
			break;
		}
		case tet_hp_cd_ftype::adiabatic: {
			temp = new neumann(*this,*fbdry(bnum));
			break;
		}

		default: {
			temp = tet_hp::getnewfaceobject(bnum,name);
			break;
		}
	}
	
	return(temp);
}
