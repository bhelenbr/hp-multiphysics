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
	
hp_vrtx_bdry* tet_hp_cd::getnewvrtxobject(int bnum, input_map &bdrydata) {
    std::string keyword,val;
    std::istringstream data;
    int type;          
    hp_vrtx_bdry *temp;  
    
		keyword =  vbdry(bnum)->idprefix + "_hp_type";
		if (!bdrydata.get(keyword,val)) {
			*gbl->log << "missing vertex type:" << keyword << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
		else {
			type = tet_hp_cd_vtype::getid(val.c_str());
		}
	
	
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
            temp = tet_hp::getnewvrtxobject(bnum,bdrydata);
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
hp_edge_bdry* tet_hp_cd::getnewedgeobject(int bnum, input_map &bdrydata) {
    std::string keyword,val;
    std::istringstream data;
    int type;          
    hp_edge_bdry *temp;  
    
	
		keyword =  ebdry(bnum)->idprefix + "_hp_type";
		if (!bdrydata.get(keyword,val)) {
			*gbl->log << "missing side type:" << keyword << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
		else {
			type = tet_hp_cd_etype::getid(val.c_str());
		}
	
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
            temp = tet_hp::getnewedgeobject(bnum,bdrydata);
            break;
        }
    }    
    gbl->ebdry_gbls(bnum) = temp->create_global_structure();
	
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

hp_face_bdry* tet_hp_cd::getnewfaceobject(int bnum, input_map &bdrydata) {
	std::string keyword,val;
	std::istringstream data;
	int type;          
	hp_face_bdry *temp;  
	

	keyword =  fbdry(bnum)->idprefix + "_hp_type";
	if (!bdrydata.get(keyword,val)) {
		*gbl->log << "missing face type:" << keyword << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
	else {
		type = tet_hp_cd_ftype::getid(val.c_str());
	}

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
			temp = tet_hp::getnewfaceobject(bnum,bdrydata);
			break;
		}
	}
	
	return(temp);
}
