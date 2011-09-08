/*
 *  getnew_bdry.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 1/12/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
#include "tet_hp_cns.h"
#include "bdry_cns.h"
#include <input_map.h>
#include <tet_boundary.h>

using namespace bdry_cns;

/** \brief Helper object for vrtx_bdry 
 *
 * \ingroup boundary
 * Contains list of all vrtx_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tet_hp_cns_vtype {
    public:
        static const int ntypes = 3;
        enum ids {unknown=-1,inflow,adiabatic,outflow};
        const static char names[ntypes][40];
        static int getid(const char *nin) {
            for(int i=0;i<ntypes;++i) 
                if (!strcmp(nin,names[i])) return(i);
            return(unknown);
        }
};

const char tet_hp_cns_vtype::names[ntypes][40] = {"inflow","adiabatic","outflow"};

hp_vrtx_bdry* tet_hp_cns::getnewvrtxobject(int bnum, input_map &bdrydata) {
    std::string keyword,val;
    std::istringstream data;
    int type;          
    hp_vrtx_bdry *temp;  
    
    keyword = vbdry(bnum)->idprefix + "_cns_type";
    if (bdrydata.get(keyword,val)) {
        type = tet_hp_cns_vtype::getid(val.c_str());
        if (type == tet_hp_cns_vtype::unknown)  {
            *gbl->log << "unknown vertex type:" << val << std::endl;
            exit(1);
        }
    }
    else {
        type = tet_hp_cns_vtype::unknown;
    }
    
    
    switch(type) {
		case tet_hp_cns_vtype::inflow: {
			temp = new inflow_pt(*this,*vbdry(bnum));
			break;
		}
		case tet_hp_cns_vtype::adiabatic: {
			temp = new adiabatic_pt(*this,*vbdry(bnum));
			break;
		}
		case tet_hp_cns_vtype::outflow: {
			temp = new outflow_pt(*this,*vbdry(bnum));
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
class tet_hp_cns_etype {
    public:
        static const int ntypes = 4;
        enum ids {unknown=-1,plain,inflow,outflow,adiabatic};
        static const char names[ntypes][40];
        static int getid(const char *nin) {
            for(int i=0;i<ntypes;++i)
                if (!strcmp(nin,names[i])) return(i);
            return(-1);
        }
};

const char tet_hp_cns_etype::names[ntypes][40] = {"plain","inflow","outflow","adiabatic"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_edge_bdry* tet_hp_cns::getnewedgeobject(int bnum, input_map &bdrydata) {
    std::string keyword,val;
    std::istringstream data;
    int type;          
    hp_edge_bdry *temp;  
    

    keyword =  ebdry(bnum)->idprefix + "_cns_type";
    if (bdrydata.get(keyword,val)) {
        type = tet_hp_cns_etype::getid(val.c_str());
        if (type == tet_hp_cns_etype::unknown)  {
            *gbl->log << "unknown edge type:" << val << std::endl;
            exit(1);
        }
    }
    else {
        type = tet_hp_cns_etype::unknown;
    }

    switch(type) {
		case tet_hp_cns_etype::plain: {
			temp = new generic_edge(*this,*ebdry(bnum));
			break;
		}
		case tet_hp_cns_etype::inflow: {
			temp = new inflow_edge(*this,*ebdry(bnum));
			break;
		}
		case tet_hp_cns_etype::outflow: {
			temp = new neumann_edge(*this,*ebdry(bnum));
			break;
		}
		case tet_hp_cns_etype::adiabatic: {
			temp = new adiabatic_edge(*this,*ebdry(bnum));
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


/** \brief Helper object for face_bdry 
 *
 * \ingroup boundary
 * Contains list of all face_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tet_hp_cns_ftype {
	public:
		static const int ntypes = 7;
		enum ids {unknown=-1,plain,inflow,outflow,adiabatic,characteristic,applied_stress,pure_slip};
		static const char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i)
				if (!strcmp(nin,names[i])) return(i);
			return(-1);
		}
};

const char tet_hp_cns_ftype::names[ntypes][40] = {"plain","inflow","outflow","adiabatic","characteristic","applied_stress","pure_slip"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_face_bdry* tet_hp_cns::getnewfaceobject(int bnum, input_map &bdrydata) {
	std::string keyword,val;
	std::istringstream data;
	int type;          
	hp_face_bdry *temp;  
	

	keyword =  fbdry(bnum)->idprefix + "_cns_type";
	if (bdrydata.get(keyword,val)) {
		type = tet_hp_cns_ftype::getid(val.c_str());
		if (type == tet_hp_cns_ftype::unknown)  {
			*gbl->log << "unknown face type:" << val << std::endl;
			exit(1);
		}
	}
	else {
		type = tet_hp_cns_ftype::unknown;
	}

	switch(type) {
		case tet_hp_cns_ftype::plain: {
			temp = new generic(*this,*fbdry(bnum));
			break;
		}
		case tet_hp_cns_ftype::inflow: {
			temp = new inflow(*this,*fbdry(bnum));
			break;
		}
		case tet_hp_cns_ftype::outflow: {
			temp = new neumann(*this,*fbdry(bnum));
			break;
		}
		case tet_hp_cns_ftype::adiabatic: {
			temp = new adiabatic(*this,*fbdry(bnum));
			break;
		}
		case tet_hp_cns_ftype::characteristic: {
			temp = new characteristic(*this,*fbdry(bnum));
			break;
		}
		case tet_hp_cns_ftype::applied_stress: {
			temp = new applied_stress(*this,*fbdry(bnum));
			break;
		}
		case tet_hp_cns_ftype::pure_slip: {
			temp = new pure_slip(*this,*fbdry(bnum));
			break;
		}
		default: {
			temp = tet_hp::getnewfaceobject(bnum,bdrydata);
			break;
		}
	}	
	gbl->fbdry_gbls(bnum) = temp->create_global_structure();

	return(temp);
}




