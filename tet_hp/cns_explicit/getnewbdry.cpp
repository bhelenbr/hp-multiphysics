/*
 *  getnew_bdry.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 1/12/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */
#include "tet_hp_cns_explicit.h"
#include "bdry_cns_explicit.h"
#include <input_map.h>
#include <tet_boundary.h>

using namespace bdry_cns_explicit;

/** \brief Helper object for vrtx_bdry 
 *
 * \ingroup boundary
 * Contains list of all vrtx_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tet_hp_cns_explicit_vtype {
    public:
        static const int ntypes = 5;
        enum ids {unknown=-1,surface_inflow,surface_periodic,surface_outflow,surface_outflow_planar,inflow};
        const static char names[ntypes][40];
        static int getid(const char *nin) {
            for(int i=0;i<ntypes;++i) 
                if (!strcmp(nin,names[i])) return(i);
            return(unknown);
        }
};

const char tet_hp_cns_explicit_vtype::names[ntypes][40] = {"surface_inflow","surface_periodic","surface_outflow","surface_outflow_planar","inflow"};

hp_vrtx_bdry* tet_hp_cns_explicit::getnewvrtxobject(int bnum, input_map &bdrydata) {
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
			type = tet_hp_cns_explicit_vtype::getid(val.c_str());
		}
    
    
    switch(type) {
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
class tet_hp_cns_explicit_etype {
    public:
        static const int ntypes = 11;
        enum ids {unknown=-1,plain,inflow,outflow,characteristic,euler,
            symmetry,applied_stress,surface,surface_slave,hybrid_surface_levelset,force_coupling};
        static const char names[ntypes][40];
        static int getid(const char *nin) {
            for(int i=0;i<ntypes;++i)
                if (!strcmp(nin,names[i])) return(i);
            return(-1);
        }
};

const char tet_hp_cns_explicit_etype::names[ntypes][40] = {"plain","inflow","outflow","characteristic","euler",
    "symmetry","applied_stress","surface","surface_slave","hybrid_surface_levelset","force_coupling"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_edge_bdry* tet_hp_cns_explicit::getnewedgeobject(int bnum, input_map &bdrydata) {
    std::string keyword,val;
    std::istringstream data;
    int type;          
    hp_edge_bdry *temp;  
    

		keyword =  ebdry(bnum)->idprefix + "_hp_type";
		if (!bdrydata.get(keyword,val)) {
			*gbl->log << "missing edge type:" << keyword << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
		else {
			type = tet_hp_cns_explicit_etype::getid(val.c_str());
		}

    switch(type) {
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
class tet_hp_cns_explicit_ftype {
	public:
		static const int ntypes = 7;
		enum ids {unknown=-1,plain,inflow,outflow,adiabatic,characteristic,applied_stress,specified_flux};
		static const char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i)
				if (!strcmp(nin,names[i])) return(i);
			return(-1);
		}
};

const char tet_hp_cns_explicit_ftype::names[ntypes][40] = {"plain","inflow","outflow","adiabatic","characteristic","applied_stress","specified_flux"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_face_bdry* tet_hp_cns_explicit::getnewfaceobject(int bnum, input_map &bdrydata) {
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
		type = tet_hp_cns_explicit_ftype::getid(val.c_str());
	}

	switch(type) {
		case tet_hp_cns_explicit_ftype::plain: {
			temp = new generic(*this,*fbdry(bnum));
			break;
		}
		case tet_hp_cns_explicit_ftype::inflow: {
			temp = new inflow(*this,*fbdry(bnum));
			break;
		}
		case tet_hp_cns_explicit_ftype::outflow: {
			temp = new neumann(*this,*fbdry(bnum));
			break;
		}
		case tet_hp_cns_explicit_ftype::adiabatic: {
			temp = new adiabatic(*this,*fbdry(bnum));
			break;
		}
		case tet_hp_cns_explicit_ftype::characteristic: {
			temp = new characteristic(*this,*fbdry(bnum));
			break;
		}
		case tet_hp_cns_explicit_ftype::applied_stress: {
			temp = new applied_stress(*this,*fbdry(bnum));
			break;
		}
		case tet_hp_cns_explicit_ftype::specified_flux: {
			temp = new specified_flux(*this,*fbdry(bnum));
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




