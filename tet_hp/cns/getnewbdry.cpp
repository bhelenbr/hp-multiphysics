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
        static const int ntypes = 5;
        enum ids {unknown=-1,surface_inflow,surface_periodic,surface_outflow,surface_outflow_planar,inflow};
        const static char names[ntypes][40];
        static int getid(const char *nin) {
            for(int i=0;i<ntypes;++i) 
                if (!strcmp(nin,names[i])) return(i);
            return(unknown);
        }
};

const char tet_hp_cns_vtype::names[ntypes][40] = {"surface_inflow","surface_periodic","surface_outflow","surface_outflow_planar","inflow"};

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
//        case tet_hp_cns_vtype::surface_inflow: {
//            temp = new surface_fixed_pt(*this,*vbdry(bnum));
//            break;
//        }
//        case tet_hp_cns_vtype::surface_periodic: {
//            temp = new surface_periodic_pt(*this,*vbdry(bnum));
//            break;
//        }
//        case tet_hp_cns_vtype::surface_outflow: {
//            temp = new surface_outflow_endpt(*this,*vbdry(bnum));
//            break;
//        }
//        case tet_hp_cns_vtype::surface_outflow_planar: {
//            temp = new surface_outflow_planar(*this,*vbdry(bnum));
//            break;
//        }
//        case tet_hp_cns_vtype::inflow: {
//            temp = new inflow_pt(*this,*vbdry(bnum));
//            break;
//        }
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

const char tet_hp_cns_etype::names[ntypes][40] = {"plain","inflow","outflow","characteristic","euler",
    "symmetry","applied_stress","surface","surface_slave","hybrid_surface_levelset","force_coupling"};

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
//        case tet_hp_cns_etype::plain: {
//            temp = new generic(*this,*ebdry(bnum));
//            break;
//        }
//        case tet_hp_cns_etype::inflow: {
//            temp = new inflow(*this,*ebdry(bnum));
//            break;
//        }
//        case tet_hp_cns_etype::outflow: {
//            temp = new neumann(*this,*ebdry(bnum));
//            break;
//        }
//        case tet_hp_cns_etype::characteristic: {
//            temp = new characteristic(*this,*ebdry(bnum));
//            break;
//        }
//        case tet_hp_cns_etype::euler: {
//            temp = new euler(*this,*ebdry(bnum));
//            break;
//        }
//        case tet_hp_cns_etype::symmetry: {
//            temp = new symmetry(*this,*ebdry(bnum));
//            break;
//        }
//        case tet_hp_cns_etype::applied_stress: {
//            temp = new applied_stress(*this,*ebdry(bnum));
//            break;
//        }
//        case tet_hp_cns_etype::surface: {
//            if (dynamic_cast<ecoupled_physics_interface *>(ebdry(bnum))) {
//                temp = new surface(*this,*ebdry(bnum));
//                dynamic_cast<ecoupled_physics_interface *>(ebdry(bnum))->physics = temp;
//            }
//            else {
//                std::cerr << "use coupled physics for surface boundary\n";
//                exit(1);
//            }
//            break;
//        }
//        case tet_hp_cns_etype::surface_slave: {
//            temp = new surface_slave(*this,*ebdry(bnum));
//            dynamic_cast<ecoupled_physics_interface *>(ebdry(bnum))->physics = temp;
//            break;
//        }
//        case tet_hp_cns_etype::hybrid_surface_levelset: {
//            temp = new hybrid_surface_levelset(*this,*ebdry(bnum));
//            break;
//        }
//		case tet_hp_cns_etype::force_coupling: {
//            temp = new force_coupling(*this,*ebdry(bnum));
//            break;
//        }
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
		static const int ntypes = 4;
		enum ids {unknown=-1,plain,inflow,outflow,applied_stress};
		static const char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i)
				if (!strcmp(nin,names[i])) return(i);
			return(-1);
		}
};

const char tet_hp_cns_ftype::names[ntypes][40] = {"plain","inflow","outflow","applied_stress"};

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
		case tet_hp_cns_ftype::applied_stress: {
			temp = new applied_stress(*this,*fbdry(bnum));
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



