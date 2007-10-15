/*
 *  mvpttobdry.cpp
 *  mblock
 *
 *  Created by helenbrk on Wed Oct 24 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_mesh.h"
#include "boundaries.h"

/** \brief Helper object for vrtx_bdry 
 *
 * \ingroup boundary
 * Contains list of all vrtx_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class vtype {
    public:
        static const int ntypes = 3;
        enum ids {plain=1,comm,prdc};
        const static char names[ntypes][40];
        static int getid(const char *nin) {
            for(int i=0;i<ntypes;++i) 
                if (!strcmp(nin,names[i])) return(i+1);
            return(-1);
        }
};

const char vtype::names[ntypes][40] = {"plain","comm","prdc"};

vrtx_bdry* tri_mesh::getnewvrtxobject(int idnum, input_map& bdrydata) {
    std::string keyword,val;
    std::istringstream data;
    ostringstream nstr;
    int type;          
    vrtx_bdry *temp;  

    nstr.str("");
    nstr << idnum << std::flush;        
    keyword = gbl->idprefix +"_v" +nstr.str() + "_type";
    if (bdrydata.get(keyword,val)) {
        type = vtype::getid(val.c_str());
        if (type < 0)  {
            *gbl->log << "unknown vertex type:" << val << std::endl;
            exit(1);
        }
    }
    else {
        type = vtype::plain;
    }
        
    switch(type) {
        case vtype::plain: {
            temp = new vrtx_bdry(idnum,*this);
            break;
        }
        case vtype::comm: {
            temp = new vcomm(idnum,*this);
            break;
        }
        case vtype::prdc: {
            temp = new vprdc(idnum,*this);
            break;
        }
        default: {
            std::cout << "unrecognizable pnt type: " <<  type << " idnum: " << idnum << std::endl;
            temp = new vrtx_bdry(idnum,*this);
            break;
        }
    } 
    
    temp->input(bdrydata);
    
    return(temp);
}


/** \brief Helper object for edge_bdry 
 *
 * \ingroup boundary
 * Contains list of all side_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class etype {
    public:
        static const int ntypes = 12;
        enum ids {plain=1, comm, partition, prdc, symbolic, coupled_symbolic, coupled_symbolic_comm, spline,
            circle, naca, ellipse};
        static const char names[ntypes][40];
        static int getid(const char *nin) {
            for(int i=0;i<ntypes;++i)
                if (!strcmp(nin,names[i])) return(i+1);
            return(-1);
        }
};

const char etype::names[ntypes][40] = {"plain", "comm", "partition", "prdc", "symbolic","coupled_symbolic","coupled_symbolic_comm", "spline",
    "circle", "naca","ellipse"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
edge_bdry* tri_mesh::getnewedgeobject(int idnum, input_map& bdrydata) {
    std::string keyword,val;
    std::istringstream data;
    ostringstream nstr;
    int type;          
    edge_bdry *temp;  
    
    type = etype::plain;

    nstr.str("");
    nstr << idnum << std::flush;        
    keyword = gbl->idprefix +"_s" +nstr.str() + "_type";
    if (bdrydata.get(keyword,val)) {
        type = etype::getid(val.c_str());
        if (type < 0)  {
            *gbl->log << "unknown edge type:" << val << std::endl;
            exit(1);
        }
    }
    else {
        type = etype::plain;
    }

    switch(type) {
        case etype::plain: {
            temp = new edge_bdry(idnum,*this);
            break;
        }
        case etype::comm: {
            temp = new ecomm(idnum,*this); 
            break;
        }
        case etype::partition: {
            temp = new epartition(idnum,*this);
            break;
        }
        case etype::prdc: {
            temp = new eprdc(idnum,*this);
            break;
        }
        case etype::symbolic: {
            temp = new eboundary_with_geometry<edge_bdry,symbolic_shape<2> >(idnum,*this);
            break;
        }
        case etype::coupled_symbolic: {
            temp = new ecoupled_physics<eboundary_with_geometry<edge_bdry,symbolic_shape<2> > >(idnum,*this);
            break;
        }
        case etype::coupled_symbolic_comm: {
            temp = new ecoupled_physics<eboundary_with_geometry<ecomm,symbolic_shape<2> > >(idnum,*this);
            break;
        }
        case etype::spline: {
     //      temp = new spline(idnum,*this);
            break;
        }

        /* SPECIAL CASES FOLLOW (DEPRECATED -- USE SYMBOLIC) */
        case etype::circle: {
            temp = new eboundary_with_geometry<edge_bdry,circle>(idnum,*this);
            break;
        }
        case etype::naca: {
            temp = new eboundary_with_geometry<edge_bdry,naca>(idnum,*this);
            break;
        }
        case etype::ellipse: {
            temp = new eboundary_with_geometry<edge_bdry,ellipse>(idnum,*this);
            break;
        }
                
        default: {
            temp = new edge_bdry(idnum,*this);
            std::cout << "unrecognizable edge type: " << idnum << "type " << type << std::endl;
            break;
        }
    }
    
    temp->input(bdrydata);
    
    return(temp);
}

    
    
