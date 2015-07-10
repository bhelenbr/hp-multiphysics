/*
 *  mvpttobdry.cpp
 *  mblock
 *
 *  Created by helenbrk on Wed Oct 24 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_mesh.h"
#include "tet_boundary.h"

/** \brief Helper object for vrtx_bdry 
 *
 * \ingroup boundary
 * Contains list of all vrtx_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class vtype {
	public:
		static const int ntypes = 4;
		enum ids {plain=1,comm,prdc,symbolic};
		const static char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i) 
				if (!strcmp(nin,names[i])) return(i+1);
			return(-1);
		}
};

const char vtype::names[ntypes][40] = {"plain","comm","prdc","symbolic"};

vrtx_bdry* tet_mesh::getnewvrtxobject(int idnum, input_map& in_map) {
	std::string keyword,typ_str;
	ostringstream nstr;
	int type;          
	vrtx_bdry *temp;  

	nstr.str("");
	nstr << idnum << std::flush;        
	keyword = gbl->idprefix +"_v" +nstr.str() + "_type";
	
	in_map.getwdefault(keyword,typ_str,std::string("plain"));
	type = vtype::getid(typ_str.c_str());
	if (type < 0)  {
		*gbl->log << "unknown vertex type:" << typ_str << std::endl;
		exit(1);
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
		case vtype::symbolic: {
			temp = new vboundary_with_geometry<vrtx_bdry,symbolic_point<tet_mesh::ND> >(idnum,*this);
			break;
		}
		default: {
			std::cout << "unrecognizable vrtx type: " <<  type << " idnum: " << idnum << std::endl;
			temp = new vrtx_bdry(idnum,*this);
			break;
		}
	} 
	
	temp->init(in_map);
	
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
		enum ids {plain=1, comm, partition, prdc, symbolic, symbolic_comm, coupled_symbolic, 
			coupled_symbolic_comm, spline,circle, naca, ellipse};
		static const char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i)
				if (!strcmp(nin,names[i])) return(i+1);
			return(-1);
		}
};

const char etype::names[ntypes][40] = {"plain", "comm", "partition", "prdc", "symbolic","symbolic_comm",
	"coupled_symbolic","coupled_symbolic_comm", "spline","circle", "naca","ellipse"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
edge_bdry* tet_mesh::getnewedgeobject(int idnum, input_map& in_map) {
	std::string keyword,typ_str;
	ostringstream nstr;
	int type;          
	edge_bdry *temp;  
	
	type = etype::plain;

	nstr.str("");
	nstr << idnum << std::flush;        
	keyword = gbl->idprefix +"_e" +nstr.str() + "_type";
	
	in_map.getwdefault(keyword,typ_str,std::string("plain"));
	type = etype::getid(typ_str.c_str());
	if (type < 0)  {
		*gbl->log << "unknown edge type:" << typ_str << std::endl;
		exit(1);
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
			temp = new eboundary_with_geometry<edge_bdry,symbolic_shape<tet_mesh::ND> >(idnum,*this);
			break;
		}
		case etype::symbolic_comm: {
			temp = new eboundary_with_geometry<ecomm,symbolic_shape<tet_mesh::ND> >(idnum,*this);
			break;
		}
		case etype::coupled_symbolic: {
			temp = new ecoupled_physics<eboundary_with_geometry<edge_bdry,symbolic_shape<tet_mesh::ND> > >(idnum,*this);
			break;
		}
		case etype::coupled_symbolic_comm: {
			temp = new ecoupled_physics<eboundary_with_geometry<ecomm,symbolic_shape<tet_mesh::ND> > >(idnum,*this);
			break;
		}
		case etype::spline: {
			temp = new spline_bdry(idnum,*this);
			break;
		}

		/* SPECIAL CASES FOLLOW (DEPRECATED -- USE SYMBOLIC) */
//        case etype::circle: {
//            temp = new eboundary_with_geometry<edge_bdry,circle>(idnum,*this);
//            break;
//        }
//        case etype::naca: {
//            temp = new eboundary_with_geometry<edge_bdry,naca>(idnum,*this);
//            break;
//        }
//        case etype::ellipse: {
//            temp = new eboundary_with_geometry<edge_bdry,ellipse>(idnum,*this);
//            break;
//        }
				
		default: {
			temp = new edge_bdry(idnum,*this);
			std::cout << "unrecognizable edge type: " << idnum << "type " << type << std::endl;
			break;
		}
	}
	
	temp->init(in_map);
	
	return(temp);
}

/** \brief Helper object for face_bdry 
 *
 * \ingroup boundary
 * Contains list of all side_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class ftype {
	public:
		static const int ntypes = 6;
		enum ids {plain=1, comm, prdc, partition,symbolic,symbolic_comm};
		static const char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i)
				if (!strcmp(nin,names[i])) return(i+1);
			return(-1);
		}
};

const char ftype::names[ntypes][40] = {"plain","comm","prdc","partition","symbolic","symbolic_comm"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
face_bdry* tet_mesh::getnewfaceobject(int idnum, input_map& inmap) {
	std::string keyword,val;
	std::istringstream data;
	ostringstream nstr;
	int type;          
	face_bdry *temp;  
	
	type = ftype::plain;

	nstr.str("");
	nstr << idnum << std::flush;        
	keyword = gbl->idprefix +"_f" +nstr.str() + "_type";
	if (inmap.get(keyword,val)) {
		type = ftype::getid(val.c_str());
		if (type < 0)  {
			*gbl->log << "unknown face type:" << val << std::endl;
			exit(1);
		}
	}
	else {
		type = ftype::plain;
	}

	switch(type) {
		case ftype::plain: {
			temp = new face_bdry(idnum,*this);
			break;
		}
		case ftype::comm: {
			temp = new fcomm(idnum,*this);
			break;
		}
		case ftype::prdc: {
			temp = new fprdc(idnum,*this);
			break;
		}
		case ftype::partition: {
			temp = new fpartition(idnum,*this);
			break;
		}
		case ftype::symbolic: {
			temp = new fboundary_with_geometry<face_bdry,symbolic_shape<tet_mesh::ND> >(idnum,*this);
			break;
		}
		case ftype::symbolic_comm: {
			temp = new fboundary_with_geometry<fcomm,symbolic_shape<tet_mesh::ND> >(idnum,*this);
			break;
		}
		default: {
			temp = new face_bdry(idnum,*this);
			std::cout << "unrecognizable face type: " << idnum << "type " << type << std::endl;
			break;
		}
	}
	
	temp->init(inmap);
	
	return(temp);
}

	
	
