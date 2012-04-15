/*
 *  mvpttobdry.cpp
 *  mblock
 *
 *  Created by helenbrk on Wed Oct 24 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_mesh.h"
#include "tri_boundary.h"

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

vrtx_bdry* tri_mesh::getnewvrtxobject(int idnum, input_map& in_map) {
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
		sim::abort(__LINE__,__FILE__,gbl->log);
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
			temp = new vboundary_with_geometry<vrtx_bdry,symbolic_point<tri_mesh::ND> >(idnum,*this);
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
		static const int ntypes = 16;
		enum ids {plain=1, comm, partition, prdc, symbolic, symbolic_comm, coupled_symbolic,
			coupled_symbolic_comm, spline, spline_comm, coupled_spline, coupled_spline_comm, circle, naca, ellipse,planar};
		static const char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i)
				if (!strcmp(nin,names[i])) return(i+1);
			return(-1);
		}
};

const char etype::names[ntypes][40] = {"plain", "comm", "partition", "prdc", "symbolic","symbolic_comm",
	"coupled_symbolic","coupled_symbolic_comm", "spline","spline_comm","coupled_spline","coupled_spline_comm","circle", "naca","ellipse","planar"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
edge_bdry* tri_mesh::getnewedgeobject(int idnum, input_map& in_map) {
	std::string keyword,typ_str;
	ostringstream nstr;
	int type;
	edge_bdry *temp;

	type = etype::plain;

	nstr.str("");
	nstr << idnum << std::flush;
	keyword = gbl->idprefix +"_s" +nstr.str() + "_type";

	in_map.getwdefault(keyword,typ_str,std::string("plain"));
	type = etype::getid(typ_str.c_str());
	if (type < 0)  {
		*gbl->log << "unknown edge type:" << typ_str << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
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
			temp = new eboundary_with_geometry<edge_bdry,symbolic_shape<tri_mesh::ND> >(idnum,*this);
			break;
		}
		case etype::symbolic_comm: {
			temp = new eboundary_with_geometry<ecomm,symbolic_shape<tri_mesh::ND> >(idnum,*this);
			break;
		}
		case etype::coupled_symbolic: {
			temp = new ecoupled_physics<eboundary_with_geometry<edge_bdry,symbolic_shape<tri_mesh::ND> > >(idnum,*this);
			break;
		}
		case etype::coupled_symbolic_comm: {
			temp = new ecoupled_physics<eboundary_with_geometry<ecomm,symbolic_shape<tri_mesh::ND> > >(idnum,*this);
			break;
		}
		case etype::spline: {
			temp = new spline_bdry(idnum,*this);
			break;
		}
		case etype::spline_comm: {
			temp = new eboundary_with_geometry<ecomm,spline_geometry>(idnum,*this);
			break;
		}
		case etype::coupled_spline: {
			temp = new ecoupled_physics<eboundary_with_geometry<edge_bdry,spline_geometry> >(idnum,*this);
			break;
		}
		case etype::coupled_spline_comm: {
			temp = new ecoupled_physics<eboundary_with_geometry<ecomm,spline_geometry> >(idnum,*this);
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
		case etype::planar: {
			temp = new eboundary_with_geometry<edge_bdry,plane>(idnum,*this);
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

	temp->init(in_map);

	return(temp);
}



