		
/*
 *  getnewibc_ins.cpp
 *  tet_hp
 *
 *  Created by Brian Helenbrook on 2/6/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp_cd.h"
#include <spline.h>

namespace ibc_cd {
	
	class zero_src : public init_bdry_cndtn {
		public:
			FLT f(int n, TinyVector<FLT,tet_mesh::ND> x,FLT time) {
				return(0.0);
			}
			void input(input_map &inmap,std::string idnty) {}
	};
	
	class unsteady_spline : public init_bdry_cndtn {
		spline<1> time_history;
		public:
			FLT f(int n, TinyVector<FLT,tet_mesh::ND> x,FLT time) {
				TinyVector<FLT,1> source;
				time_history.interpolate(time, source);
				return(source(0));
			}
			void input(input_map &inmap,std::string idnty) {
				std::string filename;
				if (!inmap.get(idnty+"_filename",filename)) {
					std::cerr << "couldn't find unsteady_spline filename" << std::endl;
					exit(1);
				}
				time_history.read(filename);
			}
	};

	class ibc_type {
		public:
			const static int ntypes = 2;
			enum ids {zero,unsteady_spline};
			const static char names[ntypes][40];
			static int getid(const char *nin) {
			int i;
			for(i=0;i<ntypes;++i) 
				if (!strcmp(nin,names[i])) return(i);
			return(-1);
			}
	};
	const char ibc_type::names[ntypes][40] = {"zero","unsteady_spline"};
}


init_bdry_cndtn *tet_hp_cd::getnewibc(std::string suffix, input_map& inmap) {
	std::string keyword,ibcname;
	init_bdry_cndtn *temp;
	int type;

	/* FIND INITIAL CONDITION TYPE */
	keyword = gbl->idprefix + "_" +suffix;
	if (!inmap.get(keyword,ibcname)) {
		keyword = suffix;
		if (!inmap.get(keyword,ibcname)) {
			*gbl->log << "couldn't find cd initial condition type " << keyword << std::endl;
		}
	}
	type = ibc_cd::ibc_type::getid(ibcname.c_str());

	switch(type) {
		case(ibc_cd::ibc_type::zero): {
			temp = new ibc_cd::zero_src;
			break;
		}
		case(ibc_cd::ibc_type::unsteady_spline): {
			temp = new ibc_cd::unsteady_spline;
			break;
		}
		default: {
			return(tet_hp::getnewibc(suffix,inmap));
		}
	}
	temp->input(inmap,keyword);
	return(temp);
}
