/*
 *  getnewibc.cpp
 *  tet_hp
 *
 *  Created by Brian Helenbrook on 2/6/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp_cns.h"
#include "bdry_cns.h"
#include <tet_boundary.h>

namespace ibc_cns {

	class freestream : public init_bdry_cndtn {
		private:
			FLT alpha, speed,perturb_amp;

		public:
			FLT f(int n, TinyVector<FLT,tet_mesh::ND> x, FLT time) {
				FLT amp = (time > 0.0 ? 0.0 : perturb_amp); 
				switch(n) {
					case(0):
						return(0.7);
					case(1):
						return(speed/340.0*sin(alpha));
					case(2):
						return(speed/340.0*cos(alpha));
					case(3):
						return(0.7);
				}
				return(0.0);
			}

			void input(input_map &blockdata,std::string idnty) {
				std::string keyword,val;
				std::istringstream data;

				keyword = idnty +"_flowspeed";
				if (!blockdata.get(keyword,speed)) 
						blockdata.getwdefault("flowspeed",speed,1.0);

				keyword = idnty +"_flowangle";
				if (!blockdata.get(keyword,alpha)) 
						blockdata.getwdefault("flowangle",alpha,0.0);  

				keyword = idnty +"_perturb_amplitude";
				if (!blockdata.get(keyword,perturb_amp)) 
						blockdata.getwdefault("perturb_amplitude",perturb_amp,0.0); 

				alpha *= M_PI/180.0;
			}
	};

	class ibc_type {
		public:
			const static int ntypes = 1;
			enum ids {freestream};
			const static char names[ntypes][40];
			static int getid(const char *nin) {
				int i;
				for(i=0;i<ntypes;++i) 
						if (!strcmp(nin,names[i])) return(i);
				return(-1);
		}
	};
	const char ibc_type::names[ntypes][40] = {"freestream"};
	
	
	class helper_type {
	public:
		const static int ntypes = 4;
		enum ids {translating_drop,parameter_changer,unsteady_body_force,force_coupling};
		const static char names[ntypes][40];
		static int getid(const char *nin) {
			int i;
			for(i=0;i<ntypes;++i) 
				if (!strcmp(nin,names[i])) return(i);
			return(-1);
		}
	};
	const char helper_type::names[ntypes][40] = {"translating_drop","parameter_changer","unsteady_body_force","force_coupling"};
	

}

tet_hp_helper *tet_hp_cns::getnewhelper(input_map& inmap) {
	std::string keyword,movername;
	int type;
	
	/* FIND INITIAL CONDITION TYPE */
	keyword = std::string(gbl->idprefix) + "_helper";
	if (!inmap.get(keyword,movername)) {
		if (!inmap.get("helper",movername)) {
			type = -1;
		}
	}
	
	type = ibc_cns::helper_type::getid(movername.c_str());
	
	switch(type) {
//		case ibc_cns::helper_type::parameter_changer: {
//			tet_hp_helper *temp = new ibc_cns::parameter_changer(*this);
//			return(temp);
//		}

		default: {
			return(tet_hp::getnewhelper(inmap));
		}
	}
}


init_bdry_cndtn *tet_hp_cns::getnewibc(std::string suffix, input_map& inmap) {
	std::string keyword,ibcname;
	init_bdry_cndtn *temp;
	int type;

	/* FIND INITIAL CONDITION TYPE */
	keyword = gbl->idprefix + "_" +suffix;
	if (!inmap.get(keyword,ibcname)) {
		keyword = suffix;
		if (!inmap.get(keyword,ibcname)) {
			*gbl->log << "couldn't find initial condition type" << std::endl;
		}
	}
	type = ibc_cns::ibc_type::getid(ibcname.c_str());

		
	switch(type) {
		case ibc_cns::ibc_type::freestream: {
			temp = new ibc_cns::freestream;
			break;
		}
		default: {
			return(tet_hp::getnewibc(suffix,inmap));
		}
	}
	temp->input(inmap,keyword);
	return(temp);
}






