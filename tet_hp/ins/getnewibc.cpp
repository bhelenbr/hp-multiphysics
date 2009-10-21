/*
 *  getnewibc_ins.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 2/6/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp_ins.h"
#include "bdry_ins.h"
#include <tet_boundary.h>

namespace ibc_ins {

	class freestream : public init_bdry_cndtn {
		private:
			FLT alpha, speed,perturb_amp;

		public:
			FLT f(int n, TinyVector<FLT,tet_mesh::ND> x, FLT time) {
				FLT amp = (time > 0.0 ? 0.0 : perturb_amp); 
				switch(n) {
					case(0):
						return(speed*cos(alpha) + amp*(0.5-x(0))*(0.5+x(0))*(0.5-x(1))*(0.5+x(1))*(0.5-x(2))*(0.5+x(2)));
						//return(speed*cos(alpha) + amp*(1.0-x(0))*(x(0))*(1.0-x(1))*(x(1))*(1.0-x(2))*(x(2)));

						//return(speed*cos(alpha) + .25*amp*exp(-x(0)*x(0)*100)*exp(-x(1)*x(1)*100)*exp(-x(2)*x(2)*100));
						//return(speed*cos(alpha) + .25*amp*exp(-(x(0)-0.5)*(x(0)-0.5)*100)*exp(-(x(1)-0.5)*(x(1)-0.5)*100)*exp(-(x(2)-0.5)*(x(2)-0.5)*100));
						//return(0);
					case(1):
						return(speed*sin(alpha));
					case(2):
						return(0.0);
						//return(16/0.08333333333/0.08333333333/0.1875/0.1875*(0.1875/2.0+x(0)-0.375)*(0.1875/2.0-x(0)+0.375)*(0.08333333333/2.0+x(1)-0.104166666665)*(0.08333333333/2.0-x(1)+0.104166666665));
					case(3):
						return(0.0);
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

}


init_bdry_cndtn *tet_hp_ins::getnewibc(std::string suffix, input_map& inmap) {
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
	type = ibc_ins::ibc_type::getid(ibcname.c_str());

		
	switch(type) {
		case ibc_ins::ibc_type::freestream: {
			temp = new ibc_ins::freestream;
			break;
		}
		default: {
			return(tet_hp::getnewibc(suffix,inmap));
		}
	}
	temp->input(inmap,keyword);
	return(temp);
}






