/*
 *  getnewibc.cpp
 *  tri_hp
 *
 *  Created by Erik Durkish on 2/7/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_swe.h"
#include <cstring>

namespace ibc_swe {

	class rossby : public init_bdry_cndtn {
		private:
			FLT amplitude;
		public:
			FLT f(int n, TinyVector<FLT,tri_mesh::ND> x, FLT time) {
				FLT A = 0.771*amplitude*amplitude;
				FLT phi = A*pow(cosh(amplitude*x(0)),-2);
				FLT dphidx = -2*amplitude*tanh(amplitude*x(0))*phi;
				FLT h = 1.0 +phi*(3.+6.*x(1)*x(1))/4.*exp(-x(1)*x(1)/2);
				switch(n) {
					case(0):
						return(h*phi*(-9.+6.*x(1)*x(1))/4.*exp(-x(1)*x(1)/2.));
					case(1):
						return(h*dphidx*2*x(1)*exp(-x(1)*x(1)/2.));
					case(2):
						return(h +0.00*exp(-x(1)*x(1)));

				}
				return(0.0);
			}

			void init(input_map &in_map,std::string idnty) {

				if (!in_map.get(idnty + "_rossby_amp",amplitude)) in_map.getwdefault("rossby_amp",amplitude,0.395);
			}
	};
	
	class flat : public init_bdry_cndtn {
	private:
	public:
		FLT f(int n, TinyVector<FLT,tri_mesh::ND> x, FLT time) {
			return(0.0);
		}
	};

	class ibc_type {
		public:
			const static int ntypes = 2;
			enum ids {rossby};
			const static char names[ntypes][40];
			static int getid(const char *nin) {
				int i;
				for(i=0;i<ntypes;++i) 
					if (!strcmp(nin,names[i])) return(i);
				return(-1);
		}
    };
    const char ibc_swe::ibc_type::names[ntypes][40] = {"rossby","flat"};


};


init_bdry_cndtn *tri_hp_swe::getnewibc(std::string name) {
	std::string ibcname, keyword;
	init_bdry_cndtn *temp;
	int type;

	type = ibc_swe::ibc_type::getid(name.c_str());
	switch(type) {
		case ibc_swe::ibc_type::rossby: {
			temp = new ibc_swe::rossby;
			break;
		}
		case ibc_swe::ibc_type::flat: {
			temp = new ibc_swe::flat;
			break;
		}

		default: {
			return(tri_hp_ins::getnewibc(name));
		}
	}
	return(temp);
}
