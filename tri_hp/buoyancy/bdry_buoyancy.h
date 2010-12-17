/*
 *  bdry_buoyancy.h
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 12/13/10.
 *  Copyright 2010 Clarkson University. All rights reserved.
 *
 */
 
#ifndef _bdry_buoyancy_h_
#define _bdry_buoyancy_h_


#include "tri_hp_buoyancy.h"
#include "../ins/bdry_ins.h"
#include <myblas.h>
#include <blitz/tinyvec-et.h>
#include <blitz/array.h>
#include <symbolic_function.h>


using namespace blitz;

namespace bdry_buoyancy {
	
	class surface : public bdry_ins::surface {
		Array<vector_function,1> fluxes;
											
		public:
			surface(tri_hp_ins &xin, edge_bdry &bin) : bdry_ins::surface(xin,bin) {mytype = "surface"; fluxes.resize(x.NV);
				Array<string,1> names(4);
				Array<int,1> dims(4);
				dims = x.ND;
				names(0) = "u";
				dims(0) = x.NV;
				names(1) = "x";
				names(2) = "xt";
				names(3) = "n";
				for(int n=0;n<x.NV;++n) {
					fluxes(n).set_arguments(4,dims,names);
				}
			}
			surface(const surface& inbdry, tri_hp_ins &xin, edge_bdry &bin)  : bdry_ins::surface(inbdry,xin,bin), fluxes(inbdry.fluxes) {}
			surface* create(tri_hp& xin, edge_bdry &bin) const {return new surface(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
			
			void init(input_map& input,void* gbl_in); 
			void element_rsdl(int sind, Array<FLT,2> lf);  // FIXME Not really compatible need to make all consistent
	};
}
#endif
