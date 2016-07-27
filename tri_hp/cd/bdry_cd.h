/*
 *  cd_bdry.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _bdry_cd_h_
#define _bdry_cd_h_

#include "tri_hp_cd.h"
#include "../hp_boundary.h"
#include "myblas.h"
#include <symbolic_function.h>

namespace bdry_cd {
	class generic : public hp_edge_bdry {
		protected:
			tri_hp_cd &x;
			FLT diff_total,conv_total,circumference;
			
		public:
			generic(tri_hp_cd &xin, edge_bdry &bin) : hp_edge_bdry(xin,bin), x(xin) {mytype = "generic";}
			generic(const generic& inbdry, tri_hp_cd &xin, edge_bdry &bin) : hp_edge_bdry(inbdry,xin,bin), x(xin) {}
			generic* create(tri_hp& xin, edge_bdry &bin) const {return new generic(*this,dynamic_cast<tri_hp_cd&>(xin),bin);}
			void output(const std::string& fname, tri_hp::filetype typ,int tlvl = 0);
	};

	class dirichlet : public generic {
		public:
			dirichlet(tri_hp_cd &xin, edge_bdry &bin) : generic(xin,bin) {mytype = "dirichlet";}
			dirichlet(const dirichlet& inbdry, tri_hp_cd &xin, edge_bdry &bin) : generic(inbdry,xin,bin) {}
			dirichlet* create(tri_hp& xin, edge_bdry &bin) const {return new dirichlet(*this,dynamic_cast<tri_hp_cd&>(xin),bin);}
			void init(input_map& inmap,void* gbl_in) {
				generic::init(inmap,gbl_in);
				type[0] = essential;
				essential_indices.push_back(0);
			}
	};

	class characteristic : public generic {
		public:
			void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx) {
				FLT vel;

#ifdef CONST_A
				vel =  (x.gbl->ax-mv(0))*norm(0) +(x.gbl->ay -mv(1))*norm(1);        
#else
				vel =  (x.gbl->a->f(0,pt,x.gbl->time)-mv(0))*norm(0) +(x.gbl->a->f(1,pt,x.gbl->time) -mv(1))*norm(1);
#endif

				if (vel > 0.0)
					flx(0) = x.gbl->rhocv*vel*u(0);
				else
					flx(0) = x.gbl->rhocv*vel*ibc->f(0, xpt, x.gbl->time);
			}
			characteristic(tri_hp_cd &xin, edge_bdry &bin) : generic(xin,bin) {mytype = "characteristic";}
			characteristic(const characteristic &inbdry, tri_hp_cd &xin, edge_bdry &bin) : generic(inbdry,xin,bin) {}
			characteristic* create(tri_hp& xin, edge_bdry &bin) const {return new characteristic(*this,dynamic_cast<tri_hp_cd&>(xin),bin);}
	};
}

#endif
