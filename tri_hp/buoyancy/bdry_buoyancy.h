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
#include <blitz/array.h>
#include <symbolic_function.h>

using namespace blitz;

namespace bdry_buoyancy {
	
	class solid_fluid : public hp_edge_bdry {
	public:
		solid_fluid(tri_hp_buoyancy &xin, edge_bdry &bin) : hp_edge_bdry(xin,bin) {
			mytype = "solid_fluid";
		}
		solid_fluid(const solid_fluid& inbdry, tri_hp_ins &xin, edge_bdry &bin)  : hp_edge_bdry(inbdry,xin,bin) {}
		solid_fluid* create(tri_hp& xin, edge_bdry &bin) const {return new solid_fluid(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
		void init(input_map& inmap);
	};

	/* Same thing as incompressible characteristic B.C. except adjusted for variable density */
	class characteristic : public bdry_ins::characteristic {
		protected:
			tri_hp_buoyancy& x;

			void flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx) {
				FLT store_rho = x.hp_ins_gbl->rho;
				x.hp_ins_gbl->rho = x.hp_buoyancy_gbl->rho_vs_T.Eval(u(2));
				bdry_ins::characteristic::flux(u,xpt,mv,norm,side_length,flx);
				flx(2) *= x.hp_buoyancy_gbl->cp;
				x.hp_ins_gbl->rho = store_rho;
			}
		public:
			characteristic(tri_hp_buoyancy &xin, edge_bdry &bin) : x(xin), bdry_ins::characteristic(xin,bin) {mytype = "characteristic";}
			characteristic(const characteristic& inbdry, tri_hp_buoyancy &xin, edge_bdry &bin) : x(xin), bdry_ins::characteristic(inbdry,xin,bin) {}
			characteristic* create(tri_hp& xin, edge_bdry &bin) const {return new characteristic(*this,dynamic_cast<tri_hp_buoyancy&>(xin),bin);}
	};
}
#endif

