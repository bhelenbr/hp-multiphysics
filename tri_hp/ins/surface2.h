//
//  surface2.h
//  tri_hp
//
//  Created by Brian Helenbrook on 12/30/14.
//
//

#ifndef _surface2_h_
#define _surface2_h_

#include <iostream>
#include "../hp_coupled_boundary.h"
#include "tri_hp_ins.h"

namespace bdry_ins {

//#define MPDEBUG
//#define DEBUG

class surface2 : public hp_deformable_bdry {
	protected:
		tri_hp_ins &x;
	
	public:
		struct global : public hp_deformable_bdry::global {
			/* FLUID PROPERTIES */
			FLT sigma,rho2,mu2,p_ext;
		} *gbl;
		
		void* create_global_structure() {return new global;}
		surface2(tri_hp_ins &xin, edge_bdry &bin) : hp_deformable_bdry(xin,bin), x(xin) {mytype = "surface2";}
		surface2(const surface2& inbdry, tri_hp_ins &xin, edge_bdry &bin)  : hp_deformable_bdry(inbdry,xin,bin), x(xin) {
			gbl = inbdry.gbl;
		};
		surface2* create(tri_hp& xin, edge_bdry &bin) const {return new surface2(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
		
		void init(input_map& input,void* gbl_in);
		void element_rsdl(int sind, Array<TinyVector<FLT,MXTM>,1> lf);
		void setup_preconditioner();
#ifndef petsc
		// void rsdl(int stage); //!< applies mass flux preconditioner for multgrid case
#endif
};

class surface_outflow2 : public hp_deformable_free_pnt {
		/* For a surface point sliding on a vertical or horizontal surface */
		/* For periodic wall have tri_mesh vertex type be comm */
	protected:
		surface2* surface;
		FLT contact_angle;
		TinyVector<FLT,tri_mesh::ND> wall_normal;
		
	public:
		surface_outflow2(tri_hp_ins &xin, vrtx_bdry &bin) : hp_deformable_free_pnt(xin,bin) {mytype = "surface_outflow2";}
		surface_outflow2(const surface_outflow2& inbdry, tri_hp_ins &xin, vrtx_bdry &bin) : hp_deformable_free_pnt(inbdry,xin,bin), contact_angle(inbdry.contact_angle) {}
		surface_outflow2* create(tri_hp& xin, vrtx_bdry &bin) const {return new surface_outflow2(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
		
		void init(input_map& inmap,void* gbl_in);
	
		/* Routine to add surface tension stress */
		void element_rsdl(Array<FLT,1> lf);
	};
}
#endif /* defined(_surface2_h_) */
