//
//  surface.h
//  tri_hp
//
//  Created by Brian Helenbrook on 12/30/14.
//
//

#ifndef _surface_h_
#define _surface_h_

#include <iostream>
#include "../hp_coupled_boundary.h"
#include "tri_hp_ins.h"

namespace bdry_ins {

//#define MPDEBUG
//#define DEBUG

class surface : public hp_coupled_bdry {
	protected:
		tri_hp_ins &x;
	
	public:
		struct surface_global {
			/* FLUID PROPERTIES */
			FLT sigma,rho2,mu2,p_ext;
#ifdef VOLUMEFLUX
            FLT area, area_now;
#endif
        };
        shared_ptr<surface_global> surface_gbl;
		surface(tri_hp_ins &xin, edge_bdry &bin) : hp_coupled_bdry(xin,bin), x(xin) {mytype = "surface";}
		surface(const surface& inbdry, tri_hp_ins &xin, edge_bdry &bin)  : hp_coupled_bdry(inbdry,xin,bin), x(xin), surface_gbl(inbdry.surface_gbl) {}
		surface* create(tri_hp& xin, edge_bdry &bin) const {return new surface(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
		
		void init(input_map& inmap);
		void element_rsdl(int sind, Array<TinyVector<FLT,MXTM>,1> lf);
		int setup_preconditioner();
#ifndef petsc
		// void rsdl(int stage); //!< applies mass flux preconditioner for multgrid case
#endif
};

class surface_outflow : public hp_deformable_free_pnt {
		/* For a surface point sliding on a vertical or horizontal surface */
		/* For periodic wall have tri_mesh vertex type be comm */
	protected:
		FLT contact_angle;
		
	public:
		surface_outflow(tri_hp_ins &xin, vrtx_bdry &bin) : hp_deformable_free_pnt(xin,bin) {mytype = "surface_outflow";}
		surface_outflow(const surface_outflow& inbdry, tri_hp_ins &xin, vrtx_bdry &bin) : hp_deformable_free_pnt(inbdry,xin,bin), contact_angle(inbdry.contact_angle) {}
		surface_outflow* create(tri_hp& xin, vrtx_bdry &bin) const {return new surface_outflow(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
		
		void init(input_map& inmap);
	
		/* Routine to add surface tension stress */
		void element_rsdl(Array<FLT,1> lf);
	};
}
#endif /* defined(_surface_h_) */
