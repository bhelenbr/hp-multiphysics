/*
 *  bdry_buoyancy.h
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 12/13/10.
 *  Copyright 2010 Clarkson University. All rights reserved.
 *
 */

#ifndef _melt_buoyancy_h_
#define _melt_buoyancy_h_


#include "tri_hp_buoyancy.h"
#include "../ins/surface2.h"
#include "../hp_coupled_boundary.h"
#include "../cd/melt_cd.h"
#include <myblas.h>
#include <blitz/array.h>
#include <symbolic_function.h>

using namespace blitz;

namespace bdry_buoyancy {
	class surface_marangoni : public bdry_ins::surface2 {
	public:
		symbolic_function<1> sigma_vs_T;
		
		surface_marangoni(tri_hp_ins &xin, edge_bdry &bin) : bdry_ins::surface2(xin,bin) { mytype = "surface_marangoni";}
		surface_marangoni(const surface_marangoni& inbdry, tri_hp_ins &xin, edge_bdry &bin)  : bdry_ins::surface2(inbdry,xin,bin), sigma_vs_T(inbdry.sigma_vs_T) {}
		surface_marangoni* create(tri_hp& xin, edge_bdry &bin) const {return new surface_marangoni(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
		void init(input_map& input,void* gbl_in);
		void element_rsdl(int sind, Array<TinyVector<FLT,MXTM>,1> lf);
	};
	
	class melt_buoyancy : public bdry_cd::melt_cd {
	public:
		tri_hp_buoyancy &x;
		
		melt_buoyancy(tri_hp_buoyancy &xin, edge_bdry &bin) : melt_cd(xin,bin), x(xin) {
			mytype = "melt_buoyancy";
		}
		melt_buoyancy(const melt_buoyancy& inbdry, tri_hp_buoyancy &xin, edge_bdry &bin)  : melt_cd(inbdry,xin,bin), x(xin) {}
		melt_buoyancy* create(tri_hp& xin, edge_bdry &bin) const {return new melt_buoyancy(*this,dynamic_cast<tri_hp_buoyancy&>(xin),bin);}
		
		void element_rsdl(int sind, Array<TinyVector<FLT,MXTM>,1> lf);
		FLT calculate_kinetic_coefficients(FLT DT,FLT sint);
		void output(const std::string& fname, tri_hp::filetype typ,int tlvl);
	};
	
	
	class melt_facet_pt2 : public hp_deformable_free_pnt {
	protected:
		tri_hp_buoyancy &x;
	public:
		melt_facet_pt2(tri_hp_buoyancy &xin, vrtx_bdry &bin) : hp_deformable_free_pnt(xin,bin), x(xin) {mytype = "melt_facet_pt";}
		melt_facet_pt2(const melt_facet_pt2& inbdry, tri_hp_buoyancy &xin, vrtx_bdry &bin) : hp_deformable_free_pnt(inbdry,xin,bin), x(xin) {}
		melt_facet_pt2* create(tri_hp& xin, vrtx_bdry &bin) const {return new melt_facet_pt2(*this,dynamic_cast<tri_hp_buoyancy&>(xin),bin);}
		
		void rsdl(int stage);
		void element_rsdl(Array<FLT,1> lf);
#ifdef petsc
		void petsc_jacobian();
#endif
	};
	
	class triple_junction : public melt_facet_pt2 {
	protected:
		FLT growth_angle;
	public:
		triple_junction(tri_hp_buoyancy &xin, vrtx_bdry &bin) : melt_facet_pt2(xin,bin) {mytype = "triple_junction";}
		triple_junction(const triple_junction& inbdry, tri_hp_buoyancy &xin, vrtx_bdry &bin) : melt_facet_pt2(inbdry,xin,bin) {}
		triple_junction* create(tri_hp& xin, vrtx_bdry &bin) const {return new triple_junction(*this,dynamic_cast<tri_hp_buoyancy&>(xin),bin);}
		
		void init(input_map& input,void* gbl_in);
		void element_rsdl(Array<FLT,1> lf);
#ifdef petsc
		void petsc_jacobian();
#endif
	};
	
	
	
}

#endif // _melt_buoyancy_h_
