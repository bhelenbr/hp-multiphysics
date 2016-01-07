/*
 *  cd_bdry.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _melt_cd_h_
#define _melt_cd_h_

#include "tri_hp_cd.h"
#include "../hp_coupled_boundary.h"
#include "myblas.h"
#include <symbolic_function.h>

//#define OLDKINETICS

// Put shared routines in this class and let buoyancy class add extra things
namespace bdry_cd {
	class melt_cd : public hp_coupled_bdry {
		
	public:
		struct global : public hp_coupled_bdry::global {
			/* PROPERTIES */
			FLT Lf, rho_s, cp_s, k_s, rho_l, cp_l, k_l;
			
			/* Kinetic Coefficients */
			FLT Krough, Ksn, A2Dn, K2Dn;
#ifdef OLDKINETICS
			FLT K2Dn_max; // Coefficients for Weinstein Kinetic model
#else
			FLT K2Dn_DT_min; // Minimum supercooling DT for nucleation
#endif
			FLT Kgt; // Gibbs Thompson curvature effect (not working)
			TinyVector<FLT,tri_mesh::ND> facetdir; // Diretion of facet
			FLT surge_time;
			
			Array<FLT,1> vdt_kinetic;
			Array<FLT,1> sdt_kinetic;
		} *gbl;
		
	public:
		melt_cd(tri_hp &xin, edge_bdry &bin) : hp_coupled_bdry(xin,bin) {mytype = "melt_cd";}
		melt_cd(const melt_cd& inbdry, tri_hp &xin, edge_bdry &bin) : hp_coupled_bdry(inbdry,xin,bin) {}
		melt_cd* create(tri_hp& xin, edge_bdry &bin) const {return new melt_cd(*this,xin,bin);}
		void* create_global_structure() {return new global;}
		
		void init(input_map& inmap, void* gbl_in);
		void output(const std::string& filename, tri_hp::filetype typ,int tlvl);
		void setup_preconditioner();
		void element_rsdl(int sind, Array<TinyVector<FLT,MXTM>,1> lf);
		FLT calculate_kinetic_coefficients(FLT DT,FLT sint);
	};
	
	class melt_facet_pt2 : public hp_deformable_free_pnt {
	public:
		melt_facet_pt2(tri_hp &xin, vrtx_bdry &bin) : hp_deformable_free_pnt(xin,bin) {mytype = "melt_facet_pt";}
		melt_facet_pt2(const melt_facet_pt2& inbdry, tri_hp &xin, vrtx_bdry &bin) : hp_deformable_free_pnt(inbdry,xin,bin) {}
		melt_facet_pt2* create(tri_hp& xin, vrtx_bdry &bin) const {return new melt_facet_pt2(*this,xin,bin);}
		
		void rsdl(int stage);
		void element_rsdl(Array<FLT,1> lf);
#ifdef petsc
		void petsc_jacobian();
#endif
	};
}
#endif

