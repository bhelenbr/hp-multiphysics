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

#ifdef ROTATE_RESIDUALS
#define SWAP_ROWS
#endif

// Put shared routines in this class and let buoyancy class add extra things
namespace bdry_cd {
	class melt_cd : public hp_deformable_bdry {
		
	public:
		struct global : public hp_deformable_bdry::global {
			/* PROPERTIES */
			FLT Lf, rho_s, cp_s, k_s, rho_l, cp_l, k_l;
			
			/* Kinetic Coefficients */
			FLT Krough, Ksn, A2Dn, K2Dn, K2Dn_max; // Coefficients for Weinstein Kinetic model
			FLT Kgt; // Gibbs Thompson curvature effect (not working)
			TinyVector<FLT,tri_mesh::ND> facetdir; // Diretion of facet
			FLT surge_time;
			
			Array<FLT,1> vdt_kinetic;
			Array<FLT,1> sdt_kinetic;
		} *gbl;
		
		public:
			melt_cd(tri_hp &xin, edge_bdry &bin) : hp_deformable_bdry(xin,bin) {mytype = "melt_cd";}
			melt_cd(const melt_cd& inbdry, tri_hp &xin, edge_bdry &bin) : hp_deformable_bdry(inbdry,xin,bin) {}
			melt_cd* create(tri_hp& xin, edge_bdry &bin) const {return new melt_cd(*this,xin,bin);}
			void* create_global_structure() {return new global;}

			void init(input_map& inmap, void* gbl_in);
			void setup_preconditioner();
		
#if defined(petsc) && defined(SWAP_ROWS)
		/* routines to swap residuals between kinetic equation and energy equation */
		void non_sparse(Array<int,1> &nnzero);
		void non_sparse_rcv(Array<int, 1> &nnzero, Array<int, 1> &nnzero_mpi);
		void petsc_make_1D_rsdl_vector(Array<double,1> res);
		void petsc_premultiply_jacobian();
#endif
	};
}
#endif

