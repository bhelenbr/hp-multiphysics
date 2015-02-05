//
//  hp_coupled_boundary.h
//  tri_hp
//
//  Created by Brian Helenbrook on 12/27/14.
//
//

#ifndef __tri_hp__hp_coupled_boundary__
#define __tri_hp__hp_coupled_boundary__

#include "hp_boundary.h"
#include <iostream>

//#define MPDEBUG
//#define DEBUG
//#define DETAILED_MINV

/* Non-deforming coupled boundary with variables on boundary */
class hp_coupled_bdry : public hp_edge_bdry {
    protected:
        int NV; //!< number of manifold variables
        bool is_master; //!< master slave relationship for b.c. pairs

	/* Stuff needed for a surface variable */
//		/** Stores vertex, side coefficients of solution */
//		struct vs {
//			Array<FLT,2> v;
//			Array<FLT,3> s;
//		} ug;
//	
//		/** Array for time history information */
//		Array<vs,1> ugbd;
//		Array<Array<FLT,3>,1> dugdt; //!< Precalculated unsteady sources at Gauss points
		Array<FLT,2> vug_frst; //!< value of solution on entry to multigrid (pnts,NV)
		Array<FLT,3> vdres; //!< Driving term for multigrid (log2p, pnts,NV)
		Array<FLT,4> sdres; //!< Driving term for multigrid (log2p, side, order, NV)
		const hp_coupled_bdry *fine, *coarse;  //!< Pointers to coarse and fine mesh objects for multigrid transfers
		
	public:
		struct global {
			
			bool is_loop;  //!< true if boundary is a loop

			/* SOLUTION STORAGE ON FIRST ENTRY TO NSTAGE */
			Array<FLT,2> vug0; //!< vertex solution on entry to multigrid (pnts,NV)
			Array<FLT,3> sug0; //!< side solution on entry to multigrid (segs,mode,NV)
			
			/* RESIDUALS */
			Array<FLT,2> vres; //!< vertex residuals (pnts,NV)
			Array<FLT,3> sres; //!< side residuals (segs,mode,NV)
			Array<FLT,2> vres0; //!< vertex multigrid driving terms (pnts,NV)
			Array<FLT,3> sres0; //!< side multigrid driving terms (segs,mode,NV)

			/* PRECONDITIONER */
			Array<FLT,3> vdt; //!< Time advancement for vertices (pnts,NV,NV)
			Array<FLT,3> sdt; //!< Time advancement for sides (segs,NV,NV)
			Array<FLT,2> meshc; //!< upwinding constants (segs,NV)
			
#ifdef DETAILED_MINV
			Array<FLT,3> ms;  //!< Side mode mass matrix (segs,NV*sm,NV*sm)
			Array<FLT,5> vms;  //!< remove vertex contribution to side mode residual (segs,NV,endpoint,mode,endpoint)
			Array<int,2> ipiv; //!< Pivot vector for Side mode mass matrix (segs,NV*sm)
#endif
			Array<FLT,1> fadd; //!< Multiplier for driving term of multigrid (NV)
			Array<FLT,2> cfl; //!< CFL #'s for each system at each order (log2p+1,NV)
			FLT adis;  //!< Dissipation constant multiplier
		} *gbl;
		
	public:
		void* create_global_structure() {return new global;}
		hp_coupled_bdry(tri_hp &xin, edge_bdry &bin) : hp_edge_bdry(xin,bin) {mytype = "hp_coupled_bdry"; is_master = base.is_frst();}
		hp_coupled_bdry(const hp_coupled_bdry& inbdry, tri_hp &xin, edge_bdry &bin)  : hp_edge_bdry(inbdry,xin,bin), NV(inbdry.NV), is_master(inbdry.is_master), gbl(inbdry.gbl) {
			fine = &inbdry;
			if (!is_master) return;
			/* default initializer for p=1 multigrid (overriddent in init for fine level) */
			vug_frst.resize(base.maxseg+1,NV);
			vdres.resize(1,base.maxseg+1,NV);
		};
		hp_coupled_bdry* create(tri_hp& xin, edge_bdry &bin) const {return new hp_coupled_bdry(*this,dynamic_cast<tri_hp&>(xin),bin);}
		void init(input_map& input,void* gbl_in);
		
		/* FOR COUPLED DYNAMIC BOUNDARIES (COMMENTED ROUTINES NEED TO BE WRITTEN) */
		// void tadvance();
		void rsdl(int stage);
		void maxres();
		void minvrt();  // TEMPO: NOT GENERAL YET
		// void update(int stage);
		void mg_restrict();
#ifdef petsc
		//void non_sparse(Array<int,1> &nnzero);
		//void petsc_jacobian();
		//void petsc_make_1D_rsdl_vector(Array<FLT,1> res);
#endif
};

/* deforming coupled boundary */
class hp_deformable_bdry : public hp_coupled_bdry {
	protected:
		Array<FLT,1> ksprg;
	public:
		hp_deformable_bdry(tri_hp &xin, edge_bdry &bin) : hp_coupled_bdry(xin,bin) {mytype = "hp_deformable_bdry";}
		hp_deformable_bdry(const hp_deformable_bdry& inbdry, tri_hp &xin, edge_bdry &bin) : hp_coupled_bdry(inbdry,xin,bin) {if (is_master) ksprg.resize(base.maxseg);}
		hp_deformable_bdry* create(tri_hp& xin, edge_bdry &bin) const {return new hp_deformable_bdry(*this,dynamic_cast<tri_hp&>(xin),bin);}
		void init(input_map& inmap,void* gbl_in);
		void tadvance();
		void rsdl(int stage);
		void update(int stage);

#ifdef petsc
		void non_sparse(Array<int,1> &nnzero);
		void non_sparse_rcv(Array<int, 1> &nnzero, Array<int, 1> &nnzero_mpi);
		void petsc_make_1D_rsdl_vector(Array<FLT,1> res);
		void petsc_jacobian();
		void petsc_matchjacobian_rcv(int phase);
		void element_jacobian(int indx, Array<FLT,2>& K);
#endif
};
#endif /* defined(__tri_hp__hp_coupled_boundary__) */


