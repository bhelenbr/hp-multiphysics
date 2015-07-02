//
//  hp_coupled_boundary.h
//  tri_hp
//
//  Created by Brian Helenbrook on 12/27/14.
//
//

#ifndef _hp_coupled_boundary_h_
#define _hp_coupled_boundary_h_

#include "hp_boundary.h"
#include <iostream>

//#define MPDEBUG
//#define DEBUG
//#define DETAILED_MINV
//#define ROTATE_RESIDUALS
//#define ONE_SIDED

/* Non-deforming coupled boundary with variables on boundary */
class hp_coupled_bdry : public hp_edge_bdry {
    public:
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
		void init(input_map& inmap,void* gbl_in);
		
		/* FOR COUPLED DYNAMIC BOUNDARIES (COMMENTED ROUTINES NEED TO BE WRITTEN) */
		// void tadvance();
		void rsdl(int stage);
		void maxres();
		void minvrt();  //FIXME: NOT GENERAL YET
		// void update(int stage);
		void mg_restrict();
		void mg_source();
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
		void element_jacobian(int indx, Array<FLT,2>& K);
		void update(int stage);

#ifdef petsc
		void non_sparse(Array<int,1> &nnzero);
		void non_sparse_rcv(Array<int, 1> &nnzero, Array<int, 1> &nnzero_mpi);
		void petsc_make_1D_rsdl_vector(Array<FLT,1> res);
		void petsc_jacobian();
		void petsc_matchjacobian_rcv(int phase);
#ifdef ROTATE_RESIDUALS
		void petsc_premultiply_jacobian();
#endif
#endif
};

class hp_deformable_fixed_pnt : public hp_vrtx_bdry {
	/* INTERSECTING BOUNDARY CONTAINING END POINT MUST HAVE GEOMETRY NOT BE DEFINED SOLELY BY MESH */
protected:
	hp_deformable_bdry *surface;
	int surfbdry;
	
public:
	hp_deformable_fixed_pnt(tri_hp &xin, vrtx_bdry &bin) : hp_vrtx_bdry(xin,bin) {mytype = "hp_deformable_fixed_pnt";}
	hp_deformable_fixed_pnt(const hp_deformable_fixed_pnt& inbdry, tri_hp &xin, vrtx_bdry &bin) : hp_vrtx_bdry(inbdry,xin,bin), surfbdry(inbdry.surfbdry) {
		if (surfbdry > -1) {
			if (!(surface = dynamic_cast<hp_deformable_bdry *>(x.hp_ebdry(base.ebdry(surfbdry))))) {
				*x.gbl->log << "something's wrong can't find surface boundary" << std::endl;
				sim::abort(__LINE__,__FILE__,x.gbl->log);
			}
		}
		else {
			surface = 0;  // Not a surface endpoint just a mesh fixed point
		}
	}
	hp_deformable_fixed_pnt* create(tri_hp& xin, vrtx_bdry &bin) const {return new hp_deformable_fixed_pnt(*this,dynamic_cast<tri_hp&>(xin),bin);}
	
	void init(input_map& inmap,void* gbl_in);
	void vdirichlet();
#ifdef petsc
	void petsc_jacobian_dirichlet();
#endif
};



class hp_deformable_free_pnt : public hp_deformable_fixed_pnt {
	/* For a surface point sliding on a vertical or horizontal surface */
	/* For periodic wall have tri_mesh vertex type be comm */
protected:
	FLT position;
	enum {vertical, horizontal, curved} wall_type;
	
public:
	hp_deformable_free_pnt(tri_hp &xin, vrtx_bdry &bin) : hp_deformable_fixed_pnt(xin,bin) {mytype = "hp_deformable_free_pnt";}
	hp_deformable_free_pnt(const hp_deformable_free_pnt& inbdry, tri_hp &xin, vrtx_bdry &bin) : hp_deformable_fixed_pnt(inbdry,xin,bin), position(inbdry.position) {}
	hp_deformable_free_pnt* create(tri_hp& xin, vrtx_bdry &bin) const {return new hp_deformable_free_pnt(*this,dynamic_cast<tri_hp&>(xin),bin);}
	
	void init(input_map& inmap,void* gbl_in);
	void element_rsdl(Array<FLT,1> lf);
	void rsdl(int stage);
	void vdirichlet() {hp_vrtx_bdry::vdirichlet();}
#ifdef petsc
	void petsc_jacobian();
	void petsc_jacobian_dirichlet() {hp_vrtx_bdry::petsc_jacobian_dirichlet();}
#endif
};


class translating_surface : public hp_deformable_bdry {
protected:
	TinyVector<FLT,tri_mesh::ND> vel;
	
public:
	translating_surface(tri_hp &xin, edge_bdry &bin) : hp_deformable_bdry(xin,bin) {mytype = "translating_surface";}
	translating_surface(const translating_surface& inbdry, tri_hp &xin, edge_bdry &bin)  : hp_deformable_bdry(inbdry,xin,bin), vel(inbdry.vel) {};
	translating_surface* create(tri_hp& xin, edge_bdry &bin) const {return new translating_surface(*this,xin,bin);}
	
	void init(input_map& inmap,void* gbl_in);
	void element_rsdl(int sind, Array<TinyVector<FLT,MXTM>,1> lf);
	void setup_preconditioner();
};

class hp_deformable_follower_pnt : public hp_deformable_free_pnt {
public:
	hp_deformable_follower_pnt(tri_hp &xin, vrtx_bdry &bin) : hp_deformable_free_pnt(xin,bin) {mytype = "hp_deformable_follower_pnt";}
	hp_deformable_follower_pnt(const hp_deformable_follower_pnt& inbdry, tri_hp &xin, vrtx_bdry &bin) : hp_deformable_free_pnt(inbdry,xin,bin) {}
	hp_deformable_follower_pnt* create(tri_hp& xin, vrtx_bdry &bin) const {return new hp_deformable_follower_pnt(*this,dynamic_cast<tri_hp&>(xin),bin);}
	
	/** This routine zeros residual so it can just follow */
	void element_rsdl(Array<FLT,1> lf) { hp_vrtx_bdry::element_rsdl(lf);}
	void rsdl(int stage);
#ifdef petsc
	void petsc_jacobian(); 
#endif
};

#endif /* defined(__hp_coupled_boundary__) */


