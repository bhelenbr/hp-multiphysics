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


/* Non-deforming coupled boundary with variables on boundary */
class hp_coupled_bdry : public hp_edge_bdry {
public:
	
	/* non shared data */
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
	Array<FLT,1> ksprg;  //!< For tangential deforming boundary equation
	
	/* shared data */
	/* Is this the right way to do this?  Can just make these a member of class and then have other objects reference these */
	/* Only problem with that is that constants become a little weird.  Would have to have pointers to doubles instead of an actual double :-( */
	/* Should distinguish those that can be changed i.e. physical parameters and those that can't i.e. NV */
	/* Could make a structure of constants, then just keep a reference sort of like gbl now except don't have to keep inheriting globals?? */
	/* That's a little weird too, what would subclasses name their structures of globals? */
	struct global {
		bool is_loop;  //!< true if boundary is a loop
		//  Different ways of treating continuity between sides of boundary in iteration
		//  traditional way is that vertex variable equations are duplicated and side modes are treated with explicit constraint
		//  one_sided - all c0 variables (vertex & side modes) are just forced to be continuous with an explicit constraint
		//  symmetric - equations are duplicated for both sides of the boundary
		//  precondition - premultiplys jacobian and residual to increase diagonal dominance before inverting
		//                 this can only be done with symmetric or one sided approach (not traditional)
		bool one_sided, symmetric, precondition;

		/* SOLUTION STORAGE ON FIRST ENTRY TO NSTAGE */
		Array<FLT,2> vug0; //!< vertex solution on entry to multigrid (pnts,NV)
		Array<FLT,3> sug0; //!< side solution on entry to multigrid (segs,mode,NV)
		
		/* RESIDUALS */
		Array<FLT,2> vres; //!< vertex residuals (pnts,NV)
		Array<FLT,3> sres; //!< side residuals (segs,mode,NV)
		Array<FLT,2> vres0; //!< vertex multigrid driving terms (pnts,NV)
		Array<FLT,3> sres0; //!< side multigrid driving terms (segs,mode,NV)
		
		/* PRECONDITIONER */
		bool field_is_coupled; //!< determines whether jacobian must be large enough to accomodate off-edge values
		Array<FLT,3> vdt; //!< Time advancement for vertices (pnts,x.NV+NV,x.NV+NV)
		Array<int,2> vpiv; //!< Pivot matrix for vdt inversion
		Array<FLT,3> sdt; //!< Time advancement for sides (segs,x.NV+NV,x.NV+NV)
		Array<int,2> spiv; //!< Pivot matrix for sdt inversion
		Array<FLT,4> sdt2; //!< Time advancement for sides (segs,sm,x.NV+NV,x.NV+NV)
		Array<int,3> spiv2; //!< Pivot matrix for sdt2 inversion
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
	void delete_global_structure() { if(shared_owner) delete gbl;}
	hp_coupled_bdry(tri_hp &xin, edge_bdry &bin) : hp_edge_bdry(xin,bin) {mytype = "hp_coupled_bdry"; is_master = base.is_frst();}
	hp_coupled_bdry(const hp_coupled_bdry& inbdry, tri_hp &xin, edge_bdry &bin)  : hp_edge_bdry(inbdry,xin,bin), NV(inbdry.NV), is_master(inbdry.is_master), gbl(inbdry.gbl) {
		fine = &inbdry;
		if (!gbl->symmetric && !is_master) return;
		/* default initializer for p=1 multigrid (overriddent in init for fine level) */
		ksprg.resize(base.maxseg);
		vug_frst.resize(base.maxseg+1,NV);
		vdres.resize(1,base.maxseg+1,NV);
	};
	hp_coupled_bdry* create(tri_hp& xin, edge_bdry &bin) const {return new hp_coupled_bdry(*this,dynamic_cast<tri_hp&>(xin),bin);}
	void init(input_map& inmap,void* gbl_in);
	
	/* FOR COUPLED DYNAMIC BOUNDARIES (COMMENTED ROUTINES NEED TO BE WRITTEN) */
	void tadvance();
	void rsdl(int stage);
	void maxres();
	void update(int stage);
	void minvrt();  //FIXME: NOT GENERAL YET
	void mg_restrict();
	void mg_source();
	void smatchsolution_snd(FLT *sdata, int bgn, int end, int stride);
	int smatchsolution_rcv(FLT *sdata, int bgn, int end, int stride);
	
	void element_jacobian(int indx, Array<FLT,2>& K);
#ifdef petsc
	void non_sparse(Array<int,1> &nnzero);
	void non_sparse_snd(Array<int, 1> &nnzero, Array<int, 1> &nnzero_mpi);
	int non_sparse_rcv(int phase, Array<int, 1> &nnzero, Array<int, 1> &nnzero_mpi);
	void petsc_jacobian();
	void petsc_matchjacobian_snd();
	int petsc_matchjacobian_rcv(int phase);
	void petsc_make_1D_rsdl_vector(Array<FLT,1> res);
	void petsc_premultiply_jacobian();
#endif
};

class hp_deformable_fixed_pnt : public multi_physics_pnt {
	/* INTERSECTING BOUNDARY CONTAINING END POINT MUST HAVE GEOMETRY NOT BE DEFINED SOLELY BY MESH */
protected:
	hp_coupled_bdry *surf;
	int surfbdry;
	
public:
	hp_deformable_fixed_pnt(tri_hp &xin, vrtx_bdry &bin) : multi_physics_pnt(xin,bin) {mytype = "hp_deformable_fixed_pnt";}
	hp_deformable_fixed_pnt(const hp_deformable_fixed_pnt& inbdry, tri_hp &xin, vrtx_bdry &bin) : multi_physics_pnt(inbdry,xin,bin), surfbdry(inbdry.surfbdry) {
		if (surfbdry > -1) {
			if (!(surf = dynamic_cast<hp_coupled_bdry *>(x.hp_ebdry(base.ebdry(surfbdry))))) {
				*x.gbl->log << base.idprefix << " something's wrong can't find surface boundary " << surfbdry << ' ' << base.ebdry(surfbdry) << ' ' << x.ebdry(base.ebdry(surfbdry))->idprefix << std::endl;
				x.tri_mesh::output("darn");
				sim::abort(__LINE__,__FILE__,x.gbl->log);
			}
		}
		else {
			surf = 0;  // Not a surface endpoint just a mesh fixed point
		}
	}
	hp_deformable_fixed_pnt* create(tri_hp& xin, vrtx_bdry &bin) const {return new hp_deformable_fixed_pnt(*this,dynamic_cast<tri_hp&>(xin),bin);}
	
	void init(input_map& inmap,void* gbl_in);
	void vdirichlet();
#ifdef petsc
	void petsc_jacobian_dirichlet();
	void setup_preconditioner();
#endif
};



class hp_deformable_free_pnt : public hp_deformable_fixed_pnt {
	/* For a surface point sliding on a vertical or horizontal surface */
	/* For periodic wall have tri_mesh vertex type be comm */
protected:
	FLT position;
	enum {vertical, horizontal, curved} wall_type;
	TinyVector<FLT,tri_mesh::ND> wall_normal;
	
public:
	hp_deformable_free_pnt(tri_hp &xin, vrtx_bdry &bin) : hp_deformable_fixed_pnt(xin,bin) {mytype = "hp_deformable_free_pnt";}
	hp_deformable_free_pnt(const hp_deformable_free_pnt& inbdry, tri_hp &xin, vrtx_bdry &bin) : hp_deformable_fixed_pnt(inbdry,xin,bin), position(inbdry.position) {}
	hp_deformable_free_pnt* create(tri_hp& xin, vrtx_bdry &bin) const {return new hp_deformable_free_pnt(*this,dynamic_cast<tri_hp&>(xin),bin);}
	
	void init(input_map& inmap,void* gbl_in);
	void element_rsdl(Array<FLT,1> lf);
	void rsdl(int stage);
	void vdirichlet() {hp_vrtx_bdry::vdirichlet();}
	void mvpttobdry(TinyVector<FLT,tri_mesh::ND> &pt);
#ifdef petsc
	void petsc_jacobian();
	void petsc_jacobian_dirichlet() {hp_vrtx_bdry::petsc_jacobian_dirichlet();}
	void setup_preconditioner() {hp_vrtx_bdry::setup_preconditioner();}
#endif
};


class translating_surface : public hp_coupled_bdry {
protected:
	TinyVector<FLT,tri_mesh::ND> vel;
	
public:
	translating_surface(tri_hp &xin, edge_bdry &bin) : hp_coupled_bdry(xin,bin) {mytype = "translating_surface";}
	translating_surface(const translating_surface& inbdry, tri_hp &xin, edge_bdry &bin)  : hp_coupled_bdry(inbdry,xin,bin), vel(inbdry.vel) {};
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
	void setup_preconditioner() {hp_deformable_fixed_pnt::setup_preconditioner();}
#endif
};

#endif /* defined(__hp_coupled_boundary__) */


