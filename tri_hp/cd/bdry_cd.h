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


// #define MELT1 /* Must be uncommented in two places */

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
			dirichlet(tri_hp_cd &xin, edge_bdry &bin) : generic(xin,bin) {
				mytype = "dirichlet";
				type[0] = essential;
				essential_indices.push_back(0);
			}
			dirichlet(const dirichlet& inbdry, tri_hp_cd &xin, edge_bdry &bin) : generic(inbdry,xin,bin) {}
			dirichlet* create(tri_hp& xin, edge_bdry &bin) const {return new dirichlet(*this,dynamic_cast<tri_hp_cd&>(xin),bin);}
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

	class melt : public generic {
		protected:
			tri_hp_cd &x;
#ifdef MELT1
			Array<TinyVector<FLT,tri_mesh::ND>,1> vug_frst;
			Array<TinyVector<FLT,tri_mesh::ND>,2> vdres; //!< Driving term for multigrid (log2p, pnts)
			Array<TinyVector<FLT,tri_mesh::ND>,3> sdres; //!< Driving term for multigrid (log2p, side, order)
			const melt *fine, *coarse;
			
		public:
			struct global {    
				bool is_loop;
				
				/* SOLUTION STORAGE ON FIRST ENTRY TO NSTAGE */
				Array<TinyVector<FLT,tri_mesh::ND>,1> vug0;
				Array<TinyVector<FLT,tri_mesh::ND>,2> sug0;
				
				/* RESIDUALS */
				Array<TinyVector<FLT,tri_mesh::ND>,1> vres;
				Array<TinyVector<FLT,tri_mesh::ND>,2> sres;
				Array<TinyVector<FLT,tri_mesh::ND>,1> vres0;
				Array<TinyVector<FLT,tri_mesh::ND>,2> sres0;
				
				/* PRECONDITIONER */
				Array<FLT,1> meshc;
				
				FLT fadd;
				FLT adis;
			} *gbl; 
			void* create_global_structure() {return new global;}
#endif
			
		public:
			melt(tri_hp_cd &xin, edge_bdry &bin) : generic(xin,bin), x(xin) {mytype = "melt";}
			melt(const melt& inbdry, tri_hp_cd &xin, edge_bdry &bin) : generic(inbdry,xin,bin), x(xin) {
#ifdef MELT1
				gbl = inbdry.gbl;
				vug_frst.resize(base.maxseg+1);
				vdres.resize(1,base.maxseg+1);
				fine = &inbdry;
#endif
			}
			melt* create(tri_hp& xin, edge_bdry &bin) const {return new melt(*this,dynamic_cast<tri_hp_cd&>(xin),bin);}
			
			void init(input_map& inmap,void* gbl_in);
			/* FOR COUPLED DYNAMIC BOUNDARIES */
			void tadvance();
			void rsdl(int stage);
#ifdef MELT1
			void element_rsdl(int sind, Array<TinyVector<FLT,MXTM>,1> lf);
			void setup_preconditioner();
			void mg_restrict(); 
			void element_jacobian(int indx, Array<FLT,2>& K);
#endif
			void vdirichlet();
			void sdirichlet(int mode);
			void update(int stage);
#ifdef petsc
			void petsc_make_1D_rsdl_vector(Array<FLT,1> res);
			/* This is wickedly complicated (doh! suffering for confusion early on about manifold b.c.'s) */
			/* element_rsdl/rsdl adds to residual for temperature, mass conservation */
			/* and then makes tangential equation residual */
			/* rsdl also applies b.c. to tangent residual vbdry->rsdl() */
			/* tangential residual goes in vres(0)/sres(0) */
			/* vdirichlet moves heat equation residual to vres(1)/sres(1) after it has been made */
			/* Then applies dirichlet b.c. to temperature */
			/* in petsc case vdirichlet rotates residual to be aligned with x,y and puts it in r_mesh rsdl */
			/* vbdry-vdirchlet zeros residuals for mesh movement in x,y directions in both places (r_gbl and vres */
			/* petsc_make_1D_rsdl_vector does the rotating instead of sdirichlet for the side modes (called from petsc_make_1D_rsdl_vector) */
			
			
			/* For jacobian */
			/* petsc_jacobian calculates how residuals change in normal way */
			/* petsc_matchjacobian_snd()  Sends T, tangent, normal */
			/* petsc_matchjacobian_rcv()	Receives T */
			/* petsc_matchjacobian_dirichlet() -> moves heat equation jacobian to normal row */
			/* then rotates normal and tangential row to x,y directions */
			
			void petsc_jacobian();
			void petsc_matchjacobian_snd();
			void petsc_matchjacobian_rcv(int phase);
			void non_sparse(Array<int,1> &nnzero);
			void non_sparse_snd(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi);
			void non_sparse_rcv(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi);
	#endif
	};
}
#endif
