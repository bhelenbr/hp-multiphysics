//
//  shock.h
//  tri_hp
//
//  Created by Brian Helenbrook on 12/30/14.
//
//

#ifndef _shock_h_
#define _shock_h_

#include <iostream>
#include "../hp_coupled_boundary.h"
#include "tri_hp_cns.h"
//#include <tri_boundary.h>

namespace bdry_cns {

//#define MPDEBUG
//#define DEBUG
//#define WAY1

class shock : public hp_coupled_bdry{
	protected:
		tri_hp_cns &x;
		Array<FLT,2> u_opp_v;
		Array<FLT,3> u_opp_e;
	
	public:
		struct global : public hp_coupled_bdry::global {
			/* FLUID PROPERTIES */
			// FLT sigma,rho2,mu2,p_ext;
		} *gbl;
		
		void* create_global_structure() {return new global;}
		void delete_global_structure() { if(shared_owner) delete gbl;}
		shock(tri_hp_cns &xin, edge_bdry &bin) : hp_coupled_bdry(xin,bin), x(xin) {mytype = "shock";}
		shock(const shock& inbdry, tri_hp_cns &xin, edge_bdry &bin)  : hp_coupled_bdry(inbdry,xin,bin), x(xin) {
			gbl = inbdry.gbl;
		};
		shock* create(tri_hp& xin, edge_bdry &bin) const {return new shock(*this,dynamic_cast<tri_hp_cns&>(xin),bin);}
		
		void init(input_map& inmap,void* gbl_in);
        void send_opposite();
		void rsdl(int stage);
		void element_rsdl(int sind, Array<TinyVector<FLT,MXTM>,1> lf);
        void element_jacobian_opp(int indx, Array<FLT,2>& K);
		int shock_mach(FLT &Mu, FLT cu, FLT cd, FLT vdiff);
#ifdef petsc
		void petsc_jacobian();
#ifdef WAY1
        void non_sparse_snd(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi);
#else
        void element_jacobian(int indx, Array<FLT,2>& K);
#endif
        int non_sparse_rcv(int phase, Array<int,1> &nnzero, Array<int,1> &nnzero_mpi);
#endif
};
}
#endif /* defined(_shock_h_) */
