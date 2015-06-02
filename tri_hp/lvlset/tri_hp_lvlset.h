/*
 *  hp_mgrid.h
 *  lvlset++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _tri_hp_lvlset_h_
#define _tri_hp_lvlset_h_

#include "../ins/tri_hp_ins.h"
#include <blocks.h>

/* DIFFERENT WAYS OF DOING LEVEL-SET */
// #define CONSERVATIVE
#define NONCONSERVATIVE
// #define LOCALIZED_WITH_DISTANCE_FUNCTION

class tri_hp_lvlset : public tri_hp_ins {
	public:
		struct global : public tri_hp_ins::global {
			/* PHYSICAL CONSTANTS */
			FLT sigma, width;
			FLT rho2, mu2;

		} *gbl;
		hp_vrtx_bdry* getnewvrtxobject(int bnum, std::string name);
		hp_edge_bdry* getnewsideobject(int bnum, std::string name);
		tri_hp_helper *getnewhelper(std::string name);

		
		/* Reinitialization Stuff */
		int reinit_iterations;
		class reinit_bc {
			public:
				tri_hp_lvlset& x;
				edge_bdry& base;
				bool active;
				TinyVector<FLT,ND> norm;
				int closer_endpt;
				FLT phi_at_endpt;
				FLT sign_phi;
				void init();
				void set_values();
				void vdirichlet();
				void sdirichlet(int mode);
				reinit_bc(tri_hp_lvlset& xin, edge_bdry &bin) : x(xin), base(bin) {}
		};
		Array<reinit_bc *,1> reinit_bdry;
		void reinitialize();
		void reinit_update();
		void reinit_minvrt();
		void reinit_rsdl(int stage);
		void reinit_setup_preconditioner();
		void reinit_element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uht,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im);
	
		/* Main functions */
		void* create_global_structure() {return new global;}
		tri_hp_lvlset* create() { return new tri_hp_lvlset(); }
		void init(input_map& inmap, void *gin);  
		void init(const multigrid_interface& in, init_purpose why=duplicate, FLT sizereduce1d=1.0);

		// functions to reinitialize phi
		void tadvance() { 
			if (!coarse_level && gbl->substep == 0 && reinit_iterations)
				reinitialize();
			tri_hp_ins::tadvance();
		}
		void calculate_unsteady_sources();
		void setup_preconditioner();
		void element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uhat,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im);

		FLT heavyside(FLT phidw) {
			FLT onemphi2 = (1-phidw*phidw);
			FLT term = 0.5;
			FLT sum = term;
			FLT r = 2.0;
			for (int i=1;i<6;++i) {
				term *= onemphi2*(r-1.)/r;
				sum += term;
				r += 2.0;
			}
			return(0.5 +phidw*sum);
			// return(0.5*(phidw +sin(M_PI*phidw)/M_PI));
		}

		FLT delta(FLT phidw) {
			return(693./512./gbl->width*pow(1.-phidw*phidw,5));
			// return(0.5/gbl->width*(1.+cos(M_PI*phidw)));
		}

		FLT heavyside_if(FLT phidw) {
			if (phidw < -1.0) return(0.0);
			if (phidw >  1.0) return(1.0);
			return(heavyside(phidw));
		}

		void heavyside_and_delta_if(FLT phidw, FLT& h, FLT& d) {
			if (phidw < -1.0) {
				h = 0.0;
				d = 0.0;
			}
			else if (phidw >  1.0) {
				h = 1.0;
				d = 0.0;
			}
			else {
				h = heavyside(phidw);
				d = delta(phidw);
			}
			return;
		}            

};
#endif
