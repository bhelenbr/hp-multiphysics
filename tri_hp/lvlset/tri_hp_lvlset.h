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
		bool reinit_flag;
		struct global : public tri_hp_ins::global {
			/* PHYSICAL CONSTANTS */
			FLT sigma, width;
			FLT rho2, mu2;

		} *gbl;
		hp_vrtx_bdry* getnewvrtxobject(int bnum, input_map &bdrydata);
		hp_edge_bdry* getnewsideobject(int bnum, input_map &bdrydata);
	public:
		// functions to reinitialize phi
		void tadvance() { 
			if (!coarse_level && gbl->substep == 0 && reinit_flag)
				reinitialize();
			tri_hp_ins::tadvance();
		}
	private:
		struct normstuff{
			bool use;
			FLT nx, ny, xloc, yloc, basephi;
		};

		void reinitialize();
		void reinit(normstuff norm0list[], normstuff norm1list[]);
		void minvrt_reinit(normstuff norm0list[], normstuff norm1list[]);
		void rsdl_reinit(int stage);
		void setup_preconditioner_reinit();
		void element_rsdl_reinit(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uht,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im);
		bool minvrt_reinit_phival(int edgenum, int segnum, FLT nx, FLT ny, FLT xloc, FLT yloc, FLT basephi, int onezero, FLT& phi);
		int minvrt_reinit_phiside(int edgenum);
		void minvrt_reinit_direction(int edgenum, int side, FLT& nx, FLT& ny, FLT& xloc, FLT &yloc, FLT &basephi);
	public:
		void setup_preconditioner();
		void element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uhat,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im);

		void* create_global_structure() {return new global;}
		tri_hp_lvlset* create() { return new tri_hp_lvlset(); }

		void init(input_map& input, void *gin);  
		void init(const multigrid_interface& in, init_purpose why=duplicate, FLT sizereduce1d=1.0);

		void calculate_unsteady_sources();

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
