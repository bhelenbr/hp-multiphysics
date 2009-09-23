/*
 *  tri_hp_cd.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _tet_hp_cd_h_
#define _tet_hp_cd_h_

#include "../tet_hp.h"
#include <blocks.h>

class tet_hp_cd : public tet_hp {
	public:
		/* THINGS SHARED BY ALL tri_hp_ins in same multigrid block */
		struct global : public tet_hp::global {
			/* STABILIZATION */
			Array<FLT,1> tau;

			/* PHYSICAL CONSTANTS */
			FLT ax,ay,az, nu;
			
			/* SOURCE FUNCTION */
			init_bdry_cndtn *src;

		} *gbl;
		
		FLT adis; // DISSIPATION CONSTANT
		
		hp_face_bdry* getnewfaceobject(int bnum, input_map &bdrydata);
		init_bdry_cndtn *getnewibc(std::string suffix, input_map& inmap);
		
	public:
		void* create_global_structure() {return new global;}
		tet_hp_cd* create() { return new tet_hp_cd(); }
		void init(input_map& input, void *gin); 
		void init(const multigrid_interface& fine, init_purpose why=duplicate, FLT sizereduce1d=1.0);
//      void length();
		void setup_preconditioner();
//		void rsdl(int stage);
		void element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uht,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im);

};
#endif
