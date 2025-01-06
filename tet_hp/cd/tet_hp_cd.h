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
		struct hp_cd_global {
			/* STABILIZATION */
			Array<FLT,1> tau;

			/* PHYSICAL CONSTANTS */
			FLT rhocv,ax,ay,az,kcond;
			
			/* SOURCE FUNCTION */
			init_bdry_cndtn *src;
        };
        shared_ptr<hp_cd_global> hp_cd_gbl;
		
		FLT adis; // DISSIPATION CONSTANT

		hp_vrtx_bdry* getnewvrtxobject(int bnum, std::string name);
		hp_edge_bdry* getnewedgeobject(int bnum, std::string name);
		hp_face_bdry* getnewfaceobject(int bnum, std::string name);

	
		init_bdry_cndtn* getnewibc(std::string ibcname);
		
	public:
		tet_hp_cd* create() { return new tet_hp_cd(); }
		void init(input_map& inmap, shared_ptr<block_global> gin); 
		void init(const multigrid_interface& fine, init_purpose why=duplicate, FLT sizereduce1d=1.0);
		void calculate_unsteady_sources();
//    void length();
//		void minvrt();
		int setup_preconditioner();
		void element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uht,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im);

};
#endif
