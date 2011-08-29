/*
 *  hp_mgrid.h
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _tet_hp_cns_h_
#define _tet_hp_cns_h_

#include "../tet_hp.h"
#include <blocks.h>

class tet_hp_cns : public tet_hp {
public:
	/* THINGS SHARED BY ALL tri_hp_ins in same multigrid block */
	struct global : public tet_hp::global {
		/* STABILIZATION */
		Array<FLT,3> tau;
		
		/* PHYSICAL CONSTANTS */
		FLT mu,kcond,R,gamma;
		TinyVector<double,tet_mesh::ND> body;
		
		/* STORAGE FOR CALCULATION OF ENERGY AND AREA */
		TinyVector<FLT,2> eanda, eanda_recv;
				
		Array<FLT,3> vpreconditioner,epreconditioner;
		
		/* SOURCE FUNCTION FOR MMS */
		//init_bdry_cndtn *src;
		
	} *gbl;
	
	
	FLT adis; // DISSIPATION CONSTANT
	
	hp_vrtx_bdry* getnewvrtxobject(int bnum, input_map &bdrydata);
	hp_edge_bdry* getnewedgeobject(int bnum, input_map &bdrydata);
	hp_face_bdry* getnewfaceobject(int bnum, input_map &bdrydata);
	init_bdry_cndtn* getnewibc(std::string suffix, input_map& inmap);
	tet_hp_helper* getnewhelper(input_map& inmap);
	
public:
	void* create_global_structure() {return new global;}
	tet_hp_cns* create() { return new tet_hp_cns(); }
	
	void init(input_map& input, void *gin); 
	void init(const multigrid_interface& in, init_purpose why=duplicate, FLT sizereduce1d=1.0);
	
	void update();
	void minvrt();
	
	//void length(); 
	void setup_preconditioner();
	void element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uhat,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im);
	void calculate_unsteady_sources();
	void calculate_preconditioner_tau_timestep(Array<double,1> pvu, FLT h, Array<FLT,2> &Pinv, Array<FLT,2> &Tau, FLT &timestep);
	//void project_new_variables();
	void switch_variables(Array<double,1> pvu, Array<double,1> &a);
	
	
};
#endif
