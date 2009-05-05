/*
 *  hp_mgrid.h
 *  swe++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */
 
#ifndef _tri_hp_swe_h_
#define _tri_hp_swe_h_

#include "../ins/tri_hp_ins.h"
#include <blocks.h>

class tri_hp_swe : public tri_hp_ins {
    public:
        struct global : public tri_hp_ins::global {
        
            /* PHYSICAL CONSTANTS */
            FLT f0, cbeta, cd, ptest;
            
            /* BATHYMETRY DATA */
            init_bdry_cndtn *bathy;
            
        } *gbl;

        init_bdry_cndtn* getnewbathy(std::string suffix, input_map& inmap);
		init_bdry_cndtn* getnewibc(std::string suffix, input_map& inmap);
		hp_edge_bdry* getnewsideobject(int bnum, input_map &bdrydata);

    public:
        void* create_global_structure() {return new global;}
        tri_hp_swe* create() { return new tri_hp_swe(); }
        
        void init(input_map& input, void *gin);  
        void init(const multigrid_interface& in, init_purpose why=duplicate, FLT sizereduce1d=1.0);

        void setup_preconditioner();
        void rsdl(int stage);
        void calculate_unsteady_sources() {tri_hp::calculate_unsteady_sources();}
};
#endif
