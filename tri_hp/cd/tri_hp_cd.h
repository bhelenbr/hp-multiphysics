/*
 *  tri_hp_cd.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */
 
#ifndef _tri_hp_cd_h_
#define _tri_hp_cd_h_

#include "../tri_hp.h"
#include <blocks.h>

class tri_hp_cd : public tri_hp {
    public:
        /* THINGS SHARED BY ALL tri_hp_ins in same multigrid block */
        struct global : public tri_hp::global {
            /* STABILIZATION */
            Array<FLT,1> tau;
                
            /* PHYSICAL CONSTANTS */
            FLT ax, ay, nu;
            
            /* SOURCE FUNCTION */
            init_bdry_cndtn *src;
            
            /* ADAPTATION LIMITS */
            FLT minlngth, maxlngth;

        } *gbl;
        
        FLT adis; // DISSIPATION CONSTANT
        
        hp_edge_bdry* getnewsideobject(int bnum, input_map &bdrydata);
        init_bdry_cndtn* getnewibc(std::string suffix, input_map& inmap);
        tri_hp_helper* getnewhelper(input_map& inmap);
        
    public:
        void* create_global_structure() {return new global;}
        tri_hp_cd* create() { return new tri_hp_cd(); }
        void init(input_map& input, void *gin); 
        void init(const multigrid_interface& in, init_purpose why=duplicate, FLT sizereduce1d=1.0);
        void length();
        void setup_preconditioner();
        void rsdl(int stage);
};
#endif
