/*
 *  hp_mgrid.h
 *  heat++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */
 
#ifndef _tri_hp_buoyancy_h_
#define _tri_hp_buoyancy_h_

#include "../ins/tri_hp_ins.h"
#include <blocks.h>
#include <mathclass.h>

class tri_hp_buoyancy : public tri_hp_ins {
    public:
        struct gbl : public tri_hp_ins::gbl {
        
            /* PHYSICAL CONSTANTS */
            FLT kcond,cp;
            symbolic_function<1> rhovsT;
            
        } *gbl_ptr;
    
    private:
        int excpt;
        
    public:
        void init(input_map& input, gbl *gin); 
        tri_hp_buoyancy* create() { return new tri_hp_buoyancy(); }
        block::ctrl setup_preconditioner(block::ctrl ctrl_message);
        block::ctrl rsdl(block::ctrl ctrl_message, int stage=sim::NSTAGE);
        void calculate_unsteady_sources(bool coarse);
};
#endif
