/*
 *  hp_mgrid.h
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */
 
#ifndef _tet_hp_ins_h_
#define _tet_hp_ins_h_

#include "../tet_hp.h"
#include <blocks.h>

class tet_hp_ins : public tet_hp {
    public:
        /* THINGS SHARED BY ALL tri_hp_ins in same multigrid block */
        struct global : public tet_hp::global {
            /* STABILIZATION */
            Array<FLT,2> tau;
                
            /* PHYSICAL CONSTANTS */
            FLT rho, mu;
            Array<FLT,1> D;
            
            /* STORAGE FOR CALCULATION OF ENERGY AND AREA */
            TinyVector<FLT,2> eanda, eanda_recv;

        } *gbl;

     
        FLT adis; // DISSIPATION CONSTANT
        
       // hp_vrtx_bdry* getnewvrtxobject(int bnum, input_map &bdrydata);
       // hp_edge_bdry* getnewedgeobject(int bnum, input_map &bdrydata);
        hp_face_bdry* getnewfaceobject(int bnum, input_map &bdrydata);
        init_bdry_cndtn* getnewibc(std::string suffix, input_map& inmap);
//        tet_hp_helper* getnewhelper(input_map& inmap);

    public:
        void* create_global_structure() {return new global;}
        tet_hp_ins* create() { return new tet_hp_ins(); }

        void init(input_map& input, void *gin); 
        void init(const multigrid_interface& in, init_purpose why=duplicate, FLT sizereduce1d=1.0);

       // void length(); dont need yet for adaption
       void setup_preconditioner();
       void rsdl(int stage);
       // void calculate_unsteady_sources();
                
};
#endif
