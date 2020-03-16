/*
 *  hp_mgrid.h
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 29 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _tri_hp_komega_h_
#define _tri_hp_komega_h_

#include "../ins/tri_hp_ins.h"
#include <blocks.h>
#include <symbolic_function.h>

#define MMS

class tri_hp_komega : public tri_hp_ins {
public:
    struct global : public tri_hp_ins::global {
        /* Physical Constants */
        FLT linf,uinf;
        
        /* Model constants */
        FLT c_mu, kinf, omginf, epslnk;
        
        /* SOURCE FUNCTION FOR MMS */
        #ifdef MMS
            init_bdry_cndtn *src;
        #endif
        
    } *gbl;


    
public:
    void* create_global_structure() {return new global;}
    void delete_global_structure() {tri_hp::delete_global_structure(); delete gbl;}
    tri_hp_komega* create() { return new tri_hp_komega(); }
    
    void init(input_map& inmap, void *gin);
    void init(const multigrid_interface& in, init_purpose why=duplicate, FLT sizereduce1d=1.0);
    
    void error_estimator();
    void element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uhat,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im);
};
#endif
