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

class tri_hp_komega : public tri_hp_ins {
public:
    struct hp_komega_global {
        /* Physical Constants */
        FLT linf,uinf;
        
        /* Model constants */
        FLT kinf, omginf, epslnk;
        FLT Clim;  // only used for Wilcox2006
        int version;
        int kmom_on, sust_on;
        
        /* SOURCE FUNCTION FOR MMS */
         init_bdry_cndtn *src;
    };
    shared_ptr<hp_komega_global> hp_komega_gbl;


    
public:
    tri_hp_komega* create() { return new tri_hp_komega(); }
    
    void init(input_map& inmap,shared_ptr<block_global> gin);
    void init(const multigrid_interface& in, init_purpose why=duplicate, FLT sizereduce1d=1.0);
    
    void error_estimator();
    void element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uhat,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im);
};
#endif
