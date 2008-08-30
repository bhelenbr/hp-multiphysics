/*
 *  initfunc.cpp
 *  planar++
 *
 *  Created by helenbrk on Wed Oct 24 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */
 
#include <blocks.h>
#include <block.h>
#include <r_tri_mesh.h>

#include "cd/tri_hp_cd.h"
#include "ins/tri_hp_ins.h"
#include "planestrain/tri_hp_ps.h"
#include "swirl/tri_hp_swirl.h"
#include "buoyancy/tri_hp_buoyancy.h"
#include "pod/pod.h"
#include "swe/tri_hp_swe.h"
#include "lvlset/tri_hp_lvlset.h"
#include "hp_boundary.h"

// #define JUST_INS

class btype {
    public:
        const static int ntypes = 12;
        enum ids {r_tri_mesh,cd,ins,ps,swirl,buoyancy,pod_ins_gen,pod_cd_gen,pod_ins_sim,pod_cd_sim,swe,lvlset};
        const static char names[ntypes][40];
        static int getid(const char *nin) {
            int i;
            for(i=0;i<ntypes;++i) 
                if (!strcmp(nin,names[i])) return(i);
            return(-1);
        }
};
const char btype::names[ntypes][40] = {"r_tri_mesh","cd","ins","ps","swirl","buoyancy",
    "pod_ins_gen","pod_cd_gen","pod_ins_sim","pod_cd_sim","swe","lvlset"};

multigrid_interface* block::getnewlevel(input_map& input) {
    std::string keyword,val,ibcname,srcname;
    std::istringstream data;
    int type;          
    
    /* FIND BLOCK TYPE */
    if (input.get(idprefix+"_type",val)) {
        type = btype::getid(val.c_str());
    }
    else {
        if (!input.get("blocktype",val)) {
            std::cerr << "couldn't find block type" << std::endl;
            exit(1);
        }
        type = btype::getid(val.c_str());
    }
            
    switch(type) {
        case btype::r_tri_mesh: {
            r_tri_mesh *temp = new r_tri_mesh();
            return(temp);
        }
#ifndef JUST_INS
        case btype::cd: {
            tri_hp_cd *temp = new tri_hp_cd();
            return(temp);
        }
#endif
        
        case btype::ins: {
            tri_hp_ins *temp = new tri_hp_ins();
            return(temp);
        }

#ifndef JUST_INS
        case btype::ps: {
            tri_hp_ps *temp = new tri_hp_ps();
            return(temp);
        }
		
		case btype::swirl: {
            tri_hp_swirl *temp = new tri_hp_swirl();
            return(temp);
        }
        
        case btype::buoyancy: {
            tri_hp_buoyancy *temp = new tri_hp_buoyancy();
            return(temp);
        }
        
        case btype::pod_ins_gen: {
            pod_generate<tri_hp_ins> *temp = new pod_generate<tri_hp_ins>();
            return(temp);
        }
        
        case btype::pod_cd_gen: {
            pod_generate<tri_hp_cd> *temp = new pod_generate<tri_hp_cd>();
            return(temp);
        }
        
        case btype::pod_ins_sim: {
            pod_simulate<tri_hp_ins> *temp = new pod_simulate<tri_hp_ins>();
            return(temp);
        }
        
        case btype::pod_cd_sim: {
            pod_simulate<tri_hp_cd> *temp = new pod_simulate<tri_hp_cd>();
            return(temp);
        }
        
        case btype::swe: {
            tri_hp_swe *temp = new tri_hp_swe();
            return(temp);
        }

        case btype::lvlset: {
            tri_hp_lvlset *temp = new tri_hp_lvlset();
            return(temp);
        }
#endif
        
        default: {
            std::cerr << "unrecognizable block type: " <<  type << std::endl;
            r_tri_mesh *temp = new r_tri_mesh();
            return(temp);
        }
    } 
        
    return(0);
}

