/*
 *  initfunc.cpp
 *  planar++
 *
 *  Created by helenbrk on Wed Oct 24 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */
 
#include "tet_hp.h"

#define CD
#define INS

#ifdef CD
#include "cd/tet_hp_cd.h"
#endif

#ifdef INS
#include "ins/tet_hp_ins.h"
#endif

class btype {
    public:
        const static int ntypes = 3;
        enum ids {plain,cd,ins};
        const static char names[ntypes][40];
        static int getid(const char *nin) {
            int i;
            for(i=0;i<ntypes;++i) 
                if (!strcmp(nin,names[i])) return(i);
            return(-1);
        }
};
const char btype::names[ntypes][40] = {"plain","cd","ins"};

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
        case btype::plain: {
            tet_hp *temp = new tet_hp();
            return(temp);
        }
#ifdef CD
        case btype::cd: {
            tet_hp_cd *temp = new tet_hp_cd();
            return(temp);
        }
#endif

#ifdef INS
        case btype::ins: {
            tet_hp_ins *temp = new tet_hp_ins();
            return(temp);
        }
#endif        
        default: {
            std::cerr << "unrecognizable block type: " <<  type << std::endl;
            tet_hp *temp = new tet_hp();
            return(temp);
        }
    } 
        
    return(0);
}

