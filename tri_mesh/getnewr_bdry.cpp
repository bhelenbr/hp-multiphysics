/*
 *  getnewr_bdry.cpp
 *  mblock
 *
 *  Created by helenbrk on Wed Oct 24 2004.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "r_tri_mesh.h"
#include "r_tri_boundary.h"

class r_stype {
    public:
        static const int ntypes = 4;
        enum ids {free=1, fixed, translating, oscillating};
        const static char names[ntypes][40];
        static int getid(const char *nin) {
            for(int i=0;i<ntypes;++i) 
                if (!strcmp(nin,names[i])) return(i+1);
            return(-1);
        }
};

const char r_stype::names[ntypes][40] = {"free", "fixed", "translating", "oscillating"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
r_side_bdry* r_tri_mesh::getnewedgeobject(int bnum, input_map& in_map) {
    std::istringstream data;
    std::map<std::string,std::string>::const_iterator mi;
    r_side_bdry *temp;  
    int type = 2;

    mi = in_map.find(ebdry(bnum)->idprefix + "_r_type");
    if (mi != in_map.end()) {
        type = r_stype::getid((*mi).second.c_str());
        if (type < 0)  {
            *gbl->log << "unknown type:" << (*mi).second << std::endl;
            exit(1);
        }
    }
    else {
        /* SOME DEFAULTS FOR VARIOUS BOUNDARY TYPES */
        if (ebdry(bnum)->mytype == "comm") {
            type = r_stype::free;
        } else if (ebdry(bnum)->mytype == "partition") {
            type = r_stype::free;
        } else if (ebdry(bnum)->mytype == "prdc") {
            type = r_stype::fixed;
            int dir;
            in_map.getwdefault(ebdry(bnum)->idprefix + "_dir",dir,0);
            if (dir == 0)
                in_map[ebdry(bnum)->idprefix+"_r_dir"] = "0 0";
            else 
                in_map[ebdry(bnum)->idprefix+"_r_dir"] = "1 1";
        }
        else {
            type = r_stype::fixed;
        }
        *gbl->log << "# using default r_type of " << r_stype::names[type-1] << " for " << ebdry(bnum)->idprefix << std::endl;
    }

    switch(type) {
        case r_stype::free: {
            temp = new r_side_bdry(*this,*ebdry(bnum));
            break;
        }
        case r_stype::fixed: {
            temp = new r_fixed(*this,*ebdry(bnum)); 
            break;
        }
        case r_stype::translating: {
            temp = new r_translating(*this,*ebdry(bnum)); 
            break;
        }
        case r_stype::oscillating: {
            temp = new r_oscillating(*this,*ebdry(bnum)); 
            break;
        }        
        default: {
            temp = new r_fixed(*this,*ebdry(bnum));
            std::cout << "Don't know this r_side_bdry type\n";
        }
    }
    
    temp->input(in_map);
    
    return(temp);
}