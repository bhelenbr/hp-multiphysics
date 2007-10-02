#include "blocks.h"
#include "block.h"
#include "r_mesh.h"


multigrid_interface* block::getnewlevel(input_map& input) {
    int type;          
    multigrid_interface *temp;  
    
    if (!input.get(idprefix+"_type",type)) input.getwdefault("blocktype",type,1);
    
    switch(type) {
        default: {
            temp = new r_mesh();
            break;
        }
    } 
        
    return(temp);
}

block* blocks::getnewblock(int idnum, input_map& input) {
    return new block(idnum);
}

