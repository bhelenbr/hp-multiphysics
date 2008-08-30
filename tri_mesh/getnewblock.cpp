#include "blocks.h"
#include "block.h"
#include "r_tri_mesh.h"


multigrid_interface* block::getnewlevel(input_map& input) {
    int type;          
    multigrid_interface *temp;  
    
    if (!input.get(idprefix+"_type",type)) input.getwdefault("type",type,1);
    
    switch(type) {
        default: {
            temp = new r_tri_mesh();
            break;
        }
    } 
        
    return(temp);
}

