#include "blocks.h"
#include "block.h"
#include "mgblock.h"
#include "r_mesh.h"

block* blocks::getnewblock(int type) {
   return(new mgrid<r_mesh>);
}
