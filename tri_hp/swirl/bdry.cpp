#include "bdry_swirl.h"
#include <myblas.h>

/*************************************************/
/* SET DIRICHLET BOUNDARY VALUES & FLUXES ********/
/* (THINGS THAT ARE INDEPENDENT OF THE SOLUTION) */
/* BUT NOT INDEPENDENT OF TIME/MESH POSITION    */
/*************************************************/

using namespace bdry_swirl;

block::ctrl symmetry::tadvance(bool coarse, block::ctrl ctrl_message) {
   int j,m,v0,sind,indx;
   TinyVector<FLT,mesh::ND> pt;
   block::ctrl state;
   
   if (ctrl_message == block::begin) excpt1 = 0;
   
   if (excpt1 == 0) {
      if ((state = hp_side_bdry::tadvance(coarse,ctrl_message)) != block::stop) return(state);
      ++excpt1;
   
      /* UPDATE BOUNDARY CONDITION VALUES */
      for(j=0;j<base.nel;++j) {
         sind = base.el(j);
         v0 = x.sd(sind).vrtx(0);
         x.ug.v(v0,0) = 0.0;
         x.ug.v(v0,2) = 0.0;
      }
      v0 = x.sd(sind).vrtx(1);
      x.ug.v(v0,0) = 0.0;
      x.ug.v(v0,2) = 0.0;

      /*******************/   
      /* SET SIDE VALUES */
      /*******************/
      for(j=0;j<base.nel;++j) {
         sind = base.el(j);
         indx = sind*x.sm0;
         for(m=0;m<basis::tri(x.log2p).sm;++m) {
            x.ug.s(sind,m,0) = 0.0;
            x.ug.s(sind,m,2) = 0.0;
         }
      }
   }
   
   return(block::stop);
}