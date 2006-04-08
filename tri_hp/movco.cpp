#include "tri_hp.h"
#include "hp_boundary.h"

block::ctrl tri_hp::mg_getfres(block::ctrl ctrl_message, Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, tri_hp *fmesh) {
   int i,j,k,m,n,tind,v0,indx,indx1;
   block::ctrl state;
   
   isfrst = true;
   coarse = true;
      
   if (ctrl_message == block::begin) excpt = 0;
         
   switch(excpt) {
      case(0): {
         if(p0 > 1) {
            --log2p;
            
            /* TRANSFER IS ON FINE MESH */
            hp_gbl->res0.v(Range(0,nvrtx-1),Range::all()) = hp_gbl->res.v(Range(0,nvrtx-1),Range::all());

            if (basis::tri(log2p).p > 1) {
               hp_gbl->res0.s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all()) = hp_gbl->res.s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all());
               
               if (basis::tri(log2p).p > 2) {
              
                  for(tind=0;tind<ntri;++tind) {
                     indx = 0;
                     indx1 = 0;
                     for(m=1;m<basis::tri(log2p).sm;++m) {
                        for(k=0;k<basis::tri(log2p).sm-m;++k) {
                           hp_gbl->res0.i(tind,indx,Range::all()) = hp_gbl->res.i(tind,indx1,Range::all());
                           ++indx;
                           ++indx1;
                        }
                        indx1 += basis::tri(log2p).p;
                     }
                  }
               }
            }
         }
         else {

            hp_gbl->res0.v(Range(0,nvrtx-1),Range::all()) = 0.0;
                  
            /* LOOP THROUGH FINE VERTICES TO CALCULATE RESIDUAL  */
            for(i=0;i<fmesh->nvrtx;++i) {
               tind = fv_to_ct(i).tri;
               for(j=0;j<3;++j) {
                  v0 = td(tind).vrtx(j);
                  for(n=0;n<NV;++n)
                     hp_gbl->res0.v(v0,n) += fv_to_ct(i).wt(j)*hp_gbl->res.v(i,n);
               }
            }
            
            /* LOOP THROUGH COARSE VERTICES   */
            /* TO CALCULATE VUG ON COARSE MESH */
            for(i=0;i<nvrtx;++i) {
               tind = cv_to_ft(i).tri;

               ug.v(i,Range::all()) = 0.0;

               for(j=0;j<3;++j) {
                  ug.v(i,Range::all()) += cv_to_ft(i).wt(j)*fmesh->ug.v(fmesh->td(tind).vrtx(j),Range::all());
               }

               vug_frst(i,Range::all()) = ug.v(i,Range::all());
            }
         }
         ++excpt;
      }
      
      case(1): {
         if (ctrl_message != block::advance1) {
            if (mmovement == coupled_deformable && p0 == 1) { 
               block::ctrl state = r_mesh::mg_getfres(ctrl_message, fv_to_ct, cv_to_ft, fmesh);
               if (state != block::stop) return(state);
            }
            return(block::advance1);
         }
         ++excpt;
         ctrl_message = block::begin;
      }

      case(2): {
         if (ctrl_message != block::advance1) {
            state = block::stop;
            for(i=0;i<nsbd;++i)
               state &= hp_sbdry(i)->mg_getfres(ctrl_message, fv_to_ct, cv_to_ft, fmesh, i);
            if (state != block::stop) return(state);
            return(block::advance1);
         }
         ++excpt;
      }
      
   }

   return(block::stop);
}

