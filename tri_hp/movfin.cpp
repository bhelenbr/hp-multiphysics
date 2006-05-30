#include "tri_hp.h"
#include "hp_boundary.h"

block::ctrl tri_hp::mg_getcchng(block::ctrl ctrl_message,Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, tri_hp *cmesh) {
   int i,j,ind,tind;
   block::ctrl state;
   int stop = 1;
         
   if (ctrl_message == block::begin) {
      excpt = 0;
      if(cmesh == this) {
         excpt = -1;
         ++log2p;
         if (log2p == log2pmax) coarse = false;
         return(block::stop);
      }
   }
      
   switch (excpt) {
      case(0): {
         if (ctrl_message != block::advance1) {
            if (mmovement == coupled_deformable) {
               state = r_mesh::mg_getcchng(ctrl_message,fv_to_ct,cv_to_ft,cmesh);
               if (state != block::stop) return(state);
            }
            return(block::advance1);
         }
         ++excpt;
         ctrl_message = block::begin;
      }

      case(1): {
         if (ctrl_message != block::advance1) {
            /* TRANFER COUPLED BOUNDARY UPDATES */
            state = block::stop;
            for(i=0;i<nsbd;++i)
               state &= hp_sbdry(i)->mg_getcchng(ctrl_message, fv_to_ct, cv_to_ft, cmesh, i);
            if (state != block::stop) return(state);
            return(block::advance1);
         }
         ++excpt;
      }
         
      case(2): {
         /* DETERMINE CORRECTIONS ON COARSE MESH   */   
         cmesh->vug_frst(Range(0,nvrtx-1),Range::all()) -= cmesh->ug.v(Range(0,nvrtx-1),Range::all());
   
         /* LOOP THROUGH FINE VERTICES   */
         /* TO DETERMINE CHANGE IN SOLUTION */   
         for(i=0;i<nvrtx;++i) {
            tind = fv_to_ct(i).tri;

            hp_gbl->res.v(i,Range::all()) = 0.0;
            
            for(j=0;j<3;++j) {
               ind = cmesh->td(tind).vrtx(j);
               hp_gbl->res.v(i,Range::all()) -= fv_to_ct(i).wt(j)*cmesh->vug_frst(ind,Range::all());
            }
         }
         ++excpt;
         mp_phase = -1;
         ctrl_message = block::stay;
      }
         
      case(3): {
         /* SEND COMMUNICATION PACKETS  */ 
         if (ctrl_message == block::stay) {
            ++mp_phase;

            switch(mp_phase%3) {
               case(0):
                  for(i=0;i<nsbd;++i)
                     sbdry(i)->vloadbuff(boundary::partitions,(FLT *) hp_gbl->res.v.data(),0,NV-1,NV);
                     
                  for(i=0;i<nsbd;++i) 
                     sbdry(i)->comm_prepare(boundary::partitions,mp_phase/3,boundary::symmetric);
                  
                  return(block::stay);
               case(1):
                  for(i=0;i<nsbd;++i) 
                     sbdry(i)->comm_transmit(boundary::partitions,mp_phase/3,boundary::symmetric);
                  return(block::stay);
               case(2):
                  stop = 1;
                  for(i=0;i<nsbd;++i) {
                     stop &= sbdry(i)->comm_wait(boundary::partitions,mp_phase/3,boundary::symmetric);
                     sbdry(i)->vfinalrcv(boundary::partitions,mp_phase/3,boundary::symmetric,boundary::average,(FLT *) hp_gbl->res.v.data(),0,NV-1,NV);
                  }
                  return(static_cast<block::ctrl>(stop));
            }
         }
         ++excpt;
      }
      case(4): {
         /* ADD CORRECTION */
         ug.v(Range(0,nvrtx-1),Range::all()) += hp_gbl->res.v(Range(0,nvrtx-1),Range::all());     
         ++excpt;          
      }
   }
   return(block::stop);
}
   
   
   
   
      
      
