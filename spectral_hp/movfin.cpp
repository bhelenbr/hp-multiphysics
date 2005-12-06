#include "tri_hp.h"
#include "hp_boundary.h"

block::ctrl tri_hp::mg_getcchng(int excpt,Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, tri_hp *cmesh) {
   int i,j,ind,tind;
   block::ctrl state;
   int stop = 1;
   static int last_r_mesh;
   
   if(cmesh == this) {
      ++log2p;
      if (log2p == log2pmax) coarse = false;
      return(block::stop);
   }
   if (excpt == 0) last_r_mesh = 0;
   
   if (mmovement == coupled_deformable) {
      state = r_mesh::mg_getcchng(excpt,fv_to_ct,cv_to_ft,cmesh);
      if (state != block::stop) {
         last_r_mesh = excpt +1;
         return(state);
      }
   }
   excpt -= last_r_mesh;

      
   switch (excpt) {
      case(0): {
         mp_phase = -1;
   
         /* TRANFER COUPLED BOUNDARY UPDATES */
         for(i=0;i<nsbd;++i)
            hp_sbdry(i)->mg_getcchng(excpt, fv_to_ct, cv_to_ft, cmesh->hp_sbdry(i));

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
         return(block::advance);
      }
         
      case(1): {
         /* SEND COMMUNICATION PACKETS  */ 
         ++mp_phase;

         switch(mp_phase%3) {
            case(0):
               for(i=0;i<nsbd;++i)
                  if (sbdry(i)->group()&0x1) sbdry(i)->vloadbuff((FLT *) hp_gbl->res.v.data(),0,NV-1,NV);
                  
               for(i=0;i<nsbd;++i) 
                  if (sbdry(i)->group()&0x1) sbdry(i)->comm_prepare(mp_phase/3);
               
               return(block::stay);
            case(1):
               for(i=0;i<nsbd;++i) 
                  if (sbdry(i)->group()&0x1) sbdry(i)->comm_transmit(mp_phase/3);
               return(block::stay);
            case(2):
               stop = 1;
               for(i=0;i<nsbd;++i) {
                  if (sbdry(i)->group()&0x1) {
                     stop &= sbdry(i)->comm_wait(mp_phase/3);
                     sbdry(i)->vfinalrcv(mp_phase/3,(FLT *) hp_gbl->res.v.data(),0,NV-1,NV);
                  }
               }
               return(static_cast<block::ctrl>(stop));
         }
      }
      case(2): {
         /* ADD CORRECTION */
         ug.v(Range(0,nvrtx-1),Range::all()) += hp_gbl->res.v(Range(0,nvrtx-1),Range::all());               
         return(block::stop);
      }
   }
   *sim::log << "Flow control error\n" << std::endl;
   return(block::stop);
}
   
   
   
   
      
      
