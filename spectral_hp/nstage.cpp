#include "tri_hp.h"
#include "hp_boundary.h"

block::ctrl tri_hp::update(int excpt) {
   int i,m,k,n,indx,indx1;
   FLT cflalpha;
   int unshiftedexcpt;
   static int lastresidual,lastminvrt,laststage,addtostage,stage;
   block::ctrl state;
   
   if (excpt == 0) {

      /* STORE INITIAL VALUES FOR NSTAGE EXPLICIT SCHEME */
      hp_gbl->ug0.v(Range(0,nvrtx),Range::all()) = ug.v(Range(0,nvrtx),Range::all());

      if (basis::tri(log2p).p > 1) {
         hp_gbl->ug0.s(Range(0,nside),Range(0,sm0),Range::all()) = ug.s(Range(0,nside),Range::all(),Range::all());
         
         if (basis::tri(log2p).p > 2) {
            hp_gbl->ug0.i(Range(0,ntri),Range(0,im0),Range::all()) = ug.i(Range(0,ntri),Range::all(),Range::all());
         }   
      }

      /* COUPLED BOUNDARY MOVEMENT EQUATIONS */   
      for(i=0;i<nsbd;++i)
         hp_sbdry(i)->update(0);
      
      lastresidual = excpt +1;
      laststage = 0;
      stage = 0;
   }
   
   unshiftedexcpt = excpt;
   excpt -= laststage;
   
   /* CALCULATE RESIDUAL */
   if (excpt < lastresidual) {
      state = rsdl(excpt);
      if (state != block::stop) {
         lastresidual = excpt +1;
         return(state);
      }
      lastminvrt = excpt +1;
   }
   
   /* INVERT MASS MATRIX */
   if (excpt < lastminvrt) {
      state = minvrt(excpt-lastresidual);
      if(state != block::stop) {
         lastminvrt = excpt +1;
         return(state);
      }
   }

   if (excpt < lastminvrt+1) {
      cflalpha = 2.0; // cfl[log2p]*alpha[stage];
      
      ug.v(Range(0,nvrtx),Range::all()) = hp_gbl->ug0.v(Range(0,nvrtx),Range::all()) -cflalpha*hp_gbl->res.v(Range(0,nvrtx),Range::all());

      if (basis::tri(log2p).sm > 0) {
         ug.s(Range(0,nside),Range(0,basis::tri(log2p).sm),Range::all()) = hp_gbl->ug0.s(Range(0,nside),Range(0,basis::tri(log2p).sm),Range::all()) -cflalpha*hp_gbl->res.s(Range(0,nside),Range(0,basis::tri(log2p).sm),Range::all());

         if (basis::tri(log2p).im > 0) {

            for(i=0;i<ntri;++i) {
               indx = 0;
               indx1 = 0;
               for(m=1;m<basis::tri(log2p).sm;++m) {
                  for(k=0;k<basis::tri(log2p).sm-m;++k) {
                     for(n=0;n<NV;++n) {
                        ug.i(i,indx1,n) =  hp_gbl->ug0.i(i,indx1,n) -cflalpha*hp_gbl->res.i(i,indx,n);
                     }
                     ++indx; ++indx1;
                  }
                  indx1 += sm0 -basis::tri(log2p).sm;
               }
            }
         }
      }
      
      /* COUPLED BOUNDARY MOVEMENT EQUATIONS */   
      for(i=0;i<nsbd;++i) {
         hp_sbdry(i)->update(stage);
      }
      
      addtostage = 1;
      return(block::advance);
   }
   
   /* UPDATE STAGE & REPEAT */
   laststage = unshiftedexcpt;
   stage += addtostage;
   addtostage = 0;
      
   if (stage < sim::NSTAGE) {
      return(block::advance);
   }
   
   return(block::stop);
}
