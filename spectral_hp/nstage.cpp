#include "tri_hp.h"
#include "hp_boundary.h"

#define NODEBUG

block::ctrl tri_hp::update(int excpt) {
   int i,m,k,n,indx,indx1;
   FLT cflalpha;
   int unshiftedexcpt;
   static int lastresidual,lastminvrt,laststage,addtostage,stage,last_r_mesh;
   block::ctrl state;
   
   if (excpt == 0) last_r_mesh = 0;
   
   /* COUPLED MESH MOVMEMENT */
   if (mmovement == coupled_deformable) {
      state = r_mesh::update(excpt);
      if (state != block::stop) {
         last_r_mesh = excpt+1;
         return(state);
      }
   }
   excpt -= last_r_mesh;
         
   if (excpt == 0) {
   
#ifdef DEBUG 
      *sim::log << "loading ug0 last_r_mesh: " << last_r_mesh << std::endl;
#endif

      /* STORE INITIAL VALUES FOR NSTAGE EXPLICIT SCHEME */
      hp_gbl->ug0.v(Range(0,nvrtx-1),Range::all()) = ug.v(Range(0,nvrtx-1),Range::all());
      if (basis::tri(log2p).sm) {
         hp_gbl->ug0.s(Range(0,nside-1),Range(0,sm0-1),Range::all()) = ug.s(Range(0,nside-1),Range::all(),Range::all());
         if (basis::tri(log2p).im) {
            hp_gbl->ug0.i(Range(0,ntri-1),Range(0,im0-1),Range::all()) = ug.i(Range(0,ntri-1),Range::all(),Range::all());
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
#ifdef DEBUG
      *sim::log << "performing residual step " << excpt << ' ' << stage << std::endl;
#endif
      state = rsdl(excpt,stage);
      if (state != block::stop) {
         lastresidual = excpt +1;
         return(state);
      }
      lastminvrt = excpt +1;
   }
      
   /* INVERT MASS MATRIX */
   if (excpt < lastminvrt) {
#ifdef DEBUG
      *sim::log << "performing minvrt step " << excpt -lastresidual +1 << std::endl;
#endif
      state = minvrt(excpt-lastresidual +1);
      if(state != block::stop) {
         lastminvrt = excpt +2;
         return(state);
      }
   }

   /* UPDATE SOLUTION */
   if (excpt == lastminvrt-1) {
      cflalpha = sim::alpha[stage]; // cfl[log2p]*alpha[stage];

#ifdef DEBUG 
      *sim::log << "updating solution: " << stage << std::endl;
#endif
      
      ug.v(Range(0,nvrtx-1),Range::all()) = hp_gbl->ug0.v(Range(0,nvrtx-1),Range::all()) -cflalpha*hp_gbl->res.v(Range(0,nvrtx-1),Range::all());

      if (basis::tri(log2p).sm > 0) {
         ug.s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all()) = hp_gbl->ug0.s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all()) -cflalpha*hp_gbl->res.s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all());

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
   laststage = unshiftedexcpt+1;
   stage += addtostage;
   addtostage = 0;
      
   if (stage < sim::NSTAGE) {
      return(block::advance);
   }
   
   return(block::stop);
}
