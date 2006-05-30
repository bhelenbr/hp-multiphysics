#include "tri_hp.h"
#include "hp_boundary.h"

#define NODEBUG

block::ctrl tri_hp::update(block::ctrl ctrl_message) {
   int i,m,k,n,indx,indx1;
   FLT cflalpha;
   block::ctrl state;
 
   if (ctrl_message == block::begin) excpt1 = 0;

   switch (excpt1) {
      case(0): {
#ifdef CTRL_DEBUG
         *sim::log << "step 0 of tri_hp::update: ctrl_message: " << ctrl_message << " excpt1: " << excpt1 << std::endl;
#endif
         /* COUPLED MESH MOVMEMENT */
         if (ctrl_message != block::advance1) {
            if (mmovement == coupled_deformable) {
               state = r_mesh::update(ctrl_message);
               if (state != block::stop) return(state);
            }
            return(block::advance1);
         }
         ++excpt1;
      }
      
      case(1): {
#ifdef CTRL_DEBUG
         *sim::log << "step 1 of tri_hp::update: ctrl_message: " << ctrl_message << " excpt1: " << excpt1 << std::endl;
#endif

         /* STORE INITIAL VALUES FOR NSTAGE EXPLICIT SCHEME */
         hp_gbl->ug0.v(Range(0,nvrtx-1),Range::all()) = ug.v(Range(0,nvrtx-1),Range::all());
         if (basis::tri(log2p).sm) {
            hp_gbl->ug0.s(Range(0,nside-1),Range(0,sm0-1),Range::all()) = ug.s(Range(0,nside-1),Range::all(),Range::all());
            if (basis::tri(log2p).im) {
               hp_gbl->ug0.i(Range(0,ntri-1),Range(0,im0-1),Range::all()) = ug.i(Range(0,ntri-1),Range::all(),Range::all());
            }
         }
         ++excpt1;
      }
      
      case(2): {
#ifdef CTRL_DEBUG
         *sim::log << "step 2 of tri_hp::update: ctrl_message: " << ctrl_message << " excpt1: " << excpt1 << std::endl;
#endif
         for(i=0;i<nsbd;++i)
            hp_sbdry(i)->update(block::begin);
         
         hp_gbl->mover->update(block::begin,nvrtx,vrtx,vrtxbd(1));
         ++excpt1;
         stage = 0;
      }
      
      case(3): {
#ifdef CTRL_DEBUG
         *sim::log << "step 3 of tri_hp::update: ctrl_message: " << ctrl_message << " excpt1: " << excpt1 << std::endl;
#endif
         ctrl_message = block::begin;
         ++excpt1;
      }
      
      case(4): {
#ifdef CTRL_DEBUG
         *sim::log << "step 4 of tri_hp::update: ctrl_message: " << ctrl_message << " excpt1: " << excpt1 << std::endl;
#endif
         if (ctrl_message != block::advance2) {
            state = rsdl(ctrl_message,stage); 
            if (state != block::stop) return(state);
            return(block::advance2);
         }
         ++excpt1;
         ctrl_message = block::begin;
      }
      
      case(5): {
#ifdef CTRL_DEBUG
         *sim::log << "step 5 of tri_hp::update: ctrl_message: " << ctrl_message << " excpt1: " << excpt1 << std::endl;
#endif
         if (ctrl_message != block::advance1) {
            state = minvrt(ctrl_message);
            if (state != block::stop) return(state);
            return(block::advance1);
         }
         ++excpt1;
      }
     
      case(6): {
   
      
   printf("nstage: %d nvrtx: %d log2p: %d\n",stage,nvrtx,log2p);

   for(i=0;i<nvrtx;++i)
      printf("nstage: %d %e %e\n", i,hp_gbl->vprcn(i,0),hp_gbl->vprcn(i,2));

   for(i=0;i<nvrtx;++i)
      printf("v: %d %e %e %e\n",i,hp_gbl->res.v(i,0),hp_gbl->res.v(i,1),hp_gbl->res.v(i,2));
      
   for(i=0;i<nside;++i)
      for(m=0;m<basis::tri(log2p).sm;++m)
         printf("s: %d %d %e %e %e\n",i,m,hp_gbl->res.s(i,m,0),hp_gbl->res.s(i,m,1),hp_gbl->res.s(i,m,2));

   for(i=0;i<ntri;++i)
      for(m=0;m<basis::tri(log2p).im;++m)
         printf("i: %d %d %e %e %e\n",i,m,hp_gbl->res.i(i,m,0),hp_gbl->res.i(i,m,1),hp_gbl->res.i(i,m,2));

   for(i=0;i<nvrtx;++i)
      printf("ug.v: %d %e %e %e\n",i,ug.v(i,0),ug.v(i,1),ug.v(i,2));
      
   for(i=0;i<nside;++i)
      for(m=0;m<basis::tri(log2p).sm;++m)
         printf("ug.s: %d %d %e %e %e\n",i,m,ug.s(i,m,0),ug.s(i,m,1),ug.s(i,m,2));

   for(i=0;i<ntri;++i)
      for(m=0;m<basis::tri(log2p).im;++m)
         printf("ug.i: %d %d %e %e %e\n",i,m,ug.i(i,m,0),ug.i(i,m,1),ug.i(i,m,2));



//   exit(1);
         
         cflalpha = sim::alpha[stage]*hp_gbl->cfl(log2p);
#ifdef CTRL_DEBUG
         *sim::log << "step 6 of tri_hp::update: ctrl_message: " << ctrl_message << " excpt1: " << excpt1 << std::endl;
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
         ++excpt1;
         ctrl_message = block::advance;
      }
      
      case(7): {
#ifdef CTRL_DEBUG
         *sim::log << "step 7 of tri_hp::update: ctrl_message: " << ctrl_message << " excpt1: " << excpt1 << std::endl;
#endif
         if (ctrl_message != block::advance2) {
            state = block::stop;
            for(i=0;i<nsbd;++i) {
               state &= hp_sbdry(i)->update(ctrl_message);
            }
            state &= hp_gbl->mover->update(ctrl_message,nvrtx,vrtx,vrtxbd(1));
            
            if (state != block::stop) return(state);
            return(block::advance2);
         }
         ++excpt1;
      }
      
      case(8): {
#ifdef CTRL_DEBUG
         *sim::log << "step 8 of tri_hp::update: ctrl_message: " << ctrl_message << " excpt1: " << excpt1 << std::endl;
#endif
         stage += 1;
         if (stage < sim::NSTAGE) {
            excpt1 = 3;
            return(block::advance);
         }
         ++excpt1;
      }
   }
   
   return(block::stop);
}
