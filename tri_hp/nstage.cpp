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
         ctrl_message = block::begin;
      }
      
      case(2): {
#ifdef CTRL_DEBUG
         *sim::log << "step 2 of tri_hp::update: ctrl_message: " << ctrl_message << " excpt1: " << excpt1 << std::endl;
#endif
         if (ctrl_message != block::advance1) {
            state = block::stop; 
            for(i=0;i<nsbd;++i)
               state &= hp_sbdry(i)->update(ctrl_message);
            
            state &= hp_gbl->mover->update(ctrl_message,nvrtx,vrtx,vrtxbd(1));

            if (state != block::stop) return(state);
            else return(block::advance1);
         }
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
         ctrl_message = block::begin;
      }
     
      case(6): {
      
//            for(i=0;i<nvrtx;++i)
//               printf("nstage: %d %e %e %e %e %e\n", i,hp_gbl->vprcn(i,0),hp_gbl->vprcn(i,2),hp_gbl->res.v(i,0),hp_gbl->res.v(i,1),hp_gbl->res.v(i,2));      

//         for(i=0;i<nvrtx;++i)
//      printf("v: %d %e %e %e\n",i,hp_gbl->res.v(i,0),hp_gbl->res.v(i,1),hp_gbl->res.v(i,2));
//      
//   for(i=0;i<nside*b->sm;++i)
//      printf("s: %d %e %e %e\n",i,gbl->res.s[i][0],gbl->res.s[i][1],gbl->res.s[i][2]);
//
//   for(i=0;i<ntri*b->im;++i)
//      printf("i: %d %e %e %e\n",i,gbl->res.i[i][0],gbl->res.i[i][1],gbl->res.i[i][2]); 
//      
////         std::cout << hp_gbl->res.v(Range(0,nvrtx-1),Range::all());
//         std::cout << hp_gbl->res.s(Range(0,nside-1),Range::all(),Range::all());
//         std::cout << hp_gbl->res.i(Range(0,ntri-1),Range::all(),Range::all());
//         exit(1);
         
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
         ctrl_message = block::advance2;
      }
      
      case(7): {
#ifdef CTRL_DEBUG
         *sim::log << "step 7 of tri_hp::update: ctrl_message: " << ctrl_message << " excpt1: " << excpt1 << std::endl;
#endif
         if (ctrl_message != block::advance1) {
            state = block::stop;
            for(i=0;i<nsbd;++i) {
               state &= hp_sbdry(i)->update(ctrl_message);
            }
            state &= hp_gbl->mover->update(ctrl_message,nvrtx,vrtx,vrtxbd(1));
            
            if (state != block::stop) return(state);
            return(block::advance1);
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
