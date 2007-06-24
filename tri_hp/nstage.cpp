#include "tri_hp.h"
#include "hp_boundary.h"

// #define DEBUG
// #define CTRL_DEBUG

block::ctrl tri_hp::update(block::ctrl ctrl_message) {
    int i,m,k,n,indx,indx1;
    FLT cflalpha;
    block::ctrl state;
 
    if (ctrl_message == block::begin) excpt1 = 0;

    switch (excpt1) {
        case(0): {
#ifdef CTRL_DEBUG
            *sim::log << idprefix << " step 0 of tri_hp::update: ctrl_message: " << ctrl_message << " excpt1: " << excpt1 << std::endl;
#endif
            /* COUPLED MESH MOVMEMENT */
            if (ctrl_message != block::advance1) {
                if (mmovement == coupled_deformable  && log2p == 0) {
                    state = r_mesh::update(ctrl_message);
                    if (state != block::stop) return(state);
                }
                return(block::advance1);
            }
            ++excpt1;
        }
        
        case(1): {
#ifdef CTRL_DEBUG
            *sim::log << idprefix << " step 1 of tri_hp::update: ctrl_message: " << ctrl_message << " excpt1: " << excpt1 << std::endl;
#endif

            /* STORE INITIAL VALUES FOR NSTAGE EXPLICIT SCHEME */
            gbl_ptr->ug0.v(Range(0,nvrtx-1),Range::all()) = ug.v(Range(0,nvrtx-1),Range::all());
            if (basis::tri(log2p).sm) {
                gbl_ptr->ug0.s(Range(0,nside-1),Range(0,sm0-1),Range::all()) = ug.s(Range(0,nside-1),Range::all(),Range::all());
                if (basis::tri(log2p).im) {
                    gbl_ptr->ug0.i(Range(0,ntri-1),Range(0,im0-1),Range::all()) = ug.i(Range(0,ntri-1),Range::all(),Range::all());
                }
            }
            ++excpt1;
        }
        
        case(2): {
#ifdef CTRL_DEBUG
            *sim::log << idprefix << " step 2 of tri_hp::update: ctrl_message: " << ctrl_message << " excpt1: " << excpt1 << std::endl;
#endif
            mover->update(block::begin);
            
            for(i=0;i<nsbd;++i)
                hp_sbdry(i)->update(block::begin);
            
            ++excpt1;
            stage = 0;
        }
        
        case(3): {
#ifdef CTRL_DEBUG
            *sim::log << idprefix << " step 3 of tri_hp::update: ctrl_message: " << ctrl_message << " excpt1: " << excpt1 << std::endl;
#endif
            ctrl_message = block::begin;
            ++excpt1;
        }
        
        case(4): {
#ifdef CTRL_DEBUG
            *sim::log << idprefix << " step 4 of tri_hp::update: ctrl_message: " << ctrl_message << " excpt1: " << excpt1 << std::endl;
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
            *sim::log << idprefix << " step 5 of tri_hp::update: ctrl_message: " << ctrl_message << " excpt1: " << excpt1 << std::endl;
#endif
            if (ctrl_message != block::advance1) {
                state = minvrt(ctrl_message);
                if (state != block::stop) return(state);
                return(block::advance1);
            }
            ++excpt1;
        }
      
        case(6): {
 #ifdef DEBUG      
    printf("%s nstage: %d nvrtx: %d log2p: %d\n",idprefix.c_str(),stage,nvrtx,log2p);

    for(i=0;i<nvrtx;++i) {
        printf("%s nstage: %d ",idprefix.c_str(),i);
        for(n=0;n<NV;++n) {
            if (fabs(gbl_ptr->vprcn(i,n)) > 1.0e-9) printf("%8.5e ",gbl_ptr->vprcn(i,n));
            else printf("%8.5e ",0.0);
        }
        printf("\n");
    }

    for(i=0;i<nvrtx;++i) {
        printf("%s v: %d ",idprefix.c_str(),i);
        for(n=0;n<NV;++n) {
            if (fabs(gbl_ptr->res.v(i,n)) > 1.0e-9) printf("%8.5e ",gbl_ptr->res.v(i,n));
            else printf("%8.5e ",0.0);
        }
        printf("\n");
    }

    for(i=0;i<nside;++i) {
        for(m=0;m<basis::tri(log2p).sm;++m) {
            printf("%s s: %d ",idprefix.c_str(),i);
            for(n=0;n<NV;++n) {
                if (fabs(gbl_ptr->res.s(i,m,n)) > 1.0e-9) printf("%8.5e ",gbl_ptr->res.s(i,m,n));
                else printf("%8.5e ",0.0);
            }
            printf("\n");
        }
    }
    
    
    for(i=0;i<ntri;++i) {
        for(m=0;m<basis::tri(log2p).im;++m) {
            printf("%s i: %d ",idprefix.c_str(),i);
            for(n=0;n<NV;++n) {
                if (fabs(gbl_ptr->res.i(i,m,n)) > 1.0e-9) printf("%8.5e ",gbl_ptr->res.i(i,m,n));
                else printf("%8.5e ",0.0);
            }
            printf("\n");
        }
    }
    
    for(i=0;i<nvrtx;++i) {
        printf("%s ug.v: %d ",idprefix.c_str(),i);
        for(n=0;n<NV;++n) {
            if (fabs(ug.v(i,n)) > 1.0e-9) printf("%8.5e ",ug.v(i,n));
            else printf("%8.5e ",0.0);
        }
        printf("\n");
    }

    for(i=0;i<nside;++i) {
        for(m=0;m<basis::tri(log2p).sm;++m) {
            printf("%s ug.s: %d ",idprefix.c_str(),i);
            for(n=0;n<NV;++n) {
                if (fabs(ug.s(i,m,n)) > 1.0e-9) printf("%8.5e ",ug.s(i,m,n));
                else printf("%8.5e ",0.0);
            }
            printf("\n");
        }
    }
    
    
    for(i=0;i<ntri;++i) {
        for(m=0;m<basis::tri(log2p).im;++m) {
            printf("%s ug.i: %d ",idprefix.c_str(),i);
            for(n=0;n<NV;++n) {
                if (fabs(ug.i(i,m,n)) > 1.0e-9) printf("%8.5e ",ug.i(i,m,n));
                else printf("%8.5e ",0.0);
            }
            printf("\n");
        }
    }
#endif
            
            cflalpha = sim::alpha[stage]*gbl_ptr->cfl(log2p);
#ifdef CTRL_DEBUG
            *sim::log << idprefix << " step 6 of tri_hp::update: ctrl_message: " << ctrl_message << " excpt1: " << excpt1 << std::endl;
#endif
        
            ug.v(Range(0,nvrtx-1),Range::all()) = gbl_ptr->ug0.v(Range(0,nvrtx-1),Range::all()) -cflalpha*gbl_ptr->res.v(Range(0,nvrtx-1),Range::all());

            if (basis::tri(log2p).sm > 0) {
                ug.s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all()) = gbl_ptr->ug0.s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all()) -cflalpha*gbl_ptr->res.s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all());

                if (basis::tri(log2p).im > 0) {

                    for(i=0;i<ntri;++i) {
                        indx = 0;
                        indx1 = 0;
                        for(m=1;m<basis::tri(log2p).sm;++m) {
                            for(k=0;k<basis::tri(log2p).sm-m;++k) {
                                for(n=0;n<NV;++n) {
                                    ug.i(i,indx1,n) =  gbl_ptr->ug0.i(i,indx1,n) -cflalpha*gbl_ptr->res.i(i,indx,n);
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

            if (ctrl_message != block::advance2) {
                state = mover->update(ctrl_message);
                if (state != block::stop) return(state);
                return(block::advance2);
            }
            ++excpt1;
            ctrl_message = block::advance;
        }
        case(8): {
#ifdef CTRL_DEBUG
            *sim::log << idprefix << " step 8 of tri_hp::update: ctrl_message: " << ctrl_message << " excpt1: " << excpt1 << std::endl;
#endif
            if (ctrl_message != block::advance2) {
                state = block::stop;
                for(i=0;i<nsbd;++i) {
                    state &= hp_sbdry(i)->update(ctrl_message);
                }
                
                if (state != block::stop) return(state);
                return(block::advance2);
            }
#ifdef DEBUG
            exit(1);
#endif
            ++excpt1;
        }
        
        case(9): {
#ifdef CTRL_DEBUG
            *sim::log << idprefix << " step 9 of tri_hp::update: ctrl_message: " << ctrl_message << " excpt1: " << excpt1 << std::endl;
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
