#include "tri_hp.h"
#include "hp_boundary.h"
#include <myblas.h>

#define NO_CTRL_DEBUG

/************************************************/
/**********        INVERT MASS MATRIX     **********/
/************************************************/
block::ctrl tri_hp::minvrt(block::ctrl ctrl_message) {
    int i,j,k,m,n,tind,sind,v0,indx,indx1,indx2,sgn,msgn;
    TinyVector<int,3> sign,side;
    Array<FLT,2> tinv(NV,NV);
    Array<FLT,1> temp(NV);
        
    if (ctrl_message == block::begin) excpt = 0;
    else excpt += ctrl_message;

    switch(excpt) {
        case(0): {
#ifdef CTRL_DEBUG
            *sim::log << idprefix << " step 0 of tri_hp::minvrt: ctrl_message: " << ctrl_message << " excpt: " << excpt << " mode " << mode << std::endl;
#endif
            /* LOOP THROUGH SIDES */
            if (basis::tri(log2p).sm > 0) {
                indx = 0;
                for(sind = 0; sind<nside;++sind) {
                    /* SUBTRACT SIDE CONTRIBUTIONS TO VERTICES */            
                    for (k=0; k <basis::tri(log2p).sm; ++k) {
                        for (i=0; i<2; ++i) {
                            v0 = sd(sind).vrtx(i);
                            for(n=0;n<NV;++n)
                                gbl_ptr->res.v(v0,n) -= basis::tri(log2p).sfmv(i,k)*gbl_ptr->res.s(sind,k,n);
                        }
                        ++indx;
                    }
                }
                    
                if (basis::tri(log2p).im > 0) {
                    /* SUBTRACT INTERIORS */
                    indx = 0;
                    for(tind = 0; tind<ntri;++tind) {
                        indx2 = 3;
                        for (i=0; i<3; ++i) {
                            v0 = td(tind).vrtx(i);
                            for (k=0;k<basis::tri(log2p).im;++k)
                                for(n=0;n<NV;++n)
                                    gbl_ptr->res.v(v0,n) -= basis::tri(log2p).ifmb(i,k)*gbl_ptr->res.i(tind,k,n);

                            sind = td(tind).side(i);
                            sgn = td(tind).sign(i);
                            msgn = 1;
                            for (j=0;j<basis::tri(log2p).sm;++j) {
                                for (k=0;k<basis::tri(log2p).im;++k)
                                    for(n=0;n<NV;++n)
                                        gbl_ptr->res.s(sind,j,n) -= msgn*basis::tri(log2p).ifmb(indx2,k)*gbl_ptr->res.i(tind,k,n);
                                msgn *= sgn;
                                ++indx2;
                            }
                        }
                      indx += basis::tri(log2p).im;
                    }
                }
            }
            
            if (gbl_ptr->diagonal_preconditioner) {
                gbl_ptr->res.v(Range(0,nvrtx-1),Range::all()) *= gbl_ptr->vprcn(Range(0,nvrtx-1),Range::all());
            } else {
                /* ASSUMES LOWER TRIANGULAR FOR NOW */
                for(i=0;i<nvrtx;++i) {
                    DGETUS(&gbl_ptr->vprcn_ut(i,0,0), NV, NV, &gbl_ptr->res.v(i,0));
                }
            }
            /* PREPARE MESSAGE PASSING */
            mp_phase = -1;
            return(block::advance);
        }
        
        case(1): {
#ifdef CTRL_DEBUG
            *sim::log << idprefix << "step 1 of tri_hp::minvrt: ctrl_message: " << ctrl_message << " excpt: " << excpt << " mode " << mode << std::endl;
#endif
            ++mp_phase;
            switch(mp_phase%3) {
                case(0):
                    vc0load(mp_phase/3,gbl_ptr->res.v.data());
                    return(block::stay);
                case(1):
                    vmsgpass(boundary::all_phased,mp_phase/3,boundary::symmetric);
                    return(block::stay);
                case(2):
                    return(static_cast<block::ctrl>(vc0wait_rcv(mp_phase/3,gbl_ptr->res.v.data())));
            }
        }
        
        case(2): {
#ifdef CTRL_DEBUG
            *sim::log << idprefix << "step 2 of tri_hp::minvrt: ctrl_message: " << ctrl_message << " excpt: " << excpt << " mode " << mode << std::endl;
#endif
            /* APPLY VERTEX DIRICHLET B.C.'S */
            for(i=0;i<nsbd;++i)
                hp_sbdry(i)->vdirichlet();
                
            for(i=0;i<nvbd;++i)
                hp_vbdry(i)->vdirichlet2d();
                    
            if(basis::tri(log2p).sm == 0) {
                excpt = 7;
                return(block::stop);
            }
            
            /* REMOVE VERTEX CONTRIBUTION FROM SIDE MODES */
            /* SOLVE FOR SIDE MODES */
            /* PART 1 REMOVE VERTEX CONTRIBUTIONS */
            for(tind=0;tind<ntri;++tind) {
            
            if (gbl_ptr->diagonal_preconditioner) {
                    for(i=0;i<3;++i) {
                        v0 = td(tind).vrtx(i);
                        for(n=0;n<NV;++n)
                            uht(n)(i) = gbl_ptr->res.v(v0,n)*gbl_ptr->tprcn(tind,n);
                    }

                    for(i=0;i<3;++i) {
                        sind = td(tind).side(i);
                        sgn  = td(tind).sign(i);
                        for(j=0;j<3;++j) {
                            indx1 = (i+j)%3;
                            msgn = 1;
                            for(k=0;k<basis::tri(log2p).sm;++k) {
                                for(n=0;n<NV;++n)
                                    gbl_ptr->res.s(sind,k,n) -= msgn*basis::tri(log2p).vfms(j,k)*uht(n)(indx1);
                                msgn *= sgn;
                            }
                        }
                    }
                }
                else {
                    /* THIS IS TO USE A MATRIX PRECONDITIONER */
                    for(i=0;i<3;++i) {
                        v0 = td(tind).vrtx(i);
                        for(n=0;n<NV;++n)
                            uht(n)(i) = gbl_ptr->res.v(v0,n);
                    }

                    for(i=0;i<3;++i) {
                        indx = td(tind).side(i);
                        sgn  = td(tind).sign(i);
                        for(j=0;j<3;++j) {
                            indx1 = (i+j)%3;
                            msgn = 1;
                            for(k=0;k<basis::tri(log2p).sm;++k) {
                                for(n=0;n<NV;++n) {
                                    for(m=0;m<NV;++m) {
                                        gbl_ptr->res.s(indx,k,n) -= msgn*basis::tri(log2p).vfms(j,k)*gbl_ptr->tprcn_ut(tind,n,m)*uht(m)(indx1);
                                    }
                                }
                                msgn *= sgn;
                            }
                        }
                    }
                }
            }
            mode = 0;
            return(block::advance);
        }
        
        case(3): {
#ifdef CTRL_DEBUG
            *sim::log << idprefix << "step 3 of tri_hp::minvrt: ctrl_message: " << ctrl_message << " excpt: " << excpt << " mode " << mode << std::endl;
#endif
            /* SOLVE FOR SIDE MODE */
            if (gbl_ptr->diagonal_preconditioner) {
                gbl_ptr->res.s(Range(0,nside-1),mode,Range::all()) *= gbl_ptr->sprcn(Range(0,nside-1),Range::all())*basis::tri(log2p).sdiag(mode);
            }
            else {
                for(sind = 0; sind < nside; ++sind) {
                    DGETUS(&gbl_ptr->sprcn_ut(sind,0,0), NV, NV, &gbl_ptr->res.s(sind,mode,0));
                    gbl_ptr->res.s(sind,mode,Range::all()) /= basis::tri(log2p).sdiag(mode);
                }
            }
            mp_phase = -1;
            return(block::advance);
        }
        
        case(4): {
#ifdef CTRL_DEBUG
            *sim::log << idprefix << "step 4 of tri_hp::minvrt: ctrl_message: " << ctrl_message << " excpt: " << excpt << " mode " << mode << std::endl;
#endif
            ++mp_phase;
            switch(mp_phase%3) {
                case(0):
                    sc0load(gbl_ptr->res.s.data(),mode,mode,gbl_ptr->res.s.extent(secondDim));
                    return(block::stay);
                case(1):
                    smsgpass(boundary::all,0,boundary::symmetric);
                    return(block::stay);
                case(2):
                    sc0wait_rcv(gbl_ptr->res.s.data(),mode,mode,gbl_ptr->res.s.extent(secondDim));
                    return(block::advance);
            }
        }
        
        case(5): {
#ifdef CTRL_DEBUG
            *sim::log << idprefix << "step 5 of tri_hp::minvrt: ctrl_message: " << ctrl_message << " excpt: " << excpt << " mode " << mode << std::endl;
#endif
            /* APPLY DIRCHLET B.C.S TO MODE */
            for(i=0;i<nsbd;++i)
                hp_sbdry(i)->sdirichlet(mode);
                
            if (mode == basis::tri(log2p).sm-1) {
                return(block::advance);
            }
          
            /* REMOVE MODE FROM HIGHER MODES */
            for(tind=0;tind<ntri;++tind) {

            if (gbl_ptr->diagonal_preconditioner) {
                    for(i=0;i<3;++i) {
                        side(i) = td(tind).side(i);
                        sign(i) = td(tind).sign(i);
                        sgn      = (mode % 2 ? sign(i) : 1);
                        for(n=0;n<NV;++n)
                            uht(n)(i) = sgn*gbl_ptr->res.s(side(i),mode,n)*gbl_ptr->tprcn(tind,n);
                    }
                    
                    /* REMOVE MODES J,K FROM MODE I,M */
                    for(i=0;i<3;++i) {
                        msgn = ( (mode +1) % 2 ? sign(i) : 1);
                        for(m=mode+1;m<basis::tri(log2p).sm;++m) {
                            for(j=0;j<3;++j) {
                                indx = (i+j)%3;
                                for(n=0;n<NV;++n) {
                                    gbl_ptr->res.s(side(i),m,n) -= msgn*basis::tri(log2p).sfms(mode,m,j)*uht(n)(indx);
                                }
                            }
                            msgn *= sign(i);
                        }
                    }
                }
                else {
                    for(i=0;i<3;++i) {
                        side(i) = td(tind).side(i);
                        sign(i) = td(tind).sign(i);
                        sgn      = (mode % 2 ? sign(i) : 1);
                        for(n=0;n<NV;++n)
                            uht(n)(i) = sgn*gbl_ptr->res.s(side(i),mode,n);
                    }
                    
                    /* REMOVE MODES J,K FROM MODE I,M */
                    for(i=0;i<3;++i) {
                        msgn = ( (mode +1) % 2 ? sign(i) : 1);
                        for(m=mode+1;m<basis::tri(log2p).sm;++m) {
                            for(j=0;j<3;++j) {
                                indx = (i+j)%3;
                                for(n=0;n<NV;++n) {
                                    gbl_ptr->res.s(side(i),m,n) -= msgn*basis::tri(log2p).sfms(mode,m,j)*gbl_ptr->tprcn_ut(tind,n,0)*uht(0)(indx);
                                    for(k=1;k<NV;++k) {
                                        gbl_ptr->res.s(side(i),m,n) -= msgn*basis::tri(log2p).sfms(mode,m,j)*gbl_ptr->tprcn_ut(tind,n,k)*uht(k)(indx);
                                    }
                                }
                            }
                            msgn *= sign(i);
                        }
                    }
                }
            }
            
            /* GO BACK AND DO NEXT MODE */
            excpt = 2;
            ++mode;
            return(block::advance);
        }     
        
        case(6): {
#ifdef CTRL_DEBUG
            *sim::log << idprefix << "step 6 of tri_hp::minvrt: ctrl_message: " << ctrl_message << " excpt: " << excpt << " mode " << mode << std::endl;
#endif
            if (basis::tri(log2p).im == 0) return(block::advance);
            
            /* SOLVE FOR INTERIOR MODES */
            for(tind = 0; tind < ntri; ++tind) {
                DPBTRSNU2(&basis::tri(log2p).idiag(0,0),basis::tri(log2p).ibwth+1,basis::tri(log2p).im,basis::tri(log2p).ibwth,&(gbl_ptr->res.i(tind,0,0)),NV);
                restouht_bdry(tind);
                if (gbl_ptr->diagonal_preconditioner) {
                    for(k=0;k<basis::tri(log2p).im;++k) {
                        gbl_ptr->res.i(tind,k,Range::all()) /= gbl_ptr->tprcn(tind,Range::all());
                        
                        for (i=0;i<basis::tri(log2p).bm;++i)
                            for(n=0;n<NV;++n) 
                                gbl_ptr->res.i(tind,k,n) -= basis::tri(log2p).bfmi(i,k)*uht(n)(i);
                    }
                }
                else {
                    for(k=0;k<basis::tri(log2p).im;++k) {
                        /* SUBTRACT BOUNDARY MODES (bfmi is multipled by interior inverse matrix so do this after DPBSLN) */
                        for (i=0;i<basis::tri(log2p).bm;++i) {
                            for(n=0;n<NV;++n) {
                                for(m=0;m<NV;++m) {
                                    gbl_ptr->res.i(tind,k,n) -= basis::tri(log2p).bfmi(i,k)*uht(m)(i)*gbl_ptr->tprcn_ut(tind,n,m);
                                }
                            }
                        }
                        /* INVERT PRECONDITIONER (ASSUMES LOWER TRIANGULAR) */
                        DGETUS(&gbl_ptr->tprcn_ut(tind,0,0), NV, NV, &gbl_ptr->res.i(tind,k,0));
                    }
                }
            }
        }
    }
    return(block::stop);
}

void tri_hp::restouht_bdry(int tind) {
     int i,m,n,indx,cnt;
     int sign, msgn;
    
    for (i=0; i<3; ++i) {
        indx = td(tind).vrtx(i);
        for(n=0; n<NV; ++n)
            uht(n)(i) = gbl_ptr->res.v(indx,n);
    }

    cnt = 3;
    for(i=0;i<3;++i) {
        indx = td(tind).side(i);
        sign = td(tind).sign(i);
        msgn = 1;
        for (m = 0; m < basis::tri(log2p).sm; ++m) {
            for(n=0; n<NV; ++n)
                uht(n)(cnt) = msgn*gbl_ptr->res.s(indx,m,n);
            msgn *= sign;
            ++cnt;
        }
    }
      
    return;
}

block::ctrl tri_hp::setup_preconditioner(block::ctrl ctrl_message) {
    int i;
    block::ctrl state;

    /* GENERIC ROUTINE TO DO PARALLEL COMMUNICATIONS & CLEAN-UP */
    /* PHYSICS MODULES MUST SET-UP PRECONDITIONER FIRST */
    
    if (ctrl_message == block::begin) excpt = 0;
    
    switch (excpt) {
        case(0): {
#ifdef CTRL_DEBUG
            *sim::log << idprefix << "Step 0 of tri_hp::setup_preconditioner" << ctrl_message << " excpt: " << excpt << std::endl;
#endif
            if (ctrl_message != block::advance1) {
                /* SET UP TSTEP FOR ACTIVE BOUNDARIES */
                if (mmovement == coupled_deformable && log2p == 0) {
                    state = r_mesh::setup_preconditioner(ctrl_message);    
                    if (state != block::stop) return(state);
                }
                return(block::advance1);
            }
            else {
                ctrl_message = block::begin;
                ++excpt;
            }
        }
        
        case(1): {
#ifdef CTRL_DEBUG
            *sim::log << idprefix << "Step 1 of tri_hp::setup_preconditioner" << ctrl_message << " excpt: " << excpt << std::endl;
#endif
            if (ctrl_message != block::advance1) {
                /* SET UP TSTEP FOR ACTIVE BOUNDARIES */
                state = mover->setup_preconditioner(ctrl_message);    
                if (state != block::stop) return(state);
                return(block::advance1);
            }
            else {
                ctrl_message = block::begin;
                ++excpt;
            }
        }
        
        case(2): {
#ifdef CTRL_DEBUG
            *sim::log << idprefix << "Step 2 of tri_hp::setup_preconditioner" << ctrl_message << " excpt: " << excpt << std::endl;
#endif
            if (ctrl_message != block::advance1) {
                state = block::stop;
                /* SET UP TSTEP FOR ACTIVE BOUNDARIES */                
                for(i=0;i<nsbd;++i)
                    state &= hp_sbdry(i)->setup_preconditioner(ctrl_message);

                if (state != block::stop) return(state);
                return(block::advance1);
            }
            else {
                ++excpt;
            }
        }
        
        case(3): {
#ifdef CTRL_DEBUG
            *sim::log << idprefix << "Step 3 of tri_hp::setup_preconditioner" << ctrl_message << " excpt: " << excpt << std::endl;;
#endif
            mp_phase = -1;
            stage = 0;
            ++excpt;
            ctrl_message = block::stay;
        }
        
        case(4): {
#ifdef CTRL_DEBUG
            *sim::log << idprefix << "Step 4 of tri_hp::setup_preconditioner" << ctrl_message << " excpt: " << excpt << std::endl;
#endif
            if (gbl_ptr->diagonal_preconditioner) {
                if (ctrl_message == block::stay) {
                    ++mp_phase;
                    switch(mp_phase%3) {
                        case(0):
                            vc0load(mp_phase/3,gbl_ptr->vprcn.data());
                            return(block::stay);
                        case(1):
                            vmsgpass(boundary::all_phased,mp_phase/3,boundary::symmetric);
                            return(block::stay);
                        case(2):
                            return(static_cast<block::ctrl>(vc0wait_rcv(mp_phase/3,gbl_ptr->vprcn.data())));
                    }
                }
                else {
                    mp_phase = -1;
                    ++excpt;
                    ctrl_message = block::stay;
                }
            }
            else {
                /* NEED STUFF HERE FOR CONTINUITY OF MATRIX PRECONDITIONER */
                if (ctrl_message == block::stay) {
                    ++mp_phase;
                    switch(mp_phase%3) {
                        case(0):
                            vc0load(mp_phase/3,gbl_ptr->vprcn_ut.data() +stage*NV,NV);
                            return(block::stay);
                        case(1):
                            vmsgpass(boundary::all_phased,mp_phase/3,boundary::symmetric);
                            return(block::stay);
                        case(2):
                            return(static_cast<block::ctrl>(vc0wait_rcv(mp_phase/3,gbl_ptr->vprcn_ut.data()+stage*NV,NV)));
                    }
                }
                else {
                    mp_phase = -1;
                    ctrl_message = block::stay;
                    if (++stage < NV) return(ctrl_message);
                    ++excpt;
                    stage = 0;
                }
            }
        }
        
        case(5): {
#ifdef CTRL_DEBUG
            *sim::log << idprefix << "Step 5 of tri_hp::setup_preconditioner" << ctrl_message << " excpt: " << excpt << std::endl;
#endif
            if (gbl_ptr->diagonal_preconditioner) {
                if (ctrl_message == block::stay && log2p > 0) {
                    ++mp_phase;
                    switch(mp_phase%3) {
                        case(0):
                            sc0load(gbl_ptr->sprcn.data(),0,0,1);
                            return(block::stay);
                        case(1):
                            smsgpass(boundary::all,0,boundary::symmetric);
                            return(block::stay);
                        case(2):
                            sc0wait_rcv(gbl_ptr->sprcn.data(),0,0,1);
                            return(block::advance);
                    }
                }
                else {
                    ++excpt;
                }
            }
            else {
                if (ctrl_message == block::stay && log2p > 0) {
                    ++mp_phase;
                    switch(mp_phase%3) {
                        case(0):
                            sc0load(gbl_ptr->sprcn_ut.data()+stage*NV,0,0,NV);
                            return(block::stay);
                        case(1):
                            smsgpass(boundary::all,0,boundary::symmetric);
                            return(block::stay);
                        case(2):
                            sc0wait_rcv(gbl_ptr->sprcn_ut.data()+stage*NV,0,0,NV);
                            return(block::advance);
                    }
                }
                else {
                    mp_phase = -1;
                    ctrl_message = block::stay;
                    if (++stage < NV) return(ctrl_message);
                    ++excpt;
                    stage = 0;
                }
            }
        }
        
        case(6): {
#ifdef CTRL_DEBUG
            *sim::log << idprefix << "Step 6 of tri_hp::setup_preconditioner" << ctrl_message << " excpt: " << excpt << std::endl;
#endif

            if (gbl_ptr->diagonal_preconditioner) {
                /* PREINVERT PRECONDITIONER FOR VERTICES */
                gbl_ptr->vprcn(Range(0,nvrtx-1),Range::all()) = 1.0/(basis::tri(log2p).vdiag*gbl_ptr->vprcn(Range(0,nvrtx-1),Range::all()));
              
                if (basis::tri(log2p).sm > 0) {
                    /* INVERT DIAGANOL PRECONDITIONER FOR SIDES */                
                    gbl_ptr->sprcn(Range(0,nside-1),Range::all()) = 1.0/gbl_ptr->sprcn(Range(0,nside-1),Range::all());
                }
            }
            else {
                /* FACTORIZE PRECONDITIONER FOR VERTICES ASSUMES LOWER TRIANGULAR NOTHING  */
    //            for(i=0;i<nvrtx;++i)
    //                for(int n=0;n<NV;++n)
    //                    gbl_ptr->vprcn_ut(i,n,n) = 1.0/(basis::tri(log2p).vdiag*gbl_ptr->vprcn_ut(i,n,n));
    //          
    //            if (basis::tri(log2p).sm > 0) {
    //                /* INVERT DIAGANOL PRECONDITIONER FOR SIDES ASSUMES LOWER TRIANGULAR */     
    //                for(i=0;i<nside;++i)
    //                    for(int n=0;n<NV;++n)
    //                        gbl_ptr->sprcn_ut(i,n,n)= 1.0/gbl_ptr->sprcn_ut(i,n,n);
    //            }
            }
            ++excpt;
            ctrl_message = block::begin;
        }
    }
    return(block::stop);
}

    

block::ctrl tri_hp::minvrt_test(block::ctrl ctrl_message) {
    int i,j,k,m,n,tind,side,indx,indx1;
    TinyVector<int,3> v;
    block::ctrl state;
    TinyVector<FLT,mesh::ND> pt;
    FLT jcb,dtstari;
    
    if (ctrl_message == block::begin) excpt1 = 0;
    
    switch (excpt1) {
        case(0): {
            /*	SET TIME STEP TO BE 1 */
            for(tind = 0; tind < ntri; ++tind) {
                jcb = 0.25*area(tind);
                v = td(tind).vrtx;

                /* SET UP DIAGONAL PRECONDITIONER */
                dtstari = jcb*RAD((vrtx(v(0))(0) +vrtx(v(1))(0) +vrtx(v(2))(0))/3.);
                gbl_ptr->tprcn(tind,Range::all()) = dtstari;        
                
                for(i=0;i<3;++i) {
                    gbl_ptr->vprcn(v(i),Range::all())  += gbl_ptr->tprcn(tind,Range::all());
                    if (basis::tri(log2p).sm > 0) {
                        side = td(tind).side(i);
                        gbl_ptr->sprcn(side,Range::all()) += gbl_ptr->tprcn(tind,Range::all());
                    }
                }
            }
            ++excpt1;
            ctrl_message = block::begin;
        }
    
        case(1): {
            state = setup_preconditioner(ctrl_message);
            if (state != block::stop) return(state);
            ++excpt1;
        }
            
        case(2): {
            gbl_ptr->res.v(Range(0,nvrtx-1),Range::all()) = 0.0;
            gbl_ptr->res.s(Range(0,nside-1),Range::all(),Range::all()) = 0.0;
            gbl_ptr->res.i(Range(0,ntri-1),Range::all(),Range::all()) = 0.0;
            
            ug.v(Range(0,nvrtx-1),Range::all()) = 0.0;
            ug.s(Range(0,nside-1),Range::all(),Range::all()) = 0.0;
            ug.i(Range(0,ntri-1),Range::all(),Range::all()) = 0.0;

            for(tind = 0; tind<ntri;++tind) {
            
                if (td(tind).info > -1) {
                    crdtocht(tind);
                    for(n=0;n<ND;++n)
                        basis::tri(log2p).proj_bdry(&cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), &dcrd(n,1)(0,0),MXGP);

                }
                else {
                    for(n=0;n<ND;++n)
                        basis::tri(log2p).proj(vrtx(td(tind).vrtx(0))(n),vrtx(td(tind).vrtx(1))(n),vrtx(td(tind).vrtx(2))(n),&crd(n)(0,0),MXGP);

                    for(i=0;i<basis::tri(log2p).gpx;++i) {
                        for(j=0;j<basis::tri(log2p).gpn;++j) {
                            for(n=0;n<ND;++n) {
                                dcrd(n,0)(i,j) = 0.5*(vrtx(td(tind).vrtx(1))(n) -vrtx(td(tind).vrtx(0))(n));
                                dcrd(n,1)(i,j) = 0.5*(vrtx(td(tind).vrtx(2))(n) -vrtx(td(tind).vrtx(0))(n));
                            }
                        }
                    }
                }
                 
                for(n=0;n<NV;++n)
                    for(i=0;i<basis::tri(log2p).tm;++i)
                        lf(n)(i) = 0.0;

                for(i=0;i<basis::tri(log2p).gpx;++i) {
                    for(j=0;j<basis::tri(log2p).gpn;++j) {
                        pt(0) = crd(0)(i,j);
                        pt(1) = crd(1)(i,j);
                        cjcb(i,j) = dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j);
                        for(n=0;n<NV;++n)
                            res(n)(i,j) = RAD(crd(0)(i,j))*gbl_ptr->ibc->f(n,pt)*cjcb(i,j);
                    }
                }
                for(n=0;n<NV;++n)
                    basis::tri(log2p).intgrt(&lf(n)(0),&res(n)(0,0),MXGP);
                                  
                lftog(tind,gbl_ptr->res);
            }
            ctrl_message = block::begin;
            ++excpt1;
        }
        
        
        case(3): {
            state = minvrt(ctrl_message);
            if (state != block::stop) return(state);
            ++excpt1;
        }
            
        case(4): {
            /* Inversion finished */
            ug.v(Range(0,nvrtx-1),Range::all()) = gbl_ptr->res.v(Range(0,nvrtx-1),Range::all());

            if (basis::tri(log2p).sm > 0) {
                ug.s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all()) = gbl_ptr->res.s(Range(0,nside-1),Range(0,basis::tri(log2p).sm-1),Range::all());
         
                if (basis::tri(log2p).im > 0) {

                    for(i=0;i<ntri;++i) {
                        indx = 0;
                        indx1 = 0;
                        for(m=1;m<basis::tri(log2p).sm;++m) {
                            for(k=0;k<basis::tri(log2p).sm-m;++k) {
                                for(n=0;n<NV;++n) {
                                    ug.i(i,indx1,n) = gbl_ptr->res.i(i,indx,n);
                                }
                                ++indx; ++indx1;
                            }
                            indx1 += sm0 -basis::tri(log2p).sm;
                        }
                    }
                }
            }
            ++excpt1;
        }
    }
    return(block::stop);
}


