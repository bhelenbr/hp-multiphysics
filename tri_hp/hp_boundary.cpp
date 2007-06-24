/*
 *  hp_boundary.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/2/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp.h"
#include "hp_boundary.h"

#define NO_MPDEBUG

hp_vrtx_bdry* tri_hp::getnewvrtxobject(int bnum, input_map &bdrydata) {
    hp_vrtx_bdry *temp = new hp_vrtx_bdry(*this,*vbdry(bnum));    
    return(temp);
}

hp_side_bdry* tri_hp::getnewsideobject(int bnum, input_map &bdrydata) {
    hp_side_bdry *temp = new hp_side_bdry(*this,*sbdry(bnum));
    return(temp);
}

FLT hp_vrtx_bdry::dummy;

void hp_side_bdry::smatchsolution_rcv(FLT *sdata, int bgn, int end, int stride) {
    /* CAN'T USE sfinalrcv BECAUSE OF CHANGING SIGNS */
    int j,k,m,n,count,countdn,countup,offset,sind,sign;
    FLT mtchinv;
    
    if (!base.is_comm()) return;

    int matches = 1;
    int bgnsign = (bgn % 2 ? -1 : 1);
    
    /* ASSUMES REVERSE ORDERING OF SIDES */
    for(m=0;m<base.matches();++m) {    
            
        ++matches;
        
        int ebp1 = end-bgn+1;
        countdn = (base.nel-1)*ebp1*x.NV;
        countup = 0;
        for(j=0;j<base.nel;++j) {
            sign = bgnsign;
            for(k=0;k<ebp1;++k) {
                for(n=0;n<x.NV;++n)
                    base.fsndbuf(countup++) += sign*base.frcvbuf(m,countdn++);
                sign *= -1;
            }
            countdn -= 2*ebp1*x.NV;
        }
    }
    
    if (matches > 1) {
        mtchinv = 1./matches;

#ifdef MPDEBUG
        *sim::log << "side finalrcv"  << base.idnum << " " << base.is_frst() << std::endl;
#endif
        count = 0;
        for(j=0;j<base.nel;++j) {
            sind = base.el(j);
            offset = (sind*stride +bgn)*x.NV;
            for (k=bgn;k<=end;++k) {
                for(n=0;n<x.NV;++n) {
                    sdata[offset++] = base.fsndbuf(count++)*mtchinv;
#ifdef MPDEBUG
                    *sim::log << "\t" << sdata[offset-1] << std::endl;
#endif
                }
            }
        }
    }
    return;
}

void hp_side_bdry::copy_data(const hp_side_bdry &bin) {
    
    if (!curved || !x.sm0) return;
    
    for(int i=0;i<sim::nadapt; ++i)
        crvbd(i)(Range(0,base.nel-1),Range::all()) = bin.crvbd(i)(Range(0,base.nel-1),Range::all());
}

void hp_side_bdry::init(input_map& inmap,void* &gbl_in) {
    int i;
    bool coarse;
    std::string keyword;
    std::istringstream data;
    std::string filename;
    
    keyword = base.idprefix + "_curved";
    inmap.getwdefault(keyword,curved,false);

    keyword = x.idprefix + "_coarse";
    inmap.getwdefault(keyword,coarse,false);
    
    keyword = base.idprefix + "_coupled";
    inmap.getwdefault(keyword,coupled,false);
    
    if (curved && !coarse) {
        crv.resize(base.maxel,x.sm0);
        for(i=1;i<sim::nhist+1;++i)
            crvbd(i).resize(base.maxel,x.sm0);
        crvbd(0).reference(crv);
    }
    
    dxdt.resize(x.log2pmax+1,base.maxel);
        
    base.resize_buffers(base.maxel*(x.sm0+2)*x.NV);

    return;
}

void hp_side_bdry::output(std::ostream& fout, tri_hp::filetype typ,int tlvl) {
    int j,m,n;
    
    switch(typ) {
        case(tri_hp::text):
            fout << base.idprefix << " " << mytype << std::endl;
            if (curved) {
                fout << "p0: " << x.p0 << std::endl;
                
                for(j=0;j<base.nel;++j) {
                    for(m=0;m<x.sm0;++m) {
                        for(n=0;n<mesh::ND;++n)
                            fout << crvbd(tlvl)(j,m)(n) << ' ';
                        fout << std::endl;
                    }
                }
            }
            
        default:
            break;
    }
    return;
}

#define NO_RESTARTFROMOLD
void hp_side_bdry::input(ifstream& fin,tri_hp::filetype typ,int tlvl) {
    int j,m,n,pmin;
    std::string idin, mytypein;
#ifdef RESTARTFROMOLD
    char buff[100];
#endif
    
    switch(typ) {
        case(tri_hp::text):
#ifdef RESTARTFROMOLD
            fin.getline(buff,100);
            pmin = x.p0;
#else
            fin >> idin >> mytypein;
#endif
            if (curved) { 
#ifndef RESTARTFROMOLD
                fin.ignore(80,':');
                fin >>  pmin;
#endif
                pmin = x.p0;
                for(j=0;j<base.nel;++j) {
                    for(m=0;m<pmin -1;++m) {
                        for(n=0;n<mesh::ND;++n)
                            fin >> crvbd(tlvl)(j,m)(n);
                    }
                    for(m=pmin-1;m<x.sm0;++m) {
                        for(n=0;n<mesh::ND;++n)
                            crvbd(tlvl)(j,m)(n) = 0.0;
                    }
                }
            }
#ifdef RESTARTFROMOLD
            else {
                while (fin.get(buff[0])) {
                    if (buff[0] == 'B') break;
                    fin.ignore(80,'\n');
                }
            }
#endif
            break;
            
        default:
            break;
    }
}

void hp_side_bdry::curv_init(int tlvl) {
    int i,j,m,n,v0,v1,sind,info;
    TinyVector<FLT,2> pt;
    char uplo[] = "U";

    if (!curved) return;
    
    /* SKIP END VERTICES */
    for(j=1;j<base.nel;++j) {
        sind = base.el(j);
        v0 = x.sd(sind).vrtx(0);
        base.mvpttobdry(j,-1.0, x.vrtx(v0));
    }
//    v0 = x.sd(sind).vrtx(1);
//    base.mvpttobdry(base.nel-1,1.0, x.vrtx(v0));

    if (basis::tri(x.log2p).p == 1) return;
        
    /*****************************/
    /* SET UP HIGHER ORDER MODES */
    /*****************************/
    for(j=0;j<base.nel;++j) {
        sind = base.el(j);

        v0 = x.sd(sind).vrtx(0);
        v1 = x.sd(sind).vrtx(1);
        
        for(n=0;n<mesh::ND;++n) 
            basis::tri(x.log2p).proj1d(x.vrtx(v0)(n),x.vrtx(v1)(n),&x.crd(n)(0,0));
    
        for(i=0;i<basis::tri(x.log2p).gpx;++i) {
            pt(0) = x.crd(0)(0,i);
            pt(1) = x.crd(1)(0,i);
            base.mvpttobdry(j,basis::tri(x.log2p).xp(i),pt);
            x.crd(0)(0,i) -= pt(0);
            x.crd(1)(0,i) -= pt(1);
        }
        
        for(n=0;n<mesh::ND;++n) {
            basis::tri(x.log2p).intgrt1d(&x.cf(n,0),&x.crd(n)(0,0));
            DPBTRS(uplo,basis::tri(x.log2p).sm,basis::tri(x.log2p).sbwth,1,&basis::tri(x.log2p).sdiag1d(0,0),basis::tri(x.log2p).sbwth+1,&x.cf(n,2),basis::tri(x.log2p).sm,info);
        
            for(m=0;m<basis::tri(x.log2p).sm;++m)
                crvbd(tlvl)(j,m)(n) = -x.cf(n,m+2);
        }
        
        /* TEST FOR A CIRCLE 
        basis::tri(x.log2p).ptvalues1d(0.0);
        basis::tri(x.log2p).ptprobe1d(mesh::ND,&pt(0),&x.cht(0,0),MXTM);
        *sim::log << pt << ' ' << pt(0)*pt(0) +pt(1)*pt(1) << std::endl;
        */

    }
    return;
}

void hp_side_bdry::findandmovebdrypt(TinyVector<FLT,2>& xp,int &bel,FLT &psi) const {
    int sind,v0,v1,iter;
    FLT dx,dy,ol,roundoff,dpsi;
    TinyVector<FLT,2> pt;
        
    base.findbdrypt(xp,bel,psi);
    if (!curved) {
        base.side_bdry::mvpttobdry(bel,psi,xp);
        basis::tri(x.log2p).ptvalues1d(psi);
        return;
    }
    
    sind = base.el(bel);
    v0 = x.sd(sind).vrtx(0);
    v1 = x.sd(sind).vrtx(1);
    dx = x.vrtx(v1)(0) - x.vrtx(v0)(0);
    dy = x.vrtx(v1)(1) - x.vrtx(v0)(1);
    ol = 2./(dx*dx +dy*dy);
    dx *= ol;
    dy *= ol;
    
    /* FIND PSI SUCH THAT TANGENTIAL POSITION ALONG LINEAR SIDE STAYS THE SAME */
    /* THIS WAY, MULTIPLE CALLS WILL NOT GIVE DIFFERENT RESULTS */ 
    x.crdtocht1d(sind);
    
    iter = 0;
    roundoff = 10.0*EPSILON*(1.0 +(fabs(xp(0)*dx) +fabs(xp(1)*dy)));
    do {
        basis::tri(x.log2p).ptprobe1d(x.ND,pt.data(),psi,&x.cht(0,0),MXTM);
        dpsi = (pt(0) -xp(0))*dx +(pt(1) -xp(1))*dy;
        psi -= dpsi;
        if (iter++ > 100) {
            *sim::log << "#Warning: max iterations for curved side in bdry_locate type: " << base.idnum << " el: " << bel << " sind: " << sind << " loc: " << xp << " dpsi: " << dpsi << std::endl;
            break;
        }  
    } while (fabs(dpsi) > roundoff);
    xp = pt;
}

void hp_side_bdry::mvpttobdry(int bel,FLT psi,TinyVector<FLT,2> &xp) {

    /* SOLUTION IS BEING ADAPTED MUST GET INFO FROM ADAPT STORAGE */
    /* FIRST GET LINEAR APPROXIMATION TO LOCATION */
    base.side_bdry::mvpttobdry(bel, psi, xp);
    adapt_storage->findandmovebdrypt(xp,bel,psi);
    
    return;
}


block::ctrl hp_side_bdry::tadvance(bool coarse, block::ctrl ctrl_message) {
    int stage = sim::substep +sim::esdirk;  
        
    if (ctrl_message == block::begin) {
        if (x.p0 > 1 && curved) {
            if (stage) {
                /* BACK CALCULATE K TERM */
                for(int j=0;j<base.nel;++j) {
                    for(int m=0;m<basis::tri(x.log2p).sm;++m) {
                        for(int n=0;n<mesh::ND;++n)
                            crvbd(stage+1)(j,m)(n) = (crvbd(0)(j,m)(n)-crvbd(1)(j,m)(n))*sim::adirk[stage-1][stage-1];
                    }
                }
            }
            
            if (sim::substep == 0) {
                /* STORE TILDE W */
                for(int j=0;j<base.nel;++j) {
                    for(int m=0;m<basis::tri(x.log2p).sm;++m) {
                        for(int n=0;n<mesh::ND;++n)
                            crvbd(1)(j,m)(n) = crv(j,m)(n);
                    }
                }
            }
            
            /* UPDATE TILDE W */
            for (int s=0;s<stage;++s) {            
                for(int j=0;j<base.nel;++j) {
                    for(int m=0;m<basis::tri(x.log2p).sm;++m) {
                        for(int n=0;n<mesh::ND;++n) {
                            crvbd(1)(j,m)(n) += sim::adirk[stage][s]*crvbd(s+2)(j,m)(n);
                        }
                    }
                }
            }
        }
        
        calculate_unsteady_sources(coarse);
        
        /* EXTRAPOLATE SOLUTION OR MOVE TO NEXT TIME POSITION */
        if (!coupled && curved) curv_init();
    }
    return(block::stop);
}

void hp_side_bdry::calculate_unsteady_sources(bool coarse) {
    int i,j,n,sind;
    
    for(i=0;i<=x.log2pmax;++i) {
        for(j=0;j<base.nel;++j) {
            sind = base.el(j);
            x.crdtocht1d(sind,1);
            for(n=0;n<mesh::ND;++n)
                basis::tri(i).proj1d(&x.cht(n,0),&dxdt(i,j)(n,0));
        }
    }
    
    return;
}




void tri_hp::vc0load(int phase, FLT *vdata, int vrtstride) {
    int i;
            
    /* SEND COMMUNICATIONS TO ADJACENT MESHES */\
    for(i=0;i<nsbd;++i) 
        hp_sbdry(i)->vmatchsolution_snd(phase,vdata,vrtstride);
    for(i=0;i<nvbd;++i)
        hp_vbdry(i)->vmatchsolution_snd(phase,vdata,vrtstride);
    
    for(i=0;i<nsbd;++i)
        sbdry(i)->comm_prepare(boundary::all_phased,phase,boundary::symmetric);
    for(i=0;i<nvbd;++i)
        vbdry(i)->comm_prepare(boundary::all_phased,phase,boundary::symmetric);
    
    return;
}
int tri_hp::vc0wait_rcv(int phase, FLT *vdata, int vrtstride) {
    int stop = 1;
    int i;
        
    for(i=0;i<nsbd;++i) {
        stop &= sbdry(i)->comm_wait(boundary::all_phased,phase,boundary::symmetric);
    }
            
    for(i=0;i<nvbd;++i) {
        stop &= vbdry(i)->comm_wait(boundary::all_phased,phase,boundary::symmetric);
    }
        

    for(i=0;i<nsbd;++i) 
        hp_sbdry(i)->vmatchsolution_rcv(phase,vdata,vrtstride);
    for(i=0;i<nvbd;++i)
        hp_vbdry(i)->vmatchsolution_rcv(phase,vdata,vrtstride);
        
    return(stop);
}

int tri_hp::vc0rcv(int phase, FLT *vdata, int vrtstride) {
    int stop = 1,i;
    
    for(i=0;i<nsbd;++i)
        stop &= sbdry(i)->comm_nowait(boundary::all_phased,phase,boundary::symmetric);
    for(i=0;i<nvbd;++i)
        stop &= vbdry(i)->comm_nowait(boundary::all_phased,phase,boundary::symmetric);
        
    for(i=0;i<nsbd;++i) 
        hp_sbdry(i)->vmatchsolution_rcv(phase,vdata,vrtstride);
    for(i=0;i<nvbd;++i)
        hp_vbdry(i)->vmatchsolution_rcv(phase,vdata,vrtstride);
        
    return(stop);
}

void tri_hp::sc0load(FLT *sdata, int bgnmode, int endmode, int modestride) {
    int i;
        
    /* SEND COMMUNICATIONS TO ADJACENT MESHES */\
    for(i=0;i<nsbd;++i) 
        hp_sbdry(i)->smatchsolution_snd(sdata,bgnmode,endmode,modestride);
    
    for(i=0;i<nsbd;++i)
        sbdry(i)->comm_prepare(boundary::all,0,boundary::symmetric);
    
    return;
}

int tri_hp::sc0wait_rcv(FLT *sdata, int bgnmode, int endmode, int modestride) {
    int stop = 1;
    int i;
    
    for(i=0;i<nsbd;++i) {
        stop &= sbdry(i)->comm_wait(boundary::all,0,boundary::symmetric);
    }

    for(i=0;i<nsbd;++i) 
        hp_sbdry(i)->smatchsolution_rcv(sdata,bgnmode,endmode,modestride);
        
    return(stop);
}

int tri_hp::sc0rcv(FLT *sdata, int bgnmode, int endmode, int modestride) {
    int stop = 1,i;
    
    for(i=0;i<nsbd;++i)
        stop &= sbdry(i)->comm_nowait(boundary::all,0,boundary::symmetric);
        
    for(i=0;i<nsbd;++i) 
        hp_sbdry(i)->smatchsolution_rcv(sdata,bgnmode,endmode,modestride);
        
    return(stop);
}
    

block::ctrl tri_hp::matchboundaries(block::ctrl ctrl_message) {
    int i, m, n, msgn, bnum, count;
    block::ctrl state;
    
    if (ctrl_message == block::begin) excpt = 0;

#define NO_CTRL_DEBUG
#ifdef CTRL_DEBUG
    *sim::log << idprefix << " In tri_hp::matchboundares with excpt " << excpt << std::endl;
#endif
    
    switch(excpt) {
        case 0: {
            /* Match boundary vertices */
            state = mesh::matchboundaries(ctrl_message);
            if (state != block::stop) return(state);
            else {
                ++excpt;
                mp_phase = -1;
                ctrl_message = block::stay;
            }
        }
        case 1: {            
            if (ctrl_message == block::stay) {
                
                if (!sm0) {
                    excpt += 2;
                }
                else {
                    /* Match curved sides */
                    for(bnum=0;bnum<nsbd;++bnum) {
                        if (sbdry(bnum)->is_comm() && hp_sbdry(bnum)->is_curved()) {                
                            count = 0;
                            for(i=0;i<sbdry(bnum)->nel;++i) {
                                for(m=0;m<basis::tri(log2p).sm;++m) {
                                    for(n=0;n<ND;++n)
                                        sbdry(bnum)->fsndbuf(count++) = hp_sbdry(bnum)->crds(i,m,n);
                                }
                            }
                            sbdry(bnum)->sndsize() = count;
                            sbdry(bnum)->sndtype() = boundary::flt_msg;
                            sbdry(bnum)->comm_prepare(boundary::all,0,boundary::master_slave);
                        }
                    }
                }
            }
            ++excpt;
            return(block::advance);
        }
        
        case 2: {
            for(bnum=0;bnum<nsbd;++bnum) {
                if (sbdry(bnum)->is_comm() && hp_sbdry(bnum)->is_curved()) {                
                    sbdry(bnum)->comm_exchange(boundary::all,0,boundary::master_slave);
                }
            }
            ++excpt;
            return(block::advance);
        }
        case 3: {
            for(bnum=0;bnum<nsbd;++bnum) {
                if (sbdry(bnum)->is_comm() && hp_sbdry(bnum)->is_curved()) {                
                    sbdry(bnum)->comm_wait(boundary::all,0,boundary::master_slave);
                    
                    if (!sbdry(bnum)->is_frst()) {
                        count = 0;
                        for(i=sbdry(bnum)->nel-1;i>=0;--i) {
                            msgn = 1;
                            for(m=0;m<basis::tri(log2p).sm;++m) {
                                for(n=0;n<ND;++n)
                                    hp_sbdry(bnum)->crds(i,m,n) = msgn*sbdry(bnum)->frcvbuf(0,count++);
                                msgn *= -1;
                            }
                        }
                    }
                }
            }
            ++excpt;
            return(block::advance);
        }
        
        case 4: {
            mp_phase = -1;
            ++excpt;
            ctrl_message = block::stay;
        }
                
        case 5: {
            ++mp_phase;
            excpt += ctrl_message;
            switch(mp_phase%3) {
                case(0):
                    vc0load(mp_phase/3,ug.v.data());
                    return(block::stay);
                case(1):
                    vmsgpass(boundary::all_phased,mp_phase/3,boundary::symmetric);
                    return(block::stay);
                case(2):
                    return(static_cast<block::ctrl>(vc0wait_rcv(mp_phase/3,ug.v.data())));
            }
        }
        case 6: {
            mp_phase = -1;
            ++excpt;
            ctrl_message = block::stay;
        }
        case 7: {
            if (ctrl_message == block::stay) {

                if (!sm0) return(block::advance);
                
                ++mp_phase;
                switch(mp_phase%3) {
                    case(0):
                        sc0load(ug.s.data(),0,sm0-1,ug.s.extent(secondDim));
                        return(block::stay);
                    case(1):
                        smsgpass(boundary::all,0,boundary::symmetric);
                        return(block::stay);
                    case(2):
                        sc0wait_rcv(ug.s.data(),0,sm0-1,ug.s.extent(secondDim));
                        return(block::advance);
                }
            }
            else
                ++excpt;
        }
    }
    
    return(block::stop);
}

block::ctrl hp_side_bdry::findmax(block::ctrl ctrl_message, FLT (*fxy)(TinyVector<FLT,2> &x)) {
    FLT ddpsi1, ddpsi2, psil, psir;
    TinyVector<FLT,2> xp, dx, maxloc, minloc;
    FLT max,min;
    int v0, sind;
    
    if (ctrl_message == block::begin) excpt = 0;
    else ++excpt;
    
    switch (excpt) {
        case(0):
            /* CALCULATE SLOPE AT ENDPOINT & TRANSMIT TO NEXT SURFACE */
            sind = base.el(base.nel-1);
            x.crdtocht1d(sind);
            basis::tri(x.log2p).ptprobe1d(mesh::ND,&xp(0),&dx(0),1.0,&x.cht(0,0),MXTM);
            ddpsi2 = (*fxy)(dx);
            if (base.vbdry(1) >= 0) {
                x.vbdry(base.vbdry(1))->vloadbuff(boundary::manifolds,&ddpsi2,0,1,1);
                x.vbdry(base.vbdry(1))->comm_prepare(boundary::manifolds,0,boundary::master_slave);
            }
            if (base.vbdry(0) >= 0) {
                x.vbdry(base.vbdry(0))->comm_prepare(boundary::manifolds,0,boundary::master_slave);
            }
            return(block::advance);
        
        case(1):
            if (base.vbdry(1) >= 0) 
                x.vbdry(base.vbdry(1))->comm_exchange(boundary::manifolds,0,boundary::master_slave);
            if (base.vbdry(0) >= 0)
                x.vbdry(base.vbdry(0))->comm_exchange(boundary::manifolds,0,boundary::master_slave);
            return(block::advance);
        
        case(2):
            if (base.vbdry(1) >= 0) 
                x.vbdry(base.vbdry(1))->comm_wait(boundary::manifolds,0,boundary::master_slave);
            if (base.vbdry(0) >= 0) {
                x.vbdry(base.vbdry(0))->comm_wait(boundary::manifolds,0,boundary::master_slave);
                if (x.vbdry(base.vbdry(0))->is_comm()) 
                    ddpsi2 = x.vbdry(base.vbdry(0))->frcvbuf(0,0);
                else
                    ddpsi2 = 0.0;
            }
                
            max = -1.0e99;
            min = 1.0e99;
            for(int indx=0;indx<base.nel;++indx) {
                sind = base.el(indx);
                x.crdtocht1d(sind);
                basis::tri(x.log2p).ptprobe1d(mesh::ND, &xp(0), &dx(0), -1.0, &x.cht(0,0), MXTM);
                ddpsi1 = (*fxy)(dx);
                if (ddpsi1 * ddpsi2 <= 0.0) {
                    v0 = x.sd(base.el(indx)).vrtx(0);
                    if ((*fxy)(x.vrtx(v0)) > max) {
                        maxloc[0] = x.vrtx(v0)(0);
                        maxloc[1] = x.vrtx(v0)(1);
                        max = (*fxy)(x.vrtx(v0));
                    }
                    if ((*fxy)(x.vrtx(v0)) < min) {
                        minloc[0] = x.vrtx(v0)(0);
                        minloc[1] = x.vrtx(v0)(1);
                        min = (*fxy)(x.vrtx(v0));
                    }
                    *sim::log << "#LOCAL EXTREMA: " << x.vrtx(v0)(0) << ' ' << x.vrtx(v0)(1) << ' ' <<(*fxy)(x.vrtx(v0)) << std::endl;
                }
                basis::tri(x.log2p).ptprobe1d(mesh::ND, &xp(0), &dx(0), 1.0, &x.cht(0,0), MXTM);
                ddpsi2 = (*fxy)(dx);
                if (ddpsi1 *ddpsi2 <= 0.0) {
                    /* INTERIOR MAXIMUM */
                    psil = -1.0;
                    psir = 1.0;
                    while (psir-psil > 1.0e-10) {
                        basis::tri(x.log2p).ptprobe1d(mesh::ND, &xp(0), &dx(0), 0.5*(psil +psir), &x.cht(0,0), MXTM);
                        if ((*fxy)(dx)*ddpsi1 < 0.0) 
                            psir = 0.5*(psil+psir);
                        else
                            psil = 0.5*(psil+psir);
                    }
                    if ((*fxy)(xp) > max) {
                        maxloc[0] = xp[0];
                        maxloc[1] = xp[1];
                        max = (*fxy)(xp);
                    }
                    if ((*fxy)(xp) < min) {
                        minloc[0] = xp[0];
                        minloc[1] = xp[1];
                        min = (*fxy)(xp);
                    }
                    *sim::log << "#LOCAL EXTREMA: " << xp[0] << ' ' << xp[1] << ' ' << (*fxy)(xp) << std::endl;
                }  
            }
            *sim::log << "#MAX EXTREMA: " << maxloc[0] << ' ' << maxloc[1] << ' ' << max << std::endl;
            *sim::log << "#MIN EXTREMA: " << minloc[0] << ' ' << minloc[1] << ' ' << min << std::endl;
    }
    return(block::stop);
}

 void hp_side_bdry::findintercept(FLT (*fxy)(TinyVector<FLT,2> &x)) {
    FLT psil, psir;
    TinyVector<FLT,2> xp, dx;
    int v0, sind;
    FLT vl, vr;
    
    sind = base.el(0);
    x.crdtocht1d(sind);
    v0 = x.sd(sind).vrtx(0);
    vl = (*fxy)(x.vrtx(v0));

    for(int indx=0;indx<base.nel;++indx) {
        sind = base.el(indx);
        x.crdtocht1d(sind);
        v0 = x.sd(sind).vrtx(1);
        vr = (*fxy)(x.vrtx(v0));

        if (vl*vr <= 0.0) {
            /* INTERIOR INTERCEPT */
            psil = -1.0;
            psir = 1.0;
            while (psir-psil > 1.0e-10) {
                basis::tri(x.log2p).ptprobe1d(mesh::ND,&xp(0),&dx(0),0.5*(psil+psir),&x.cht(0,0),MXTM);
                if ((*fxy)(xp)*vl < 0.0) 
                    psir = 0.5*(psil+psir);
                else
                    psil = 0.5*(psil+psir);
            }
            *sim::log << "#INTERSECTION: " << xp[0] << ' ' << xp[1] << ' ' << (*fxy)(xp) << std::endl;
        }
        vl = vr; 
    }
    
    return;
}
