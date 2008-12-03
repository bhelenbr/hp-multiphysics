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
#include <blitz/tinyvec-et.h>
#include <libbinio/binwrap.h>
#include <myblas.h>

//#define MPDEBUG

hp_vrtx_bdry* tri_hp::getnewvrtxobject(int bnum, input_map &bdrydata) {
    hp_vrtx_bdry *temp = new hp_vrtx_bdry(*this,*vbdry(bnum));  
    gbl->vbdry_gbls(bnum) = temp->create_global_structure();
    return(temp);
}

hp_edge_bdry* tri_hp::getnewsideobject(int bnum, input_map &bdrydata) {
    hp_edge_bdry *temp = new hp_edge_bdry(*this,*ebdry(bnum));
    gbl->ebdry_gbls(bnum) = temp->create_global_structure();
    return(temp);
}

void hp_edge_bdry::smatchsolution_rcv(FLT *sdata, int bgn, int end, int stride) {
    /* CAN'T USE sfinalrcv BECAUSE OF CHANGING SIGNS */
    int j,k,m,n,count,countdn,countup,offset,sind,sign;
    FLT mtchinv;
    
    if (!base.is_comm()) return;

    int matches = 1;
    int bgnsign = (bgn % 2 ? -1 : 1);
    
    /* ASSUMES REVERSE ORDERING OF SIDES */
    for(m=0;m<base.nmatches();++m) {    
            
        ++matches;
        
        int ebp1 = end-bgn+1;
        countdn = (base.nseg-1)*ebp1*x.NV;
        countup = 0;
        for(j=0;j<base.nseg;++j) {
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
        *x.gbl->log << "side finalrcv"  << base.idnum << " " << base.is_frst() << std::endl;
#endif
        count = 0;
        for(j=0;j<base.nseg;++j) {
            sind = base.seg(j);
            offset = (sind*stride +bgn)*x.NV;
            for (k=bgn;k<=end;++k) {
                for(n=0;n<x.NV;++n) {
                    sdata[offset++] = base.fsndbuf(count++)*mtchinv;
#ifdef MPDEBUG
                    *x.gbl->log << "\t" << sdata[offset-1] << std::endl;
#endif
                }
            }
        }
    }
    return;
}

void hp_edge_bdry::copy(const hp_edge_bdry &bin) {
    
    if (!curved || !x.sm0) return;
    
    for(int i=0;i<sim::nadapt; ++i)
        crvbd(i)(Range(0,base.nseg-1),Range::all()) = bin.crvbd(i)(Range(0,base.nseg-1),Range::all());
}

void hp_edge_bdry::init(input_map& inmap,void* gbl_in) {
    int i;
    std::string keyword;
    std::istringstream data;
    std::string filename;
	
	if (inmap.find(base.idprefix +"_ibc") != inmap.end()) {
		ibc = x.getnewibc(base.idprefix+"_ibc",inmap);
	}
    
    keyword = base.idprefix + "_curved";
    inmap.getwdefault(keyword,curved,false);

    keyword = base.idprefix + "_coupled";
    inmap.getwdefault(keyword,coupled,false);
    
    if (curved && !x.coarse_level) {
        crv.resize(base.maxseg,x.sm0);
        for(i=1;i<sim::nhist+1;++i)
            crvbd(i).resize(base.maxseg,x.sm0);
        crvbd(0).reference(crv);
    }
    
    dxdt.resize(x.log2pmax+1,base.maxseg);
        
    base.resize_buffers(base.maxseg*(x.sm0+2)*x.NV);

    return;
}

void hp_edge_bdry::output(std::ostream& fout, tri_hp::filetype typ,int tlvl) {
    int j,m,n;

    switch(typ) {
        case(tri_hp::text):
            fout << base.idprefix << " " << mytype << std::endl;
            if (curved) {
                fout << "p0: " << x.p0 << std::endl;
                
                for(j=0;j<base.nseg;++j) {
                    for(m=0;m<x.sm0;++m) {
                        for(n=0;n<tri_mesh::ND;++n)
                            fout << crvbd(tlvl)(j,m)(n) << ' ';
                        fout << std::endl;
                    }
                }
            }
			break;
			
		case(tri_hp::binary):
			if (curved) {
				binowstream bout(&fout);
				bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::BigEndian)),1);
				bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::FloatIEEE)),1);
				bout.writeInt(x.p0,sizeof(int));
                
                for(j=0;j<base.nseg;++j) {
                    for(m=0;m<x.sm0;++m) {
                        for(n=0;n<tri_mesh::ND;++n)
                            bout.writeFloat(crvbd(tlvl)(j,m)(n),binio::Double);
                    }
                }
			}
            break;
			
        default:
            break;
    }
    return;
}

void hp_edge_bdry::input(ifstream& fin,tri_hp::filetype typ,int tlvl) {
    int j,m,n,pmin;
    std::string idin, mytypein;    
    switch(typ) {
        case(tri_hp::text):
            fin >> idin >> mytypein;
            if (curved) { 
                fin.ignore(80,':');
                fin >>  pmin;
                pmin = x.p0;
                for(j=0;j<base.nseg;++j) {
                    for(m=0;m<pmin -1;++m) {
                        for(n=0;n<tri_mesh::ND;++n)
                            fin >> crvbd(tlvl)(j,m)(n);
                    }
                    for(m=pmin-1;m<x.sm0;++m) {
                        for(n=0;n<tri_mesh::ND;++n)
                            crvbd(tlvl)(j,m)(n) = 0.0;
                    }
                }
            }
            break;
		case(tri_hp::binary):
			if (curved) {
				biniwstream bin(&fin);
				
				/* HEADER INFORMATION */
				bin.setFlag(binio::BigEndian,bin.readInt(1));
				bin.setFlag(binio::FloatIEEE,bin.readInt(1));
                
				pmin = bin.readInt(sizeof(int));
                for(j=0;j<base.nseg;++j) {
                    for(m=0;m<pmin-1;++m) {
                        for(n=0;n<tri_mesh::ND;++n)
                            crvbd(tlvl)(j,m)(n) = bin.readFloat(binio::Double);
                    }
	                for(m=pmin-1;m<x.sm0;++m) {
                        for(n=0;n<tri_mesh::ND;++n)
                            crvbd(tlvl)(j,m)(n) = 0.0;
                    }				
                }
			}
            break;
            
        default:
            break;
    }
}

void hp_edge_bdry::curv_init(int tlvl) {
    int i,j,m,n,v0,v1,sind,info;
    TinyVector<FLT,2> pt;
    char uplo[] = "U";

    if (!curved) return;
    
    /* SKIP END VERTICES */
    for(j=1;j<base.nseg;++j) {
        sind = base.seg(j);
        v0 = x.seg(sind).pnt(0);
        base.mvpttobdry(j,-1.0, x.pnts(v0));
    }
//    v0 = x.seg(sind).pnt(1);
//    base.mvpttobdry(base.nseg-1,1.0, x.pnts(v0));

    if (basis::tri(x.log2p).p == 1) return;
        
    /*****************************/
    /* SET UP HIGHER ORDER MODES */
    /*****************************/
    for(j=0;j<base.nseg;++j) {
        sind = base.seg(j);

        v0 = x.seg(sind).pnt(0);
        v1 = x.seg(sind).pnt(1);
        
        for(n=0;n<tri_mesh::ND;++n) 
            basis::tri(x.log2p).proj1d(x.pnts(v0)(n),x.pnts(v1)(n),&x.crd(n)(0,0));
    
        for(i=0;i<basis::tri(x.log2p).gpx;++i) {
            pt(0) = x.crd(0)(0,i);
            pt(1) = x.crd(1)(0,i);
            base.mvpttobdry(j,basis::tri(x.log2p).xp(i),pt);
            x.crd(0)(0,i) -= pt(0);
            x.crd(1)(0,i) -= pt(1);
        }
        
        for(n=0;n<tri_mesh::ND;++n) {
            basis::tri(x.log2p).intgrt1d(&x.cf(n,0),&x.crd(n)(0,0));
            DPBTRS(uplo,basis::tri(x.log2p).sm,basis::tri(x.log2p).sbwth,1,&basis::tri(x.log2p).sdiag1d(0,0),basis::tri(x.log2p).sbwth+1,&x.cf(n,2),basis::tri(x.log2p).sm,info);
        
            for(m=0;m<basis::tri(x.log2p).sm;++m)
                crvbd(tlvl)(j,m)(n) = -x.cf(n,m+2);
        }
        
        /* TEST FOR A CIRCLE 
        basis::tri(x.log2p).ptvalues1d(0.0);
        basis::tri(x.log2p).ptprobe1d(tri_mesh::ND,&pt(0),&x.cht(0,0),MXTM);
        *gbl->log << pt << ' ' << pt(0)*pt(0) +pt(1)*pt(1) << std::endl;
        */

    }
    return;
}

void hp_edge_bdry::findandmovebdrypt(TinyVector<FLT,2>& xp,int &bel,FLT &psi) const {
    int sind,v0,v1,iter;
    FLT dx,dy,ol,roundoff,dpsi;
    TinyVector<FLT,2> pt;
        
    base.findbdrypt(xp,bel,psi);
    if (!curved) {
        base.edge_bdry::mvpttobdry(bel,psi,xp);
        basis::tri(x.log2p).ptvalues1d(psi);
        return;
    }
    
    sind = base.seg(bel);
    v0 = x.seg(sind).pnt(0);
    v1 = x.seg(sind).pnt(1);
    dx = x.pnts(v1)(0) - x.pnts(v0)(0);
    dy = x.pnts(v1)(1) - x.pnts(v0)(1);
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
            *x.gbl->log << "#Warning: max iterations for curved side in bdry_locate type: " << base.idnum << " seg: " << bel << " sind: " << sind << " loc: " << xp << " dpsi: " << dpsi << std::endl;
            break;
        }  
    } while (fabs(dpsi) > roundoff);
    xp = pt;
}

void hp_edge_bdry::mvpttobdry(int bel,FLT psi,TinyVector<FLT,2> &xp) {

    /* SOLUTION IS BEING ADAPTED MUST GET INFO FROM ADAPT STORAGE */
    /* FIRST GET LINEAR APPROXIMATION TO LOCATION */
    base.edge_bdry::mvpttobdry(bel, psi, xp);
    adapt_storage->findandmovebdrypt(xp,bel,psi);
    
    return;
}


void hp_edge_bdry::tadvance() {
    int stage = x.gbl->substep +x.gbl->esdirk;  
        
    if (x.p0 > 1 && curved) {
        if (stage) {
            /* BACK CALCULATE K TERM */
            for(int j=0;j<base.nseg;++j) {
                for(int m=0;m<basis::tri(x.log2p).sm;++m) {
                    for(int n=0;n<tri_mesh::ND;++n)
                        crvbd(stage+1)(j,m)(n) = (crvbd(0)(j,m)(n)-crvbd(1)(j,m)(n))*x.gbl->adirk[stage-1][stage-1];
                }
            }
        }
        
        if (x.gbl->substep == 0) {
            /* STORE TILDE W */
            for(int j=0;j<base.nseg;++j) {
                for(int m=0;m<basis::tri(x.log2p).sm;++m) {
                    for(int n=0;n<tri_mesh::ND;++n)
                        crvbd(1)(j,m)(n) = crv(j,m)(n);
                }
            }
        }
        
        /* UPDATE TILDE W */
        for (int s=0;s<stage;++s) {            
            for(int j=0;j<base.nseg;++j) {
                for(int m=0;m<basis::tri(x.log2p).sm;++m) {
                    for(int n=0;n<tri_mesh::ND;++n) {
                        crvbd(1)(j,m)(n) += x.gbl->adirk[stage][s]*crvbd(s+2)(j,m)(n);
                    }
                }
            }
        }

        /* EXTRAPOLATE GUESS? */
        if (stage && x.gbl->dti > 0.0) {
            FLT constant =  x.gbl->cdirk[x.gbl->substep];
            crvbd(0)(Range(0,base.nseg-1),Range(0,basis::tri(x.log2p).sm-1)) += constant*crvbd(stage+1)(Range(0,base.nseg-1),Range(0,basis::tri(x.log2p).sm-1));
        }
    }

    calculate_unsteady_sources();
        
    /* EXTRAPOLATE SOLUTION OR MOVE TO NEXT TIME POSITION */
    if (!coupled && curved) curv_init();
        
    return;
}

void hp_edge_bdry::calculate_unsteady_sources() {
    int i,j,n,sind;
    
    for(i=0;i<=x.log2pmax;++i) {
        for(j=0;j<base.nseg;++j) {
            sind = base.seg(j);
            x.crdtocht1d(sind,1);
            for(n=0;n<tri_mesh::ND;++n)
                basis::tri(i).proj1d(&x.cht(n,0),&dxdt(i,j)(n,0));
        }
    }
    
    return;
}

void tri_hp::vc0load(int phase, FLT *pdata, int vrtstride) {
    int i;
            
    /* SEND COMMUNICATIONS TO ADJACENT MESHES */\
    for(i=0;i<nebd;++i) 
        hp_ebdry(i)->pmatchsolution_snd(phase,pdata,vrtstride);
    for(i=0;i<nvbd;++i)
        hp_vbdry(i)->pmatchsolution_snd(phase,pdata,vrtstride);
    
    for(i=0;i<nebd;++i)
        ebdry(i)->comm_prepare(boundary::all_phased,phase,boundary::symmetric);
    for(i=0;i<nvbd;++i)
        vbdry(i)->comm_prepare(boundary::all_phased,phase,boundary::symmetric);
    
    return;
}
int tri_hp::vc0wait_rcv(int phase, FLT *pdata, int vrtstride) {
    int stop = 1;
    int i;
        
    for(i=0;i<nebd;++i) {
        stop &= ebdry(i)->comm_wait(boundary::all_phased,phase,boundary::symmetric);
    }
            
    for(i=0;i<nvbd;++i) {
        stop &= vbdry(i)->comm_wait(boundary::all_phased,phase,boundary::symmetric);
    }
        

    for(i=0;i<nebd;++i) 
        hp_ebdry(i)->pmatchsolution_rcv(phase,pdata,vrtstride);
    for(i=0;i<nvbd;++i)
        hp_vbdry(i)->pmatchsolution_rcv(phase,pdata,vrtstride);
        
    return(stop);
}

int tri_hp::vc0rcv(int phase, FLT *pdata, int vrtstride) {
    int stop = 1,i;
    
    for(i=0;i<nebd;++i)
        stop &= ebdry(i)->comm_nowait(boundary::all_phased,phase,boundary::symmetric);
    for(i=0;i<nvbd;++i)
        stop &= vbdry(i)->comm_nowait(boundary::all_phased,phase,boundary::symmetric);
        
    for(i=0;i<nebd;++i) 
        hp_ebdry(i)->pmatchsolution_rcv(phase,pdata,vrtstride);
    for(i=0;i<nvbd;++i)
        hp_vbdry(i)->pmatchsolution_rcv(phase,pdata,vrtstride);
        
    return(stop);
}

void tri_hp::sc0load(FLT *sdata, int bgnmode, int endmode, int modestride) {
    int i;
        
    /* SEND COMMUNICATIONS TO ADJACENT MESHES */\
    for(i=0;i<nebd;++i) 
        hp_ebdry(i)->smatchsolution_snd(sdata,bgnmode,endmode,modestride);
    
    for(i=0;i<nebd;++i)
        ebdry(i)->comm_prepare(boundary::all,0,boundary::symmetric);
    
    return;
}

int tri_hp::sc0wait_rcv(FLT *sdata, int bgnmode, int endmode, int modestride) {
    int stop = 1;
    int i;
    
    for(i=0;i<nebd;++i) {
        stop &= ebdry(i)->comm_wait(boundary::all,0,boundary::symmetric);
    }

    for(i=0;i<nebd;++i) 
        hp_ebdry(i)->smatchsolution_rcv(sdata,bgnmode,endmode,modestride);
        
    return(stop);
}

int tri_hp::sc0rcv(FLT *sdata, int bgnmode, int endmode, int modestride) {
    int stop = 1,i;
    
    for(i=0;i<nebd;++i)
        stop &= ebdry(i)->comm_nowait(boundary::all,0,boundary::symmetric);
        
    for(i=0;i<nebd;++i) 
        hp_ebdry(i)->smatchsolution_rcv(sdata,bgnmode,endmode,modestride);
        
    return(stop);
}
    

void tri_hp::matchboundaries() {
    int i, m, n, msgn, bnum, count;
    int last_phase, mp_phase;
        
    /* Match boundary vertices */
    tri_mesh::matchboundaries();
        
    for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
        vc0load(mp_phase,ug.v.data());
        pmsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
        last_phase = true;
        last_phase &= vc0wait_rcv(mp_phase,ug.v.data());
    }
    
    if (!sm0) return;
    
    sc0load(ug.s.data(),0,sm0-1,ug.s.extent(secondDim));
    smsgpass(boundary::all,0,boundary::symmetric);
    sc0wait_rcv(ug.s.data(),0,sm0-1,ug.s.extent(secondDim));    

    /* Match curved sides */
    for(bnum=0;bnum<nebd;++bnum) {
        if (ebdry(bnum)->is_comm() && hp_ebdry(bnum)->is_curved()) {                
            count = 0;
            for(i=0;i<ebdry(bnum)->nseg;++i) {
                for(m=0;m<basis::tri(log2p).sm;++m) {
                    for(n=0;n<ND;++n)
                        ebdry(bnum)->fsndbuf(count++) = hp_ebdry(bnum)->crds(i,m,n);
                }
            }
            ebdry(bnum)->sndsize() = count;
            ebdry(bnum)->sndtype() = boundary::flt_msg;
            ebdry(bnum)->comm_prepare(boundary::all,0,boundary::master_slave);
        }
    }
    
    for(bnum=0;bnum<nebd;++bnum) {
        if (ebdry(bnum)->is_comm() && hp_ebdry(bnum)->is_curved()) {                
            ebdry(bnum)->comm_exchange(boundary::all,0,boundary::master_slave);
        }
    }
    
    
    for(bnum=0;bnum<nebd;++bnum) {
        if (ebdry(bnum)->is_comm() && hp_ebdry(bnum)->is_curved()) {                
            ebdry(bnum)->comm_wait(boundary::all,0,boundary::master_slave);
            
            if (!ebdry(bnum)->is_frst()) {
                count = 0;
                for(i=ebdry(bnum)->nseg-1;i>=0;--i) {
                    msgn = 1;
                    for(m=0;m<basis::tri(log2p).sm;++m) {
                        for(n=0;n<ND;++n)
                            hp_ebdry(bnum)->crds(i,m,n) = msgn*ebdry(bnum)->frcvbuf(0,count++);
                        msgn *= -1;
                    }
                }
            }
        }
    }
    
    return;
}

void hp_edge_bdry::findmax(FLT (*fxy)(TinyVector<FLT,2> &x)) {
    FLT ddpsi1, ddpsi2, psil, psir;
    TinyVector<FLT,2> xp, dx, maxloc, minloc;
    FLT max,min;
    int v0, sind;
    

    /* CALCULATE SLOPE AT ENDPOINT & TRANSMIT TO NEXT SURFACE */
    sind = base.seg(base.nseg-1);
    x.crdtocht1d(sind);
    basis::tri(x.log2p).ptprobe1d(tri_mesh::ND,&xp(0),&dx(0),1.0,&x.cht(0,0),MXTM);
    ddpsi2 = (*fxy)(dx);
    if (base.vbdry(1) >= 0) {
        x.vbdry(base.vbdry(1))->vloadbuff(boundary::manifolds,&ddpsi2,0,1,1);
        x.vbdry(base.vbdry(1))->comm_prepare(boundary::manifolds,0,boundary::master_slave);
    }
    if (base.vbdry(0) >= 0) {
        x.vbdry(base.vbdry(0))->comm_prepare(boundary::manifolds,0,boundary::master_slave);
    }
        

    if (base.vbdry(1) >= 0) 
        x.vbdry(base.vbdry(1))->comm_exchange(boundary::manifolds,0,boundary::master_slave);
    if (base.vbdry(0) >= 0)
        x.vbdry(base.vbdry(0))->comm_exchange(boundary::manifolds,0,boundary::master_slave);
        
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
    for(int indx=0;indx<base.nseg;++indx) {
        sind = base.seg(indx);
        x.crdtocht1d(sind);
        basis::tri(x.log2p).ptprobe1d(tri_mesh::ND, &xp(0), &dx(0), -1.0, &x.cht(0,0), MXTM);
        ddpsi1 = (*fxy)(dx);
        if (ddpsi1 * ddpsi2 <= 0.0) {
            v0 = x.seg(base.seg(indx)).pnt(0);
            if ((*fxy)(x.pnts(v0)) > max) {
                maxloc[0] = x.pnts(v0)(0);
                maxloc[1] = x.pnts(v0)(1);
                max = (*fxy)(x.pnts(v0));
            }
            if ((*fxy)(x.pnts(v0)) < min) {
                minloc[0] = x.pnts(v0)(0);
                minloc[1] = x.pnts(v0)(1);
                min = (*fxy)(x.pnts(v0));
            }
            *x.gbl->log << "#LOCAL EXTREMA: " << x.pnts(v0)(0) << ' ' << x.pnts(v0)(1) << ' ' <<(*fxy)(x.pnts(v0)) << std::endl;
        }
        basis::tri(x.log2p).ptprobe1d(tri_mesh::ND, &xp(0), &dx(0), 1.0, &x.cht(0,0), MXTM);
        ddpsi2 = (*fxy)(dx);
        if (ddpsi1 *ddpsi2 <= 0.0) {
            /* INTERIOR MAXIMUM */
            psil = -1.0;
            psir = 1.0;
            while (psir-psil > 1.0e-10) {
                basis::tri(x.log2p).ptprobe1d(tri_mesh::ND, &xp(0), &dx(0), 0.5*(psil +psir), &x.cht(0,0), MXTM);
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
            *x.gbl->log << "#LOCAL EXTREMA: " << xp[0] << ' ' << xp[1] << ' ' << (*fxy)(xp) << std::endl;
        }  
    }
    *x.gbl->log << "#MAX EXTREMA: " << maxloc[0] << ' ' << maxloc[1] << ' ' << max << std::endl;
    *x.gbl->log << "#MIN EXTREMA: " << minloc[0] << ' ' << minloc[1] << ' ' << min << std::endl;

    return;
}

 void hp_edge_bdry::findintercept(FLT (*fxy)(TinyVector<FLT,2> &x)) {
    FLT psil, psir;
    TinyVector<FLT,2> xp, dx;
    int v0, sind;
    FLT vl, vr;
    
    sind = base.seg(0);
    x.crdtocht1d(sind);
    v0 = x.seg(sind).pnt(0);
    vl = (*fxy)(x.pnts(v0));

    for(int indx=0;indx<base.nseg;++indx) {
        sind = base.seg(indx);
        x.crdtocht1d(sind);
        v0 = x.seg(sind).pnt(1);
        vr = (*fxy)(x.pnts(v0));

        if (vl*vr <= 0.0) {
            /* INTERIOR INTERCEPT */
            psil = -1.0;
            psir = 1.0;
            while (psir-psil > 1.0e-10) {
                basis::tri(x.log2p).ptprobe1d(tri_mesh::ND,&xp(0),&dx(0),0.5*(psil+psir),&x.cht(0,0),MXTM);
                if ((*fxy)(xp)*vl < 0.0) 
                    psir = 0.5*(psil+psir);
                else
                    psil = 0.5*(psil+psir);
            }
            *x.gbl->log << "#INTERSECTION: " << xp[0] << ' ' << xp[1] << ' ' << (*fxy)(xp) << std::endl;
        }
        vl = vr; 
    }
    
    return;
}
