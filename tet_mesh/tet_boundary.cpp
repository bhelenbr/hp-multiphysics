/*
 *  boundaries.cpp
 *  mesh3d
 *
 *  Created by Michael Brazell on 7/31/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_mesh.h"
#include "tet_boundary.h"
#include <assert.h>
#include <float.h>

/********************/
/* VERTEX FUNCTIONS */
/********************/
    
/* GENERIC VERTEX COMMUNICATIONS */
void vcomm::ploadbuff(boundary::groups grp,FLT *base,int bgn,int end, int stride) {
    int i,offset;

    if (!in_group(grp)) return;

    sndsize()=end-bgn+1;
    sndtype()=flt_msg;
    
    /* LOAD SEND BUFFER */    
    offset = pnt*stride +bgn;
    for (i=0;i<end-bgn+1;++i) 
        fsndbuf(i) = base[offset+i];
}

void vcomm::pfinalrcv(boundary::groups grp, int phi, comm_type type, operation op, FLT *base, int bgn, int end, int stride) {
    int i,offset;
        
    bool reload = comm_finish(grp,phi,type,op);
    if (!reload) return;
    
    offset = pnt*stride +bgn;
    for(i=0;i<end-bgn+1;++i) 
        base[offset++] = fsndbuf(i);
        
    return;
}
 
/**************************************/
/* GENERIC FUNCTIONS FOR EDGES        */
/**************************************/
void edge_bdry::alloc(int n) {
    maxseg = n;
    seg.resize(Range(-1,n));
    prev.resize(Range(-1,n));
    next.resize(Range(-1,n));
    prev = 0;
    next = 0;
}
      
void edge_bdry::copy(const edge_bdry& bin) {
    int i;
        
    if (!maxseg) alloc(bin.maxseg);
    else assert(bin.nseg <= maxseg);
    vbdry = bin.vbdry;
    
    nseg = bin.nseg;
    
    for(i=0;i<nseg;++i)
        seg(i) = bin.seg(i);
        
    return;
}

void edge_bdry::mvpttobdry(int indx, FLT psi, TinyVector<FLT,tet_mesh::ND> &pt) {
    /* FOR A LINEAR SIDE */
    int n;
    
    for (n=0;n<tet_mesh::ND;++n)
        pt(n) = 0.5*((1. -psi)*x.pnts(x.seg(seg(indx)).pnt(0))(n) +(1.+psi)*x.pnts(x.seg(seg(indx)).pnt(1))(n));
    
    return;
}

void edge_bdry::findbdrypt(const TinyVector<FLT,tet_mesh::ND> xpt, int &sidloc, FLT &psiloc) const {
    int k,sind,p0,p1,sidlocprev;
    FLT dx,dy,ol,psi,normdist;
    FLT psiprev,normdistprev;
    FLT mindist = 1.0e32;
        
    if (x.seg(seg(0)).pnt(0) == x.seg(seg(nseg-1)).pnt(1)) {
        /* BOUNDARY IS A LOOP */
        sind = seg(nseg-1);
        p0 = x.seg(sind).pnt(0);
        p1 = x.seg(sind).pnt(1);
        dx = x.pnts(p1)(0) - x.pnts(p0)(0);
        dy = x.pnts(p1)(1) - x.pnts(p0)(1);
        ol = 2./(dx*dx +dy*dy);
        psi = ol*((xpt(0) -x.pnts(p0)(0))*dx +(xpt(1) -x.pnts(p0)(1))*dy) -1.;
        normdist = dx*(xpt(1)-x.pnts(p0)(1))-dy*(xpt(0)-x.pnts(p1)(0));
        normdist *= sqrt(ol/2.);
        psiprev = psi;
        normdistprev = normdist;
        sidlocprev = nseg-1;
    } 
    else {
        psiprev = -1.0;
    }
    
    for(k=0;k<nseg;++k) {
        sind = seg(k);
        p0 = x.seg(sind).pnt(0);
        p1 = x.seg(sind).pnt(1);
        dx = x.pnts(p1)(0) - x.pnts(p0)(0);
        dy = x.pnts(p1)(1) - x.pnts(p0)(1);
        ol = 2./(dx*dx +dy*dy);
        psi = ol*((xpt(0) -x.pnts(p0)(0))*dx +(xpt(1) -x.pnts(p0)(1))*dy) -1.;
        normdist = dx*(xpt(1)-x.pnts(p0)(1))-dy*(xpt(0)-x.pnts(p1)(0));
        normdist *= sqrt(ol/2.);
        
        if (psi <= -1.0 && psiprev >= 1.0) {
            /* PREVIOUS & THIS SIDE ARE POTENTIAL MATCHES */
            if (fabs(normdist) < mindist) {
                mindist = fabs(normdist);
                sidloc = k;
                psiloc = -1.0;
            }
            if (fabs(normdistprev) < mindist) {
                mindist = fabs(normdistprev);
                sidloc = sidlocprev;
                psiloc = 1.0;
            }
        }
        else if (psi >= -1.0 && psi <= 1.0) {
            /* POTENTIAL SIDE */
            if (fabs(normdist) < mindist) {
                mindist = fabs(normdist);
                sidloc = k;
                psiloc = psi;
            }
        }
        psiprev = psi;
        normdistprev = normdist;
        sidlocprev = k;
    }
    
    return;
}



void edge_bdry::mgconnect(Array<tet_mesh::transfer,1> &cnnct,tet_mesh& tgt, int bnum) {
    int j,k,sind,tind,p0,sidloc;
    FLT psiloc;
        
    for(k=1;k<nseg;++k) {
        p0 = x.seg(seg(k)).pnt(0);
        tgt.ebdry(bnum)->findbdrypt(x.pnts(p0), sidloc, psiloc);
        sind = tgt.ebdry(bnum)->seg(sidloc);
        tind = tgt.seg(sind).tet;                            
        cnnct(p0).tet = tind;
        for (j=0;j<3;++j) 
            if (tgt.tri(tind).seg(j) == sind) break;
        assert(j < 3);
        cnnct(p0).wt(j) = 0.0;
        cnnct(p0).wt((j+1)%3) = 0.5*(1.-psiloc);
        cnnct(p0).wt((j+2)%3) = 0.5*(1.+psiloc);
    } 
    
    return;
}

/* SWAP ELEMENTS IN LIST */
void edge_bdry::swap(int s1, int s2) {
    int ind;
        
    /* FIXME NOT SURE HOW TO SWAP S VALUES */    
    ind = seg(s1);
    seg(s1) = seg(s2);
    seg(s2) = ind;
    
    /* PREV/NEXT IMPLEMENTATION */
    ind = prev(s1);
    prev(s1) = prev(s2);
    prev(s2) = ind;
    
    ind = next(s1);
    next(s1) = next(s2);
    next(s2) = ind;
    
    int prev1 = prev(s1);
    int next1 = next(s1);
    int prev2 = prev(s2);
    int next2 = next(s2);
    
    if (prev1 == s1) prev1 = s2;
    if (next1 == s1) next1 = s2;
    if (prev2 == s2) prev2 = s1;
    if (next2 == s2) next2 = s1;
    
    next(prev1) = s1;
    prev(next1) = s1;
    
    next(prev2) = s2;
    prev(next2) = s2;

    return;
}


/* REORDERS BOUNDARIES TO BE SEQUENTIAL */
/* USES gbl->i1wk & gbl->i2wk AS WORK ARRAYS */
void edge_bdry::reorder() {
    int i,count,total,sind,minp,first;

    total = nseg;
    
    /* DON'T ASSUME wk INITIALIZED TO -1 */
    for(i=0;i<nseg;++i) {
        sind = seg(i);
        x.gbl->i1wk(x.seg(sind).pnt(0)) = -1;
        x.gbl->i2wk(x.seg(sind).pnt(1)) = -1;
    }
    
    /* STORE SIDE INDICES BY VERTEX NUMBER */
    for(i=0; i < nseg; ++i) {
        sind = seg(i);
        x.gbl->i1wk(x.seg(sind).pnt(1)) = i;
        x.gbl->i2wk(x.seg(sind).pnt(0)) = i;
    }

    /* FIND FIRST SIDE */    
    first = -1;
    for(i=0;i<nseg;++i) {
        sind = seg(i);
        if (x.gbl->i1wk(x.seg(sind).pnt(0)) == -1) {
            first = i;
            break;
        }
    }
    
    /* SPECIAL CONSTRAINT IF LOOP */
    /* THIS IS TO ELIMINATE ANY INDEFINITENESS ABOUT SIDE ORDERING FOR LOOP */
    if (first < 0) {
        minp = x.npnt;
        for(i=0;i<nseg;++i) {
            sind = seg(i);
            if (x.seg(sind).pnt(0) < minp) {
                first = i;
                minp = x.seg(sind).pnt(0);
            }
        }
    }
    
    /* SWAP FIRST SIDE */
    count = 0;
    swap(count,first);
    x.gbl->i1wk(x.seg(seg(first)).pnt(1)) = first;
    x.gbl->i2wk(x.seg(seg(first)).pnt(0)) = first;
    x.gbl->i1wk(x.seg(seg(count)).pnt(1)) = count;
    x.gbl->i2wk(x.seg(seg(count)).pnt(0)) = -1;  // TO MAKE SURE LOOP STOPS

    /* REORDER LIST */
    while ((first = x.gbl->i2wk(x.seg(seg(count++)).pnt(1))) >= 0) {
        swap(count,first);
        x.gbl->i1wk(x.seg(seg(first)).pnt(1)) = first;
        x.gbl->i2wk(x.seg(seg(first)).pnt(0)) = first;
        x.gbl->i1wk(x.seg(seg(count)).pnt(1)) = count;
        x.gbl->i2wk(x.seg(seg(count)).pnt(0)) = count;
    }
    
    /* RESET gbl->i1wk TO -1 */
    for(i=0; i <total; ++i) {
        sind = seg(i);
        x.gbl->i1wk(x.seg(sind).pnt(1)) = -1;
    }
    
    if (count < total) {
        ++x.nebd;
        x.ebdry.resizeAndPreserve(x.nebd);
        x.ebdry(x.nebd-1) = create(x);
        x.ebdry(x.nebd-1)->copy(*this);
        nseg = count;

        for(i=0;i<total-nseg;++i)
            x.ebdry(x.nebd-1)->swap(i,i+nseg);
        x.ebdry(x.nebd-1)->nseg = total-nseg;
        *x.gbl->log << "#creating new boundary: " << idnum << " num: " << x.ebdry(x.nebd-1)->nseg << std::endl;
    }
    
    for (i=0;i<nseg;++i) {
        prev(i) = i-1;
        next(i) = i+1;
    }
    next(nseg-1) = -1;
    
    return;
}

void ecomm::ploadbuff(boundary::groups grp,FLT *base,int bgn,int end, int stride) {
    int j,k,count,sind,offset;
    
    if (!((1<<grp)&groupmask)) return;

    count = 0;
    for(j=0;j<nseg;++j) {
        sind = seg(j);
        offset = x.seg(sind).pnt(0)*stride;
        for (k=bgn;k<=end;++k) {
            fsndbuf(count++) = base[offset+k];
        }
    }
    offset = x.seg(sind).pnt(1)*stride;
    for (k=bgn;k<=end;++k) 
        fsndbuf(count++) = base[offset+k]; 
        
    sndsize() = count;
    sndtype() = boundary::flt_msg;
}

void ecomm::pfinalrcv(boundary::groups grp, int phi, comm_type type, operation op, FLT *base,int bgn,int end, int stride) {
    int j,k,m,count,countdn,countup,offset,sind;
    int matches = 1;
    FLT mtchinv;
    /* ASSUMES REVERSE ORDERING OF SIDES */
    /* WON'T WORK IN 3D */
        
    if (!((1<<grp)&groupmask)) return;
    
    switch(type) {
        case(slave_master): {
            if (!first) return;
        }
        
        case(master_slave): {
            if (first || phase(grp)(0) != phi) return;
            
#ifdef MPDEBUG
            *x.gbl->log << "finalrcv"  << idnum << " " << is_frst() << std::endl;
#endif
            int ebp1 = end-bgn+1;
            countdn = nseg*ebp1;
            for(j=0;j<nseg;++j) {
                sind = seg(j);
                offset = x.seg(sind).pnt(0)*stride +bgn;
                for(k=0;k<ebp1;++k) {
                    base[offset+k] = frcvbuf(0,countdn +k);
#ifdef MPDEBUG
                    *x.gbl->log << "\t" << base[offset+k] << std::endl;
#endif
                }
                countdn -= ebp1;
            }
            offset = x.seg(sind).pnt(1)*stride +bgn;
            for(k=0;k<ebp1;++k) {
                base[offset+k] = frcvbuf(0,countdn+k);
#ifdef MPDEBUG
                *x.gbl->log << "\t" << base[offset+k] << std::endl;
#endif
            }
            return;
        }
            
        default: {
            switch(op) {
                case(average):
    
                    /* RELOAD FROM BUFFER */
                    /* ELIMINATES V/S/F COUPLING IN ONE PHASE */
                    /* FINALRCV SHOULD BE CALLED F,S,V ORDER (V HAS FINAL AUTHORITY) */    
                    for(m=0;m<nmatch;++m) {    
                        if (phase(grp)(m) != phi) continue;
                        
                        ++matches;
                        
                        int ebp1 = end-bgn+1;
                        countdn = nseg*ebp1;
                        countup = 0;
                        for(j=0;j<nseg+1;++j) {
                            for(k=0;k<ebp1;++k)
                                fsndbuf(countup +k) += frcvbuf(m,countdn +k);
                            countup += ebp1;
                            countdn -= ebp1;
                        }
                    }

                    if (matches > 1) {
                        mtchinv = 1./matches;

#ifdef MPDEBUG
                        *x.gbl->log << "finalrcv"  << idnum << " " << is_frst() << std::endl;
#endif
                        count = 0;
                        for(j=0;j<nseg;++j) {
                            sind = seg(j);
                            offset = x.seg(sind).pnt(0)*stride;
                            for (k=bgn;k<=end;++k) {
                                base[offset+k] = fsndbuf(count++)*mtchinv;
#ifdef MPDEBUG
                                *x.gbl->log << "\t" << base[offset+k] << std::endl;
#endif
                            }

                        }
                        offset = x.seg(sind).pnt(1)*stride;
                        for (k=bgn;k<=end;++k) {
                            base[offset+k] = fsndbuf(count++)*mtchinv;
#ifdef MPDEBUG
                            *x.gbl->log << "\t" << base[offset+k] << std::endl;
#endif
                        }
                    }
                    return;
                case(sum): 
                     matches = 1;
    
                    /* RELOAD FROM BUFFER */
                    /* ELIMINATES V/S/F COUPLING IN ONE PHASE */
                    /* FINALRCV SHOULD BE CALLED F,S,V ORDER (V HAS FINAL AUTHORITY) */    
                    for(m=0;m<nmatch;++m) {    
                        if (phase(grp)(m) != phi) continue;
                        
                        ++matches;
                        
                        int ebp1 = end-bgn+1;
                        countdn = nseg*ebp1;
                        countup = 0;
                        for(j=0;j<nseg+1;++j) {
                            for(k=0;k<ebp1;++k)
                                fsndbuf(countup +k) += frcvbuf(m,countdn +k);
                            countup += ebp1;
                            countdn -= ebp1;
                        }
                    }

                    if (matches > 1) {
#ifdef MPDEBUG
                        *x.gbl->log << "finalrcv"  << idnum << " " << is_frst() << std::endl;
#endif
                        count = 0;
                        for(j=0;j<nseg;++j) {
                            sind = seg(j);
                            offset = x.seg(sind).pnt(0)*stride;
                            for (k=bgn;k<=end;++k) {
                                base[offset+k] = fsndbuf(count++);
#ifdef MPDEBUG
                                *x.gbl->log << "\t" << base[offset+k] << std::endl;
#endif
                            }

                        }
                        offset = x.seg(sind).pnt(1)*stride;
                        for (k=bgn;k<=end;++k) {
                            base[offset+k] = fsndbuf(count++);
#ifdef MPDEBUG
                            *x.gbl->log << "\t" << base[offset+k] << std::endl;
#endif
                        }
                    }
                    return;
                
                case(maximum):                    
                    matches = 1;
    
                    /* RELOAD FROM BUFFER */
                    /* ELIMINATES V/S/F COUPLING IN ONE PHASE */
                    /* FINALRCV SHOULD BE CALLED F,S,V ORDER (V HAS FINAL AUTHORITY) */    
                    for(m=0;m<nmatch;++m) {    
                        if (phase(grp)(m) != phi) continue;
                        
                        ++matches;
                        
                        int ebp1 = end-bgn+1;
                        countdn = nseg*ebp1;
                        countup = 0;
                        for(j=0;j<nseg+1;++j) {
                            for(k=0;k<ebp1;++k)
                                fsndbuf(countup +k) = MAX(fsndbuf(countup+k),frcvbuf(m,countdn +k));
                            countup += ebp1;
                            countdn -= ebp1;
                        }
                    }

                    if (matches > 1) {
#ifdef MPDEBUG
                        *x.gbl->log << "finalrcv"  << idnum << " " << is_frst() << std::endl;
#endif
                        count = 0;
                        for(j=0;j<nseg;++j) {
                            sind = seg(j);
                            offset = x.seg(sind).pnt(0)*stride;
                            for (k=bgn;k<=end;++k) {
                                base[offset+k] = fsndbuf(count++);
#ifdef MPDEBUG
                                *x.gbl->log << "\t" << base[offset+k] << std::endl;
#endif
                            }

                        }
                        offset = x.seg(sind).pnt(1)*stride;
                        for (k=bgn;k<=end;++k) {
                            base[offset+k] = fsndbuf(count++);
#ifdef MPDEBUG
                            *x.gbl->log << "\t" << base[offset+k] << std::endl;
#endif
                        }
                    }
                    return;
                
                default: 
                    *x.gbl->log << "replacement with symmetric sending?" << std::endl;
                    exit(1);
            }
            break;
        }
    }
}
    
void ecomm::sloadbuff(boundary::groups grp,FLT *base,int bgn,int end, int stride) {
    int j,k,count,sind,offset;
    
    if (!((1<<grp)&groupmask)) return;

    count = 0;
    for(j=0;j<nseg;++j) {
        sind = seg(j);
        offset = sind*stride;
        for (k=bgn;k<=end;++k) {
            fsndbuf(count++) = base[offset+k];
        }
    }
        
    sndsize() = count;
    sndtype() = boundary::flt_msg;
}

void ecomm::sfinalrcv(boundary::groups grp, int phi, comm_type type, operation op, FLT *base,int bgn,int end, int stride) {
    int j,k,count,offset,sind;
    /* ASSUMES REVERSE ORDERING OF SIDES */
    /* WON'T WORK IN 3D */
    
    bool reload = comm_finish(grp,phi,type,op);
    if (!reload) return;
   
     count = 0;
    for(j=0;j<nseg;++j) {
        sind = seg(j);
        offset = sind*stride;
        for (k=bgn;k<=end;++k) {
            base[offset+k] = fsndbuf(count++);
        }
    }
    return;
}

void epartition::mgconnect(Array<tet_mesh::transfer,1> &cnnct,tet_mesh& tgt, int bnum) {
    int i,j,k,p0;
    
 
    /* BOUNDARY IS AN INTERNAL PARTITION BOUNDARY */
    /* MAKE SURE ENDPOINTS ARE OK */
    i = x.seg(seg(0)).pnt(0);
    if (cnnct(i).tet < 0) {
        tgt.otree.nearpt(x.pnts(i).data(),p0);
        cnnct(i).tet=tgt.pnt(p0).tet;
        for(j=0;j<4;++j) {
            cnnct(i).wt(j) = 0.0;
            if (tgt.tet(cnnct(i).tet).pnt(j) == p0) cnnct(i).wt(j) = 1.0;
        }
    }
    i = x.seg(seg(nseg-1)).pnt(1);
    if (cnnct(i).tet < 0) {
        tgt.otree.nearpt(x.pnts(i).data(),p0);
        cnnct(i).tet=tgt.pnt(p0).tet;
        for(j=0;j<3;++j) {
            cnnct(i).wt(j) = 0.0;
            if (tgt.tet(cnnct(i).tet).pnt(j) == p0) cnnct(i).wt(j) = 1.0;
        }
    }
    
    if (first) {
        sndsize() = 0;
        sndtype() = int_msg;
        for(k=1;k<nseg;++k) {
            p0 = x.seg(seg(k)).pnt(0);
            if (cnnct(p0).tet > 0) {
                isndbuf(sndsize()++) = -1;
            }
            else {
                isndbuf(sndsize()++) = +1;
                cnnct(p0).tet = 0;
                for(j=0;j<3;++j)
                    cnnct(p0).wt(j) = 0.0;
            }
        }
    }
    
    comm_prepare(boundary::all,0,slave_master); 
    comm_exchange(boundary::all,0,slave_master);
    comm_wait(boundary::all,0,slave_master);
    
    if (!first) {
        i = 0;
        for(k=nseg-1;k>0;--k) {
            p0 = x.seg(seg(k)).pnt(1);
            if (ircvbuf(0,i) < 0) {
                cnnct(p0).tet = 0;
                for(j=0;j<3;++j)
                    cnnct(p0).wt(j) = 0.0;
            }
        }
    }                    
}

/********************/
/* FACE FUNCTIONS */
/********************/
void face_bdry::alloc(int mxsize) {
     
     /* SIDE INFO */
     maxpst = mxsize;
     seg.resize(Range(-1,maxpst));

     /* VERTEX INFO */                                
     pnt.resize(Range(-1,maxpst));
     
     /* TRI INFO */ 
     tri.resize(Range(-1,maxpst));
          
     return;
}

void face_bdry::copy(const face_bdry& tgt) {
    int i;
        
    if (!maxpst) {
        alloc(tgt.maxpst);
    }
    else {
        /* CHECK IF BIG ENOUGH */
        if (tgt.nseg > maxpst) {
            *x.gbl->log << "face_bdry is too big to copy" << std::endl;
            exit(1);
        }
    }
    
    /* COPY VERTEX INFO OVER */
    npnt = tgt.npnt;
    for(i=0;i<npnt;++i)
        pnt(i) = tgt.pnt(i);
            
    /* COPY SIDE INFORMATION */
    nseg = tgt.nseg;
    for(i=0;i<nseg;++i) {
        seg(i) = tgt.seg(i);
    }

    nebd = tgt.nebd;
    ebdry = tgt.ebdry;
    
    /* COPY ELEMENT DATA */
    ntri = tgt.ntri;
    for(i=0;i<ntri;++i) {
        tri(i) = tgt.tri(i);
    }
        
    otree.copy(tgt.otree);
    otree.change_vptr((FLT (*)[tet_mesh::ND]) tgt.x.pnts(0).data() );
    
    return;  
}

void face_bdry::mgconnect(Array<tet_mesh::transfer,1> &cnnct,tet_mesh& tgt, int bnum) {
    int j,k,tind,ttind,p0,p1,facloc;
    FLT r, s;
    
    // vertices of each face on a tet
    TinyMatrix<int,4,3> vf1;
    vf1(0,0)=1, vf1(0,1)=2, vf1(0,2)=3;
    vf1(1,0)=0, vf1(1,1)=3, vf1(1,2)=2;
    vf1(2,0)=0, vf1(2,1)=1, vf1(2,2)=3;
    vf1(3,0)=0, vf1(3,1)=2, vf1(3,2)=1;
        
    for(j=0;j<npnt;++j) {
        p0 = pnt(j).gindx;

        FLT dist = tgt.otree.nearpt(x.pnts(p0).data(),p1);  //FIXME  
        if (dist > EPSILON*10.) {
            *x.gbl->log << "findbdrypt is not general enough for arbitrary points yet\n";
            exit(1);  // FIXME
        }
        ttind = tgt.pnt(p0).tet;
        cnnct(p0).tet = ttind;
        for (k=0;k<4;++k) {
            if (tgt.tet(ttind).pnt(k) == p1) 
                cnnct(p0).wt(k) = 1.0;
            else 
                cnnct(p0).wt(k) = 0.0;
        }
        
         /* THIS IS THE RIGHT WAY, BUT FINDBDRYPT ISN'T WORKING YET */
//        tgt.fbdry(bnum)->findbdrypt(x.pnts(p0), facloc, r, s);
//        tind = tgt.fbdry(bnum)->tri(facloc).gindx;
//        ttind = tgt.tri(tind).tet(0);                             
//        cnnct(p0).tet = ttind;
//        for (k=0;k<4;++k) 
//            if (tgt.tet(ttind).tri(k) == tind) break;
//        assert(k < 3);
//        cnnct(p0).wt(k) = 0.0;
//        
//        // Assumes Definition of face is compatible with definition on tet
//        cnnct(p0).wt(vf1(k,0)) = 0.5*(1.+s);
//        cnnct(p0).wt(vf1(k,1)) = 0.5*(-r -s);
//        cnnct(p0).wt(vf1(k,2)) = 0.5*(1.+r);
    } 
    
    return;
}

void face_bdry::findbdrypt(const TinyVector<FLT,tet_mesh::ND> xpt, int &facloc, FLT &r, FLT &s) const {
    if (1) {
        *x.gbl->log << "findbdrypt is not general enough for arbitrary points yet\n";
        exit(1);  // FIXME
    }

    return;
}

void face_bdry::mvpttobdry(int indx, FLT r, FLT s, TinyVector<FLT,tet_mesh::ND> &pt) {
    /* FOR A LINEAR SIDE */
    int n;
    
    for (n=0;n<tet_mesh::ND;++n)
        pt(n) = 0.5*((1.+s)*x.pnts(x.tri(tri(indx).gindx).pnt(0))(n) +(-r-s)*x.pnts(x.tri(tri(indx).gindx).pnt(1))(n) +(1.+r)*x.pnts(x.tri(tri(indx).gindx).pnt(2))(n));
    
    return;
}

void fcomm::ploadbuff(boundary::groups grp,FLT *base,int bgn,int end, int stride) {
    int j,k,count,offset;
    
    if (!in_group(grp)) return;

    count = 0;
    for(j=0;j<npnt;++j) {
        offset = pnt(j).gindx*stride;
        for (k=bgn;k<=end;++k) {
            fsndbuf(count++) = base[offset+k];
        }
    }
        
    sndsize() = count;
    sndtype() = boundary::flt_msg;
}

void fcomm::pfinalrcv(boundary::groups grp, int phi, comm_type type, operation op, FLT *base,int bgn,int end, int stride) {
    int j,k,count,offset;
    /* ASSUMES REVERSE ORDERING OF SIDES */
    /* WON'T WORK IN 3D */
        
    bool reload = comm_finish(grp,phi,type,op);
    if (!reload) return;
  
    count = 0;
    for(j=0;j<npnt;++j) {
        offset = pnt(j).gindx*stride;
        for (k=bgn;k<=end;++k) {
            base[offset+k] = fsndbuf(count++);
        }
    }  
    return;
}

void fcomm::sloadbuff(boundary::groups grp,FLT *base,int bgn,int end, int stride) {
    int j,k,count,sind,offset;
    
    if (!in_group(grp)) return;

    count = 0;
    for(j=0;j<nseg;++j) {
        sind = seg(j).gindx;
        offset = sind*stride;
        for (k=bgn;k<=end;++k) {
            fsndbuf(count++) = base[offset+k];
        }
    }
        
    sndsize() = count;
    sndtype() = boundary::flt_msg;
}

void fcomm::sfinalrcv(boundary::groups grp, int phi, comm_type type, operation op, FLT *base,int bgn,int end, int stride) {
    int j,k,count,offset,sind;
    
    bool reload = comm_finish(grp,phi,type,op);
    if (!reload) return;
    
    count = 0;
    for(j=0;j<nseg;++j) {
        sind = seg(j).gindx;
        offset = sind*stride;
        for (k=bgn;k<=end;++k) {
            base[offset+k] = fsndbuf(count++);
        }

    }
    return;
}

void fcomm::tloadbuff(boundary::groups grp,FLT *base,int bgn,int end, int stride) {
    int j,k,count,tind,offset;
    
    if (!in_group(grp)) return;

    count = 0;
    for(j=0;j<ntri;++j) {
        tind = tri(j).gindx;
        offset = tind*stride;
        for (k=bgn;k<=end;++k) {
            fsndbuf(count++) = base[offset+k];
        }
    }
        
    sndsize() = count;
    sndtype() = boundary::flt_msg;
}

void fcomm::tfinalrcv(boundary::groups grp, int phi, comm_type type, operation op, FLT *base,int bgn,int end, int stride) {
    int j,k,count,offset,tind;
    /* ASSUMES REVERSE ORDERING OF SIDES */
    /* WON'T WORK IN 3D */
    
    bool reload = comm_finish(grp,phi,type,op);
    if (!reload) return;
 
    count = 0;
    for(j=0;j<ntri;++j) {
        tind = tri(j).gindx;
        offset = tind*stride;
        for (k=bgn;k<=end;++k) {
             base[offset+k] = fsndbuf(count++);
        }
    }
    return;
}

void fcomm::match_numbering(int step) {
    switch (step) {
        case(1): {
            if (is_frst()) {
                return;
            }
            else {
                /* Slave receives master boundaries vertex positions */
                int count = 0;
                for (int i=0;i<npnt;++i) {
                    TinyVector<FLT,tet_mesh::ND> mpnt;
                    for (int n=0;n<tet_mesh::ND;++n)
                        mpnt(n) = fsndbuf(count++);
                        
                    FLT dist = x.otree.nearpt(mpnt.data(),pnt(i).gindx);
                    if (dist > 10.*EPSILON) {
                        *x.gbl->log << "Matching face numbering error: " << dist << ' ' << mpnt << '\n';
                    }
                }
            }
            return;
        }
        
        case(2): {
            if (is_frst()) {
             /* Master loads integer data */
                int count = 0;
                for (int i=0;i<nseg;++i) {
                    isndbuf(count++) = seg(i).pnt(0);
                    isndbuf(count++) = seg(i).pnt(1);
                }
                
                for (int i=0;i<ntri;++i) {
                    for(int j=0;j<3;++j) {
                        isndbuf(count++) = tri(i).pnt(j);
                    }
                }
                sndsize() = count;
                sndtype() = boundary::int_msg;
            }
            else {
                sndsize() = 2*nseg +3*ntri;
                sndtype() = boundary::int_msg;
            }
            return;
        }
        
        case(3): {
            if (is_frst()) return;
            
            int count = 0;
            for (int i=0;i<nseg;++i) {
                seg(i).pnt(0) = isndbuf(count++);
                seg(i).pnt(1) = isndbuf(count++);
                *x.gbl->log << seg(i).pnt(0) << ' ' << seg(i).pnt(1) << std::endl;
            }
            gbltolclside();
            
            /* Triangles are reverse oriented */
            for (int i=0;i<ntri;++i) {
                tri(i).pnt(0) = isndbuf(count++);
                tri(i).pnt(2) = isndbuf(count++);
                tri(i).pnt(1) = isndbuf(count++);
            }
            gbltolcltri();
        }
    }
    
    return;
}
        
            
            




