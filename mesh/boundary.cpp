#include "boundaries.h"
#include <assert.h>
#include <float.h>

/********************/
/* VERTEX FUNCTIONS */
/********************/
    
/* GENERIC VERTEX COMMUNICATIONS */
void vcomm::vloadbuff(boundary::groups grp,FLT *base,int bgn,int end, int stride) {
    int i,offset;

    if (!((1<<grp)&groupmask)) return;

    sndsize()=end-bgn+1;
    sndtype()=flt_msg;
    
    /* LOAD SEND BUFFER */    
    offset = v0*stride +bgn;
    for (i=0;i<end-bgn+1;++i) 
        fsndbuf(i) = base[offset+i];
}

void vcomm::vfinalrcv(boundary::groups grp, int phi, comm_type type, operation op, FLT *base, int bgn, int end, int stride) {
    int i,m,offset;
    int matches = 1;
    
    
    if (!((1<<grp)&groupmask)) return;

    switch(type) {
        case(slave_master): {
            if (!first) return;
        }
        case(master_slave): {
            if (first || phase(grp)(0) != phi) return;
            
            offset = v0*stride +bgn;
            for(i=0;i<end-bgn+1;++i) 
                base[offset++] = frcvbuf(0,i);
                
            return;
        }
            
        default: {
            switch(op) {
                case(average):
                    for(m=0;m<nmatch;++m) {
                        if (phase(grp)(m) != phi) continue;
                        ++matches;
                        
                        for(i=0;i<end-bgn+1;++i) 
                            fsndbuf(i) += frcvbuf(m,i);
                    }
                    
                    if (matches > 1) {
                        offset = v0*stride +bgn;
                        for(i=0;i<end-bgn+1;++i) 
                            base[offset++] = fsndbuf(i)/matches;
                    }
                    return;
                case(sum): 
                    for(m=0;m<nmatch;++m) {
                        if (phase(grp)(m) != phi) continue;
                        ++matches;
                                                
                        for(i=0;i<end-bgn+1;++i) 
                            fsndbuf(i) += frcvbuf(m,i);
                    }
                    
                    if (matches > 1) {
                        offset = v0*stride +bgn;
                        for(i=0;i<end-bgn+1;++i) 
                            base[offset++] = fsndbuf(i);
                    }
                    return;
                case(maximum): 
                    for(m=0;m<nmatch;++m) {
                        if (phase(grp)(m) != phi) continue;
                        ++matches;
                                                
                        for(i=0;i<end-bgn+1;++i) 
                            fsndbuf(i) = MAX(fsndbuf(i),frcvbuf(m,i));
                    }
                    
                    if (matches > 1) {
                        offset = v0*stride +bgn;
                        for(i=0;i<end-bgn+1;++i) 
                            base[offset++] = fsndbuf(i);
                    }
                    return;
                default: 
                    *sim::log << "replacement with symmetric sending?" << std::endl;
                    exit(1);
            }
            break;
        }
    }
}

/**************************************/
/* GENERIC FUNCTIONS FOR SIDES          */
/**************************************/
void side_bdry::alloc(int n) {
    maxel = n;
    el.resize(n);
}
      
void side_bdry::copy(const side_bdry& bin) {
    int i;
    
        
    if (!maxel) alloc(bin.maxel);
	else assert(bin.nel <= maxel);
    vbdry = bin.vbdry;
    
    nel = bin.nel;
    
    for(i=0;i<nel;++i)
        el(i) = bin.el(i);
        
    return;
}

void side_bdry::mvpttobdry(int indx, FLT psi, TinyVector<FLT,mesh::ND> &pt) {
    /* FOR A LINEAR SIDE */
    int n;
    
    for (n=0;n<mesh::ND;++n)
        pt(n) = 0.5*((1. -psi)*x.vrtx(x.sd(el(indx)).vrtx(0))(n) +(1.+psi)*x.vrtx(x.sd(el(indx)).vrtx(1))(n));
    
    return;
}

void side_bdry::findbdrypt(const TinyVector<FLT,2> xpt, int &sidloc, FLT &psiloc) const {
    int k,sind,v0,v1,sidlocprev;
    FLT dx,dy,ol,psi,normdist;
    FLT psiprev,normdistprev;
    FLT mindist = 1.0e32;
        
    if (x.sd(el(0)).vrtx(0) == x.sd(el(nel-1)).vrtx(1)) {
        /* BOUNDARY IS A LOOP */
        sind = el(nel-1);
        v0 = x.sd(sind).vrtx(0);
        v1 = x.sd(sind).vrtx(1);
        dx = x.vrtx(v1)(0) - x.vrtx(v0)(0);
        dy = x.vrtx(v1)(1) - x.vrtx(v0)(1);
        ol = 2./(dx*dx +dy*dy);
        psi = ol*((xpt(0) -x.vrtx(v0)(0))*dx +(xpt(1) -x.vrtx(v0)(1))*dy) -1.;
        normdist = dx*(xpt(1)-x.vrtx(v0)(1))-dy*(xpt(0)-x.vrtx(v1)(0));
        normdist *= sqrt(ol/2.);
        psiprev = psi;
        normdistprev = normdist;
        sidlocprev = nel-1;
    } 
    else {
        psiprev = -1.0;
    }
    
    for(k=0;k<nel;++k) {
        sind = el(k);
        v0 = x.sd(sind).vrtx(0);
        v1 = x.sd(sind).vrtx(1);
        dx = x.vrtx(v1)(0) - x.vrtx(v0)(0);
        dy = x.vrtx(v1)(1) - x.vrtx(v0)(1);
        ol = 2./(dx*dx +dy*dy);
        psi = ol*((xpt(0) -x.vrtx(v0)(0))*dx +(xpt(1) -x.vrtx(v0)(1))*dy) -1.;
        normdist = dx*(xpt(1)-x.vrtx(v0)(1))-dy*(xpt(0)-x.vrtx(v1)(0));
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



void side_bdry::mgconnect(Array<mesh::transfer,1> &cnnct, mesh& tgt, int bnum) {
    int j,k,sind,tind,v0,sidloc;
    FLT psiloc;
        
    for(k=1;k<nel;++k) {
        v0 = x.sd(el(k)).vrtx(0);
        tgt.sbdry(bnum)->findbdrypt(x.vrtx(v0), sidloc, psiloc);
        sind = tgt.sbdry(bnum)->el(sidloc);
        tind = tgt.sd(sind).tri(0);                            
        cnnct(v0).tri = tind;
        for (j=0;j<3;++j) 
            if (tgt.td(tind).side(j) == sind) break;
        assert(j < 3);
        cnnct(v0).wt(j) = 0.0;
        cnnct(v0).wt((j+1)%3) = 0.5*(1.-psiloc);
        cnnct(v0).wt((j+2)%3) = 0.5*(1.+psiloc);
    } 
    
    return;
}

/* SWAP ELEMENTS IN LIST */
void side_bdry::swap(int s1, int s2) {
    int ind;
    
    /* TEMPORARY NOT SURE HOW TO SWAP S VALUES */    
    ind = el(s1);
    el(s1) = el(s2);
    el(s2) = ind;

    return;
}


/* REORDERS BOUNDARIES TO BE SEQUENTIAL */
/* USES gbl_ptr->intwk & gbl_ptr->i2wk AS WORK ARRAYS */
void side_bdry::reorder() {
    int i,count,total,sind,minv,first;

    total = nel;
    
    /* DON'T ASSUME wk INITIALIZED TO -1 */
    for(i=0;i<nel;++i) {
        sind = el(i);
        x.gbl_ptr->intwk(x.sd(sind).vrtx(0)) = -1;
        x.gbl_ptr->i2wk(x.sd(sind).vrtx(1)) = -1;
    }
    
    /* STORE SIDE INDICES BY VERTEX NUMBER */
    for(i=0; i < nel; ++i) {
        sind = el(i);
        x.gbl_ptr->intwk(x.sd(sind).vrtx(1)) = i;
        x.gbl_ptr->i2wk(x.sd(sind).vrtx(0)) = i;
    }

    /* FIND FIRST SIDE */    
    first = -1;
    for(i=0;i<nel;++i) {
        sind = el(i);
        if (x.gbl_ptr->intwk(x.sd(sind).vrtx(0)) == -1) {
            first = i;
            break;
        }
    }
    
    /* SPECIAL CONSTRAINT IF LOOP */
    /* THIS IS TO ELIMINATE ANY INDEFINITENESS ABOUT SIDE ORDERING FOR LOOP */
    if (first < 0) {
        minv = x.nvrtx;
        for(i=0;i<nel;++i) {
            sind = el(i);
            if (x.sd(sind).vrtx(1) < minv) {
                first = i;
                minv = x.sd(sind).vrtx(1);
            }
        }
    }
    
    /* SWAP FIRST SIDE */
    count = 0;
    swap(count,first);
    x.gbl_ptr->intwk(x.sd(el(first)).vrtx(1)) = first;
    x.gbl_ptr->i2wk(x.sd(el(first)).vrtx(0)) = first;
    x.gbl_ptr->intwk(x.sd(el(count)).vrtx(1)) = count;
    x.gbl_ptr->i2wk(x.sd(el(count)).vrtx(0)) = -1;  // TO MAKE SURE LOOP STOPS

    /* REORDER LIST */
    while ((first = x.gbl_ptr->i2wk(x.sd(el(count++)).vrtx(1))) >= 0) {
        swap(count,first);
        x.gbl_ptr->intwk(x.sd(el(first)).vrtx(1)) = first;
        x.gbl_ptr->i2wk(x.sd(el(first)).vrtx(0)) = first;
        x.gbl_ptr->intwk(x.sd(el(count)).vrtx(1)) = count;
        x.gbl_ptr->i2wk(x.sd(el(count)).vrtx(0)) = count;
    }
    
    /* RESET gbl_ptr->intwk TO -1 */
    for(i=0; i <total; ++i) {
        sind = el(i);
        x.gbl_ptr->intwk(x.sd(sind).vrtx(1)) = -1;
    }
    
    if (count < total) {
        ++x.nsbd;
        x.sbdry.resizeAndPreserve(x.nsbd);
        x.sbdry(x.nsbd-1) = create(x);
        x.sbdry(x.nsbd-1)->copy(*this);
        nel = count;

        for(i=0;i<total-nel;++i)
            x.sbdry(x.nsbd-1)->swap(i,i+nel);
        x.sbdry(x.nsbd-1)->nel = total-nel;
        *sim::log << "#creating new boundary: " << idnum << " num: " << x.sbdry(x.nsbd-1)->nel << std::endl;
        return;
    }
    
    return;
}

void scomm::vloadbuff(boundary::groups grp,FLT *base,int bgn,int end, int stride) {
    int j,k,count,sind,offset;
    
    if (!((1<<grp)&groupmask)) return;

    count = 0;
    for(j=0;j<nel;++j) {
        sind = el(j);
        offset = x.sd(sind).vrtx(0)*stride;
        for (k=bgn;k<=end;++k) {
            fsndbuf(count++) = base[offset+k];
        }
    }
    offset = x.sd(sind).vrtx(1)*stride;
    for (k=bgn;k<=end;++k) 
        fsndbuf(count++) = base[offset+k]; 
        
    sndsize() = count;
    sndtype() = boundary::flt_msg;
}

void scomm::vfinalrcv(boundary::groups grp, int phi, comm_type type, operation op, FLT *base,int bgn,int end, int stride) {
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
            *sim::log << "finalrcv"  << idnum << " " << is_frst() << std::endl;
#endif
            int ebp1 = end-bgn+1;
            countdn = nel*ebp1;
            for(j=0;j<nel;++j) {
                sind = el(j);
                offset = x.sd(sind).vrtx(0)*stride +bgn;
                for(k=0;k<ebp1;++k) {
                    base[offset+k] = frcvbuf(0,countdn +k);
#ifdef MPDEBUG
                    *sim::log << "\t" << base[offset+k] << std::endl;
#endif
                }
                countdn -= ebp1;
            }
            offset = x.sd(sind).vrtx(1)*stride +bgn;
            for(k=0;k<ebp1;++k) {
                base[offset+k] = frcvbuf(0,countdn+k);
#ifdef MPDEBUG
                *sim::log << "\t" << base[offset+k] << std::endl;
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
                        countdn = nel*ebp1;
                        countup = 0;
                        for(j=0;j<nel+1;++j) {
                            for(k=0;k<ebp1;++k)
                                fsndbuf(countup +k) += frcvbuf(m,countdn +k);
                            countup += ebp1;
                            countdn -= ebp1;
                        }
                    }

                    if (matches > 1) {
                        mtchinv = 1./matches;

#ifdef MPDEBUG
                        *sim::log << "finalrcv"  << idnum << " " << is_frst() << std::endl;
#endif
                        count = 0;
                        for(j=0;j<nel;++j) {
                            sind = el(j);
                            offset = x.sd(sind).vrtx(0)*stride;
                            for (k=bgn;k<=end;++k) {
                                base[offset+k] = fsndbuf(count++)*mtchinv;
#ifdef MPDEBUG
                                *sim::log << "\t" << base[offset+k] << std::endl;
#endif
                            }

                        }
                        offset = x.sd(sind).vrtx(1)*stride;
                        for (k=bgn;k<=end;++k) {
                            base[offset+k] = fsndbuf(count++)*mtchinv;
#ifdef MPDEBUG
                            *sim::log << "\t" << base[offset+k] << std::endl;
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
                        countdn = nel*ebp1;
                        countup = 0;
                        for(j=0;j<nel+1;++j) {
                            for(k=0;k<ebp1;++k)
                                fsndbuf(countup +k) += frcvbuf(m,countdn +k);
                            countup += ebp1;
                            countdn -= ebp1;
                        }
                    }

                    if (matches > 1) {
#ifdef MPDEBUG
                        *sim::log << "finalrcv"  << idnum << " " << is_frst() << std::endl;
#endif
                        count = 0;
                        for(j=0;j<nel;++j) {
                            sind = el(j);
                            offset = x.sd(sind).vrtx(0)*stride;
                            for (k=bgn;k<=end;++k) {
                                base[offset+k] = fsndbuf(count++);
#ifdef MPDEBUG
                                *sim::log << "\t" << base[offset+k] << std::endl;
#endif
                            }

                        }
                        offset = x.sd(sind).vrtx(1)*stride;
                        for (k=bgn;k<=end;++k) {
                            base[offset+k] = fsndbuf(count++);
#ifdef MPDEBUG
                            *sim::log << "\t" << base[offset+k] << std::endl;
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
                        countdn = nel*ebp1;
                        countup = 0;
                        for(j=0;j<nel+1;++j) {
                            for(k=0;k<ebp1;++k)
                                fsndbuf(countup +k) = MAX(fsndbuf(countup+k),frcvbuf(m,countdn +k));
                            countup += ebp1;
                            countdn -= ebp1;
                        }
                    }

                    if (matches > 1) {
#ifdef MPDEBUG
                        *sim::log << "finalrcv"  << idnum << " " << is_frst() << std::endl;
#endif
                        count = 0;
                        for(j=0;j<nel;++j) {
                            sind = el(j);
                            offset = x.sd(sind).vrtx(0)*stride;
                            for (k=bgn;k<=end;++k) {
                                base[offset+k] = fsndbuf(count++);
#ifdef MPDEBUG
                                *sim::log << "\t" << base[offset+k] << std::endl;
#endif
                            }

                        }
                        offset = x.sd(sind).vrtx(1)*stride;
                        for (k=bgn;k<=end;++k) {
                            base[offset+k] = fsndbuf(count++);
#ifdef MPDEBUG
                            *sim::log << "\t" << base[offset+k] << std::endl;
#endif
                        }
                    }
                    return;
                
                default: 
                    *sim::log << "replacement with symmetric sending?" << std::endl;
                    exit(1);
            }
            break;
        }
    }
}
    
void scomm::sloadbuff(boundary::groups grp,FLT *base,int bgn,int end, int stride) {
    int j,k,count,sind,offset;
    
    if (!((1<<grp)&groupmask)) return;

    count = 0;
    for(j=0;j<nel;++j) {
        sind = el(j);
        offset = sind*stride;
        for (k=bgn;k<=end;++k) {
            fsndbuf(count++) = base[offset+k];
        }
    }
        
    sndsize() = count;
    sndtype() = boundary::flt_msg;
}

void scomm::sfinalrcv(boundary::groups grp, int phi, comm_type type, operation op, FLT *base,int bgn,int end, int stride) {
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
            *sim::log << "finalrcv"  << idnum << " " << is_frst() << std::endl;
#endif
            int ebp1 = end-bgn+1;
            countdn = (nel-1)*ebp1;
            countup = 0;
            for(j=0;j<nel;++j) {
                sind = el(j);
                offset = sind*stride +bgn;
                for (k=0;k<ebp1;++k) {
                    base[offset+k] = frcvbuf(0,countdn +k);
#ifdef MPDEBUG
                    *sim::log << "\t" << base[offset+k] << std::endl;
#endif
                }
                countdn -= ebp1;
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
                        countdn = (nel-1)*ebp1;
                        countup = 0;
                        for(j=0;j<nel;++j) {
                            for(k=0;k<ebp1;++k)
                                fsndbuf(countup +k) += frcvbuf(m,countdn +k);
                            countup += ebp1;
                            countdn -= ebp1;
                        }
                    }
                    
                    if (matches > 1) {
                        mtchinv = 1./matches;

#ifdef MPDEBUG
                        *sim::log << "finalrcv"  << idnum << " " << is_frst() << std::endl;
#endif
                        count = 0;
                        for(j=0;j<nel;++j) {
                            sind = el(j);
                            offset = sind*stride;
                            for (k=bgn;k<=end;++k) {
                                base[offset+k] = fsndbuf(count++)*mtchinv;
#ifdef MPDEBUG
                                *sim::log << "\t" << base[offset+k] << std::endl;
#endif
                            }

                        }
                    }
                    return;
                case(sum): 
                    for(m=0;m<nmatch;++m) {    
                        if (phase(grp)(m) != phi) continue;
                        
                        ++matches;
                        
                        int ebp1 = end-bgn+1;
                        countdn = (nel-1)*ebp1;
                        countup = 0;
                        for(j=0;j<nel;++j) {
                            for(k=0;k<ebp1;++k)
                                fsndbuf(countup +k) += frcvbuf(m,countdn +k);
                            countup += ebp1;
                            countdn -= ebp1;
                        }
                    }
                    
                    if (matches > 1) {
                        mtchinv = 1./matches;

#ifdef MPDEBUG
                        *sim::log << "finalrcv"  << idnum << " " << is_frst() << std::endl;
#endif
                        count = 0;
                        for(j=0;j<nel;++j) {
                            sind = el(j);
                            offset = sind*stride;
                            for (k=bgn;k<=end;++k) {
                                base[offset+k] = fsndbuf(count++);
#ifdef MPDEBUG
                                *sim::log << "\t" << base[offset+k] << std::endl;
#endif
                            }

                        }
                    }
                    return;
                case(maximum):                    
                    for(m=0;m<nmatch;++m) {    
                        if (phase(grp)(m) != phi) continue;
                        
                        ++matches;
                        
                        int ebp1 = end-bgn+1;
                        countdn = (nel-1)*ebp1;
                        countup = 0;
                        for(j=0;j<nel;++j) {
                            for(k=0;k<ebp1;++k)
                                fsndbuf(countup +k) = MAX(fsndbuf(countup+k),frcvbuf(m,countdn +k));
                            countup += ebp1;
                            countdn -= ebp1;
                        }
                    }
                    
                    if (matches > 1) {
#ifdef MPDEBUG
                        *sim::log << "finalrcv"  << idnum << " " << is_frst() << std::endl;
#endif
                        count = 0;
                        for(j=0;j<nel;++j) {
                            sind = el(j);
                            offset = sind*stride;
                            for (k=bgn;k<=end;++k) {
                                base[offset+k] = fsndbuf(count++);
#ifdef MPDEBUG
                                *sim::log << "\t" << base[offset+k] << std::endl;
#endif
                            }

                        }
                    }
                    return;
                    
                case(replace):
                    *sim::log << "Should only call replace with master_slave messages\n" << std::endl;
                    exit(1);
            }
        }
    }
}


void spartition::mgconnect(Array<mesh::transfer,1> &cnnct, mesh& tgt, int bnum) {
    int i,j,k,v0;
    
 
    /* BOUNDARY IS AN INTERNAL PARTITION BOUNDARY */
    /* MAKE SURE ENDPOINTS ARE OK */
    i = x.sd(el(0)).vrtx(0);
    if (cnnct(i).tri < 0) {
        tgt.qtree.nearpt(x.vrtx(i).data(),v0);
        cnnct(i).tri=tgt.vd(v0).tri;
        for(j=0;j<3;++j) {
            cnnct(i).wt(j) = 0.0;
            if (tgt.td(cnnct(i).tri).vrtx(j) == v0) cnnct(i).wt(j) = 1.0;
        }
    }
    i = x.sd(el(nel-1)).vrtx(1);
    if (cnnct(i).tri < 0) {
        tgt.qtree.nearpt(x.vrtx(i).data(),v0);
        cnnct(i).tri=tgt.vd(v0).tri;
        for(j=0;j<3;++j) {
            cnnct(i).wt(j) = 0.0;
            if (tgt.td(cnnct(i).tri).vrtx(j) == v0) cnnct(i).wt(j) = 1.0;
        }
    }
    
    if (first) {
        sndsize() = 0;
        sndtype() = int_msg;
        for(k=1;k<nel;++k) {
            v0 = x.sd(el(k)).vrtx(0);
            if (cnnct(v0).tri > 0) {
                isndbuf(sndsize()++) = -1;
            }
            else {
                isndbuf(sndsize()++) = +1;
                cnnct(v0).tri = 0;
                for(j=0;j<3;++j)
                    cnnct(v0).wt(j) = 0.0;
            }
        }
    }
    
    comm_prepare(boundary::all,0,slave_master); 
    comm_exchange(boundary::all,0,slave_master);
    comm_wait(boundary::all,0,slave_master);
    
    if (!first) {
        i = 0;
        for(k=nel-1;k>0;--k) {
            v0 = x.sd(el(k)).vrtx(1);
            if (ircvbuf(0,i) < 0) {
                cnnct(v0).tri = 0;
                for(j=0;j<3;++j)
                    cnnct(v0).wt(j) = 0.0;
            }
        }
    }                    
}

void curved_analytic_interface::mvpttobdry(TinyVector<FLT,mesh::ND> &pt) {
    int iter,n;
    FLT mag, delt_dist;
        
    /* FOR AN ANALYTIC SURFACE */
    iter = 0;
    do {
        mag = 0.0;
        for(n=0;n<mesh::ND;++n)
            mag += pow(dhgt(n,pt.data()),2);
        mag = sqrt(mag);
        delt_dist = -hgt(pt.data())/mag;
        for(n=0;n<mesh::ND;++n)
            pt(n) += delt_dist*dhgt(n,pt.data())/mag;
        if (++iter > 100) {
            *sim::log << "curved iterations exceeded curved boundary " << pt(0) << ' ' << pt(1) << '\n';
            exit(1);
        }
    } while (fabs(delt_dist) > 10.*EPSILON);
    
    return;
}

