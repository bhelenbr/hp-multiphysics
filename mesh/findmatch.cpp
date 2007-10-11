#include "mesh.h"
#include "boundaries.h"
#include <stdlib.h>
#include <new>

int tri_mesh::comm_entity_size() {
    int i,tsize,nvcomm,nscomm;
    
    /*	VERTEX INFO */    
    nvcomm = 0;
    for(i=0;i<nvbd;++i) 
        if (vbdry(i)->is_comm()) ++nvcomm;
        
    tsize = 1 + 2*nvcomm; // bdry number & id
    
    /* SIDE INFO */
    nscomm = 0;
    for(i=0;i<nsbd;++i)
        if (sbdry(i)->is_comm()) ++nscomm; 
    
    // tsize += 1 +4*nscomm; // bdry number, id, v0id, v1id
    tsize += 1 +2*nscomm; // bdry number, id

    tsize += 1;  // nfcomm = 0
    
    return(tsize);
}

int tri_mesh::comm_entity_list(Array<int,1>& list) {
    int i,nvcomm,nscomm,tsize;
    
    /* MAKE 1D PACKED LIST OF ALL INFORMATION ON COMMUNICATION BOUNDARIES */
    tsize = 0;

    /*	VERTEX INFO */    
    nvcomm = 0;
    for(i=0;i<nvbd;++i) 
        if (vbdry(i)->is_comm()) ++nvcomm;
        
    list(tsize++) = nvcomm;
    
    for(i=0;i<nvbd;++i) {
        if (vbdry(i)->is_comm()) {
            list(tsize++) = i;
            list(tsize++) = vbdry(i)->idnum;
        }
    }
    
    /* SIDE INFO */
    nscomm = 0;
    for(i=0;i<nsbd;++i)
        if (sbdry(i)->is_comm()) ++nscomm;
        
    list(tsize++) = nscomm;
    
    for(i=0;i<nsbd;++i) {
        if (sbdry(i)->is_comm()) {
            list(tsize++) = i;
            list(tsize++) = sbdry(i)->idnum;
#ifdef SKIP
            v0 = sd(sbdry(i)->el(0)).vrtx(0);
            v0id = -1;
            for(j=0;j<nvbd;++j) {
                if (vbdry(j)->v0 == v0) {
                    v0id = vbdry(j)->idnum;
                    break;
                }
            }
            list(tsize++) = v0id;
            v0 = sd(sbdry(i)->el(sbdry(i)->nel-1)]).vrtx(1);
            v0id = -1;
            for(j=0;j<nvbd;++j) {
                if (vbdry(j)->v0 == v0) {
                    v0id = vbdry(j)->idnum;
                    break;
                }
            }
            list(tsize++) = v0id;
#endif
        }
    }
    
    /* FACE BOUNDARIES */
    list(tsize++) = 0;
    
    return(tsize);
}

boundary* tri_mesh::getvbdry(int num) {return vbdry(num);}
boundary* tri_mesh::getsbdry(int num) {return sbdry(num);}
        
        
void tri_mesh::vmsgload(boundary::groups group, int phase, boundary::comm_type type, FLT *base,int bgn, int end, int stride) {
    int i;
        
    /* SEND COMMUNICATIONS TO ADJACENT MESHES */\
    for(i=0;i<nsbd;++i) 
        sbdry(i)->vloadbuff(group,base,bgn,end,stride);
    for(i=0;i<nvbd;++i)
        vbdry(i)->vloadbuff(group,base,bgn,end,stride);
    
    for(i=0;i<nsbd;++i)
        sbdry(i)->comm_prepare(group,phase,type);
    for(i=0;i<nvbd;++i)
        vbdry(i)->comm_prepare(group,phase,type);
    
    return;
}

void tri_mesh::vmsgpass(boundary::groups group, int phase, boundary::comm_type type) {

    for(int i=0;i<nsbd;++i) 
        sbdry(i)->comm_exchange(group,phase,type);
    for(int i=0;i<nvbd;++i) 
        vbdry(i)->comm_exchange(group,phase,type);

    return;
}

int tri_mesh::vmsgwait_rcv(boundary::groups group, int phase, boundary::comm_type type, boundary::operation op, FLT *base,int bgn, int end, int stride) {
    int stop = 1;
    int i;
    
    for(i=0;i<nsbd;++i)
        stop &= sbdry(i)->comm_wait(group,phase,type);
    for(i=0;i<nvbd;++i)
        stop &= vbdry(i)->comm_wait(group,phase,type);
        

    for(i=0;i<nsbd;++i) 
        sbdry(i)->vfinalrcv(group,phase,type,op,base,bgn,end,stride);
    for(i=0;i<nvbd;++i)
        vbdry(i)->vfinalrcv(group,phase,type,op,base,bgn,end,stride);
        
    return(stop);
}

int tri_mesh::vmsgrcv(boundary::groups group, int phase, boundary::comm_type type, boundary::operation op, FLT *base,int bgn, int end, int stride) {
    int stop = 1,i;
    
    for(i=0;i<nsbd;++i)
        stop &= sbdry(i)->comm_nowait(group,phase,type);
    for(i=0;i<nvbd;++i)
        stop &= vbdry(i)->comm_nowait(group,phase,type);
        
    for(i=0;i<nsbd;++i) 
        sbdry(i)->vfinalrcv(group,phase,type,op,base,bgn,end,stride);
    for(i=0;i<nvbd;++i)
        vbdry(i)->vfinalrcv(group,phase,type,op,base,bgn,end,stride);
        
    return(stop);
}


void tri_mesh::smsgload(boundary::groups group, int phase, boundary::comm_type type, FLT *base,int bgn, int end, int stride) {
    int i;
        
    /* SEND COMMUNICATIONS TO ADJACENT MESHES */\
    for(i=0;i<nsbd;++i) 
        sbdry(i)->sloadbuff(group,base,bgn,end,stride);
    
    for(i=0;i<nsbd;++i)
        sbdry(i)->comm_prepare(group,phase,type);
    
    return;
}

void tri_mesh::smsgpass(boundary::groups group, int phase, boundary::comm_type type) {

    for(int i=0;i<nsbd;++i) 
        sbdry(i)->comm_exchange(group,phase,type);

    return;
}

int tri_mesh::smsgwait_rcv(boundary::groups group, int phase, boundary::comm_type type, boundary::operation op, FLT *base,int bgn, int end, int stride) {
    int stop = 1;
    int i;
    
    for(i=0;i<nsbd;++i)
        stop &= sbdry(i)->comm_wait(group,phase,type);

    for(i=0;i<nsbd;++i) 
        sbdry(i)->sfinalrcv(group,phase,type,op,base,bgn,end,stride);
        
    return(stop);
}

int tri_mesh::smsgrcv(boundary::groups group,int phase, boundary::comm_type type, boundary::operation op, FLT *base,int bgn, int end, int stride) {
    int stop = 1,i;
    
    for(i=0;i<nsbd;++i)
        stop &= sbdry(i)->comm_nowait(group,phase,type);
        
    for(i=0;i<nsbd;++i) 
        sbdry(i)->sfinalrcv(group,phase,type,op,base,bgn,end,stride);
        
    return(stop);
}

void tri_mesh::matchboundaries() {
    int last_phase;
    int mp_phase;
        
    for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
    
        /* LOAD POSITIONS INTO BUFFERS */
        for(int i=0;i<nvbd;++i)
            vbdry(i)->loadpositions();
        for(int i=0;i<nsbd;++i) 
            sbdry(i)->loadpositions();
            
        /* FIRST PART OF SENDING, POST ALL RECEIVES */
        for(int i=0;i<nsbd;++i)
            sbdry(i)->comm_prepare(boundary::all_phased,mp_phase,boundary::master_slave);
        for(int i=0;i<nvbd;++i)
            vbdry(i)->comm_prepare(boundary::all_phased,mp_phase,boundary::master_slave);
                            
                
        /* SECOND PART */
        vmsgpass(boundary::all_phased,mp_phase,boundary::master_slave);
      
        /* FINAL PART OF SENDING */
        last_phase = true;
        for(int i=0;i<nsbd;++i)
            last_phase &= sbdry(i)->comm_wait(boundary::all_phased,mp_phase,boundary::master_slave);
                            
        for(int i=0;i<nvbd;++i)
            last_phase &= vbdry(i)->comm_wait(boundary::all_phased,mp_phase,boundary::master_slave);
                            
        for(int i=0;i<nsbd;++i)
            sbdry(i)->rcvpositions(mp_phase);
        for(int i=0;i<nvbd;++i)
            vbdry(i)->rcvpositions(mp_phase);
    }

    return;
}

#ifdef METIS

extern "C" void METIS_PartMeshNodal(int *ne, int *nn, int *elmnts, int *etype, int *numflag, int *nparts, int *edgecut,
int *epart, int *npart);

void tri_mesh::setpartition(int nparts) {
    int i,n;
    int etype = 1;
    int numflag = 0;
    int edgecut;
    
    /* CREATE ISOLATED TVRTX ARRAY */
    Array<TinyVector<int,3>,1> tvrtx(maxvst);    
    for(i=0;i<ntri;++i)
        for(n=0;n<3;++n)
            tvrtx(i)(n) = td(i).vrtx(n);
            
    METIS_PartMeshNodal(&ntri, &nvrtx, &tvrtx(0)(0), &etype, &numflag, &nparts, &edgecut,&(gbl->intwk(0)),&(gbl->i2wk(0)));
    
    for(i=0;i<ntri;++i)
        td(i).info = gbl->intwk(i);
        
    gbl->intwk = -1;

    return;
}
#endif

/* WHEN THIS ROUTINE EXITS */
/* THE OLD MESH STORES  */
/* vd(vind).info = new vrtx index or -1 */
/* THE NEW MESH STORES */
/* td(tind).info = old tri index */ 
void tri_mesh::partition(class tri_mesh& xin, int npart) {
    int i,j,n,tind,sind,v0,indx;
    Array<int,2> bcntr(xin.nsbd +5,3); 
    int bnum,bel,match;
    
    for(i=0;i<xin.nvrtx;++i)
        xin.vd(i).info = -1;

    ntri = 0;
    for(i=0;i<xin.ntri;++i) {
        if (xin.td(i).info == npart) {
            ++ntri;
            for(n=0;n<3;++n)
                xin.vd(xin.td(i).vrtx(n)).info = npart;
        }
    }
    
    *gbl->log << "New mesh with " << ntri << " of " << xin.ntri << " tris\n";
    
    if (!initialized) {
        maxvst = static_cast<int>(static_cast<FLT>(ntri*xin.maxvst)/xin.ntri);
        allocate(maxvst);
    }
    else if (3*ntri > maxvst) {
        *gbl->log << "mesh is too small" << std::endl;
        exit(1);
    }

    nvrtx = 0;
    for(i=0;i<xin.nvrtx;++i) {
        if (xin.vd(i).info == npart) {
            for(n=0;n<ND;++n)
                vrtx(nvrtx)(n) = xin.vrtx(i)(n);
            xin.vd(i).info = nvrtx;
            ++nvrtx;
        }
        else
            xin.vd(i).info = -1;
    }

    ntri = 0;
    for(i=0;i<xin.ntri;++i) {
        if (xin.td(i).info == npart) {
            for(n=0;n<3;++n)
                td(ntri).vrtx(n) = xin.vd(xin.td(i).vrtx(n)).info;
            td(ntri).info = i;
            ++ntri;
        }
    }

    createsideinfo();
    
    nsbd = 0;
    for(i=0;i<nside;++i) {
        sd(i).info = -1;
        if (sd(i).tri(1) < 0) {
            tind = td(sd(i).tri(0)).info;

            v0 = sd(i).vrtx(0);
            for(n=0;n<3;++n)
                if (xin.vd(xin.td(tind).vrtx(n)).info == v0) break;
            if (n==3) *gbl->log << "error in partitioning\n";
            n = (n+2)%3;

            indx = xin.td(tind).tri(n);
            if (indx < 0) {
                /* BOUNDARY SIDE */
                bnum = getbdrynum(indx);
                bel = getbdryel(indx);

                for (j = 0; j <nsbd;++j) {
                    if (bnum == bcntr(j,0)) {
                        ++bcntr(j,1);
                        sd(i).info = j;
                        goto next1;
                    }
                }
                /* NEW SIDE */
                sd(i).info = nsbd;
                bcntr(nsbd,0) = bnum;
                bcntr(nsbd++,1) = 1;
            }
            else {
                /* PARTITION SIDE */
                match = xin.td(indx).info;
                if (match < npart) bnum = (match<<16) + (npart << 24);
                else bnum = (npart<<16) + (match << 24);
                for (j = 0; j <nsbd;++j) {
                    if (bcntr(j,0) == -bnum) {
                        ++bcntr(j,1);
                        sd(i).info = j;
                        goto next1;
                    }
                }
                /* NEW SIDE */
                sd(i).info = nsbd;
                bcntr(nsbd,0) = -bnum;
                bcntr(nsbd++,1) = 1;
            }
        }
        next1: continue;
    }

    sbdry.resize(nsbd);
    for(i=0;i<nsbd;++i) {
        if (bcntr(i,0) < 0) 
            sbdry(i) = new spartition(abs(bcntr(i,0)),*this);
        else
            sbdry(i) = xin.sbdry(bcntr(i,0))->create(*this);
        sbdry(i)->alloc(static_cast<int>(bcntr(i,1)*2));
        sbdry(i)->nel = 0;
    }         
    
    
    for(i=0;i<nside;++i) {
        if (sd(i).info > -1) 
            sbdry(sd(i).info)->el(sbdry(sd(i).info)->nel++) = i;
    }

    for(i=0;i<nsbd;++i) {
        /* CREATES NEW BOUNDARY FOR DISCONNECTED SEGMENTS OF SAME TYPE */
        sbdry(i)->reorder();
    }
    
    /* MOVE VERTEX BOUNDARY INFO */
    nvbd = 0;
    for(i=0;i<xin.nvbd;++i) 
        if (xin.vd(xin.vbdry(i)->v0).info > -1)
            ++nvbd;
    vbdry.resize(nvbd+2*nsbd);
    
    nvbd = 0;
    for(i=0;i<xin.nvbd;++i) {
        if (xin.vd(xin.vbdry(i)->v0).info > -1) {
            vbdry(nvbd) = xin.vbdry(i)->create(*this);
            vbdry(nvbd)->alloc(4);
            vbdry(nvbd)->v0 = xin.vd(xin.vbdry(i)->v0).info;
            ++nvbd;
        }
    }
    
    
    /* CREATE COMMUNICATION ENDPOINT BOUNDARIES */
    for(i=0;i<nsbd;++i) {
        if (sbdry(i)->mytype == "partition") {
            sind = sbdry(i)->el(0);
            v0 = sd(sind).vrtx(0);
            for(j=0;j<nvbd;++j)
                if (vbdry(j)->v0 == v0) goto nextv0;
         
            /* ENDPOINT IS NEW NEED TO DEFINE BOUNDARY */
            tind = td(sd(sind).tri(0)).info;
            for(n=0;n<3;++n)
                if (xin.vd(xin.td(tind).vrtx(n)).info == v0) break;
            vbdry(nvbd) = new vcomm(xin.td(tind).vrtx(n),*this);
            vbdry(nvbd)->alloc(4);
            vbdry(nvbd)->v0 = v0;
            ++nvbd;
        
            nextv0:
            sind = sbdry(i)->el(sbdry(i)->nel-1);
            v0 = sd(sind).vrtx(1);
            for(j=0;j<nvbd;++j)
                if (vbdry(j)->v0 == v0) goto nextv1;
         
            /* NEW ENDPOINT */
            tind = td(sd(sind).tri(0)).info;
            for(n=0;n<3;++n)
                if (xin.vd(xin.td(tind).vrtx(n)).info == v0) break;
            vbdry(nvbd) = new vcomm(xin.td(tind).vrtx(n),*this);
            vbdry(nvbd)->alloc(4);
            vbdry(nvbd)->v0 = v0;
            ++nvbd;

            nextv1:
            continue;
        }
    }
    
    bdrylabel();  // CHANGES STRI / TTRI ON BOUNDARIES TO POINT TO GROUP/ELEMENT

    createttri();
    createvtri();
    cnt_nbor();
    FLT xmin[ND], xmax[ND];
    for(n=0;n<ND;++n) {
        xmin[n] = xin.qtree.xmin(n);
        xmax[n] = xin.qtree.xmax(n);
    }
    treeinit(xmin,xmax);

    initialized = 1;

    return;
}
                
    
    
    
            

