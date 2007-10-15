#include "tri_mesh.h"
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
    for(i=0;i<nebd;++i)
        if (ebdry(i)->is_comm()) ++nscomm; 
    
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
    for(i=0;i<nebd;++i)
        if (ebdry(i)->is_comm()) ++nscomm;
        
    list(tsize++) = nscomm;
    
    for(i=0;i<nebd;++i) {
        if (ebdry(i)->is_comm()) {
            list(tsize++) = i;
            list(tsize++) = ebdry(i)->idnum;
#ifdef SKIP
            p0 = seg(ebdry(i)->el(0)).pnt(0);
            v0id = -1;
            for(j=0;j<nvbd;++j) {
                if (vbdry(j)->p0 == p0) {
                    v0id = vbdry(j)->idnum;
                    break;
                }
            }
            list(tsize++) = v0id;
            p0 = seg(ebdry(i)->el(ebdry(i)->nel-1)]).pnt(1);
            v0id = -1;
            for(j=0;j<nvbd;++j) {
                if (vbdry(j)->p0 == p0) {
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
boundary* tri_mesh::getebdry(int num) {return ebdry(num);}
        
        
void tri_mesh::pmsgload(boundary::groups group, int phase, boundary::comm_type type, FLT *base,int bgn, int end, int stride) {
    int i;
        
    /* SEND COMMUNICATIONS TO ADJACENT MESHES */\
    for(i=0;i<nebd;++i) 
        ebdry(i)->vloadbuff(group,base,bgn,end,stride);
    for(i=0;i<nvbd;++i)
        vbdry(i)->vloadbuff(group,base,bgn,end,stride);
    
    for(i=0;i<nebd;++i)
        ebdry(i)->comm_prepare(group,phase,type);
    for(i=0;i<nvbd;++i)
        vbdry(i)->comm_prepare(group,phase,type);
    
    return;
}

void tri_mesh::pmsgpass(boundary::groups group, int phase, boundary::comm_type type) {

    for(int i=0;i<nebd;++i) 
        ebdry(i)->comm_exchange(group,phase,type);
    for(int i=0;i<nvbd;++i) 
        vbdry(i)->comm_exchange(group,phase,type);

    return;
}

int tri_mesh::pmsgwait_rcv(boundary::groups group, int phase, boundary::comm_type type, boundary::operation op, FLT *base,int bgn, int end, int stride) {
    int stop = 1;
    int i;
    
    for(i=0;i<nebd;++i)
        stop &= ebdry(i)->comm_wait(group,phase,type);
    for(i=0;i<nvbd;++i)
        stop &= vbdry(i)->comm_wait(group,phase,type);
        

    for(i=0;i<nebd;++i) 
        ebdry(i)->vfinalrcv(group,phase,type,op,base,bgn,end,stride);
    for(i=0;i<nvbd;++i)
        vbdry(i)->vfinalrcv(group,phase,type,op,base,bgn,end,stride);
        
    return(stop);
}

int tri_mesh::pmsgrcv(boundary::groups group, int phase, boundary::comm_type type, boundary::operation op, FLT *base,int bgn, int end, int stride) {
    int stop = 1,i;
    
    for(i=0;i<nebd;++i)
        stop &= ebdry(i)->comm_nowait(group,phase,type);
    for(i=0;i<nvbd;++i)
        stop &= vbdry(i)->comm_nowait(group,phase,type);
        
    for(i=0;i<nebd;++i) 
        ebdry(i)->vfinalrcv(group,phase,type,op,base,bgn,end,stride);
    for(i=0;i<nvbd;++i)
        vbdry(i)->vfinalrcv(group,phase,type,op,base,bgn,end,stride);
        
    return(stop);
}


void tri_mesh::smsgload(boundary::groups group, int phase, boundary::comm_type type, FLT *base,int bgn, int end, int stride) {
    int i;
        
    /* SEND COMMUNICATIONS TO ADJACENT MESHES */\
    for(i=0;i<nebd;++i) 
        ebdry(i)->sloadbuff(group,base,bgn,end,stride);
    
    for(i=0;i<nebd;++i)
        ebdry(i)->comm_prepare(group,phase,type);
    
    return;
}

void tri_mesh::smsgpass(boundary::groups group, int phase, boundary::comm_type type) {

    for(int i=0;i<nebd;++i) 
        ebdry(i)->comm_exchange(group,phase,type);

    return;
}

int tri_mesh::smsgwait_rcv(boundary::groups group, int phase, boundary::comm_type type, boundary::operation op, FLT *base,int bgn, int end, int stride) {
    int stop = 1;
    int i;
    
    for(i=0;i<nebd;++i)
        stop &= ebdry(i)->comm_wait(group,phase,type);

    for(i=0;i<nebd;++i) 
        ebdry(i)->sfinalrcv(group,phase,type,op,base,bgn,end,stride);
        
    return(stop);
}

int tri_mesh::smsgrcv(boundary::groups group,int phase, boundary::comm_type type, boundary::operation op, FLT *base,int bgn, int end, int stride) {
    int stop = 1,i;
    
    for(i=0;i<nebd;++i)
        stop &= ebdry(i)->comm_nowait(group,phase,type);
        
    for(i=0;i<nebd;++i) 
        ebdry(i)->sfinalrcv(group,phase,type,op,base,bgn,end,stride);
        
    return(stop);
}

void tri_mesh::matchboundaries() {
    int last_phase;
    int mp_phase;
        
    for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
    
        /* LOAD POSITIONS INTO BUFFERS */
        for(int i=0;i<nvbd;++i)
            vbdry(i)->loadpositions();
        for(int i=0;i<nebd;++i) 
            ebdry(i)->loadpositions();
            
        /* FIRST PART OF SENDING, POST ALL RECEIVES */
        for(int i=0;i<nebd;++i)
            ebdry(i)->comm_prepare(boundary::all_phased,mp_phase,boundary::master_slave);
        for(int i=0;i<nvbd;++i)
            vbdry(i)->comm_prepare(boundary::all_phased,mp_phase,boundary::master_slave);
                            
                
        /* SECOND PART */
        pmsgpass(boundary::all_phased,mp_phase,boundary::master_slave);
      
        /* FINAL PART OF SENDING */
        last_phase = true;
        for(int i=0;i<nebd;++i)
            last_phase &= ebdry(i)->comm_wait(boundary::all_phased,mp_phase,boundary::master_slave);
                            
        for(int i=0;i<nvbd;++i)
            last_phase &= vbdry(i)->comm_wait(boundary::all_phased,mp_phase,boundary::master_slave);
                            
        for(int i=0;i<nebd;++i)
            ebdry(i)->rcvpositions(mp_phase);
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
    Array<TinyVector<int,3>,1> tvrtx(maxpst);    
    for(i=0;i<ntri;++i)
        for(n=0;n<3;++n)
            tvrtx(i)(n) = tri(i).pnt(n);
            
    METIS_PartMeshNodal(&ntri, &npnt, &tvrtx(0)(0), &etype, &numflag, &nparts, &edgecut,&(gbl->intwk(0)),&(gbl->i2wk(0)));
    
    for(i=0;i<ntri;++i)
        tri(i).info = gbl->intwk(i);
        
    gbl->intwk = -1;

    return;
}
#endif

/* WHEN THIS ROUTINE EXITS */
/* THE OLD MESH STORES  */
/* pnt(pind).info = new pnt index or -1 */
/* THE NEW MESH STORES */
/* tri(tind).info = old tri index */ 
void tri_mesh::partition(class tri_mesh& xin, int npart) {
    int i,j,n,tind,sind,p0,indx;
    Array<int,2> bcntr(xin.nebd +5,3); 
    int bnum,bel,match;
    
    for(i=0;i<xin.npnt;++i)
        xin.pnt(i).info = -1;

    ntri = 0;
    for(i=0;i<xin.ntri;++i) {
        if (xin.tri(i).info == npart) {
            ++ntri;
            for(n=0;n<3;++n)
                xin.pnt(xin.tri(i).pnt(n)).info = npart;
        }
    }
    
    *gbl->log << "New mesh with " << ntri << " of " << xin.ntri << " tris\n";
    
    if (!initialized) {
        maxpst = static_cast<int>(static_cast<FLT>(ntri*xin.maxpst)/xin.ntri);
        allocate(maxpst);
    }
    else if (3*ntri > maxpst) {
        *gbl->log << "mesh is too small" << std::endl;
        exit(1);
    }

    npnt = 0;
    for(i=0;i<xin.npnt;++i) {
        if (xin.pnt(i).info == npart) {
            for(n=0;n<ND;++n)
                pnts(npnt)(n) = xin.pnts(i)(n);
            xin.pnt(i).info = npnt;
            ++npnt;
        }
        else
            xin.pnt(i).info = -1;
    }

    ntri = 0;
    for(i=0;i<xin.ntri;++i) {
        if (xin.tri(i).info == npart) {
            for(n=0;n<3;++n)
                tri(ntri).pnt(n) = xin.pnt(xin.tri(i).pnt(n)).info;
            tri(ntri).info = i;
            ++ntri;
        }
    }

    createseg();
    
    nebd = 0;
    for(i=0;i<nseg;++i) {
        seg(i).info = -1;
        if (seg(i).tri(1) < 0) {
            tind = tri(seg(i).tri(0)).info;

            p0 = seg(i).pnt(0);
            for(n=0;n<3;++n)
                if (xin.pnt(xin.tri(tind).pnt(n)).info == p0) break;
            if (n==3) *gbl->log << "error in partitioning\n";
            n = (n+2)%3;

            indx = xin.tri(tind).tri(n);
            if (indx < 0) {
                /* BOUNDARY SIDE */
                bnum = getbdrynum(indx);
                bel = getbdryel(indx);

                for (j = 0; j <nebd;++j) {
                    if (bnum == bcntr(j,0)) {
                        ++bcntr(j,1);
                        seg(i).info = j;
                        goto next1;
                    }
                }
                /* NEW SIDE */
                seg(i).info = nebd;
                bcntr(nebd,0) = bnum;
                bcntr(nebd++,1) = 1;
            }
            else {
                /* PARTITION SIDE */
                match = xin.tri(indx).info;
                if (match < npart) bnum = (match<<16) + (npart << 24);
                else bnum = (npart<<16) + (match << 24);
                for (j = 0; j <nebd;++j) {
                    if (bcntr(j,0) == -bnum) {
                        ++bcntr(j,1);
                        seg(i).info = j;
                        goto next1;
                    }
                }
                /* NEW SIDE */
                seg(i).info = nebd;
                bcntr(nebd,0) = -bnum;
                bcntr(nebd++,1) = 1;
            }
        }
        next1: continue;
    }

    ebdry.resize(nebd);
    for(i=0;i<nebd;++i) {
        if (bcntr(i,0) < 0) 
            ebdry(i) = new epartition(abs(bcntr(i,0)),*this);
        else
            ebdry(i) = xin.ebdry(bcntr(i,0))->create(*this);
        ebdry(i)->alloc(static_cast<int>(bcntr(i,1)*2));
        ebdry(i)->nel = 0;
    }         
    
    
    for(i=0;i<nseg;++i) {
        if (seg(i).info > -1) 
            ebdry(seg(i).info)->el(ebdry(seg(i).info)->nel++) = i;
    }

    for(i=0;i<nebd;++i) {
        /* CREATES NEW BOUNDARY FOR DISCONNECTED SEGMENTS OF SAME TYPE */
        ebdry(i)->reorder();
    }
    
    /* MOVE VERTEX BOUNDARY INFO */
    nvbd = 0;
    for(i=0;i<xin.nvbd;++i) 
        if (xin.pnt(xin.vbdry(i)->p0).info > -1)
            ++nvbd;
    vbdry.resize(nvbd+2*nebd);
    
    nvbd = 0;
    for(i=0;i<xin.nvbd;++i) {
        if (xin.pnt(xin.vbdry(i)->p0).info > -1) {
            vbdry(nvbd) = xin.vbdry(i)->create(*this);
            vbdry(nvbd)->alloc(4);
            vbdry(nvbd)->p0 = xin.pnt(xin.vbdry(i)->p0).info;
            ++nvbd;
        }
    }
    
    
    /* CREATE COMMUNICATION ENDPOINT BOUNDARIES */
    for(i=0;i<nebd;++i) {
        if (ebdry(i)->mytype == "partition") {
            sind = ebdry(i)->el(0);
            p0 = seg(sind).pnt(0);
            for(j=0;j<nvbd;++j)
                if (vbdry(j)->p0 == p0) goto nextv0;
         
            /* ENDPOINT IS NEW NEED TO DEFINE BOUNDARY */
            tind = tri(seg(sind).tri(0)).info;
            for(n=0;n<3;++n)
                if (xin.pnt(xin.tri(tind).pnt(n)).info == p0) break;
            vbdry(nvbd) = new vcomm(xin.tri(tind).pnt(n),*this);
            vbdry(nvbd)->alloc(4);
            vbdry(nvbd)->p0 = p0;
            ++nvbd;
        
            nextv0:
            sind = ebdry(i)->el(ebdry(i)->nel-1);
            p0 = seg(sind).pnt(1);
            for(j=0;j<nvbd;++j)
                if (vbdry(j)->p0 == p0) goto nextv1;
         
            /* NEW ENDPOINT */
            tind = tri(seg(sind).tri(0)).info;
            for(n=0;n<3;++n)
                if (xin.pnt(xin.tri(tind).pnt(n)).info == p0) break;
            vbdry(nvbd) = new vcomm(xin.tri(tind).pnt(n),*this);
            vbdry(nvbd)->alloc(4);
            vbdry(nvbd)->p0 = p0;
            ++nvbd;

            nextv1:
            continue;
        }
    }
    
    bdrylabel();  // CHANGES STRI / TTRI ON BOUNDARIES TO POINT TO GROUP/ELEMENT

    createtritri();
    createpnttri();
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
                
    
    
    
            

