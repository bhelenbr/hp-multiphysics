#include "tet_mesh.h"
#include <tet_boundary.h>
#include <stdlib.h>
#include <new>

int tet_mesh::comm_entity_size() {
	int i,tsize,nvcomm,necomm,nfcomm;
	
	/*    VERTEX INFO */    
	nvcomm = 0;
	for(i=0;i<nvbd;++i) 
		if (vbdry(i)->is_comm()) ++nvcomm;
		
	tsize = 1 + 2*nvcomm; // bdry number & id
	
	/* SIDE INFO */
	necomm = 0;
	for(i=0;i<nebd;++i)
		if (ebdry(i)->is_comm()) ++necomm; 
	
	tsize += 1 +2*necomm; // bdry number & id

	/* SIDE INFO */
	nfcomm = 0;
	for(i=0;i<nfbd;++i)
		if (fbdry(i)->is_comm()) ++nfcomm;     
	
	tsize += 1 +2*nfcomm;  //  bdry number & id
	
	return(tsize);
}

int tet_mesh::comm_entity_list(Array<int,1>& list) {
	int i,nvcomm,necomm,nfcomm,tsize;
	
	/* MAKE 1D PACKED LIST OF ALL INFORMATION ON COMMUNICATION BOUNDARIES */
	tsize = 0;

	/*    VERTEX INFO */    
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
	
	/* EDGE INFO */
	necomm = 0;
	for(i=0;i<nebd;++i)
		if (ebdry(i)->is_comm()) ++necomm;
		
	list(tsize++) = necomm;
	
	for(i=0;i<nebd;++i) {
		if (ebdry(i)->is_comm()) {
			list(tsize++) = i;
			list(tsize++) = ebdry(i)->idnum;
		}
	}
	
	/* FACE INFO */
	nfcomm = 0;
	for(i=0;i<nfbd;++i)
		if (fbdry(i)->is_comm()) ++nfcomm;
		
	list(tsize++) = nfcomm;
	
	for(i=0;i<nfbd;++i) {
		if (fbdry(i)->is_comm()) {
			list(tsize++) = i;
			list(tsize++) = fbdry(i)->idnum;
		}
	}
	
	return(tsize);
}

boundary* tet_mesh::getvbdry(int num) {return vbdry(num);}
boundary* tet_mesh::getebdry(int num) {return ebdry(num);}
boundary* tet_mesh::getfbdry(int num) {return fbdry(num);}

		
		
void tet_mesh::pmsgload(boundary::groups group, int phase, boundary::comm_type type, FLT *base,int bgn, int end, int stride) {
	int i;
		
	/* SEND COMMUNICATIONS TO ADJACENT MESHES */
	for(i=0;i<nfbd;++i) 
		fbdry(i)->ploadbuff(group,base,bgn,end,stride);    
	for(i=0;i<nebd;++i) 
		ebdry(i)->ploadbuff(group,base,bgn,end,stride);
	for(i=0;i<nvbd;++i)
		vbdry(i)->ploadbuff(group,base,bgn,end,stride);

	for(i=0;i<nfbd;++i)
		fbdry(i)->comm_prepare(group,phase,type);   
	for(i=0;i<nebd;++i)
		ebdry(i)->comm_prepare(group,phase,type);
	for(i=0;i<nvbd;++i)
		vbdry(i)->comm_prepare(group,phase,type);
	
	return;
}

void tet_mesh::pmsgpass(boundary::groups group, int phase, boundary::comm_type type) {

	for(int i=0;i<nfbd;++i) 
		fbdry(i)->comm_exchange(group,phase,type);
	for(int i=0;i<nebd;++i) 
		ebdry(i)->comm_exchange(group,phase,type);
	for(int i=0;i<nvbd;++i) 
		vbdry(i)->comm_exchange(group,phase,type);

	return;
}

int tet_mesh::pmsgwait_rcv(boundary::groups group, int phase, boundary::comm_type type, boundary::operation op, FLT *base,int bgn, int end, int stride) {
	int stop = 1;
	int i;
	
	for(i=0;i<nfbd;++i)
		stop &= fbdry(i)->comm_wait(group,phase,type);    
	for(i=0;i<nebd;++i)
		stop &= ebdry(i)->comm_wait(group,phase,type);
	for(i=0;i<nvbd;++i)
		stop &= vbdry(i)->comm_wait(group,phase,type);
		
	for(i=0;i<nfbd;++i) 
		fbdry(i)->pfinalrcv(group,phase,type,op,base,bgn,end,stride);
	for(i=0;i<nebd;++i) 
		ebdry(i)->pfinalrcv(group,phase,type,op,base,bgn,end,stride);
	for(i=0;i<nvbd;++i)
		vbdry(i)->pfinalrcv(group,phase,type,op,base,bgn,end,stride);
		
	return(stop);
}

int tet_mesh::pmsgrcv(boundary::groups group, int phase, boundary::comm_type type, boundary::operation op, FLT *base,int bgn, int end, int stride) {
	int stop = 1,i;

	for(i=0;i<nfbd;++i)
		stop &= fbdry(i)->comm_nowait(group,phase,type);    
	for(i=0;i<nebd;++i)
		stop &= ebdry(i)->comm_nowait(group,phase,type);
	for(i=0;i<nvbd;++i)
		stop &= vbdry(i)->comm_nowait(group,phase,type);
		
	for(i=0;i<nfbd;++i) 
		fbdry(i)->pfinalrcv(group,phase,type,op,base,bgn,end,stride);
	for(i=0;i<nebd;++i) 
		ebdry(i)->pfinalrcv(group,phase,type,op,base,bgn,end,stride);
	for(i=0;i<nvbd;++i)
		vbdry(i)->pfinalrcv(group,phase,type,op,base,bgn,end,stride);
		
	return(stop);
}


void tet_mesh::smsgload(boundary::groups group, int phase, boundary::comm_type type, FLT *base,int bgn, int end, int stride) {
	int i;
		
	/* SEND COMMUNICATIONS TO ADJACENT MESHES */
	for(i=0;i<nfbd;++i) 
		fbdry(i)->sloadbuff(group,base,bgn,end,stride);   
	for(i=0;i<nebd;++i) 
		ebdry(i)->sloadbuff(group,base,bgn,end,stride);

	for(i=0;i<nfbd;++i)
		fbdry(i)->comm_prepare(group,phase,type);    
	for(i=0;i<nebd;++i)
		ebdry(i)->comm_prepare(group,phase,type);
	
	return;
}

void tet_mesh::smsgpass(boundary::groups group, int phase, boundary::comm_type type) {

	for(int i=0;i<nfbd;++i) 
		fbdry(i)->comm_exchange(group,phase,type);
	for(int i=0;i<nebd;++i) 
		ebdry(i)->comm_exchange(group,phase,type);

	return;
}

int tet_mesh::smsgwait_rcv(boundary::groups group, int phase, boundary::comm_type type, boundary::operation op, FLT *base,int bgn, int end, int stride) {
	int stop = 1;
	int i;

	for(i=0;i<nfbd;++i)
		stop &= fbdry(i)->comm_wait(group,phase,type);   
	for(i=0;i<nebd;++i)
		stop &= ebdry(i)->comm_wait(group,phase,type);

	for(i=0;i<nfbd;++i) 
		fbdry(i)->sfinalrcv(group,phase,type,op,base,bgn,end,stride);
	for(i=0;i<nebd;++i) 
		ebdry(i)->sfinalrcv(group,phase,type,op,base,bgn,end,stride);
		
	return(stop);
}

int tet_mesh::smsgrcv(boundary::groups group,int phase, boundary::comm_type type, boundary::operation op, FLT *base,int bgn, int end, int stride) {
	int stop = 1,i;

	for(i=0;i<nfbd;++i)
		stop &= fbdry(i)->comm_nowait(group,phase,type);    
	for(i=0;i<nebd;++i)
		stop &= ebdry(i)->comm_nowait(group,phase,type);

	for(i=0;i<nfbd;++i) 
		fbdry(i)->sfinalrcv(group,phase,type,op,base,bgn,end,stride);        
	for(i=0;i<nebd;++i) 
		ebdry(i)->sfinalrcv(group,phase,type,op,base,bgn,end,stride);
		
	return(stop);
}


void tet_mesh::tmsgload(boundary::groups group, int phase, boundary::comm_type type, FLT *base,int bgn, int end, int stride) {
	int i;
		
	/* SEND COMMUNICATIONS TO ADJACENT MESHES */
	for(i=0;i<nfbd;++i) 
		fbdry(i)->tloadbuff(group,base,bgn,end,stride);  

	for(i=0;i<nfbd;++i)
		fbdry(i)->comm_prepare(group,phase,type);    
	
	return;
}

void tet_mesh::tmsgpass(boundary::groups group, int phase, boundary::comm_type type) {

	for(int i=0;i<nfbd;++i) 
		fbdry(i)->comm_exchange(group,phase,type);

	return;
}

int tet_mesh::tmsgwait_rcv(boundary::groups group, int phase, boundary::comm_type type, boundary::operation op, FLT *base,int bgn, int end, int stride) {
	int stop = 1;
	int i;

	for(i=0;i<nfbd;++i)
		stop &= fbdry(i)->comm_wait(group,phase,type);   

	for(i=0;i<nfbd;++i) 
		fbdry(i)->tfinalrcv(group,phase,type,op,base,bgn,end,stride);
		
	return(stop);
}

int tet_mesh::tmsgrcv(boundary::groups group,int phase, boundary::comm_type type, boundary::operation op, FLT *base,int bgn, int end, int stride) {
	int stop = 1,i;

	for(i=0;i<nfbd;++i)
		stop &= fbdry(i)->comm_nowait(group,phase,type);    

	for(i=0;i<nfbd;++i) 
		fbdry(i)->tfinalrcv(group,phase,type,op,base,bgn,end,stride);        
		
	return(stop);
}

void tet_mesh::matchboundaries() {
	int last_phase;
	int mp_phase;
		
	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
	
		/* LOAD POSITIONS INTO BUFFERS */
		for(int i=0;i<nvbd;++i)
			vbdry(i)->loadpositions();
		for(int i=0;i<nebd;++i) 
			ebdry(i)->loadpositions();
		for(int i=0;i<nfbd;++i) 
			fbdry(i)->loadpositions();
			
		/* FIRST PART OF SENDING, POST ALL RECEIVES */
		for(int i=0;i<nfbd;++i)
			fbdry(i)->comm_prepare(boundary::all_phased,mp_phase,boundary::master_slave);
		for(int i=0;i<nebd;++i)
			ebdry(i)->comm_prepare(boundary::all_phased,mp_phase,boundary::master_slave);
		for(int i=0;i<nvbd;++i)
			vbdry(i)->comm_prepare(boundary::all_phased,mp_phase,boundary::master_slave);
							
				
		/* SECOND PART */
		pmsgpass(boundary::all_phased,mp_phase,boundary::master_slave);

		/* FINAL PART OF SENDING */
		last_phase = true;
		for(int i=0;i<nfbd;++i)
			last_phase &= fbdry(i)->comm_wait(boundary::all_phased,mp_phase,boundary::master_slave);
		for(int i=0;i<nebd;++i)
			last_phase &= ebdry(i)->comm_wait(boundary::all_phased,mp_phase,boundary::master_slave);
		for(int i=0;i<nvbd;++i)
			last_phase &= vbdry(i)->comm_wait(boundary::all_phased,mp_phase,boundary::master_slave);

		for(int i=0;i<nfbd;++i)
			fbdry(i)->rcvpositions(mp_phase);                            
		for(int i=0;i<nebd;++i)
			ebdry(i)->rcvpositions(mp_phase);
		for(int i=0;i<nvbd;++i)
			vbdry(i)->rcvpositions(mp_phase);
	}

	return;
}


void tet_mesh::create_unique_numbering() {
	int sind,pnt0,phase;
	bool last_phase;
	
	/* vinfo will store unique numbers of each vertex */
	/* Find maximum number of vertices on any block */
	int maxpnts, sndpnt = npnt;
	sim::blks.allreduce(&sndpnt,&maxpnts,1,blocks::int_msg,blocks::max);

	for (int i=0;i<npnt;++i) 
		pnt(i).info = maxpnts*gbl->idnum +i;
	

	/* Find minimum along all communication boundaries */
	for(last_phase = false, phase = 0; !last_phase; ++phase) {
		for (int i=0;i<nfbd;++i) {
			if (!fbdry(i)->is_comm()) continue;
			
			for (int j=0;j<fbdry(i)->npnt;++j) {
				pnt0 = fbdry(i)->pnt(j).gindx;
				fbdry(i)->isndbuf(j) = pnt(pnt0).info;
			}
			fbdry(i)->sndsize() = fbdry(i)->npnt;
			fbdry(i)->sndtype() = boundary::int_msg;
		}

		for (int i=0;i<nebd;++i) {
			if (!ebdry(i)->is_comm()) continue;
			
			for (int j=0;j<ebdry(i)->nseg;++j) {
				sind = ebdry(i)->seg(j).gindx;
				pnt0 = seg(sind).pnt(0);
				ebdry(i)->isndbuf(j) = pnt(pnt0).info;
			}
			pnt0 = seg(sind).pnt(1);
			ebdry(i)->isndbuf(ebdry(i)->nseg) = pnt(pnt0).info;
			ebdry(i)->sndsize() = ebdry(i)->nseg+1;
			ebdry(i)->sndtype() = boundary::int_msg;
		}

		for (int i=0;i<nvbd;++i) {
			if (!vbdry(i)->is_comm()) continue;

			vbdry(i)->isndbuf(0) = pnt(vbdry(i)->pnt).info;
			vbdry(i)->sndsize() = 1;
			vbdry(i)->sndtype() = boundary::int_msg;
		}

		for(int i=0;i<nfbd;++i)
			fbdry(i)->comm_prepare(boundary::all_phased,phase,boundary::symmetric);   
		for(int i=0;i<nebd;++i)
			ebdry(i)->comm_prepare(boundary::all_phased,phase,boundary::symmetric);
		for(int i=0;i<nvbd;++i)
			vbdry(i)->comm_prepare(boundary::all_phased,phase,boundary::symmetric);
		
		for(int i=0;i<nfbd;++i) 
			fbdry(i)->comm_exchange(boundary::all_phased,phase,boundary::symmetric);
		for(int i=0;i<nebd;++i) 
			ebdry(i)->comm_exchange(boundary::all_phased,phase,boundary::symmetric);
		for(int i=0;i<nvbd;++i) 
			vbdry(i)->comm_exchange(boundary::all_phased,phase,boundary::symmetric);

		last_phase = true;
		for(int i=0;i<nfbd;++i) {
			last_phase &= fbdry(i)->comm_wait(boundary::all_phased,phase,boundary::symmetric); 
			fbdry(i)->comm_finish(boundary::all_phased,phase,boundary::symmetric,boundary::minimum);
		}
		for(int i=0;i<nebd;++i) {
			last_phase &= ebdry(i)->comm_wait(boundary::all_phased,phase,boundary::symmetric);
			ebdry(i)->comm_finish(boundary::all_phased,phase,boundary::symmetric,boundary::minimum);
		}
		for(int i=0;i<nvbd;++i) {
			last_phase &= vbdry(i)->comm_wait(boundary::all_phased,phase,boundary::symmetric);
			vbdry(i)->comm_finish(boundary::all_phased,phase,boundary::symmetric,boundary::minimum);
		}
		
		/* Find minimum along all communication boundaries */
		for (int i=0;i<nfbd;++i) {
			if (!fbdry(i)->is_comm()) continue;
			
			for (int j=0;j<fbdry(i)->npnt;++j) {
				pnt0 = fbdry(i)->pnt(j).gindx;
				pnt(pnt0).info = fbdry(i)->isndbuf(j);
			}
		}

		for (int i=0;i<nebd;++i) {
			if (!ebdry(i)->is_comm()) continue;
			
			for (int j=0;j<ebdry(i)->nseg;++j) {
				sind = ebdry(i)->seg(j).gindx;
				pnt0 = seg(sind).pnt(0);
				pnt(pnt0).info = ebdry(i)->isndbuf(j);
			}
			pnt0 = seg(sind).pnt(1);
			pnt(pnt0).info = ebdry(i)->isndbuf(ebdry(i)->nseg);
		}

		for (int i=0;i<nvbd;++i) {
			if (!vbdry(i)->is_comm()) continue;

			pnt(vbdry(i)->pnt).info = vbdry(i)->isndbuf(0);
		}
	}

	return;
}

void tet_mesh::match_bdry_numbering() {

	/* LOAD POSITIONS INTO BUFFERS */
	for(int i=0;i<nfbd;++i) 
		fbdry(i)->ploadbuff(boundary::all,&(pnts(0)(0)),0,ND-1,ND);
		
	/* FIRST PART OF SENDING, POST ALL RECEIVES */
	for(int i=0;i<nfbd;++i)
		fbdry(i)->comm_prepare(boundary::all,0,boundary::master_slave);
						
	/* SECOND PART */
	tmsgpass(boundary::all,0,boundary::master_slave);

	/* FINAL PART OF SENDING */
	for(int i=0;i<nfbd;++i)
		fbdry(i)->comm_wait(boundary::all,0,boundary::master_slave);
	
	/* FINAL PART OF SENDING */
	for(int i=0;i<nfbd;++i)
		fbdry(i)->comm_finish(boundary::all,0,boundary::master_slave,boundary::replace);

	/* Slave receives vertex positions & matches vertex numbering */
	for(int i=0;i<nfbd;++i)
		fbdry(i)->match_numbering(1);

	/* Create global numbering system */
	create_unique_numbering();
	 
	/* Reorder Side boundaries so direction is the same */
	for(int i=0;i<nebd;++i) 
		ebdry(i)->match_numbering(1);

	/* FIRST PART OF SENDING, POST ALL RECEIVES */
	for(int i=0;i<nebd;++i)
		ebdry(i)->comm_prepare(boundary::all,0,boundary::master_slave);
						
	/* SECOND PART */
	 for(int i=0;i<nebd;++i) 
		ebdry(i)->comm_exchange(boundary::all,0,boundary::master_slave);

	/* FINAL PART OF SENDING */
	for(int i=0;i<nebd;++i)
		ebdry(i)->comm_wait(boundary::all,0,boundary::master_slave);
	
	/* FINAL PART OF SENDING */
	for(int i=0;i<nebd;++i)
		ebdry(i)->comm_finish(boundary::all,0,boundary::master_slave,boundary::replace);

	for(int i=0;i<nebd;++i) 
		ebdry(i)->match_numbering(2);		  

	/* Redefine tets based on global numbering system */
	reorient_tets(true);
	match_all();
		
	/* Master loads integer data (seg, tri definitions) */
	for(int i=0;i<nfbd;++i)
		fbdry(i)->match_numbering(2);

	/* FIRST PART OF SENDING, POST ALL RECEIVES */
	for(int i=0;i<nfbd;++i)
		fbdry(i)->comm_prepare(boundary::all,0,boundary::master_slave);
						
	/* SECOND PART */
	tmsgpass(boundary::all,0,boundary::master_slave);

	/* FINAL PART OF SENDING */
	for(int i=0;i<nfbd;++i) {
		fbdry(i)->comm_wait(boundary::all,0,boundary::master_slave);
		fbdry(i)->comm_finish(boundary::all,0,boundary::master_slave,boundary::replace);
	}

	/* Slaves receive seg, tri definitions */
	for(int i=0;i<nfbd;++i)
		fbdry(i)->match_numbering(3);

	for(int i = 0; i < nfbd; ++i){
		if (!fbdry(i)->is_frst() && fbdry(i)->is_comm()) {
			fbdry(i)->match_tri_and_seg(); 
			fbdry(i)->create_tri_tri();
			fbdry(i)->create_pnt_nnbor();
			fbdry(i)->create_pnt_tri();
		}
	}  
	
#ifdef TEST_MATCH
	for (int i=0;i<nfbd;++i) {
		if (!fbdry(i)->is_comm()) continue;
		
		*gbl->log << fbdry(i)->idprefix << "\n";
		*gbl->log << "Reordered pnts\n";
		for (int j=0;j<fbdry(i)->npnt;++j) {
			*gbl->log << pnts(fbdry(i)->pnt(j).gindx) << std::endl;
		}
	
		if (fbdry(i)->is_frst()) {
			*gbl->log << "Reordered Side definitions\n";
			for (int j=0;j<fbdry(i)->nseg;++j) {
				*gbl->log << fbdry(i)->seg(j).pnt(0) << ' ' << fbdry(i)->seg(j).pnt(1) << std::endl;
				*gbl->log << fbdry(i)->seg(j).tri(0) << ' ' << fbdry(i)->seg(j).tri(1) << std::endl;
				*gbl->log << pnts(seg(fbdry(i)->seg(j).gindx).pnt(0)) << ' ' << pnts(seg(fbdry(i)->seg(j).gindx).pnt(0)) << std::endl;
			}
			
			
			*gbl->log << "Reordered Tri definitions\n";
			for (int j=0;j<fbdry(i)->ntri;++j) {
				*gbl->log << fbdry(i)->tri(j).pnt(0) << ' ' << fbdry(i)->tri(j).pnt(1) << ' ' << fbdry(i)->tri(j).pnt(2) << std::endl;
				*gbl->log << fbdry(i)->tri(j).seg(0) << ' ' << fbdry(i)->tri(j).seg(1) << ' ' << fbdry(i)->tri(j).seg(2) << std::endl;
				*gbl->log << fbdry(i)->tri(j).tri(0) << ' ' << fbdry(i)->tri(j).tri(1) << ' ' << fbdry(i)->tri(j).tri(2) << std::endl;
				*gbl->log << fbdry(i)->tri(j).sgn(0) << ' ' << fbdry(i)->tri(j).sgn(1) << ' ' << fbdry(i)->tri(j).sgn(2) << std::endl;
				
				*gbl->log << pnts(tri(fbdry(i)->tri(j).gindx).pnt(0)) << ' ' << pnts(tri(fbdry(i)->tri(j).gindx).pnt(1)) << ' ' << pnts(tri(fbdry(i)->tri(j).gindx).pnt(2)) << std::endl;
			}   
		}
		else {
			*gbl->log << "Reordered Side definitions\n";
			for (int j=0;j<fbdry(i)->nseg;++j) {
				*gbl->log << fbdry(i)->seg(j).pnt(0) << ' ' << fbdry(i)->seg(j).pnt(1) << std::endl;
				*gbl->log << fbdry(i)->seg(j).tri(1) << ' ' << fbdry(i)->seg(j).tri(0) << std::endl;
				*gbl->log << pnts(seg(fbdry(i)->seg(j).gindx).pnt(0)) << ' ' << pnts(seg(fbdry(i)->seg(j).gindx).pnt(0)) << std::endl;
			}
			
			*gbl->log << "Reordered Tri definitions\n";
			for (int j=0;j<fbdry(i)->ntri;++j) {
				*gbl->log << fbdry(i)->tri(j).pnt(0) << ' ' << fbdry(i)->tri(j).pnt(2) << ' ' << fbdry(i)->tri(j).pnt(1) << std::endl;
				*gbl->log << fbdry(i)->tri(j).seg(0) << ' ' << fbdry(i)->tri(j).seg(2) << ' ' << fbdry(i)->tri(j).seg(1) << std::endl;
				*gbl->log << fbdry(i)->tri(j).tri(0) << ' ' << fbdry(i)->tri(j).tri(2) << ' ' << fbdry(i)->tri(j).tri(1) << std::endl;
				*gbl->log << -fbdry(i)->tri(j).sgn(0) << ' ' << -fbdry(i)->tri(j).sgn(2) << ' ' << -fbdry(i)->tri(j).sgn(1) << std::endl;
				
				*gbl->log << pnts(tri(fbdry(i)->tri(j).gindx).pnt(0)) << ' ' << pnts(tri(fbdry(i)->tri(j).gindx).pnt(2)) << ' ' << pnts(tri(fbdry(i)->tri(j).gindx).pnt(1)) << std::endl;
			}   
		}
	}
#endif

	setinfo();
	
}


#ifdef METIS

extern "C" void METIS_PartMeshNodal(int *ne, int *nn, int *elmnts, int *etype, int *numflag, int *nparts, int *edgecut,
									int *epart, int *npart);

void tet_mesh::setpartition(int nparts) {
	int i,n;
	int etype = 2; /* 1 = triangle, 2 = tetrahedral*/
	int numflag = 0; /* C-style numbering */
	int edgecut;
	
	/* CREATE ISOLATED TVRTX ARRAY */
	Array<TinyVector<int,4>,1> tvrtx(maxvst);
	for(i=0;i<ntet;++i)
		for(n=0;n<4;++n)
			tvrtx(i)(n) = tet(i).pnt(n);
	
	METIS_PartMeshNodal(&ntet, &npnt, &tvrtx(0)(0), &etype, &numflag, &nparts, &edgecut,&(gbl->i1wk(0)),&(gbl->i2wk(0)));
	
	for(i=0;i<ntet;++i)
		tet(i).info = gbl->i1wk(i);
	
	gbl->i1wk = -1;
	
	return;
}
#endif

/* WHEN THIS ROUTINE EXITS */
/* THE OLD MESH STORES  */
/* pnt(pind).info = new pnt index or -1 */
/* THE NEW MESH STORES */
/* pnt(sind).info = old pnt index */
/* seg(sind).info = old side index */
/* tri(tind).info = old tri index */
void tet_mesh::partition(class tet_mesh& xin, int npart) {
	int i,j,n,tind,indx,lcl;
	Array<int,2> bcntr(xin.nfbd +10,2);
	int bnum,bel,match;
	
	/* TO CREATE UNIQUE FACE NUMBERS */
	int maxfnum = 0;
	for(i=0;i<xin.nfbd;++i) {
		maxfnum = MAX(maxfnum,xin.fbdry(i)->idnum);
	}
	++maxfnum;
	
	/* TO CREATE UNIQUE EDGE NUMBERS */
	int maxenum = 0;
	for(i=0;i<xin.nebd;++i) {
		maxenum = MAX(maxenum,xin.ebdry(i)->idnum);
	}
	++maxenum;
	
	/* TO CREATE UNIQUE VERTEX NUMBERS */
	int maxvnum = 0;
	for(i=0;i<xin.nvbd;++i) {
		maxvnum = MAX(maxvnum,xin.vbdry(i)->idnum);
	}
	++maxvnum;
	
	for(i=0;i<xin.npnt;++i)
		xin.pnt(i).info = -1;
	
	ntet = 0;
	for(i=0;i<xin.ntet;++i) {
		if (xin.tet(i).info == npart) {
			++ntet;
			for(n=0;n<4;++n) {
				xin.pnt(xin.tet(i).pnt(n)).info = npart;
			}
		}
	}
	if (!initialized) {
		maxvst = static_cast<int>(static_cast<FLT>(ntet*xin.maxvst)/xin.ntet);
		allocate(maxvst);
	}
	else if (4*ntet > maxvst) {//FIX ME 3->4??
		*gbl->log << "mesh is too small" << std::endl;
		exit(1);
	}
	
	ostringstream nstr;
	nstr << "b" << npart << std::flush;
	gbl->idprefix = nstr.str();
	gbl->idnum = npart;
	nstr.clear();
	
	/* FILL IN PNT ARRAY */
	npnt = 0;
	for(i=0;i<xin.npnt;++i) {
		if (xin.pnt(i).info == npart) {
			for(n=0;n<ND;++n)
				pnts(npnt)(n) = xin.pnts(i)(n);
			xin.pnt(i).info = npnt;
			++npnt;
		}
	}
	
	/* FILL IN TET ARRAY */
	ntet = 0;
	for(i=0;i<xin.ntet;++i) {
		if (xin.tet(i).info == npart) {
			for(n=0;n<4;++n)
				tet(ntet).pnt(n) = xin.pnt(xin.tet(i).pnt(n)).info;
			tet(ntet).info = i;
			++ntet;
		}
	}
	
	create_from_tet();
	
	nfbd = 0;
	for(i=0;i<ntri;++i) {
		
		/* tet number in old mesh */
		tind = tet(tri(i).tet(0)).info;
		
		lcl=0;
		for(n=0;n<4;++n)  // CHECK JUST USING LOCAL TIND & LOCAL PNT??? 
			lcl += xin.pnt(xin.tet(tind).pnt(n)).info;

		for(n=0;n<3;++n)
			lcl -= tri(i).pnt(n);
		
		for(n=0;n<4;++n)
			if(xin.pnt(xin.tet(tind).pnt(n)).info == lcl) break;
		
		assert(n < 4);
			
		tri(i).info = xin.tet(tind).tri(n);
				
		
		if (tri(i).tet(1) < 0) {		
			indx = xin.tet(tind).tet(n);

			if (indx < 0) {
				/* BOUNDARY TRI */
				bnum = getbdrynum(indx);
				
				for (j = 0; j <nfbd;++j) {
					if (bnum == bcntr(j,0)) {
						++bcntr(j,1);
						tri(i).tet(1) = numatbdry(j,0);
						goto next1;
					}
				}
				/* NEW TRI */
				tri(i).tet(1) = numatbdry(nfbd, 0);
				bcntr(nfbd,0) = bnum;
				bcntr(nfbd++,1) = 1;
			}
			else {
				/* PARTITION FACE */
				match = xin.tet(indx).info;
				if (match < npart) bnum = (match<<16) + (npart<<24);
				else bnum = (match<<24) + (npart<<16);
				for (j = 0; j <nfbd;++j) {
					if (bcntr(j,0) == -bnum) {
						++bcntr(j,1);
						tri(i).tet(1) = numatbdry(j,match);
						goto next1;
					}
				}
				/* NEW FACE */
				tri(i).tet(1) = numatbdry(nfbd,match);
				bcntr(nfbd,0) = -bnum;
				bcntr(nfbd++,1) = 1;
			}
		}
		next1: continue;
	}
	
	/* Make new boundaries and print boundary information */
	fbdry.resize(nfbd);
	for(i=0;i<nfbd;++i) {		
		if (bcntr(i,0) < 0) {
			fbdry(i) = new fpartition(abs(bcntr(i,0)),*this);
		}
		else {
			fbdry(i) = xin.fbdry(bcntr(i,0))->create(*this);
		}
		fbdry(i)->alloc(static_cast<int>(bcntr(i,1)*3));
		fbdry(i)->ntri = 0;
	}
		
	for(i=0;i<ntri;++i) {
		if (tri(i).tet(1) < 0) {
			bnum = getbdrynum(tri(i).tet(1));
			fbdry(bnum)->tri(fbdry(bnum)->ntri++).gindx = i;
		}
	}

	for(i = 0; i < nfbd; ++i) {
		for(tind = 0; tind < fbdry(i)->ntri;++tind) {
			int tri_gindx = fbdry(i)->tri(tind).gindx;
			for(j=0;j<3;++j)
				fbdry(i)->tri(tind).pnt(j) = tri(tri_gindx).pnt(j);		
		}
	}

	for(int i = 0; i < nfbd; ++i) 
		fbdry(i)->create_from_tri(); 
	
	/* stopped here */
//	for(i=0;i<nebd;++i) {
//		/* CREATES NEW BOUNDARY FOR DISCONNECTED SEGMENTS OF SAME TYPE */
//		ebdry(i)->reorder();
//	}
//	
//	/* MOVE VERTEX BOUNDARY INFO */
//	nvbd = 0;
//	for(i=0;i<xin.nvbd;++i)
//		if (xin.pnt(xin.vbdry(i)->pnt).info > -1)
//			++nvbd;
//	vbdry.resize(nvbd+2*nebd);
//	
//	nvbd = 0;
//	for(i=0;i<xin.nvbd;++i) {
//		if (xin.pnt(xin.vbdry(i)->pnt).info > -1) {
//			vbdry(nvbd) = xin.vbdry(i)->create(*this);
//			vbdry(nvbd)->alloc(4);
//			vbdry(nvbd)->pnt = xin.pnt(xin.vbdry(i)->pnt).info;
//			++nvbd;
//		}
//	}
//	
	/* CREATE COMMUNICATION ENDPOINT BOUNDARIES */
	for(i=0;i<nfbd;++i) {
		if (fbdry(i)->mytype == "partition") {
			/* Now that all independent fbdry sides are determined give new numbers to partitions */
			tind = fbdry(i)->tri(0).gindx;
			match = getbdryel(tri(tind).tet(1));
			if (npart < match) {
				int newid = tri(fbdry(i)->tri(0).gindx).info +maxfnum;
				ostringstream nstr;
				nstr << "b" << npart << "_f" << newid << std::flush;
				fbdry(i)->idprefix = nstr.str();
				fbdry(i)->idnum = newid;
				nstr.clear();
				std::cout << fbdry(i)->idprefix << "_type: partition\n";
			}
			else {
				int newid = tri(fbdry(i)->tri(fbdry(i)->ntri-1).gindx).info +maxfnum;
				ostringstream nstr;
				nstr << "b" << npart << "_f" << newid << std::flush;
				fbdry(i)->idprefix = nstr.str();
				fbdry(i)->idnum = newid;
				nstr.clear();	
				std::cout << fbdry(i)->idprefix << "_type: partition\n";
			}
//			
//			p0 = seg(sind).pnt(0);
//			for(j=0;j<nvbd;++j)
//				if (vbdry(j)->pnt == p0) goto nextv0;
//			
//			/* ENDPOINT IS NEW NEED TO DEFINE BOUNDARY */
//			tind = tri(seg(sind).tri(0)).info;
//			for(n=0;n<3;++n)
//				if (xin.pnt(xin.tri(tind).pnt(n)).info == p0) break;
//			vbdry(nvbd) = new vcomm(xin.tri(tind).pnt(n) +maxvnum,*this);
//			vbdry(nvbd)->alloc(4);
//			vbdry(nvbd)->pnt = p0;
//			std::cout << vbdry(nvbd)->idprefix << "_type: comm\n";
//			++nvbd;
//			
//		nextv0:
//			sind = ebdry(i)->seg(ebdry(i)->nseg-1);
//			p0 = seg(sind).pnt(1);
//			for(j=0;j<nvbd;++j)
//				if (vbdry(j)->pnt == p0) goto nextv1;
//			
//			/* NEW ENDPOINT */
//			tind = tri(seg(sind).tri(0)).info;
//			for(n=0;n<3;++n)
//				if (xin.pnt(xin.tri(tind).pnt(n)).info == p0) break;
//			vbdry(nvbd) = new vcomm(xin.tri(tind).pnt(n)+maxvnum,*this);
//			vbdry(nvbd)->alloc(4);
//			vbdry(nvbd)->pnt = p0;
//			std::cout << vbdry(nvbd)->idprefix << "_type: comm\n";
//			++nvbd;
//			
//		nextv1:
//			continue;
		}
	}
	
	bdrylabel();  // CHANGES STRI / TTRI ON BOUNDARIES TO POINT TO GROUP/ELEMENT
	
//	createtritri();
//	createpnttri();
//	cnt_nbor();
	FLT xmin[ND], xmax[ND];
	for(n=0;n<ND;++n) {
		xmin[n] = xin.otree.xmin(n);
		xmax[n] = xin.otree.xmax(n);
	}
	treeinit(xmin,xmax);
	
	initialized = 1;
	
	return;
}




