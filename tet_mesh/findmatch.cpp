#include "tet_mesh.h"
#include "tet_boundary.h"

#ifdef USING_MADLIB
#include "MAdLibInterface.h"
#endif

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

	/* ALL BOUNDARIES LOAD POSITIONS INTO BUFFERS */
	for(int i=0;i<nfbd;++i) 
		fbdry(i)->ploadbuff(boundary::all,&(pnts(0)(0)),0,ND-1,ND);
		
	/* FIRST PART OF SENDING, POST ALL RECEIVES, ALL MASTERS SEND TO SLAVES */
	for(int i=0;i<nfbd;++i)
		fbdry(i)->comm_prepare(boundary::all,0,boundary::master_slave);
						
	/* SECOND PART */
	tmsgpass(boundary::all,0,boundary::master_slave);

	/* FINAL PART OF SENDING */
	for(int i=0;i<nfbd;++i)
		fbdry(i)->comm_wait(boundary::all,0,boundary::master_slave);
	
	/* FINAL PART OF SENDING, REPLACE SLAVES BUFFER WITH MASTER BUFFER */
	for(int i=0;i<nfbd;++i)
		fbdry(i)->comm_finish(boundary::all,0,boundary::master_slave,boundary::replace);

	/* Slave receives vertex positions & matches vertex numbering */
	for(int i=0;i<nfbd;++i)
		fbdry(i)->match_numbering(1);

	
	/* ALL BOUNDARIES LOAD POSITIONS INTO BUFFERS */
	for(int i=0;i<nebd;++i) 
		ebdry(i)->ploadbuff(boundary::all,&(pnts(0)(0)),0,ND-1,ND);
	
	/* FIRST PART OF SENDING, POST ALL RECEIVES, ALL MASTERS SEND TO SLAVES */
	for(int i=0;i<nebd;++i)
		ebdry(i)->comm_prepare(boundary::all,0,boundary::master_slave);
	
	/* SECOND PART */
	for(int i=0;i<nebd;++i) 
		ebdry(i)->comm_exchange(boundary::all,0,boundary::master_slave);
	
	/* FINAL PART OF SENDING */
	for(int i=0;i<nebd;++i)
		ebdry(i)->comm_wait(boundary::all,0,boundary::master_slave);
	
	/* FINAL PART OF SENDING, REPLACE SLAVES BUFFER WITH MASTER BUFFER */
	for(int i=0;i<nebd;++i)
		ebdry(i)->comm_finish(boundary::all,0,boundary::master_slave,boundary::replace);
	
	/* Slave receives vertex positions & matches vertex numbering */
	for(int i=0;i<nebd;++i)
		ebdry(i)->match_numbering(0);
	
	/* Create global numbering system */
	create_unique_numbering();
	 
	
	
	
//	// dont think this needs to be done anymore
//	/* Reorder Side boundaries so direction is the same */
//	for(int i=0;i<nebd;++i) 
//		ebdry(i)->match_numbering(1);
//
//	/* FIRST PART OF SENDING, POST ALL RECEIVES */
//	for(int i=0;i<nebd;++i)
//		ebdry(i)->comm_prepare(boundary::all,0,boundary::master_slave);
//						
//	/* SECOND PART */
//	 for(int i=0;i<nebd;++i) 
//		ebdry(i)->comm_exchange(boundary::all,0,boundary::master_slave);
//
//	/* FINAL PART OF SENDING */
//	for(int i=0;i<nebd;++i)
//		ebdry(i)->comm_wait(boundary::all,0,boundary::master_slave);
//	
//	/* FINAL PART OF SENDING */
//	for(int i=0;i<nebd;++i)
//		ebdry(i)->comm_finish(boundary::all,0,boundary::master_slave,boundary::replace);
//	
//	for(int i=0;i<nebd;++i) 
//		ebdry(i)->match_numbering(2);		  

	
	
	
	
	
	
	
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
	
//#define TEST_MATCH
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
				*gbl->log << pnts(seg(fbdry(i)->seg(j).gindx).pnt(0)) << ' ' << pnts(seg(fbdry(i)->seg(j).gindx).pnt(1)) << std::endl;
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
				*gbl->log << pnts(seg(fbdry(i)->seg(j).gindx).pnt(0)) << ' ' << pnts(seg(fbdry(i)->seg(j).gindx).pnt(1)) << std::endl;
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

extern "C" void METIS_PartMeshDual(int *ne, int *nn, int *elmnts, int *etype, int *numflag, int *nparts, int *edgecut,
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
	
	METIS_PartMeshDual(&ntet, &npnt, &tvrtx(0)(0), &etype, &numflag, &nparts, &edgecut,&(gbl->i1wk(0)),&(gbl->i2wk(0)));
	
	for(i=0;i<ntet;++i)
		tet(i).info = gbl->i1wk(i);
	
	gbl->i1wk = -1;
	
	return;
}
#endif


void tet_mesh::partition(class tet_mesh& xin, int npart, int nparts) {
	int i,j,n,tind,indx,lcl,p0,p1,sind,egindx,gindx,newid,find;
	Array<int,2> bcntr(xin.nfbd +40,2);
	Array<int,2> ecntr(xin.nfbd+xin.nebd +40,4);
	int number_problem_edges = 0;
	Array<int,2> problem_edges(number_problem_edges,2);
	
	//TinyVector<int,2> a,b;//don't need?
	int bnum,match;
	
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
		maxvst = static_cast<int> (static_cast<FLT> (ntet) * static_cast<FLT> (xin.maxvst) / static_cast<FLT> (xin.ntet));
		allocate(maxvst*2);
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
	
	/* fill in the rest of the mesh data structure */
	create_from_tet();// messes up pnt,seg,tri.info
	
	/* find and create face boundaries */
	nfbd = 0;
	bcntr = -1;
	
	for(i=0;i<ntri;++i) {
		
		/* inside tet connected to tri */
		tind = tri(i).tet(0);
		
		/* tet.info = gbl index of tet*/
		int gbl_tind = tet(tind).info;		
		
		/* find gbl index of original tri in xin mesh */
		for(n=0;n<4;++n){	
			tri(i).info = xin.tet(gbl_tind).tri(n);
			lcl = 0;
			for (j=0; j<3; ++j) {
				lcl += xin.pnt(xin.tri(tri(i).info).pnt(j)).info;
				lcl -= tri(i).pnt(j);
			}
			
			if (lcl == 0) break;
			
		}
		assert(n < 4);
		
		/* if tri.tet(1) < 0 then it is a boundary tet in new mesh */
		if (tri(i).tet(1) < 0) {		
			
			/* find out if it is an original boundary or a new partition */
			indx = xin.tet(gbl_tind).tet(n);
			
			if (indx < 0) {
				/* BOUNDARY TRI */
				bnum = getbdrynum(indx);
				
				for (j = 0; j <nfbd;++j) {
					if (bnum == bcntr(j,0)) {
						++bcntr(j,1);
						tri(i).tet(1) = -j-1;// keep track of bnum and shift by one
						goto next1;
					}
				}
				/* FIRST TRI */
				tri(i).tet(1) = -j-1;
				bcntr(nfbd,0) = bnum;
				bcntr(nfbd++,1) = 1;
			}
			else {
				/* PARTITION FACE */
				/* find partition number of connected tet */
				match = xin.tet(indx).info;
				bnum = match+10;
				
				for (j = 0; j <nfbd;++j) {
					if (bcntr(j,0) == -bnum) {
						++bcntr(j,1);
						tri(i).tet(1) = -j-1;
						goto next1;
					}
				}
				/* FIRST FACE */
				tri(i).tet(1) = -j-1;
				bcntr(nfbd,0) = -bnum;
				bcntr(nfbd++,1) = 1;
			}
		}
	next1: continue;
	}
	
	/* Make new face boundaries and print boundary information */
	fbdry.resize(nfbd);
	for(i=0;i<nfbd;++i) {		
		if (bcntr(i,0) < 0) {
			/* create new partition with wacky idnum */
			fbdry(i) = new fpartition(maxfnum+abs(bcntr(i,0)),*this);
		}
		else {
			/* reuse face boundary in original mesh */
			fbdry(i) = xin.fbdry(bcntr(i,0))->create(*this);
		}
		fbdry(i)->alloc(static_cast<int>(bcntr(i,1)*3));
		fbdry(i)->ntri = 0;
	}
	
	/* create global indx pointer to mesh */
	for(i=0;i<ntri;++i) {
		if (tri(i).tet(1) < 0) {
			bnum = -tri(i).tet(1)-1;
			fbdry(bnum)->tri(fbdry(bnum)->ntri++).gindx = i;
		}
	}
	/* create triangle data on face boundaries */
	for(i = 0; i < nfbd; ++i) {
		for(tind = 0; tind < fbdry(i)->ntri;++tind) {
			int tri_gindx = fbdry(i)->tri(tind).gindx;
			for(j=0;j<3;++j)
				fbdry(i)->tri(tind).pnt(j) = tri(tri_gindx).pnt(j);		
		}
		/* fill in mesh data structure on face boundaries */
		fbdry(i)->create_from_tri();
	}
	
	/* find face boundaries that are disconnected and separate them */
	for(i = 0; i < nfbd; ++i) 
		fbdry(i)->pull_apart_face_boundaries();	
		
	/* find and create edge boundaries */
	nebd = 0;	
	for(i = 0; i < nseg; ++i)
		seg(i).info = -1; 
	
	Array<int,1> edge_partitions(nparts); /* list of 1's and 0's to keep track if a segment touches a partition */	
	Array<int,2> edge_boundary_list(0,nparts); /* each unique combination of partitions touching an edge boundary  */
	Array<int,2> edge_boundary_faces(0,2); /* two faces touching an edge boundary  */
	Array<int,2> edge_boundary_count(0,2); /* keeps track of boundary number and number of segs on edge boundary */

	for(i = 0; i < nfbd; ++i){
		for(int sind = 0; sind < fbdry(i)->nseg;++sind){
			if (fbdry(i)->seg(sind).tri(1) < 0){
				find = fbdry(i)->seg(sind).tri(0); /* inside fbdry triangle connected to seg */
				gindx = fbdry(i)->seg(sind).gindx; /* global index of segment */
				
				tind = fbdry(i)->tri(find).gindx; /* global index of triangle */
				
				/* find global index to old xin mesh but don't store yet */
				for(n=0;n<3;++n){					
					egindx = xin.tri(tri(tind).info).seg(n); 
					p0 = xin.pnt(xin.seg(egindx).pnt(0)).info;
					p1 = xin.pnt(xin.seg(egindx).pnt(1)).info;

					if (p0 == seg(gindx).pnt(0) && p1 == seg(gindx).pnt(1)) break;
									
					if (p1 == seg(gindx).pnt(0) && p0 == seg(gindx).pnt(1)) break;
					
				}
				assert(n < 3);
								
				/* predefined edge boundary */
				if (xin.seg(egindx).info != -1) {
					
					/* find index of xin edge boundary */
					int xin_index = -1;
					for (j = 0; j < xin.nebd; ++j) {
						if (xin.ebdry(j)->idnum == xin.seg(egindx).info) {
							xin_index = j;
							break;
						}
					}
					
					/* if it already exists then mark and add to count */
					for (n = 0; n < nebd; ++n) {
						if (edge_boundary_count(n,0) == xin_index) {
							seg(gindx).info = n;
							++edge_boundary_count(n,1);
							goto nextseg;							
						}
					}
					
					/* new edge boundary */					
					++nebd;
					edge_boundary_list.resizeAndPreserve(nebd,nparts);
					edge_boundary_count.resizeAndPreserve(nebd,2);
					edge_boundary_faces.resizeAndPreserve(nebd,2);
					
					/* put some random stuff in here */
					for (j=0; j<nparts; ++j) 
						edge_boundary_list(nebd-1,j) = 50;
					edge_boundary_faces(nebd-1,0) = nfbd+1;
					edge_boundary_faces(nebd-1,1) = nfbd+1;					
					
					seg(gindx).info = nebd-1;
					edge_boundary_count(nebd-1,1) = 1; /* first edge found */	
					
					/* find predefined edge boundary and keep track of number */
					for (j = 0; j < xin.nebd; ++j) {
						if (xin.ebdry(j)->idnum == xin.seg(egindx).info) 
							edge_boundary_count(nebd-1,0) = j;
						 
					}	
					
					goto nextseg;
				}	
				
				/* find all tets connected to edge */
				xin.ring(egindx);

				/* number of tets connected to edge */
				int nnbor = xin.seg(egindx).nnbor;

				/* if an edge is attached to a partition insert a 1 */
				edge_partitions = 0;
				for(n = 0; n < nnbor; ++n)    
					edge_partitions(xin.tet(xin.gbl->i2wk(n)).info) = 1;
								
				/* keep list of boundary numbers negative for partition */
				for(n=0;n<nebd;++n){
					/* check to see if edge boundary already exists */
					/* compare partitions */
					int diff_boundary_list = 0;
					for (j=0; j<nparts; ++j) 
						diff_boundary_list += abs(edge_boundary_list(n,j) - edge_partitions(j));
					
					/* compare face boundaries */
					// temp fix me: uncomment these two lines and get rid of -100 thing
					//diff_boundary_list += abs(edge_boundary_faces(n,0) - seg(gindx).info);
					//diff_boundary_list += abs(edge_boundary_faces(n,1) - i); // i is the face boundary

					if(diff_boundary_list == 0 && seg(gindx).info == -100){//temp fix me
					//if(diff_boundary_list == 0 ){//temp fix me
						
						/* found a partition edge set to a unique negative number */
						seg(gindx).info = n; 
						edge_boundary_count(n,0) = -n-1; 
						
						/* add seg to number of segs on edge boundary */
						++edge_boundary_count(n,1);					
						
						goto nextseg;
					}
				}
				
				/* edge found first time keep go to next seg */
				if (seg(gindx).info == -1) {
					seg(gindx).info = i;
					seg(gindx).info = -100;//temp fix me 
					goto nextseg;
				}
				
				
				/* new edge boundary */
				++nebd;
				edge_boundary_list.resizeAndPreserve(nebd,nparts);
				edge_boundary_count.resizeAndPreserve(nebd,2);
				edge_boundary_faces.resizeAndPreserve(nebd,2);
				
				for (j=0; j<nparts; ++j) 
					edge_boundary_list(nebd-1,j) = edge_partitions(j);

				/* new partition give a unique negative number */
				edge_boundary_faces(nebd-1,0) = seg(gindx).info;
				edge_boundary_faces(nebd-1,1) = i;
				seg(gindx).info = nebd-1;
				edge_boundary_count(nebd-1,0) = -nebd; 
				edge_boundary_count(nebd-1,1) = 1; /* first edge found */	

				
			}
			nextseg:;
		}
	}
	
	/* allocate edge boundaries */
	ebdry.resize(nebd);
	for(i=0;i<nebd;++i) {
		if (edge_boundary_count(i,0) < 0) {
			ebdry(i) = new epartition(maxenum+abs(edge_boundary_count(i,0)),*this);
		}
		else {
			ebdry(i) = xin.ebdry(edge_boundary_count(i,0))->create(*this);
		}
		ebdry(i)->alloc(static_cast<int>(5*edge_boundary_count(i,1)));
		ebdry(i)->nseg = 0;
	}
	
	/* create global index pointer for edge boundaries */
	for(i = 0; i < nfbd; ++i){
		for(int sind = 0; sind < fbdry(i)->nseg;++sind){
			if (fbdry(i)->seg(sind).tri(1) < 0){
				tind = fbdry(i)->seg(sind).tri(0);
				gindx = fbdry(i)->seg(sind).gindx;	
				int bdryindx = seg(gindx).info;
				
				if (bdryindx != -1) 
					ebdry(bdryindx)->seg(ebdry(bdryindx)->nseg++).gindx = gindx;
				
				seg(gindx).info = -1;
	
			}	
		}
	}
	
	/* now find pointer to seg from old mesh and store in seg.info 
	 used later to define partition numbering */
	for(sind = 0; sind < nseg; ++sind){
		tind=seg(sind).tet;
		for(n=0;n<6;++n){
			int eind = tet(tind).seg(n);			
			if(eind == sind){
				seg(sind).info=xin.tet(tet(tind).info).seg(n);
				break;
			}
		}
		assert(n<6);

		bool pointsmatch = false;
		p0 = xin.pnt(xin.seg(seg(sind).info).pnt(0)).info;
		p1 = xin.pnt(xin.seg(seg(sind).info).pnt(1)).info;
		if (p0 == seg(sind).pnt(0) && p1 == seg(sind).pnt(1)) pointsmatch = true;		
		if (p1 == seg(sind).pnt(0) && p0 == seg(sind).pnt(1)) pointsmatch = true;
		if (pointsmatch==false) {
			*gbl->log << "points dont match for edge" << endl;
			exit(2);
		}
	}

	for(i=0;i<nebd;++i) {
		
		/* put zeros in because next prev is not setup yet and need to do a swap */
		for(j=0;j<ebdry(i)->nseg;++j){
			ebdry(i)->seg(j).next = 0;
			ebdry(i)->seg(j).prev = 0;			
		}			
		
		/* start each boundary with minimum xin gbl index 
		   so loops start with the same edge on each partition*/
		int minseg = seg(ebdry(i)->seg(0).gindx).info;
		int minsegindx = 0;
		for(j=1;j<ebdry(i)->nseg;++j) {
			if (seg(ebdry(i)->seg(j).gindx).info < minseg) {
				minseg = seg(ebdry(i)->seg(j).gindx).info;
				minsegindx = j;			
			}
		}
		
		/* swap so that minsegindx is first */
		ebdry(i)->swap(minsegindx,0);
		
		/* set up all the next and prev for edge boundary */
		ebdry(i)->setup_next_prev();

		/* reorder edge boundaries to make it easy to find vertex boundaries 
		 also can separate disconnected edge boundaries */
		ebdry(i)->reorder();

	}

	
	/* MOVE VERTEX BOUNDARY INFO */
	nvbd = 0;
	for(i=0;i<xin.nvbd;++i)
		if (xin.pnt(xin.vbdry(i)->pnt).info > -1)
			++nvbd;
	
	vbdry.resize(nvbd+2*nebd);
	
	nvbd = 0;
	for(i=0;i<xin.nvbd;++i) {
		if (xin.pnt(xin.vbdry(i)->pnt).info > -1) {
			vbdry(nvbd) = xin.vbdry(i)->create(*this);
			vbdry(nvbd)->alloc(4);
			vbdry(nvbd)->pnt = xin.pnt(xin.vbdry(i)->pnt).info;
			++nvbd;
		}
	}
	
	/* CREATE COMMUNICATION FACE BOUNDARIES */
	for(i=0;i<nfbd;++i) {
		if (fbdry(i)->mytype == "partition") {
			/* Now that all independent fbdry sides are determined give new numbers to partitions */
			int minface = tri(fbdry(i)->tri(0).gindx).info;
			for(j=1;j<fbdry(i)->ntri;++j)
				minface=MIN(minface,tri(fbdry(i)->tri(j).gindx).info);
			newid = minface+maxfnum;
			ostringstream nstr;
			nstr << "b" << npart << "_f" << newid << std::flush;
			fbdry(i)->idprefix = nstr.str();
			fbdry(i)->idnum = newid;
			nstr.clear();	
			std::cout << fbdry(i)->idprefix << "_type: comm\n";
		}
	}
	
	/* CREATE COMMUNICATION EDGE BOUNDARIES */
	for(i=0;i<nebd;++i) {
		if (ebdry(i)->mytype == "partition") {
			/* Now that all independent ebdry sides are determined give new numbers to partitions */
			int minseg = seg(ebdry(i)->seg(0).gindx).info;
			for(j=1;j<ebdry(i)->nseg;++j)
				minseg=MIN(minseg,seg(ebdry(i)->seg(j).gindx).info);
			newid = minseg + maxenum;
			ostringstream nstr;
			nstr << "b" << npart << "_e" << newid << std::flush;
			ebdry(i)->idprefix = nstr.str();
			ebdry(i)->idnum = newid;
			nstr.clear();
			std::cout << ebdry(i)->idprefix << "_type: comm\n";
		}
	}
	
	/* CREATE COMMUNICATION VERTEX BOUNDARIES */
	for(i=0;i<nebd;++i) {
		sind = ebdry(i)->seg(0).gindx;	
		p0 = seg(sind).pnt(0);
		
		for(j=0;j<nvbd;++j)
			if (vbdry(j)->pnt == p0) goto nextv0;
				
		/* ENDPOINT IS NEW NEED TO DEFINE BOUNDARY */
		tind = tet(seg(sind).tet).info;
		for(n=0;n<4;++n)
			if (xin.pnt(xin.tet(tind).pnt(n)).info == p0) break;
		vbdry(nvbd) = new vcomm(xin.tet(tind).pnt(n) +maxvnum,*this);
		vbdry(nvbd)->alloc(4);
		vbdry(nvbd)->pnt = p0;
		std::cout << vbdry(nvbd)->idprefix << "_type: comm\n";
		++nvbd;
		
	nextv0:
		
		sind = ebdry(i)->seg(ebdry(i)->nseg-1).gindx;
		p0 = seg(sind).pnt(1);
		
		for(j=0;j<nvbd;++j)
			if (vbdry(j)->pnt == p0) goto nextv1;
		
		/* NEW ENDPOINT */
		tind = tet(seg(sind).tet).info;
		for(n=0;n<4;++n)
			if (xin.pnt(xin.tet(tind).pnt(n)).info == p0) break;
		vbdry(nvbd) = new vcomm(xin.tet(tind).pnt(n)+maxvnum,*this);
		vbdry(nvbd)->alloc(4);
		vbdry(nvbd)->pnt = p0;
		std::cout << vbdry(nvbd)->idprefix << "_type: comm\n";
		++nvbd;
		
	nextv1:
		continue;
		
	}
	
#define findlone

#ifdef findlone 
	/* find lone pnts and edges that may need to be set to dirichlet */
	Array<int,1> seginfo(xin.nseg);
	
	for(int i=0;i<xin.nseg;++i)
		seginfo(i) = -1;
	
	for(sind = 0; sind < nseg; ++sind)
		seginfo(seg(sind).info) = sind;
	
	for (int i=0; i<npnt; ++i) 
		pnt(i).info = 0;	
	
	for (int i=0; i<nseg; ++i) 
		seg(i).info = 0;
	
	for (int i=0; i < nfbd; ++i) {
		if (fbdry(i)->mytype != "partition") {
			for (int j=0; j<fbdry(i)->npnt; ++j) 
				pnt(fbdry(i)->pnt(j).gindx).info = 1;
			
			for (int j=0; j<fbdry(i)->nseg; ++j) 
				seg(fbdry(i)->seg(j).gindx).info = 1;
			
		}
	}

	for (int k=0; k < xin.nfbd; ++k) {
		for (int n=0; n< xin.fbdry(k)->npnt; ++n) {
			p0 = xin.pnt(xin.fbdry(k)->pnt(n).gindx).info;
			if (p0 != -1 && pnt(p0).info == 0) {
				pnt(p0).info = 1;
				std::cout << "# found a lone pnt idnum: " << xin.fbdry(k)->pnt(n).gindx+maxvnum << " point: " << p0 << endl;
				std::cout << "# should be set to same BC as face boundary idnum: " << xin.fbdry(k)->idnum << endl;
				
//					vbdry.resizeAndPreserve(nvbd+1);
//					vbdry(nvbd) = new vcomm(xin.fbdry(k)->pnt(n).gindx+maxvnum,*this);
//					vbdry(nvbd)->alloc(4);
//					vbdry(nvbd)->pnt = p0;
//					std::cout << vbdry(nvbd)->idprefix << "_cd_type: dirichlet\n";
//					++nvbd;
				
			}
		}
		
		/* not sure if this is correct */
		for (int n=0; n< xin.fbdry(k)->nseg; ++n) {
			p0 = seginfo(xin.fbdry(k)->seg(n).gindx);
			if (p0 != -1 && seg(p0).info == 0) {
				seg(p0).info = 1;
				std::cout << "# found a lone seg idnum: " << xin.fbdry(k)->seg(n).gindx+maxenum << " seg:" <<  p0 << endl;//<< xin.fbdry(k)->pnt(n).gindx+maxvnum << " point: " << p0 << endl;
				std::cout << "# should be set to same BC as face boundary idnum: " << xin.fbdry(k)->idnum << endl;	
//					ebdry.resizeAndPreserve(nebd+1);				
//					ebdry(nebd) = new epartition(maxenum+xin.fbdry(k)->seg(n).gindx,*this);
//					ebdry(nebd)->alloc(static_cast<int>(5));
//					ebdry(nebd)->nseg = 0;
//					ebdry(nebd)->seg(0).gindx = xin.seg(xin.fbdry(k)->seg(n).gindx).info;
//					ebdry(nebd)->setup_next_prev();
//					ebdry(nebd)->reorder();
//					std::cout << ebdry(nebd)->idprefix << "_cd_type: dirichlet\n";
//					++nebd;
				
			}
		}
	}
	
#endif

	
	/* call match_all because reorder doesnt account for sign change */
	match_all();
	for(int i = 0; i < nfbd; ++i){
		fbdry(i)->match_all(); 
	}
	
	bdrylabel();  // CHANGES STRI / TTRI ON BOUNDARIES TO POINT TO GROUP/ELEMENT
	
	FLT xmin[ND], xmax[ND];
	for(n=0;n<ND;++n) {
		xmin[n] = xin.otree.xmin(n);
		xmax[n] = xin.otree.xmax(n);
	}
	treeinit(xmin,xmax);
	
	initialized = 1;
	
	return;
}

void tet_mesh::setup_partition(int nparts, Array<int,2> & boundary_partitions, Array<int,1> & commbdrys) {
	
	int ncommbdry = 0;
	
	boundary_partitions.resize(ncommbdry,nparts);
	commbdrys.resize(3);
	commbdrys = 0;
	
	Array<int,1> temp(nparts);

	/* fill in all infos with boundary index */
	for(int i = 0; i < npnt; ++i)
		pnt(i).info = -1;
	
	for(int i = 0; i < nvbd; ++i)
		pnt(vbdry(i)->pnt).info = i;
	
	for(int i = 0; i < nseg; ++i)
		seg(i).info = -1;
	
	for(int i = 0; i < nebd; ++i)
		for(int j = 0; j < ebdry(i)->nseg; ++j)
			seg(ebdry(i)->seg(j).gindx).info = i;
	
	for(int i = 0; i < ntri; ++i)
		tri(i).info = -1;
	
	for(int i = 0; i < nfbd; ++i)
		for(int j = 0; j < fbdry(i)->ntri; ++j)
			tri(fbdry(i)->tri(j).gindx).info = i;
	
	/* find communication face boundaries */
	for(int i = 0; i < ntri; ++i){

		if(tri(i).tet(1) > -1){			
			
			/* tag partitions with 1's */
			temp = 0;
			for(int j = 0; j < 2; ++j) {
				temp(tet(tri(i).tet(j)).info) = 1;
			}
			
			/* if sum is greater than one then it is a communication boundary */
			int lcl = 0;
			for(int j = 0; j < nparts; ++j)
				lcl += temp(j);
			
			if(lcl > 1) {	
				
				/* check if boundary already defined */
				int match = -1;
				for(int j = 0; j < ncommbdry; ++j){
					match = 0;
					for(int k = 0; k < nparts; ++k)
						match += abs(temp(k)-boundary_partitions(j,k));
					
					/* boundary found tag with unique number */
					if(match == 0) {
						tri(i).info = j+nfbd;
						break;
					}				
				}
				
				/* new face boundary */
				if(match != 0){
					/* tag with unique number */
					tri(i).info = ncommbdry+nfbd;
					/* keep track of total boundaries */
					++ncommbdry;
					/* keep track of face boundaries */
					++commbdrys(2);

					boundary_partitions.resizeAndPreserve(ncommbdry,nparts);
					/* store 0,1's */
					for(int k = 0; k < nparts; ++k)
						boundary_partitions(ncommbdry-1,k) = temp(k);				
				}
			}
		}		
	}
	
		
	/* find communication edge boundaries */
	for(int i = 0; i < nseg; ++i){
		
		/* search all tets surrounding edge */
		int nbor = seg(i).nnbor; 		
		ring(i);
		
		/* tag partitions with 1's */
		temp = 0;
		for(int j = 0; j < nbor; ++j){
			temp(tet(gbl->i2wk(j)).info) = 1;
		}
		
		/* if sum is greater than one then it is a communication boundary */
		int lcl = 0;
		for(int j = 0; j < nparts; ++j)
			lcl += temp(j);
		
		if(lcl > 1) {
			/* check if boundary already defined */
			int match = -1;
			for(int j = 0; j < ncommbdry; ++j){
				match = 0;
				for(int k = 0; k < nparts; ++k)
					match += abs(temp(k)-boundary_partitions(j,k));
				
				/* boundary found tag with unique number */
				if(match == 0) {
					seg(i).info = j+nebd;
					break;
				}				
			}

			/* new edge boundary */
			if(match != 0){
				/* tag with unique number */
				seg(i).info = ncommbdry+nebd;
				/* keep track of total communication boundaries */
				++ncommbdry;
				/* keep track of total edge boundaries */
				++commbdrys(1);
				boundary_partitions.resizeAndPreserve(ncommbdry,nparts);
				for(int k = 0; k < nparts; ++k)
					boundary_partitions(ncommbdry-1,k) = temp(k);				
			}
		}
	}
	
	/* find communication vertex boundaries */
	for(int i = 0; i < npnt; ++i){
		
		/* search tets surrounding vertex */
		int nbor = pnt(i).nnbor; 		
		vertexball(i);
		
		/* tag partitions with 1's */
		temp = 0;
		for(int j = 0; j < nbor; ++j)
			temp(tet(gbl->i2wk(j)).info) = 1;
		
		/* if sum is greater than one then it is a communication boundary */
		int lcl = 0;
		for(int j = 0; j < nparts; ++j)
			lcl += temp(j);
		
		if(lcl > 1) {
			/* check if boundary already defined */
			int match = -1;
			for(int j = 0; j < ncommbdry; ++j){
				match = 0;
				for(int k = 0; k < nparts; ++k)
					match += abs(temp(k)-boundary_partitions(j,k));

				/* boundary found tag with unique number */
				if(match == 0) {
					pnt(i).info = j+nvbd;
					break;
				}				
			}
			
			/* new vertex boundary */
			if(match != 0){
				/* tag with unique number */
				pnt(i).info = ncommbdry+nvbd;
				/* keep track of total communication boundaries */
				++ncommbdry;
				/* keep track of total vertex boundaries */
				++commbdrys(0);
				boundary_partitions.resizeAndPreserve(ncommbdry,nparts);
				for(int k = 0; k < nparts; ++k)
					boundary_partitions(ncommbdry-1,k) = temp(k);				
			}
		}
	}


	/* find comm pnts that are on more than one edge boundaries */
	Array<int,1> findpnts(npnt); findpnts = -1;	
	
//	Array<int,1> storeinfo(npnt);
//	for(int i = 0; i < npnt; ++i)
//		storeinfo(i) = pnt(i).info;

	for(int i = 0; i < nseg; ++i){
		int seginfo = seg(i).info;
		if(seginfo >= nebd+commbdrys(2)){
			for(int j = 0; j < 2; ++j) {
				int p0 = seg(i).pnt(j);
				if(pnt(p0).info >= nvbd+commbdrys(2) && pnt(p0).info < nvbd+commbdrys(1)+commbdrys(2)) {
					if(findpnts(p0) != seginfo && findpnts(p0) > -1 ) {
						findpnts(p0) = -2;
						//pnt(p0).info = ncommbdry+nvbd;
						pnt(p0).info = -2;
					}
					if(findpnts(p0) == -1) {
						findpnts(p0) = seginfo;
					}
				}
			}			
		}
	}
	
//	for(int i = 0; i < npnt; ++i)
//		if(storeinfo(i) != pnt(i).info)
//			cout << "change to pnt info " << i << ' ' << storeinfo(i) << ' ' << pnt(i).info << endl;

	
			   
	return;
}

void tet_mesh::partition2(class tet_mesh& xin, int npart, int nparts, Array<int,2> boundary_list, Array<int,1> commbdrys) {
	
	/* TO CREATE UNIQUE FACE NUMBERS */
	int maxfnum = 0;
	for(int i=0;i<xin.nfbd;++i) {
		maxfnum = MAX(maxfnum,xin.fbdry(i)->idnum);
	} ++maxfnum;
	
	/* TO CREATE UNIQUE EDGE NUMBERS */
	int maxenum = 0;
	for(int i=0;i<xin.nebd;++i) {
		maxenum = MAX(maxenum,xin.ebdry(i)->idnum);
	} ++maxenum;
	
	/* TO CREATE UNIQUE VERTEX NUMBERS */
	int maxvnum = 0;
	for(int i=0;i<xin.nvbd;++i) {
		maxvnum = MAX(maxvnum,xin.vbdry(i)->idnum);
	} ++maxvnum;

	/* work array intwk(xin pnt index) = local pnt index */
	Array<int,1> intwk(xin.npnt); intwk = -1;
	
	/* count tets in partition and tag vertex */
	ntet = 0;
	for(int i=0;i<xin.ntet;++i) {
		if (xin.tet(i).info == npart) {
			++ntet;
			for(int n=0;n<4;++n) {
				intwk(xin.tet(i).pnt(n)) = npart;
			}
		}
	}
	
	/* allocate mesh */
	if (!initialized) {
		maxvst = static_cast<int> (static_cast<FLT> (ntet) * static_cast<FLT> (xin.maxvst) / static_cast<FLT> (xin.ntet));
		allocate(maxvst*2);
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
	for(int i=0;i<xin.npnt;++i) {
		if (intwk(i) == npart) {
			for(int n=0;n<ND;++n)
				pnts(npnt)(n) = xin.pnts(i)(n);
			intwk(i) = npnt;
			++npnt;
		}
	}
	
	/* FILL IN TET ARRAY */
	ntet = 0;
	for(int i=0;i<xin.ntet;++i) {
		if (xin.tet(i).info == npart) {
			for(int n=0;n<4;++n)
				tet(ntet).pnt(n) = intwk(xin.tet(i).pnt(n));
			tet(ntet).info = i;
			++ntet;
		}
	}
	
	/**********************************/
	/*             FACES              */
	/**********************************/
	
	nfbd = 0; ntri = 0;
	Array<int,2> face_boundaries(nfbd,2);

	/* find all face boundaries */
	for(int i = 0; i < xin.ntri; ++i) {
		if(xin.tri(i).info != -1){
			
			/* find out if boundary tri is connected to this partition */
			bool partitionfound = false;
			
			if(xin.tri(i).tet(1) > -1) {
				if(xin.tet(xin.tri(i).tet(0)).info == npart || xin.tet(xin.tri(i).tet(1)).info == npart ) 
					partitionfound = true;
			}
			else {
				if(xin.tet(xin.tri(i).tet(0)).info == npart ){
					partitionfound = true;
				}
			}

			/* tri is located on this partition */
			if(partitionfound == true) {
				
				/* triangle vertex points */
				for(int n = 0; n < 3; ++n)
					tri(ntri).pnt(n) = intwk(xin.tri(i).pnt(n));
				
				/* store global index from xin mesh */
				tri(ntri).info = i;
				++ntri;
				
				
				/* keep track of unique tag and count number of tri's on each boundary */
				bool match_found = false;
				for(int j = 0; j < nfbd; ++j){
					if(xin.tri(i).info == face_boundaries(j,0)) {
						++face_boundaries(j,1);
						match_found = true;
						break;
					}					   
				}
				
				if(match_found == false) {
					face_boundaries.resizeAndPreserve(++nfbd,2);
					face_boundaries(nfbd-1,0) = xin.tri(i).info;
					face_boundaries(nfbd-1,1) = 1;
				}
			}
		}			
	}
		
	/* Make new face boundaries and print boundary information */
	fbdry.resize(nfbd);
	for(int i=0;i<nfbd;++i) {		
		if (face_boundaries(i,0) >= xin.nfbd) {
			/* create new partition */
			fbdry(i) = new fpartition(maxfnum+face_boundaries(i,0),*this);
		}
		else {
			/* reuse face boundary in original mesh */
			fbdry(i) = xin.fbdry(face_boundaries(i,0))->create(*this);
		}
		fbdry(i)->alloc(static_cast<int>(face_boundaries(i,1)*3));
		fbdry(i)->ntri = 0;
	}
	
	/* fill in global index for face boundary */
	for(int i = 0; i < ntri; ++i){
		for(int j = 0; j < nfbd; ++j){
			if(xin.tri(tri(i).info).info == face_boundaries(j,0)){								
				fbdry(j)->tri(fbdry(j)->ntri).gindx = i;				
				++fbdry(j)->ntri;
				break;
			}
		}
	}
	
	/* CREATE COMMUNICATION FACE BOUNDARIES */
	for(int i=0;i<nfbd;++i) {
		if (fbdry(i)->mytype == "partition") {
			/* now that all independent fbdry sides are determined give new numbers to partitions */
			int minface = tri(fbdry(i)->tri(0).gindx).info;
			for(int j=1;j<fbdry(i)->ntri;++j)
				minface=MIN(minface,tri(fbdry(i)->tri(j).gindx).info);
			int newid = minface+maxfnum;
			ostringstream nstr;
			nstr << "b" << npart << "_f" << newid << std::flush;
			fbdry(i)->idprefix = nstr.str();
			fbdry(i)->idnum = newid;
			nstr.clear();	
			std::cout << fbdry(i)->idprefix << "_type: comm\n";
		}
	}
	
	/**********************************/
	/*             EDGES              */
	/**********************************/
	
	nebd = 0; nseg = 0;
	Array<int,2> edge_boundaries(nebd,2);

	/* find all edge boundaries */
	for(int i = 0; i < xin.nseg; ++i) {
		/* edge boundary index */
		int seginfo = xin.seg(i).info;
		
		/* boundary edge */
		if(seginfo != -1){
			
			bool partitionfound = false;

			/* predefined non-partition edge */
			if(seginfo < xin.nebd) {
				
				/* search if on partition */
				int nbor = xin.seg(i).nnbor; 		
				xin.ring(i);				
				
				for(int j = 0; j < nbor; ++j){
					if(npart == xin.tet(xin.gbl->i2wk(j)).info) {
						partitionfound = true;
						break;
					}
				}
			}
		
			/* only accept partition segments that have combinations other than face boundaries */
			if(seginfo >= xin.nebd+commbdrys(2)){
				if(boundary_list(seginfo-xin.nebd,npart) == 1) {
					partitionfound = true;
				}				
			}
		
			/* segment is on this partition */
			if(partitionfound == true) {
				/* fill in seg.pnt info */
				for(int n = 0; n < 2; ++n)
					seg(nseg).pnt(n) = intwk(xin.seg(i).pnt(n));
				
				/* store index from xin mesh */
				seg(nseg).info = i;
				
				/* add segment to mesh */
				++nseg;
				
				/* check if edge boundary is already created */
				bool match_found = false;
				for(int j = 0; j < nebd; ++j){
					/* edge boundary found */
					if(seginfo == edge_boundaries(j,0)) {
						/* keep track of total segments on edge boundary */
						++edge_boundaries(j,1);
						match_found = true;
						break;
					}					   
				}
				
				/* if edge boundary not found then create it */
				if(match_found == false) {
					edge_boundaries.resizeAndPreserve(++nebd,2);
					/* store boundary index */
					edge_boundaries(nebd-1,0) = seginfo;
					/* first segment on edge boundary */
					edge_boundaries(nebd-1,1) = 1;
				}
			}
			
		}
					
	}
	
	/* Make new edge boundaries and print boundary information */
	ebdry.resize(nebd);
	for(int i=0;i<nebd;++i) {		
		if (edge_boundaries(i,0) >= xin.nebd) {
			/* create new partition with unique idnum */
			ebdry(i) = new epartition(maxenum+maxfnum+edge_boundaries(i,0),*this);
		}
		else {
			/* reuse face boundary in original mesh */
			ebdry(i) = xin.ebdry(edge_boundaries(i,0))->create(*this);
		}
		ebdry(i)->alloc(static_cast<int>(edge_boundaries(i,1)*3+5));
		ebdry(i)->nseg = 0;
	}
	
	/* now that edge boundary is allocated fill in global index */
	for(int i = 0; i < nseg; ++i){
		for(int j = 0; j < nebd; ++j){
			if(xin.seg(seg(i).info).info == edge_boundaries(j,0)){
				ebdry(j)->seg(ebdry(j)->nseg).gindx = i;
				++ebdry(j)->nseg;
				break;
			}
		}
	}
	
	/* fill in edge boundary stuff and reorder/separate disconnected boundaries */
	/* don't mess with predifined plain edges because gmsh will recombine them*/
	for(int i=0;i<nebd;++i) {
		if(ebdry(i)->mytype == "partition"){
			ebdry(i)->setup_next_prev();
			ebdry(i)->reorder();
		}
	}
	
	/* set intwk so intwkseg(xin seg index) = local index */	
	Array<int,1> intwkseg(xin.nseg); intwkseg = -1;
	for(int i = 0; i < nseg; ++i) {
		intwkseg(seg(i).info) = i;
	}
	
	/* tag edge communication boundaries 
	 > -1 comm, -1 dirichlet */
	Array<int,1> tagseg(nseg); tagseg = -1;
	for(int i=0;i<nebd;++i) {
		if (ebdry(i)->mytype == "partition") {
			for(int j = 0; j < ebdry(i)->nseg; ++j) {
				tagseg(ebdry(i)->seg(j).gindx) = i;
			}
		}
	}

	/* check for lone segments that need non communication boundaries */
	for(int i = 0; i < xin.nfbd; ++i) {
		for(int j = 0; j < xin.fbdry(i)->nseg; ++j) {
			int gindx = xin.fbdry(i)->seg(j).gindx;
			int e0 = intwkseg(gindx);
			
			/* edge not in mesh */
			if(e0 == -1) {					
				continue;
			}
			
			if(tagseg(e0) > -1) {
				std::cout << "# found lone segment on face boundary: " <<  xin.fbdry(i)->idnum << endl;
				
				int lclindex = -1;
				for(int k = 0; k < ebdry(tagseg(e0))->nseg; ++k){
					if(e0 == ebdry(tagseg(e0))->seg(k).gindx){
						lclindex = k;
						break;
					}
				}
				if(lclindex == -1) cout << "we got problems" << endl;
				
				if(ebdry(tagseg(e0))->nseg == 1) {
					std::cout << ebdry(tagseg(e0))->idprefix << "_type: plain" << endl;
					std::cout << ebdry(tagseg(e0))->idprefix << "_cns_type: inflow" << endl;
					continue;
				}
				
				/* copy edge boundary */
				ebdry.resizeAndPreserve(nebd+1);
				ebdry(nebd) = ebdry(tagseg(e0))->create(*this);
				ebdry(nebd)->copy(*ebdry(tagseg(e0)));	
				/* put lone edge in 0 index */
				ebdry(nebd)->swap(0,lclindex);
				/* set number of segs to 1 */
				ebdry(nebd)->nseg = 1;
				/* give unique number */
				ebdry(nebd)->idnum = gindx+maxenum;
				/* change type */
				ebdry(nebd)->mytype = "plain";
				
				/* output stuff */
				ostringstream nstr;
				nstr << "b" << npart << "_e" << ebdry(nebd)->idnum << std::flush;
				ebdry(nebd)->idprefix = nstr.str();
				nstr.clear();
				
				std::cout << ebdry(nebd)->idprefix << "_type: plain" << endl;
				std::cout << ebdry(nebd)->idprefix << "_cns_type: inflow" << endl;
				
				++nebd;

				/* swap edge with last spot */
				ebdry(tagseg(e0))->swap(lclindex,ebdry(tagseg(e0))->nseg-1);
				/* subtract one off nseg */
				--ebdry(tagseg(e0))->nseg;
				/* now reorder */
				ebdry(tagseg(e0))->setup_next_prev();
				ebdry(tagseg(e0))->reorder();
				tagseg(e0) = -1;
				
			}	
		}		
	}
	
	/* CREATE COMMUNICATION EDGE BOUNDARIES */
	for(int i=0;i<nebd;++i) {
		if (ebdry(i)->mytype == "partition") {
			/* Now that all independent ebdry sides are determined give new numbers to partitions */
			int minseg = seg(ebdry(i)->seg(0).gindx).info;
			for(int j=1;j<ebdry(i)->nseg;++j)
				minseg=MIN(minseg,seg(ebdry(i)->seg(j).gindx).info);
			int newid = minseg + maxenum;
			ostringstream nstr;
			nstr << "b" << npart << "_e" << newid << std::flush;
			ebdry(i)->idprefix = nstr.str();
			ebdry(i)->idnum = newid;
			nstr.clear();
			std::cout << ebdry(i)->idprefix << "_type: comm\n";
		}
	}
		
	/**********************************/
	/*           VERTICES             */
	/**********************************/

	/* count vertex boundaries */
	nvbd = 0;
	for(int i=0;i<xin.nvbd;++i){
		if (intwk(xin.vbdry(i)->pnt) > -1){
			++nvbd;
		}
	}

	/* MOVE VERTEX BOUNDARY INFO */
	vbdry.resize(nvbd);
	nvbd = 0;
	for(int i=0;i<xin.nvbd;++i) {
		if (intwk(xin.vbdry(i)->pnt) > -1) {
			vbdry(nvbd) = xin.vbdry(i)->create(*this);
			vbdry(nvbd)->alloc(4);
			vbdry(nvbd)->pnt = intwk(xin.vbdry(i)->pnt);
			++nvbd;
		}
	}

	/* use tagpnt to store vertex boundary index 
	   used later in lone point finder  
	   > -1 comm, -1 vertex/edge/face dirichlet, -2 neither */
	Array<int,1> tagpnt(npnt); tagpnt = -2;	
	
	/* tag pnts on dirichlet boundaries with -1 */
	for(int i = 0; i < nfbd; ++i) {
		if(fbdry(i)->mytype != "partition") {
			for(int j = 0; j < fbdry(i)->ntri; ++j){
				int gindx = fbdry(i)->tri(j).gindx;
				for(int k = 0; k < 3; ++k){
					tagpnt(tri(gindx).pnt(k)) = -1;
				}
			}
		}
	}

	for(int i = 0; i < nebd; ++i) {
		if(ebdry(i)->mytype != "partition") {
			for(int j = 0; j < ebdry(i)->nseg; ++j){
				int gindx = ebdry(i)->seg(j).gindx;
				for(int k = 0; k < 2; ++k){
					tagpnt(seg(gindx).pnt(k)) = -1;
				}
			}
		}
	}
	for(int i = 0; i < nvbd; ++i) {
		tagpnt(vbdry(i)->pnt) = -1;
	}

	/* find comm vertex boundaries and tag with > -1 */
	for(int i = 0; i < xin.npnt; ++i) {
		int pntinfo = xin.pnt(i).info;
		
		/* experimental: try and find partition edges that share a vertex */
		if(pntinfo == -2) {
			
			xin.vertexball(i);
			for(int j = 0; j < xin.pnt(i).nnbor; ++j) {
				if(xin.tet(xin.gbl->i2wk(j)).info == npart) {
				
					vbdry.resizeAndPreserve(nvbd+1);
					tagpnt(intwk(i)) = nvbd;
					vbdry(nvbd) = new vcomm(i+maxvnum,*this);
					vbdry(nvbd)->alloc(4);
					vbdry(nvbd)->pnt = intwk(i);
					vbdry(nvbd)->mytype = "comm";
					++nvbd;
					break;
				}
			}
			continue;
		}
		
		/* check if point is a unique combination of partitions 
		 other than face and edge communication boundaries*/
		if(pntinfo >= xin.nvbd+commbdrys(1)+commbdrys(2)){
			if(boundary_list(pntinfo-xin.nvbd,npart) == 1) {

				vbdry.resizeAndPreserve(nvbd+1);
				tagpnt(intwk(i)) = nvbd;
				vbdry(nvbd) = new vcomm(i+maxvnum,*this);
				vbdry(nvbd)->alloc(4);
				vbdry(nvbd)->pnt = intwk(i);
				vbdry(nvbd)->mytype = "comm";
				++nvbd;
			}				
		}
	}
	
	/* check for lone points that need non communication boundaries */
	for(int i = 0; i < xin.nfbd; ++i) {
		for(int j = 0; j < xin.fbdry(i)->npnt; ++j) {
			int p0 = intwk(xin.fbdry(i)->pnt(j).gindx);
			
			if(p0 == -1) continue;

			if(tagpnt(p0) > -1) {
				vbdry(tagpnt(p0))->mytype = "plain";
				std::cout << "# found lone point type 1 on face boundary: " <<  xin.fbdry(i)->idnum << endl;
				std::cout << vbdry(tagpnt(p0))->idprefix << "_type: plain" << endl;
				std::cout << vbdry(tagpnt(p0))->idprefix << "_cns_type: inflow" << endl;
				tagpnt(p0) = -1;
			}	
			else if(tagpnt(p0) < -1) {
				cout <<"# found lone point point type 2 on face boundary: " << xin.fbdry(i)->idnum << endl;
				vbdry.resizeAndPreserve(nvbd+1);
				vbdry(nvbd)->mytype = "plain";
				vbdry(nvbd) = new vcomm(xin.fbdry(i)->pnt(j).gindx+maxvnum,*this);
				vbdry(nvbd)->alloc(4);
				vbdry(nvbd)->pnt = p0;
				std::cout << vbdry(nvbd)->idprefix << "_type: plain\n";
				std::cout << vbdry(nvbd)->idprefix << "_cns_type: inflow\n";
				tagpnt(p0) = -1;
				++nvbd;
			}
		}
	}

	/* pnt( partition pnt index) = xin pnt index */
	for(int i = 0; i < xin.npnt; ++i)
		pnt(intwk(i)).info = i;
	
	/* CREATE COMMUNICATION VERTEX BOUNDARIES */
	for(int i=0;i<nebd;++i) {
		if(ebdry(i)->mytype == "partition"){
			int sind = ebdry(i)->seg(0).gindx;	

			int p0 = seg(sind).pnt(0);
			
			for(int j=0;j<nvbd;++j)
				if (vbdry(j)->pnt == p0) goto nextv0;
			
			/* ENDPOINT IS NEW NEED TO DEFINE BOUNDARY */
			vbdry.resizeAndPreserve(nvbd+1);
			vbdry(nvbd) = new vcomm(pnt(p0).info+maxvnum,*this);
			vbdry(nvbd)->alloc(4);
			tagpnt(p0) = nvbd;
			vbdry(nvbd)->pnt = p0;
			vbdry(nvbd)->mytype = "comm";
			++nvbd;

		nextv0:

			sind = ebdry(i)->seg(ebdry(i)->nseg-1).gindx;
			p0 = seg(sind).pnt(1);

			for(int j=0;j<nvbd;++j)
				if (vbdry(j)->pnt == p0) goto nextv1;
			
			/* NEW ENDPOINT */
			vbdry.resizeAndPreserve(nvbd+1);
			vbdry(nvbd) = new vcomm(pnt(p0).info+maxvnum,*this);
			vbdry(nvbd)->alloc(4);
			tagpnt(p0) = nvbd;
			vbdry(nvbd)->pnt = p0;
			vbdry(nvbd)->mytype = "comm";
			++nvbd;

		nextv1:
			
			continue;
		}
		
	}

	/* find pnts that are on two edge boundaries */
	Array<int,1> findpnts(npnt); findpnts = -1;	
	for(int i = 0; i < nebd; ++i){
		if(ebdry(i)->mytype == "partition") {
			int p0 = seg(ebdry(i)->seg(0).gindx).pnt(0);
			
			/* vertex already found and not on this ebdry */
			if (findpnts(p0) > -1 && findpnts(p0) != i){
				findpnts(p0) = -2;
			}

			/* first time vertex found */
			if(findpnts(p0) == -1)
				findpnts(p0) = i;
				
			for(int sind = 0; sind < ebdry(i)->nseg; ++sind){
				p0 = seg(ebdry(i)->seg(sind).gindx).pnt(1);
				/* vertex already found and not on this ebdry */
				if (findpnts(p0) > -1 && findpnts(p0) != i){
					findpnts(p0) = -2;
				}
				
				/* first time vertex found */
				if(findpnts(p0) == -1)
					findpnts(p0) = i;
			}	
		}
	}

	for(int i = 0; i < nvbd; ++i)
		findpnts(vbdry(i)->pnt) = -1;
	
	for(int i = 0; i < npnt; ++i) {
		if(findpnts(i) == -2) {
			cout << "HEY LOOK HERE make new vertex comm: " << pnt(i).info+maxvnum << endl;
			cout << "code is gonna crash with this mesh" << endl;
			//exit(1);
//			vbdry.resizeAndPreserve(nvbd+1);
//			vbdry(nvbd) = new vcomm(pnt(i).info+maxvnum,*this);
//			vbdry(nvbd)->alloc(4);
//			tagpnt(i) = nvbd;
//			vbdry(nvbd)->pnt = i;
//			std::cout << vbdry(nvbd)->idprefix << "_type: comm\n";
//			++nvbd;
		}
	}
	
	/* check for lone points that need non communication boundaries */
	for(int i = 0; i < xin.nfbd; ++i) {
		for(int j = 0; j < xin.fbdry(i)->npnt; ++j) {
			int p0 = intwk(xin.fbdry(i)->pnt(j).gindx);
			
			if(p0 == -1) continue;
			
			if(tagpnt(p0) > -1) {
				vbdry(tagpnt(p0))->mytype = "plain";				
				std::cout << "# found lone point type 3 on face boundary: " <<  xin.fbdry(i)->idnum << endl;
				std::cout << vbdry(tagpnt(p0))->idprefix << "_type: plain" << endl;
				std::cout << vbdry(tagpnt(p0))->idprefix << "_cns_type: inflow" << endl;
				tagpnt(p0) = -1;
			}	
		}
	}

	/* list any vertex boundaries that are communication */
	for(int i = 0; i < nvbd; ++i){
		if(vbdry(i)->mytype == "comm"){
			std::cout << vbdry(i)->idprefix << "_type: comm" << endl;

		}
	}
			
	
	MAd::pGModel MAdModel = NULL;
	GM_create(&MAdModel,"theModel");
	MAdLibInterface::exportToMAdModel(this, MAdModel);
	MAd::pMesh MAdMesh = M_new(MAdModel);
	MAdLibInterface::exportToMAdMesh(this, MAdMesh);	
	MAdLibInterface::importFromMAdMesh(MAdMesh,this);

	delete MAdMesh;
	delete MAdModel;
	

	return;
}
