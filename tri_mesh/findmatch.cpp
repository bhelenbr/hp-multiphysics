#include "tri_mesh.h"
#include "tri_boundary.h"
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
			p0 = seg(ebdry(i)->seg(0)).pnt(0);
			v0id = -1;
			for(j=0;j<nvbd;++j) {
				if (vbdry(j)->pnt == p0) {
					v0id = vbdry(j)->idnum;
					break;
				}
			}
			list(tsize++) = v0id;
			p0 = seg(ebdry(i)->seg(ebdry(i)->nseg-1)]).pnt(1);
			v0id = -1;
			for(j=0;j<nvbd;++j) {
				if (vbdry(j)->pnt == p0) {
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

#ifdef METIS4
extern "C" void METIS_PartMeshDual(int *ne, int *nn, int *elmnts, int *etype, int *numflag, int *nparts, int *edgecut, int *epart, int *npart);

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
	
	METIS_PartMeshDual(&ntri, &npnt, &tvrtx(0)(0), &etype, &numflag, &nparts, &edgecut,&(gbl->intwk(0)),&(gbl->i2wk(0)));
	
	for(i=0;i<ntri;++i)
		tri(i).info = gbl->intwk(i);
	
	gbl->intwk = -1;
	
	/* FIX METIS SO THAT WE DON'T HAVE ANY STRANDED TRIANGLES */
	for(int i=0;i<ntri;++i) {
		for(int j=0;j<3;++j) {
			int tind = tri(i).tri(j);
			if (tind > -1) 
				if (tri(i).info == tri(tind).info) goto next;
		}
		
		std::cout << "#reassigning problem triangle " << i << '\n';
		for(int j=0;j<3;++j) {
			int tind = tri(i).tri(j);
			if (tind > -1) {
				tri(i).info = tri(tind).info;
				break;
			}
		}	
		
	next:;
	}
	
	
	return;
}
#endif

#ifdef METIS
#include <metis.h>

void tri_mesh::setpartition(int nparts) {
	idx_t edgecut;
	idx_t ncommon = 2;
	idx_t ntri_idx = ntri;
 	idx_t npnt_idx = npnt;
	idx_t nparts_idx = nparts;	
	Array<idx_t,1> iwk1(maxpst), iwk2(maxpst), eptr(maxpst), eind(3*maxpst);

	for(int i=0;i<ntri;++i) {
		eptr(i) = 3*i;
		eind(eptr(i)) = tri(i).pnt(0);
		eind(eptr(i)+1) = tri(i).pnt(1);	
		eind(eptr(i)+2) = tri(i).pnt(2);	
	}
	eptr(ntri) = 3*ntri;
		
	int err = METIS_PartMeshDual(&ntri_idx, &npnt_idx, eptr.data(), eind.data(), NULL, NULL, &ncommon, &nparts_idx, NULL, NULL, &edgecut,iwk1.data(),iwk2.data());
	if (err != METIS_OK) {
		*gbl->log << "METIS partitioning error" << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
	for(int i=0;i<ntri;++i)
		tri(i).info = iwk1(i);
	
	/* FIX METIS SO THAT WE DON'T HAVE ANY STRANDED TRIANGLES */
	for(int i=0;i<ntri;++i) {
		for(int j=0;j<3;++j) {
			int tind = tri(i).tri(j);
			if (tind > -1) 
				if (tri(i).info == tri(tind).info) goto next;
		}
		
		std::cout << "#reassigning problem triangle " << i << '\n';
		for(int j=0;j<3;++j) {
			int tind = tri(i).tri(j);
			if (tind > -1) {
				tri(i).info = tri(tind).info;
				break;
			}
		}	
		
		next:;
	}
	return;
}
#endif

/* WHEN THIS ROUTINE EXITS */
/* THE OLD MESH STORES  */
/* pnt(pind).info = new pnt index or -1 */
/* THE NEW MESH STORES */
/* pnt(sind).info = old pnt index */
void tri_mesh::partition(class tri_mesh& xin, int npart) {
	int i,j,n,tind,sind,p0,indx;
	Array<int,2> bcntr(xin.nebd +10,2);
	int bnum,match;
	
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
		
	ntri = 0;
	for(i=0;i<xin.ntri;++i) {
		if (xin.tri(i).info == npart) {
			++ntri;
			for(n=0;n<3;++n) {
				xin.pnt(xin.tri(i).pnt(n)).info = npart;
			}
		}
	}
	
	if (!initialized) {
		maxpst = static_cast<int>(static_cast<FLT>(ntri*xin.maxpst)/xin.ntri);
		allocate(maxpst);
	}
	else if (3*ntri > maxpst) {
		*gbl->log << "mesh is too small" << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
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

	/* FILL IN TRI ARRAY */
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
		/* tri number in old mesh */
		tind = tri(seg(i).tri(0)).info;
		
		/* find first point of side */
		p0 = seg(i).pnt(0);
		for(n=0;n<3;++n)
			if (xin.pnt(xin.tri(tind).pnt(n)).info == p0) break;
		if (n==3) *gbl->log << "error in partitioning\n";
		
		// Rotate 2 to get side number
		n = (n+2)%3;
		
		/* store seg number in old mesh */
		seg(i).info = xin.tri(tind).seg(n);
			
		if (seg(i).tri(1) < 0) {
			indx = xin.tri(tind).tri(n);
			if (indx < 0) {
				/* BOUNDARY SIDE */
				bnum = getbdrynum(indx);

				for (j = 0; j <nebd;++j) {
					if (bnum == bcntr(j,0)) {
						++bcntr(j,1);
						seg(i).tri(1) = trinumatbdry(j,0);
						goto next1;
					}
				}
				/* NEW SIDE */
				seg(i).tri(1) = trinumatbdry(nebd, 0);
				bcntr(nebd,0) = bnum;
				bcntr(nebd++,1) = 1;
			}
			else {
				/* PARTITION SIDE */
				match = xin.tri(indx).info;
				if (match < npart) bnum = (match<<16) + (npart<<24);
				else bnum = (match<<24) + (npart<<16);
				for (j = 0; j <nebd;++j) {
					if (bcntr(j,0) == -bnum) {
						++bcntr(j,1);
						seg(i).tri(1) = trinumatbdry(j,match);
						goto next1;
					}
				}
				/* NEW SIDE */
				seg(i).tri(1) = trinumatbdry(nebd,match);
				bcntr(nebd,0) = -bnum;
				bcntr(nebd++,1) = 1;
			}
		}
		next1: continue;
	}

	/* Make new boundaries and print boundary information */
	ebdry.resize(nebd);
	for(i=0;i<nebd;++i) {
		if (bcntr(i,0) < 0) {
			ebdry(i) = new epartition(abs(bcntr(i,0)),*this);
		}
		else {
			ebdry(i) = xin.ebdry(bcntr(i,0))->create(*this);
		}
		ebdry(i)->alloc(static_cast<int>(bcntr(i,1)*2));
		ebdry(i)->nseg = 0;
	}


	for(i=0;i<nseg;++i) {
		if (seg(i).tri(1) < 0) {
			bnum = getbdrynum(seg(i).tri(1));
			ebdry(bnum)->seg(ebdry(bnum)->nseg++) = i;
		}
	}

	for(i=0;i<nebd;++i) {
		/* CREATES NEW BOUNDARY FOR DISCONNECTED SEGMENTS OF SAME TYPE */
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

	/* CREATE COMMUNICATION ENDPOINT BOUNDARIES */
	for(i=0;i<nebd;++i) {
		if (ebdry(i)->mytype == "partition") {			
			/* Now that all independent ebdry sides are determined give new numbers to partitions */
			sind = ebdry(i)->seg(0);
			match = getbdryseg(seg(sind).tri(1));
			if (npart < match) {
				int newid = seg(ebdry(i)->seg(0)).info +maxenum;
				ostringstream nstr;
				nstr << "b" << npart << "_s" << newid << std::flush;
				ebdry(i)->idprefix = nstr.str();
				ebdry(i)->idnum = newid;
				nstr.clear();
				std::cout << ebdry(i)->idprefix << "_type: partition\n";
			}
			else {
				int newid = seg(ebdry(i)->seg(ebdry(i)->nseg-1)).info +maxenum;
				ostringstream nstr;
				nstr << "b" << npart << "_s" << newid << std::flush;
				ebdry(i)->idprefix = nstr.str();
				ebdry(i)->idnum = newid;
				nstr.clear();	
				std::cout << ebdry(i)->idprefix << "_type: partition\n";
			}
			
			p0 = seg(sind).pnt(0);
			for(j=0;j<nvbd;++j)
				if (vbdry(j)->pnt == p0) goto nextv0;

			/* ENDPOINT IS NEW NEED TO DEFINE BOUNDARY */
			tind = tri(seg(sind).tri(0)).info;
			for(n=0;n<3;++n)
				if (xin.pnt(xin.tri(tind).pnt(n)).info == p0) break;
			vbdry(nvbd) = new vcomm(xin.tri(tind).pnt(n) +maxvnum,*this);
			vbdry(nvbd)->alloc(4);
			vbdry(nvbd)->pnt = p0;
			std::cout << vbdry(nvbd)->idprefix << "_type: comm\n";
			std::cout << vbdry(nvbd)->idprefix << "_group: 0 1 2\n";  // Put this boundary in the partition group
			++nvbd;

			nextv0:
			sind = ebdry(i)->seg(ebdry(i)->nseg-1);
			p0 = seg(sind).pnt(1);
			for(j=0;j<nvbd;++j)
				if (vbdry(j)->pnt == p0) goto nextv1;

			/* NEW ENDPOINT */
			tind = tri(seg(sind).tri(0)).info;
			for(n=0;n<3;++n)
				if (xin.pnt(xin.tri(tind).pnt(n)).info == p0) break;
			vbdry(nvbd) = new vcomm(xin.tri(tind).pnt(n)+maxvnum,*this);
			vbdry(nvbd)->alloc(4);
			vbdry(nvbd)->pnt = p0;
			std::cout << vbdry(nvbd)->idprefix << "_type: comm\n";
			std::cout << vbdry(nvbd)->idprefix << "_group: 0 1 2\n";  // Put this boundary in the partition group
			++nvbd;

			nextv1:
			continue;
		}
	}

	bdrylabel();  // CHANGES STRI / TTRI ON BOUNDARIES TO POINT TO GROUP/ELEMENT

	createtritri();
	createpnttri();
	cnt_nbor();
	TinyVector<FLT,ND> xmin, xmax;
	for(n=0;n<ND;++n) {
		xmin[n] = xin.qtree.xmin(n);
		xmax[n] = xin.qtree.xmax(n);
	}
	treeinit(xmin,xmax);

	initialized = 1;

	return;
}






