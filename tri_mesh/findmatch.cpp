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

#ifdef MPISRC
#include <parmetis.h>
#endif

void tri_mesh::setpartition(int nparts, Array<int,1> part_list) {
		
	/* Determine how many triangles are in each block so we can construct graph */
	const int nblock = sim::blks.nblock;
	Array<int,1> ntri_array_snd(nblock),ntri_array(nblock);
	
	ntri_array_snd = 0;
	ntri_array_snd(gbl->idnum) = ntri;
	sim::blks.allreduce(ntri_array_snd.data(),ntri_array.data(),nblock,blocks::int_msg,blocks::sum);
	
	/* Create pmetis graph */
	Array<idx_t,1> vrtxdist(nblock+1);
	vrtxdist(0) = 0;
	for(int i = 0; i < nblock; ++i)
		vrtxdist(i+1) = vrtxdist(i) +ntri_array(i);
	const int voffset = vrtxdist(gbl->idnum);
	
	/* Send triangle numbers along communication boundaries */
	for (int i=0;i<nebd;++i) {
		if (ebdry(i)->is_comm()) {
			int count = 0;
			for(int j=ebdry(i)->nseg-1;j>=0;--j) {
				int sind = ebdry(i)->seg(j);
				ebdry(i)->isndbuf(count++) = seg(sind).tri(0)+voffset;
			}
			ebdry(i)->sndsize() = count;
			ebdry(i)->sndtype() = boundary::int_msg;
		}
	}
	
	for(int i=0;i<nebd;++i)
		ebdry(i)->comm_prepare(boundary::all,0,boundary::symmetric);
	
	for(int i=0;i<nebd;++i)
		ebdry(i)->comm_exchange(boundary::all,0,boundary::symmetric);
	
	for(int i=0;i<nebd;++i)
		ebdry(i)->comm_wait(boundary::all,0,boundary::symmetric);
	
	/* Input parameters for ParMetis */
	idx_t ncon = 1;
	idx_t wgtflag = 3;
	idx_t numflag = 0;
	idx_t edgecut;
	real_t tpwgts[nparts];
	real_t ubvec[ncon];
	idx_t nparts_idx = nparts;
	Array<idx_t,1> part(ntri), adjwgt(6*ntri), xadj(ntri+1), adjncy(6*ntri), vwgt(ntri);
	idx_t options[3];
	options[0] = 0;
	for(int i=0;i<nparts;++i)
		tpwgts[i] = 1.0/nparts;
	ubvec[0] = 1.05;
	
	int count = 0;
	for(int i=0;i<ntri;++i) {
		xadj(i) = count;
		vwgt(i) = 1;
		for(int j=0;j<3;++j) {
			if (tri(i).tri(j) > -1) {
				adjncy(count) = tri(i).tri(j) +voffset;
				adjwgt(count) = 1;
				count++;
			}
			else {
				int b = getbdrynum(tri(i).tri(j));
				int seg = getbdryseg(tri(i).tri(j));
				if (ebdry(b)->is_comm()) {
					adjncy(count) = ebdry(b)->ircvbuf(0,seg);
					adjwgt(count) = 1000;  // Fails if partition goes parallel to physical comm boundary! */
					++count;
				}
			}
		}
	}
	xadj(ntri) = count;
	
#ifdef MPISRC
	MPI_Comm comm;
	MPI_Comm_dup(MPI_COMM_WORLD,&comm);
	
	int err = ParMETIS_V3_PartKway(vrtxdist.data(),xadj.data(),adjncy.data(),vwgt.data(), adjwgt.data(), &wgtflag, &numflag, &ncon, &nparts_idx, tpwgts, ubvec, options, &edgecut, part.data(),&comm);
	
	if (err != METIS_OK) {
		*gbl->log << "ParMETIS partitioning error" << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
#endif
	
	/* Determine how many partitions on each block */
	Array<int,1> part_list_snd(nparts*nblock+2);
	part_list_snd = 0;
	for(int i=0;i<ntri;++i) {
		part_list_snd(part(i) +gbl->idnum*nparts) = 1;
	}
	
	sim::blks.allreduce(part_list_snd.data(),part_list.data(),nblock*nparts,blocks::int_msg,blocks::sum);
		
	Array<int,1> part_count(nblock+1);
	part_count(0) = 0;
	for(int i=0;i<nblock;++i) {
		part_count(i+1) = part_count(i);
		for(int j=0;j<nparts;++j) {
			part_count(i+1) += part_list(i*nparts+j);
		}
	}
	
	int myparts = part_count(gbl->idnum);
	
	Array<int,1> new_part_numbers(nparts);
	for(int i=0;i<nparts;++i) {
		if (part_list(gbl->idnum*nparts+i))
			new_part_numbers(i) = myparts++;
		else
			new_part_numbers(i) = -1;
	}
	
	for(int i=0;i<ntri;++i) {
		tri(i).info = new_part_numbers(part(i));
	}
	
	/* Calculate global max of edge and vertex numbers */
	/* TO CREATE UNIQUE EDGE NUMBERS */
	int maxenum_snd = 0;
	for(int i=0;i<nebd;++i) {
		maxenum_snd = MAX(maxenum_snd,ebdry(i)->idnum);
	}
	int maxenum = maxenum_snd;
#ifdef MPISRC
	sim::blks.allreduce(&maxenum_snd, &maxenum, 1, blocks::int_msg, blocks::max);
#endif
	
	/* TO CREATE UNIQUE VERTEX NUMBERS */
	int maxvnum_snd = 0;
	for(int i=0;i<nvbd;++i) {
		maxvnum_snd = MAX(maxvnum_snd,vbdry(i)->idnum);
	}
	int maxvnum = maxvnum_snd;
#ifdef MPISRC
	sim::blks.allreduce(&maxvnum_snd, &maxvnum, 1, blocks::int_msg, blocks::max);
#endif
	
	part_list(nblock*nparts) = maxenum;
	part_list(nblock*nparts+1) = maxvnum;
	
	return;
}

void tri_mesh::subpartition(int& nparts) {
	idx_t ncon = 1;
	idx_t objval;
	idx_t ntri_idx = ntri;
	idx_t npnt_idx = npnt;
	idx_t nparts_idx = nparts;
	Array<idx_t,1> iwk1(ntri), adjwgt(6*ntri), xadj(ntri+1), adjncy(6*ntri), vwgt(ntri);
	idx_t options[METIS_NOPTIONS];
	
	METIS_SetDefaultOptions(options);
	options[METIS_OPTION_NITER] = 1000;
	options[METIS_OPTION_CONTIG] = 1;
	options[METIS_OPTION_CTYPE] = METIS_CTYPE_SHEM;
	
	int ndomains = 0;
	int count = 0;
	for(int i=0;i<ntri;++i) {
		xadj(i) = count;
		int mark = tri(i).info;
		vwgt(i) = 1+4*mark;
		ndomains = MAX(mark,ndomains);
		for(int j=0;j<3;++j) {
			if (tri(i).tri(j) > -1) {
				adjncy(count) = tri(i).tri(j);
				if (mark == tri(tri(i).tri(j)).info) {
					adjwgt(count) = 1;
				}
				else {
					adjwgt(count) = 1;
				}
				++count;
			}
		}
	}
	
	xadj(ntri) = count;
	++ndomains;
	//int err = METIS_PartGraphRecursive(&ntri_idx,&ncon, xadj.data(), adjncy.data(), vwgt.data(), NULL, adjwgt.data(), &nparts_idx, NULL, NULL, options, &objval, iwk1.data());
	
	int err = METIS_PartGraphKway(&ntri_idx,&ncon, xadj.data(), adjncy.data(), vwgt.data(), NULL, adjwgt.data(), &nparts_idx, NULL, NULL, options, &objval, iwk1.data());
	if (err != METIS_OK) {
		*gbl->log << "METIS partitioning error" << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
	
	/* Do an intersection of partitions */
	std::map<int,int> part_numbers;
	int nsubparts = 0;
	for(int i=0;i<ntri;++i) {
		int p1 = iwk1(i);
		int p2 = tri(i).info;
		if (part_numbers.find(p1*ndomains+p2) == part_numbers.end()) {
			part_numbers[p1*ndomains+p2] = nsubparts;
			tri(i).info = nsubparts;
			++nsubparts;
			
		}
		else {
			tri(i).info = part_numbers[p1*ndomains+p2];
		}
	}
	
	/* output number of parititions for each physical domain */
	for(int i=0;i<ndomains;++i) {
		int nsubdomain = 0;
		for (int j=0;j<nparts;++j) {
			if (part_numbers.find(j*ndomains+i) != part_numbers.end())
				++nsubdomain;
		}
		std::cout << 'b' << i << "_npartitions: " << nsubdomain << std::endl;
	}
	
	nparts = nsubparts;
	
	
	return;
}
#else
void tri_mesh::setpartition(int npart) {
	*gbl->log << "need metis to partition" << std::endl;
	sim::abort(__LINE__,__FILE__,gbl->log);
}
void tri_mesh::setpartition(int nparts,blitz::Array<int,1> part_list) {
	*gbl->log << "need metis to partition" << std::endl;
	sim::abort(__LINE__,__FILE__,gbl->log);
}
void tri_mesh::subpartition(int& nparts) {
	*gbl->log << "need metis to partition" << std::endl;
	sim::abort(__LINE__,__FILE__,gbl->log);
}
#endif

/* WHEN THIS ROUTINE EXITS */
/* THE OLD MESH STORES  */
/* pnt(pind).info = new pnt index or -1 */
/* THE NEW MESH STORES */
/* pnt(sind).info = old pnt index */
void tri_mesh::partition(multigrid_interface& xmesh, int npart, int maxenum, int maxvnum) {
	tri_mesh& xin(dynamic_cast<tri_mesh&>(xmesh));
	
	for(int i=0;i<xin.npnt;++i)
		xin.pnt(i).info = -1;
		
	ntri = 0;
	for(int i=0;i<xin.ntri;++i) {
		if (xin.tri(i).info == npart) {
			++ntri;
			for(int n=0;n<3;++n) {
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
	else {
		for(int i=0;i<nebd;++i)
			delete(ebdry(i));
		
		for(int i=0;i<nvbd;++i)
			delete(vbdry(i));
	}

	/* FILL IN PNT ARRAY */
	npnt = 0;
	for(int i=0;i<xin.npnt;++i) {
		if (xin.pnt(i).info == npart) {
			for(int n=0;n<ND;++n)
				pnts(npnt)(n) = xin.pnts(i)(n);
			xin.pnt(i).info = npnt;
			pnt(npnt).info = i;
			lngth(npnt) = xin.lngth(i);
			++npnt;
		}
	}

	/* FILL IN TRI ARRAY */
	ntri = 0;
	for(int i=0;i<xin.ntri;++i) {
		if (xin.tri(i).info == npart) {
			for(int n=0;n<3;++n)
				tri(ntri).pnt(n) = xin.pnt(xin.tri(i).pnt(n)).info;
			tri(ntri).info = i;
			++ntri;
		}
	}
	
	createseg();

	nebd = 0;
	std::vector<TinyVector<int,2> > bcntr;
	//Array<int,2> bcntr(xin.nebd +10,2);
	for(int i=0;i<nseg;++i) {
		/* tri number in old mesh */
		int tind = tri(seg(i).tri(0)).info;
		
		/* find first point of side */
		int p0 = seg(i).pnt(0);
		int n;
		for(n=0;n<3;++n)
			if (xin.pnt(xin.tri(tind).pnt(n)).info == p0) break;
		if (n==3) *gbl->log << "error in partitioning\n";
		
		// Rotate 2 to get side number
		n = (n+2)%3;
		
		/* store seg number in old mesh */
		seg(i).info = xin.tri(tind).seg(n);
		
		int bnum, match;
		if (seg(i).tri(1) < 0) {
			int indx = xin.tri(tind).tri(n);
			if (indx < 0) {
				/* PHYSICAL BOUNDARY */
				bnum = getbdrynum(indx);

				for (int j = 0; j <nebd;++j) {
					if (bnum == bcntr[j][0]) {
						++bcntr[j][1];
						seg(i).tri(1) = trinumatbdry(j,0);
						goto next1;
					}
				}
				/* NEW SIDE */
				seg(i).tri(1) = trinumatbdry(nebd, 0);
				bcntr.push_back(TinyVector<int,2>(bnum,1));
				++nebd;
			}
			else {
				/* PARTITION SIDE */
				match = xin.tri(indx).info;
				for (int j = 0; j <nebd;++j) {
					if (bcntr[j][0] == -(match+1)) {
						++bcntr[j][1];
						seg(i).tri(1) = trinumatbdry(j,match);
						goto next1;
					}
				}
				/* NEW SIDE */
				seg(i).tri(1) = trinumatbdry(nebd,match);
				bcntr.push_back(TinyVector<int,2>(-(match+1),1));
				++nebd;
			}
		}
		next1: continue;
	}
	
	/* MOVE VERTEX BOUNDARY INFO */
	nvbd = 0;
	for(int i=0;i<xin.nvbd;++i)
		if (xin.pnt(xin.vbdry(i)->pnt).info > -1)
			++nvbd;
	vbdry.resize(nvbd+2*nebd);
	
	nvbd = 0;
	for(int i=0;i<xin.nvbd;++i) {
		if (xin.pnt(xin.vbdry(i)->pnt).info > -1) {
			vbdry(nvbd) = xin.vbdry(i)->create(*this);
			vbdry(nvbd)->alloc(4);
			vbdry(nvbd)->pnt = xin.pnt(xin.vbdry(i)->pnt).info;
			ostringstream nstr;
			nstr << "b" << npart << "_v" << xin.vbdry(i)->idnum << std::flush;
			vbdry(nvbd)->idprefix = nstr.str();
			*gbl->log << nstr.str() << "_matches: " << xin.vbdry(i)->idprefix << std::endl;
			nstr.clear();
			++nvbd;
		}
	}

	/* Make new boundaries and print boundary information */
	ebdry.resize(nebd);
	for(int i=0;i<nebd;++i) {
		if (bcntr[i][0] < 0) {
			ebdry(i) = new epartition(abs(bcntr[i][0])-1,*this);
		}
		else {
			ebdry(i) = xin.ebdry(bcntr[i][0])->create(*this);
			if (!ebdry(i)->is_comm()) {
				/* Non communication boundaries do not need unique numbers */
				ostringstream nstr;
				nstr << "b" << npart << "_s" << xin.ebdry(bcntr[i][0])->idnum << std::flush;
				*gbl->log << nstr.str() << "_matches: " << ebdry(i)->idprefix << std::endl;
				ebdry(i)->idprefix = nstr.str();
				nstr.clear();
			}
		}
		ebdry(i)->alloc(static_cast<int>(bcntr[i][1]*2));
		ebdry(i)->nseg = 0;
	}


	for(int i=0;i<nseg;++i) {
		if (seg(i).tri(1) < 0) {
			int bnum = getbdrynum(seg(i).tri(1));
			ebdry(bnum)->seg(ebdry(bnum)->nseg++) = i;
		}
	}

	
	for(int i=0;i<nebd;++i) {
		/* CREATES NEW BOUNDARY FOR DISCONNECTED SEGMENTS OF SAME TYPE */
		ebdry(i)->reorder();
	}
	
#ifndef MPISRC
	/* TO CREATE UNIQUE EDGE NUMBERS */
	if (maxenum + maxvnum == 0) {
		for(int i=0;i<xin.nebd;++i) {
			maxenum = MAX(maxenum,xin.ebdry(i)->idnum);
		}
		
		/* TO CREATE UNIQUE VERTEX NUMBERS */
		for(int i=0;i<xin.nvbd;++i) {
			maxvnum = MAX(maxvnum,xin.vbdry(i)->idnum);
		}
	}
#else
	/* Load current definitions of new boundary conditions */
	std::map<std::pair<int,int>,int> physical_sides;
	std::map<std::pair<int,int>,int> partition_sides;
	std::map<std::pair<int,int>,int> physical_vertices;
	std::map<std::pair<int,int>,int> partition_vertices;
	

	sim::blks.begin_use_shared_memory();
	std::cout << "maxvnum " << maxvnum << std::endl;

	/* unpack current boundary information into usable form */
	int count = 0;
	int nentry = sim::blks.int_shared_memory(count++);
	for(int i=0;i<nentry;++i) {
		int type = sim::blks.int_shared_memory(count++);
		int identity1 = sim::blks.int_shared_memory(count++);
		int identity2 = sim::blks.int_shared_memory(count++);
		int idnum = sim::blks.int_shared_memory(count++);
		std::cout << 'b' << xin.gbl->idnum << " p" << npart << " type:" << type << ' ' << identity1 << ' ' << identity2 << ' ' << idnum << std::endl;
		switch(type) {
			case(0): {
				physical_sides[std::pair<int,int>(identity1,identity2)] = idnum;
				maxenum = MAX(idnum,maxenum);
				break;
			}
			case(1): {
				partition_sides[std::pair<int,int>(identity1,identity2)] = idnum;
				maxenum = MAX(idnum,maxenum);
				break;
			}
			case(2): {
				physical_vertices[std::pair<int,int>(identity1,identity2)] = idnum;
				maxvnum = max(idnum,maxvnum);
				break;
			}
			case(3): {
				partition_vertices[std::pair<int,int>(identity1,identity2)] = idnum;
				maxvnum = max(idnum,maxvnum);
				break;
			}
		}
	}
#endif
	
	std::cout << "maxvnum " << maxvnum << std::endl;
	
#ifdef MPISRC
	/* CREATE COMMUNICATION ENDPOINT BOUNDARIES FOR PHYSICAL COMM BOUNDARIES */
	for(int i=0;i<nebd;++i) {
		if (ebdry(i)->is_comm() && ebdry(i)->mytype != "partition") {
			
			/* Find idnum */
			int identity1 = ebdry(i)->idnum;
			int is_frst = xin.ebdry(getbdrynum(xin.seg(seg(ebdry(i)->seg(0)).info).tri(1)))->is_frst();
			int identity2;
			if (is_frst) {
				identity2 = getbdryseg(xin.seg(seg(ebdry(i)->seg(0)).info).tri(1));
			}
			else {
				int t1 = xin.seg(seg(ebdry(i)->seg(ebdry(i)->nseg-1)).info).tri(1);
				identity2 = xin.ebdry(getbdrynum(t1))->nseg-1 -getbdryseg(t1);
			}
			std::pair<int,int> bid(identity1,identity2);
			int newid;
			if (physical_sides.find(bid) != physical_sides.end())
				newid = physical_sides[bid];
			else {
				newid = ++maxenum;
				sim::blks.int_shared_memory(0)++;
				sim::blks.int_shared_memory(count++) = 0;
				sim::blks.int_shared_memory(count++) = identity1;
				sim::blks.int_shared_memory(count++) = identity2;
				sim::blks.int_shared_memory(count++) = newid;
				std::cout << 'b' << xin.gbl->idnum << " p" << npart << " type:" << 0 << ' ' << identity1 << ' ' << identity2 << ' ' << newid << std::endl;
			}
			ostringstream nstr;
			nstr << "b" << npart << "_s" << newid << std::flush;
			*gbl->log << nstr.str() << "_matches: " << ebdry(i)->idprefix << std::endl;
			ebdry(i)->idprefix = nstr.str();
			ebdry(i)->idnum = newid;
			nstr.str("");
			nstr.clear();
			

			/* Create vertex endpoints */
			int sind = ebdry(i)->seg(0);
			int p0 = seg(sind).pnt(0);
			int n, tind;
			for(int j=0;j<nvbd;++j)
				if (vbdry(j)->pnt == p0) goto nextv0;
			
			/* ENDPOINT IS NEW NEED TO DEFINE BOUNDARY */
			tind = tri(seg(sind).tri(0)).info;
			for(n=0;n<3;++n)
				if (xin.pnt(xin.tri(tind).pnt(n)).info == p0) break;
			
			if (is_frst) {
				identity2 = getbdryseg(xin.seg(seg(sind).info).tri(1));
			}
			else {
				int t1 = xin.seg(seg(sind).info).tri(1);
				identity2 = xin.ebdry(getbdrynum(t1))->nseg -getbdryseg(t1);
			}
			
			bid = std::pair<int,int>(identity1,identity2);
			if (physical_vertices.find(bid) != physical_vertices.end())
				newid = physical_vertices[bid];
			else {
				newid = ++maxvnum;
				sim::blks.int_shared_memory(0)++;
				sim::blks.int_shared_memory(count++) = 2;
				sim::blks.int_shared_memory(count++) = identity1;
				sim::blks.int_shared_memory(count++) = identity2;
				sim::blks.int_shared_memory(count++) = newid;
				std::cout << 'b' << xin.gbl->idnum << " p" << npart << " type:" << 2 << ' ' << identity1 << ' ' << identity2 << ' ' << newid << std::endl;
			}
			vbdry(nvbd) = new vcomm(newid,*this);
			vbdry(nvbd)->alloc(4);
			vbdry(nvbd)->pnt = p0;
			nstr << "b" << npart << "_v" << newid << std::flush;
			vbdry(nvbd)->idprefix = nstr.str();
			*gbl->log << vbdry(nvbd)->idprefix << "_type: comm\n";
			*gbl->log << vbdry(nvbd)->idprefix << "_group: 0 1 2\n";  // Put this boundary in the partition group
			nstr.str("");
			nstr.clear();
			++nvbd;
			
		nextv0:
			sind = ebdry(i)->seg(ebdry(i)->nseg-1);
			p0 = seg(sind).pnt(1);
			for(int j=0;j<nvbd;++j)
				if (vbdry(j)->pnt == p0) goto nextv1;
			
			/* NEW ENDPOINT */
			tind = tri(seg(sind).tri(0)).info;
			for(n=0;n<3;++n)
				if (xin.pnt(xin.tri(tind).pnt(n)).info == p0) break;
			
			if (is_frst) {
				identity2 = getbdryseg(xin.seg(seg(sind).info).tri(1)) +1;
			}
			else {
				int t1 = xin.seg(seg(sind).info).tri(1);
				identity2 = xin.ebdry(getbdrynum(t1))->nseg-1 -getbdryseg(t1);
			}
			
			bid = std::pair<int,int>(identity1,identity2);
			if (physical_vertices.find(bid) != physical_vertices.end())
				newid = physical_vertices[bid];
			else {
				newid = ++maxvnum;
				sim::blks.int_shared_memory(0)++;
				sim::blks.int_shared_memory(count++) = 2;
				sim::blks.int_shared_memory(count++) = identity1;
				sim::blks.int_shared_memory(count++) = identity2;
				sim::blks.int_shared_memory(count++) = newid;
				std::cout << 'b' << xin.gbl->idnum << " p" << npart << " type:" << 2 << ' ' << identity1 << ' ' << identity2 << ' ' << newid << std::endl;
			}
			vbdry(nvbd) = new vcomm(newid,*this);
			vbdry(nvbd)->alloc(4);
			vbdry(nvbd)->pnt = p0;
			nstr << "b" << npart << "_v" << newid << std::flush;
			vbdry(nvbd)->idprefix = nstr.str();
			*gbl->log << vbdry(nvbd)->idprefix << "_type: comm\n";
			*gbl->log << vbdry(nvbd)->idprefix << "_group: 0 1 2\n";  // Put this boundary in the partition group
			nstr.str("");
			nstr.clear();
			++nvbd;
			
		nextv1:
			continue;
		}
	}
#endif
	
	/* RENUMBER PARTITION BOUNDARIES AND CREATE COMMUNICATION ENDPOINT BOUNDARIES */
	for(int i=0;i<nebd;++i) {
		if (ebdry(i)->mytype == "partition") {
			/* Now that all independent ebdry sides are determined give new numbers to partitions */
			int sind = ebdry(i)->seg(0);
			int match = getbdryseg(seg(sind).tri(1));  // for partitions match is stored in location for bdry seg above
			int identity1 = xin.gbl->idnum;
			int identity2;
			int newid;
			if (npart < match) {
				identity2 = seg(ebdry(i)->seg(0)).info;
				newid = seg(ebdry(i)->seg(0)).info +maxenum+1;

			}
			else {
				identity2 = seg(ebdry(i)->seg(ebdry(i)->nseg-1)).info;
				newid = seg(ebdry(i)->seg(ebdry(i)->nseg-1)).info +maxenum+1;
			}
#ifdef MPISRC
			std::pair<int,int> bid(identity1,identity2);
			if (partition_sides.find(bid) != partition_sides.end())
				newid = partition_sides[bid];
			else {
				newid = ++maxenum;
				sim::blks.int_shared_memory(0)++;
				sim::blks.int_shared_memory(count++) = 1;
				sim::blks.int_shared_memory(count++) = identity1;
				sim::blks.int_shared_memory(count++) = identity2;
				sim::blks.int_shared_memory(count++) = newid;
				std::cout << 'b' << xin.gbl->idnum << " p" << npart << " type:" << 1 << ' ' << identity1 << ' ' << identity2 << ' ' << newid << std::endl;
			}
#endif
			ostringstream nstr;
			nstr << "b" << npart << "_s" << newid << std::flush;
			ebdry(i)->idprefix = nstr.str();
			ebdry(i)->idnum = newid;
			nstr.str("");
			nstr.clear();
			*gbl->log << ebdry(i)->idprefix << "_type: partition\n";
			

			
			/* Create vertex endpoints */
			int p0 = seg(sind).pnt(0);
			int n, tind;
			for(int j=0;j<nvbd;++j)
				if (vbdry(j)->pnt == p0) goto nextv2;
			
			/* ENDPOINT IS NEW NEED TO DEFINE BOUNDARY */
			tind = tri(seg(sind).tri(0)).info;
			for(n=0;n<3;++n)
				if (xin.pnt(xin.tri(tind).pnt(n)).info == p0) break;
			
			identity2 = xin.tri(tind).pnt(n);
			newid = maxvnum +identity2 +1;
#ifdef MPISRC
			bid = std::pair<int,int>(identity1,identity2);
			if (partition_vertices.find(bid) != partition_vertices.end())
				newid = partition_vertices[bid];
			else {
				newid = ++maxvnum;
				sim::blks.int_shared_memory(0)++;
				sim::blks.int_shared_memory(count++) = 3;
				sim::blks.int_shared_memory(count++) = identity1;
				sim::blks.int_shared_memory(count++) = identity2;
				sim::blks.int_shared_memory(count++) = newid;
				std::cout << 'b' << xin.gbl->idnum << " p" << npart << " type:" << 3 << ' ' << identity1 << ' ' << identity2 << ' ' << newid << std::endl;
			}
#endif
			vbdry(nvbd) = new vcomm(newid,*this);
			vbdry(nvbd)->alloc(4);
			vbdry(nvbd)->pnt = p0;
			nstr << "b" << npart << "_v" << newid << std::flush;
			vbdry(nvbd)->idprefix = nstr.str();
			*gbl->log << vbdry(nvbd)->idprefix << "_type: comm\n";
			*gbl->log << vbdry(nvbd)->idprefix << "_group: 0 1 2\n";  // Put this boundary in the partition group
			nstr.str("");
			nstr.clear();
			++nvbd;
			
		nextv2:
			sind = ebdry(i)->seg(ebdry(i)->nseg-1);
			p0 = seg(sind).pnt(1);
			for(int j=0;j<nvbd;++j)
				if (vbdry(j)->pnt == p0) goto nextv3;
			
			/* NEW ENDPOINT */
			tind = tri(seg(sind).tri(0)).info;
			for(n=0;n<3;++n)
				if (xin.pnt(xin.tri(tind).pnt(n)).info == p0) break;
			
			identity2 = xin.tri(tind).pnt(n);
			newid = maxvnum +identity2 +1;
#ifdef MPISRC
			bid = std::pair<int,int>(identity1,identity2);
			if (partition_vertices.find(bid) != partition_vertices.end())
				newid = partition_vertices[bid];
			else {
				newid = ++maxvnum;
				sim::blks.int_shared_memory(0)++;
				sim::blks.int_shared_memory(count++) = 3;
				sim::blks.int_shared_memory(count++) = identity1;
				sim::blks.int_shared_memory(count++) = identity2;
				sim::blks.int_shared_memory(count++) = newid;
				std::cout << 'b' << xin.gbl->idnum << " p" << npart << " type:" << 3 << ' ' << identity1 << ' ' << identity2 << ' ' << newid << std::endl;
			}
#endif
			vbdry(nvbd) = new vcomm(newid,*this);
			vbdry(nvbd)->alloc(4);
			vbdry(nvbd)->pnt = p0;
			nstr << "b" << npart << "_v" << newid << std::flush;
			vbdry(nvbd)->idprefix = nstr.str();
			*gbl->log << vbdry(nvbd)->idprefix << "_type: comm\n";
			*gbl->log << vbdry(nvbd)->idprefix << "_group: 0 1 2\n";  // Put this boundary in the partition group
			nstr.str("");
			nstr.clear();
			++nvbd;
			
		nextv3:
			continue;
		}
	}
	
#ifdef MPISRC
	sim::blks.end_use_shared_memory();
#endif

	bdrylabel();  // CHANGES STRI / TTRI ON BOUNDARIES TO POINT TO GROUP/ELEMENT
	
	/* FIND ENDPOINT MATCHES */
	for(int i=0;i<nvbd;++i) {
		/* Find two connecting boundary sides */
		for(int j=0;j<nebd;++j) {
			if (seg(ebdry(j)->seg(0)).pnt(0) == vbdry(i)->pnt) {
				vbdry(i)->ebdry(1) = j;
				ebdry(j)->vbdry(0) = i;
			}
			if (seg(ebdry(j)->seg(ebdry(j)->nseg-1)).pnt(1) == vbdry(i)->pnt) {
				vbdry(i)->ebdry(0) = j;
				ebdry(j)->vbdry(1) = i;
			}
		}
	}

	createtritri();
	createpnttri();
	cnt_nbor();
	TinyVector<FLT,ND> xmin, xmax;
	for(int n=0;n<ND;++n) {
		xmin[n] = xin.qtree.xmin(n);
		xmax[n] = xin.qtree.xmax(n);
	}
	treeinit(xmin,xmax);

	initialized = 1;

	return;
}






