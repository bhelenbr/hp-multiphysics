/*
 *  pod.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 7/26/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include <myblas.h>
#include <libbinio/binfile.h>

struct bd_str {
	int nmodes;
	int multiplicity;
};

template<class BASE> void pod_simulate<BASE>::init(input_map& input, void *gin) {
	std::string filename,keyword,linebuff;
	std::ostringstream nstr;
	std::istringstream instr;
	int i;

	/* Initialize base class */
	BASE::init(input,gin);

	input.getwdefault(BASE::gbl->idprefix + "_groups",pod_id,0);

	nstr.str("");
	nstr << "pod" << pod_id << "_nmodes";
	if (!input.get(nstr.str(),nmodes)) input.getwdefault("nmodes",nmodes,5); 
	nstr.clear();

	vsi ugstore;
	ugstore.v.reference(BASE::ugbd(0).v);
	ugstore.s.reference(BASE::ugbd(0).s);
	ugstore.i.reference(BASE::ugbd(0).i);

	modes.resize(nmodes);
	for(i=0;i<nmodes;++i) {
		nstr.str("");
		nstr << i << std::flush;
		filename = "mode" +nstr.str() +"_" +BASE::gbl->idprefix;
		nstr.clear();
		modes(i).v.resize(BASE::maxpst,BASE::NV);
		modes(i).s.resize(BASE::maxpst,BASE::sm0,BASE::NV);
		modes(i).i.resize(BASE::maxpst,BASE::im0,BASE::NV);
		BASE::ugbd(0).v.reference(modes(i).v);
		BASE::ugbd(0).s.reference(modes(i).s);
		BASE::ugbd(0).i.reference(modes(i).i);
		BASE::input(filename, BASE::binary);
	}
	BASE::ugbd(0).v.reference(ugstore.v);
	BASE::ugbd(0).s.reference(ugstore.s);
	BASE::ugbd(0).i.reference(ugstore.i);

	pod_ebdry.resize(BASE::nebd);
	/* Count how many boundary modes there are so we can size arrays before initializing boundaries */
	/* For each mesh block that is part of this pod block need to know
	/* # of pod boundaries, ids, and # of pod modes for each unique id */

	/* First need to know what block # I am and how many total blocks there are in this pod group */
	int localid = sim::blks.allreduce_local_id(pod_id,BASE::gbl->idnum);
	int npodblk = sim::blks.allreduce_nmember(pod_id);
	Array<int,1> binfo(npodblk),binfo_recv(npodblk);
	binfo = 0;

	for (int i=0;i<BASE::nebd;++i) {
		/* Not going to initialize until I can resize coeffs & rsdls arrays to accomodate boundary modes */
		pod_ebdry(i) = new pod_sim_edge_bdry<BASE>(*this,*BASE::ebdry(i));

		keyword = pod_ebdry(i)->base.idprefix +"_pod";
		input.getwdefault(keyword,pod_ebdry(i)->active,false);
		if (!pod_ebdry(i)->active) {
			pod_ebdry(i)->nmodes = 0;
			continue;
		}
		binfo(localid)++;

		keyword = pod_ebdry(i)->base.idprefix + "_pod_id";
		input.getwdefault(keyword,pod_ebdry(i)->pod_id,pod_ebdry(i)->base.idnum);

		nstr.str("");
		nstr << "bdry_pod" << pod_ebdry(i)->pod_id << "_nmodes";
		if (!input.get(nstr.str(),pod_ebdry(i)->nmodes)) input.getwdefault("bdry_nmodes",pod_ebdry(i)->nmodes,nmodes); 
		nstr.clear();
	}

	/* Send number of pod boundaries belonging to each block of pod group */
	sim::blks.allreduce(binfo.data(), binfo_recv.data(), npodblk,blocks::int_msg, blocks::sum, pod_id);


	/* Second thing is to pass id, nmodes for each active boundary */
	/* Count # of pod_boundaries before this block and total # of pod boundaries in pod group */
	int nbefore = 0;
	for (int i=0;i<localid;++i) {
		nbefore += binfo_recv(i);
	}

	int ntotal = nbefore;
	for (int i=localid;i<npodblk;++i) {
		ntotal += binfo_recv(i);
	}

	binfo.resize(ntotal*2);
	binfo_recv.resize(ntotal*2);
	nbefore *= 2;
	binfo = 0;
	for (int i=0;i<BASE::nebd;++i) {
		if (!pod_ebdry(i)->active) continue;

		binfo(nbefore++) = pod_ebdry(i)->pod_id;
		binfo(nbefore++) = pod_ebdry(i)->nmodes;
	}
	sim::blks.allreduce(binfo.data(), binfo_recv.data(), ntotal*2,blocks::int_msg, blocks::sum, pod_id);

	/* Now make a map from pod_id to number of modes & multiplicity */
	std::map<int,bd_str> pod_bdry_map;
	tmodes = nmodes;
	for (int i=0;i<2*ntotal;i+=2) {
		if (pod_bdry_map.find(binfo_recv(i)) != pod_bdry_map.end()) {
			++pod_bdry_map[binfo_recv(i)].multiplicity;
		}
		else {
			pod_bdry_map[binfo_recv(i)] = bd_str();
			pod_bdry_map[binfo_recv(i)].multiplicity = 1;
			pod_bdry_map[binfo_recv(i)].nmodes = binfo_recv(i+1);
			tmodes += binfo_recv(i+1);
		}
	}
	*BASE::gbl->log << "#There are " << tmodes << " total modes on pod block " << pod_id << std::endl;

	coeffs.resize(tmodes);
	rsdls.resize(tmodes);
	rsdls0.resize(tmodes);
	rsdls_recv.resize(tmodes);
	multiplicity.resize(tmodes);
	jacobian.resize(tmodes,tmodes);
	ipiv.resize(tmodes);

	/* Count total number of boundary modes */
	/* and make map be an accrual of previous modes */
	multiplicity = 1.0;
	int bindex = nmodes;
	for (std::map<int,bd_str>::iterator mi = pod_bdry_map.begin(); mi != pod_bdry_map.end(); ++mi) {
		int n = mi->second.nmodes;
		mi->second.nmodes = bindex;
		multiplicity(Range(bindex,bindex+n-1)) = mi->second.multiplicity;
		bindex += n;
	}


	int restartfile;
	input.getwdefault("restart",restartfile,1);
	nstr.str("");
	nstr << restartfile << std::flush;
	filename = "coeff" +nstr.str() +"_" +BASE::gbl->idprefix +".bin";
	binifstream bin;
	bin.open(filename.c_str());
	if (bin.error()) {
		*BASE::gbl->log << "couldn't open coefficient input file " << filename << std::endl;
		exit(1);
	}
	bin.setFlag(binio::BigEndian,bin.readInt(1));
	bin.setFlag(binio::FloatIEEE,bin.readInt(1));

	/* CONSTRUCT INITIAL SOLUTION DESCRIPTION */
	BASE::ug.v(Range(0,BASE::npnt-1)) = 0.;
	BASE::ug.s(Range(0,BASE::nseg-1)) = 0.;
	BASE::ug.i(Range(0,BASE::ntri-1)) = 0.;

	for (int l=0;l<nmodes;++l) {
		coeffs(l) = bin.readFloat(binio::Double); 
		BASE::ug.v(Range(0,BASE::npnt-1)) += coeffs(l)*modes(l).v(Range(0,BASE::npnt-1));
		BASE::ug.s(Range(0,BASE::nseg-1)) += coeffs(l)*modes(l).s(Range(0,BASE::nseg-1));
		BASE::ug.i(Range(0,BASE::ntri-1)) += coeffs(l)*modes(l).i(Range(0,BASE::ntri-1));
	}
	bin.close();

	/* Let boundary conditions load to and from coeff/rsdls vectors */
	/* Then initialize them */
	for (int i=0;i<BASE::nebd;++i) {
		if (!pod_ebdry(i)->active) continue;
		pod_ebdry(i)->bindex = pod_bdry_map[pod_ebdry(i)->pod_id].nmodes;
		pod_ebdry(i)->init(input);
	}	

	// *BASE::gbl->log << multiplicity << std::endl;


	return;
}

template<class BASE> void pod_simulate<BASE>::rsdl(int stage) {

	BASE::rsdl(stage);

	rsdls = 0.0;  // Need to do this, because not every block will touch every rsdl

	/* APPLY VERTEX DIRICHLET B.C.'S */
	for(int i=0;i<BASE::nebd;++i)
		BASE::hp_ebdry(i)->vdirichlet();

	for(int i=0;i<BASE::nvbd;++i)
		BASE::hp_vbdry(i)->vdirichlet2d();

	/* APPLY DIRCHLET B.C.S TO MODE */
	for(int i=0;i<BASE::nebd;++i)
		for(int sm=0;sm<basis::tri(BASE::log2p)->sm();++sm)
			BASE::hp_ebdry(i)->sdirichlet(sm);

	for (int k = 0; k < nmodes; ++k) {
		for(int i=0; i<BASE::npnt;++i)
			for(int n=0;n<BASE::NV;++n)
				rsdls(k) += modes(k).v(i,n)*BASE::gbl->res.v(i,n);

		for(int i=0; i<BASE::nseg;++i)
			for(int sm=0;sm<basis::tri(BASE::log2p)->sm();++sm)
				for(int n=0;n<BASE::NV;++n)
					rsdls(k) += modes(k).s(i,sm,n)*BASE::gbl->res.s(i,sm,n);

		for(int i=0; i<BASE::ntri;++i)
			for(int im=0;im<basis::tri(BASE::log2p)->im();++im)
				for(int n=0;n<BASE::NV;++n)
					rsdls(k) += modes(k).i(i,im,n)*BASE::gbl->res.i(i,im,n);

	}

	/* FORM RESIDUALS FOR SIDE MODES */
	for (int i=0;i<BASE::nebd;++i)
		pod_ebdry(i)->rsdl();

	/* COMMUNICATE POD BDRY RESDIUALS */
	for(int i=0;i<BASE::nebd;++i)
		pod_ebdry(i)->loadbuff(rsdls_recv);

	for(int i=0;i<BASE::nebd;++i)
		if (pod_ebdry(i)->active) BASE::ebdry(i)->comm_prepare(boundary::all,0,boundary::symmetric);

	for(int i=0;i<BASE::nebd;++i) 
		if (pod_ebdry(i)->active) BASE::ebdry(i)->comm_exchange(boundary::all,0,boundary::symmetric);

	for(int i=0;i<BASE::nebd;++i)
		if (pod_ebdry(i)->active) BASE::ebdry(i)->comm_wait(boundary::all,0,boundary::symmetric);

	for(int i=0;i<BASE::nebd;++i) 
		pod_ebdry(i)->finalrcv(rsdls_recv);

	sim::blks.allreduce(rsdls.data(),rsdls_recv.data(),tmodes,blocks::flt_msg,blocks::sum,pod_id);	

	return;
}


template<class BASE> void pod_simulate<BASE>::setup_preconditioner() {

	BASE::setup_preconditioner();
	rsdl(BASE::gbl->nstage);

	/* STORE BASELINE IN LAST COLUMN */
	jacobian(Range(0,tmodes-1),tmodes-1) = rsdls_recv;
	BASE::gbl->ug0.v(Range(0,BASE::npnt-1)) = BASE::ug.v(Range(0,BASE::npnt-1));
	BASE::gbl->ug0.s(Range(0,BASE::nseg-1)) = BASE::ug.s(Range(0,BASE::nseg-1));
	BASE::gbl->ug0.i(Range(0,BASE::ntri-1)) = BASE::ug.i(Range(0,BASE::ntri-1));

	for (int modeloop = 0; modeloop < nmodes; ++modeloop) {
		/* PERTURB EACH COEFFICIENT */
		FLT delta = 1.0e-2*coeffs(modeloop) +1.0e-12;
		BASE::gbl->res.v(Range(0,BASE::npnt-1)) = delta*modes(modeloop).v(Range(0,BASE::npnt-1));
		BASE::gbl->res.s(Range(0,BASE::nseg-1)) = delta*modes(modeloop).s(Range(0,BASE::nseg-1));
		BASE::gbl->res.i(Range(0,BASE::ntri-1)) = delta*modes(modeloop).i(Range(0,BASE::ntri-1));

		/* APPLY VERTEX DIRICHLET B.C.'S */
		for(int i=0;i<BASE::nebd;++i)
			BASE::hp_ebdry(i)->vdirichlet();

		for(int i=0;i<BASE::nvbd;++i)
			BASE::hp_vbdry(i)->vdirichlet2d();

		/* APPLY DIRCHLET B.C.S TO MODE */
		for(int i=0;i<BASE::nebd;++i)
			for(int sm=0;sm<basis::tri(BASE::log2p)->sm();++sm)
				BASE::hp_ebdry(i)->sdirichlet(sm);

		BASE::ug.v(Range(0,BASE::npnt-1)) = BASE::gbl->ug0.v(Range(0,BASE::npnt-1)) +BASE::gbl->res.v(Range(0,BASE::npnt-1));
		BASE::ug.s(Range(0,BASE::nseg-1)) = BASE::gbl->ug0.s(Range(0,BASE::nseg-1)) +BASE::gbl->res.s(Range(0,BASE::nseg-1));
		BASE::ug.i(Range(0,BASE::ntri-1)) = BASE::gbl->ug0.i(Range(0,BASE::ntri-1)) +BASE::gbl->res.i(Range(0,BASE::ntri-1));

		rsdl(BASE::gbl->nstage);

		/* STORE IN ROW */
		jacobian(Range(0,tmodes-1),modeloop) = (rsdls_recv -jacobian(Range(0,tmodes-1),tmodes-1))/delta;
	}

	/* CREATE JACOBIAN FOR BOUNDARY MODES */
	for (int modeloop = nmodes; modeloop < tmodes; ++modeloop) {

		/* ZERO COEFFICIENT */
		BASE::gbl->res.v(Range(0,BASE::npnt-1)) = 0.0;
		BASE::gbl->res.s(Range(0,BASE::nseg-1)) = 0.0; 
		BASE::gbl->res.i(Range(0,BASE::ntri-1)) = 0.0;

		for (int bind=0;bind<BASE::nebd;++bind) {		
			pod_ebdry(bind)->addto2Dsolution(BASE::gbl->res,modeloop,1.0e-4);
		}

		/* APPLY VERTEX DIRICHLET B.C.'S */
		for(int i=0;i<BASE::nebd;++i)
			BASE::hp_ebdry(i)->vdirichlet();

		for(int i=0;i<BASE::nvbd;++i)
			BASE::hp_vbdry(i)->vdirichlet2d();

		BASE::ug.v(Range(0,BASE::npnt-1)) = BASE::gbl->ug0.v(Range(0,BASE::npnt-1)) +BASE::gbl->res.v(Range(0,BASE::npnt-1));
		BASE::ug.s(Range(0,BASE::nseg-1)) = BASE::gbl->ug0.s(Range(0,BASE::nseg-1)) +BASE::gbl->res.s(Range(0,BASE::nseg-1));
		BASE::ug.i(Range(0,BASE::ntri-1)) = BASE::gbl->ug0.i(Range(0,BASE::ntri-1)) +BASE::gbl->res.i(Range(0,BASE::ntri-1));

		rsdl(BASE::gbl->nstage);

		/* STORE IN ROW */
		jacobian(Range(0,tmodes-1),modeloop) = (rsdls_recv -jacobian(Range(0,tmodes-1),tmodes-1))/1.0e-4;
	}

	/* RESTORE UG & COEFF VECTOR */
	BASE::ug.v(Range(0,BASE::npnt-1)) = BASE::gbl->ug0.v(Range(0,BASE::npnt-1));
	BASE::ug.s(Range(0,BASE::nseg-1)) = BASE::gbl->ug0.s(Range(0,BASE::nseg-1));
	BASE::ug.i(Range(0,BASE::ntri-1)) = BASE::gbl->ug0.i(Range(0,BASE::ntri-1));

	/* FACTORIZE PRECONDITIONER */
	int info;
	GETRF(tmodes,tmodes,jacobian.data(),tmodes,ipiv.data(),info);
	if (info != 0) {
		*BASE::gbl->log << "DGETRF FAILED FOR POD JACOBIAN " << info << std::endl;
		exit(1);
	}
	return;
}

template<class BASE> void pod_simulate<BASE>::update() {
	char trans[] = "T";
	int info;

	rsdl(BASE::gbl->nstage);

	GETRS(trans,tmodes,1,jacobian.data(),tmodes,ipiv.data(),rsdls_recv.data(),tmodes,info);
	if (info != 0) {
		*BASE::gbl->log << "DGETRS FAILED FOR POD UPDATE " << info << std::endl;
		exit(1);
	}
	/* Need to fix pod boundaries so coefficients are equal */
	/* store rsdls_recv to compare to after boundary comm */
	rsdls0 = rsdls_recv;	
	rsdls = rsdls_recv;

	/* COMMUNICATE POD BDRY CORRECTIONS */
	for(int i=0;i<BASE::nebd;++i)
		pod_ebdry(i)->loadbuff(rsdls);

	for(int i=0;i<BASE::nebd;++i)
		if (pod_ebdry(i)->active) BASE::ebdry(i)->comm_prepare(boundary::all,0,boundary::symmetric);

	for(int i=0;i<BASE::nebd;++i) 
		if (pod_ebdry(i)->active) BASE::ebdry(i)->comm_exchange(boundary::all,0,boundary::symmetric);

	for(int i=0;i<BASE::nebd;++i)
		if (pod_ebdry(i)->active) BASE::ebdry(i)->comm_wait(boundary::all,0,boundary::symmetric);

	for(int i=0;i<BASE::nebd;++i) 
		pod_ebdry(i)->finalrcv(rsdls);

	rsdls -= rsdls0;  // rsdls is now the correction to apply at the boundaries rsdls = rsdls -rsdls0

	/* Send correction to other blocks */
	sim::blks.allreduce(rsdls.data(),rsdls_recv.data(),tmodes,blocks::flt_msg,blocks::sum,pod_id);

	/* Problem is that for pod boundaries spanning more than 1 blocks, the correction gets added mulitple times */
	rsdls_recv /= multiplicity;
	rsdls_recv += rsdls0;

	coeffs -= rsdls_recv;

	BASE::gbl->res.v(Range(0,BASE::npnt-1)) = 0.0;
	BASE::gbl->res.s(Range(0,BASE::nseg-1)) = 0.0;
	BASE::gbl->res.i(Range(0,BASE::ntri-1)) = 0.0;

	for (int m=0;m<nmodes;++m) {
		BASE::gbl->res.v(Range(0,BASE::npnt-1)) += rsdls_recv(m)*modes(m).v(Range(0,BASE::npnt-1));
		BASE::gbl->res.s(Range(0,BASE::nseg-1)) += rsdls_recv(m)*modes(m).s(Range(0,BASE::nseg-1));
		BASE::gbl->res.i(Range(0,BASE::ntri-1)) += rsdls_recv(m)*modes(m).i(Range(0,BASE::ntri-1));
	}

	for (int m=nmodes;m<tmodes;++m) {
		for (int bind=0;bind<BASE::nebd;++bind) {		
			pod_ebdry(bind)->addto2Dsolution(BASE::gbl->res,m,rsdls_recv(m));
		}
	}

	/* APPLY VERTEX DIRICHLET B.C.'S */
	for(int i=0;i<BASE::nebd;++i)
		BASE::hp_ebdry(i)->vdirichlet();

	for(int i=0;i<BASE::nvbd;++i)
		BASE::hp_vbdry(i)->vdirichlet2d();

	/* APPLY DIRCHLET B.C.S TO MODE */
	for(int i=0;i<BASE::nebd;++i)
		for(int sm=0;sm<basis::tri(BASE::log2p)->sm();++sm)
			BASE::hp_ebdry(i)->sdirichlet(sm);


	BASE::ug.v(Range(0,BASE::npnt-1)) -= BASE::gbl->res.v(Range(0,BASE::npnt-1));
	BASE::ug.s(Range(0,BASE::nseg-1)) -= BASE::gbl->res.s(Range(0,BASE::nseg-1));
	BASE::ug.i(Range(0,BASE::ntri-1)) -= BASE::gbl->res.i(Range(0,BASE::ntri-1));

	return;
}

template<class BASE>void pod_simulate<BASE>::sc0load() {
	for(int i=0;i<BASE::nebd;++i)
		BASE::ebdry(i)->comm_prepare(boundary::all,0,boundary::symmetric);
	return;
}

template<class BASE> int pod_simulate<BASE>::sc0wait_rcv() {
	int stop = 1;
	for(int i=0;i<BASE::nebd;++i) {
		stop &= BASE::ebdry(i)->comm_wait(boundary::all,0,boundary::symmetric);
	}
	return(stop);
}


template<class BASE> FLT pod_simulate<BASE>::maxres() {
    int i;
    FLT mxr;

    mxr = 0.0;

    for(i=0;i<nmodes;++i)
		mxr = MAX(fabs(rsdls_recv(i)),mxr);

    *BASE::gbl->log << ' ' << mxr << ' ';

    return(mxr);
}

template<class BASE> void pod_sim_edge_bdry<BASE>::init(input_map& input) {
    std::string filename,keyword,linebuff;
    std::ostringstream nstr;
    std::istringstream instr;
    int i;

	modes.resize(nmodes);
    for(i=0;i<nmodes;++i) {
		nstr.str("");
		nstr << i << std::flush;
		filename = "mode" +nstr.str() +"_" +base.idprefix;
		nstr.clear();
		modes(i).v.resize(base.maxseg+1,x.NV);
		modes(i).s.resize(base.maxseg,x.sm0,x.NV);

		/* INPUT 1D MODE */
		filename = "mode" +nstr.str() + "_" + base.idprefix +".bin";
		binifstream bin;
		bin.open(filename.c_str());
		if (bin.error()) {
			*x.gbl->log << "couldn't open input file " << filename << std::endl;
			exit(1);
		}
		bin.setFlag(binio::BigEndian,bin.readInt(1));
		bin.setFlag(binio::FloatIEEE,bin.readInt(1));

		for (int bsind=0;bsind<base.nseg;++bsind)
			for (int n=0;n<x.NV;++n)
				modes(i).v(bsind,n) = bin.readFloat(binio::Double);
		for (int n=0;n<x.NV;++n)
			modes(i).v(base.nseg,n) = bin.readFloat(binio::Double);

		for (int bsind=0;bsind<base.nseg;++bsind) {
			for (int m=0;m<x.sm0;++m)
				for (int n=0;n<x.NV;++n)
					modes(i).s(bsind,m,n) = bin.readFloat(binio::Double);
		}
		bin.close();
    }

    int initfile;
    input.getwdefault("initfile",initfile,1);
    nstr.str("");
    nstr << initfile << std::flush;
    filename = "coeff" +nstr.str() +"_" +base.idprefix +".bin";
    binifstream bin;
    bin.open(filename.c_str());
    if (bin.error()) {
		*x.gbl->log << "couldn't open coefficient input file " << filename << std::endl;
		exit(1);
    }
	bin.setFlag(binio::BigEndian,bin.readInt(1));
	bin.setFlag(binio::FloatIEEE,bin.readInt(1));

    /* CONSTRUCT INITIAL SOLUTION DESCRIPTION */
    for (int l=0;l<nmodes;++l) {
		x.coeffs(bindex+l) = bin.readFloat(binio::Double); 
	}
	bin.close();

	addto2Dsolution(x.ug);

    return;
}


template<class BASE> void pod_sim_edge_bdry<BASE>::loadbuff(Array<FLT,1>& sdata) {
	if (!active) return;

	for(int i=0;i<nmodes;++i) {
		base.fsndbuf(i) = sdata(bindex+i);
	}
	base.sndsize() = nmodes;
	base.sndtype() = boundary::flt_msg;
}

template<class BASE> void pod_sim_edge_bdry<BASE>::finalrcv(Array<FLT,1>& sdata) {
	int j;

	if (!active) return;

    bool reload = base.comm_finish(boundary::all,0,boundary::symmetric,boundary::average);
	if (!reload) return;

    for(j=0;j<nmodes;++j) {
		sdata(bindex+j) = base.fsndbuf(j);
    }
	return;
}

template<class BASE> void pod_sim_edge_bdry<BASE>::addto2Dsolution(struct tri_hp::vsi ug) {

	if (!active) return;

	int sind, v0;

	for (int l=0;l<nmodes;++l) {
		int bsind = 0;
		do {
			sind = base.seg(bsind);
			v0 = x.seg(sind).pnt(0);
			ug.v(v0) += x.coeffs(bindex+l)*modes(l).v(bsind);
		} while (++bsind < base.nseg);
		v0 = x.seg(sind).pnt(1);
		ug.v(v0) += x.coeffs(bindex+l)*modes(l).v(base.nseg);

		for (int bsind=0;bsind<base.nseg;++bsind) {
			sind = base.seg(bsind);
			ug.s(sind) += x.coeffs(bindex+l)*modes(l).s(bsind);
		}
	}
	return;
}

template<class BASE> void pod_sim_edge_bdry<BASE>::addto2Dsolution(struct tri_hp::vsi ug, int mode, FLT coeff) {

	if (!active) return;

	int lmode = mode -bindex;
	if (lmode < 0 || lmode >= nmodes) return;

	int sind, v0;

	int bsind = 0;
	do {
		sind = base.seg(bsind);
		v0 = x.seg(sind).pnt(0);
		ug.v(v0) += coeff*modes(lmode).v(bsind);
	} while(++bsind < base.nseg);
	v0 = x.seg(sind).pnt(1);
	ug.v(v0) += coeff*modes(lmode).v(base.nseg);

	for (int bsind=0;bsind<base.nseg;++bsind) {
		sind = base.seg(bsind);
		ug.s(sind) += coeff*modes(lmode).s(bsind);
	}

	return;
}




template<class BASE> void pod_sim_edge_bdry<BASE>::rsdl() {

	if (!active) return;

	int sind,v0;

	for (int k = 0; k < nmodes; ++k) {		
		int bsind = 0;
		do {
			sind = base.seg(bsind);
			v0 = x.seg(sind).pnt(0);
			for (int n=0;n<x.NV;++n)
				x.rsdls(bindex+k) += modes(k).v(bsind,n)*x.gbl->res.v(v0,n);
		} while (++bsind < base.nseg);
		v0 = x.seg(sind).pnt(1);
		for (int n=0;n<x.NV;++n)
			x.rsdls(bindex+k) += modes(k).v(base.nseg,n)*x.gbl->res.v(v0,n);

		for (int bsind=0;bsind<base.nseg;++bsind) {
			sind = base.seg(bsind);
			for(int sm=0;sm<basis::tri(x.log2p)->sm();++sm)
				for(int n=0;n<x.NV;++n)
					x.rsdls(bindex+k) += modes(k).s(bsind,sm,n)*x.gbl->res.s(sind,sm,n);
		}		
    }

    return;
}





