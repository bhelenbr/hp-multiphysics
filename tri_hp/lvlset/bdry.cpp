#include "bdry_lvlset.h"
#include <myblas.h>

//#define MPDEBUG

/*************************************************/
/* SET DIRICHLET BOUNDARY VALUES & FLUXES ********/
/* (THINGS THAT ARE INDEPENDENT OF THE SOLUTION) */
/* BUT NOT INDEPENDENT OF TIME/MESH POSITION     */
/*************************************************/
using namespace bdry_lvlset;

/** \brief Helper object for vrtx_bdry 
 *
 * \ingroup boundary
 * Contains list of all vrtx_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tri_hp_lvlset_vtype {
	public:
		static const int ntypes = 1;
		enum ids {unknown=-1,hybrid_point};
		const static char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i) 
				if (!strcmp(nin,names[i])) return(i);
			return(unknown);
		}
};

const char tri_hp_lvlset_vtype::names[ntypes][40] = {"hybrid_point"};

hp_vrtx_bdry* tri_hp_lvlset::getnewvrtxobject(int bnum, input_map &bdrydata) {
	std::string keyword,val;
	std::istringstream data;
	int type;          
	hp_vrtx_bdry *temp;  

	keyword = vbdry(bnum)->idprefix + "_lvlset_type";
	if (bdrydata.get(keyword,val)) {
		type = tri_hp_lvlset_vtype::getid(val.c_str());
		if (type == tri_hp_lvlset_vtype::unknown)  {
			*gbl->log << "unknown vertex type:" << val << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
	}
	else {
		type = tri_hp_lvlset_vtype::unknown;
	}


	switch(type) {
		case tri_hp_lvlset_vtype::hybrid_point: {
			temp = new hybrid_pt(*this,*vbdry(bnum));
			break;
		}
		default: {
			temp = tri_hp::getnewvrtxobject(bnum,bdrydata);
			break;
		}
	} 
	gbl->vbdry_gbls(bnum) = temp->create_global_structure();
	return(temp);
}



class tri_hp_lvlset_stype {
	public:
		static const int ntypes = 5;
		enum ids {unknown=-1,inflow,outflow,characteristic,euler,hybrid};
		static const char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i)
				if (!strcmp(nin,names[i])) return(i);
			return(-1);
		}
};

const char tri_hp_lvlset_stype::names[ntypes][40] = {"inflow","outflow","characteristic","euler","hybrid"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_edge_bdry* tri_hp_lvlset::getnewsideobject(int bnum, input_map &bdrydata) {
	std::string keyword,val;
	std::istringstream data;
	int type;          
	hp_edge_bdry *temp;  

	if (bdrydata.get(ebdry(bnum)->idprefix + "_lvlset_type",val)) {
		type = tri_hp_lvlset_stype::getid(val.c_str());
		if (type == tri_hp_lvlset_stype::unknown)  {
			*gbl->log << "unknown side type:" << val << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
	}
	else {
		type = tri_hp_lvlset_stype::unknown;
	}

	switch(type) {
		case tri_hp_lvlset_stype::inflow: {
			temp = new characteristic<bdry_ins::inflow>(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_lvlset_stype::outflow: {
			temp = new characteristic<bdry_ins::neumann>(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_lvlset_stype::characteristic: {
			temp = new characteristic<bdry_ins::characteristic>(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_lvlset_stype::euler: {
			temp = new characteristic<bdry_ins::euler>(*this,*ebdry(bnum));
			break;
		}
		case tri_hp_lvlset_stype::hybrid: {
			temp = new hybrid(*this,*ebdry(bnum));
			break;
		}
		default: {
			temp = tri_hp_ins::getnewsideobject(bnum,bdrydata);
			break;
		}
	}    
	return(temp);
}

void hybrid::pmatchsolution_snd(int phase, FLT *pdata, int stride) {
	int j,k,count,sind,offset;

	stride*=4;


	count = 0;
	for(j=0;j<base.nseg;++j) {
		sind = base.seg(j);
		offset = x.seg(sind).pnt(0)*stride;
		for (k=0;k<x.ND;++k) {
			base.fsndbuf(count++) = pdata[offset+k];
		}
		base.fsndbuf(count++) = pdata[offset+x.NV-1];
	}
	offset = x.seg(sind).pnt(1)*stride;
	for (k=0;k<x.ND;++k) 
		base.fsndbuf(count++) = pdata[offset+k];
	base.fsndbuf(count++) = pdata[offset+x.NV-1];

	base.sndsize() = count;
	base.sndtype() = boundary::flt_msg;
}

void hybrid::pmatchsolution_rcv(int phi, FLT *pdata, int stride) {
	int j,k,m,count,countdn,countup,offset,sind;
	int matches = 1;
	FLT mtchinv;

	stride *= 4; // u,v,phi,p;
	int grp = 1;  // ALL_PHASED GROUP

	/* ASSUMES REVERSE ORDERING OF SIDES */
	/* WON'T WORK IN 3D */
	/* RELOAD FROM BUFFER */
	/* ELIMINATES V/S/F COUPLING IN ONE PHASE */
	/* FINALRCV SHOULD BE CALLED F,S,V ORDER (V HAS FINAL AUTHORITY) */    
	for(m=0;m<base.nmatches();++m) {    
		if (base.msg_phase(grp,m) != phi) continue;

		++matches;

		int ebp1 = 3;  // u,v,p but not phi
		countdn = base.nseg*ebp1;
		countup = 0;
		for(j=0;j<base.nseg+1;++j) {
			for(k=0;k<ebp1;++k) {
				base.fsndbuf(countup +k) += base.frcvbuf(m,countdn +k);
			}
			countup += ebp1;
			countdn -= ebp1;
		}
	}

	if (matches > 1) {
		mtchinv = 1./matches;

#ifdef MPDEBUG
		*x.gbl->log << "finalrcv"  << base.idnum << " " << base.is_frst() << std::endl;
#endif
		count = 0;
		for(j=0;j<base.nseg;++j) {
			sind = base.seg(j);
			offset = x.seg(sind).pnt(0)*stride;
			for (k=0;k<x.ND;++k) {
				pdata[offset+k] = base.fsndbuf(count++)*mtchinv;
#ifdef MPDEBUG
				*x.gbl->log << "\t" << pdata[offset+k] << std::endl;
#endif
			}
			pdata[offset+x.NV-1] = base.fsndbuf(count++)*mtchinv;
#ifdef MPDEBUG
			*x.gbl->log << "\t" << pdata[offset+x.NV-1] << std::endl;
#endif
		}
		offset = x.seg(sind).pnt(1)*stride;
		for (k=0;k<x.ND;++k) {
			pdata[offset+k] = base.fsndbuf(count++)*mtchinv;
#ifdef MPDEBUG
			*x.gbl->log << "\t" << pdata[offset+k] << std::endl;
#endif
		}
		pdata[offset+x.NV-1] = base.fsndbuf(count++)*mtchinv;
#ifdef MPDEBUG
		*x.gbl->log << "\t" << pdata[offset+x.NV-1] << std::endl;
#endif
	}
	return;
}

void hybrid::smatchsolution_snd(FLT *sdata, int bgn, int end, int stride) {
	int j,k,n,countup,offset;

	if (!base.is_comm()) return;

#ifdef MPDEBUG
	*x.gbl->log << base.idprefix << " In surface_snd"  << base.idnum << " " << base.is_frst() << std::endl;
#endif

	countup = 0;
	for(j=0;j<base.nseg;++j) {
		offset = base.seg(j)*stride*x.NV;
		for(k=bgn;k<=end;++k) {
			for(n=0;n<x.ND;++n) {
				base.fsndbuf(countup++) = sdata[offset +k*x.NV +n];
#ifdef MPDEBUG
				*x.gbl->log << "\t" << sdata[offset +k*x.NV +n] << std::endl;
#endif
			}
			base.fsndbuf(countup++) = sdata[offset +k*x.NV +x.NV-1];
#ifdef MPDEBUG
			*x.gbl->log << "\t" << sdata[offset +k*x.NV +x.NV-1] << std::endl;
#endif
		}
	}
	base.sndsize() = countup;
	base.sndtype() = boundary::flt_msg;
	return;
}


void hybrid::smatchsolution_rcv(FLT *sdata, int bgn, int end, int stride) {

	if (!base.is_comm()) return;

	/* CAN'T USE sfinalrcv BECAUSE OF CHANGING SIGNS */
	int j,k,m,n,count,countdn,countup,offset,sind,sign;
	FLT mtchinv;

	/* ASSUMES REVERSE ORDERING OF SIDES */
	/* WON'T WORK IN 3D */

	int matches = 1;

	int bgnsign = (bgn % 2 ? -1 : 1);

	/* RELOAD FROM BUFFER */
	/* ELIMINATES V/S/F COUPLING IN ONE PHASE */
	/* FINALRCV SHOULD BE CALLED F,S,V ORDER (V HAS FINAL AUTHORITY) */    
	for(m=0;m<base.nmatches();++m) {            
		++matches;

		int ebp1 = end-bgn+1;

		countdn = (base.nseg-1)*ebp1*(x.NV-1);
		countup = 0;
		for(j=0;j<base.nseg;++j) {
			sign = bgnsign;
			for(k=0;k<ebp1;++k) {
				for(n=0;n<x.NV-1;++n) {
					base.fsndbuf(countup++) += sign*base.frcvbuf(m,countdn++);
				}
				sign *= -1;
			}
			countdn -= 2*ebp1*(x.NV-1);
		}
	}

	if (matches > 1) {
		mtchinv = 1./matches;

#ifdef MPDEBUG
		*x.gbl->log << base.idprefix << " In surface_rcv"  << base.idnum << " " << base.is_frst() << std::endl;
#endif
		count = 0;
		for(j=0;j<base.nseg;++j) {
			sind = base.seg(j);
			offset = sind*stride*x.NV;
			for (k=bgn;k<=end;++k) {
				for(n=0;n<x.ND;++n) {
					sdata[offset +k*x.NV +n] = base.fsndbuf(count++)*mtchinv;
#ifdef MPDEBUG
					*x.gbl->log << "\t" << sdata[offset +k*x.NV +n] << std::endl;
#endif
				}
				sdata[offset +k*x.NV +x.NV-1] = base.fsndbuf(count++)*mtchinv;
#ifdef MPDEBUG
				*x.gbl->log << "\t" << sdata[offset +k*x.NV +x.NV-1] << std::endl;
#endif				
			}

		}
	}
	return;
}

void hybrid_pt::pmatchsolution_snd(int phase, FLT *pdata, int stride) {
    int count,offset;

	if (!base.is_comm()) return;

	offset = base.pnt*4;
	count = 0;
    base.fsndbuf(count++) = pdata[offset];
	base.fsndbuf(count++) = pdata[offset+1];
	base.fsndbuf(count++) = pdata[offset+3];

    base.sndsize() = count;
    base.sndtype() = boundary::flt_msg;
}

void hybrid_pt::pmatchsolution_rcv(int phi, FLT *pdata, int stride) {
	int m,matches,k,offset,count;
	FLT mtchinv;

	if (!base.is_comm()) return;

	stride *= 4; // u,v,phi,p;
	int grp = 1;  // ALL_PHASED GROUP
	matches = 1;

	for(m=0;m<base.nmatches();++m) {    
		if (base.msg_phase(grp,m) != phi) continue;
		++matches;
		for(k=0;k<3;++k) {
			base.fsndbuf(k) += base.frcvbuf(m,k);
		}
	}

	if (matches > 1) {
		mtchinv = 1./matches;

#ifdef MPDEBUG
		*x.gbl->log << "finalrcv"  << base.idnum << " " << base.is_frst() << std::endl;
#endif

		offset = base.pnt*stride;
		count = 0;
		for (k=0;k<x.ND;++k) {
			pdata[offset+k] = base.fsndbuf(count++)*mtchinv;
#ifdef MPDEBUG
			*x.gbl->log << "\t" << pdata[offset+k] << std::endl;
#endif
		}
		pdata[offset+x.NV-1] = base.fsndbuf(count++)*mtchinv;
#ifdef MPDEBUG
		*x.gbl->log << "\t" << pdata[offset+x.NV-1] << std::endl;
#endif
	}
}


void hybrid_pt::update(int stage) {

	if (stage == -1) return;

	base.sndsize() = 4;
	base.sndtype() = boundary::flt_msg;
	base.fsndbuf(0) = 0.0;
	base.fsndbuf(1) = 0.0;
	base.fsndbuf(2) = 1.0;
	base.fsndbuf(3) = 0.0;

	int sind = x.ebdry(base.ebdry(1))->seg(0);
	int v0 = x.seg(sind).pnt(0);
	int v1 = x.seg(sind).pnt(1);

	TinyVector<FLT,2> nrm, vel;
	FLT normvel;

	nrm(0) =  (x.pnts(v1)(1) -x.pnts(v0)(1));
	nrm(1) = -(x.pnts(v1)(0) -x.pnts(v0)(0));

	vel(0) = 0.5*(x.ug.v(v0,0)-(x.gbl->bd(0)*(x.pnts(v0)(0) -x.vrtxbd(1)(v0)(0))) +
				  x.ug.v(v1,0)-(x.gbl->bd(0)*(x.pnts(v1)(0) -x.vrtxbd(1)(v1)(0))));
	vel(1) = 0.5*(x.ug.v(v0,1)-(x.gbl->bd(0)*(x.pnts(v0)(1) -x.vrtxbd(1)(v0)(1))) +
				  x.ug.v(v1,1)-(x.gbl->bd(0)*(x.pnts(v1)(1) -x.vrtxbd(1)(v1)(1))));

	/* normvel is defined positive outward */
	normvel = vel(0)*nrm(0)+vel(1)*nrm(1);

	if (normvel > 0.0)  { 
		/* flow going out level-set changed from zero */
		int p0 = x.seg(x.ebdry(base.ebdry(1))->seg(0)).pnt(1);
		int p1 = x.seg(x.ebdry(base.ebdry(0))->seg(x.ebdry(base.ebdry(0))->nseg-1)).pnt(0);

		if (x.ug.v(base.pnt,2)*x.ug.v(p0,2) < 0.0) {
			/* MOVE TO INTERSECTION POINT */
			TinyVector<FLT,2> dx = x.pnts(p0)-x.pnts(base.pnt);
			dx *= (0.0-x.ug.v(base.pnt,2))/(x.ug.v(p0,2)-x.ug.v(base.pnt,2));
			x.pnts(base.pnt)(1) += dx(1);
			x.ug.v(base.pnt,2) = 0.0;
		}
		else {
			/* MOVE TO INTERSECTION POINT */
			TinyVector<FLT,2> dx = x.pnts(p1)-x.pnts(base.pnt);
			dx *= (0.0-x.ug.v(base.pnt,2))/(x.ug.v(p1,2)-x.ug.v(base.pnt,2));
			x.pnts(base.pnt)(1) += dx(1);
			x.ug.v(base.pnt,2) = 0.0;	
		}
		base.fsndbuf(2) = -1.0;
	}
	base.fsndbuf(3) = x.pnts(base.pnt)(1);

	if (!base.is_comm()) return;

	base.comm_prepare(boundary::all,0,boundary::symmetric);
	base.comm_exchange(boundary::all,0,boundary::symmetric);
	base.comm_wait(boundary::all,0,boundary::symmetric);

	for(int m=0;m<base.nmatches();++m) {
		for(int i=0;i<4;++i) 
			base.fsndbuf(i) += base.frcvbuf(m,i);
	}

	if (base.fsndbuf(0)*base.fsndbuf(2) > 0.0) {
		*x.gbl->log << "uh-oh opposite characteristics at hybrid point" << std::endl;
		*x.gbl->log << "local "  << base.idprefix << ' ' << base.fsndbuf(0) << "remote " << base.fsndbuf(2) << std::endl;
	}

	if (base.fsndbuf(0) > 0.0) {
		x.pnts(base.pnt)(1) = base.fsndbuf(3);
	}
	else {
		x.pnts(base.pnt)(1) = base.fsndbuf(1);
	}
}
