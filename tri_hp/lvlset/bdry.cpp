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
			return(tri_hp::getnewvrtxobject(bnum,bdrydata));
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
			temp = new characteristic<bdry_ins::generic>(*this,*ebdry(bnum));
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
			return(tri_hp_ins::getnewsideobject(bnum,bdrydata));
			break;
		}
	}    
	gbl->ebdry_gbls(bnum) = temp->create_global_structure();

	return(temp);
}


/* THIS CLASS IS FOR UNSTEADY FROZEN FLOW STUFF */
class reinitialize_flow : public tri_hp_helper {
	public:
		reinitialize_flow(tri_hp &xin) :tri_hp_helper(xin) {}
		void update(int stage) {
			const int npnt = x.npnt, nseg = x.nseg, ntri = x.ntri;
			const int sm0 = x.sm0, im0 = x.im0;
			
			/* This is for an unsteady frozen flow */
			/* Store current value of level-set */
			x.ugbd(1).v(Range(0,npnt-1),0) = x.ug.v(Range(0,npnt-1),2);
			if (x.sm0) x.ugbd(1).s(Range(0,nseg-1),Range(0,sm0-1),0) = x.ug.s(Range(0,nseg-1),Range(0,sm0-1),2);
			if (x.im0) x.ugbd(1).i(Range(0,ntri-1),Range(0,im0-1),0) = x.ug.i(Range(0,ntri-1),Range(0,im0-1),2);
			x.tobasis(x.gbl->ibc);
			x.ug.v(Range(0,npnt-1),2) = x.ugbd(1).v(Range(0,npnt-1),0); 
			if (x.sm0) x.ug.s(Range(0,nseg-1),Range(0,sm0-1),2) = x.ugbd(1).s(Range(0,nseg-1),Range(0,sm0-1),0);
			if (x.im0) x.ug.i(Range(0,ntri-1),Range(0,im0-1),2) = x.ugbd(1).i(Range(0,ntri-1),Range(0,im0-1),0);
		}
};

class tri_hp_lvlset_helper_type {
	public:
		const static int ntypes = 1;
		enum ids {reinitialize};
		const static char names[ntypes][40];
		static int getid(const char *nin) {
			int i;
			for(i=0;i<ntypes;++i) 
				if (!strcmp(nin,names[i])) return(i);
			return(-1);
		}
};
const char tri_hp_lvlset_helper_type::names[ntypes][40] = {"reinitialize"};


tri_hp_helper *tri_hp_lvlset::getnewhelper(input_map& inmap) {
	std::string movername;
	int type;
	
	/* FIND INITIAL CONDITION TYPE */
	if (!inmap.get(gbl->idprefix + "_helper",movername))
		inmap.getwdefault("tri_hp_helper",movername,std::string("default"));
	
	type = tri_hp_lvlset_helper_type::getid(movername.c_str());
	
	switch(type) {
		case tri_hp_lvlset_helper_type::reinitialize: {
			tri_hp_helper *temp = new reinitialize_flow(*this);
			return(temp);
		}
		default: {
			tri_hp_helper *temp = new tri_hp_helper(*this);
			return(temp);
		}
	}
}


void hybrid::pmatchsolution_snd(int phase, FLT *pdata, int stride) {
	int j,k,count,sind(0),offset;

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
	int j,k,m,count,countdn,countup,offset,sind(0);
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
	
	int sendsize(8);
	base.sndsize() = sendsize;
	base.sndtype() = boundary::flt_msg;
	// see /ins/bdry.cpp hybrid_pt::update
	// info from moving-mesh
	base.fsndbuf(0) = 0.0;
	base.fsndbuf(1) = 0.0;
	base.fsndbuf(2) = 0.0;
	// (3) = flow direction
	// (4) = x location
	// (5) = y location
	base.fsndbuf(3) = 1.0;
	base.fsndbuf(4) = 0.0;
	base.fsndbuf(5) = 0.0;
	base.fsndbuf(6) = 0.0;
	base.fsndbuf(7) = 0.0;

	int sind = x.ebdry(base.ebdry(1))->seg(0);
	int v0 = x.seg(sind).pnt(0);
	int v1 = x.seg(sind).pnt(1);

	TinyVector<FLT,2> nrm, vel;
	FLT normvel;

	/* Normal to hybrid boundary */
	nrm(0) =  (x.pnts(v1)(1) -x.pnts(v0)(1));
	nrm(1) = -(x.pnts(v1)(0) -x.pnts(v0)(0));

	vel(0) = x.ug.v(v0,0)-(x.gbl->bd(0)*(x.pnts(v0)(0) -x.vrtxbd(1)(v0)(0)));
	vel(1) = x.ug.v(v0,1)-(x.gbl->bd(0)*(x.pnts(v0)(1) -x.vrtxbd(1)(v0)(1)));
#ifdef MESH_REF_VEL
	vel -= x.gbl->mesh_ref_vel;
#endif

	/* normvel is defined positive outward */
	normvel = vel(0)*nrm(0)+vel(1)*nrm(1);

	if (normvel > 0.0)  { 
		/* flow going out level-set changed from zero */
		int p0 = x.seg(x.ebdry(base.ebdry(1))->seg(0)).pnt(1);
		int p1 = x.seg(x.ebdry(base.ebdry(0))->seg(x.ebdry(base.ebdry(0))->nseg-1)).pnt(0);

		/* MOVE TO INTERSECTION POINT */
		if (x.ug.v(base.pnt,2)*x.ug.v(p0,2) < 0.0) {
			/* Between base.pnt and point above */
			TinyVector<FLT,2> dx = x.pnts(p0)-x.pnts(base.pnt);
			dx *= (0.0-x.ug.v(base.pnt,2))/(x.ug.v(p0,2)-x.ug.v(base.pnt,2));
			x.pnts(base.pnt) += dx;
			x.ug.v(base.pnt,2) = 0.0;
		}
		else {
			/* Between base.pnt and point below */
			TinyVector<FLT,2> dx = x.pnts(p1)-x.pnts(base.pnt);
			dx *= (0.0-x.ug.v(base.pnt,2))/(x.ug.v(p1,2)-x.ug.v(base.pnt,2));
			x.pnts(base.pnt) += dx;
			x.ug.v(base.pnt,2) = 0.0;	
		}
		base.fsndbuf(3) = -1.0;
	}
	// move the approximation to the boundary if it an analytically defined boundary
	// x.ebdry(base.ebdry(1))->mvpttobdry(0,-1.0,x.pnts(base.pnt));
	base.fsndbuf(4) = x.pnts(base.pnt)(0);
	base.fsndbuf(5) = x.pnts(base.pnt)(1);

	base.comm_prepare(boundary::all,0,boundary::symmetric);
	base.comm_exchange(boundary::all,0,boundary::symmetric);
	base.comm_wait(boundary::all,0,boundary::symmetric);

	// add information from other points
	for(int m=0;m<base.nmatches();++m) {
		for(int i=0;i<sendsize;++i) 
			base.fsndbuf(i) += base.frcvbuf(m,i);
	}

	if (base.fsndbuf(0)*base.fsndbuf(3) > 0.0) {
		*x.gbl->log << "uh-oh opposite characteristics at hybrid point" << std::endl;
		*x.gbl->log << "local "  << base.idprefix << ' ' << base.fsndbuf(0) << "remote " << base.fsndbuf(3) << std::endl;
	}
	
	if (base.fsndbuf(0) < 0.0) {
		// flow is into lvl-set
		// Set point to location sent from moving mesh
		x.pnts(base.pnt)(0) = base.fsndbuf(1);
		x.pnts(base.pnt)(1) = base.fsndbuf(2);
		
		// Load tangent sent from moving mesh 
		// this will be used to determine slope of levelset along hybrid inflow boundaries
		// Had to do this check because hybrid update occurs before hybrid_pt update.  With
		// Multigrid tang got messed up by coarse log2p levels 
		if (x.log2p == x.log2pmax) {
			tang(0) = base.fsndbuf(6);
			tang(1) = base.fsndbuf(7);
			tangset = true;
		}
	}
}

void hybrid::update(int stage) {
	int sind(0),tind(0),v0(0);
	bool flag(false);
	int v1,v2;
	FLT xloc, yloc;
	FLT x2,y2,velx,vely,nrmx,nrmy;
	FLT normux, normuy, lengthu;
	TinyVector<FLT,tri_mesh::ND> tang;
	
	//  return;  // Uncomment this to simply freeze values on incoming boundaries

	if (stage == -1 || x.coarse_flag) return;

	// Find hybrid point //
	if (base.vbdry(0) > -1) {
		sind = base.seg(0);
		flag = true;
	}
	else if (base.vbdry(1) > -1) {
		sind = base.seg(base.nseg-1);
		flag = false;
	}
	else {
		*x.gbl->log << "Neither was a hybrid_pt?  What gives?" << std::endl;
		sim::abort(__LINE__,__FILE__,x.gbl->log);
	}
	
	hybrid_pt *hpt;
	if ((hpt = dynamic_cast<hybrid_pt *>(x.hp_vbdry(base.vbdry(1-flag)))) == NULL) {
		*x.gbl->log << "Neither was a hybrid_pt?  What gives? " << hpt << std::endl;
		sim::abort(__LINE__,__FILE__,x.gbl->log);
	}

	tind = x.seg(sind).tri(0);
	v0 = x.seg(sind).pnt(1-flag);
	v2 = x.seg(sind).pnt(flag);
	
	/* Load surface tangent from hybrid_pt */
	tang(0) = hpt->tang(0);
	tang(1) = hpt->tang(1);
	bool tangset = hpt->tangset;
	if (!tangset) return;
	
	/* normal to levelset */
	lengthu = sqrt(tang(0)*tang(0)+tang(1)*tang(1));
	normux = tang(1)/lengthu;
	normuy = -tang(0)/lengthu;
	
	/* get sign right? */
	TinyVector<FLT,tri_mesh::ND> dx = x.pnts(v2)-x.pnts(v0);
	if ((dx(0)*normux +dx(1)*normuy)*x.ug.v(v2,2) < 0.0) {
		normux *= -1;
		normuy *= -1;
	}	
	
	/* Location of hybrid point */
	xloc = x.pnts(v0)(0);
	yloc = x.pnts(v0)(1);

	/* Re-evaluate phi at inflow points */
	for(int j=0;j<base.nseg;++j) {
		/* Calculate normal velocity to side */
		// normvel is defined positive outward //
		sind = base.seg(j);
		v1 = x.seg(sind).pnt(0);
		v2 = x.seg(sind).pnt(1);
		nrmx =  (x.pnts(v2)(1) -x.pnts(v1)(1));
		nrmy = -(x.pnts(v2)(0) -x.pnts(v1)(0));
		velx = 0.5*(x.ug.v(v0,0)-(x.gbl->bd(0)*(x.pnts(v0)(0) -x.vrtxbd(1)(v0)(0))) +
				  x.ug.v(v1,0)-(x.gbl->bd(0)*(x.pnts(v1)(0) -x.vrtxbd(1)(v1)(0))));
		vely = 0.5*(x.ug.v(v0,1)-(x.gbl->bd(0)*(x.pnts(v0)(1) -x.vrtxbd(1)(v0)(1))) +
				  x.ug.v(v1,1)-(x.gbl->bd(0)*(x.pnts(v1)(1) -x.vrtxbd(1)(v1)(1))));
#ifdef MESH_REF_VEL
		velx -= x.gbl->mesh_ref_vel(0);
		vely -= x.gbl->mesh_ref_vel(1);
#endif
		FLT normvel = velx*nrmx+vely*nrmy;
		
		if (normvel < 0.0 ) {
			v2 = x.seg(sind).pnt(flag);
			x2 = x.pnts(v2)(0);
			y2 = x.pnts(v2)(1);

			// use these lines to apply phi along the boundary as the normal distance from the two-phase boundary tangent
			FLT temp = (x2-xloc)*normux + (y2-yloc)*normuy;
			if (fabs(x.ug.v(v2,2) - temp) > 0.1) {
				*x.gbl->log << "big change " << x2 << ' ' << y2 << std::endl;
				*x.gbl->log << "big change " << xloc << ' ' << yloc << std::endl;
				*x.gbl->log << "big change " << normux << ' ' << normuy << std::endl;
				x.output("bigchange",tri_hp::tecplot);
				// exit(1);
			}
			else {
				x.ug.v(v2,2) = temp;
			}
			// OR use this line to apply phi along the boundary as the vertical distance from the hybrid point
			// x.ug.v(v2,2) = y2 - yloc;
		}
	}

}
