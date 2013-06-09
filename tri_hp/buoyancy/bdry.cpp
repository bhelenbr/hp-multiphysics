/*
 *  bdry.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 12/13/10.
 *  Copyright 2010 Clarkson University. All rights reserved.
 *
 */

#include "bdry_buoyancy.h"
#include <tri_boundary.h>

using namespace bdry_buoyancy;
#include "bdry_buoyancy.h"
#include <myblas.h>

//#define MPDEBUG

//#define DEBUG

using namespace bdry_buoyancy;

void solid_fluid::init(input_map& inmap,void* gbl_in) {
	std::string keyword,val;
	std::istringstream data;
	std::string filename;
	
	/* Load in the additional flux for radiation if solid/gas boundary */
	keyword = base.idprefix + "_hp_typelist";
	if (inmap.find(base.idprefix+"_hp_typelist") == inmap.end()) {
		inmap[base.idprefix+"_hp_typelist"] = "0 0 1 1";
	}
	if (inmap.find(base.idprefix+"_flux2") == inmap.end()) {
		inmap[base.idprefix+"_flux2"] = "0.0)";
	}
	if (inmap.find(base.idprefix+"_flux3") == inmap.end()) {
		inmap[base.idprefix+"_flux3"] = "rho*(u0*n0 +u1*n1)";
	}
	
	symbolic::init(inmap,gbl_in);
	
	/* Now that we have loaded flux, make pressure neumann */
	essential_indices.clear();
	for(int n=0;n<x.ND;++n)
		essential_indices.push_back(n);
	type(Range(0,x.ND-1)) = essential;
	type(Range(x.ND,x.NV-1)) = natural;
		
	return;
}

void solid_fluid::smatchsolution_snd(FLT *sdata, int bgn, int end, int stride) {
	int j,k,n,countup,offset;
	
	if (!base.is_comm()) return;
	
#ifdef MPDEBUG
	*x.gbl->log << base.idprefix << " In melt::smatchsolution_snd"  << base.idnum << " " << base.is_frst() << std::endl;
#endif
	
	/* This boundary should always be first so send going up */
	countup = 0;
	for(j=0;j<base.nseg;++j) {
		offset = base.seg(j)*stride*x.NV;
		for(k=bgn;k<=end;++k) {
			for(n=x.ND;n<x.ND+1;++n) {
				base.fsndbuf(countup++) = sdata[offset +k*x.NV +n];
#ifdef MPDEBUG
				*x.gbl->log << "\t" << sdata[offset +k*x.NV +n] << std::endl;
#endif
			}
		}
	}
	base.sndsize() = countup;
	base.sndtype() = boundary::flt_msg;
	return;
}

void solid_fluid::smatchsolution_rcv(FLT *sdata, int bgn, int end, int stride) {
	
	if (!base.is_comm()) return;
	
	/* CAN'T USE sfinalrcv BECAUSE OF CHANGING SIGNS */
	int j,k,m,n,count,countup,offset,sind;
	FLT mtchinv;
	
	/* OPPOSING BOUNDARY SENDS BACKWARDS SO THIS CAN GO UP */
	
	int matches = 1;
	
	/* RELOAD FROM BUFFER */
	/* ELIMINATES V/S/F COUPLING IN ONE PHASE */
	/* FINALRCV SHOULD BE CALLED F,S,V ORDER (V HAS FINAL AUTHORITY) */
	for(m=0;m<base.nmatches();++m) {
		++matches;
		
		int ebp1 = end-bgn+1;
		
		countup = 0;
		for(j=0;j<base.nseg;++j) {
			for(k=0;k<ebp1;++k) {
				for(n=0;n<1;++n) {
					base.fsndbuf(countup) += base.frcvbuf(m,countup);
					++countup;
				}
			}
		}
	}
	
	if (matches > 1) {
		mtchinv = 1./matches;
		
#ifdef MPDEBUG
		*x.gbl->log << base.idprefix << "melt::smatchsolution_rcv"  << base.idnum << " " << base.is_frst() << std::endl;
#endif
		count = 0;
		for(j=0;j<base.nseg;++j) {
			sind = base.seg(j);
			offset = sind*stride*x.NV;
			for (k=bgn;k<=end;++k) {
				for(n=x.ND;n<x.ND+1;++n) {
					sdata[offset +k*x.NV +n] = base.fsndbuf(count++)*mtchinv;
#ifdef MPDEBUG
					*x.gbl->log << "\t" << sdata[offset +k*x.NV +n] << std::endl;
#endif
				}
			}
			
		}
	}
	return;
}

#ifdef petsc
void solid_fluid::petsc_matchjacobian_snd() {
	
	int vdofs;
	if (x.mmovement != x.coupled_deformable)
		vdofs = x.NV;
	else
		vdofs = x.NV+x.ND;
	
	int row,sind=-2;
	
	/* I am cheating here and sending floats and int's together */
	/* Send Jacobian entries */
	base.sndsize() = 0;
	base.sndtype() = boundary::flt_msg;
	base.fsndbuf(base.sndsize()++) = x.jacobian_start +FLT_EPSILON;
	
	/* First send number of entries for each vertex row */
	/* then append column numbers & values */
	for(int i=0;i<base.nseg;++i) {
		sind = base.seg(i);
		row = x.seg(sind).pnt(0)*vdofs+2;
		/* attach diagonal column # to allow continuity enforcement */
		base.fsndbuf(base.sndsize()++) = row +FLT_EPSILON;
		/* Send Temperature Equation */
		base.fsndbuf(base.sndsize()++) = x.J._cpt(row+1) -x.J._cpt(row) +FLT_EPSILON;
#ifdef MPDEBUG
		*x.gbl->log << "sending " << x.J._cpt(row+1) -x.J._cpt(row) << " jacobian entries for vertex " << row/vdofs << " and variable " << n << std::endl;
#endif
		for (int col=x.J._cpt(row);col<x.J._cpt(row+1);++col) {
#ifdef MPDEBUG
			*x.gbl->log << x.J._col(col) << ' ';
#endif
			base.fsndbuf(base.sndsize()++) = x.J._col(col) +FLT_EPSILON;
			base.fsndbuf(base.sndsize()++) = x.J._val(col);
		}
#ifdef MPDEBUG
		*x.gbl->log << std::endl;
#endif
		
		/* Send Side Information */
		row = x.npnt*vdofs +sind*x.NV*x.sm0 +2; // TEMPERATURE ROW
		/* attach diagonal column # to allow continuity enforcement */
		base.fsndbuf(base.sndsize()++) = row +FLT_EPSILON;
		for(int mode=0;mode<x.sm0;++mode) {
			base.fsndbuf(base.sndsize()++) = x.J._cpt(row+1) -x.J._cpt(row) +FLT_EPSILON;
#ifdef MPDEBUG
			*x.gbl->log << "sending " << x.J._cpt(row+1) -x.J._cpt(row) << " jacobian entries for vertex " << row/vdofs << " and variable " << n << std::endl;
#endif
			for (int col=x.J._cpt(row);col<x.J._cpt(row+1);++col) {
#ifdef MPDEBUG
				*x.gbl->log << x.J._col(col) << ' ';
#endif
				base.fsndbuf(base.sndsize()++) = x.J._col(col) +FLT_EPSILON;
				base.fsndbuf(base.sndsize()++) = x.J._val(col);
			}
			
#ifdef MPDEBUG
			*x.gbl->log << std::endl;
#endif
			row += x.NV;
		}
	}
	
	/* LAST POINT */
	row = x.seg(sind).pnt(1)*vdofs+2;
	/* attach diagonal # to allow continuity enforcement */
	base.fsndbuf(base.sndsize()++) = row +FLT_EPSILON;
	base.fsndbuf(base.sndsize()++) = x.J._cpt(row+1) -x.J._cpt(row) +FLT_EPSILON;
#ifdef MPDEBUG
	*x.gbl->log << "sending " << x.J._cpt(row+1) -x.J._cpt(row) << " jacobian entries for vertex " << row/vdofs << " and variable " << n << std::endl;
#endif
	for (int col=x.J._cpt(row);col<x.J._cpt(row+1);++col) {
#ifdef MPDEBUG
		*x.gbl->log << x.J._col(col) << ' ';
#endif
		base.fsndbuf(base.sndsize()++) = x.J._col(col) +FLT_EPSILON;
		base.fsndbuf(base.sndsize()++) = x.J._val(col);
	}
#ifdef MPDEBUG
	*x.gbl->log << std::endl;
#endif
}


void solid_fluid::petsc_matchjacobian_rcv(int phase) {
	const int ND = x.ND;
	
	if (base.matchphase(boundary::all_phased,0) != phase) return;
	
	int count = 0;
	int Jstart_mpi = static_cast<int>(base.frcvbuf(0, count++));
	
	sparse_row_major *pJ_mpi;
	if (base.is_local(0)) {
		pJ_mpi = &x.J;
		Jstart_mpi = 0;
	}
	else {
		pJ_mpi = &x.J_mpi;
	}
	
	int vdofs;
	if (x.mmovement != x.coupled_deformable)
		vdofs = x.NV;
	else
		vdofs = x.NV+x.ND;
	
	int sm = basis::tri(x.log2p)->sm();
	
	
	/* Now do stuff for communication boundaries */
	int row;
	
	std::vector<int> c0vars;
	c0vars.push_back(ND);
	
	/* Now Receive Information */
	for (int i=base.nseg-1;i>=0;--i) {
		int sind = base.seg(i);
		int rowbase = x.seg(sind).pnt(1)*vdofs;
		int row_mpi = static_cast<int>(base.frcvbuf(0,count++)) +Jstart_mpi;
		
		for(std::vector<int>::iterator it=c0vars.begin();it!=c0vars.end();++it) {
			row = rowbase + *it;
			int ncol = static_cast<int>(base.frcvbuf(0,count++));
#ifdef MPDEBUG
			*x.gbl->log << "receiving " << ncol << " jacobian entries for vertex " << row/vdofs << " and variable " << *it << std::endl;
#endif
			for (int k = 0;k<ncol;++k) {
				int col = static_cast<int>(base.frcvbuf(0,count++)) +Jstart_mpi;
#ifdef MPDEBUG
				*x.gbl->log  << col << ' ';
#endif
				FLT val = base.frcvbuf(0,count++);
				(*pJ_mpi).add_values(row,col,val);
			}
#ifdef MPDEBUG
			*x.gbl->log << std::endl;
#endif
			/* Shift all entries for this vertex.  Remote variables on solid are T */
			for(int n=0;n<1;++n) {
				FLT dval = x.J_mpi(row,row_mpi +n);
				(*pJ_mpi)(row,row_mpi +n) = 0.0;
				x.J(row,rowbase +c0vars[n]) += dval;
			}
			x.J.multiply_row(row,0.5);
			x.J_mpi.multiply_row(row,0.5);
		}
		
		/* Now receive side Jacobian information */
		row = x.npnt*vdofs +sind*x.NV*x.sm0;
		row_mpi = static_cast<int>(base.frcvbuf(0,count++)) +Jstart_mpi;
		
		int mcnt = x.ND;
		int sgn = 1;
		for(int mode=0;mode<x.sm0;++mode) {
			int ncol = static_cast<int>(base.frcvbuf(0,count++));
#ifdef MPDEBUG
			*x.gbl->log << "receiving " << ncol << " jacobian entries for side " << sind << " and variable " << 2 << std::endl;
#endif
			for (int k = 0;k<ncol;++k) {
				int col = static_cast<int>(base.frcvbuf(0,count++)) +Jstart_mpi;
#ifdef MPDEBUG
				*x.gbl->log  << col << ' ';
#endif
				FLT val = sgn*base.frcvbuf(0,count++);
				(*pJ_mpi).add_values(row+mcnt,col,val);
			}
#ifdef MPDEBUG
			*x.gbl->log << std::endl;
#endif
			
			/* Shift all modes in equation */
			int mcnt_mpi = 0;
			int sgn_mpi = 1;
			for(int mode_mpi=0;mode_mpi<x.sm0;++mode_mpi) {
				FLT dval = x.J_mpi(row+mcnt,row_mpi+mcnt_mpi);
				(*pJ_mpi)(row+mcnt,row_mpi+mcnt_mpi) = 0.0;
				x.J(row+mcnt,row+x.ND+mode_mpi*x.NV) += sgn_mpi*dval;
				sgn_mpi *= -1;
				mcnt_mpi += 1; // Only temperature on solid block
			}
			x.J.multiply_row(row+mcnt,0.5);
			x.J_mpi.multiply_row(row+mcnt,0.5);
			mcnt += x.NV;
			sgn *= -1;
		}
	}
	
	int sind = base.seg(0);
	int rowbase = x.seg(sind).pnt(0)*vdofs;
	int row_mpi = static_cast<int>(base.frcvbuf(0,count++)) +Jstart_mpi;
	
	for(std::vector<int>::iterator it=c0vars.begin();it!=c0vars.end();++it) {
		row = rowbase + *it;
		int ncol = static_cast<int>(base.frcvbuf(0,count++));
#ifdef MPDEBUG
		*x.gbl->log << "receiving " << ncol << " jacobian entries for vertex " << row/vdofs << " and variable " << *it << std::endl;
#endif
		for (int k = 0;k<ncol;++k) {
			int col = static_cast<int>(base.frcvbuf(0,count++)) +Jstart_mpi;
#ifdef MPDEBUG
			*x.gbl->log << col << ' ';
#endif
			FLT val = base.frcvbuf(0,count++);
			(*pJ_mpi).add_values(row,col,val);
		}
#ifdef MPDEBUG
		*x.gbl->log << std::endl;
#endif
		/* Shift all entries for this vertex.  Remote variables on solid are T, x, y */
		for(int n=0;n<1;++n) {
			FLT dval = x.J_mpi(row,row_mpi +n);
			(*pJ_mpi)(row,row_mpi +n) = 0.0;
			x.J(row,rowbase +c0vars[n]) += dval;
		}
		x.J.multiply_row(row,0.5);
		x.J_mpi.multiply_row(row,0.5);
	}
}

void solid_fluid::non_sparse_snd(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi) {
	const int sm=basis::tri(x.log2p)->sm();
	const int NV = x.NV;
	const int ND = tri_mesh::ND;
	
	int vdofs;
	if (x.mmovement != tri_hp::coupled_deformable)
		vdofs = NV;
	else
		vdofs = ND+NV;
	
	int begin_seg = x.npnt*vdofs;
	
	
	/* INTERFACE TEMPERATURE EQUATION JACOBIAN WILL BE EXCHANGED */
	/* IT WILL THEN BE REPLACED BY A DIRCHLET CONDITION AND MOVED TO X/Y EQUATION */
	std::vector<int> c0vars;
	c0vars.push_back(ND);
	for(int n=x.NV;n<vdofs;++n) {
		c0vars.push_back(n);
	}
	
	/* Going to send all jacobian entries,  Diagonal entries for matching DOF's will be merged together not individual */
	/* Send number of non-zeros to matches */
	base.sndsize() = 0;
	base.sndtype() = boundary::int_msg;
	for (int i=0;i<base.nseg;++i) {
		int sind = base.seg(i);
		int pind = x.seg(sind).pnt(0)*vdofs;
		for(std::vector<int>::iterator it=c0vars.begin();it!=c0vars.end();++it)
			base.isndbuf(base.sndsize()++) = nnzero(pind +*it);
	}
	int sind = base.seg(base.nseg-1);
	int pind = x.seg(sind).pnt(1)*vdofs;
	for(std::vector<int>::iterator it=c0vars.begin();it!=c0vars.end();++it)
		base.isndbuf(base.sndsize()++) = nnzero(pind +*it);
	
	/* Last thing to send is nnzero for edges */
	if (sm) {
		for (int i=0;i<base.nseg;++i) {
			int sind = base.seg(i);
			for (int m=0;m<sm;++m) {
				base.isndbuf(base.sndsize()++) = nnzero(begin_seg +sind*sm*NV +m*NV +2);
			}
		}
	}

	return;
}

void solid_fluid::non_sparse_rcv(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi) {
	
	if (!base.is_comm()) return;
	
	const int sm=basis::tri(x.log2p)->sm();
	const int NV = x.NV;
	const int ND = tri_mesh::ND;
	
	int vdofs;
	if (x.mmovement != tri_hp::coupled_deformable)
		vdofs = NV;
	else
		vdofs = ND+NV;
	
	int begin_seg = x.npnt*vdofs;
	
	std::vector<int> c0vars;
	c0vars.push_back(ND);
	for(int n=x.NV;n<vdofs;++n) {
		c0vars.push_back(n);
	}
	
	int count = 0;
	for (int i=base.nseg-1;i>=0;--i) {
		int sind = base.seg(i);
		int pind = x.seg(sind).pnt(1)*vdofs;
		int nentry = base.ircvbuf(0,count);
		for(std::vector<int>::iterator it=c0vars.begin();it!=c0vars.end();++it) {
			// T received from conv-diffusive side;
			nnzero_mpi(pind+*it) += nentry;
			++count; // Skip sizes sent for x & y rows (must match T size)
		}
		
	}
	int sind = base.seg(0);
	int pind = x.seg(sind).pnt(0)*vdofs;
	int nentry = base.ircvbuf(0,count);
	for(std::vector<int>::iterator it=c0vars.begin();it!=c0vars.end();++it) {
		// T received from conv-diffusive side;
		nnzero_mpi(pind+*it) += nentry;
		++count; // Skip sizes sent for x & y rows (must match T size)
	}
	
	
	/* Now add to side degrees of freedom */
	if (sm) {
		int toadd = base.ircvbuf(0,count++);
		for (int i=0;i<base.nseg;++i) {
			int sind = base.seg(i);
			for (int mode=0;mode<sm;++mode) {
				for(int n=x.ND;n<x.ND+1;++n) {
					nnzero_mpi(begin_seg+sind*NV*sm +mode*NV +n) += toadd;
				}
			}
		}
	}
}

#endif


/** \brief Helper object for vrtx_bdry 
 *
 * \ingroup boundary
 * Contains list of all vrtx_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tri_hp_buoyancy_vtype {
public:
	static const int ntypes = 3;
	enum ids {unknown=-1,melt_end,melt_inflow,melt_facet_pt};
	const static char names[ntypes][40];
	static int getid(const char *nin) {
		for(int i=0;i<ntypes;++i) 
			if (!strcmp(nin,names[i])) return(i);
		return(unknown);
	}
};

const char tri_hp_buoyancy_vtype::names[ntypes][40] = {"melt_end","melt_inflow","melt_facet_pt"};

hp_vrtx_bdry* tri_hp_buoyancy::getnewvrtxobject(int bnum, input_map &bdrydata) {
	std::string keyword,val;
	std::istringstream data;
	int type;          
	hp_vrtx_bdry *temp;  
	
	keyword = vbdry(bnum)->idprefix + "_buoyancy_type";
	if (bdrydata.get(keyword,val)) {
		type = tri_hp_buoyancy_vtype::getid(val.c_str());
		if (type == tri_hp_buoyancy_vtype::unknown)  {
			*gbl->log << "unknown vertex type:" << val << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
	}
	else {
		type = tri_hp_buoyancy_vtype::unknown;
	}
	
	
	switch(type) {
		case tri_hp_buoyancy_vtype::melt_end: {
			temp = new melt_end_pt(*this,*vbdry(bnum));
			break;
		}
		case tri_hp_buoyancy_vtype::melt_inflow: {
			temp = new melt_inflow_pt(*this,*vbdry(bnum));
			break;
		}
		case tri_hp_buoyancy_vtype::melt_facet_pt: {
			temp = new melt_facet_pt(*this,*vbdry(bnum));
			break;
		}
		default: {
			return(tri_hp_ins::getnewvrtxobject(bnum,bdrydata));
		}
	} 
	gbl->vbdry_gbls(bnum) = temp->create_global_structure();
	return(temp);
}



class tri_hp_buoyancy_stype {
	public:
		static const int ntypes = 5;
		enum ids {unknown=-1,surface,melt,melt_kinetics,kellerman,solid_fluid};
		static const char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i)
				if (!strcmp(nin,names[i])) return(i);
			return(-1);
		}
};

const char tri_hp_buoyancy_stype::names[ntypes][40] = {"surface","melt","melt_kinetics","kellerman","solid_fluid"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_edge_bdry* tri_hp_buoyancy::getnewsideobject(int bnum, input_map &bdrydata) {
	std::string keyword,val;
	std::istringstream data;
	int type;          
	hp_edge_bdry *temp;  
	
	if (bdrydata.get(ebdry(bnum)->idprefix + "_buoyancy_type",val)) {
		type = tri_hp_buoyancy_stype::getid(val.c_str());
		if (type == tri_hp_buoyancy_stype::unknown)  {
			*gbl->log << "unknown side type:" << val << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
	}
	else {
		type = tri_hp_buoyancy_stype::unknown;
	}
	
	switch(type) {
		case tri_hp_buoyancy_stype::surface: {
			if (dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))) {
				temp = new surface(*this,*ebdry(bnum));
				dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))->physics = temp;
			}
			else {
				std::cerr << "use coupled physics for surface boundary" << std::endl;
				sim::abort(__LINE__,__FILE__,&std::cerr);
				assert(0);
			}
			break;
		}
		case tri_hp_buoyancy_stype::melt: {
			if (dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))) {
				temp = new melt(*this,*ebdry(bnum));
				dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))->physics = temp;
			}
			else {
				std::cerr << "use coupled physics for melt boundary" << std::endl;
				sim::abort(__LINE__,__FILE__,&std::cerr);
				assert(0);
			}
			break;
		}
		case tri_hp_buoyancy_stype::melt_kinetics: {
			if (dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))) {
				temp = new melt_kinetics(*this,*ebdry(bnum));
				dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))->physics = temp;
			}
			else {
				std::cerr << "use coupled physics for melt_kinetics boundary" << std::endl;
				sim::abort(__LINE__,__FILE__,&std::cerr);
				assert(0);
			}
			break;
		}
		case tri_hp_buoyancy_stype::kellerman: {
			if (dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))) {
				temp = new kellerman(*this,*ebdry(bnum));
				dynamic_cast<ecoupled_physics_ptr *>(ebdry(bnum))->physics = temp;
			}
			else {
				std::cerr << "use coupled physics for kellerman boundary" << std::endl;
				sim::abort(__LINE__,__FILE__,&std::cerr);
				assert(0);
			}
			break;
		}
		case tri_hp_buoyancy_stype::solid_fluid: {
			if (dynamic_cast<eboundary_with_geometry<ecomm,symbolic_shape<tri_mesh::ND> > *>(ebdry(bnum))) {
				temp = new solid_fluid(*this,*ebdry(bnum));
			}
			else {
				std::cerr << "use symbolic_comm for solid_liquid boundary" << std::endl;
				sim::abort(__LINE__,__FILE__,&std::cerr);
				assert(0);
			}
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

