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

void surface::init(input_map& inmap,void* gbl_in) {
	bdry_ins::surface::init(inmap,gbl_in);
	
	input_map zeromap;
	zeromap["zero"] = "0.0";

	std::ostringstream nstr;
	for (int n=0;n<x.NV;++n) {
		nstr.str("");
		nstr << base.idprefix << "_flux" << n << std::flush;
		if (inmap.find(nstr.str()) != inmap.end()) {
			fluxes(n).init(inmap,nstr.str());
		}
		else {
			fluxes(n).init(zeromap,"zero");
		}
	}
	
	return;
}



/* Free-Surface that Allows Radiative Heat Flux & Marangoni Effects */
void surface::element_rsdl(int indx, Array<FLT,2> lf) {
	int i,n,sind,v0,v1;
	TinyVector<FLT,tri_mesh::ND> norm, rp;
	Array<FLT,1> ubar(x.NV),flx(x.NV);
	FLT jcb;
	Array<TinyVector<FLT,MXGP>,1> u(x.NV);
	TinyMatrix<FLT,tri_mesh::ND,MXGP> crd, dcrd, mvel;
	TinyMatrix<FLT,9,MXGP> res;
	Array<FLT,1> axpt(tri_mesh::ND), amv(tri_mesh::ND), anorm(tri_mesh::ND),au(x.NV);

	sind = base.seg(indx);
	v0 = x.seg(sind).pnt(0);
	v1 = x.seg(sind).pnt(1);
	
	x.crdtocht1d(sind);
	for(n=0;n<tri_mesh::ND;++n)
		basis::tri(x.log2p)->proj1d(&x.cht(n,0),&crd(n,0),&dcrd(n,0));
	
	for(n=0;n<x.NV;++n)
		basis::tri(x.log2p)->proj1d(&x.uht(n)(0),&u(n)(0));    
	
	for(i=0;i<basis::tri(x.log2p)->gpx();++i) {
		norm(0) =  dcrd(1,i);
		norm(1) = -dcrd(0,i);
		jcb = sqrt(norm(0)*norm(0) +norm(1)*norm(1));
		
		/* RELATIVE VELOCITY STORED IN MVEL(N)*/
		for(n=0;n<tri_mesh::ND;++n) {
			mvel(n,i) = u(n)(i) -(x.gbl->bd(0)*(crd(n,i) -dxdt(x.log2p,indx)(n,i)));
		}
	
		/* Evaluate Fluxes */
		axpt(0) = crd(0,i); axpt(1) = crd(1,i);
		amv(0) = mvel(0,i)-u(0)(i); amv(1) = mvel(1,i)-u(1)(i);
		anorm(0)= norm(0); anorm(1) = norm(1);
		for(int n=0;n<x.NV;++n)
			au(n) = u(n)(i);
		for(int n=0;n<x.NV;++n) {
			flx(n) = fluxes(n).Eval(au,axpt,amv,anorm,x.gbl->time);
		}
			
		/* TANGENTIAL SPACING */                
		res(0,i) = -ksprg(indx)*jcb;
		/* NORMAL FLUX */
		res(1,i) = -RAD(crd(0,i))*(mvel(0,i)*norm(0) +mvel(1,i)*norm(1)) +flx(x.NV-1);     
		/* UPWINDING BASED ON TANGENTIAL VELOCITY */
		res(2,i) = -res(1,i)*(-norm(1)*mvel(0,i) +norm(0)*mvel(1,i))/jcb*gbl->meshc(indx);
		
		/* SURFACE TENSION SOURCE TERM X-DIRECTION */ 
		res(4,i) = +RAD(crd(0,i))*((x.gbl->rho -gbl->rho2)*x.gbl->g*crd(1,i) +gbl->p_ext)*norm(0) +flx(0)*jcb;
#ifdef AXISYMMETRIC
		res(4,i) += gbl->sigma*jcb;
#endif
		/* AND INTEGRATION BY PARTS TERM */
		res(5,i) = +RAD(crd(0,i))*gbl->sigma*norm(1)/jcb;
		
		
		/* SURFACE TENSION SOURCE TERM Y-DIRECTION */
		res(6,i) = +RAD(crd(0,i))*((x.gbl->rho -gbl->rho2)*x.gbl->g*crd(1,i) +gbl->p_ext)*norm(1) +flx(1)*jcb;                
		/* AND INTEGRATION BY PARTS TERM */
		res(7,i) = -RAD(crd(0,i))*gbl->sigma*norm(0)/jcb;
		
		/* Auxiliary Fluxes Here */
		res(8,i) = flx(2)*jcb;
	}
	
	lf = 0.0;
	/* INTEGRATE & STORE SURFACE TENSION SOURCE TERM */
	basis::tri(x.log2p)->intgrt1d(&lf(0,0),&res(4,0));
	basis::tri(x.log2p)->intgrtx1d(&lf(0,0),&res(5,0));
	basis::tri(x.log2p)->intgrt1d(&lf(1,0),&res(6,0));
	basis::tri(x.log2p)->intgrtx1d(&lf(1,0),&res(7,0));
	
	/* Heat Flux Source Term */
	basis::tri(x.log2p)->intgrt1d(&lf(2,0),&res(8,0));

	
	/* INTEGRATE & STORE MESH MOVEMENT RESIDUALS */                    
	basis::tri(x.log2p)->intgrtx1d(&lf(x.NV,0),&res(0,0));
	basis::tri(x.log2p)->intgrt1d(&lf(x.NV+1,0),&res(1,0));
	basis::tri(x.log2p)->intgrtx1d(&lf(x.NV+1,0),&res(2,0));
	
#ifndef petsc
	/* mass flux preconditioning */
	for(int m=0;m<basis::tri(x.log2p)->sm()+2;++m)
		lf(x.NV-1,m) = -x.gbl->rho*lf(x.NV+1,m); 
#ifndef INERTIALESS
	for (n=0;n<x.NV-1;++n) 
		ubar(n) = 0.5*(x.uht(n)(0) +x.uht(n)(1));
	
	for (n=0;n<x.NV-1;++n) {
		lf(n,0) -= x.uht(n)(0)*(x.gbl->rho -gbl->rho2)*lf(x.NV+1,0);
		lf(n,1) -= x.uht(n)(1)*(x.gbl->rho -gbl->rho2)*lf(x.NV+1,1);
		for(int m=0;m<basis::tri(x.log2p)->sm();++m)
			lf(n,m+2) -= ubar(n)*(x.gbl->rho -gbl->rho2)*lf(x.NV+1,m+2);
	}
#endif
#endif
	
#ifdef DROP
	basis::tri(x.log2p)->intgrt1d(&lf(x.NV+2,0),&res(3,0));
#endif
	
	return;
}

/** \brief Helper object for vrtx_bdry 
 *
 * \ingroup boundary
 * Contains list of all vrtx_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tri_hp_buoyancy_vtype {
public:
	static const int ntypes = 2;
	enum ids {unknown=-1,melt_end,melt_inflow};
	const static char names[ntypes][40];
	static int getid(const char *nin) {
		for(int i=0;i<ntypes;++i) 
			if (!strcmp(nin,names[i])) return(i);
		return(unknown);
	}
};

const char tri_hp_buoyancy_vtype::names[ntypes][40] = {"melt_end","melt_inflow"};

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
		default: {
			return(tri_hp_ins::getnewvrtxobject(bnum,bdrydata));
		}
	} 
	gbl->vbdry_gbls(bnum) = temp->create_global_structure();
	return(temp);
}



class tri_hp_buoyancy_stype {
	public:
		static const int ntypes = 2;
		enum ids {unknown=-1,surface,melt};
		static const char names[ntypes][40];
		static int getid(const char *nin) {
			for(int i=0;i<ntypes;++i)
				if (!strcmp(nin,names[i])) return(i);
			return(-1);
		}
};

const char tri_hp_buoyancy_stype::names[ntypes][40] = {"surface","melt"};

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
				std::cerr << "use coupled physics for surface boundary" << std::endl;
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

