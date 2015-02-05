//
//  surface2.h
//  tri_hp
//
//  Created by Brian Helenbrook on 12/30/14.
//
//

#ifndef __tri_hp__surface2__
#define __tri_hp__surface2__

#include <iostream>
#include "../hp_coupled_boundary.h"
#include "tri_hp_ins.h"

//#define MPDEBUG
//#define DEBUG

class surface2 : public hp_deformable_bdry {
protected:
	tri_hp_ins &x;
public:
	struct global : public hp_deformable_bdry::global {
		/* FLUID PROPERTIES */
		FLT sigma,rho2,mu2,p_ext;
	} *gbl;
	
	void* create_global_structure() {return new global;}
	surface2(tri_hp_ins &xin, edge_bdry &bin) : hp_deformable_bdry(xin,bin), x(xin) {mytype = "surface2";}
	surface2(const surface2& inbdry, tri_hp_ins &xin, edge_bdry &bin)  : hp_deformable_bdry(inbdry,xin,bin), x(xin) {
		gbl = inbdry.gbl;
	};
	surface2* create(tri_hp& xin, edge_bdry &bin) const {return new surface2(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
	
	void init(input_map& input,void* gbl_in);
	void element_rsdl(int sind, Array<TinyVector<FLT,MXTM>,1> lf);
	void setup_preconditioner();
#ifndef petsc
	void rsdl(int stage); //!< applies mass flux preconditioner for multgrid case
#endif
};

class surface_fixed_pt2 : public hp_vrtx_bdry {
	/* INTERSECTING BOUNDARY CONTAINING END POINT MUST HAVE GEOMETRY NOT BE DEFINED SOLELY BY MESH */
protected:
	tri_hp_ins &x;
	surface2 *surf;
	int surfbdry;
	bool fix_norm;
	
public:
	surface_fixed_pt2(tri_hp_ins &xin, vrtx_bdry &bin) : hp_vrtx_bdry(xin,bin), x(xin), fix_norm(1) {mytype = "surface_fixed_pt2";}
	surface_fixed_pt2(const surface_fixed_pt2& inbdry, tri_hp_ins &xin, vrtx_bdry &bin) : hp_vrtx_bdry(inbdry,xin,bin), x(xin),
	surfbdry(inbdry.surfbdry),fix_norm(inbdry.fix_norm) {
		if (!(surf = dynamic_cast<surface2 *>(x.hp_ebdry(base.ebdry(surfbdry))))) {
			*x.gbl->log << "something's wrong can't find surface boundary" << std::endl;
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}
	}
	surface_fixed_pt2* create(tri_hp& xin, vrtx_bdry &bin) const {return new surface_fixed_pt2(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
	
	void init(input_map& inmap,void* gbl_in) {
		std::string keyword,val;
		std::istringstream data;
		std::string filename;
		
		hp_vrtx_bdry::init(inmap,gbl_in);
		
		if ((surf = dynamic_cast<surface2 *>(x.hp_ebdry(base.ebdry(0))))) {
			surfbdry = 0;
		}
		else if ((surf = dynamic_cast<surface2 *>(x.hp_ebdry(base.ebdry(1))))) {
			surfbdry = 1;
		}
		else {
			*x.gbl->log << "something's wrong neither side is a surface boundary" << std::endl;
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}
		fix_norm = 1;  // THIS IS ONLY TO BE CHANGE BY INHERITED CLASSES
	}
	
	void rsdl(int stage) {
		if (surfbdry == 0) {
			/* SET TANGENT RESIDUAL TO ZERO */
			surf->gbl->vres(x.ebdry(base.ebdry(0))->nseg,0) = 0.0;
			if (fix_norm) {
#ifndef petsc
				/* POST-REMOVE ADDED MASS FLUX TERM FOR FIXED POINT */
				x.gbl->res.v(base.pnt,x.NV-1) += surf->gbl->vres(x.ebdry(base.ebdry(0))->nseg,1)*x.gbl->rho;
#endif
				/* AND ZERO RESIDUAL */
				surf->gbl->vres(x.ebdry(base.ebdry(0))->nseg,1) = 0.0;
			}
		}
		else {
			/* SET TANGENT RESIDUAL TO ZERO */
			surf->gbl->vres(0,0) = 0.0;
			if (fix_norm) {
#ifndef petsc
				/* POST-REMOVE ADDED MASS FLUX TERM FOR FIXED POINT */
				x.gbl->res.v(base.pnt,x.NV-1) += surf->gbl->vres(0,1)*x.gbl->rho;
#endif
				/* AND ZERO RESIDUAL */
				surf->gbl->vres(0,1) = 0.0;
			}
		}
		return;
	}
	
	void vdirichlet() {
		if (surfbdry == 0) {
			for(int n=0;n<=fix_norm;++n)
				surf->gbl->vres(x.ebdry(base.ebdry(0))->nseg,n) = 0.0;
		}
		else {
			for(int n=0;n<=fix_norm;++n)
				surf->gbl->vres(0,n) = 0.0;
		}
	}
	
	void mvpttobdry(TinyVector<FLT,tri_mesh::ND> &pt) {
		if (surfbdry == 0) {
			x.ebdry(base.ebdry(1))->mvpttobdry(0,-1.0,pt);
		}
		else {
			x.ebdry(base.ebdry(0))->mvpttobdry(x.ebdry(base.ebdry(0))->nseg-1,1.0,pt);
		}
	}
	
#ifdef petsc
	void petsc_jacobian_dirichlet() {
		/* BOTH X & Y ARE FIXED */
		if (fix_norm) {
			Array<int,1> rows(tri_mesh::ND);
			for(int n=0;n<tri_mesh::ND;++n)
				rows(n) = (x.NV+tri_mesh::ND)*base.pnt +x.NV +n;
#ifdef MY_SPARSE
			x.J.zero_rows(tri_mesh::ND,rows);
			x.J_mpi.zero_rows(tri_mesh::ND,rows);
			x.J.set_diag(tri_mesh::ND,rows,1.0);
#else
			MatZeroRows(x.petsc_J,tri_mesh::ND,rows.data(),1.0);
#endif
		}
	}
#endif
};

class surface_outflow2 : public surface_fixed_pt2 {
	/* For a surface point sliding on a vertical or horizontal surface */
	/* For periodic wall have tri_mesh vertex type be comm */
protected:
	FLT position;
	enum {vertical, horizontal, angled, curved} wall_type;
	enum {prdc, free_angle, fixed_angle} contact_type;
	FLT contact_angle;  // ONLY USED FOR FIXED_ANGLE CONTACT TYPE
	TinyVector<FLT,tri_mesh::ND> wall_normal;
	
public:
	surface_outflow2(tri_hp_ins &xin, vrtx_bdry &bin) : surface_fixed_pt2(xin,bin) {mytype = "surface_outflow2";}
	surface_outflow2(const surface_outflow2& inbdry, tri_hp_ins &xin, vrtx_bdry &bin) : surface_fixed_pt2(inbdry,xin,bin), position(inbdry.position),
	wall_type(inbdry.wall_type), contact_type(inbdry.contact_type), contact_angle(inbdry.contact_angle), wall_normal(inbdry.wall_normal) {}
	surface_outflow2* create(tri_hp& xin, vrtx_bdry &bin) const {return new surface_outflow2(*this,dynamic_cast<tri_hp_ins&>(xin),bin);}
	
	void init(input_map& inmap,void* gbl_in) {
		surface_fixed_pt2::init(inmap,gbl_in);
		fix_norm = 0;
		
		std::string input;
		if (!inmap.get(base.idprefix + "_contact_type",input))
			contact_type = prdc;
		else if (input == "free_angle")
			contact_type = free_angle;
		else if (input == "fixed_angle") {
			contact_type = fixed_angle;
			inmap.getwdefault(base.idprefix +"_contact_angle",contact_angle,90.0);
			contact_angle *= M_PI/180.0;
		}
		else {
			*x.gbl->log << "Unrecognized contact type" << std::endl;
		}
		
		if (inmap.get(base.idprefix + "_wall_type",input)) {
			if (input == "vertical") {
				wall_type = vertical;
				position = x.pnts(base.pnt)(0);
			}
			else if (input == "horizontal") {
				wall_type = horizontal;
				position = x.pnts(base.pnt)(1);
			}
			else if (input == "angled")
				wall_type = angled;
			else if (input == "curved")
				wall_type = curved;
			else {
				*x.gbl->log << "Unrecognized wall type" << std::endl;
			}
		}
		else {
			wall_type = vertical;
			position = x.pnts(base.pnt)(0);
		}
		
		/* Find tangent to wall and use to constrain motion */
		int bnumwall = base.ebdry(1-surfbdry);
		TinyVector<FLT,tri_mesh::ND> rp;
		if (surfbdry == 0) {
			int sindwall = x.ebdry(bnumwall)->seg(0);
			x.crdtocht1d(sindwall);
			basis::tri(x.log2p)->ptprobe1d(2,&rp(0),&wall_normal(0),-1.0,&x.cht(0,0),MXTM);
		}
		else {
			int sindwall = x.ebdry(bnumwall)->seg(x.ebdry(bnumwall)->nseg-1);
			x.crdtocht1d(sindwall);
			basis::tri(x.log2p)->ptprobe1d(2,&rp(0),&wall_normal(0),1.0,&x.cht(0,0),MXTM);
		}
		FLT length = sqrt(wall_normal(0)*wall_normal(0) +wall_normal(1)*wall_normal(1));
		wall_normal /= length;
		FLT temp = wall_normal(0);
		wall_normal(0) = wall_normal(1);
		wall_normal(1) = -temp;
	}
	
	void mvpttobdry(TinyVector<FLT,tri_mesh::ND> &pt) {
		switch(wall_type) {
			case(vertical):
				pt(0) = position;
				break;
			case(horizontal):
				pt(1) = position;
				break;
			case(angled):case(curved):
				surface_fixed_pt2::mvpttobdry(pt);
				break;
		}
	}
	
	/* Routine to add surface tension stress */
	/* also zero's tangent residual in no petsc */
	void rsdl(int stage);
	
#ifdef petsc
	void petsc_jacobian();
	void petsc_jacobian_dirichlet() {}
#endif
};

#endif /* defined(__tri_hp__surface2__) */
