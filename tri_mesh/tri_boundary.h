#ifndef _tri_boundary_h_
#define _tri_boundary_h_

#include "tri_mesh.h"
#include "boundary.h"
#include <iostream>
#include <input_map.h>
#include <string>
#include <sstream>
#include <fstream>
#include "mappings.h"
#include "mapped_mesh.h"

/** \brief Specialization for a communiation vertex
 *
 * \ingroup boundary
 */
class vcomm : public comm_bdry<vrtx_bdry,tri_mesh> {
	public:
		vcomm(int inid, tri_mesh& xin) : comm_bdry<vrtx_bdry,tri_mesh>(inid,xin) {mytype="comm";}
		vcomm(const vcomm &inbdry, tri_mesh& xin) : comm_bdry<vrtx_bdry,tri_mesh>(inbdry,xin) {}

		vcomm* create(tri_mesh& xin) const {return new vcomm(*this,xin);}

		/** Generic routine to load buffers from array */
		void vloadbuff(boundary::groups group,FLT *base,int bgn,int end, int stride);
		/** Generic routine to receive into array */
		void vfinalrcv(boundary::groups group,int phase, comm_type type, operation op, FLT *base,int bgn,int end, int stride);
};

/** \brief Specialization for a communiation edge
 *
 * \ingroup boundary
 */
class ecomm : public comm_bdry<edge_bdry,tri_mesh> {
	public:
		/* CONSTRUCTOR */
		ecomm(int inid, tri_mesh& xin) : comm_bdry<edge_bdry,tri_mesh>(inid,xin) {mytype="comm";}
		ecomm(const ecomm &inbdry, tri_mesh& xin) : comm_bdry<edge_bdry,tri_mesh>(inbdry,xin) {}

		ecomm* create(tri_mesh& xin) const {return new ecomm(*this,xin);}

		/* GENERIC COMMUNICATIONS */
		void vloadbuff(boundary::groups group,FLT *base,int bgn,int end, int stride);
		void vfinalrcv(boundary::groups group,int phase, comm_type type, operation op, FLT *base,int bgn,int end, int stride);
		void sloadbuff(boundary::groups group,FLT *base,int bgn,int end, int stride);
		void sfinalrcv(boundary::groups group,int phase, comm_type type, operation op, FLT *base,int bgn,int end, int stride);
};

class emapped_comm : public ecomm {
    mapped_mesh& x;
    public:
        /* CONSTRUCTOR */
        emapped_comm(int inid, tri_mesh& xin) : ecomm(inid,xin), x(dynamic_cast<mapped_mesh& >(xin)) {mytype="emapped_comm";}
        emapped_comm(const emapped_comm &inbdry, tri_mesh& xin) : ecomm(inbdry,xin), x(inbdry.x) {}
        emapped_comm* create(tri_mesh& xin) const {return new emapped_comm(*this,xin);}
        void loadpositions() {vloadbuff(all,&(x.mapped_pnts(0)(0)),0,tri_mesh::ND-1,tri_mesh::ND);}
        void rcvpositions(int phase);
        void mvpttobdry(int nseg,FLT psi, TinyVector<FLT,tri_mesh::ND> &pt);  /* Move point to psi location in physical space */
};

class vmapped_comm : public vcomm {
    mapped_mesh& x;
    public:
        /* CONSTRUCTOR */
        vmapped_comm(int inid, tri_mesh& xin) : vcomm(inid,xin), x(dynamic_cast<mapped_mesh& >(xin)) {mytype="vmapped_comm";}
        vmapped_comm(const vmapped_comm &inbdry, tri_mesh& xin) : vcomm(inbdry,xin), x(inbdry.x) {}
        vmapped_comm* create(tri_mesh& xin) const {return new vmapped_comm(*this,xin);}
        void loadpositions() {vloadbuff(all,&(x.mapped_pnts(0)(0)),0,tri_mesh::ND-1,tri_mesh::ND);}
        void rcvpositions(int phase);
};



/** \brief Specialization for a parition edge
 *
 * \ingroup boundary
 * This is specifically for partitions that are not smooth boundaries
 * as typically created by metis.  Because they are not smooth it
 * changes what happens between different levels of multigrid. 
 * It also changes what happens during mesh adaptation 
 */
class epartition : public ecomm {
public:
	tri_mesh remote_halo;
	int npnt_h; /**< number of points in halo */
	int nseg_h; /**< number of interior segments in halo */
	int ntri_h; /**< number of triangles in halo */
	Array<int,1> pnt_h; /**< points in halo */
	Array<int,1> seg_h; /**< interior segments in halo */
	Array<int,1> sgn_h; /**< direction of segments relative to remote mesh */
	Array<int,1> tri_h; /**< triangles in halo */

public:
	/* CONSTRUCTOR */
	epartition(int inid, tri_mesh& xin) : ecomm(inid,xin), npnt_h(0), nseg_h(0), ntri_h(0) {add_to_group(boundary::partitions); mytype="partition";}
	epartition(const epartition &inbdry, tri_mesh& xin) : ecomm(inbdry,xin), npnt_h(0), nseg_h(0), ntri_h(0) {}
	epartition* create(tri_mesh& xin) const {return new epartition(*this,xin);}
	void alloc(int size);
	void copy(const edge_bdry& bin);
	void mgconnect(Array<tri_mesh::transfer,1> &cnnct, tri_mesh& tgt, int bnum);
	void mgconnect1(Array<tri_mesh::transfer,1> &cnnct, tri_mesh& tgt, int bnum);
	void calculate_halo();
	void receive_halo();
};


/** \brief Template to make a periodic vertex or edge
 *
 * \ingroup boundary
 * can do periodic in x or y (dir = 0/1)
 */
template<class BASE> class prdc_template : public BASE {
	protected:
		int dir;
	public:
		/* CONSTRUCTOR */
		prdc_template(int idin, tri_mesh &xin) : BASE(idin,xin), dir(0) {BASE::mytype="prdc";}
		prdc_template(const prdc_template<BASE> &inbdry, tri_mesh &xin) : BASE(inbdry,xin), dir(inbdry.dir) {}

		prdc_template<BASE>* create(tri_mesh& xin) const {return(new prdc_template<BASE>(*this,xin));}

		int& setdir() {return(dir);}
		void init(input_map& inmap) {
			std::string keyword;
			std::map<std::string,std::string>::const_iterator mi;
			std::istringstream data;

			BASE::init(inmap);

			keyword = BASE::idprefix + "_dir";
			mi = inmap.find(keyword);
			if (mi != inmap.end()) {
				data.str(mi->second);
				data >> dir;
				data.clear();
			}
		}

		/* SEND/RCV POINT POSITION */
		void loadpositions() { BASE::vloadbuff(BASE::all,&(BASE::x.pnts(0)(0)),1-dir,1-dir +tri_mesh::ND-2,tri_mesh::ND); }
		void rcvpositions(int phase) { BASE::vfinalrcv(BASE::all_phased,phase,BASE::master_slave,boundary::replace,&(BASE::x.pnts(0)(0)),1-dir,1-dir +tri_mesh::ND-2,tri_mesh::ND); }
};

/** \ingroup boundary */
//@{
typedef prdc_template<vcomm> vprdc;  /**< Periodic vertex point */
typedef prdc_template<ecomm> eprdc;  /**< Periodic edge boundary */
//@}

/** \brief Template to make a boundary and give it geometry
 *
 * \ingroup boundary
 * BASE is the boundary object
 * GEOM object should be of type "geometry"
 */
template<class BASE,class GEOM> class eboundary_with_geometry : public BASE, public rigid_movement_interface2D {
	public:
		GEOM geometry_object;

	public:
		eboundary_with_geometry(int inid, tri_mesh &xin) : BASE(inid,xin) {BASE::mytype=BASE::mytype+"analytic";}
		eboundary_with_geometry(const eboundary_with_geometry<BASE,GEOM> &inbdry, tri_mesh &xin) : BASE(inbdry,xin), rigid_movement_interface2D(inbdry), geometry_object(inbdry.geometry_object) {}
		eboundary_with_geometry* create(tri_mesh& xin) const {return(new eboundary_with_geometry<BASE,GEOM>(*this,xin));}

		void init(input_map& inmap) {
			BASE::init(inmap);
			rigid_movement_interface2D::init(inmap,BASE::idprefix);
			geometry_object.init(inmap,BASE::idprefix,*BASE::x.gbl->log);
		}

		void mvpttobdry(int nseg,FLT psi, TinyVector<FLT,tri_mesh::ND> &pt) {
			to_geometry_frame(pt);
			if (geometry_object.mvpttobdry(pt,BASE::x.gbl->time)) {
				*BASE::x.gbl->log << "trouble moving point to geometry for side " <<BASE::idprefix << std::endl;
				BASE::x.output("error");
				sim::abort(__LINE__,__FILE__,BASE::x.gbl->log);
			};
			to_physical_frame(pt);
			return;
		}
		
		void bdry_normal(int indx,FLT psi, TinyVector<FLT,tri_mesh::ND> &norm) {
			TinyVector<FLT,tri_mesh::ND> pt;
			for (int n=0;n<tri_mesh::ND;++n)
				pt(n) = 0.5*((1. -psi)*BASE::x.pnts(BASE::x.seg(BASE::seg(indx)).pnt(0))(n) +(1.+psi)*BASE::x.pnts(BASE::x.seg(BASE::seg(indx)).pnt(1))(n));
			to_geometry_frame(pt);
			geometry_object.bdry_normal(pt,BASE::x.gbl->time,norm);
			rotate_vector(norm);
			return;
		}
};

template<class BASE,class GEOM> class vboundary_with_geometry : public BASE {
	public:
		GEOM geometry_object;

	public:
		vboundary_with_geometry(int inid, tri_mesh &xin) : BASE(inid,xin) {BASE::mytype=BASE::mytype+"analytic";}
		vboundary_with_geometry(const vboundary_with_geometry<BASE,GEOM> &inbdry, tri_mesh &xin) : BASE(inbdry,xin), geometry_object(inbdry.geometry_object) {}
		vboundary_with_geometry* create(tri_mesh& xin) const {return(new vboundary_with_geometry<BASE,GEOM>(*this,xin));}

		void init(input_map& inmap) {
			BASE::init(inmap);
			geometry_object.init(inmap,BASE::idprefix,*BASE::x.gbl->log);
		}

		void mvpttobdry(TinyVector<FLT,tri_mesh::ND> &pt) {
			geometry_object.mvpttobdry(pt,BASE::x.gbl->time);
			return;
		}
};

/* Boundary objects that redirect to a physics class for obtaining geoemtry information */
class vcoupled_physics_ptr {
	public:
		vgeometry_interface<tri_mesh::ND> *physics;
};

class ecoupled_physics_ptr {
	public:
		egeometry_interface<tri_mesh::ND> *physics;
};

class fcoupled_physics_ptr {
	public:
		fgeometry_interface<tri_mesh::ND> *physics;
};

/* Class which uses geometry information for tstep < 0, then redirects to a physics class */
template<class BASE> class ecoupled_physics : public ecoupled_physics_ptr, public BASE {
  public:
		ecoupled_physics(int inid, tri_mesh &xin) : BASE(inid,xin) {BASE::mytype=BASE::mytype+"coupled";}
		ecoupled_physics(const ecoupled_physics<BASE> &inbdry, tri_mesh &xin) : BASE(inbdry,xin) {}
		ecoupled_physics* create(tri_mesh& xin) const {return(new ecoupled_physics<BASE>(*this,xin));}

		virtual void mvpttobdry(int nseg,FLT psi, TinyVector<FLT,tri_mesh::ND> &pt) {
			if (BASE::x.gbl->tstep < 0) BASE::mvpttobdry(nseg,psi,pt);
			else physics->mvpttobdry(nseg,psi,pt);
			return;
		}
};

#include <spline.h>

class spline_bdry : public edge_bdry, public rigid_movement_interface2D {
    public:
        spline<tri_mesh::ND> my_spline;
        Array<FLT,1> s;  // STORE S COORDINATE OF BOUNDARY POINTS (NOT WORKING)?
        FLT smin, smax; // LIMITS FOR BOUNDARY
        FLT scale, norm_dist;

		spline_bdry(int inid, tri_mesh &xin) : edge_bdry(inid,xin) {mytype="spline";}
		spline_bdry(const spline_bdry &inbdry, tri_mesh &xin) : edge_bdry(inbdry,xin), rigid_movement_interface2D(inbdry), my_spline(inbdry.my_spline), smin(inbdry.smin), smax(inbdry.smax), scale(inbdry.scale), norm_dist(inbdry.norm_dist) {}
		spline_bdry* create(tri_mesh& xin) const {return(new spline_bdry(*this,xin));}

		/* TEMPORARY INPUT/OUTPUTING/INIT NEEDS TO BE STRAIGHTENED OUT */
		void alloc(int n) {edge_bdry::alloc(n); s.resize(n+1);}
		void init(input_map& inmap) {
			edge_bdry::init(inmap);
			rigid_movement_interface2D::init(inmap,idprefix);
			std::string line;
			if (!inmap.get(edge_bdry::idprefix+"_filename",line)) {
				*x.gbl->log << "Couldn't fine spline file name in input file\n";
				sim::abort(__LINE__,__FILE__,x.gbl->log);
			}
			my_spline.read(line);

			inmap.getlinewdefault(edge_bdry::idprefix+"_s_limits",line,"0 1");
			std::istringstream data(line);
			data >> smin >> smax;
			data.clear();
            
            inmap.getwdefault(edge_bdry::idprefix+"_scale",scale,1.0);
            inmap.getwdefault(edge_bdry::idprefix+"_norm_dist",norm_dist,0.0);
		}

		void mvpttobdry(int seg_ind,FLT psi, TinyVector<FLT,tri_mesh::ND> &pt) {
			int sind = seg(seg_ind);
			int p0 = x.seg(sind).pnt(0);
			int p1 = x.seg(sind).pnt(1);
            TinyVector<FLT,tri_mesh::ND> pt0(x.pnts(p0));
            TinyVector<FLT,tri_mesh::ND> pt1(x.pnts(p1));
            
            to_geometry_frame(pt);
            to_geometry_frame(pt0);
            to_geometry_frame(pt1);
            pt /= scale;
            pt0 /= scale;
            pt1 /= scale;
            
			/* METHOD 1 */
			FLT sloc,sloc0,sloc1;
			my_spline.find(sloc0,pt0);
			my_spline.find(sloc1,pt1);
			/* FOR LOOPS */
            if (sloc0 >= smax) sloc0 = smin;
            if (sloc0 > sloc1) sloc1 = smax;
            
			sloc = 0.5*((1-psi)*sloc0 +(1+psi)*sloc1);
			my_spline.offset(sloc, norm_dist, pt);
            
//            *x.gbl->log << sloc0 << ' ' << pt0 << ' ' << sloc1 << ' ' << pt1 << ' ' << sloc << ' ' << pt << std::endl;

//            TinyVector<FLT,tri_mesh::ND> dx = pt1 -pt0;
//            FLT l2 = dx(0)*dx(0) +dx(1)*dx(1);
//            FLT ds,psinew;
//            int iter;
//            for (iter = 0; iter < 100; ++iter) {
//                psinew = 2*((pt(0)-pt0(0))*dx(0) +(pt(1)-pt0(1))*dx(1))/l2 -1.0;
//                ds = -(psinew-psi)*(sloc1-sloc0)/2.0;
//                sloc += 0.5*ds;
//                my_spline.interpolate(sloc,pt);
//                if (fabs(psinew-psi) < 1.0e-8) break;
//            }
//            if (iter > 99) {
//                *x.gbl->log << "too many spline iterations at " << pt << " with final change of " << ds << '\n';
//            }
            pt *= scale;
			to_physical_frame(pt);
			return;
		}
};

#endif


