#include "tri_mesh.h"

#include <iostream>
#include <input_map.h>
#include <string>
#include <sstream>
#include <fstream>

//#define MPDEBUG

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

/** \brief Specialization for a parition edge 
 *
 * \ingroup boundary
 * This is specifically for partitions that are not smooth boundaries
 * as typically created by metis.  Because they are not smooth it 
 * changes what happens between different levels of multigrid.
 */
class epartition : public ecomm {
    public:
        /* CONSTRUCTOR */
        epartition(int inid, tri_mesh& xin) : ecomm(inid,xin) {groupmask = 1;mytype="partition";}
        epartition(const epartition &inbdry, tri_mesh& xin) : ecomm(inbdry,xin) {}

        epartition* create(tri_mesh& xin) const {return new epartition(*this,xin);}
        void mgconnect(Array<tri_mesh::transfer,1> &cnnct, tri_mesh& tgt, int bnum);
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
        void output(std::ostream& fout) {
            BASE::output(fout);
            fout << BASE::idprefix << "_dir" << ": " << dir << std::endl;  
        }
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
template<class BASE,class GEOM> class eboundary_with_geometry : public BASE {
    public:
        GEOM geometry_object;
        
    public: 
        eboundary_with_geometry(int inid, tri_mesh &xin) : BASE(inid,xin) {BASE::mytype=BASE::mytype+"analytic";}
        eboundary_with_geometry(const eboundary_with_geometry<BASE,GEOM> &inbdry, tri_mesh &xin) : BASE(inbdry,xin), geometry_object(inbdry.geometry_object) {}
        eboundary_with_geometry* create(tri_mesh& xin) const {return(new eboundary_with_geometry<BASE,GEOM>(*this,xin));}

        void output(std::ostream& fout) {
            BASE::output(fout);
            geometry_object.output(fout,BASE::idprefix);
        }
        void init(input_map& inmap) {
            BASE::init(inmap);
            geometry_object.init(inmap,BASE::idprefix,*BASE::x.gbl->log);
        }
        
        void mvpttobdry(int nseg,FLT psi, TinyVector<FLT,tri_mesh::ND> &pt) {
            geometry_object.mvpttobdry(pt,BASE::x.gbl->time);
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

        void output(std::ostream& fout) {
            BASE::output(fout);
            geometry_object.output(fout,BASE::idprefix);
        }
        void init(input_map& inmap) {
            BASE::init(inmap);
            geometry_object.init(inmap,BASE::idprefix,*BASE::x.gbl->log);
        }
        
        void mvpttobdry(TinyVector<FLT,tri_mesh::ND> &pt) {
            geometry_object.mvpttobdry(pt,BASE::x.gbl->time);
            return;
        }
};


/** \brief Interface & template to make a boundary that can be coupled to some other geometry object after tstep = 0
 *
 * \ingroup boundary
 * geometry in class is for initial condition, then
 * Physics object must provide geometry
 * BASE is the boundary object 
 */
class ecoupled_physics_interface : public geometry<tri_mesh::ND> {
    public:
        egeometry_interface<2> *physics;
};

template<class BASE> class ecoupled_physics : public ecoupled_physics_interface, public BASE {
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
#include <blitz/tinyvec-et.h>

class spline_bdry : public edge_bdry {
    spline<2> my_spline;
    Array<FLT,1> s;  // STORE S COORDINATE OF BOUNDARY POINTS (NOT WORKING)?
    FLT smin, smax; // LIMITS FOR BOUNDARY
    
    public: 
        spline_bdry(int inid, tri_mesh &xin) : edge_bdry(inid,xin) {mytype="spline";}
        spline_bdry(const spline_bdry &inbdry, tri_mesh &xin) : edge_bdry(inbdry,xin), my_spline(inbdry.my_spline), smin(inbdry.smin), smax(inbdry.smax) {}
        spline_bdry* create(tri_mesh& xin) const {return(new spline_bdry(*this,xin));}
        
        /* TEMPORARY INPUT/OUTPUTING/INIT NEEDS TO BE STRAIGHTENED OUT */
        void alloc(int n) {edge_bdry::alloc(n); s.resize(n+1);}
        void output(std::ostream& fout) {
            edge_bdry::output(fout);
            fout << s(Range(0,edge_bdry::nseg-1)) << std::endl;
        }
        void init(input_map& inmap) {
            edge_bdry::init(inmap);
			
			std::string line;
			if (!inmap.get(edge_bdry::idprefix+"_filename",line)) {
				*x.gbl->log << "Couldn't fine spline file name in input file\n";
				exit(1);
			}
            my_spline.read(line);
			
            inmap.getlinewdefault(edge_bdry::idprefix+"_s_limits",line,"0 1");
            std::istringstream data(line);
            data >> smin >> smax;
            data.clear();     
        }
        
        void input(std::istream& fin, tri_mesh::filetype type) {
            edge_bdry::input(fin,type);
			for(int i=0;i<nseg+1;++i) {
				fin >> s(i);
			}
        }
        
        void mvpttobdry(int nseg,FLT psi, TinyVector<FLT,tri_mesh::ND> &pt) {
            /* TEMPORARY THIS IS A HACK UNTIL I GET PARAMETRIC BOUNDARIES WORKING BETTER */
            int sind = seg(nseg);
            int p0 = x.seg(sind).pnt(0);
            int p1 = x.seg(sind).pnt(1);
            float sloc;
						
            /* METHOD 1 */
            FLT sloc0,sloc1;
            my_spline.find(sloc0,x.pnts(p0));
            /* FOR LOOPS */
            if (sloc0 > smax) sloc0 = smin;
            
            my_spline.find(sloc1,x.pnts(p1));
            /* FOR LOOPS */
            if (sloc1 < smin) sloc1 = smax;
       
            sloc = 0.5*((1-psi)*sloc0 +(1+psi)*sloc1);
            my_spline.interpolate(sloc,pt);
			            
            TinyVector<FLT,tri_mesh::ND> dx = x.pnts(p1) -x.pnts(p0);
            FLT l2 = dx(0)*dx(0) +dx(1)*dx(1);
            FLT ds,psinew;
            int iter;
            for (iter = 0; iter < 100; ++iter) {
                psinew = 2*((pt(0)-x.pnts(p0)(0))*dx(0) +(pt(1)-x.pnts(p0)(1))*dx(1))/l2 -1.0;
                ds = -(psinew-psi)*(sloc1-sloc0)/2.0;
                sloc += 0.5*ds;
                my_spline.interpolate(sloc,pt);
                if (fabs(psinew-psi) < 1.0e-8) break;
            }
            if (iter > 99) {
                *x.gbl->log << "too many spline iterations\n";
            }

            return;
        }
};



