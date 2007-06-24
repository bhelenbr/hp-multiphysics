#ifndef _mesh_h_
#define _mesh_h_

#include <math.h>
#include <quadtree.h>
#include <iostream>
#include <float.h>
#include <utilities.h>
#include <iostream>
#include <input_map.h>
#include <string>
#include <sstream>
#include <blitz/array.h>
#include "blocks.h"
#include "boundary.h"

#ifdef SINGLE
#define FLT float
#define EPSILON FLT_EPSILON
#else
#ifndef FLT
#define FLT double
#define EPSILON DBL_EPSILON
#endif
#endif

using namespace blitz;

class side_bdry;
class vrtx_bdry;

/** This is an unstructured triangular mesh class which has adaptation and 
parallel communication capabilities */

class mesh {

    /***************/
    /* DATA          */
    /***************/
    public: 
        std::string idprefix;
        int maxvst;
        static const int ND = 2;
                
        /* VERTEX DATA */
        int nvrtx;  /**< total number of vertices in the mesh */
        Array<TinyVector<FLT,ND>,1> vrtx; /**< physical location of the vertices in the mesh */
        Array<FLT,1> vlngth; /**< target mesh resolution in the neighborhood of each vertex */
        /** data structure for integer data at each vertex */
        struct vstruct {
            int tri; /**< pointer to a triangle connected to this vertex */
            int nnbor; /**< number of vertices sharing this vertex */
            int info; /**< General purpose (mostly for adaptation) */
        };
        Array<vstruct,1> vd;
        class quadtree<ND> qtree;
    
        /* VERTEX BOUNDARY INFO */
        int nvbd;
        Array<vrtx_bdry *,1> vbdry;
        vrtx_bdry* getnewvrtxobject(int idnum, input_map& bdrydata);

        /* SIDE DATA */        
        int nside;
        struct sstruct {
            TinyVector<int,2> vrtx;
            TinyVector<int,2> tri;
            int info;
        };
        Array<sstruct,1> sd;
        
        /* SIDE BOUNDARY INFO */
        int nsbd;
        Array<side_bdry *,1> sbdry;
        side_bdry* getnewsideobject(int idnum, input_map& bdrydata);

        /* TRIANGLE DATA */        
        int ntri;
        struct tstruct {
            TinyVector<int,3> vrtx;
            TinyVector<int,3> side;
            TinyVector<int,3> sign;
            TinyVector<int,3> tri;
            int info;
        };
        Array<tstruct,1> td;
        
        /* SOME WORK VARIABLES */
        /* ANY ROUTINE THAT USES i1wk SHOULD RESET IT TO -1 */
        /* THIS IS ONLY SHARED BETWEEN MESH OBJECTS */
        /* BECAUSE IT MUST BE KEPT AT -1 */
        static Array<int,1> i1wk;

        /* THESE ARE WORK VARIABLES. */
        /* THEY POINT TO GLOBALLY SHARED MEMORY AND THUS ARE */
        /* NOT GUARANTEED TO BE STABLE AFTER CALLS TO OTHER BLOCKS */
        static Array<FLT,1> fscr1;
        static Array<int,1> i2wk, i2wk_lst1, i2wk_lst2, i2wk_lst3;

        int initialized;
        static int maxsrch;
        
        /**************/
        /*  INTERFACE */
        /**************/
        /* INITIALIZATION & ALLOCATION */
        mesh() : idprefix(""), nvbd(0), nsbd(0), initialized(0)  {}
        void allocate(int mxsize);
        void allocate_duplicate(FLT sizereduce1d,const class mesh& xmesh);
        void reload_scratch_pointers();
        size_t needed_scratch_size();
        void copy(const mesh& tgt);
        virtual ~mesh();
        
        /* INPUT/OUTPUT MESH (MAY MODIFY VINFO/SINFO/TINFO) */
        enum filetype {easymesh, gambit, tecplot, grid, text, binary, BRep, mavriplis, boundary, vlength, debug_adapt, datatank};
        void input(const std::string &filename, filetype ftype,  FLT grwfac, input_map& bdrymap);
        int output(const std::string &filename, filetype ftype = grid) const;
        void bdry_output(const std::string &filename) const;
        virtual void setinfo();  // FOR EASYMESH OUTPUT (NOT USED)

        /* MESH MODIFICATION UTILTIES */
        void coarsen_substructured(const class mesh &tgt,int p);
	    void symmetrize();
        void append(const mesh &z);
        void shift(TinyVector<FLT,ND>& s);
        void scale(TinyVector<FLT,ND>& s);
        int smooth_cofa(int niter);
        void refineby2(const class mesh& xmesh);
        void settrim();
        void initvlngth();
        block::ctrl adapt(block::ctrl ctrl_message, FLT tolsize);
        int coarsen(FLT factor, const class mesh& xmesh);
        block::ctrl coarsen2(block::ctrl ctrl_message, FLT factor, const class mesh& inmesh, FLT size_reduce = 1.0);
        void coarsen3();

        /* UTILITIES FOR PARALLEL COMPUTATIONS */
#ifdef METIS
        void setpartition(int nparts); 
#endif 
        void partition(class mesh& xmesh, int npart);
        int comm_entity_size();
        int comm_entity_list(Array<int,1>& list);
        void vmsgload(boundary::groups group, int phase, boundary::comm_type type, FLT *base, int bgn, int end, int stride);
        void vmsgpass(boundary::groups group,int phase, boundary::comm_type type);
        int vmsgwait_rcv(boundary::groups group, int phase, boundary::comm_type type, boundary::operation op, FLT *base,int bgn, int end, int stride);
        int vmsgrcv(boundary::groups group, int phase, boundary::comm_type type, boundary::operation, FLT *base,int bgn, int end, int stride);
        void smsgload(boundary::groups group, int phase, boundary::comm_type type, FLT *base,int bgn, int end, int stride);
        void smsgpass(boundary::groups group, int phase, boundary::comm_type type);
        int smsgwait_rcv(boundary::groups group,int phase, boundary::comm_type type, boundary::operation op, FLT *base,int bgn, int end, int stride);
        int smsgrcv(boundary::groups group,int phase, boundary::comm_type type,  boundary::operation op, FLT *base,int bgn, int end, int stride);
        block::ctrl matchboundaries(block::ctrl ctrl_message);
        
        
        /* UTILITIES FOR INTERPOLATION BETWEEN MESHES */
        struct transfer {
            int tri;
            TinyVector<FLT,3> wt;
        };
        block::ctrl mgconnect(block::ctrl ctrl_message, Array<transfer,1> &cnnct, const class mesh& tgt);
        void testconnect(const std::string &fname,Array<transfer,1> &cnnct, mesh *cmesh);
        
        /* SOME DEGUGGING FUNCTIONS */
        void checkintegrity();
        void checkintwk() const;
        
    protected:
    
        /*******************/
        /* INTERNAL FUNCTIONS */
        /*******************/  
        void cnt_nbor(void);
        void bdrylabel(void);
        void createsideinfo(void);
        void createtdstri(void);
        void createttri(void);
        void createvtri(void);
        void treeinit();
        void treeinit(FLT x1[ND], FLT x2[ND]);

        /* MESH MODIFICATION FUNCTIONS */
        /* 4 binary digits for vertex, side, tri */
        const static int VSPEC = 0x4, VDLTE = 0x2, VTOUC=0x1;
        const static int SDLTE = 0x10*0x2, STOUC=0x10*0x1;
        const static int TSRCH = 0x100*0x4, TDLTE = 0x100*0x2, TTOUC=0x100*0x1;
        void setup_for_adapt();
        void swap(FLT swaptol = EPSILON);
        void bdry_yaber(FLT tolsize);
        void bdry_yaber1();
        void yaber(FLT tolsize);
        void bdry_rebay(FLT tolsize);
        void bdry_rebay1();
        void rebay(FLT tolsize);
        void cleanup_after_adapt();
     
        /* TO CREATE AN INITIAL TRIANGUlATION */
        void triangulate(int nside);
        void addtri(int v0,int v1,int v2,int sind,int dir);

        /* TO INSERT A POINT */
        int insert(const TinyVector<FLT,ND> &x);
        int insert(int vnum, int tnum);
        void bdry_insert(int vnum, int sind, int endpt = 0);
        int findtri(TinyVector<FLT,ND> x, int vnear) const;

        /* FOR COARSENING A SIDE */
        void collapse(int sind, int endpt);

        void dltvrtx(int vind);
        virtual void movevdata(int frm, int to) {}
        virtual void movevdata_bdry(int bnum,int bel,int endpt) {}
        virtual void updatevdata(int v) {}
        virtual void updatevdata_bdry(int bnum,int bel,int endpt) {}
        void dltsd(int sind);
        virtual void movesdata(int frm, int to) {}
        virtual void movesdata_bdry(int bnum,int bel) {}
        virtual void updatesdata(int s) {}
        virtual void updatesdata_bdry(int bnum,int bel) {}
        void dlttri(int tind);
        virtual void movetdata(int frm, int to) {}
        virtual void updatetdata(int t) {}
        
        /* ORDERED LIST FUNCTIONS */
        void putinlst(int sind);
        void tkoutlst(int sind);
        
        /* TO EDGE SWAP FOR IMPROVED QUALITY */
        int swap(int sind, FLT tol = 0.0);  //SWAPS A SINGLE SIDE 
        void swap(int nswp, int *swp, FLT tol = 0.0); //SWAPS A LIST OF SIDES

        /* PRIMITIVE FUNCTIONS */
    public:
        inline FLT distance(int v0, int v1) const {
            FLT d = 0.0;
            for(int n = 0; n<ND;++n)
                d += pow(vrtx(v0)(n) -vrtx(v1)(n),2);
            return(sqrt(d));
        }
        inline FLT distance2(int v0, int v1) const {
            FLT d = 0.0;
            for(int n = 0; n<ND;++n)
                d += pow(vrtx(v0)(n) -vrtx(v1)(n),2);
            return(d);
        }  
        FLT incircle(int tind, const TinyVector<FLT,ND> &x) const;
        FLT insidecircle(int sind, const TinyVector<FLT,ND> &x) const;
        FLT area(int sind, int vind) const;
        FLT area(int v0, int v1, int v2) const;
        FLT area(int tind) const;
        FLT minangle(int v0, int v1, int v2) const;
        FLT angle(int v0, int v1, int v2) const;
        FLT circumradius(int tind) const;
        void circumcenter(int tind, TinyVector<FLT,2> &x) const;
        FLT inscribedradius(int tind) const;
        FLT aspect(int tind) const;
        FLT intri(int tind, const TinyVector<FLT,2> &x) const;
        void getwgts(TinyVector<FLT,3> &wgt) const;
        /* tri numbers at boundary point to el and group */
        int getbdrynum(int trinum) const { return((-trinum>>16) -1);}
        int getbdryel(int trinum) const { return(-trinum&0xFFFF);}
        int trinumatbdry(int bnum, int bel) const { return(-(((bnum+1)<<16) +bel));}
    private:
        int excpt, excpt1, mp_phase;
};

class vgeometry_interface {
    public:
        virtual void mvpttobdry(TinyVector<FLT,mesh::ND> &pt) {}
        virtual ~vgeometry_interface() {}
};

class vgeometry_pointer {
    public:
        vgeometry_interface *solution_data;
};

class sgeometry_interface {
    public:
        virtual void mvpttobdry(int nel, FLT psi, TinyVector<FLT,mesh::ND> &pt) {}
        virtual ~sgeometry_interface() {}
};

class sgeometry_pointer {
    public:
        sgeometry_interface *solution_data;
};


/** \brief Specialization for a vertex 
 *
 * \ingroup boundary
 * Specialization of communication routines and 
 * and storage for a vertex boundary 
 */
class vrtx_bdry : public boundary, public vgeometry_interface {
    public:
        mesh &x;
        TinyVector<int,2> sbdry;
        int v0;
        
        /* CONSTRUCTOR */
        vrtx_bdry(int intype, mesh &xin) : boundary(intype), x(xin) {idprefix = x.idprefix +"_v" +idprefix; mytype="plain";}
        vrtx_bdry(const vrtx_bdry &inbdry, mesh &xin) : boundary(inbdry.idnum), x(xin)  {idprefix = inbdry.idprefix; mytype = inbdry.mytype; sbdry = inbdry.sbdry;}

        /* OTHER USEFUL STUFF */
        virtual vrtx_bdry* create(mesh &xin) const { return(new vrtx_bdry(*this,xin));}
        virtual void copy(const vrtx_bdry& bin) {
            v0 = bin.v0;
            sbdry = bin.sbdry;
        }
        virtual void vloadbuff(boundary::groups group, FLT *base, int bgn, int end, int stride) {}
        virtual void vfinalrcv(boundary::groups group, int phase, comm_type type, operation op, FLT *base, int bgn, int end, int stride) {}

        virtual void mvpttobdry(TinyVector<FLT,mesh::ND> &pt) {}
        virtual void loadpositions() {vloadbuff(all,&(x.vrtx(0)(0)),0,mesh::ND-1,mesh::ND);}
        virtual void rcvpositions(int phase) {vfinalrcv(all_phased,phase,master_slave,replace,&(x.vrtx(0)(0)),0,mesh::ND-1,mesh::ND);}
};


/** \brief Specialization for a side 
 *
 * \ingroup boundary
 * Specialization of communication routines and 
 * and storage for a side boundary 
 */
class side_bdry : public boundary, public sgeometry_interface {
    public:
        mesh &x;
        TinyVector<int,2> vbdry;
        int maxel;
        int nel;
        Array<int,1> el;
        
        /* CONSTRUCTOR */
        side_bdry(int inid, mesh &xin) : boundary(inid), x(xin), maxel(0)  {idprefix = x.idprefix +"_s" +idprefix; mytype="plain"; vbdry = -1;}
        side_bdry(const side_bdry &inbdry, mesh &xin) : boundary(inbdry.idnum), x(xin), maxel(0)  {idprefix = inbdry.idprefix; mytype = inbdry.mytype; vbdry = inbdry.vbdry;}
        
        /* BASIC B.C. STUFF */
        void alloc(int n);
        virtual side_bdry* create(mesh &xin) const {
            return(new side_bdry(*this,xin));
        }
        virtual void copy(const side_bdry& bin);
        
        /* ADDITIONAL STUFF FOR SIDES */
        virtual void swap(int s1, int s2);
        virtual void reorder();
        virtual block::ctrl mgconnect(block::ctrl ctrl_message, Array<mesh::transfer,1> &cnnct, const class mesh& tgt, int bnum);
        virtual void mvpttobdry(int nel, FLT psi, TinyVector<FLT,mesh::ND> &pt);
        virtual void findbdrypt(const TinyVector<FLT,2> xpt, int &sidloc, FLT &psiloc) const;
        
        /* DEFAULT SENDING FOR SIDE VERTICES */
        virtual void vloadbuff(boundary::groups group,FLT *base,int bgn,int end, int stride) {}
        virtual void vfinalrcv(boundary::groups group,int phase, comm_type type, operation op, FLT *base,int bgn,int end, int stride) {}
        virtual void sloadbuff(boundary::groups group,FLT *base,int bgn,int end, int stride) {}
        virtual void sfinalrcv(boundary::groups group,int phase, comm_type type, operation op, FLT *base,int bgn,int end, int stride) {}
        virtual void loadpositions() {vloadbuff(all,&(x.vrtx(0)(0)),0,mesh::ND-1,mesh::ND);}
        virtual void rcvpositions(int phase) {vfinalrcv(all_phased,phase,master_slave,replace,&(x.vrtx(0)(0)),0,mesh::ND-1,mesh::ND);}
};

#endif
