#ifndef _tri_mesh_h_
#define _tri_mesh_h_

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

class vrtx_bdry;
class edge_bdry;

/** This is an unstructured triangular mesh class which has adaptation and 
parallel communication capabilities */

class tri_mesh : public multigrid_interface {

    /***************/
    /* DATA        */
    /***************/
    public: 
        int maxpst;  /**< maximum length of point, segment, or tri lists (all same for convenience) */
        static const int ND = 2;  /**< 2 eimensional */
                
        /** @name Variables describing mesh points */
        //@{
        int npnt;  /**< Total number of points in the mesh */
        Array<TinyVector<FLT,ND>,1> pnts; /**< Physical location of the points in the mesh */
        Array<FLT,1> lngth; /**< Target mesh resolution in the neighborhood of each point */
        /** Data structure for integer data at each point */
        struct pntstruct {
            int tri; /**< Pointer to a triangle connected to this point */
            int nnbor; /**< Number of segments connected to this point */
            int info; /**< General purpose (mostly for adaptation) */
        };
        Array<pntstruct,1> pnt;  /**< Array of point integer data */
        class quadtree<ND> qtree; /**< Quadtree for finding points */
        //@}
   
        /** @name Variables describing vertex (0D) boundary conditions */
        //@{
        int nvbd;  /**< number of vertex boundaries */
        Array<vrtx_bdry *,1> vbdry; /**< Array of vertex boundary objects */
        vrtx_bdry* getnewvrtxobject(int idnum, input_map& bdrydata);  /**< function for obtaining different vertex boundary objects */
        //@}

        /** @name Variables describing mesh segments */
        //@{    
        int nseg; /**< number of segments in mesh */
        /** Data structure for each segment */
        struct segstruct {
            TinyVector<int,2> pnt; /**< enpoints */
            TinyVector<int,2> tri; /**< adjacent triangles */
            int info; /**< General purpose (mostly for adaptation) */
        };
        Array<segstruct,1> seg; /**< Array of segment data */
        //@}

        /** @name Variables describing edge (1D) boundary conditions */
        //@{
        int nebd; /**< number of edge boundaries */
        Array<edge_bdry *,1> ebdry; /**< array of edge boundary objects */
        edge_bdry* getnewedgeobject(int idnum, input_map& bdrydata); /**< function for obtaining different edge boundary objects */
        //@}

        /** @name Variables describing triangles */
        //@{          
        int ntri; /**< Number of triangles in mesh */
        /** data structure for each triangle */
        struct tristruct {
            TinyVector<int,3> pnt; /**< triangle points */
            TinyVector<int,3> seg; /**< triangle segments */
            TinyVector<int,3> sgn; /**< sign is positive if counter-clockwise is the same direction as segment definition */
            TinyVector<int,3> tri; /**< adjacent triangles */
            int info; /**< General purpose (mostly for adaptation) */
        };
        Array<tristruct,1> tri;
        //@}

        /** /struct For information shared between meshes not used simultaneously (multigrid levels) */
        struct global : public block_global {            
            Array<int,1> intwk; /**< Integer work array, any routine that uses intwk should reset it to -1 */
            Array<FLT,1> fltwk; /**< Floating point work array */
            int nlst; /**< variable to keep track of number of entities in list */
            Array<int,1> i2wk, i2wk_lst1, i2wk_lst2, i2wk_lst3; /**< Arrays for storing lists */
            int maxsrch; /**< integer describing maximum number of triangles to search before giving up */
        } *gbl;
        
        bool initialized; /**< Initialization flag */
        
        /**************/
        /*  INTERFACE */
        /**************/
        /** Constructor */
        tri_mesh() : nvbd(0), nebd(0), gbl(0), initialized(0)  {}
        /** Routine to allocate shared variables */
        void* create_global_structure() {return new global;}
        /** Routine to initialize with using information in map and shared resource in gbl_in */
        void init(input_map& input, void *gbl_in);
        /** Routine to initialze from another mesh with option of increasing or decreasing storage (compatible with block.h) */
        void init(const multigrid_interface& mgin, FLT sizereduce1d = 1.0);
        /** Routine to copy */
        void copy(const tri_mesh& tgt);

        /** Innput/Output file types */
        enum filetype {easymesh, gambit, tecplot, grid, text, binary, BRep, mavriplis, boundary, vlength, debug_adapt, datatank};
        
        /** Input mesh */
        void input(const std::string &filename, filetype ftype,  FLT grwfac, input_map &input);
        /** Virtual routine so inheritors can set up info after input/adaptation */
        virtual void setinfo();  
        
        /** Outputs solution in various filetypes */
        void output(const std::string &filename, filetype ftype = grid) const;

         /** @name adapt adaptation routines */
        //@{  
        void initlngth();  /**< Set target mesh resolution based on current mesh */
        virtual void length() {} /**< Virtual so inheritors can set target resolution before adaptation */
        void adapt(); /**< Adapt mesh */
        //@}

        /* Destructor (deletes boundary objects) */
        virtual ~tri_mesh();
        
         /** @name Mesh modification utilties */
        //@{   
	    void symmetrize();  /**< Creates a symmetric mesh about y = 0 */
        void cut(); /**< Cut's mesh based on indicator array inf fltwk (Positive / Negative at each point) and aligns cut edge along 0 */
        void trim(); /**< Starting from boundaries deletes triangles based on sign of fltwk */
        void append(const tri_mesh &z); /**< Appends mesh */
        void shift(TinyVector<FLT,ND>& s); /**< Tranlates mesh */
        void scale(TinyVector<FLT,ND>& s); /**< Scales mesh */
        int smooth_cofa(int niter); /**< Does a center of area smoothing */
        void refineby2(const class tri_mesh& xmesh); /**< Refines by 2 */
        void coarsen_substructured(const class tri_mesh &tgt,int p); /**< Coarsens mesh that was output by hp FEM method */ 
        //@}

         /** @name Routines for parallel computations */
        //@{ 
#ifdef METIS
        void setpartition(int nparts);  /**< Set partition of mesh (in tri(i).info) */
#endif 
        void partition(class tri_mesh& xmesh, int npart); /**< Creates a partition from xmesh */
        int comm_entity_size(); /**< Returns size of list of communication entities (for blocks.h) */
        int comm_entity_list(Array<int,1>& list);  /**< Returns list of communication entities */
        class boundary* getvbdry(int num); /**< Returns pointer to vertex boundary (for blocks.h) */
        class boundary* getebdry(int num); /**< Returns pointer to edge boundary (for blocks.h) */
        void pmsgload(boundary::groups group, int phase, boundary::comm_type type, FLT *base, int bgn, int end, int stride); /**< Loads message buffers by point index from base (vertex and edge boundaries) */
        void pmsgpass(boundary::groups group,int phase, boundary::comm_type type); /**< Sends vertex and edge boundary messages */
        int pmsgwait_rcv(boundary::groups group, int phase, boundary::comm_type type, boundary::operation op, FLT *base,int bgn, int end, int stride); /**< Receives vertex and edge boundary messages */
        int pmsgrcv(boundary::groups group, int phase, boundary::comm_type type, boundary::operation, FLT *base,int bgn, int end, int stride);/**< Receives without waiting (Don't use) */
        void smsgload(boundary::groups group, int phase, boundary::comm_type type, FLT *base,int bgn, int end, int stride);/**< Loads message buffers by segment index from base (edge boundaries) */
        void smsgpass(boundary::groups group, int phase, boundary::comm_type type);/**< Sends edge boundary messages */
        int smsgwait_rcv(boundary::groups group,int phase, boundary::comm_type type, boundary::operation op, FLT *base,int bgn, int end, int stride); /**< Receives edge boundary messages */
        int smsgrcv(boundary::groups group,int phase, boundary::comm_type type,  boundary::operation op, FLT *base,int bgn, int end, int stride); /**< Receives without waiting (Don't use) */
        void matchboundaries();  /**< Matches location of boundary points */
        //@}
        
         /** @name Routines for multigrid */
        //@{ 
        int coarse_level; /**< Integer telling which multigrid level I am */
        multigrid_interface *fine, *coarse; /**< Pointers to fine and coarse objects */
        void coarsen(FLT factor, const class tri_mesh& xmesh); /**< Coarsens by triangulating then inserting */
        void coarsen2(FLT factor, const class tri_mesh& inmesh, FLT size_reduce = 1.0);  /**< Coarsens by using yaber */
        void coarsen3(); /**< Coarsens based on marks stored in pnt().info */

        /** Structure for mesh connection data */
        struct transfer {
            int tri; /**< Tri containing this point */
            TinyVector<FLT,3> wt; /**< Interpolation weights */
        };
        Array<transfer,1> fcnnct, ccnnct; /**< Arrays for connection data for each point */
        void connect(multigrid_interface& tgt);  /**< Set-up connections  to fine mesh compatible with block.h */
        void mgconnect(tri_mesh &tgt, Array<transfer,1> &cnnct); /**< Utility called by connect */
        void testconnect(const std::string &fname,Array<transfer,1> &cnnct, tri_mesh *cmesh); /**< Tests by outputing grids with interpolated point locations */
        //@}
        
        /* SOME DEGUGGING FUNCTIONS */
        void checkintegrity(); /**< Checks data arrays for compatibility */
        void checkintwk() const; /**< Makes sure intwk is all -1 */
        
    protected:
    
        /** @name Setup routines */
        //@{
        void allocate(int mxsize);  /**< Allocates memory */
        void cnt_nbor(void); /**< Fills in pnt().nnbor */
        void bdrylabel(void); /**< Makes seg().tri and tri().tri on boundary have pointer to boundary group/element */
        void createseg(void); /**< Creates all segment information from list of triangle points and also tri().seg/sgn (if necessary) */
        void createsegtri(void); /**< Creates seg.tri() and tri().seg connections */
        void createtritri(void); /**< Creates tri().tri data */
        void createpnttri(void); /**< Creates pnt().tri data */
        void treeinit(); /**< Initialize quadtree data */
        void treeinit(FLT x1[ND], FLT x2[ND]); /** Initialize quadtree data with specified domain */
        //@}
        
        /** @name Mesh modification functions */
        //@{
         /** 4 binary digits  used for each pnt, seg, tri adapt data 
         *  Typically special,deleted,touched,searched 
         */
        const static int PSPEC = 0x4, PDLTE = 0x2, PTOUC=0x1;
        const static int SDLTE = 0x10*0x2, STOUC=0x10*0x1;
        const static int TSRCH = 0x100*0x4, TDLTE = 0x100*0x2, TTOUC=0x100*0x1;
        void setup_for_adapt(); /**< Set all flags */

        void triangulate(int nseg); /**< Creates an initial triangulation */
        void addtri(int p0,int p1,int p2,int sind,int dir); /**< Utility for creating triangles used by triangulate */
        
        void swap(FLT swaptol = EPSILON); /**< Edge swap */
        int swap(int sind, FLT tol = 0.0);  /**< Swaps a single segment */
        void swap(int nswp, int *swp, FLT tol = 0.0); /**< Swaps a list of segments */
        
        void bdry_yaber(FLT tolsize); /**< Coarsen edge boundaries */
        void bdry_yaber1(); /**< Coarsen slave edge boundaries */
        void yaber(FLT tolsize); /**< Coarsen mesh (Rebay backwards) */
        void collapse(int sind, int endpt); /**< Removes by collapsing segment sind to endpt */

        void bdry_rebay(FLT tolsize); /**< Refine edges */
        void bdry_rebay1(); /**< Refine slave edges */
        void rebay(FLT tolsize); /**< Refine mesh using Rebay point placement */
        int insert(const TinyVector<FLT,ND> &x);  /**< Inserts a point */
        int insert(int pnum, int tnum);  /**< Inserts point at pnum into tnum */
        void bdry_insert(int pnum, int sind, int endpt = 0); /**< Inserts a boundary point in segment sind if endpt is 0 makes left old and right new seg */
        int findtri(TinyVector<FLT,ND> x, int pnear); /**< Locate triangle containing point with initial seed of pnear */
        
        void cleanup_after_adapt(); /**< Clean up and move data etc.. */
        void dltpnt(int pind); /**< Removes leftover point references */
        virtual void movepdata(int frm, int to) {} /**< Virtual routine so inheritors can automatically have their point data moved */
        virtual void movepdata_bdry(int bnum,int bel,int endpt) {} /**< Virtual routine so inheritors can automatically have their boundary point data moved */ 
        virtual void updatepdata(int v) {} /**< Virtual routine so inheritors can automatically have new point data updated */
        virtual void updatepdata_bdry(int bnum,int bel,int endpt) {} /**< Virtual routine so inheritors can automatically have new bondary point data updated */
        void dltseg(int sind); /**< Removes leftover segment references */
        virtual void movesdata(int frm, int to) {} /**< Virtual routine so inheritors can automatically have their segment data moved */
        virtual void movesdata_bdry(int bnum,int bel) {} /**< Virtual routine so inheritors can automatically have their boundry segment data moved  */
        virtual void updatesdata(int s) {} /**< Virtual routine so inheritors can update segment data */
        virtual void updatesdata_bdry(int bnum,int bel) {} /**< Virtual routine so inheritors can update boundary segment data */
        void dlttri(int tind); /**< Removes leftover triangle references */
        virtual void movetdata(int frm, int to) {} /**< Virtual routine so inheritors can automatically have their triangle data moved */
        virtual void updatetdata(int t) {}  /**< Virtual routine so inheritors can automatically have new triangles updated */
        
        void putinlst(int sind);  /**< For inserting into adaptation priority queues */
        void tkoutlst(int sind);  /**< For removing from adaption priority queues */ 
        //@}
        
        /** @name Some primitive functions */
        //@{
        /** Calculate distance between two points */
        inline FLT distance(int p0, int p1) const {
            FLT d = 0.0;
            for(int n = 0; n<ND;++n)
                d += pow(pnts(p0)(n) -pnts(p1)(n),2);
            return(sqrt(d));
        }
        /** Calculate the distance squared between two points */
        inline FLT distance2(int p0, int p1) const {
            FLT d = 0.0;
            for(int n = 0; n<ND;++n)
                d += pow(pnts(p0)(n) -pnts(p1)(n),2);
            return(d);
        }  
        FLT incircle(int tind, const TinyVector<FLT,ND> &x) const; /**< Calculate whether point is in circumscribed circle on tind*/
        FLT insegcircle(int sind, const TinyVector<FLT,ND> &x) const; /**< Calculate whether point is in in circumscribed circle of sind */
        FLT area(int sind, int pind) const;  /**< Calculates area of triangle formed by sind and pind */
        FLT area(int p0, int p1, int p2) const; /**< Calculates area of triangle formed by p0,p1,p2 */
        FLT area(int tind) const; /**< Calculates area of triangle */
        FLT minangle(int p0, int p1, int p2) const;  /**< Calculates minimum angle associated with triangle */
        FLT angle(int p0, int p1, int p2) const;  /**< Calculates angle point 1 is the vertex */
        FLT circumradius(int tind) const; /**< Calcultes circumbscribed radius */
        void circumcenter(int tind, TinyVector<FLT,2> &x) const;  /**< Calculates center of circumscribed circle */
        FLT inscribedradius(int tind) const; /**< Calculates inscribed radius */
        FLT aspect(int tind) const;  /**< Calculates aspect ratio basied on perimeter to area ratio normalized so 1 is equilateral */
        FLT intri(int tind, const TinyVector<FLT,2> &x);  /**< Determine whether a triangle contains a point (0.0 or negative) */
        TinyVector<FLT,3> tri_wgt; /**< Used in intri for searching (normalized value returned by getwgts after successful search) */
        void getwgts(TinyVector<FLT,3> &wgt) const; /**< Returns weighting for point interpolation within last triangle find by intri */

        int getbdrynum(int trinum) const { return((-trinum>>16) -1);}  /**< Uses info in seg.tri or tri.tri to determine boundary object number */
        int getbdryel(int trinum) const { return(-trinum&0xFFFF);}  /**< Uses info in seg.tri or tri.tri to determine boundary element */
        int trinumatbdry(int bnum, int bel) const { return(-(((bnum+1)<<16) +bel));} /**< Combines bnum & bel into 1 integer for storage in boundary of seg.tri or tri.tri */
        //@}

};


/** \brief Specialization for a vertex 
 *
 * \ingroup boundary
 * Specialization of communication routines and 
 * and storage for a vertex boundary 
 */
class vrtx_bdry : public boundary, vgeometry_interface<tri_mesh::ND> {
    public:
        tri_mesh &x;
        TinyVector<int,2> ebdry;
        int p0;
        
        /* CONSTRUCTOR */
        vrtx_bdry(int intype, tri_mesh &xin) : boundary(intype), x(xin) {idprefix = x.gbl->idprefix +"_v" +idprefix; mytype="plain";}
        vrtx_bdry(const vrtx_bdry &inbdry, tri_mesh &xin) : boundary(inbdry.idnum), x(xin)  {idprefix = inbdry.idprefix; mytype = inbdry.mytype; ebdry = inbdry.ebdry;}

        /* OTHER USEFUL STUFF */
        virtual vrtx_bdry* create(tri_mesh &xin) const { return(new vrtx_bdry(*this,xin));}
        virtual void copy(const vrtx_bdry& bin) {
            p0 = bin.p0;
            ebdry = bin.ebdry;
        }
        virtual void vloadbuff(boundary::groups group, FLT *base, int bgn, int end, int stride) {}
        virtual void vfinalrcv(boundary::groups group, int phase, comm_type type, operation op, FLT *base, int bgn, int end, int stride) {}

        virtual void mvpttobdry(TinyVector<FLT,tri_mesh::ND> &pt) {}
        virtual void loadpositions() {vloadbuff(all,&(x.pnts(0)(0)),0,tri_mesh::ND-1,tri_mesh::ND);}
        virtual void rcvpositions(int phase) {vfinalrcv(all_phased,phase,master_slave,replace,&(x.pnts(0)(0)),0,tri_mesh::ND-1,tri_mesh::ND);}
};


/** \brief Specialization for an edge
 *
 * \ingroup boundary
 * Specialization of communication routines and 
 * and storage for an edge boundary 
 */
class edge_bdry : public boundary, egeometry_interface<tri_mesh::ND> {
    public:
        tri_mesh &x;
        TinyVector<int,2> vbdry;
        int maxel;
        int nel;
        Array<int,1> el;
        
        /* CONSTRUCTOR */
        edge_bdry(int inid, tri_mesh &xin) : boundary(inid), x(xin), maxel(0)  {idprefix = x.gbl->idprefix +"_s" +idprefix; mytype="plain"; vbdry = -1;}
        edge_bdry(const edge_bdry &inbdry, tri_mesh &xin) : boundary(inbdry.idnum), x(xin), maxel(0)  {idprefix = inbdry.idprefix; mytype = inbdry.mytype; vbdry = inbdry.vbdry;}
        
        /* BASIC B.C. STUFF */
        void alloc(int n);
        virtual edge_bdry* create(tri_mesh &xin) const {
            return(new edge_bdry(*this,xin));
        }
        virtual void copy(const edge_bdry& bin);
        
        /* ADDITIONAL STUFF FOR SIDES */
        virtual void swap(int s1, int s2);
        virtual void reorder();
        virtual void mgconnect(Array<tri_mesh::transfer,1> &cnnct, tri_mesh& tgt, int bnum);
        virtual void mvpttobdry(int nel, FLT psi, TinyVector<FLT,tri_mesh::ND> &pt);
        virtual void findbdrypt(const TinyVector<FLT,2> xpt, int &sidloc, FLT &psiloc) const;
        
        /* DEFAULT SENDING FOR SIDE POINTS */
        virtual void vloadbuff(boundary::groups group,FLT *base,int bgn,int end, int stride) {}
        virtual void vfinalrcv(boundary::groups group,int phase, comm_type type, operation op, FLT *base,int bgn,int end, int stride) {}
        virtual void sloadbuff(boundary::groups group,FLT *base,int bgn,int end, int stride) {}
        virtual void sfinalrcv(boundary::groups group,int phase, comm_type type, operation op, FLT *base,int bgn,int end, int stride) {}
        virtual void loadpositions() {vloadbuff(all,&(x.pnts(0)(0)),0,tri_mesh::ND-1,tri_mesh::ND);}
        virtual void rcvpositions(int phase) {vfinalrcv(all_phased,phase,master_slave,replace,&(x.pnts(0)(0)),0,tri_mesh::ND-1,tri_mesh::ND);}
};

#endif
