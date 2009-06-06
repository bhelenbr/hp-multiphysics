#ifndef _tet_mesh_h_
#define _tet_mesh_h_

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
#include <blocks.h>

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
class face_bdry;

/** This is an unstructured tetrahedral mesh class which has adaptation and 
parallel communication capabilities */

class tet_mesh : public multigrid_interface {

	/***************/
	/*     DATA    */
	/***************/
	public: 
		int maxvst; /**< maximum length of point, segment, tri, or tet lists (all same for convenience) */
		static const int ND = 3;  /**< 3 dimensional */
				
		/** @name Variables describing mesh points */
		//@{
		int npnt;  /**< This is the total number of points in the mesh */
		Array<TinyVector<FLT,ND>,1> pnts;  /**< Physical location of the points in the mesh */
		Array<FLT,1> lngth; /**< Target mesh resolution in the neighborhood of each point */
		/** Data structure for integer data at each point */
		struct pntstruct {
			int tet;   /**< Pointer to a tet connected to this point */
			int seg;  /**< Pointer to a seg connected to this point */
			int nnbor; /**< Number of tets connected to this point */
			int nspk; /**< Number of segs connected to this point */
			int info; /**< General purpose (mostly for adaptation) */
		};
		Array<pntstruct,1> pnt; /**< Array of point integer data */
		class quadtree<ND> otree; /**< Octtree for finding points */
		//@}
	
		/* VERTEX BOUNDARY INFO */
		/** @name Variables describing vertex (0D) boundary conditions */
		//@{
		int nvbd; /**< number of vertex boundaries */
		Array<vrtx_bdry *,1> vbdry; /**< Array of vertex boundary objects */
		vrtx_bdry* getnewvrtxobject(int idnum, input_map& bdrydata); /**< function for obtaining different vertex boundary objects */
		//@}
		
		/* SEGMENT DATA */        
		/** @name Variables describing mesh segments */
		//@{
		int nseg;    /**< number of segments in mesh */
//        int nspk; /**< number of spokes for a given vertex */
		/** Data structure for each segment */
		struct segstruct {
			TinyVector<int,2> pnt;  /**< two points on a segment */
			int tet;   /**< one tet connected to a segment */
			int nnbor;    /**< number of neighboring tets */
			int info;   /**< General purpose (mostly for adaptation) */
		};
		Array<segstruct,1> seg; /**< Array of segment data */
		//@}
		
		/* EDGE BOUNDARY INFO */
		/** @name Variables describing edge (1D) boundary conditions */
		//@{
		int nebd; /**< number of edge boundaries */
		Array<edge_bdry *,1> ebdry; /**< array of edge boundary objects */
		edge_bdry* getnewedgeobject(int idnum, input_map& bdrydata); /**< function for obtaining different edge boundary objects */
		//@}
		
		/* TRIANGLE DATA */  
		/** @name Variables describing triangles */
		//@{       
		int ntri;            /**< Number of triangles in mesh */
		struct tristruct {
			TinyVector<int,3> pnt;   /**< triangle points */
			TinyVector<int,3> seg;  /**< triangle segments */
			TinyVector<int,2> tet;  /**< 2 tets sharing common tri */
			TinyVector<int,3> sgn;  /**< sign convention for each segment on tri */
			int info; /**< General purpose (mostly for adaptation) */
		};
		Array<tristruct,1> tri; /**< Array of triangle data */
		
		/* FACE BOUNDARY INFO */
		/** @name Variables describing face (2D) boundary conditions */
		//@{
		int nfbd; /**< number of face boundaries */
		Array<face_bdry *,1> fbdry;  /**< array of face boundary objects */
		face_bdry* getnewfaceobject(int idnum, input_map& bdrydata); /**< function for obtaining different face boundary objects */
		int getbdrynum(int tetnum) const { return((-tetnum>>16) -1);}  /**< Uses info in seg.tri or tri.tri to determine boundary object number */
		int getbdryseg(int tetnum) const { return(-tetnum&0xFFFF);}  /**< Uses info in seg.tri or tri.tri to determine boundary element */
		int tetnumatbdry(int bnum, int bel) const { return(-(((bnum+1)<<16) +bel));} /**< Combines bnum & bel into 1 integer for storage in boundary of seg.tri or tri.tri */
		//@}
		
		/* TETRAHEDRAL DATA */    
		/** @name Variables describing tetrahedrals */
		//@{     
		int ntet; /**< Number of tetrahedrals in mesh*/
		struct tstruct {
			TinyVector<int,4> pnt;       /**< tetrahedral points*/  
			TinyVector<int,6> seg;     /**< tetrahedral segments*/ 
			TinyVector<int,4> tri;     /**< tetrahedral triangles*/ 
			TinyVector<int,4> rot;     /**< face orientation + ccw - cw*/
			TinyVector<int,4> tet;     /**< four tet's connected to a tet*/  
			TinyVector<int,6> sgn;     /**< sign convention for segment on tet*/ 
			int info; /**< General purpose (mostly for adaptation) */
			FLT vol; /**< volume of a linear tetrahedral */
		};
		Array<tstruct,1> tet; /**< Array of tetrahedral data */
		//@}
		
		/** /struct For information shared between meshes not used simultaneously (multigrid levels) */
		struct global : public block_global {            
			Array<int,1> i1wk; /**< Integer work array, any routine that uses i1wk should reset it to -1 */
			Array<FLT,1> fltwk; /**< Floating point work array */
			int nlst; /**< variable to keep track of number of entities in list */
			Array<int,1> i2wk,i2wk_lst1, i2wk_lst2, i2wk_lst3; /**< Arrays for storing lists */
			int maxsrch; /**< integer describing maximum number of triangles to search before giving up */
		} *gbl;

		bool initialized;

		/**************/
		/*  INTERFACE */
		/**************/
		/* INITIALIZATION & ALLOCATION */
		tet_mesh() : nvbd(0), nebd(0), nfbd(0), initialized(0)  {}
		/** Routine to allocate shared variables */
		void* create_global_structure() {return new global;}
		/** Routine to initialize with using information in map and shared resource in gbl_in */
		void init(input_map& input, void *gbl_in);
		/** Routine to initialze from another mesh with option of increasing or decreasing storage (compatible with block.h) */
		void init(const multigrid_interface& mgin, init_purpose why=duplicate, FLT sizereduce1d=1.0);
		/** Routine to copy */
		void copy(const tet_mesh& tgt);
		
		/** Innput/Output file types */
		enum filetype {easymesh, baker, gambit, tecplot, grid, text, binary, BRep, mavriplis, boundary, vlength, debug_adapt, datatank};
		
		/** Input mesh */
		void input(const std::string &filename, filetype ftype,  FLT grwfac, input_map &input);
		/** Virtual routine so inheritors can set up info after input/adaptation */
		virtual void setinfo();  

		/** Outputs solution in various filetypes */
		void output(const std::string &outname,block::output_purpose why) {output(outname,grid);}
		void output(const std::string &filename, filetype ftype = grid) const;

		/** @name adapt adaptation routines */
		//@{  
		void initlngth();  /**< Set target mesh resolution based on current mesh */
		virtual void length() {} /**< Virtual so inheritors can set target resolution before adaptation */
//        void adapt(); /**< Adapt mesh */
		//@}

		/* Destructor (deletes boundary objects) */
		virtual ~tet_mesh();
		
		void test();
		void vertexball(int vind);
		//void spokes(int vind);

		
		/** @name Mesh modification utilties */
		//@{   
//        void symmetrize();  /**< Creates a symmetric mesh about y = 0 */
//        void cut(); /**< Cut's mesh based on indicator array inf fltwk (Positive / Negative at each point) and aligns cut edge along 0 */
//        void trim(); /**< Starting from boundaries deletes triangles based on sign of fltwk */
//        void append(const tet_mesh &z); /**< Appends mesh */
		void shift(TinyVector<FLT,ND>& s); /**< Tranlates mesh */
		void scale(TinyVector<FLT,ND>& s); /**< Scales mesh */
		int smooth_cofa(int niter); /**< Does a center of area smoothing */
		void refineby2(const class tet_mesh& xmesh); /**< Refines by 2 */
//        void coarsen_substructured(const class tet_mesh &tgt,int p); /**< Coarsens mesh that was output by hp FEM method */ 
		//@}

//         /** @name Routines for parallel computations */
//        //@{ 
//#ifdef METIS
//        void setpartition(int nparts);  /**< Set partition of mesh (in tri(i).info) */
//#endif 
//        void partition(class tet_mesh& xmesh, int npart); /**< Creates a partition from xmesh */
		int comm_entity_size(); /**< Returns size of list of communication entities (for blocks.h) */
		int comm_entity_list(Array<int,1>& list);  /**< Returns list of communication entities */
		class boundary* getvbdry(int num); /**< Returns pointer to vertex boundary (for blocks.h) */
		class boundary* getebdry(int num); /**< Returns pointer to edge boundary (for blocks.h) */
		class boundary* getfbdry(int num); /**< Returns pointer to edge boundary (for blocks.h) */
		void create_unique_numbering();
		void match_bdry_numbering();
		void pmsgload(boundary::groups group, int phase, boundary::comm_type type, FLT *base, int bgn, int end, int stride); /**< Loads message buffers by point index from base (vertex, edge, and face boundaries) */
		void pmsgpass(boundary::groups group,int phase, boundary::comm_type type); /**< Sends vertex and edge face point boundary messages */
		int pmsgwait_rcv(boundary::groups group, int phase, boundary::comm_type type, boundary::operation op, FLT *base,int bgn, int end, int stride); /**< Receives all point messages (vertex, edge, and face boundaries) */
		int pmsgrcv(boundary::groups group, int phase, boundary::comm_type type, boundary::operation, FLT *base,int bgn, int end, int stride);/**< Receives without waiting (Don't use) */
		void smsgload(boundary::groups group, int phase, boundary::comm_type type, FLT *base,int bgn, int end, int stride);/**< Loads message buffers by segment index from base (edge & face boundaries) */
		void smsgpass(boundary::groups group, int phase, boundary::comm_type type);/**< Sends edge & face segment boundary messages */
		int smsgwait_rcv(boundary::groups group,int phase, boundary::comm_type type, boundary::operation op, FLT *base,int bgn, int end, int stride); /**< Receives segment messages */
		int smsgrcv(boundary::groups group,int phase, boundary::comm_type type,  boundary::operation op, FLT *base,int bgn, int end, int stride); /**< Receives without waiting (Don't use) */
		void tmsgload(boundary::groups group, int phase, boundary::comm_type type, FLT *base,int bgn, int end, int stride);/**< Loads message buffers by triangle index from base (face boundaries) */
		void tmsgpass(boundary::groups group, int phase, boundary::comm_type type);/**< Sends triangle boundary messages */
		int tmsgwait_rcv(boundary::groups group, int phase, boundary::comm_type type, boundary::operation op, FLT *base,int bgn, int end, int stride);/**< Receives triangle messages */
		int tmsgrcv(boundary::groups group,int phase, boundary::comm_type type, boundary::operation op, FLT *base,int bgn, int end, int stride); /**< Receives without waiting (Don't use) */        
		void matchboundaries();  /**< Matches location of boundary points */
		//@}
		
		/** @name Routines for multigrid */
		//@{ 
		int coarse_level; /**< Integer telling which multigrid level I am */
		multigrid_interface *fine, *coarse; /**< Pointers to fine and coarse objects */
//        void coarsen(FLT factor, const class tet_mesh& xmesh); /**< Coarsens by triangulating then inserting */
//        void coarsen2(FLT factor, const class tet_mesh& inmesh, FLT size_reduce = 1.0);  /**< Coarsens by using yaber */
//        void coarsen3(); /**< Coarsens based on marks stored in pnt().info */

		/** Structure for mesh connection data */
		struct transfer {
			int tet; /**< Tri containing this point */
			TinyVector<FLT,4> wt; /**< Interpolation weights */
		};
		Array<transfer,1> fcnnct, ccnnct; /**< Arrays for connection data for each point */
		void connect(multigrid_interface& tgt);  /**< Set-up connections  to fine mesh compatible with block.h */
		void mgconnect(tet_mesh &tgt, Array<transfer,1> &cnnct); /**< Utility called by connect */
		void testconnect(const std::string &fname,Array<transfer,1> &cnnct, tet_mesh *cmesh); /**< Tests by outputing grids with interpolated point locations */
		//@}

		/* SOME DEGUGGING FUNCTIONS */
		void checkintegrity(); /**< Checks data arrays for compatibility */
		void checki1wk() const; /**< Makes sure i1wk is all -1 */

		
	protected:
		/** @name Setup routines */
		//@{
		void allocate(int mxsize);  /**< Allocates memory */
		void cnt_nbor(void); /**< Fills in pnt().nnbor */
		void bdrylabel(void); /**< Makes seg().tri and tri().tri on boundary have pointer to boundary group/element */
		void createseg(void); /**< Creates all segment information from list of triangle points and also tri().seg/sgn (if necessary) */
		void createsegtri(void); /**< Creates seg.tri() and tri().seg connections */
		void createtritri(void); /**< Creates tri().tri data */
		void treeinit(); /**< Initialize octtree data */
		void treeinit(FLT x1[ND], FLT x2[ND]); /** Initialize octtree data with specified domain */
		//@}
			
		/*******************/
		/* INTERNAL FUNCTIONS */
		/*******************/  
		void fixvertexinfo(void);  /**< reorients vertex ordering, and calculates volume */
		void feedinvertexinfo(void); /**< reorients vertex ordering based on pnt.info and calculates volume */
		void vertexnnbor(void);  /**< counts number of neighboring tet's connected to a point pnt().nnbor*/
		void createedgeinfo(void);  /**<  */
		void edgeinfo(void);
		void createfaceinfo(void);
		void faceinfo(void);
		void morefaceinfo(void);
		void createtetinfo(void);
		void createtdstri(void);
		void createttet(void);
		void createvtet(void);
		void ring(int eind);
		void createfaceorientation(void);

//        /** @name Mesh modification functions */
//        //@{
//         /** 4 binary digits  used for each pnt, seg, tri adapt data 
//         *  Typically special,deleted,touched,searched 
//         */
//        const static int PSPEC = 0x4, PDLTE = 0x2, PTOUC=0x1;
//        const static int SDLTE = 0x10*0x2, STOUC=0x10*0x1;
//        const static int TSRCH = 0x100*0x4, TDLTE = 0x100*0x2, TTOUC=0x100*0x1;
//        void setup_for_adapt(); /**< Set all flags */
//
//        void triangulate(int nseg); /**< Creates an initial triangulation */
//        void addtri(int p0,int p1,int p2,int sind,int dir); /**< Utility for creating triangles used by triangulate */
//        
//        void swap(FLT swaptol = EPSILON); /**< Edge swap */
//        int swap(int sind, FLT tol = 0.0);  /**< Swaps a single segment */
//        void swap(int nswp, int *swp, FLT tol = 0.0); /**< Swaps a list of segments */
//        
//        void bdry_yaber(FLT tolsize); /**< Coarsen edge boundaries */
//        void bdry_yaber1(); /**< Coarsen slave edge boundaries */
//        void yaber(FLT tolsize); /**< Coarsen mesh (Rebay backwards) */
//        void collapse(int sind, int endpt); /**< Removes by collapsing segment sind to endpt */
//
//        void bdry_rebay(FLT tolsize); /**< Refine edges */
//        void bdry_rebay1(); /**< Refine slave edges */
//        void rebay(FLT tolsize); /**< Refine mesh using Rebay point placement */
//        int insert(const TinyVector<FLT,ND> &x);  /**< Inserts a point */
//        int insert(int pnum, int tnum);  /**< Inserts point at pnum into tnum */
//        void bdry_insert(int pnum, int sind, int endpt = 0); /**< Inserts a boundary point in segment sind if endpt is 0 makes left old and right new seg */
//        int findtri(TinyVector<FLT,ND> x, int pnear); /**< Locate triangle containing point with initial seed of pnear */
//        
//        void cleanup_after_adapt(); /**< Clean up and move data etc.. */
//        void dltpnt(int pind); /**< Removes leftover point references */
//        virtual void movepdata(int frm, int to) {} /**< Virtual routine so inheritors can automatically have their point data moved */
//        virtual void movepdata_bdry(int bnum,int bel,int endpt) {} /**< Virtual routine so inheritors can automatically have their boundary point data moved */ 
//        virtual void updatepdata(int v) {} /**< Virtual routine so inheritors can automatically have new point data updated */
//        virtual void updatepdata_bdry(int bnum,int bel,int endpt) {} /**< Virtual routine so inheritors can automatically have new bondary point data updated */
//        void dltseg(int sind); /**< Removes leftover segment references */
//        virtual void movesdata(int frm, int to) {} /**< Virtual routine so inheritors can automatically have their segment data moved */
//        virtual void movesdata_bdry(int bnum,int bel) {} /**< Virtual routine so inheritors can automatically have their boundry segment data moved  */
//        virtual void updatesdata(int s) {} /**< Virtual routine so inheritors can update segment data */
//        virtual void updatesdata_bdry(int bnum,int bel) {} /**< Virtual routine so inheritors can update boundary segment data */
//        void dlttri(int tind); /**< Removes leftover triangle references */
//        virtual void movetdata(int frm, int to) {} /**< Virtual routine so inheritors can automatically have their triangle data moved */
//        virtual void updatetdata(int t) {}  /**< Virtual routine so inheritors can automatically have new triangles updated */
//        
//        void putinlst(int sind);  /**< For inserting into adaptation priority queues */
//        void tkoutlst(int sind);  /**< For removing from adaption priority queues */ 
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
		bool findtet(const TinyVector<FLT,3> xp, int seedvrtx, int &tind);
		bool findtet(const TinyVector<FLT,3> xp, int &tind);
		FLT intet(int tind, const TinyVector<FLT,ND> &x);  /**< Determine whether a triangle contains a point (0.0 or negative) */
		TinyVector<FLT,3> tet_wgt; /**< Used in intri for searching (normalized value returned by getwgts after successful search) */
		void getwgts(TinyVector<FLT,ND+1> &wgt) const; /**< Returns weighting for point interpolation within last triangle find by intri */
		//@}
};

/** \brief Specialization for a vertex 
 *
 * \ingroup boundary
 * Specialization of communication routines and 
 * and storage for a vertex boundary 
 */
class vrtx_bdry : public boundary {
	public:
		tet_mesh &x;
		TinyVector<int,2> ebdry;
		int pnt;
		
		/* CONSTRUCTOR */
		vrtx_bdry(int intype, tet_mesh &xin) : boundary(intype), x(xin) {idprefix = x.gbl->idprefix +"_v" +idprefix; mytype="plain";}
		vrtx_bdry(const vrtx_bdry &inbdry, tet_mesh &xin) : boundary(inbdry.idnum), x(xin)  {idprefix = inbdry.idprefix; mytype = inbdry.mytype; ebdry = inbdry.ebdry;}

		/* OTHER USEFUL STUFF */
		virtual vrtx_bdry* create(tet_mesh &xin) const { return(new vrtx_bdry(*this,xin));}
		virtual void copy(const vrtx_bdry& bin) {
			pnt = bin.pnt;
			ebdry = bin.ebdry;
		}
		
		/* INPUT/OUTPUT NOT USED YET // FIXME */
		virtual void input(istream &fin,tet_mesh::filetype type = tet_mesh::grid) {}
		virtual void output(ostream &fin,tet_mesh::filetype type = tet_mesh::grid) const {}
		
		
		virtual void ploadbuff(boundary::groups group, FLT *base, int bgn, int end, int stride) {}
		virtual void pfinalrcv(boundary::groups group, int phase, comm_type type, operation op, FLT *base, int bgn, int end, int stride) {}

		virtual void mvpttobdry(TinyVector<FLT,tet_mesh::ND> &pt) {}
		virtual void loadpositions() {ploadbuff(all,&(x.pnts(0)(0)),0,tet_mesh::ND-1,tet_mesh::ND);}
		virtual void rcvpositions(int phase) {pfinalrcv(all_phased,phase,master_slave,replace,&(x.pnts(0)(0)),0,tet_mesh::ND-1,tet_mesh::ND);}
};


/** \brief Specialization for a edge 
 *
 * \ingroup boundary
 * Specialization of communication routines and 
 * and storage for a side boundary 
 */
class edge_bdry : public boundary {
	public:
		tet_mesh &x;
		  TinyVector<int,2> vbdry;  // FIXME
		int maxseg;
		int nseg;
		  struct segstruct {
				int next; /**< Not used except in adaptation (kept ordered) */
				int prev; /**< Not used except in adaptation (kept ordered) */
			int gindx;  /**< global index of side */
			int info; /**< General purpose (mostly for adaptation) */
		};
		  Array<segstruct,1> seg;

		/* CONSTRUCTOR */
		edge_bdry(int inid, tet_mesh &xin) : boundary(inid), x(xin), maxseg(0)  {idprefix = x.gbl->idprefix +"_e" +idprefix; mytype="plain"; vbdry = -1;}
		edge_bdry(const edge_bdry &inbdry, tet_mesh &xin) : boundary(inbdry.idnum), x(xin), maxseg(0)  {idprefix = inbdry.idprefix; mytype = inbdry.mytype; vbdry = inbdry.vbdry;}
		
		/* BASIC B.C. STUFF */
		void alloc(int n);
		virtual edge_bdry* create(tet_mesh &xin) const {
			return(new edge_bdry(*this,xin));
		}
		virtual void copy(const edge_bdry& bin);
		
		/* INPUT/OUTPUT NOT USED YET // FIXME */
		virtual void input(istream &fin,tet_mesh::filetype type = tet_mesh::grid) {
			if (type == tet_mesh::grid) {
				for(int j=0;j<nseg;++j) {
					fin.ignore(80,':');
					fin >> seg(j).gindx;
					fin.ignore(80,'\n');
				}
			}
		}
		virtual void output(ostream &fin,tet_mesh::filetype type = tet_mesh::grid) const {}

		
		/* ADDITIONAL STUFF FOR EDGES */
		  virtual void match_numbering(int step) {}
		virtual void swap(int s1, int s2);
		  virtual void setup_next_prev();
		virtual void reorder();
		virtual void mgconnect(Array<tet_mesh::transfer,1> &cnnct,tet_mesh& tgt, int bnum);
		virtual void mvpttobdry(int nseg, FLT psi, TinyVector<FLT,tet_mesh::ND> &pt);
		virtual void findbdrypt(const TinyVector<FLT,tet_mesh::ND> xpt, int &sidloc, FLT &psiloc) const;
		
		/* DEFAULT SENDING FOR SIDE VERTICES */
		virtual void ploadbuff(boundary::groups group,FLT *base,int bgn,int end, int stride) {}
		virtual void pfinalrcv(boundary::groups group,int phase, comm_type type, operation op, FLT *base,int bgn,int end, int stride) {}
		virtual void sloadbuff(boundary::groups group,FLT *base,int bgn,int end, int stride) {}
		virtual void sfinalrcv(boundary::groups group,int phase, comm_type type, operation op, FLT *base,int bgn,int end, int stride) {}
		virtual void loadpositions() {ploadbuff(all,&(x.pnts(0)(0)),0,tet_mesh::ND-1,tet_mesh::ND);}
		virtual void rcvpositions(int phase) {pfinalrcv(all_phased,phase,master_slave,replace,&(x.pnts(0)(0)),0,tet_mesh::ND-1,tet_mesh::ND);}
};

/** \brief Specialization for a face 
 *
 * \ingroup boundary
 * Specialization of communication routines and 
 * and storage for a side boundary 
 */
class face_bdry : public boundary {
	/***************/
	/* DATA        */
	/***************/
	public: 
		tet_mesh &x;
		int maxpst;  /**< maximum length of point, segment, or tri lists (all same for convenience) */
				
		/** @name Variables describing mesh points */
		//@{
		int npnt;  /**< Total number of points in the mesh */
		/** Data structure for integer data at each point */
		struct pntstruct {
			int tri; /**< Pointer to a triangle connected to this point */
			int nnbor; /**< Number of segments connected to this point */
			int gindx;
			int info; /**< General purpose (mostly for adaptation) */
		};
		Array<pntstruct,1> pnt;  /**< Array of point integer data */
		class quadtree<tet_mesh::ND> otree; /**< Octtree for finding points */  //FIXME
		//@}

		/** @name Variables describing mesh segments */
		//@{    
		int nseg; /**< number of segments in mesh */
		/** Data structure for each segment */
		struct segstruct {
			TinyVector<int,2> pnt; /**< enpoints */
			TinyVector<int,2> tri; /**< adjacent triangles */
			int gindx;
			int info; /**< General purpose (mostly for adaptation) */
		};
		Array<segstruct,1> seg; /**< Array of segment data */
		//@}

		/** @name Variables describing edges (1D) of face on this block */
		//@{
		int nebd; /**< number of edge boundaries */  // FIXME
		Array<int,1> ebdry; /**< array of edge numbers */ // FIXME
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
			int gindx;
			int info; /**< General purpose (mostly for adaptation) */
		};
		Array<tristruct,1> tri;
		//@}
				
		/**************/
		/*  INTERFACE */
		/**************/
		/* INITIALIZATION & ALLOCATION */
		face_bdry(int inid, tet_mesh &xin) : boundary(inid), x(xin) {idprefix = x.gbl->idprefix +"_f" +idprefix; mytype="plain";}
		face_bdry(const face_bdry &inbdry, tet_mesh &xin) : boundary(inbdry.idnum), x(xin) {idprefix = inbdry.idprefix; mytype = inbdry.mytype;}
		void alloc(int mxsize);      
		virtual face_bdry* create(tet_mesh &xin) const {
			return(new face_bdry(*this,xin));
		}
		virtual void copy(const face_bdry& bin);

		/* INPUT/OUTPUT // FIX ME: NOT USED YET */
		virtual void input(istream &fin, tet_mesh::filetype type = tet_mesh::grid) {}
		virtual void output(ostream &fout, tet_mesh::filetype type = tet_mesh::grid) const {}
		
		/*********************/
		/* SOME UTITILITIES */
		/*******************/  
		virtual void match_numbering(int step) {}
		void vertexcircle(int vind);
		void cnt_nbor(void);
		void createsideinfo(void);
		void createtdstri(void);
		void createttri(void);
		void createvtri(void);
		void gbltolclvrtx(void);
		void gbltolclside(void);
		void gbltolcltri(void);
		void allinfo(void);
		void treeinit();
		void treeinit(FLT x1[tet_mesh::ND], FLT x2[tet_mesh::ND]);
		virtual void mgconnect(Array<tet_mesh::transfer,1> &cnnct,tet_mesh& tgt, int bnum);
		virtual void findbdrypt(const TinyVector<FLT,tet_mesh::ND> xpt, int &facloc, FLT &r, FLT &s) const;
		virtual void mvpttobdry(int ntri, FLT r, FLT s, TinyVector<FLT,tet_mesh::ND> &pt);
		
//        /* MESH MODIFICATION FUNCTIONS */
//        /* 4 binary digits for vertex, side, tri */
//        const static int VSPEC = 0x4, VDLTE = 0x2, VTOUC=0x1;
//        const static int SDLTE = 0x10*0x2, STOUC=0x10*0x1;
//        const static int TSRCH = 0x100*0x4, TDLTE = 0x100*0x2, TTOUC=0x100*0x1;
//        void setup_for_adapt();
//        void swap(FLT swaptol = EPSILON);
//        void bdry_yaber(FLT tolsize);
//        void bdry_yaber1();
//        void yaber(FLT tolsize);
//        void bdry_rebay(FLT tolsize);
//        void bdry_rebay1();
//        void rebay(FLT tolsize);
//        void cleanup_after_adapt();
//     
//        /* TO CREATE AN INITIAL TRIANGUlATION */
//        void triangulate(int nseg);
//        void addtri(int v0,int v1,int v2,int sind,int dir);
//
//        /* TO INSERT A POINT */
//        int insert(const TinyVector<FLT,ND> &x);
//        int insert(int vnum, int tnum);
//        void bdry_insert(int vnum, int sind, int endpt = 0);
//        int findtri(TinyVector<FLT,ND> x, int vnear) const;
//
//        /* FOR COARSENING A SIDE */
//        void collapse(int sind, int endpt);
//
//        void dltvrtx(int vind);
//        virtual void movevdata(int frm, int to) {}
//        virtual void movevdata_bdry(int bnum,int bel,int endpt) {}
//        virtual void updatevdata(int v) {}
//        virtual void updatevdata_bdry(int bnum,int bel,int endpt) {}
//        void dltsd(int sind);
//        virtual void movesdata(int frm, int to) {}
//        virtual void movesdata_bdry(int bnum,int bel) {}
//        virtual void updatesdata(int s) {}
//        virtual void updatesdata_bdry(int bnum,int bel) {}
//        void dlttri(int tind);
//        virtual void movetdata(int frm, int to) {}
//        virtual void updatetdata(int t) {}
//        
//        /* ORDERED LIST FUNCTIONS */
//        void putinlst(int sind);
//        void tkoutlst(int sind);
//        
//        /* TO EDGE SWAP FOR IMPROVED QUALITY */
//        int swap(int sind, FLT tol = 0.0);  //SWAPS A SINGLE SIDE 
//        void swap(int nswp, int *swp, FLT tol = 0.0); //SWAPS A LIST OF SIDES
//
//        /* PRIMITIVE FUNCTIONS */
//    public:
//        inline FLT distance(int v0, int v1) const {
//            FLT d = 0.0;
//            for(int n = 0; n<ND;++n)
//                d += pow(pnts(v0)(n) -pnts(v1)(n),2);
//            return(sqrt(d));
//        }
//        inline FLT distance2(int v0, int v1) const {
//            FLT d = 0.0;
//            for(int n = 0; n<ND;++n)
//                d += pow(pnts(v0)(n) -pnts(v1)(n),2);
//            return(d);
//        }  
//        FLT incircle(int tind, const TinyVector<FLT,ND> &x) const;
//        FLT insidecircle(int sind, const TinyVector<FLT,ND> &x) const;
//        FLT area(int sind, int vind) const;
//        FLT area(int v0, int v1, int v2) const;
//        FLT area(int tind) const;
//        FLT minangle(int v0, int v1, int v2) const;
//        FLT angle(int v0, int v1, int v2) const;
//        FLT circumradius(int tind) const;
//        void circumcenter(int tind, TinyVector<FLT,2> &x) const;
//        FLT inscribedradius(int tind) const;
//        FLT aspect(int tind) const;
//        FLT inttri(int tind, const TinyVector<FLT,tet_mesh::ND> &x) const;
//        void getwgts(TinyVector<FLT,tet_mesh::ND+1> &wgt) const;
//        /* tri numbers at boundary point to el and group */
//        int getbdrynum(int trinum) const { return((-trinum>>16) -1);}
//        int getbdryel(int trinum) const { return(-trinum&0xFFFF);}
//        int trinumatbdry(int bnum, int bel) const { return(-(((bnum+1)<<16) +bel));}
		
		/* DEFAULT SENDING FOR SIDE VERTICES */
		virtual void ploadbuff(boundary::groups group,FLT *base,int bgn,int end, int stride) {}
		virtual void pfinalrcv(boundary::groups group,int phase, comm_type type, operation op, FLT *base,int bgn,int end, int stride) {}
		virtual void sloadbuff(boundary::groups group,FLT *base,int bgn,int end, int stride) {}
		virtual void sfinalrcv(boundary::groups group,int phase, comm_type type, operation op, FLT *base,int bgn,int end, int stride) {}
		virtual void tloadbuff(boundary::groups group,FLT *base,int bgn,int end, int stride) {}
		virtual void tfinalrcv(boundary::groups group,int phase, comm_type type, operation op, FLT *base,int bgn,int end, int stride) {} 
		virtual void loadpositions() {ploadbuff(all,&(x.pnts(0)(0)),0,tet_mesh::ND-1,tet_mesh::ND);}
		virtual void rcvpositions(int phase) {pfinalrcv(all_phased,phase,master_slave,replace,&(x.pnts(0)(0)),0,tet_mesh::ND-1,tet_mesh::ND);}



};
#endif
