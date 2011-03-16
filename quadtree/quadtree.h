/*
 *  quadtree.h
 *  mblock
 *
 *  Created by helenbrk on Thu Aug 09 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#ifndef _quadtree_h_
#define _quadtree_h_

#ifdef SINGLE
#define EPSILON FLT_EPSILON
#define FLT float
#else
#ifndef FLT
#define EPSILON DBL_EPSILON
#define FLT double
#endif
#endif

template<int ND> class quadtree;

template<int ND> class box {
    private:
        friend class quadtree<ND>;
        int num; 			// NUMBER OF POINTS
        class box<ND> *prnt;	// POINTER TO PARENT QUAD
        int pind; 			// DAUGHTER NUMBER OF PARENT
        FLT xmin[ND], xmax[ND];
        union {
            int node[(1<<ND)];     	// STORES EITHER NODES OR POINTER TO DAUGHTERS
            class box<ND> *dghtr[(1<<ND)];
        };
        
    public:
        box() {}
        box(box *p, int i, FLT x1[ND], FLT x2[ND]) : num(0), prnt(p), pind(i) {
            for(int n=0;n<ND;++n) {
                xmin[n] = x1[n];
                xmax[n] = x2[n];
            }
        }
            
};

template<int ND> class quadtree {
    private:
        int maxvrtx;
        FLT (*vrtx)[ND];
        class box<ND> *base;
        class box<ND> **indx;
        int size;
        int current;

/*		THIS IS USED FOR SEARCHING */        
        class box<ND> **srchlst;
        int maxsrch;

    public:
        quadtree() : maxvrtx(0), vrtx(0), size(0), current(0) {};
        void copy(const class quadtree& tgt);
        void allocate(FLT (*v)[ND],int mxv);
        inline void init(FLT (*v)[ND], int mxv, FLT x1[ND], FLT x2[ND]) { allocate(v,mxv); init(x1,x2);}
        void init(); // RESETS WITH SAME AREA
        void init(FLT x1[ND], FLT x2[ND]);  // RESETS box WITH NEW AREA
        void reinit(); // ERASES TREE AND REINSERTS POINTS THUS REMOVING UNUSED boxS
        
        inline void change_vptr(FLT (*v)[ND]) { vrtx = v;}
        inline FLT xmin(int i) const {return(base[0].xmin[i]);}
        inline FLT xmax(int i) const {return(base[0].xmax[i]);}
              
        void addpt(int v0, class box<ND> *start = 0);
        
        FLT nearpt(int const v0, int& pt) const;
        FLT nearpt(FLT const x[ND], int& pt) const;
        
        void dltpt(int v0);
        
        void movept(int from, int to);
        
        void update(int bgn, int end);
        void update(int v0);
        
        enum FILETYPE {text,tecplot};
        void output(const char *filename, FILETYPE type=tecplot);
};
#endif

