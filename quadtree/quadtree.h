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

#include <blitz/array.h>
using namespace blitz;

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
	TinyVector<FLT,ND> xmin, xmax;
	union {
		int node[(1<<ND)];     	// STORES EITHER NODES OR POINTER TO DAUGHTERS
		class box<ND> *dghtr[(1<<ND)];
	};
	
public:
	box() {}
	box(box *p, int i, TinyVector<FLT,ND> x1, TinyVector<FLT,ND> x2) : num(0), prnt(p), pind(i), xmin(x1), xmax(x2) {}
};

template<int ND> class quadtree {
private:
	int maxvrtx;
	Array<TinyVector<FLT,ND>,1> vrtx;
	class box<ND> *base;
	class box<ND> **indx;
	int size;
	int current;
	
	/*		THIS IS USED FOR SEARCHING */
	class box<ND> **srchlst;
	int maxsrch;
	
public:
	quadtree() : maxvrtx(0), vrtx(0), size(0), current(0) {};
	quadtree(const class quadtree& tgt) {
		size = 0;
		copy(tgt);
	}
	void copy(const class quadtree& tgt);
	void allocate(Array<TinyVector<FLT,ND>,1> v,int mxv);
	inline void init(Array<TinyVector<FLT,ND>,1> v, int mxv, TinyVector<FLT,ND> x1, TinyVector<FLT,ND> x2) { allocate(v,mxv); init(x1,x2);}
	void init(); // RESETS WITH SAME AREA
	void init(TinyVector<FLT,ND> x1, TinyVector<FLT,ND> x2);  // RESETS box WITH NEW AREA
	void reinit(); // ERASES TREE AND REINSERTS POINTS THUS REMOVING UNUSED boxS
	
	inline void change_vptr(Array<TinyVector<FLT,ND>,1> v) {vrtx.reference(v);}
	inline FLT xmin(int i) const {return(base[0].xmin[i]);}
	inline FLT xmax(int i) const {return(base[0].xmax[i]);}
	
	void addpt(int v0, class box<ND> *start = 0);
	
	FLT nearpt(int const v0, int& pt) const;
	FLT nearpt(const TinyVector<FLT,ND> x, int& pt) const;
	
	void dltpt(int v0);
	
	void movept(int from, int to);
	
	void update(int bgn, int end);
	void update(int v0);
	
	enum FILETYPE {text,tecplot};
	void output(const char *filename, FILETYPE type=tecplot);
	~quadtree() {
		delete []srchlst;
		delete []indx;
		delete []base;
	}
};
#endif

