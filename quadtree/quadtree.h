/*
 *  quadtree.h
 *  mblock
 *
 *  Created by helenbrk on Thu Aug 09 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#include<cstdio>
#ifndef _quadtree_h_
#define _quadtree_h_

#ifndef FLT
#define FLT double
#endif

#if (FLT == double)
#define EPSILON DBL_EPSILON
#else
#define EPSILON FLT_EPSILON
#endif

#ifndef ND
#define ND 2
#endif

class quadtree {
   private:
      int maxvrtx;
      FLT (*vrtx)[ND];
      class quad *base;
      class quad **indx;
      int size;
      int current;

   public:
      quadtree() : size(0) {};
      quadtree(const quadtree& copy) : size(0) {*this = copy;}
      quadtree& operator=(const quadtree& copy);
      void change_vptr(FLT (*v)[ND]) { vrtx = v;}
      void init(FLT (*v)[ND], int mxv, FLT x1, FLT y1, FLT x2, FLT y2);
      void init(FLT x1, FLT y1, FLT x2, FLT y2);
      void addpt(int v0, class quad *start = NULL);
      FLT nearpt(int v0, int& pt) const;
      FLT nearpt(FLT x, FLT y, int& pt) const;
      void dltpt(int v0);
      void movept(int from, int to);
      void output(void);
};

class quad {
   private:
      friend class quadtree;
      int num; 			// NUMBER OF POINTS
      class quad *prnt;	// POINTER TO PARENT QUAD
      int pind; 			// DAUGHTER NUMBER OF PARENT
      FLT xmin, ymin;
      FLT xmax, ymax;
      union {
         int node[4];    	// STORES EITHER NODES OR POINTER TO DAUGHTERS
         class quad *dghtr[4];
      };
      
   public:
      quad() {}
      quad(quad *p, int i, FLT x1, FLT y1, FLT x2, FLT y2) : 
         num(0), prnt(p), pind(i), xmin(x1), ymin(y1), xmax(x2), ymax(y2) {}
         
};
#endif