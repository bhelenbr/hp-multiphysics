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

#ifdef SINGLE
#define EPSILON FLT_EPSILON
#define FLT float
#else
#define EPSILON DBL_EPSILON
#define FLT double
#endif

#ifndef ND
#define ND 2
#endif

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

class quadtree {
   private:
      int maxvrtx;
      FLT (*vrtx)[ND];
      class quad *base;
      class quad **indx;
      int size;
      int current;

/*		THIS IS USED BY ALL QUADS FOR SEARCHING */      
      static class quad **srchlst;
      static int maxsrch;

   public:
      quadtree() : size(0) {};
      void copy(const class quadtree& tgt);
      void allocate(FLT (*v)[ND],int mxv);
      inline void init(FLT (*v)[ND], int mxv, FLT x1, FLT y1, FLT x2, FLT y2) { allocate(v,mxv); init(x1,y1,x2,y2);}
      void init(); // RESETS WITH SAME AREA
      void init(FLT x1, FLT y1, FLT x2, FLT y2);  // RESETS QUAD WITH NEW AREA
      void reinit(); // ERASES TREE AND REINSERTS POINTS THUS REMOVING UNUSED QUADS
      
      inline void change_vptr(FLT (*v)[ND]) { vrtx = v;}
      inline FLT xmin() {return(base[0].xmin);}
      inline FLT ymin() {return(base[0].ymin);}
      inline FLT xmax() {return(base[0].xmax);}
      inline FLT ymax() {return(base[0].ymax);} 
           
      void addpt(int v0, class quad *start = NULL);
      
      FLT nearpt(int v0, int& pt) const;
      FLT nearpt(FLT x, FLT y, int& pt) const;
      
      void dltpt(int v0);
      
      void movept(int from, int to);
      
      void update(int bgn, int end);
      void update(int v0);
      
      enum FILETYPE {text,tecplot};
      void output(char *filename, FILETYPE type=tecplot);
};
#endif

