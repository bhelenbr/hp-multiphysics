/*
 *  quad.cpp
 *  mblock
 *
 *  Created by helenbrk on Wed Aug 08 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */
#include"quadtree.h"
#include<utilities.h>
#include<float.h>
#include<cstdio>
#include<assert.h>
#include<string.h>
#include<math.h>

/* STATIC VARIABLES FOR SEARCHING */ 
class quad **quadtree::srchlst;
int quadtree::maxsrch = 0;
  
void quadtree::allocate(FLT (*v)[ND], int mxv) {
   int i;
   
   vrtx = v;
   maxvrtx = mxv;
   if (size != 0) {
      delete [] base;
      delete [] indx;
   }
   
/*	TAKES APPROXIMATELY 1:1 QUADS TO VERTICES WITH 4 NODES / QUAD */
/*	+10 IS FACTOR OF SAFETY FOR SMALL mxv */
   size = (int) 1.1*mxv +10; 
   base = new class quad[size];
   indx = new class quad*[maxvrtx];
   for(i=0;i<maxvrtx;++i)
      indx[i] = NULL;

   base[0] = quad(NULL,0,0.0,0.0,1.0,1.0);
   current = 1;
   
   if (maxsrch == 0) {
      srchlst = new class quad*[size];
      maxsrch = size;
   }
   else {
      if (size > maxsrch) { 
         printf("#Warning: Better to allocate largest quadtree to smallest\n");
         delete []srchlst;
         srchlst = new class quad*[size];
         maxsrch = size;
      }
   }
}

void quadtree::init(FLT xmin, FLT ymin, FLT xmax, FLT ymax) {

   if (xmin >= xmax || ymin >= ymax) {
      printf("quadtree initialization error: zero domain area %f %f %f %f\n",xmin,xmax,ymin,ymax);
      exit(1);
   }

//   base[0] = quad(NULL,0,xmin,ymin,xmax,ymax);

   base[0].num = 0; base[0].prnt = NULL; base[0].pind = 0; base[0].xmin = xmin; base[0].ymin = ymin; base[0].xmax = xmax; base[0].ymax = ymax;

   for(int i = 0; i<maxvrtx; ++i)
      indx[i] = NULL;
      
   current = 1;
}

void quadtree::init() {

//   base[0] = quad(NULL,0,base[0].xmin,base[0].ymin,base[0].xmax,base[0].ymax);
   base[0].num = 0; base[0].prnt = NULL; base[0].pind = 0;
   for(int i = 0; i<maxvrtx; ++i)
      indx[i] = NULL;
      
   current = 1;
}

void quadtree::copy(const class quadtree& tgt) {
   int i,j;
      
   vrtx = tgt.vrtx;
   if (size == 0) {
      maxvrtx = tgt.maxvrtx;
      size = tgt.size;
      base = new class quad[size];
      indx = new class quad*[maxvrtx];
   }
   else {
      assert(size >= tgt.current);
   }

   if (maxvrtx < tgt.maxvrtx) {
      for(i=0;i<maxvrtx;++i) {
         if (tgt.indx[i] == NULL) indx[i] = NULL;
         else indx[i] = base + (tgt.indx[i] -tgt.base);
      }
   
      for(i=maxvrtx;i<tgt.maxvrtx;++i)
         if (tgt.indx[i] != NULL) {
            printf("quadtree is too small for copy\n");
            exit(1);
      }
   }
   else {
      for(i=0;i<tgt.maxvrtx;++i) {
         if (tgt.indx[i] == NULL) indx[i] = NULL;
         else indx[i] = base + (tgt.indx[i] -tgt.base);
      }
      for(i=tgt.maxvrtx;i<maxvrtx;++i)
         indx[i] = NULL;
   }

   current = tgt.current;
   
   for(i=0;i<current;++i) {
      base[i].num = tgt.base[i].num;
      base[i].prnt = base +(tgt.base[i].prnt -tgt.base);
      base[i].pind = tgt.base[i].pind;
      base[i].xmin = tgt.base[i].xmin;
      base[i].ymin = tgt.base[i].ymin;
      base[i].xmax = tgt.base[i].xmax;
      base[i].ymax = tgt.base[i].ymax;
      if (base[i].num > 0) {
         for(j=0;j<4;++j)
            base[i].node[j] = tgt.base[i].node[j];
      }
      else {
         for(j=0;j<4;++j)
            base[i].dghtr[j] = base +(tgt.base[i].dghtr[j] -tgt.base);
      }
   }
   
   base[0].prnt = NULL;
   
   return;
}


void quadtree::reinit() {
   int i;
      
   assert(size != 0); 
   
   base[0].num = 0;
   current = 1;
   
   for(i=0;i<maxvrtx;++i)
      if (indx[i] != NULL) addpt(i);
   
   return;
}



void quadtree::addpt(int v0, class quad* start) {
   int i,j;  // DON'T MAKE THESE STATIC SCREWS UP RECURSION
   class quad *qpt;  // SAME
   int store[5];
   FLT avoid0,xshift,yshift,x1,y1,dx,dy;
   
   if (start == NULL) 
      qpt = base;
   else
      qpt = start;

/*	LOOP TO BOTTOM OF TREE */
   while (qpt->num < 0) {
      xshift = vrtx[v0][0] -0.5*(qpt->xmax +qpt->xmin);
      avoid0 = MAX(fabs(xshift), 10.*EPSILON);
      i = ((int) ((1+10.*EPSILON)*xshift/avoid0) +1)/2;
      yshift = vrtx[v0][1] -0.5*(qpt->ymax +qpt->ymin);
      avoid0 = MAX(fabs(yshift), 10.*EPSILON);
      i = 2*i +((int) ((1+10.*EPSILON)*yshift/avoid0) +1)/2;
      qpt = qpt->dghtr[i];
   }
   
   if (qpt->num < 4) {
/*		QUAD CAN ACCEPT NEW POINT */
      qpt->node[qpt->num++] = v0;
      indx[v0] = qpt;
      return;
   }

/*	QUAD IS FULL SUBDIVIDE */
   if (current +4 >= size) {
      printf("Need to allocate bigger quadtree %d\n",size);
      output("quad_error");
      exit(1);
   }
   for (i=0;i<4;++i)
      store[i] = qpt->node[i];
   store[4] = v0;
   
   qpt->num = -1;
   dx = (qpt->xmax -qpt->xmin)*0.5;
   dy = (qpt->ymax -qpt->ymin)*0.5;
      
   for (i=0;i<2;++i) {
      for(j=0;j<2;++j) {
         x1 = qpt->xmin +dx*i;
         y1 = qpt->ymin +dy*j;
         base[current] = quad(qpt,2*i +j,x1,y1,x1+dx,y1+dy);
         qpt->dghtr[2*i +j] = base +current;
         ++current;
      }
   }
   
   for(i=0;i<5;++i)
      addpt(store[i],qpt);
      
   return;
}

/*	FIND CLOSEST POINT TO A POINT IN THE ARRAY */
FLT quadtree::nearpt(int v0, int& pt) const {
   int i,nsrch,exclude;
   class quad *topbox;
   FLT dist,dx,dy,xh,xl,yh,yl,mindist;
   class quad *qpt;

   xl = base[0].xmin;
   xh = base[0].xmax;
   yl = base[0].ymin;
   yh = base[0].ymax;
   mindist = (xh -xl)*(xh -xl) +(yh-yl)*(yh-yl);
   nsrch = 0;
   
   if (indx[v0] == NULL) 
      srchlst[nsrch++] = base;
   else
      srchlst[nsrch++] = indx[v0];
      
   topbox = srchlst[0];	 

   for(;;) {
      while (nsrch > 0) {
         qpt = srchlst[--nsrch];
         if (qpt->xmin > xh || qpt->xmax < xl || qpt->ymin > yh || qpt->ymax < yl) continue;
         if (qpt->num < 0) {
/*				ADD DAUGHTERS TO LIST */
            for(i=0;i<4;++i)
               srchlst[nsrch++] = qpt->dghtr[i];
            continue;
         }
         else {
/*				CHECK POINTS */
            for(i=0;i<qpt->num;++i) {
               if (qpt->node[i] == v0) continue;
               dx = vrtx[v0][0] -vrtx[qpt->node[i]][0];
               dy = vrtx[v0][1] -vrtx[qpt->node[i]][1];
               dist = dx*dx + dy*dy;
               if (dist < mindist) {
                  pt = qpt->node[i];
                  mindist = dist;
                  dist = sqrt(dist);
                  xh = vrtx[v0][0] +dist;
                  xl = vrtx[v0][0] -dist;
                  yh = vrtx[v0][1] +dist;
                  yl = vrtx[v0][1] -dist;
               }
            }
         }
      }
      if (topbox->xmin < xl && topbox->xmax > xh && 
            topbox->ymin < yl && topbox->ymax > yh) return(mindist);
      
      exclude = topbox->pind;
      topbox = topbox->prnt;

/*		CHECK IF WE ARE AT THE TOP LEVEL */
      if (topbox == NULL) return(mindist);
      
      for(i=0;i<4;++i)
         if (i != exclude) srchlst[nsrch++] = topbox->dghtr[i];
   }

   return(mindist);
}

/* FIND CLOSEST POINT TO AN ARBITRARY POINT */
FLT quadtree::nearpt(FLT x, FLT y, int& pt) const {
   int i,nsrch,exclude;
   class quad *topbox;
   FLT dist,dx,dy,xh,xl,yh,yl,mindist;
   class quad *qpt;

   xl = base[0].xmin;
   xh = base[0].xmax;
   yl = base[0].ymin;
   yh = base[0].ymax;
   mindist = (xh -xl)*(xh -xl) +(yh-yl)*(yh-yl);
   
   nsrch = 0;
   srchlst[nsrch++] = base;
   topbox = srchlst[0];	 

   for(;;) {
      while (nsrch > 0) {
         qpt = srchlst[--nsrch];
         if (qpt->xmin > xh || qpt->xmax < xl || qpt->ymin > yh || qpt->ymax < yl) continue;
         if (qpt->num < 0) {
/*				ADD DAUGHTERS TO LIST */
            for(i=0;i<4;++i)
               srchlst[nsrch++] = qpt->dghtr[i];
            continue;
         }
         else {
/*				CHECK POINTS */
            for(i=0;i<qpt->num;++i) {
               dx = x -vrtx[qpt->node[i]][0];
               dy = y -vrtx[qpt->node[i]][1];
               dist = dx*dx + dy*dy;
               if (dist < mindist) {
                  pt = qpt->node[i];
                  mindist = dist;
                  dist = sqrt(mindist);
                  xh = x +dist;
                  xl = x -dist;
                  yh = y +dist;
                  yl = y -dist;
               }
            }
         }
      }
      if (topbox->xmin < xl && topbox->xmax > xh && 
            topbox->ymin < yl && topbox->ymax > yh) return(mindist);
      
      exclude = topbox->pind;
      topbox = topbox->prnt;

/*		CHECK IF WE ARE AT THE TOP LEVEL */
      if (topbox == NULL) return(mindist);
      
      for(i=0;i<4;++i)
         if (i != exclude) srchlst[nsrch++] = topbox->dghtr[i];
   }

   return(mindist);
}

void quadtree::dltpt(int v0) {
   int ind, i;
   class quad *qpt;
   
   qpt = indx[v0];
   
   if (qpt == NULL) {
      printf("error: deleting point which isn't in tree %d\n",v0);
      exit(1);
   }
   
   indx[v0] = NULL;

   ind = -1;
   for(i=0;i<qpt->num;++i) {
      if (qpt->node[i] == v0) {
         ind = i;
         break;
      }
   }
   
   if (ind == -1) {
      printf("Couldn't find node to remove: %d\n",v0);
      exit(1);
   }
   
   for(i=ind+1;i<qpt->num;++i)
      qpt->node[i-1] = qpt->node[i];
   
   --qpt->num;
   
   return;
}

void quadtree::output(char *filename, FILETYPE type) {
   int i,nsrch;
   char fnmapp[100];
   class quad *qpt;
   FILE *out;
   
   strcpy(fnmapp,filename);
   
   switch (type) {
      case(text):
         strcat(fnmapp,".text");
         out = fopen(fnmapp,"w");
         if (out == NULL ) {
            printf("couldn't open output file %s for output\n",fnmapp);
            exit(1);
         }
         nsrch = 0;
         srchlst[nsrch++] = base;
         
         while (nsrch > 0) {
            qpt = srchlst[--nsrch];
            if (qpt->num < 0) {
      /*			ADD DAUGHTERS TO LIST */
               for(i=0;i<4;++i)
                  srchlst[nsrch++] = qpt->dghtr[i];
               continue;
            }
            else {
      /*			PRINT POINTS */
               fprintf(out,"box: %f %f %f %f\n",qpt->xmin,qpt->xmax,qpt->ymin,qpt->ymax);
               for(i=0;i<qpt->num;++i) 
                  fprintf(out,"%d %f %f\n",qpt->node[i],vrtx[qpt->node[i]][0],vrtx[qpt->node[i]][1]);
            }
         }
         break;
      case(tecplot):
         strcat(fnmapp,".dat");
         out = fopen(fnmapp,"w");
         if (out == NULL ) {
            printf("couldn't open output file %s for output\n",fnmapp);
            exit(1);
         }
         nsrch = 0;
         srchlst[nsrch++] = base;
         
         while (nsrch > 0) {
            qpt = srchlst[--nsrch];
            if (qpt->num < 0) {
      /*			ADD DAUGHTERS TO LIST */
               for(i=0;i<4;++i)
                  srchlst[nsrch++] = qpt->dghtr[i];
               continue;
            }
            else {
   /*				PRINT BOX & POINTS */
               if (qpt->num > 0) {
                  fprintf(out,"ZONE N=4, E=1, F=FEPOINT, ET=QUADRILATERAL, C=RED\n");
                  fprintf(out,"%f  %f\n",qpt->xmin,qpt->ymin);
                  fprintf(out,"%f  %f\n",qpt->xmax,qpt->ymin);
                  fprintf(out,"%f  %f\n",qpt->xmax,qpt->ymax);
                  fprintf(out,"%f  %f\n",qpt->xmin,qpt->ymax);
                  fprintf(out,"%d %d %d %d\n",1,2,3,4);
                  
                  fprintf(out,"ZONE I=%d, C=BLACK\n",qpt->num);
                  for(i=0;i<qpt->num;++i) 
                     fprintf(out,"%f  %f\n",vrtx[qpt->node[i]][0],vrtx[qpt->node[i]][1]);
               }
            }
         }
         break;
   }

   return;
}

void quadtree::movept(int from, int to) {
   int i;
   class quad *qpt;
   
   qpt = indx[from];
   
   if (qpt == NULL) {
      printf("error: moving point which isn't in tree %d\n",from);
      exit(1);
   }
   
   if (indx[to] != NULL) {
      printf("error: moving to a point which exists %d\n",to);
      exit(1);
   }
   
   indx[to] = indx[from];
   indx[from] = NULL;

   for(i=0;i<qpt->num;++i) {
      if (qpt->node[i] == from) {
         qpt->node[i] = to;
         break;
      }
   }
   assert(i != qpt->num);

   return;
}

void quadtree::update(int bgn, int end) {
   int i;
   class quad *qpt;
   
   for(i=bgn;i<end;++i) {
      qpt = indx[i];
      if (qpt == NULL) continue;
      if (vrtx[i][0] < qpt->xmin || vrtx[i][0] > qpt->xmax || vrtx[i][1] < qpt->ymin || vrtx[i][1] > qpt->ymax) update(i);
   }
   
   return;
}

void quadtree::update(int v0) {
   int i,nsrch,exclude;
   FLT x,y;
   class quad *topbox;
   class quad *qpt;

  	x = vrtx[v0][0];
   y = vrtx[v0][1];
      
   nsrch = 0;
   srchlst[nsrch++] = indx[v0];
   topbox = srchlst[0];	

   for(;;) {
      while (nsrch > 0) {
         qpt = srchlst[--nsrch];
 
         if (qpt->num < 0) {
/*				ADD DAUGHTERS TO LIST */
            for(i=0;i<4;++i)
               srchlst[nsrch++] = qpt->dghtr[i];
            continue;
         }
         else {
            if (qpt->xmin <= x && qpt->xmax >= x && qpt->ymin <= y && qpt->ymax >= y) {
               dltpt(v0);
               addpt(v0, qpt);
               return;
            }
         }
      }
      
      exclude = topbox->pind;
      topbox = topbox->prnt;

/*		CHECK IF WE ARE AT THE TOP LEVEL */
      if (topbox == NULL) {
         printf("#Warning: point is outside of quadtree domain\n");
         printf("#(%f,%f) (%f,%f): (%f, %f)\n",base[0].xmin,base[0].ymin,base[0].xmax,base[0].ymax,x,y);
         addpt(v0);
         return;
      }
      
      for(i=0;i<4;++i)
         if (i != exclude) srchlst[nsrch++] = topbox->dghtr[i];
   }

   return;
}

