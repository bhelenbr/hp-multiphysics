/*
 *  main.cpp
 *  quadtree
 *
 *  Created by helenbrk on Wed Oct 31 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include"quadtree.h"

int main() {
   int i;
   
   class quadtree testtree;
   FLT x1[2] = {0.0,0.0},x2[2] ={1.0,1.0};
   FLT vtest[11][2] = {{0.1,0.1},{0.2,0.4},{.9,0.2},{.125,0.7},{0.51,0.53},{0.85,0.15},{0.25,0.85},{0.99,0.99},{0.001,.3},{0.01,0.3},{0.2,0.1}};
   FLT dist;
   
   testtree.init(vtest, 20, x1, x2);
   printf("%f %f %f %f\n",testtree.xmin(0),testtree.xmax(0),testtree.xmin(1),testtree.xmax(1));
   for(i=0;i<11;++i) {
      printf("%d\n",i);
      testtree.addpt(i);
   }

   quadtree copytree(testtree);
      
   dist = testtree.nearpt(9, i); 
   printf("%d %f\n",i,dist);
   dist = copytree.nearpt(9, i); 
   printf("%d %f\n",i,dist);
   
   
   vtest[0][0] = 0.9;
   vtest[0][1] = 0.9;
   testtree.update(0,10);
   
   printf("update ok\n");
   
   testtree.reinit();
   
   dist = testtree.nearpt(0, i); 
   printf("%d %f\n",i,dist);
   
   testtree.output("out");
   testtree.output("out",quadtree::text);
   
   return(0);
}