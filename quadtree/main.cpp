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
   FLT vtest[11][2] = {{0.1,0.1},{0.2,0.4},{.9,0.2},{.125,0.7},{0.5,0.5},{0.85,0.15},{0.25,0.85},{0.99,0.99},{0.001,.3},{0.01,0.3},{0.2,0.1}};
   FLT dist;
   
   testtree.init(vtest, 10, 0.0, 0.0, 1.0, 1.0);
   for(i=0;i<11;++i) {
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
   
   testtree.reinit();
   
   dist = testtree.nearpt(0, i); 
   printf("%d %f\n",i,dist);
   
   testtree.output("out");
   
   
   return(0);
}