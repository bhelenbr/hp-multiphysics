/*
 *  main.cpp
 *  quadtree
 *
 *  Created by helenbrk on Wed Oct 31 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include"quadtree.h"
#include<cstdio>

int main() {
   int i;
   
   class quadtree<2> testtree;
   FLT x1[2] = {0.0,0.0},x2[2] ={1.0,1.0};
   FLT vtest[11][2] = {{0.1,0.1},{0.2,0.4},{.9,0.2},{.125,0.7},{0.51,0.53},{0.85,0.15},{0.25,0.85},{0.99,0.99},{0.001,.3},{0.01,0.3},{0.2,0.1}};
   FLT dist;
   
   testtree.init(vtest, 20, x1, x2);
   printf("%f %f %f %f\n",testtree.xmin(0),testtree.xmax(0),testtree.xmin(1),testtree.xmax(1));
   for(i=0;i<11;++i) {
      printf("%d\n",i);
      testtree.addpt(i);
   }

   quadtree<2> copytree(testtree);
      
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
   testtree.output("out",quadtree<2>::text);
   
   class quadtree<3> octtree;
   FLT y1[3] = {0.0,0.0,0.0},y2[3] ={1.0,1.0,1.0};
   FLT wtest[9][3] = {{0.1,0.1,0.1},{0.1,0.1,0.9},{0.1,0.9,0.1},{0.1,0.9,0.9},{0.9,0.1,0.1},{0.9,0.1,0.9},{0.9,0.9,0.1},{0.9,0.9,0.9}, {0.9,0.85,0.9}};
   octtree.init(wtest, 20, y1, y2);
   for(i=0;i<9;++i) {
      printf("%d\n",i);
      octtree.addpt(i);
   }

   quadtree<3> copytree3(octtree);
      
   dist = octtree.nearpt(8, i); 
   printf("%d %f\n",i,dist);
   dist = copytree3.nearpt(8, i); 
   printf("%d %f\n",i,dist);
   
   
   wtest[0][0] = 0.9;
   wtest[0][1] = 0.9;
   wtest[0][2] = 0.87;
   octtree.update(0,9);
   
   printf("update ok\n");
   
   octtree.reinit();
   
   dist = octtree.nearpt(0, i); 
   printf("%d %f\n",i,dist);
   
   octtree.output("out3");
   octtree.output("out3",quadtree<3>::text);
   
   return(0);
}