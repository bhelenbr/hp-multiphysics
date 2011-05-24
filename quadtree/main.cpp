/*
 *  main.cpp
 *  quadtree
 *
 *  Created by helenbrk on Wed Oct 31 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "quadtree.h"
#include <cstdio>
#include <utilities.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <blitz/array.h>

/* TEST OR FINDUNIQUE */

#define TEST

int main() {
	int i,j,n,pt;


#ifdef TEST
	class quadtree<2> testtree,testtree2;
	TinyVector<FLT,2> x1(0.0,0.0),x2(1.0,1.0);
	FLT vtest[17][2] = {{0.1,0.1},{0.2,0.4},{.9,0.2},{.125,0.7},{0.51,0.53},{0.85,0.15},{0.25,0.85},{0.99,0.99},{0.001,.3},{0.01,0.3},{0.2,0.1},
	{2.0,2.0},{0.9,0.9},{0.9,0.9},{0.9,0.9},{0.9,0.9},{0.9,0.9}};
	blitz::Array<blitz::TinyVector<FLT,2>,1> vtest2(17);
	FLT dist;

	for(i=0;i<17;++i) {
		vtest2(i)(0) = vtest[i][0];
		vtest2(i)(1) = vtest[i][1];
	}
	
	testtree.init(vtest2, 2000, x1, x2);
	printf("%f %f %f %f\n",testtree.xmin(0),testtree.xmax(0),testtree.xmin(1),testtree.xmax(1));
	for(i=0;i<11;++i) {
		printf("%d\n",i);
		testtree.addpt(i);
	}
	
//	/* Test what happens for point outside domain */
//	printf("point outside test\n");
//	testtree.addpt(11);
//	
//	/* Test what happens for same point 5 times */
//	printf("5 times test\n");
//	for(i=12;i<17;++i)
//		testtree.addpt(i);	

	quadtree<2> copytree(testtree);
			
	dist = testtree.nearpt(9, i); 
	printf("%d %f\n",i,dist);
	dist = copytree.nearpt(9, i); 
	printf("%d %f\n",i,dist);
	
	vtest2(0)[0] = 0.9;
	vtest2(0)[1] = 0.9;
	testtree.update(0,10);

	printf("update ok\n");
	
	testtree.reinit();
	
	dist = testtree.nearpt(vtest2(0), i); 
	printf("%d %f\n",i,dist);  
	
	dist = testtree.nearpt(0, i); 
	printf("%d %f\n",i,dist);
	
	testtree.output("out");
	testtree.output("out",quadtree<2>::text);
	
	/* What happens if I insert point twice? */
	testtree.addpt(0);
	
	
//	class quadtree<3> octtree;
//	FLT y1[3] = {0.0,0.0,0.0},y2[3] ={1.0,1.0,1.0};
//	FLT wtest[9][3] = {{0.1,0.1,0.1},{0.1,0.1,0.9},{0.1,0.9,0.1},{0.1,0.9,0.9},{0.9,0.1,0.1},{0.9,0.1,0.9},{0.9,0.9,0.1},{0.9,0.9,0.9}, {0.9,0.85,0.9}};
//	octtree.init(wtest, 20, y1, y2);
//	for(i=0;i<9;++i) {
//			printf("%d\n",i);
//			octtree.addpt(i);
//	}
//
//	quadtree<3> copytree3(octtree);
//			
//	dist = octtree.nearpt(8, i); 
//	printf("%d %f\n",i,dist);
//	dist = copytree3.nearpt(8, i); 
//	printf("%d %f\n",i,dist);
//	
//	
//	wtest[0][0] = 0.9;
//	wtest[0][1] = 0.9;
//	wtest[0][2] = 0.87;
//	octtree.update(0,9);
//	
//	printf("update ok\n");
//	
//	octtree.reinit();
//	
//	dist = octtree.nearpt(0, i); 
//	printf("%d %f\n",i,dist);
//	
//	octtree.output("out3");
//	octtree.output("out3",quadtree<3>::text);
//
//	/* What happens if I insert 4 points outside of domain? */
//	quadtree<2> bigtree;
//	bigtree.init(vtest,100000, x1, x2);
//		for(i=0;i<16;++i) {
//				printf("%d\n",i);
//				bigtree.addpt(i);
//		}
//	bigtree.output("huh",quadtree<2>::text);


#endif

#ifdef FINDUNIQUE
	int nvrts = 556732;

	class quadtree<3> dups;
	FLT (*vrts)[3];
	FLT xmax[3],xmin[3];
	FLT dist1;
	int tvrtx[3];
	int *ipoint;

	ipoint = new int[nvrts];
	if (ipoint == 0) printf("shit\n");
	vrts = (FLT (*)[3]) malloc(nvrts*3*sizeof(double));
	if (vrts == NULL) printf("shit\n");

	dups.allocate(vrts,nvrts);
	scanf("%lf %lf %lf\n",&vrts[0][0],&vrts[0][1],&vrts[0][2]);
	
	for(n=0;n<3;++n) {
			xmax[n] = vrts[0][n];
			xmin[n] = vrts[0][n];
	}
	
	for(i=1;i<nvrts;++i) {
			scanf("%lf %lf %lf",&vrts[i][0],&vrts[i][1],&vrts[i][2]);
			for(n=0;n<3;++n) {
					xmax[n] = MAX(vrts[i][n],xmax[n]);
					xmin[n] = MIN(vrts[i][n],xmin[n]);
			}
	}
	
	for(i=0;i<nvrts;++i)
			for(n=0;n<3;++n)
					vrts[i][n] -= xmax[n];

	for(n=0;n<3;++n) {
			xmin[n] -= xmax[n];
			xmax[n] = 0.0;
	}

	dups.init(vrts,nvrts,xmin,xmax);
	dups.addpt(0);
	printf("%17.10e %17.10e %17.10e\n",vrts[0][0],vrts[0][1],vrts[0][2]);
	ipoint[0] = 0;
	int count = 1;

	for(i=0;i<nvrts;i+=4) {
			for(j=0;j<3;++j) {
					dist1 = dups.nearpt(vrts[i+j], pt);
					if (fabs((vrts[i+j][0]-vrts[pt][0])/xmin[0]) > 1.0e-4
							 || fabs((vrts[i+j][1]-vrts[pt][1])/xmin[1]) > 1.0e-4
							 || fabs((vrts[i+j][2]-vrts[pt][2])/xmin[2]) > 1.0e-4) {
							printf("%17.10e %17.10e %17.10e\n",vrts[i+j][0],vrts[i+j][1],vrts[i+j][2]);
							dups.addpt(i+j);
							ipoint[i+j] = count;
							pt = count++;
					}
					else {
							pt = ipoint[pt];
					}
					tvrtx[j] = pt;
			}
			printf("# %d %d %d\n",tvrtx[0]+1,tvrtx[1]+1,tvrtx[2]+1);

	}
#endif


	return(0);
}
