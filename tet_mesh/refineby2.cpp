/*
 *  refineby2.cpp
 *  mesh
 *
 *  Created by Brian Helenbrook on Tue Jun 11 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */
#include "tet_mesh.h"
#include <utilities.h>
#include <assert.h>
//#include <blitz/tinyvec-et.h>

void tet_mesh::refineby2(const class tet_mesh& inmesh) {
	int n,sind,ind,find,p0,p1;
	TinyVector<FLT,ND> xpt;
	TinyVector<FLT,ND> edge_length;
	int ijind[3][3][3];
	
	/* INPUT MESH MUST HAVE GROWTH FACTOR OF 4 */
	/* BECAUSE OF gbl->intwk USAGE */
	if (!initialized) {
		/* VERTEX STORAGE ALLOCATION */
		init(inmesh,duplicate,0.5);
	}
	
	this->copy(inmesh);
					
	/* CALCULATE LOCATION OF NEW INTERIOR POINTS */
	for(sind=0;sind<nseg;++sind) {            
		p0 = seg(sind).pnt(0);
		p1 = seg(sind).pnt(1);
		
		/* MIDPOINT */
		for(n=0;n<ND;++n)
			xpt(n) = 0.5*(pnts(p0)(n) +pnts(p1)(n));
					
		/* INSERT POINT */
		for(n=0;n<ND;++n)
			pnts(npnt)(n) = xpt(n);
			
		otree.addpt(npnt);
		++npnt;
	}

	/* ADD NEW TETS TO MESH */
	ind=ntet;
	for(int tind=0;tind<ntet;++tind){                
		 /* VERTICES */
		ijind[0][0][2] = tet(tind).pnt(0);
		ijind[0][2][0] = tet(tind).pnt(1);
		ijind[0][0][0] = tet(tind).pnt(2);
		ijind[2][0][0] = tet(tind).pnt(3);

		/* EDGES */
		ijind[1][0][0] = npnt +tet(tind).seg(0);
		ijind[1][1][0] = npnt +tet(tind).seg(1);
		ijind[0][1][0] = npnt +tet(tind).seg(2);
		ijind[0][1][1] = npnt +tet(tind).seg(3);
		ijind[0][0][1] = npnt +tet(tind).seg(4);
		ijind[1][0][1] = npnt +tet(tind).seg(5);

		p0 = npnt +tet(tind).seg(2);
		p1 = npnt +tet(tind).seg(5);
		edge_length(0) = sqrt(pow(pnts(p0)(0)-pnts(p1)(0),2)+pow(pnts(p0)(1)-pnts(p1)(1),2)+pow(pnts(p0)(2)-pnts(p1)(2),2));
		
		p0 = npnt +tet(tind).seg(1);
		p1 = npnt +tet(tind).seg(4);
		edge_length(1) = sqrt(pow(pnts(p0)(0)-pnts(p1)(0),2)+pow(pnts(p0)(1)-pnts(p1)(1),2)+pow(pnts(p0)(2)-pnts(p1)(2),2));
		
		p0 = npnt +tet(tind).seg(0);
		p1 = npnt +tet(tind).seg(3);
		edge_length(2) = sqrt(pow(pnts(p0)(0)-pnts(p1)(0),2)+pow(pnts(p0)(1)-pnts(p1)(1),2)+pow(pnts(p0)(2)-pnts(p1)(2),2));
	
		if(edge_length(0) < edge_length(1)) {
			sind = 0;
		}
		else{
			sind = 1;
		}
		
		if(edge_length(2) < edge_length(sind)) {
			sind = 2;
		}
		
		switch(sind) {
			case(0):
				// repeated seg2 seg5 010 101 
				//(e4,e5,e0,e2)
				tet(++ind).pnt(0)=ijind[0][0][1],tet(ind).pnt(1)=ijind[1][0][1],tet(ind).pnt(2)=ijind[1][0][0],tet(ind).pnt(3)=ijind[0][1][0];
				//(e4,e5,e2,e3)
				tet(++ind).pnt(0)=ijind[0][0][1],tet(ind).pnt(1)=ijind[1][0][1],tet(ind).pnt(2)=ijind[0][1][0],tet(ind).pnt(3)=ijind[0][1][1];
				//(e5,e2,e1,e0)
				tet(++ind).pnt(0)=ijind[1][0][1],tet(ind).pnt(1)=ijind[0][1][0],tet(ind).pnt(2)=ijind[1][1][0],tet(ind).pnt(3)=ijind[1][0][0];
				//(e5,e2,e3,e1)
				tet(++ind).pnt(0)=ijind[1][0][1],tet(ind).pnt(1)=ijind[0][1][0],tet(ind).pnt(2)=ijind[0][1][1],tet(ind).pnt(3)=ijind[1][1][0];
				
			case(1):
				// repeated seg1 seg4 110 001
				//(e4,e5,e0,e1)
				tet(++ind).pnt(0)=ijind[0][0][1],tet(ind).pnt(1)=ijind[1][0][1],tet(ind).pnt(2)=ijind[1][0][0],tet(ind).pnt(3)=ijind[1][1][0];
				//(e4,e5,e1,e3)
				tet(++ind).pnt(0)=ijind[0][0][1],tet(ind).pnt(1)=ijind[1][0][1],tet(ind).pnt(2)=ijind[1][1][0],tet(ind).pnt(3)=ijind[0][1][1];
				//(e5,e4,e1,e0)
				tet(++ind).pnt(0)=ijind[1][0][1],tet(ind).pnt(1)=ijind[0][0][1],tet(ind).pnt(2)=ijind[1][1][0],tet(ind).pnt(3)=ijind[1][0][0];
				//(e5,e4,e3,e1)
				tet(++ind).pnt(0)=ijind[1][0][1],tet(ind).pnt(1)=ijind[0][0][1],tet(ind).pnt(2)=ijind[0][1][1],tet(ind).pnt(3)=ijind[1][1][0];
				
			case(2):
				// repeated seg0 seg3 100 011
				//(e4,e5,e0,e3)
				tet(++ind).pnt(0)=ijind[0][0][1],tet(ind).pnt(1)=ijind[1][0][1],tet(ind).pnt(2)=ijind[1][0][0],tet(ind).pnt(3)=ijind[0][1][1];
				//(e4,e5,e0,e3)
				tet(++ind).pnt(0)=ijind[0][0][1],tet(ind).pnt(1)=ijind[1][0][1],tet(ind).pnt(2)=ijind[1][0][0],tet(ind).pnt(3)=ijind[0][1][1];
				//(e5,e3,e1,e0)
				tet(++ind).pnt(0)=ijind[1][0][1],tet(ind).pnt(1)=ijind[0][1][1],tet(ind).pnt(2)=ijind[1][1][0],tet(ind).pnt(3)=ijind[1][0][0];
				//(e5,e0,e3,e1)
				tet(++ind).pnt(0)=ijind[1][0][1],tet(ind).pnt(1)=ijind[1][0][0],tet(ind).pnt(2)=ijind[0][1][1],tet(ind).pnt(3)=ijind[1][1][0];
				
		}
		

		
		/* single tet */
		tet(++ind).pnt(0)=ijind[0][0][0],tet(ind).pnt(1)=ijind[0][0][1],tet(ind).pnt(2)=ijind[1][0][0],tet(ind).pnt(3)=ijind[0][1][0];
	
		/* single tet */
		tet(++ind).pnt(0)=ijind[0][0][1],tet(ind).pnt(1)=ijind[1][0][1],tet(ind).pnt(2)=ijind[0][1][1],tet(ind).pnt(3)=ijind[0][0][2];
		
		/* single tet */
		tet(++ind).pnt(0)=ijind[0][1][0],tet(ind).pnt(1)=ijind[1][1][0],tet(ind).pnt(2)= ijind[0][2][0],tet(ind).pnt(3)=ijind[0][1][1];

		/* single tet */
		tet(tind).pnt(0)=ijind[1][0][0],tet(tind).pnt(1)=ijind[2][0][0],tet(tind).pnt(2)=ijind[1][1][0],tet(tind).pnt(3)=ijind[1][0][1];
	}
	
	
		/* REFINE BY 2 FACE BOUNDARIES */
	for(int i = 0; i < nfbd; ++i){
		ind=fbdry(i)->ntri;
		for(int j = 0; j < fbdry(i)->ntri; ++j){
			find = fbdry(i)->tri(j).gindx;
			fbdry(i)->tri(j).pnt(0)=npnt+tri(find).seg(0),fbdry(i)->tri(j).pnt(1)=npnt+tri(find).seg(1),fbdry(i)->tri(j).pnt(2)=npnt+tri(find).seg(2);
			fbdry(i)->tri(++ind).pnt(0)=tri(find).pnt(0),fbdry(i)->tri(ind).pnt(1)=npnt+tri(find).seg(2),fbdry(i)->tri(ind).pnt(2)=npnt+tri(find).seg(1);
			fbdry(i)->tri(++ind).pnt(0)=npnt+tri(find).seg(2),fbdry(i)->tri(ind).pnt(1)=tri(find).pnt(1),fbdry(i)->tri(ind).pnt(2)=npnt+tri(find).seg(0);
			fbdry(i)->tri(++ind).pnt(0)=npnt+tri(find).seg(1),fbdry(i)->tri(ind).pnt(1)=npnt+tri(find).seg(0),fbdry(i)->tri(ind).pnt(2)=tri(find).pnt(2);
		}
	}
	
					
	for(int i = 0; i < nebd; ++i) {
		ind=ebdry(i)->nseg;
		for(int j = 0; j < ebdry(i)->nseg; ++j) {
			sind = ebdry(i)->seg(j).gindx;
		}
		ebdry(i)->setup_next_prev();
		ebdry(i)->reorder();
	}
	
	ntet=ind;
	
	/* UPDATE DATA STRUCTURES */
	reorient_tets();
	create_from_tet();
	
	for(int i = 0; i < nfbd; ++i){
		fbdry(i)->create_from_tri();
	}
	
	return;
}

