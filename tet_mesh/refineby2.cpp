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
#include <blitz/tinyvec-et.h>

void tet_mesh::refineby2(const class tet_mesh& inmesh) {
    int i,j,n,sind,ind,find,tind,p0,p1,count,pnear,err,initialsidenumber;
    bool found;
    TinyVector<FLT,ND> xpt;
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


		/* five tet */
		tet(++ind).pnt(0)=ijind[0][0][0],tet(ind).pnt(1)=ijind[0][0][1],tet(ind).pnt(2)=ijind[1][0][0],tet(ind).pnt(3)=ijind[0][1][0];
		tet(++ind).pnt(0)=ijind[0][0][1],tet(ind).pnt(1)=ijind[1][0][1],tet(ind).pnt(2)=ijind[1][0][0],tet(ind).pnt(3)=ijind[0][1][0];
		tet(++ind).pnt(0)=ijind[0][0][1],tet(ind).pnt(1)=ijind[1][0][1],tet(ind).pnt(2)=ijind[0][1][0],tet(ind).pnt(3)=ijind[0][1][1];
		tet(++ind).pnt(0)=ijind[1][0][1],tet(ind).pnt(1)=ijind[0][1][0],tet(ind).pnt(2)=ijind[1][1][0],tet(ind).pnt(3)=ijind[1][0][0];
		tet(++ind).pnt(0)=ijind[1][0][1],tet(ind).pnt(1)=ijind[0][1][0],tet(ind).pnt(2)=ijind[0][1][1],tet(ind).pnt(3)=ijind[1][1][0];

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
	
					
	for(int i = 0; i < nebd; ++i){
		ind=ebdry(i)->nseg;
		for(int j = 0; j < ebdry(i)->nseg; ++j){
			sind = ebdry(i)->seg(j).gindx;

      }
		ebdry(i)->setup_next_prev();
		ebdry(i)->reorder();
	}
	
	npnt=count;
	ntet=ind;
	
	/* UPDATE DATA STRUCTURES */
	fixvertexinfo();
	createedgeinfo();
	createfaceinfo();
	morefaceinfo();
	createtetinfo();
	vertexnnbor();	
	for(int i = 0; i < nfbd; ++i){
		fbdry(i)->gbltolclvrtx();
		fbdry(i)->createsideinfo(); 
		fbdry(i)->createttri();
		fbdry(i)->gbltolcltri();
		fbdry(i)->gbltolclside();        
		fbdry(i)->cnt_nbor();
		fbdry(i)->createvtri();    
	}
			
    
    return;
}
