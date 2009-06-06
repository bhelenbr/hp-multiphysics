/*
 *  vertexball.cpp
 *  mesh3d
 *
 *  Created by Michael Brazell on 8/20/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_mesh.h"
#include <utilities.h>
#include <float.h>
#include <math.h>


void tet_mesh::vertexball(int vind){
	int i,j,k,ind,tind,tind2;
	int nbor = pnt(vind).nnbor;        

	ind = 0;
	// known tet connected to vertex
	gbl->i2wk(ind) = pnt(vind).tet;
	gbl->i1wk(gbl->i2wk(ind)) = 0;
	
	for(i = 0; i < nbor; ++i) {    
		tind = gbl->i2wk(i);  
		for(j = 0; j < 4; ++j) {            
			tind2 = tet(tind).tet(j);
			if (tind2 == -1)
				goto NEXTFACE;
			if (gbl->i1wk(tind2) < 0) {            
				for(k = 0; k < 4; ++k) {
					if(tet(tind2).pnt(k) == vind) {
						gbl->i2wk(++ind) = tind2; // connected tet found
						gbl->i1wk(tind2) = 0;
						goto NEXTFACE;                            
					}
				}
			}
			NEXTFACE:;        
		}
	}
	
	for(i = 0; i < nbor; ++i) {    
		gbl->i1wk(gbl->i2wk(i))=-1; // reset i1wk to -1
	}
	return;
}

//void tet_mesh::spokes(int vind){
//    int tind,sind;
//    int ind = 0;
//    
//    for(int i=0; i < pnt(vind).nnbor; ++i){
//        tind = gbl->i2wk(i);
//        for(int j=0; j < 6; ++j){
//            sind=tet(tind).seg(j);
//            if(gbl->i1wk(sind) < 0){
//                for(int k=0; k < 2; ++k){
//                    if(vind == seg(sind).pnt(k)){
//                        gbl->i3wk(ind++)=sind;
//                        gbl->i1wk(sind) = 1;
//                    }
//                }
//            }
//            else {
//                ++gbl->i1wk(sind);
//            }
//        }
//    }
//    
//    nspk=ind;
//    for(int i = 0; i < nspk; ++i) {    
//        gbl->i2wk(i)=gbl->i1wk(gbl->i3wk(i));//number of times edge is hit
//        gbl->i1wk(gbl->i3wk(i))=-1; // reset i1wk to -1
//    }    
//
//    return;
//}

void tet_mesh::ring(int eind){
	int i,j,k,ind,tind,tind2;
	int nbor = seg(eind).nnbor;        

	ind = 0;
	// known tet connected to edge
	gbl->i2wk(ind) = seg(eind).tet;
	gbl->i1wk(gbl->i2wk(ind)) = 0;
	
	for(i = 0; i < nbor; ++i) {    
		tind = gbl->i2wk(i);  
		for(j = 0; j < 4; ++j) {            
			tind2 = tet(tind).tet(j);
			if (tind2 == -1)
				goto NEXTFACE;
			if (gbl->i1wk(tind2) < 0) {            
				for(k = 0; k < 6; ++k) {
					if(tet(tind2).seg(k) == eind) {
						gbl->i2wk(++ind) = tind2; // connected tet found
						gbl->i1wk(tind2) = 0;
						goto NEXTFACE;                            
					}
				}
			}
			NEXTFACE:;        
		}
	}
	
	for(i = 0; i < nbor; ++i) {    
		gbl->i1wk(gbl->i2wk(i))=-1; // reset i1wk to -1
	}
	return;
}




