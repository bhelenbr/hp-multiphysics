/*
 *  vertexcircle.cpp
 *  mesh3d
 *
 *  Created by Michael Brazell on 8/27/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_mesh.h"
#include <utilities.h>
#include <float.h>
#include <math.h>


void face_bdry::vertexcircle(int vind){
    int i,j,k,ind,tind,tind2;
    int nbor = pnt(vind).nnbor;
            
    ind = 0;
    // known tri connected to vertex
    x.gbl->i2wk(ind) = pnt(vind).tri;
    x.gbl->i1wk(x.gbl->i2wk(ind)) = 0;
    
    for(i = 0; i < nbor; ++i) {    
        tind = x.gbl->i2wk(i);  
        for(j = 0; j < 3; ++j) {            
            tind2 = tri(tind).tri(j);
            if (tind2 == -1)
                goto NEXTSIDE;
            if (x.gbl->i1wk(tind2) < 0) {            
                for(k = 0; k < 3; ++k) {
                    if(tri(tind2).pnt(k) == vind) {
                        x.gbl->i2wk(++ind) = tind2; // connected tri found
                        x.gbl->i1wk(tind2) = 0;
                        goto NEXTSIDE;                            
                    }
                }
            }    
            NEXTSIDE:;
        }
    }
    
    for(i = 0; i < nbor; ++i) {    
        x.gbl->i1wk(x.gbl->i2wk(i))=-1; // reset i1wk to -1
    }
    return;
}








