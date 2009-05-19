/*
 *  tobasis.cpp
 *  planar++
 *
 *  Created by helenbrk on Fri Oct 12 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include"tet_hp.h"
#include<myblas.h>

 void tet_hp::tobasis(init_bdry_cndtn *ibc, int tlvl) {
   int tind,i,j,k,m,n,v0,v1,eind,find;
   TinyVector<FLT,3> pt;
   int stridey = MXGP;
   int stridex = MXGP*MXGP; 
      
        
   /* LOOP THROUGH VERTICES */
   for(i=0;i<npnt;++i)
      for(n=0;n<NV;++n)
         ugbd(tlvl).v(i,n) = ibc->f(n,pnts(i),gbl->time);

   if (basis::tet(log2p).em == 0) return;

   /* LOOP THROUGH EDGES */   
	for(eind = 0; eind < nseg; ++eind) {
            
		v0 = seg(eind).pnt(0);
		v1 = seg(eind).pnt(1);
      
		if (seg(eind).info < 0) {
			for(n=0;n<ND;++n)
				basis::tet(log2p).proj1d(pnts(v0)(n),pnts(v1)(n),&crd1d(n)(0));
		}
		else {
			crdtocht1d(eind,tlvl);
			for(n=0;n<ND;++n)
				basis::tet(log2p).proj1d(&cht(n)(0),&crd1d(n)(0));
		}
      
		for(n=0;n<NV;++n)
			basis::tet(log2p).proj1d(ugbd(tlvl).v(v0,n),ugbd(tlvl).v(v1,n),&u1d(n)(0));

		for(i=0;i<basis::tet(log2p).gpx; ++i) {
			pt(0) = crd1d(0)(i);
			pt(1) = crd1d(1)(i);
			pt(2) = crd1d(2)(i);
			for(n=0;n<NV;++n){
				// difference between actual function and linear 
				u1d(n)(i) -= ibc->f(n,pt,gbl->time);
			}
		}
            
		for(n=0;n<NV;++n)
			basis::tet(log2p).intgrt1d(&lf(n)(0),&u1d(n)(0));
		 
		for(n=0;n<NV;++n) {
			for(m=0;m<basis::tet(log2p).em;++m){ 
				ugbd(tlvl).e(eind,m,n) = -lf(n)(2+m)*basis::tet(log2p).diag1d(m);
			}
		}
   }
   
   if (basis::tet(log2p).fm == 0) return;
   
   /* LOOP THROUGH FACES */
   for(find = 0; find < ntri; ++find) {
		ugtouht2d_bdry(find,tlvl);
		for(n=0;n<NV;++n)
			basis::tet(log2p).proj2d_bdry(&uht(n)(0),&u2d(n)(0)(0),stridey);

		if (tri(find).info < 0) {
			for(n=0;n<ND;++n)
				basis::tet(log2p).proj2d(vrtxbd(tlvl)(tri(find).pnt(0))(n),vrtxbd(tlvl)(tri(find).pnt(1))(n),vrtxbd(tlvl)(tri(find).pnt(2))(n),&crd2d(n)(0)(0),stridey);
				
		}
		else {
			crdtocht2d(find,tlvl);
			for(n=0;n<ND;++n)
				basis::tet(log2p).proj2d_bdry(&cht(n)(0),&crd2d(n)(0)(0),stridey);
		}
         
      for (i=0; i < basis::tet(log2p).gpx; ++i ) {
         for (j=0; j < basis::tet(log2p).gpy; ++j ) {
            pt(0) = crd2d(0)(i)(j);
            pt(1) = crd2d(1)(i)(j);
            pt(2) = crd2d(2)(i)(j);
            for(n=0;n<NV;++n) {
               u2d(n)(i)(j) -= ibc->f(n,pt,gbl->time);
            }
         }
      }
	  
                     
      for(n=0;n<NV;++n) {
         basis::tet(log2p).intgrt2d(&lf(n)(0),&u2d(n)(0)(0),stridey);
         for(i=0;i<basis::tet(log2p).fm;++i){
            ugbd(tlvl).f(find,i,n) = -lf(n)(3+3*basis::tet(log2p).em+i)*basis::tet(log2p).diag2d(i);
         }
      }
   }
  
   if (basis::tet(log2p).im == 0) return;
   
   /* LOOP THROUGH TETS */
   for(tind = 0; tind < ntet; ++tind) {
		ugtouht_bdry(tind,tlvl);
		for(n=0;n<NV;++n)
         basis::tet(log2p).proj_bdry(&uht(n)(0),&u(n)(0)(0)(0),stridex,stridey);
         
		if (tet(tind).info < 0) {
			for(n=0;n<ND;++n)
				basis::tet(log2p).proj(vrtxbd(tlvl)(tet(tind).pnt(0))(n),vrtxbd(tlvl)(tet(tind).pnt(1))(n),vrtxbd(tlvl)(tet(tind).pnt(2))(n),vrtxbd(tlvl)(tet(tind).pnt(3))(n),&crd(n)(0)(0)(0),stridex,stridey);
		}
		else {
			crdtocht(tind,tlvl);
			for(n=0;n<ND;++n)
				basis::tet(log2p).proj_bdry(&cht(n)(0),&crd(n)(0)(0)(0),stridex,stridey);
		}
         
      for (i=0; i < basis::tet(log2p).gpx; ++i ) {
         for (j=0; j < basis::tet(log2p).gpy; ++j ) {
            for (k=0; k < basis::tet(log2p).gpz; ++k ) {
               pt(0) = crd(0)(i)(j)(k);
               pt(1) = crd(1)(i)(j)(k);
               pt(2) = crd(2)(i)(j)(k);
               for(n=0;n<NV;++n)
                  u(n)(i)(j)(k) -= ibc->f(n,pt,gbl->time);
            }
         }
      }
                     
      for(n=0;n<NV;++n) {
         basis::tet(log2p).intgrt(&lf(n)(0),&u(n)(0)(0)(0),stridex,stridey);
         for(i=0;i<basis::tet(log2p).im;++i) {
            ugbd(tlvl).i(tind,i,n) = -lf(n)(basis::tet(log2p).bm+i)*basis::tet(log2p).diag3d(i);
          }
      }
   }
	
   return;
}   

