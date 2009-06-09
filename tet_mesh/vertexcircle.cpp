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
#include <blitz/tinyvec-et.h>



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

void face_bdry::findbdrypt(const TinyVector<FLT,tet_mesh::ND> xpt, int &facloc, FLT &r, FLT &s) {
	FLT normdist, minnormdist = 1.0e99;
	
	/* Need to fix this so it does less repetitive work FIXME */
	for (int i=0;i<npnt;++i)
		pnt(i).info = 0;
		
	for (int i=0;i<nseg;++i)
		seg(i).info = 0;

	/* For now going to search the whole boundary until I fix quadtree FIXME */
	for (int tind=0;tind<ntri;++tind) {
	
		TinyVector<int,3> v;
		for (int i=0;i<3;++i)
			v(i) = pnt(tri(tind).pnt(i)).gindx;
			
		TinyVector<FLT,3> dx0 = x.pnts(v(2))-x.pnts(v(1));
		TinyVector<FLT,3> dx1 = x.pnts(v(0))-x.pnts(v(2));
		TinyVector<FLT,3> dx2 = x.pnts(v(1))-x.pnts(v(0));
		
		TinyVector<FLT,3> norm = cross(dx0,dx2);
		
		/* Calculate contra-variant vectors */
		TinyVector<FLT,3> nx0 = cross(dx0,norm);
		TinyVector<FLT,3> nx1 = cross(dx1,norm);
		TinyVector<FLT,3> nx2 = cross(dx2,norm);
		
		/* Calculate areas */
		TinyVector<FLT,3> tri_wgt;
		tri_wgt(0) = dot(nx0,xpt-x.pnts(v(1)));
		tri_wgt(1) = dot(nx1,xpt-x.pnts(v(2)));
		tri_wgt(2) = dot(nx2,xpt-x.pnts(v(0)));
		int casenum = (tri_wgt(2) > -EPSILON ? 4 : 0) +(tri_wgt(1) > -EPSILON ? 2 : 0) +(tri_wgt(0) > -EPSILON ? 1 : 0);
		
		switch (casenum) {
			case(1): {
				/* Vertex 0 may be possibillity */
				int pnum = tri(tind).pnt(0);
				if ((++pnt(pnum).info) == pnt(pnum).nnbor) {
					/* This vertex has been tagged by all surrounding triangles */
					TinyVector<FLT,tet_mesh::ND> dx = xpt - x.pnts(v(0));
					normdist = sqrt(dot(dx,dx));
					if (normdist < minnormdist) {
						facloc = tind;
						minnormdist = normdist;
						r = -1.0;
						s = 1.0;
						*x.gbl->log << xpt << " case " << casenum << ' ' << tri(facloc).gindx << ' ' << r << ' ' << s << '\n';
					}
				}
				break;
			}
			case(2): {
				/* Vertex 1 may be possibillity */
				int pnum = tri(tind).pnt(1);
				if ((++pnt(pnum).info) == pnt(pnum).nnbor) {
					/* This vertex has been tagged by all surrounding triangles */
					TinyVector<FLT,tet_mesh::ND> dx = xpt - x.pnts(v(1));
					normdist = sqrt(dot(dx,dx));
					if (normdist < minnormdist) {
						facloc = tind;
						minnormdist = normdist;
						r = -1.0;
						s = -1.0;
						*x.gbl->log << xpt << " case " << casenum << ' ' << tri(facloc).gindx << ' ' << r << ' ' << s << '\n';
					}
				}
				break;
			}
			
			case(4): {
				/* Vertex 2 may be possibillity */
				int pnum = tri(tind).pnt(2);
				if ((++pnt(pnum).info) == pnt(pnum).nnbor) {
					/* This vertex has been tagged by all surrounding triangles */
					TinyVector<FLT,tet_mesh::ND> dx = xpt - x.pnts(v(2));
					normdist = sqrt(dot(dx,dx));
					if (normdist < minnormdist) {
						facloc = tind;
						minnormdist = normdist;
						r = 1.0;
						s = -1.0;
						*x.gbl->log << xpt << " case " << casenum << ' ' << tri(facloc).gindx << ' ' << r << ' ' << s << '\n';
					}
				}
				break;
			} 
			
			case(6): {
				/* Side 0 negative, 1 & 2 positive */
				/* Negative edge is possible location */
				int sind = tri(tind).seg(0);
				seg(sind).info++;
				
				if (seg(sind).info == 2) {
					/* side has been possible two times */
					/* Find distance to edge */
					TinyVector<FLT,tet_mesh::ND> dx = xpt - x.pnts(v(1));
					/* Subtract component of vector along edge */
					FLT psi = dot(dx,dx0)/dot(dx0,dx0);
					TinyVector<FLT,tet_mesh::ND> dxn = dx -psi*dx0;
					normdist = sqrt(dot(dxn,dxn));
					if (normdist < minnormdist) {
						facloc = tind;
						minnormdist = normdist;
						r = 2.*psi -1.;
						s = 0.0;
						*x.gbl->log << xpt << " case " << casenum << ' ' << tri(facloc).gindx << ' ' << r << ' ' << s << '\n';
					}
				}
				break;
			}
			
			case(5): {
				/* Side 1 negative, 0 & 2 positive */
				/* Negative edge is possible location */
				int sind = tri(tind).seg(1);
				seg(sind).info++;
				
				if (seg(sind).info == 2) {
					/* side has been possible two times */
					/* Find distance to edge */
					TinyVector<FLT,tet_mesh::ND> dx = xpt - x.pnts(v(2));
					/* Subtract component of vector along edge */
					FLT psi = dot(dx,dx1)/dot(dx1,dx1);
					TinyVector<FLT,tet_mesh::ND> dxn = dx -psi*dx1;
					normdist = sqrt(dot(dxn,dxn));
					if (normdist < minnormdist) {
						facloc = tind;
						minnormdist = normdist;
						s = 2.*psi -1.;
						r = -s;
						*x.gbl->log << xpt << " case " << casenum << ' ' << tri(facloc).gindx << ' ' << r << ' ' << s << '\n';
					}
				}
				break;
			}
			
			
			case(3): {
				/* Side 2 negative, 0 & 1 positive */
				/* Negative edge is possible location */
				int sind = tri(tind).seg(2);
				seg(sind).info++;
				
				if (seg(sind).info == 2) {
					/* side has been possible two times */
					/* Find distance to edge */
					TinyVector<FLT,tet_mesh::ND> dx = xpt - x.pnts(v(0));
					/* Subtract component of vector along edge */
					FLT psi = dot(dx,dx2)/dot(dx2,dx2);
					TinyVector<FLT,tet_mesh::ND> dxn = dx -psi*dx2;
					normdist = sqrt(dot(dx,dx));
					if (normdist < minnormdist) {
						facloc = tind;
						minnormdist = normdist;
						r = 2.*psi -1.;
						s = 0.0;
						*x.gbl->log << xpt << " case " << casenum << ' ' << tri(facloc).gindx << ' ' << r << ' ' << s << '\n';
					}
				}
				break;
			}
			case(7): {
				/* All positive */
				/* Inside triangle */
				norm /= sqrt(dot(norm,norm));
				TinyVector<FLT,3> dx = xpt - x.pnts(v(1));
				normdist = fabs(dot(dx,norm));
				if (normdist < minnormdist) {
					facloc = tind;
					minnormdist = normdist;
					r = 2.*dot(dx,nx2)/dot(dx0,nx2) -1.;
					s = -2.*dot(dx,nx0)/dot(dx2,nx0) -1.;
//					*x.gbl->log << xpt << " case " << casenum << ' ' << tri(facloc).gindx << ' ' << r << ' ' << s << '\n';
				}
			}
		}
	}
	
	return;
}

double face_bdry::intri(int tind,const TinyVector<FLT,3>& pt) {
	
	TinyVector<int,3> v;
	for (int i=0;i<3;++i)
		v(i) = pnt(tri(tind).pnt(i)).gindx;
		
	TinyVector<FLT,3> dx0 = x.pnts(v(2))-x.pnts(v(1));
	TinyVector<FLT,3> dx1 = x.pnts(v(0))-x.pnts(v(2));
	TinyVector<FLT,3> dx2 = x.pnts(v(1))-x.pnts(v(0));
	
	norm = cross(dx2,dx0);
	
	/* Calculate contra-variants vectors */
	dx0 = cross(norm,dx0);
	dx1 = cross(norm,dx1);
	dx2 = cross(norm,dx2);
	
	/* Calculate areas */
	tri_wgt(0) = dot(dx0,pt-x.pnts(v(1)));
	tri_wgt(1) = dot(dx1,pt-x.pnts(v(2)));
	tri_wgt(2) = dot(dx2,pt-x.pnts(v(0)));	
	
	return(fabs(tri_wgt(0)) +fabs(tri_wgt(1)) +fabs(tri_wgt(2)) - (tri_wgt(0) +tri_wgt(1) +tri_wgt(2)));
}

/* RETURNS WEIGHTS FROM INTRI FUNCTION */
void face_bdry::getwgts(TinyVector<FLT,3> &wt) const {
	int i;
	FLT sum;

	sum = tri_wgt(0) +tri_wgt(1) +tri_wgt(2);
	sum *= sum;
	
	for(i=0;i<3;++i)
		wt(i) = tri_wgt(i)/sum;

	return;
}


