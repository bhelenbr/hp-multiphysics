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

//	/* TO TEST FINDBDRYPT */
//  Using gmsh mesh spheres_mod.msh
//	int i;
//	for (i=0;i<nfbd;++i) 
//		if (fbdry(i)->idnum == 48) 
//			break;
//			
//	for (FLT phi=0.0;phi<M_PI/2;phi+=M_PI/40.0) {
//		FLT rad = 1.5;
//		int tind;
//		FLT r,s;
//		FLT angle = M_PI/3.114;
//		TinyVector<FLT,3> xpt(rad*sin(phi)*cos(angle),rad*sin(phi)*sin(angle),rad*cos(phi));
//		std::cout << xpt << ' ';
//		fbdry(i)->findbdrypt(xpt, tind, r, s);
//		*gbl->log << xpt << ' ' << phi << ' ' << tind << ' ' << r << ' ' << s << std::endl;
//		if (tind > -1) {
//			fbdry(i)->mvpttobdry(tind, r, s, xpt);	
//			std::cout << xpt << '\n';
//		}
//	}
//	exit(1);


void face_bdry::findbdrypt(const TinyVector<FLT,tet_mesh::ND> xpt, int &facloc, FLT &r, FLT &s) {
	FLT normdist, minnormdist = 1.0e99;
	facloc = -1;

	/* For now I am going to search the whole boundary for the best possible answer */
	for (int i=0;i<npnt;++i)
		pnt(i).info = 0;
		
	for (int i=0;i<nseg;++i)
		seg(i).info = 0;
		
		
	for (int tind=0;tind<ntri;++tind) {	
		TinyVector<int,3> v;
		for (int i=0;i<3;++i) {
			v(i) = pnt(tri(tind).pnt(i)).gindx;
			if (tri(tind).tri(i) < 0) ++seg(tri(tind).seg(i)).info;
		}
			
		TinyVector<FLT,3> dx0 = x.pnts(v(2))-x.pnts(v(1));
		TinyVector<FLT,3> dx1 = x.pnts(v(0))-x.pnts(v(2));
		TinyVector<FLT,3> dx2 = x.pnts(v(1))-x.pnts(v(0));
		
		/* Check Position along side */
		TinyMatrix<FLT,3,2> psi;
		psi(0,0) = dot(dx0,xpt-x.pnts(v(1)));
		psi(0,1) = dot(dx0,xpt-x.pnts(v(2)));
		psi(1,0) = dot(dx1,xpt-x.pnts(v(2)));
		psi(1,1) = dot(dx1,xpt-x.pnts(v(0)));
		psi(2,0) = dot(dx2,xpt-x.pnts(v(0)));
		psi(2,1) = dot(dx2,xpt-x.pnts(v(1)));	
		
		if (psi(1,1) > -EPSILON && psi(2,0) < EPSILON) {
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
				}
			}
		}
		else if (psi(2,1) > -EPSILON && psi(0,0) < EPSILON) {
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
				}
			}
		}
		else if (psi(0,1) > -EPSILON && psi(1,0) < EPSILON) {
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
				}
			}
		}

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
			
		if (psi(0,0)*psi(0,1) < EPSILON && tri_wgt(0) < EPSILON) {
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
					s = -1.0;
				}
			}
		}
		else if (psi(1,0)*psi(1,1) < EPSILON && tri_wgt(1) < EPSILON) {
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
				}
			}
		}
		else if (psi(2,0)*psi(2,1) < EPSILON && tri_wgt(2) < EPSILON) {
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
					s = 1. -2.*psi;
					r = -1.0;
				}
			}
		}
			
		if ((fabs(tri_wgt(0)) +fabs(tri_wgt(1)) +fabs(tri_wgt(2)) -(tri_wgt(0) +tri_wgt(1) +tri_wgt(2))) < EPSILON) {
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
			}
		}
	}
	
	if (facloc < 0) {
		*x.gbl->log << idprefix << " couldn't find " << xpt << std::endl;
		x.output(idprefix +"error");
		for (int i=0;i<npnt;++i)
			std::cout << i << ' ' << pnt(i).gindx << ' ' << pnt(i).info << ' ' << pnt(i).nnbor << std::endl;
		for (int i=0;i<nseg;++i)
			std::cout << i << ' ' << seg(i).gindx << ' ' << seg(i).info << std::endl;
	
		exit(1);
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


