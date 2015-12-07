/*
 *  setupbdry.cpp
 *  mesh3d
 *
 *  Created by Michael Brazell on 8/8/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#include "tet_mesh.h"
#include <float.h>

/* this assumes that the global vertex numbers are in tri.pnt */
void face_bdry::load_gbl_tri_pnt_from_mesh() {

	for(int i = 0; i < ntri; ++i) {
		int tind = tri(i).gindx;
		for(int j = 0; j < 3; ++j) {
			tri(i).pnt(j) = x.tri(tind).pnt(j);
		}
	}

	return;
}

void face_bdry::create_from_gbl_tri_pnt() {
	
	/* Create pnt.gindx and tri.pnt */
	npnt = 0;
	for(int i = 0; i < ntri; ++i) {
		for(int j = 0; j < 3; ++j) {
			int v0 = tri(i).pnt(j);
			if (x.gbl->i1wk(v0) < 0) {
				x.gbl->i1wk(v0) = npnt;
				tri(i).pnt(j) = npnt;
				pnt(npnt).gindx = v0;
				++npnt;
			}
			else {
				tri(i).pnt(j) = x.gbl->i1wk(v0);
			}
		}
	}
	
	for(int i = 0; i < npnt; ++i) {
		x.gbl->i1wk(pnt(i).gindx) = -1;
	}
	
	/* Fill in gbl indices of triangles */
	/* Triangles are oriented the same as in the mesh */
	create_tri_gindx();
	
	/* Create seg definitions */
	int v1,v2,vout;
	
	for(int i=0;i<npnt;++i)
		pnt(i).info = -1;
	
	nseg = 0;
	for(int tind=0;tind<ntri;++tind) {
		int gtind = tri(tind).gindx;
		tri(tind).sgn = x.tri(gtind).sgn;
		vout = tri(tind).pnt(0);
		v1 = tri(tind).pnt(1);
		v2 = tri(tind).pnt(2);
		for(int j=0;j<3;++j) {
			assert(x.tri(gtind).pnt(j) == pnt(tri(tind).pnt(j)).gindx);
			int gsind = x.tri(gtind).seg(j);
			if (x.gbl->i1wk(gsind) == -1) {
				tri(tind).seg(j) = nseg;
				seg(nseg).gindx = gsind;
				int sign = tri(tind).sgn(j);
				seg(nseg).pnt((1-sign)/2) = v1;
				seg(nseg).pnt((1+sign)/2) = v2;
				seg(nseg).tri((1-sign)/2) = tind;
				seg(nseg).tri((1+sign)/2) = -1;
				x.gbl->i1wk(gsind) = nseg++;
			}
			else {
				int sind = x.gbl->i1wk(gsind);
				tri(tind).seg(j) = sind;
				int sign = tri(tind).sgn(j);
				
				if (seg(sind).tri((1-sign)/2) == -1)
					seg(sind).tri((1-sign)/2) = tind;
				else {
					*x.gbl->log << "different orientation triangles on face boundary? " << tind << ' ' << gtind << std::endl;
					sim::abort(__LINE__, __FILE__, x.gbl->log);
				}
			}
			int temp = vout;
			vout = v1;
			v1 = v2;
			v2 = temp;
		}
	}
	
	for(int sind=0;sind<nseg;++sind) {
		x.gbl->i1wk(seg(sind).gindx) = -1;
	}
	create_pntnnbor_tritri_pnttri();

	return;
}

void face_bdry::create_seg_gindx() {
	int i,j,lcl2,eind,sind;    
	long lcl0,lcl1;
	TinyVector<int,2> v,a;
	
	for(sind=0;sind<nseg; ++sind) {
		v(0)=pnt(seg(sind).pnt(0)).gindx;
		v(1)=pnt(seg(sind).pnt(1)).gindx;
		lcl0=v(0)+v(1);
		lcl1=(v(0)+1)*(v(1)+1);
		x.vertexball(v(0));
		for(i=0; i < x.pnt(v(0)).nnbor; ++i) {
			for(j = 0; j < 6; ++j) {
				eind=x.tet(x.gbl->i2wk(i)).seg(j);        
				a(0)=x.seg(eind).pnt(0);    
				a(1)=x.seg(eind).pnt(1);
				lcl2=abs(lcl0-a(0)-a(1));
				lcl2+=abs(lcl1-(a(0)+1)*(a(1)+1));
				if (lcl2 == 0) {
					seg(sind).gindx=eind;
					goto NEXTSIDE;
				}
			}        
		}
		*x.gbl->log << "Trouble matching face boundary seg to global side definitions " << idprefix << ' ' << sind << ' ' << seg(sind).pnt(0) << ' ' << seg(sind).pnt(1) << ' ' << v(0) << ' ' << v(1) << std::endl;
NEXTSIDE:;
	}
	
	return;
}

void face_bdry::create_tri_gindx() {
	int i,j,tind,lcl2,find;    
	long lcl0,lcl1;
	TinyVector<int,3> v,a;
	
	for(tind=0;tind<ntri; ++tind) {
		v(0)=pnt(tri(tind).pnt(0)).gindx;
		v(1)=pnt(tri(tind).pnt(1)).gindx;
		v(2)=pnt(tri(tind).pnt(2)).gindx;
		lcl0=v(0)+v(1)+v(2);
		lcl1=(v(0)+1)*(v(1)+1)*(v(2)+1);
		x.vertexball(v(0));
		for(i = 0; i < x.pnt(v(0)).nnbor; ++i) {
			for(j = 0; j < 4; ++j) {
				find=x.tet(x.gbl->i2wk(i)).tri(j);        
				a(0)=x.tri(find).pnt(0);    
				a(1)=x.tri(find).pnt(1);
				a(2)=x.tri(find).pnt(2);
				lcl2=abs(lcl0-a(0)-a(1)-a(2));
				lcl2+=abs(lcl1-(a(0)+1)*(a(1)+1)*(a(2)+1));
				if (lcl2 == 0) {
					tri(tind).gindx=find;
					goto NEXTTRI;
				}
			}
		}
		*x.gbl->log << "Trouble matching face boundary tri to global tri definitions " << idprefix << ' ' << tind << ' ' << v(0) << ' ' << v(1) << ' ' << v(2) << std::endl;
NEXTTRI:;
		if (v(0) == a(0) && v(1) == a(1) && v(2)==a(2)) continue;
		
		*x.gbl->log << "Rotation trouble matching face boundary tri to global tri definitions " << idprefix << ' ' << tind << ' ' << v(0) << ' ' << v(1) << ' ' << v(2) << ' ' << a(0) << ' ' << a(1) << ' ' << a(2) << std::endl;
	}
	
	return;
}

/* fills in all info after loading minimal mesh data (grid) */
void face_bdry::create_from_gindx() {
	int ind;
	
	/* fill in pnt data */
	for(int i = 0; i < npnt; ++i) {
		x.gbl->i1wk(pnt(i).gindx)=i;
	}
	
	for(int i = 0; i < nseg; ++i) {
		ind = seg(i).gindx;
		seg(i).pnt(0)=x.gbl->i1wk(x.seg(ind).pnt(0));
		seg(i).pnt(1)=x.gbl->i1wk(x.seg(ind).pnt(1));
		seg(i).tri = -1;
	}
	
	for(int i = 0; i < ntri; ++i) {
		ind = tri(i).gindx;
		tri(i).pnt(0)=x.gbl->i1wk(x.tri(ind).pnt(0));
		tri(i).pnt(1)=x.gbl->i1wk(x.tri(ind).pnt(1));
		tri(i).pnt(2)=x.gbl->i1wk(x.tri(ind).pnt(2));
	}
	
	/* Reset i1wk */
	for(int i = 0; i < npnt; ++i) {
		x.gbl->i1wk(pnt(i).gindx) = -1;
	}

	
	/* Fill in tri.seg and seg.tri data */
	for(int i = 0; i < nseg; ++i) {
		x.gbl->i1wk(seg(i).gindx)=i;
	}
	
	for(int i = 0; i < ntri; ++i) {
		int tind = tri(i).gindx;
		tri(i).sgn = x.tri(tind).sgn;
		for (int j=0;j<3;++j) {
			int sind = x.gbl->i1wk(x.tri(tind).seg(j));
			tri(i).seg(j)=sind;
			int sign = tri(i).sgn(j);
			if (seg(sind).tri((1-sign)/2) == -1)
				seg(sind).tri((1-sign)/2) = i;
			else {
				*x.gbl->log << "tri's not consistently defined cw or ccw on face " << std::endl;
				sim::abort(__LINE__, __FILE__, x.gbl->log);
			}
		}
	}
	for(int i = 0; i < nseg; ++i) {
		x.gbl->i1wk(seg(i).gindx) = -1;
	}

	create_pntnnbor_tritri_pnttri();

	return;
}

void face_bdry::create_pntnnbor_tritri_pnttri() {
	
	/* pnt.nnbor */
	for (int i=0;i<npnt;++i)
		pnt(i).nnbor = 0;
	
	for(int i=0;i<nseg;++i) {
		++pnt(seg(i).pnt(0)).nnbor;
		++pnt(seg(i).pnt(1)).nnbor;
	}

	/* tri.tri */
	for(int tind=0;tind<ntri;++tind) {
		for(int j=0;j<3;++j) {
			int sind = tri(tind).seg(j);
			int flip = (1 +tri(tind).sgn(j))/2;
			tri(tind).tri(j) = seg(sind).tri(flip);
		}
	}
	
	/* pnt.tri */
	for (int tind=0;tind<ntri;++tind)
		for(int i=0;i<3;++i)
			pnt(tri(tind).pnt(i)).tri = tind;
	
	return;
}

/* An octree which contains only the face points */
void face_bdry::treeinit() {
	int i,n,v0;
	FLT x1[tet_mesh::ND], x2[tet_mesh::ND], dx;
	
	for(n=0;n<tet_mesh::ND;++n)    {
		x1[n] = x.pnts(pnt(0).gindx)(n);
		x2[n] = x.pnts(pnt(0).gindx)(n);
	}

	for (i=0;i<npnt;++i) {
		v0 = pnt(i).gindx;
		//cout << "idnum" << idnum <<" v0 in treeinit "<<v0 << ' '<< i<< '/'<< npnt<<  endl;
		for(n=0;n<tet_mesh::ND;++n) {
			x1[n] = MIN(x1[n],x.pnts(v0)(n));
			x2[n] = MAX(x2[n],x.pnts(v0)(n));
		}
	}
	
	for(n=0;n<tet_mesh::ND;++n) {
		dx = MAX(x2[n]-x1[n],100.0*EPSILON);
		x1[n] -= 0.25*dx;
		x2[n] += 0.25*dx;
	}
	
	//treeinit(x1,x2);

	return;
}

void face_bdry::treeinit(FLT x1[tet_mesh::ND], FLT x2[tet_mesh::ND]) {
	
	otree.init(x1,x2);
	
	for(int i=0;i<npnt;++i) 
		otree.addpt(pnt(i).gindx);
	
	return;
}

void face_bdry::checkintegrity() {
	int i,j,sind,dir;
	bool abort = false;
	
	for(i=0;i<npnt;++i) {
		int tind = pnt(i).tri;
		for(j=0;j<3;++j)
			if (tri(tind).pnt(j) == i) goto next;
		
		*x.gbl->log << "tri.pnt is out of whack for " << idprefix << ' ' <<  i << tind << std::endl;
		sim::abort(__LINE__,__FILE__,x.gbl->log);
		
	next: continue;
	}
	
	for(i=0;i<nseg;++i) {
		for(j=0;j<2;++j) {
			if (pnt(seg(i).pnt(j)).gindx != x.seg(seg(i).gindx).pnt(j)) {
				*x.gbl->log << "failed segment gindx check for " << idprefix << " sind " << i << ' ' << pnt(seg(i).pnt(j)).gindx << ' ' << x.seg(seg(i).gindx).pnt(j) << std::endl;
				*x.gbl->log << "seg(i).gindx " << seg(i).gindx << " seg(i).pnt " << seg(i).pnt << std::endl;
				abort = true;
			}
		}
	}

	
	for(i=0;i<ntri;++i) {
		
		//if (tri(i).info < 0) continue;
		//if (area(i) < 0.0) *x.gbl->log << "negative area" << i << std::endl;
		
		for(j=0;j<3;++j) {
			sind = tri(i).seg(j);
			dir = (1-tri(i).sgn(j))/2;
			
			if (seg(sind).pnt(dir) != tri(i).pnt((j+1)%3) && seg(sind).pnt(1-dir) != tri(i).pnt((j+2)%3)) {
				*x.gbl->log << "failed pnt check for " << idprefix << " tind " << i << " sind " << sind << std::endl;
				*x.gbl->log << seg(sind).pnt << ' ' << tri(i).pnt << std::endl;
				abort = true;
			}
			
			if (seg(sind).tri(dir) != i) {
				*x.gbl->log << "failed segment check for " << idprefix << " tind " << i << " sind " << sind << ' ' << seg(sind).tri << ' ' << dir << std::endl;
				abort = true;
			}
			
			if (tri(i).tri(j) != seg(sind).tri(1-dir)) {
				*x.gbl->log << "failed ttri check for " << idprefix << " tind " << i << " sind " << sind << std::endl;
				abort = true;
			}
			
			if (pnt(tri(i).pnt(j)).gindx != x.tri(tri(i).gindx).pnt(j)) {
				*x.gbl->log << "failed tri gindx check for " << idprefix << " sind " << i << ' ' << pnt(tri(i).pnt(j)).gindx << ' ' << x.tri(tri(i).gindx).pnt(j) << std::endl;
				abort = true;
			}
		}
	}
	

	
	return;
}

void edge_bdry::checkintegrity() {
	bool abort = false;
	
	for(int i=0;i<nseg;++i) {
		if (seg(i).next > -1) {
			if (seg(i).next != i+1) {
				*x.gbl->log << "sides are out of order for " << idprefix << ' ' << i << ' ' << seg(i).next << std::endl;
				abort = true;
			}
			if (seg(seg(i).next).prev != i) {
				*x.gbl->log << "next prev is out of whack for " << idprefix << ' ' << i << ' ' << seg(i).next << ' ' << seg(seg(i).next).prev << std::endl;
				abort = true;
			}
			else {
				if (x.seg(seg(i).gindx).pnt(1) != x.seg(seg(seg(i).next).gindx).pnt(0)) {
					*x.gbl->log << "something funny in definition of edge boundary for " << idprefix << ' ' << i << ' ' << x.seg(seg(i).gindx).pnt(1) << ' ' << x.seg(seg(seg(i).next).gindx).pnt(0) << std::endl;
					abort = true;
				}
			}
		}
	}
	
	if (abort)
		sim::abort(__LINE__, __FILE__, x.gbl->log);
	
	return;
}
