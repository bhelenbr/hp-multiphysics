/*
 *  setupbdry.cpp
 *  mesh3d
 *
 *  Created by Michael Brazell on 8/8/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#include "tet_mesh.h"
#include <utilities.h>
#include <float.h>

/* only run once else it gets screwed up */
/* this assumes that the global vertex numbers are in tri.pnt */
void face_bdry::create_gbl_pnt_from_tri(void) {

	for(int i = 0; i < ntri; ++i) {
		int tind = tri(i).gindx;
		for(int j = 0; j < 3; ++j) {
			tri(i).pnt(j) = x.tri(tind).pnt(j);
		}
	}

	return;
}


void face_bdry::create_tri_pnt_and_pnt_gindx_from_gbltris(void) {
	int nvrt,vrt;    
	Array<int,1> vinfo(x.npnt);
	
	vinfo = -1;
		
	nvrt = 0;
	for(int i = 0; i < ntri; ++i) {
		for(int j = 0; j < 3; ++j) {
			vrt = tri(i).pnt(j);
			if (vinfo(vrt) < 0) {
				vinfo(vrt) = nvrt;
				tri(i).pnt(j) = nvrt;
				pnt(nvrt).gindx = vrt;
				++nvrt;                
			}
			else {
				tri(i).pnt(j) = vinfo(vrt);
			}
							
		}		
	}	
	
	// store total number of vertex on face bdry
	npnt = nvrt;
	
	return;
}

void face_bdry::create_seg_gindx(void) {
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

void face_bdry::create_tri_gindx(void) {
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
	}
	
	return;
}

/* fills in all info after loading minimal mesh data (grid) */
void face_bdry::convert_gbl_to_lcl(void) {
	int ind;
		
	for(int i = 0; i < npnt; ++i) {
		x.gbl->i2wk(pnt(i).gindx)=i;
	}
	
	for(int i = 0; i < nseg; ++i) {
		ind = seg(i).gindx;
		seg(i).pnt(0)=x.gbl->i2wk(x.seg(ind).pnt(0));    
		seg(i).pnt(1)=x.gbl->i2wk(x.seg(ind).pnt(1));    
	}
	
	for(int i = 0; i < ntri; ++i) {
		ind = tri(i).gindx;
		tri(i).pnt(0)=x.gbl->i2wk(x.tri(ind).pnt(0));    
		tri(i).pnt(1)=x.gbl->i2wk(x.tri(ind).pnt(1));    
		tri(i).pnt(2)=x.gbl->i2wk(x.tri(ind).pnt(2));
	}
	
	for(int i = 0; i < nseg; ++i) {
		x.gbl->i2wk(seg(i).gindx)=i;
	}
	for(int i = 0; i < ntri; ++i) {
		ind = tri(i).gindx;
		tri(i).seg(0)=x.gbl->i2wk(x.tri(ind).seg(0));
		tri(i).seg(1)=x.gbl->i2wk(x.tri(ind).seg(1));
		tri(i).seg(2)=x.gbl->i2wk(x.tri(ind).seg(2));
		tri(i).sgn(0)=x.tri(ind).sgn(0);
		tri(i).sgn(1)=x.tri(ind).sgn(1);
		tri(i).sgn(2)=x.tri(ind).sgn(2);
	}
	
	for(int i = 0; i < nseg; ++i)
		seg(i).tri = -1;
		
	/* tri's sharing a side */
	for(int i = 0; i < ntri; ++i) {
		for(int j = 0; j < 3; ++j) {
			ind = tri(i).seg(j);
			if (seg(ind).tri(0) < 0) {
				seg(ind).tri(0) = i;
			}
			else {
				seg(ind).tri(1) = i;
			}			
		}
	}
	
	/* 3 tri's connected to each tri */
	for(int i = 0; i < ntri; ++i) {
		for(int j = 0; j < 3; ++j) {
			ind = tri(i).seg(j);
			if (seg(ind).tri(0) == i) {
				tri(i).tri(j) = seg(ind).tri(1);
			}
			else if (seg(ind).tri(1) == i) {
				tri(i).tri(j) = seg(ind).tri(0);
			}
		}
	}

	return;
}


/* CREATE SIDELIST FROM TRIANGLE VERTEX LIST */
/* USES VINFO TO STORE FIRST SIND FROM VERTEX */
/* USES SINFO TO STORE NEXT SIND FROM SIND */
/* TVRTX MUST BE COUNTERCLOCKWISE ORDERED */
void face_bdry::create_seg_from_tri(void) {
	int i,j,tind,v1,v2,vout,temp,minv,maxv,order,sind=-1,sindprev=-1;
	
	for(i=0;i<npnt;++i)
		pnt(i).info = -1;
		
	nseg = 0;
	for(tind=0;tind<ntri;++tind) {
		vout = tri(tind).pnt(0);
		v1 = tri(tind).pnt(1);
		v2 = tri(tind).pnt(2);
		for(j=0;j<3;++j) {
			/* CHECK IF SIDE HAS BEEN CREATED ALREADY */
			if (v2 > v1) {
				minv = v1;
				maxv = v2;
				order = 0;
			}
			else {
				minv = v2;
				maxv = v1;
				order = 1;
			}
			
			sind = pnt(minv).info;
			while (sind >= 0) {
				if (maxv == seg(sind).pnt(order)) {
					if (seg(sind).tri(1) >= 0) {
						*x.gbl->log << "Error: side " << sind << " has been matched with Triangle " << tind << " 3 times" << std::endl;  
						*x.gbl->log << "tets "<< x.tri(tri(tind).gindx).tet(0) << ' ' << x.tri(tri(tind).gindx).tet(1)<< " pnts "<< x.pnts(x.tri(tri(tind).gindx).pnt(0)) << ' ' << x.pnts(x.tri(tri(tind).gindx).pnt(1)) << ' ' << x.pnts(x.tri(tri(tind).gindx).pnt(2))<<std::endl;  

						x.output("error",tet_mesh::easymesh);
						x.output("error",tet_mesh::grid);
						exit(1);
					}
					else {
						seg(sind).tri(1) = tind;
						tri(tind).seg(j) = sind;
						tri(tind).sgn(j) = -1;
						goto NEXTTRISIDE;
					}
				}
				sindprev = sind;
				sind = seg(sind).info;
			}
			/* NEW SIDE */
			seg(nseg).pnt(0) = v1;
			seg(nseg).pnt(1) = v2;
			seg(nseg).tri(0) = tind;
			seg(nseg).tri(1) = -1;
			tri(tind).seg(j) = nseg;
			tri(tind).sgn(j) = 1;
			seg(nseg).info = -1;
			if (pnt(minv).info < 0)
				pnt(minv).info = nseg;
			else 
				seg(sindprev).info = nseg;
			++nseg;
NEXTTRISIDE:
			temp = vout;
			vout = v1;
			v1 = v2;
			v2 = temp;
		}
	}

	return;
}

void face_bdry::match_tri_and_seg(void) {
	int i,j,tind,v1,v2,vout,temp,minv,maxv,order,sind,sindprev;
	
	for(i=0;i<npnt;++i)
		pnt(i).info = -1;
	
	
	// pnt(minv).info points to side
	// seg(sind).info points to next side connected to minv
	for(i=0;i<nseg;++i) {
		v1 = seg(i).pnt(0);//fix me temp switch these to be consistent?
		v2 = seg(i).pnt(1);
		minv = (v1 < v2 ? v1 : v2);
		seg(i).info = -1;
		if ((sind = pnt(minv).info) < 0)
			pnt(minv).info = i;
		else {
			while (sind >= 0) {
				sindprev = sind;
				sind = seg(sind).info;
			}
			seg(sindprev).info = i;
		}
	}

	for(i=0;i<nseg;++i) {
		seg(i).tri(0) = -1;
		seg(i).tri(1) = -1;
	}
		

	for(tind=0;tind<ntri;++tind) {
		vout = tri(tind).pnt(0);
		v1 = tri(tind).pnt(1);
		v2 = tri(tind).pnt(2);
		for(j=0;j<3;++j) {
			/* CHECK IF SIDE HAS BEEN CREATED ALREADY */
			if (v2 > v1) {
				minv = v1;
				maxv = v2;
				order = 0;
			}
			else {
				minv = v2;
				maxv = v1;
				order = 1;
			}
			sind = pnt(minv).info;
			while (sind >= 0) {
				if (maxv == seg(sind).pnt(1)) {
					seg(sind).tri(order) = tind;
					tri(tind).seg(j) = sind;
					tri(tind).sgn(j) = 1 -2*order;
					goto NEXTTRISIDE;
				}
				if (maxv == seg(sind).pnt(0)) {
					seg(sind).tri(1-order) = tind;
					tri(tind).seg(j) = sind;
					tri(tind).sgn(j) = 2*order -1;
					goto NEXTTRISIDE;
				}
				sind = seg(sind).info;
			}
			*x.gbl->log << "didn't match side: " << v1 << ' ' << v2 << ' ' << pnt(v1).gindx << ' ' << pnt(v2).gindx << std::endl;
			for(i=0;i<npnt;++i)
				*x.gbl->log << "pnt " << i << ' ' << pnt(i).gindx;
			
			
			for(i=0;i<nseg;++i)
				*x.gbl->log << "seg " << i << ": " << seg(i).pnt(0) << ' ' << seg(i).pnt(1) << ' ' << seg(i).gindx << x.seg(seg(i).gindx).pnt(0) << ' ' << x.seg(seg(i).gindx).pnt(1) << std::endl;
		
			
			for(i=0;i<ntri;++i) {
				*x.gbl->log << "tri " << i << ": " << tri(i).pnt(0) << ' ' << tri(i).pnt(1) << ' ' << tri(i).pnt(2) << ' ' << tri(i).gindx << ' ' << x.tri(tri(i).gindx).pnt(0) << ' ' << x.tri(tri(i).gindx).pnt(1) << ' ' << x.tri(tri(i).gindx).pnt(2) << std::endl;
			}
			
			exit(1);
			
NEXTTRISIDE:
			temp = vout;
			vout = v1;
			v1 = v2;
			v2 = temp;
		}
	}

	return;
}


void face_bdry::create_pnt_tri(void) {
	int i,tind;
	
	/* THIS ALLOWS US TO GET TO LOCAL HIGHER ENTITIES FROM VERTEX NUMBER */
	for (tind=0;tind<ntri;++tind)
		for(i=0;i<3;++i)
			pnt(tri(tind).pnt(i)).tri = tind;
	
	return;
}

/* CALCULATE NUMBER OF NEIGHBORS TO EACH CELL */
void face_bdry::create_pnt_nnbor(void) {
	int i;
	
	for (i=0;i<npnt;++i)
		pnt(i).nnbor = 0;
	
	for(i=0;i<nseg;++i) {
		++pnt(seg(i).pnt(0)).nnbor;
		++pnt(seg(i).pnt(1)).nnbor;
	}

	return;
}

/* CREATES TRIANGLE TO TRIANGLE POINTER */
void face_bdry::create_tri_tri(void) {
	int tind,sind,j,flip;
	
	for(tind=0;tind<ntri;++tind) {
		for(j=0;j<3;++j) {
			sind = tri(tind).seg(j);
			flip = (1 +tri(tind).sgn(j))/2;
			tri(tind).tri(j) = seg(sind).tri(flip);
		}
	}
	
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
