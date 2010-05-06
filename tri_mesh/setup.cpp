#include "tri_mesh.h"
#include <utilities.h>
#include <float.h>

/* CREATE SIDELIST FROM TRIANGLE VERTEX LIST */
/* USES PINFO TO STORE FIRST SIND FROM VERTEX */
/* USES SINFO TO STORE NEXT SIND FROM SIND */
/* TRI.PNT MUST BE COUNTERCLOCKWISE ORDERED */
void tri_mesh::createseg(void) {
	int i,j,tind,p1,p2,vout,temp,minp,maxp,order,sind;
	int sindprev = 0; // To avoid may be used uninitialized warning

	for(i=0;i<npnt;++i)
		pnt(i).info = -1;

	nseg = 0;
	for(tind=0;tind<ntri;++tind) {
		vout = tri(tind).pnt(0);
		p1 = tri(tind).pnt(1);
		p2 = tri(tind).pnt(2);
		for(j=0;j<3;++j) {
			/* CHECK IF SIDE HAS BEEN CREATED ALREADY */
			if (p2 > p1) {
				minp = p1;
				maxp = p2;
				order = 0;
			}
			else {
				minp = p2;
				maxp = p1;
				order = 1;
			}

			sind = pnt(minp).info;
			while (sind >= 0) {
				if (maxp == seg(sind).pnt(order)) {
					if (seg(sind).tri(1) >= 0) {
						*gbl->log << "Error: seg " << sind << "has been matched with Triangle" << tind << "3 times" << std::endl;                        exit(1);
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
			seg(nseg).pnt(0) = p1;
			seg(nseg).pnt(1) = p2;
			seg(nseg).tri(0) = tind;
			seg(nseg).tri(1) = -1;
			tri(tind).seg(j) = nseg;
			tri(tind).sgn(j) = 1;
			seg(nseg).info = -1;
			if (pnt(minp).info < 0)
				pnt(minp).info = nseg;
			else
				seg(sindprev).info = nseg;
			++nseg;
NEXTTRISIDE:
			temp = vout;
			vout = p1;
			p1 = p2;
			p2 = temp;
		}
	}

	return;
}

void tri_mesh::createsegtri(void) {
	int i,j,tind,p1,p2,vout,temp,minp,maxp,order,sind,sindprev;

	for(i=0;i<npnt;++i)
		pnt(i).info = -1;

	for(i=0;i<nseg;++i) {
		p1 = seg(i).pnt(0);
		p2 = seg(i).pnt(1);
		minp = (p1 < p2 ? p1 : p2);
		sind = pnt(minp).info;
		while (sind >= 0) {
			sindprev = sind;
			sind = seg(sind).info;
		}
		if (pnt(minp).info < 0)
			pnt(minp).info = i;
		else
			seg(sindprev).info = i;
		seg(i).info = -1;
	}

	for(i=0;i<nseg;++i)
		seg(i).tri(1) = -1;

	for(tind=0;tind<ntri;++tind) {
		vout = tri(tind).pnt(0);
		p1 = tri(tind).pnt(1);
		p2 = tri(tind).pnt(2);
		for(j=0;j<3;++j) {
			/* CHECK IF SIDE HAS BEEN CREATED ALREADY */
			if (p2 > p1) {
				minp = p1;
				maxp = p2;
				order = 0;
			}
			else {
				minp = p2;
				maxp = p1;
				order = 1;
			}
			sind = pnt(minp).info;
			while (sind >= 0) {
				if (maxp == seg(sind).pnt(1)) {
					seg(sind).tri(order) = tind;
					tri(tind).seg(j) = sind;
					tri(tind).sgn(j) = 1 -2*order;
					goto NEXTTRISIDE;
				}
				if (maxp == seg(sind).pnt(0)) {
					seg(sind).tri(1-order) = tind;
					tri(tind).seg(j) = sind;
					tri(tind).sgn(j) = 2*order -1;
					goto NEXTTRISIDE;
				}
				sind = seg(sind).info;
			}
			*gbl->log << "didn't match seg: " << p1 << p2 << std::endl;
			exit(1);

NEXTTRISIDE:
			temp = vout;
			vout = p1;
			p1 = p2;
			p2 = temp;
		}
	}

	return;
}


void tri_mesh::createpnttri(void) {
	int i,tind;

	/* THIS ALLOWS US TO GET TO LOCAL HIGHER ENTITIES FROM VERTEX NUMBER */
	for (tind=0;tind<ntri;++tind)
		for(i=0;i<3;++i)
			pnt(tri(tind).pnt(i)).tri = tind;

	return;
}

/* CALCULATE NUMBER OF NEIGHBORS TO EACH CELL */
void tri_mesh::cnt_nbor(void) {
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
void tri_mesh::createtritri(void) {
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

void tri_mesh::treeinit() {
	int i,j,n,sind,p0;
	FLT x1[ND], x2[ND], dx;

	for(n=0;n<ND;++n)	{
		x1[n] = pnts(0)(n);
		x2[n] = pnts(0)(n);
	}


	for (i=0;i<nebd;++i) {
		for(j=0;j<ebdry(i)->nseg;++j) {
			sind = ebdry(i)->seg(j);
			p0 = seg(sind).pnt(0);
			for(n=0;n<ND;++n) {
				x1[n] = MIN(x1[n],pnts(p0)(n));
				x2[n] = MAX(x2[n],pnts(p0)(n));
			}
		}
	}

//	std::cout << gbl->idprefix << " tree init bl: " << x1[0] << ' ' <<  x1[1] << std::endl;
//	std::cout << gbl->idprefix << " tree init tr: " << x2[0] << ' ' <<  x2[1] << std::endl;

	for(n=0;n<ND;++n) {
		dx = MAX(x2[n]-x1[n],100.0*EPSILON);
		x1[n] -= 0.25*dx;
		x2[n] += 0.25*dx;
	}

	treeinit(x1,x2);

//	std::cout << gbl->idprefix << " tree init bl: " << x1[0] << ' ' <<  x1[1] << std::endl;
//	std::cout << gbl->idprefix << " tree init tr: " << x2[0] << ' ' <<  x2[1] << std::endl;

	return;
}

void tri_mesh::treeinit(FLT x1[ND], FLT x2[ND]) {

	qtree.init(x1,x2);

	for(int i=0;i<npnt;++i)
		qtree.addpt(i);

	return;
}

/* FIX STRI TTRI TO POINT TO GROUP/SIDE ON BOUNDARY */
void tri_mesh::bdrylabel() {
	int i,j,k,sind,tind;

	for(i=0;i<nebd;++i) {
		for(j=0;j<ebdry(i)->nseg;++j) {
			sind = ebdry(i)->seg(j);
			seg(sind).tri(1) = trinumatbdry(i,j);
			tind = seg(sind).tri(0);
			for(k=0;k<3;++k)
				if (tri(tind).seg(k) == sind) break;

			tri(tind).tri(k) = seg(sind).tri(1);
		}
	}

	return;
}

void tri_mesh::initlngth() {
	int i,j,p0,p1;
	FLT l;

	for(i=0;i<npnt;++i)
		lngth(i) = 0.0;

	for(i=0;i<nseg;++i) {
		p0 = seg(i).pnt(0);
		p1 = seg(i).pnt(1);
		l = distance(seg(i).pnt(0),seg(i).pnt(1));
		lngth(p0) += l;
		lngth(p1) += l;
	}

	for(i=0;i<npnt;++i)
		lngth(i) /= pnt(i).nnbor;


	for(i=0;i<nebd;++i) {
		for(j=0;j<ebdry(i)->nseg;++j) {
			p0 = seg(ebdry(i)->seg(j)).pnt(0);
			lngth(p0) = 1.0e32;
		}
	}

	for(i=0;i<nebd;++i) {
		for(j=0;j<ebdry(i)->nseg;++j) {
			p0 = seg(ebdry(i)->seg(j)).pnt(0);
			p1 = seg(ebdry(i)->seg(j)).pnt(1);
			l = distance(p0,p1);
			lngth(p0) = MIN(l,lngth(p0));
			lngth(p1) = MIN(l,lngth(p1));
		}
	}
	
	smooth_lngth(2);

	return;
}

