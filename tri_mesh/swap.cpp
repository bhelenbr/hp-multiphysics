/*
 *  edgeswap.cpp
 *  mblock
 *
 *  Created by helenbrk on Mon Sep 17 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#include "tri_mesh.h"
#include<assert.h>
#include<float.h>

//#define DEBUG_ADAPT

#ifdef DEBUG_ADAPT
extern int adapt_count;
static std::string adapt_file;
#endif

int tri_mesh::swap(int sind, FLT tol) {
	int j,t1,t2,s1p,s2p,snum1,snum2,tind,sind1,p0,p1,vt1,vt2,dir;

	t1 = seg(sind).tri(0);
	t2 = seg(sind).tri(1);

	if (t1 < 0 || t2 < 0) return(0);

	for(snum1 = 0; snum1 < 3; ++snum1)
		if (tri(t1).seg(snum1) == sind) break;

	for(snum2 = 0; snum2 < 3; ++snum2)
		if (tri(t2).seg(snum2) == sind) break;

	p0 = seg(sind).pnt(0);
	p1 = seg(sind).pnt(1);
	vt1 = tri(t1).pnt(snum1);
	vt2 = tri(t2).pnt(snum2);

	if (MIN(minangle(p0,p1,vt1),minangle(p1,p0,vt2)) >
		MIN(minangle(vt2,vt1,p0),minangle(vt1,vt2,p1)) -tol -10.0*EPSILON) return(0);

	/* MARK TOUCHED */
	tri(sind).info |= STOUC;
	tri(t1).info |= TTOUC;
	tri(t2).info |= TTOUC;

	/* SWAP SIDE */
	seg(sind).pnt(0) = vt2;
	seg(sind).pnt(1) = vt1;

	/* KEEP 2/3 POINTS IN SAME SPOT */
	tri(t1).pnt((snum1 +2)%3) = vt2;
	tri(t2).pnt((snum2 +2)%3) = vt1;

	s1p = (snum1 +1)%3;
	s2p = (snum2 +1)%3;

	/* FIX 2 CHANGED EXTERIOR SIDES */
	tri(t1).seg(snum1) = tri(t2).seg(s2p);
	tri(t1).sgn(snum1) = tri(t2).sgn(s2p);
	tri(t1).tri(snum1) = tri(t2).tri(s2p);

	tri(t2).seg(snum2) = tri(t1).seg(s1p);
	tri(t2).sgn(snum2) = tri(t1).sgn(s1p);
	tri(t2).tri(snum2) = tri(t1).tri(s1p);

	/* FIX STRI/TTRI FOR 2 CHANGED EXTERIOR SIDES */
	sind1 = tri(t1).seg(snum1);
	dir = (1 -tri(t1).sgn(snum1))/2;
	seg(sind1).tri(dir) = t1;
	tind = tri(t1).tri(snum1);
	if (tind > -1) {
		for(j=0;j<3;++j)
			if (tri(tind).seg(j) == sind1) break;
		tri(tind).tri(j) = t1;
	}

	sind1 = tri(t2).seg(snum2);
	dir = (1 -tri(t2).sgn(snum2))/2;
	seg(sind1).tri(dir) = t2;
	tind = tri(t2).tri(snum2);
	if (tind > -1) {
		for(j=0;j<3;++j)
			if (tri(tind).seg(j) == sind1) break;
		tri(tind).tri(j) = t2;
	}

	pnt(p0).tri = t1;
	pnt(p1).tri = t2;

	tri(t1).tri(s1p) = t2;
	tri(t2).tri(s2p) = t1;

	tri(t1).seg(s1p) = sind;
	tri(t2).seg(s2p) = sind;
	tri(t1).sgn(s1p) =  1;
	tri(t2).sgn(s2p) = -1;

#ifdef DEBUG_ADAPT
		std::ostringstream nstr;
		nstr << adapt_count++ << std::flush;
		adapt_file = "adapt" +nstr.str();
		nstr.str("");
		output(adapt_file,debug_adapt);
#endif

	return(1);
}


void tri_mesh::swap(int nswp, int *swp, FLT tol) {
	int i,flag;

	do {
		flag = 0;
		for(i=0;i<nswp;++i)
			flag += swap(swp[i],tol);
	} while(flag > 0);

	return;
}

void tri_mesh::swap(FLT tol) {
	int nswap,i;

	/* PERFORM EDGE SWAPPING */
	do {
		nswap = 0;
		for(i=0;i<nseg;++i) {
			nswap += swap(i,tol);
		}
		*gbl->log << "#Swap cycle finished: " << nswap << " sides swapped" << std::endl;
	} while(nswap > 0);

	return;
}

