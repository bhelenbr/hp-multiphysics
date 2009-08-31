/*
 *  insert.cpp
 *  mblock
 *
 *  Created by helenbrk on Fri Aug 31 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#include "tri_mesh.h"
#include <math.h>
#include <float.h>
#include <assert.h>
#include <stdlib.h>

int tri_mesh::insert(const TinyVector<FLT,ND> &x) {
	int n,tind,pnear,err;

	for(n=0;n<ND;++n)
		pnts(npnt)(n) = x(n);

	qtree.addpt(npnt);
	qtree.nearpt(npnt,pnear);

	/* FIND TRIANGLE CONTAINING POINT */
	bool isfound = findtri(x,pnear,tind);
	if (!isfound) {
		std::cerr << "couldn't find triangle for point: " << x(0) << ' ' << x(1) << " pnear: " << pnear << std::endl;
		std::cerr << "maxsrch: " << gbl->maxsrch << "vtri: " << pnt(pnear).tri << std::endl;
		output("error");
		exit(1);
	}
	if (npnt >= maxpst) {
		*gbl->log << "need to use larger growth factor: too many vertices" << std::endl;
		output("error");
		exit(1);
	}
	err = insert(npnt,tind);
	npnt += 1 -err;

	return(err);
}

int tri_mesh::insert(int pnum, int tnum) {
	int ntdel, nskeep, nsdel;
	int i,j,tind,tin,tnxt,p0,sstart,dir,rstrt;
	int sind,sind1,snum;

	/* SEARCH SURROUNDING FOR NONDELAUNEY TRIANGLES */
	/* SIDES ON BOUNDARY OF HOLE (SKEEP) */
	/* SIDES IN INTERIOR OF HOLE (SDEL) (STARTS AT HALFWAY THROUGH LIST) */
	ntdel = 0;
	nskeep = 0;
	nsdel = 0;

	gbl->i2wk_lst1(ntdel++) = tnum;
	tri(tnum).info |= TSRCH;

	tind = tnum;
	p0 = tri(tnum).pnt(0);
	/* FIRST SIDE DONE SEPARATELY SO KNOW WHEN TO STOP */
	snum = 2;
	tnxt = tri(tind).tri(snum);
	sind = tri(tind).seg(snum);
	dir = (1+tri(tind).sgn(snum))/2;
	sstart = sind;
	rstrt = 1;

	if (tnxt < 0) {
		gbl->i2wk_lst2(nskeep++) = sind;
		gbl->i2wk_lst2(nskeep++) = dir;
		p0 = seg(sind).pnt(dir);
		snum = (snum+1)%3;
	}
	else if (incircle(tnxt,pnts(pnum)) > 0.0) {
		if (dir > 0) gbl->i2wk_lst3(nsdel++) = sind;
		rstrt = 2;
		gbl->i2wk_lst1(ntdel++) = tnxt;
		tri(tnxt).info |= TSRCH;
		tind = tnxt;
		for(snum=0;snum<3;++snum) {
			if (tri(tind).pnt(snum) == p0) break;
		}
		snum = (snum+2)%3;
	}
	else {
		gbl->i2wk_lst2(nskeep++) = sind;
		gbl->i2wk_lst2(nskeep++) = dir;
		p0 = seg(sind).pnt(dir);
		snum = (snum+1)%3;
	}
	tnxt = tri(tind).tri(snum);
	sind = tri(tind).seg(snum);
	dir = (1+tri(tind).sgn(snum))/2;

	/* GO COUNTER-CLOCKWISE AROUND POINTS OF HOLE */
	/* IF START SIDE IS IN THE INTERIOR MUST HIT IT TWICE */
	for(j=0;j<rstrt;++j) {
		do  {
			if (tnxt < 0) {
				gbl->i2wk_lst2(nskeep++) = sind;
				gbl->i2wk_lst2(nskeep++) = dir;
				p0 = seg(sind).pnt(dir);
				snum = (snum+1)%3;
			}
			else if (tri(tnxt).info&TSRCH) {
				if (dir > 0) gbl->i2wk_lst3(nsdel++) = sind;
				tind = tnxt;
				for(snum=0;snum<3;++snum) {
					if (tri(tind).pnt(snum) == p0) break;
				}
				snum = (snum+2)%3;
			}
			else if (incircle(tnxt,pnts(pnum)) > 0.0) {
				if (dir > 0) gbl->i2wk_lst3(nsdel++) = sind;
				gbl->i2wk_lst1(ntdel++) = tnxt;
				tri(tnxt).info |= TSRCH;
				tind = tnxt;
				for(snum=0;snum<3;++snum) {
					if (tri(tind).pnt(snum) == p0) break;
				}
				snum = (snum+2)%3;
			}
			else {
				gbl->i2wk_lst2(nskeep++) = sind;
				gbl->i2wk_lst2(nskeep++) = dir;
				p0 = seg(sind).pnt(dir);
				snum = (snum+1)%3;
			}
			tnxt = tri(tind).tri(snum);
			sind = tri(tind).seg(snum);
			dir = (1+tri(tind).sgn(snum))/2;
		} while (sind != sstart);
	}

	/* RESET TSRCH FLAGS */
	for(i=0;i<ntdel;++i)
		tri(gbl->i2wk_lst1(i)).info &= ~TSRCH;

	nskeep = nskeep >> 1;

	/*	CHECK THAT WE AREN'T INSERTING POINT VERY CLOSE TO BOUNDARY */
//    for(i=0;i<nskeep;++i) {
//        sind = gbl->i2wk_lst2(i);
//        if(fabs(minangle(pnum, seg(sind).pnt(0) , seg(sind).pnt(1))) < 6.0*M_PI/180.0) {
//            *gbl->log << "#Warning: inserting close to boundary" << std::endl;
//        }
//    }

	/* APPEND NEW INDICES TO LIST OF INDICES TO USE FOR NEW SIDES & TRIS*/
	for(i=nsdel;i<nskeep;++i)
		gbl->i2wk_lst3(i) = nseg +(i-nsdel);

	for(i=ntdel;i<nskeep;++i)
		gbl->i2wk_lst1(i) = ntri +(i-ntdel);

	/* PERIODIC TRIANGLE */
	gbl->i2wk_lst1(nskeep) = gbl->i2wk_lst1(0);

	ntri += 2;
	nseg += 3;

	if (ntri > maxpst || nseg > maxpst) {
		*gbl->log << "need to use bigger growth factor: too many sides/tris" << nseg << ntri << std::endl;
		output("error");
		exit(1);
	}

	for(i=0;i<nskeep;++i) {
		tind = gbl->i2wk_lst1(i);
		tri(tind).info |= TTOUC;

		tnxt = gbl->i2wk_lst1(i+1);
		sind = gbl->i2wk_lst2(i<<1);
		dir = gbl->i2wk_lst2((i<<1) +1);

		/* CREATE NEW INFO */
		p0 = seg(sind).pnt(1-dir);
		tri(tind).pnt(0) = p0;
		pnt(p0).tri = tind;
		p0 = seg(sind).pnt(dir);
		tri(tind).pnt(1) = p0;
		pnt(p0).tri = tind;
		tri(tind).pnt(2) = pnum;
		pnt(pnum).tri = tind;

		/* SIDE 2 INFO */
		tri(tind).seg(2) = sind;
		tri(tind).sgn(2) = -1 +2*dir;

		seg(sind).tri(1-dir) = tind;
		tin = seg(sind).tri(dir);
		tri(tind).tri(2) = tin;
		if (tin > -1) {
			for(j=0;j<3;++j) {
				if (tri(tin).seg(j) == sind) {
					tri(tin).tri(j) = tind;
					break;
				}
			}
		}

		/* CREATE SIDE 0 */
		p0 = seg(sind).pnt(dir);
		sind1 = gbl->i2wk_lst3(i);
		seg(sind1).tri(0) = tind;
		seg(sind1).pnt(0) = p0;
		seg(sind1).pnt(1) = pnum;
		tri(sind1).info |= STOUC;


		tri(tind).seg(0) = sind1;
		tri(tind).sgn(0) = 1;
		tri(tind).tri(0) = tnxt;

		seg(sind1).tri(1) = tnxt;
		tri(tnxt).seg(1) = sind1;
		tri(tnxt).sgn(1) = -1;
		tri(tnxt).tri(1) = tind;
	}

	gbl->i2wk_lst1(-1) = ntdel;
	gbl->i2wk_lst2(-1) = nskeep;
	gbl->i2wk_lst3(-1) = nsdel;

	return(0);
}

void tri_mesh::bdry_insert(int pnum, int sind, int endpt) {
	int ntdel, nsdel, nskeep;
	int sstart;
	int i,j,tin,tind,tnxt,p0,dir;
	int snum, sind1;

	/* ADD POINT TO QUADTREE */
	qtree.addpt(pnum);
	tri(pnum).info |= PTOUC;

	/* SEARCH SURROUNDING FOR NONDELAUNEY TRIANGLES */
	/* SIDES ON BOUNDARY OF HOLE (SKEEP) */
	/* SIDES IN INTERIOR OF HOLE (SDEL) (STARTS AT HALFWAY THROUGH LIST) */
	ntdel = 0;
	nskeep = 0;
	nsdel = 0;

	tind = seg(sind).tri(0);
	gbl->i2wk_lst1(ntdel++) = tind;
	tri(tind).info |= TSRCH;
	for(snum=0;snum<3;++snum)
		if (tri(tind).seg(snum) == sind) break;

	sstart = sind;
	p0 = seg(sind).pnt(1);
	snum = (snum+1)%3;
	tnxt = tri(tind).tri(snum);
	sind = tri(tind).seg(snum);
	dir = (1+tri(tind).sgn(snum))/2;

	/* GO COUNTER-CLOCKWISE AROUND POINTS OF HOLE */
	do  {
		if (tnxt < 0) {
			gbl->i2wk_lst2(nskeep++) = sind;
			gbl->i2wk_lst2(nskeep++) = dir;
			p0 = seg(sind).pnt(dir);
			snum = (snum+1)%3;
		}
		else if (tri(tnxt).info&TSRCH) {
			if (dir > 0) gbl->i2wk_lst3(nsdel++) = sind;
			tind = tnxt;
			for(snum=0;snum<3;++snum) {
				if (tri(tind).pnt(snum) == p0) break;
			}
			snum = (snum+2)%3;
		}
		else if (incircle(tnxt,pnts(pnum)) > 0.0) {
			if (dir > 0) gbl->i2wk_lst3(nsdel++) = sind;
			gbl->i2wk_lst1(ntdel++) = tnxt;
			tri(tnxt).info |= TSRCH;
			tind = tnxt;
			for(snum=0;snum<3;++snum) {
				if (tri(tind).pnt(snum) == p0) break;
			}
			snum = (snum+2)%3;
		}
		else {
			gbl->i2wk_lst2(nskeep++) = sind;
			gbl->i2wk_lst2(nskeep++) = dir;
			p0 = seg(sind).pnt(dir);
			snum = (snum+1)%3;
		}
		tnxt = tri(tind).tri(snum);
		sind = tri(tind).seg(snum);
		dir = (1+tri(tind).sgn(snum))/2;
	} while (sind != sstart);

	/* RESET TSRCH FLAGS */
	for(i=0;i<ntdel;++i)
		tri(gbl->i2wk_lst1(i)).info &= ~TSRCH;

	/* ALTER OLD BOUNDARY SIDE & CREATE NEW SIDE */
	seg(nseg).pnt(endpt) = pnum;
	seg(nseg).pnt(1-endpt) = seg(sind).pnt(1-endpt);
	seg(sind).pnt(1-endpt) = pnum;
	tri(sind).info |= STOUC;
	tri(nseg).info |= STOUC;

	/* ADD NEW SIDE TO BOUNDARY GROUP */
	/* NEED TO REORDER WHEN FINISHED */
	i = getbdrynum(seg(sind).tri(1));
	seg(nseg).tri(1) = trinumatbdry(i,ebdry(i)->nseg);
	if (ebdry(i)->nseg >= ebdry(i)->maxseg) {
		output("error");
		*gbl->log << "need to use bigger growth factor (too many boundary sides)" << std::endl;
		exit(1);
	}
	ebdry(i)->seg(ebdry(i)->nseg++) = nseg;
	++nseg;

	nskeep = nskeep>>1;

	/* APPEND NEW INDICES TO LIST OF INDICES TO USE FOR NEW SIDES & TRIS*/
	for(i=nsdel;i<nskeep;++i)
		gbl->i2wk_lst3(i) = nseg +(i-nsdel);

	for(i=ntdel;i<nskeep;++i)
		gbl->i2wk_lst1(i) = ntri +(i-ntdel);

	++ntri;
	++nseg;

	if (ntri > maxpst || nseg > maxpst) {
		*gbl->log << "need to use bigger growth factor: too many sides/tris:" << nseg << ntri << std::endl;
		output("error");
		exit(1);
	}

	/* FIX FIRST AND LAST TRIANGLE BOUNDARY SIDE */
	if (endpt) {
		tind = gbl->i2wk_lst1(0);
		tri(tind).seg(1) = sind;
		tri(tind).sgn(1) = 1;
		tri(tind).tri(1) = seg(sind).tri(1);
		seg(sind).tri(0) = tind;

		tind = gbl->i2wk_lst1(nskeep-1);
		tri(tind).seg(0) = nseg-2;
		tri(tind).sgn(0) = 1;
		tri(tind).tri(0) = seg(nseg-2).tri(1);
		seg(nseg-2).tri(0) = tind;

	}
	else {
		tind = gbl->i2wk_lst1(0);
		tri(tind).seg(1) = nseg-2;
		tri(tind).sgn(1) = 1;
		tri(tind).tri(1) = seg(nseg-2).tri(1);
		seg(nseg-2).tri(0) = tind;

		tind = gbl->i2wk_lst1(nskeep-1);
		tri(tind).seg(0) = sind;
		tri(tind).sgn(0) = 1;
		tri(tind).tri(0) = seg(sind).tri(1);
		seg(sind).tri(0) = tind;
	}


	for(i=0;i<nskeep-1;++i) {
		tind = gbl->i2wk_lst1(i);
		tri(tind).info |= TTOUC;

		tnxt = gbl->i2wk_lst1(i+1);
		sind = gbl->i2wk_lst2(i<<1);
		dir = gbl->i2wk_lst2((i<<1) +1);

		/* CREATE NEW INFO */
		p0 = seg(sind).pnt(1-dir);
		tri(tind).pnt(0) = p0;
		pnt(p0).tri = tind;
		p0 = seg(sind).pnt(dir);
		tri(tind).pnt(1) = p0;
		pnt(p0).tri = tind;
		tri(tind).pnt(2) = pnum;
		pnt(pnum).tri = tind;

		/* SIDE 2 INFO */
		tri(tind).seg(2) = sind;
		tri(tind).sgn(2) = -1 +2*dir;

		seg(sind).tri(1-dir) = tind;
		tin = seg(sind).tri(dir);
		tri(tind).tri(2) = tin;
		if (tin > -1) {
			for(j=0;j<3;++j) {
				if (tri(tin).seg(j) == sind) {
					tri(tin).tri(j) = tind;
					break;
				}
			}
		}

		/* CREATE SIDE 0 */
		p0 = seg(sind).pnt(dir);
		sind1 = gbl->i2wk_lst3(i);
		seg(sind1).tri(0) = tind;
		seg(sind1).pnt(0) = p0;
		seg(sind1).pnt(1) = pnum;
		tri(sind1).info |= STOUC;


		tri(tind).seg(0) = sind1;
		tri(tind).sgn(0) = 1;
		tri(tind).tri(0) = tnxt;

		seg(sind1).tri(1) = tnxt;
		tri(tnxt).seg(1) = sind1;
		tri(tnxt).sgn(1) = -1;
		tri(tnxt).tri(1) = tind;
	}

	/* LAST TRIANGLE */
	i = nskeep-1;
	tind = gbl->i2wk_lst1(i);
	tri(tind).info |= TTOUC;

	sind = gbl->i2wk_lst2(i<<1);
	dir = gbl->i2wk_lst2((i<<1) +1);

	/* CREATE NEW INFO */
	p0 = seg(sind).pnt(1-dir);
	tri(tind).pnt(0) = p0;
	pnt(p0).tri = tind;
	p0 = seg(sind).pnt(dir);
	tri(tind).pnt(1) = p0;
	pnt(p0).tri = tind;
	tri(tind).pnt(2) = pnum;
	pnt(pnum).tri = tind;

	/* SIDE 2 INFO */
	tri(tind).seg(2) = sind;
	tri(tind).sgn(2) = -1 +2*dir;

	seg(sind).tri(1-dir) = tind;
	tin = seg(sind).tri(dir);
	tri(tind).tri(2) = tin;
	if (tin > -1) {
		for(j=0;j<3;++j) {
			if (tri(tin).seg(j) == sind) {
				tri(tin).tri(j) = tind;
				break;
			}
		}
	}

	return;
}

bool tri_mesh::findtri(const TinyVector<FLT,ND> x, int pnear, int& tind) {
	int i,j,vn,dir,stoptri,tin;
	int ntdel;
	int tclose,nsurround;
	FLT minclosest;
	int p0,p1,p2;
	bool found = true;
	FLT dx0,dy0,dx1,dy1,dx2,dy2;
	TinyVector<FLT,3> a;

	/* TSRCH = 0x100*0x4 */
#if ((-1)&(0x100*0x4))
#define ISSRCH(A) (!((A)&TSRCH))
#define SETSRCH(A) A&=(~TSRCH)
#define CLRSRCH(A) A|=(TSRCH)
#else
#define ISSRCH(A) (((A)&TSRCH))
#define SETSRCH(A) A|=(TSRCH)
#define CLRSRCH(A) A&=(~TSRCH)
#endif

	/* HERE WE USE gbl->intwk & gbl->i2wk THIS MUST BE -1 BEFORE USING */
	tind = pnt(pnear).tri;
	stoptri = tind;
	dir = 1;
	ntdel = 0;
	do {
		if (intri(tind,x) < area(tind)*10.*EPSILON) goto FOUND;
		SETSRCH(gbl->intwk(tind));
		gbl->i2wk(ntdel++) = tind;

		for(vn=0;vn<3;++vn)
			if (tri(tind).pnt(vn) == pnear) break;

		tind = tri(tind).tri((vn +dir)%3);
		if (tind < 0) {
			if (dir > 1) break;
			/* REVERSE DIRECTION AND GO BACK TO START */
			++dir;
			tind = pnt(pnear).tri;
			for(vn=0;vn<3;++vn)
				if (tri(tind).pnt(vn) == pnear) break;

			tind = tri(tind).tri((vn +dir)%3);
			if (tind < 0) break;
			stoptri = -1;
		}
	} while(tind != stoptri);

	nsurround = ntdel;

	/* DIDN'T FIND TRIANGLE */
	/* NEED TO SEARCH SURROUNDING TRIANGLES */
	for(i=0;i<ntdel;++i) {
		tin = gbl->i2wk(i);
		for(j=0;j<3;++j) {
			tind = tri(tin).tri(j);
			if (tind < 0) continue;
			if (ISSRCH(gbl->intwk(tind))) continue;
			SETSRCH(gbl->intwk(tind));
			gbl->i2wk(ntdel++) = tind;
			if (intri(tind,x) < area(tind)*10.*EPSILON) goto FOUND;
		}
		if (ntdel >= gbl->maxsrch-4) break;
	}
	//    std::cerr << "couldn't find tri for point " << x[0] << ' ' << x[1] << ' ' << pnear << std::endl;
	minclosest = -1.0e16;
	tclose = -1;
	for (i=0;i<ntdel;++i) {
		tind = gbl->i2wk(i);

		p0 = tri(tind).pnt(0);
		p1 = tri(tind).pnt(1);
		p2 = tri(tind).pnt(2);

		dx0 =  (x(0) -pnts(p0)(0));
		dy0 =  (x(1) -pnts(p0)(1));
		dx1 =  (x(0) -pnts(p1)(0));
		dy1 =  (x(1) -pnts(p1)(1));
		dx2 =  (x(0) -pnts(p2)(0));
		dy2 =  (x(1) -pnts(p2)(1));

		a(0) = (dy2*dx1 -dx2*dy1);
		a(1) = (dy0*dx2 -dx0*dy2);
		a(2) = (dy1*dx0 -dx1*dy0);

		/* FIND NEGATIVE SIDE */
		/* CHECK IF 2 SIDES POSITIVE & 1 NEGATIVE */
		if (a(0)*a(1)*a(2) < 0) {
			a(0) /= distance(p2,p1);
			a(1) /= distance(p0,p2);
			a(2) /= distance(p1,p0);

			for (int j=0;j<3;++j) {
				if (a(j) < 0.0 && a(j) > minclosest) {
					minclosest = a(j);
					tclose = tind;
					break;
				}
			}
		}
	}
	if (tclose < 0) {
		*gbl->log << "Major Trouble in Findtri " << x << ' ' << pnear << '\n';
		exit(1);
	}

	intri(tclose,x);
	tind = tclose;
	found = false;

FOUND:
	/* RESET gbl->intwkW1 TO -1 */
	for(i=0;i<ntdel;++i) {
		CLRSRCH(gbl->intwk(gbl->i2wk(i)));
	}

	return(found);
}

bool tri_mesh::findtri(TinyVector<FLT,ND> x, int& tind) {
	int i,j,tin;
	int ntdel;
	int tclose;
	FLT minclosest;
	int p0,p1,p2;
	FLT dx0,dy0,dx1,dy1,dx2,dy2;
	TinyVector<FLT,3> a;
	bool found = true;

	ntdel = 0;
	if (intri(tind,x) < area(tind)*10.*EPSILON) return(true);
	SETSRCH(gbl->intwk(tind));
	gbl->i2wk(ntdel++) = tind;

	/* SEARCH SURROUNDING TRIANGLES */
	for(i=0;i<ntdel;++i) {
		tin = gbl->i2wk(i);
		for(j=0;j<3;++j) {
			tind = tri(tin).tri(j);
			if (tind < 0) continue;
			if (ISSRCH(gbl->intwk(tind))) continue;
			SETSRCH(gbl->intwk(tind));
			gbl->i2wk(ntdel++) = tind;
			if (intri(tind,x) < area(tind)*10.*EPSILON) goto FOUND2;
		}
		if (ntdel >= gbl->maxsrch-4) break;
	}
	//    std::cerr << "couldn't find tri for point " << x[0] << ' ' << x[1] << ' ' << pnear << std::endl;
	minclosest = -1.0e16;
	tclose = -1;
	for (i=0;i<ntdel;++i) {
		tind = gbl->i2wk(i);

		p0 = tri(tind).pnt(0);
		p1 = tri(tind).pnt(1);
		p2 = tri(tind).pnt(2);

		dx0 =  (x(0) -pnts(p0)(0));
		dy0 =  (x(1) -pnts(p0)(1));
		dx1 =  (x(0) -pnts(p1)(0));
		dy1 =  (x(1) -pnts(p1)(1));
		dx2 =  (x(0) -pnts(p2)(0));
		dy2 =  (x(1) -pnts(p2)(1));

		a(0) = (dy2*dx1 -dx2*dy1);
		a(1) = (dy0*dx2 -dx0*dy2);
		a(2) = (dy1*dx0 -dx1*dy0);

		/* FIND NEGATIVE SIDE */
		/* CHECK IF 2 SIDES POSITIVE & 1 NEGATIVE */
		if (a(0)*a(1)*a(2) < 0) {
			a(0) /= distance(p2,p1);
			a(1) /= distance(p0,p2);
			a(2) /= distance(p1,p0);

			for (int j=0;j<3;++j) {
				if (a(j) < 0.0 && a(j) > minclosest) {
					minclosest = a(j);
					tclose = tind;
					break;
				}
			}
		}
	}
	if (tclose < 0) {
		*gbl->log << "Major Trouble in Findtri " << x << ' ' << tind << '\n';
		exit(1);
	}

	intri(tclose,x);
	tind = tclose;
	found = false;

FOUND2:
	/* RESET gbl->intwkW1 TO -1 */
	for(i=0;i<ntdel;++i) {
		CLRSRCH(gbl->intwk(gbl->i2wk(i)));
	}

	return(found);
}


