/*
 *  srebay.cpp
 *  mblock
 *
 *  Created by helenbrk on Wed Sep 12 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */
#include "tri_mesh.h"
#include <assert.h>
#include <math.h>

#define REBAY
//#define DEBUG_ADAPT
#define VERBOSE

#ifdef DEBUG_ADAPT
int adapt_count = 0;
static std::string adapt_file;
#endif

/*  fscr1 = ratio of triangle size to target
	seg(i).info = triangle refinement queue
	pnt(tind).info = back reference from triangle into queue
*/

void tri_mesh::rebay(FLT tolsize) {
	int i,j,n,tind,tfind,p0,p1,p2,pnear,nsnew,ntnew,snum,intrcnt,err;
	TinyVector<FLT,ND> xpt,dx,xmid,xdif,rn;
	TinyVector<FLT,3> wt;
	FLT maxvl;
	FLT norm,p,q,s1sq,s2sq,rad1,rad2,rs;
	FLT densty,cirrad,arg,dist,rsign;

	/* TO ADJUST FOR CIRCUMSCRIBED RADIUS */
	tolsize /= sqrt(3.);

	/* COUNTER TO SEE HOW MANY POINTS INSERTED */
	intrcnt = 0;

	/* CLASSIFY TRIANGLES AS ACCEPTED OR UNACCEPTED */
	gbl->nlst = 0;
	for(i=0;i<ntri;++i) {
		if (tri(i).info&TDLTE) continue;
		maxvl = lngth(tri(i).pnt(0));
		maxvl = MAX(maxvl,lngth(tri(i).pnt(1)));
		maxvl = MAX(maxvl,lngth(tri(i).pnt(2)));
		gbl->fltwk(i) = circumradius(i)/maxvl;
		if (gbl->fltwk(i) > tolsize) putinlst(i);
	}

	/* BEGIN REFINEMENT ALGORITHM */
	while (gbl->nlst > 0) {
		for(i=gbl->nlst-1;i>=0;--i) {
			
			/* FIND TRIANGLE FACES ON BOUNDARY FIRST */
			for(j=0;j<3;++j) {
				tind = tri(seg(i).info).tri(j);
				if (tind < 0)  {
					snum = j;
					tind = seg(i).info;
					goto TFOUND;
				}
			}
			
			/* OR FIND TRIANGLES WITH 2 FACES ON BOUNDARY OF ACCEPTED REGIONS */
			int naccept = 0;
			for(j=0;j<3;++j) {
				tind = tri(seg(i).info).tri(j);
				if (pnt(tind).info == -1)  {
					++naccept;
					snum = j;
					tind = seg(i).info;
					if (naccept > 1) goto TFOUND;
				}
			}
		}
					
		for(i=gbl->nlst-1;i>=0;--i) {	
			/* FIND TRIANGLES WITH 1 FACE ON BOUNDARY OF ACCEPTED REGIONS */
			for(j=0;j<3;++j) {
				tind = tri(seg(i).info).tri(j);
				if (pnt(tind).info == -1)  {
					snum = j;
					tind = seg(i).info;
					goto TFOUND;
				}
			}
		}
		*gbl->log << "Didn't find triangle???" << std::endl;
		tind = seg(gbl->nlst-1).info;
		snum = 0;

		TFOUND:

		/* CHECK THAT NOT ALL SIDES ARE ACCEPTED */
		tfind = 0;
		for (j=0;j<3;++j)
			if (tri(tind).tri(j) < 0 || pnt(tri(tind).tri(j)).info == -1) ++tfind;

		if (tfind == 3) {
			tkoutlst(tind);
			continue;
		}

		if (npnt > maxpst -2) {
			*gbl->log << "too many vertices" << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
		}

		/* THIS IS TRIANGLE BASED REFINEMENT */
		/*	FIND ACCEPTED EDGE ON TRIANGLE & BUILD OFF OF IT */
		p0 = tri(tind).pnt((snum+1)%3);
		p1 = tri(tind).pnt((snum+2)%3);
		p2 = tri(tind).pnt(snum);

#ifdef REBAY
		/* USE REBAY'S ALGORITHM FOR INSERT POINT */
		p = 0.0;
		for(n=0;n<ND;++n) {
			dx(n) = pnts(p0)(n) -pnts(p1)(n);
			p += pow(dx(n),2);
		}
		p = 0.5*sqrt(p);

		s1sq = 0.0;
		for(n=0;n<ND;++n) {
			dx(n) = pnts(p2)(n) -pnts(p1)(n);
			s1sq += pow(dx(n),2);
		}
		s1sq *= 0.25;

		s2sq = 0.0;
		for(n=0;n<ND;++n) {
			dx(n) = pnts(p2)(n) -pnts(p0)(n);
			s2sq += pow(dx(n),2);
		}
		s2sq *= 0.25;
		circumcenter(tind,xpt);
		if (p*p > s1sq +s2sq) goto INSRT;

		q = 0.0;
		for(n=0;n<ND;++n) {
			xmid(n) = .5*(pnts(p0)(n) +pnts(p1)(n));
			dx(n) = xpt(n) -xmid(n);
			q += pow(dx(n),2);
		}
		q = sqrt(q);
		if (q < p) goto INSRT;

		densty = 0.5*(lngth(p0)  +lngth(p1))/sqrt(3.0);
		rad1 = MAX(densty,p);
		rad2 = .5*(p*p  +q*q)/q;
		cirrad = MIN(rad1,rad2);
		arg = fabs(cirrad*cirrad  -p*p);
		dist = cirrad  +sqrt(arg);
		rs = 0.0;
		for(n=0;n<ND;++n) {
			dx(n) = pnts(p1)(n) -pnts(p0)(n);
			rs += pow(dx(n),2);
		}
		rs = 1./sqrt(rs);
		rn(0) =  dx(1)*rs;
		rn(1) = -dx(0)*rs;
		for(n=0;n<ND;++n)
			xdif(n) = pnts(p2)(n)  -xmid(n);
		rsign = 1.;
		if (xdif(0)*rn(0) +xdif(1)*rn(1) < 0.) rsign = -1.;
		for(n=0;n<ND;++n)
			xpt(n) = xmid(n) +rsign*dist*rn(n);
#elif defined(MIDPOINT)
			/* MIDPOINT RULE (VERY SIMPLE) */
		for(n=0;n<ND;++n)
			xpt(n) = 0.5*(pnts(p1)(n) +pnts(p2)(n));
#else
		densty = 0.5*(lngth(p0)  +lngth(p1))/sqrt(3.0);
		rs = 0.0;
		for(n=0;n<ND;++n) {
			xmid(n) = .5*(pnts(p0)(n) +pnts(p1)(n));
			dx(n) = pnts(p2)(n) -xmid(n);
			rs += pow(dx(n),2);
		}
		rs = sqrt(rs);
		if (densty > rs) {
			*gbl->log << "just checking\n";
			densty = 0.0;
			for(n=0;n<ND;++n) {
				densty += pow(pnts(p1)(n) -pnts(p0)(n),2);
			}
			densty = sqrt(densty);
		}
		rs = densty/rs;

		xpt = xmid +rs*dx;
#endif

INSRT:
		/* INSERT POINT */
		for(n=0;n<ND;++n)
			pnts(npnt)(n) = xpt(n);

#ifdef DEBUG_ADAPT
		*gbl->log << "Inserting interior segment " << intrcnt << ' ' << adapt_count << ' ';
		for(n=0;n<ND;++n)
			*gbl->log << pnts(npnt)(n) << ' ';
		*gbl->log << std::endl;
#endif

		dist = qtree.nearpt(pnts(npnt),pnear);
		norm = 0.0;
		for (n=0;n<ND;++n)
			norm += fabs(pnts(npnt)(n));
		if (dist < 100.0*EPSILON*norm) {
#ifdef VERBOSE
			*gbl->log << "#Point to close to insert " << dist << std::endl;
			*gbl->log << pnts(p0) << std::endl;
			*gbl->log << pnts(p1) << std::endl;
			*gbl->log << pnts(p2) << std::endl;
			*gbl->log << pnts(npnt) << std::endl;
			*gbl->log << tri(tind).pnt <<  snum << std::endl;
			for(j=0;j<3;++j) {
				if (tri(tind).tri(j) < 0 || pnt(tri(tind).tri(j)).info == -1) {
					*gbl->log << "segment " << j << " is accepted\n";
				}
			}
#endif
			tkoutlst(tind);
			continue;
		}

		bool found = findtri(xpt,pnear,tfind);
		if (!found) {
#ifdef VERBOSE
			*gbl->log << "#Warning: Trying to insert outside domain " << std::endl;
			*gbl->log << pnts(p0) << std::endl;
			*gbl->log << pnts(p1) << std::endl;
			*gbl->log << pnts(p2) << std::endl;
			*gbl->log << pnts(npnt) << std::endl;
			*gbl->log << tri(tind).pnt << snum << std::endl;
			for(j=0;j<3;++j) {
				if (tri(tind).tri(j) < 0 || pnt(tri(tind).tri(j)).info == -1) {
					*gbl->log << "segment " << j << " is accepted\n";
				}
			}
#endif
			tkoutlst(tind);
			continue;
		}

		getwgts(wt);
		lngth(npnt) = 0.0;
		for(i=0;i<3;++i)
			lngth(npnt) += wt(i)*lngth(tri(tfind).pnt(i));

		err = insert(npnt,tfind);
		if (!err) {
			/* ADD POINT TO QUADTREE */
			tri(npnt).info |= PTOUC;
			qtree.addpt(npnt);
			nsnew = gbl->i2wk_lst3(-1) +3;
			ntnew = gbl->i2wk_lst1(-1) +2;
			++intrcnt;
		}
		else {
#ifdef VERBOSE
			*gbl->log << "#Warning: Makes Bad Triangle " << std::endl;
			*gbl->log << pnts(p0) << std::endl;
			*gbl->log << pnts(p1) << std::endl;
			*gbl->log << pnts(p2) << std::endl;
			*gbl->log << pnts(npnt) << std::endl;
			*gbl->log << tri(tind).pnt << snum << std::endl;
			for(j=0;j<3;++j) {
				if (tri(tind).tri(j) < 0 || pnt(tri(tind).tri(j)).info == -1) {
					*gbl->log << "segment " << j << " is accepted\n";
				}
			}
#endif
			tkoutlst(tind);
			continue;
		}
		++npnt;

		for(i=0;i<gbl->i2wk_lst1(-1);++i)
			if (pnt(gbl->i2wk_lst1(i)).info > -1) tkoutlst(gbl->i2wk_lst1(i));

		if (pnt(tind).info > -1) tkoutlst(tind);


		for(i=0;i<ntnew;++i) {
			tind = gbl->i2wk_lst1(i);
			maxvl = lngth(tri(tind).pnt(0));
			maxvl = MAX(maxvl,lngth(tri(tind).pnt(1)));
			maxvl = MAX(maxvl,lngth(tri(tind).pnt(2)));
			gbl->fltwk(tind) = circumradius(tind)/maxvl;
			if (gbl->fltwk(tind) > tolsize) putinlst(tind);
		}
#ifdef DEBUG_ADAPT
		std::ostringstream nstr;
		nstr << adapt_count++ << std::flush;
		adapt_file = "adapt" +nstr.str();
		nstr.str("");
		output(adapt_file,debug_adapt);
#endif

	}

	*gbl->log << "#Rebay finished: new interior points " << intrcnt << std::endl;

	return;
}

void tri_mesh::bdry_rebay(FLT tolsize) {
	int sind,p0,count,psifxpt,prev,sind_prev;
	int bseg = 0; // To avoid may be used uninitialized warning
	FLT psi;
	TinyVector<FLT,tri_mesh::ND> endpt;

	/* REFINE BOUNDARY SIDES */
	for(int bnum=0;bnum<nebd;++bnum) {
		count = 0;

		if (!ebdry(bnum)->is_frst()) {
			ebdry(bnum)->sndtype() = boundary::int_msg;
			ebdry(bnum)->sndsize() = 0;
			ebdry(bnum)->comm_prepare(boundary::all,0,boundary::master_slave);
			continue;
		}

		/* CHECK THAT ENDPOINTS ARE ON CURVE JUST TO BE SURE */
		sind = ebdry(bnum)->seg(0);
		p0 = seg(sind).pnt(0);
		endpt = pnts(p0);
		ebdry(bnum)->mvpttobdry(0,-1.0,endpt);
		endpt -= pnts(p0);
		if (fabs(endpt(0)) +fabs(endpt(1)) > FLT_EPSILON*(fabs(qtree.xmax(0)-qtree.xmin(0)) +fabs(qtree.xmax(1)-qtree.xmin(1)))) {
			*gbl->log << "first endpoint of boundary " << ebdry(bnum)->idprefix << " does not seem to be on curve\n";
			*gbl->log << "current" << pnts(p0)(0) << ',' << pnts(p0)(1) << " projected "<< (endpt+pnts(p0))(0) << ',' << (endpt+pnts(p0))(1) << '\n';
		}

		sind = ebdry(bnum)->seg(ebdry(bnum)->nseg-1);
		p0 = seg(sind).pnt(1);
		endpt = pnts(p0);
		ebdry(bnum)->mvpttobdry(ebdry(bnum)->nseg-1,1.0,endpt);
		endpt -= pnts(p0);
		if (fabs(endpt(0)) +fabs(endpt(1)) > FLT_EPSILON*(fabs(qtree.xmax(0)-qtree.xmin(0)) +fabs(qtree.xmax(1)-qtree.xmin(1)))) {
			*gbl->log << "last endpoint of boundary " << ebdry(bnum)->idprefix << " does not seem to be on curve\n";
			*gbl->log << "current" << pnts(p0)(0)  << ',' << pnts(p0)(1) << " projected " <<  (endpt+pnts(p0))(0) << ',' << (endpt+pnts(p0))(1) << '\n';
		}

		gbl->nlst = 0;
		for(int indx=0;indx<ebdry(bnum)->nseg;++indx) {
			sind = ebdry(bnum)->seg(indx);
			if (tri(sind).info&SDLTE) continue;
			gbl->fltwk(sind) = distance(seg(sind).pnt(0),seg(sind).pnt(1))/MAX(lngth(seg(sind).pnt(0)),lngth(seg(sind).pnt(1)));
			if (gbl->fltwk(sind) > tolsize) putinlst(sind);
		}

		/* SKIP FIRST SPOT SO CAN SEND LENGTH FIRST */
		ebdry(bnum)->sndsize() = 1;
		while (gbl->nlst > 0) {

			for(int i=gbl->nlst-1;i>=0;--i) {
				sind = seg(i).info;
				bseg = getbdryseg(seg(sind).tri(1));
				prev = ebdry(bnum)->prev(bseg);

				/* Check if it is the first edge. Ok to start there */
				if (prev < 0) break;

				/* Check that previous edge is not in list */
				sind_prev = ebdry(bnum)->seg(prev);
				if (pnt(sind_prev).info == -1)  break;
			}

			/* INSERT POINT TO MATCH LENGTH FUNCTION L(1+psi)/2 = lbar +psi dl/2 */
			FLT L = distance(seg(sind).pnt(0),seg(sind).pnt(1));
			FLT lbar = 0.5*(lngth(seg(sind).pnt(0)) +lngth(seg(sind).pnt(1)));
			FLT dl = lngth(seg(sind).pnt(1)) -lngth(seg(sind).pnt(0));
			psi = (lbar -L/2)/(L/2 -dl/2);
			if (psi > 0.0) psi = 0.0;
			
			psifxpt = static_cast<int>((1<<16)*psi);
			psi = psifxpt/static_cast<FLT>(1<<16);
			pnts(npnt) = 0.5*((1.-psi)*pnts(seg(sind).pnt(0)) +(1.+psi)*pnts(seg(sind).pnt(1)));
			ebdry(bnum)->mvpttobdry(bseg,psi,pnts(npnt));
			lngth(npnt) = 0.5*((1.-psi)*lngth(seg(sind).pnt(0)) +(1.+psi)*lngth(seg(sind).pnt(1)));

#ifdef DEBUG_ADAPT
			*gbl->log << "Inserting boundary segment " << count << ' ' << adapt_count << ' ' << npnt << ' ';
			for(int n=0;n<ND;++n)
				*gbl->log << pnts(npnt)(n) << ' ';
			*gbl->log << " seg " << bseg << " psi " << psi << std::endl;
#endif
			/* INSERT POINT */
			bdry_insert(npnt,sind);
			++npnt;
			++count;

			/* UPDATE NEXT & PREV POINTERS */
			ebdry(bnum)->prev(ebdry(bnum)->next(bseg)) = ebdry(bnum)->nseg-1;
			ebdry(bnum)->next(ebdry(bnum)->nseg-1) = ebdry(bnum)->next(bseg);
			ebdry(bnum)->next(bseg) = ebdry(bnum)->nseg-1;
			ebdry(bnum)->prev(ebdry(bnum)->nseg-1) = bseg;

			/* STORE ADAPTATION INFO FOR COMMUNICATION BOUNDARIES */
			ebdry(bnum)->isndbuf(ebdry(bnum)->sndsize()++) = bseg;
			ebdry(bnum)->isndbuf(ebdry(bnum)->sndsize()++) = psifxpt;

			/* UPDATE MODIFIED SIDE */
			tkoutlst(sind);
			gbl->fltwk(sind) = distance(seg(sind).pnt(0),seg(sind).pnt(1))/MAX(lngth(seg(sind).pnt(0)),lngth(seg(sind).pnt(1)));
			if (gbl->fltwk(sind) > tolsize) putinlst(sind);

			/* UPDATE NEW BOUNDARY SIDE */
			sind = ebdry(bnum)->seg(ebdry(bnum)->nseg -1);
			gbl->fltwk(sind) = distance(seg(sind).pnt(0),seg(sind).pnt(1))/MAX(lngth(seg(sind).pnt(0)),lngth(seg(sind).pnt(1)));
			if (gbl->fltwk(sind) > tolsize) putinlst(sind);
#ifdef DEBUG_ADAPT
			std::ostringstream nstr;
			nstr << adapt_count++ << std::flush;
			adapt_file = "adapt" +nstr.str();
			nstr.str("");
			output(adapt_file,debug_adapt);
#endif
		}
		ebdry(bnum)->isndbuf(0) = ebdry(bnum)->sndsize();
		ebdry(bnum)->sndtype() = boundary::int_msg;
		ebdry(bnum)->comm_prepare(boundary::all,0,boundary::master_slave);
		*gbl->log << "#Boundary refinement finished, " << ebdry(bnum)->idnum << ' ' << count << " sides added" << std::endl;
	}

	return;
}

void tri_mesh::bdry_rebay1() {
	int i,sind;
	int bseg, sndsize;
	int nel_bgn, psifxpt;
	FLT psi;

	/* REFINE MATCHING BOUNDARIES */
	for(int bnum=0;bnum<nebd;++bnum) {

		ebdry(bnum)->comm_wait(boundary::all,0,boundary::master_slave);

		if (ebdry(bnum)->is_frst() || !ebdry(bnum)->is_comm()) continue;

		nel_bgn = ebdry(bnum)->nseg;

		sndsize = ebdry(bnum)->ircvbuf(0,0);
		for(i=1;i<sndsize;i+=2) {
			bseg = ebdry(bnum)->ircvbuf(0,i);
			if (bseg < nel_bgn) bseg = nel_bgn -1 -bseg;

			sind = ebdry(bnum)->seg(bseg);
			psifxpt = -ebdry(bnum)->ircvbuf(0,i+1);
			psi = psifxpt/static_cast<FLT>(1<<16);
			pnts(npnt) = 0.5*((1.-psi)*pnts(seg(sind).pnt(0)) +(1.+psi)*pnts(seg(sind).pnt(1)));
			ebdry(bnum)->mvpttobdry(bseg,psi,pnts(npnt));
			lngth(npnt) = 0.5*((1.-psi)*lngth(seg(sind).pnt(0)) +(1.+psi)*lngth(seg(sind).pnt(1)));
#ifdef DEBUG_ADAPT
			*gbl->log << "Inserting boundary segment " << i << ' ' << adapt_count << ' ';
			for(int n=0;n<ND;++n)
				*gbl->log << pnts(npnt)(n) << ' ';
			*gbl->log << " seg " << bseg << " psi " << psi << std::endl;
#endif
			bdry_insert(npnt,sind,1);
			++npnt;

			/* UPDATE NEXT & PREV POINTERS */
			ebdry(bnum)->next(ebdry(bnum)->prev(bseg)) = ebdry(bnum)->nseg-1;
			ebdry(bnum)->prev(ebdry(bnum)->nseg-1) = ebdry(bnum)->prev(bseg);
			ebdry(bnum)->prev(bseg) = ebdry(bnum)->nseg-1;
			ebdry(bnum)->next(ebdry(bnum)->nseg-1) = bseg;

#ifdef DEBUG_ADAPT
			std::ostringstream nstr;
			nstr << adapt_count++ << std::flush;
			adapt_file = "adapt" +nstr.str();
			nstr.str("");
			output(adapt_file,grid);
#endif
		}
		*gbl->log << "#Slave boundary refinement finished, " << ebdry(bnum)->idnum << ' ' << (sndsize-1)/2 << " sides added" << std::endl;

	}

	return;
}


