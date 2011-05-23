/*
 *  yabers.cpp
 *  mblock
 *
 *  Created by helenbrk on Fri Sep 14 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#include "tri_mesh.h"
#include <utilities.h>
#include <assert.h>
#include <blitz/tinyvec-et.h>
#include <vector>

/* THIS IS SUPPOSED TO DO THE REVERSE OF THE REBAY ROUTINE I HOPE */
/* THUS THE NAME YABER -> REBAY */

//#define DEBUG_ADAPT

#ifdef DEBUG_ADAPT
int adapt_count;
static std::string adapt_file;
#endif

void tri_mesh::yaber(FLT tolsize) {
	int i,j,tind,sind,sind1,p0,cnt,endpt,sum;
	FLT x,y,a,asum,dx,dy,l0,l1;
	int ntsrnd, nssrnd, nperim;
	int vn,pnear,prev,tind1,stoptri,dir;
	int p1;
	int snum;
	FLT sratio;
	TinyVector<int,3> badside;
	FLT minvl;

	/* TO ADJUST FOR INSCRIBED RADIUS */
	tolsize *= 6./sqrt(3.);

	/* SET UP gbl->fltwk */
	gbl->nlst = 0;
	for(i=0;i<ntri;++i) {
		if (tri(i).info&TDLTE) continue;
		minvl = lngth(tri(i).pnt(0));
		minvl = MIN(minvl,lngth(tri(i).pnt(1)));
		minvl = MIN(minvl,lngth(tri(i).pnt(2)));
		gbl->fltwk(i) = minvl/inscribedradius(i);
		if (gbl->fltwk(i) > tolsize) putinlst(i);
	}

	/* MARK BOUNDARY VERTEX POINTS */
	/* THESE SHOULD NOT BE DELETED */
	for(i=0;i<nebd;++i) {
		for(j=0;j<ebdry(i)->nseg;++j) {
			sind = ebdry(i)->seg(j);
			p0 = seg(sind).pnt(0);
			tri(p0).info |= PSPEC;
		}
	}

	cnt = 0;
	/* BEGIN COARSENING ALGORITHM */
	while (gbl->nlst > 0) {
		for(i=gbl->nlst-1;i>=0;--i) {  // START WITH LARGEST TGT TO ACTUAL RATIO
			for(j=0;j<3;++j) {
				tind = tri(seg(i).info).tri(j);
				if (tind < 0 || pnt(tind).info == -1)  goto TFOUND;
			}
		}

		TFOUND:
		tind = seg(i).info;
		tkoutlst(tind);

		/* FIND SIDE ON TRIANGLE WITH LARGEST TARGET TO LENGTH RATIO */
		minvl = 0.0;
		sind1 = -1;
		for(j=0;j<3;++j) {
			sind = tri(tind).seg(j);
			if (seg(sind).tri(1) < 0) {
				badside(j) = false;
				continue;
			}
			sratio = MIN(lngth(seg(sind).pnt(0)),lngth(seg(sind).pnt(1)))/distance(seg(sind).pnt(0),seg(sind).pnt(1));
			if (sratio > tolsize*sqrt(3.)/6.) {
				badside(j) = true;
			}
			else {
				badside(j) = false;
			}
			if (sratio > minvl) {
				sind1 = sind;
				snum = j;
				minvl =sratio;
			}
		}
		if (sind1 < 0) {
			*gbl->log << "Trying to coarsen triangle with all three sides on boundaries\n";
			continue;
		}
		sind = sind1;

		/* REMOVE TRIANGLES THAT WILL BE DELETED */
		tind = seg(sind).tri(0);
		if (pnt(tind).info > -1) tkoutlst(tind);
		tind = seg(sind).tri(1);
		if (pnt(tind).info > -1) tkoutlst(tind);

		/* DON'T DELETE BOUNDARY POINT */
		sum = (tri(seg(sind).pnt(0)).info&PSPEC) +(tri(seg(sind).pnt(1)).info&PSPEC);
		if (sum > 0) {
			if (sum > PSPEC) {
				*gbl->log << "Trying to Delete edge with two endpoints on boundary" << pnts(seg(sind).pnt(0)) << std::endl;
				continue;
			}
			if (tri(seg(sind).pnt(0)).info&PSPEC) endpt = 1;
			else endpt = 0;
		}
		else {
			/* TRY TO FIND DIRECTION OF ACCEPTED TRIS AND DELETE AWAY FROM THAT DIRECTION */
//            if (pnt(tri(tind).tri((snum+1)%3)).info == -1  && pnt(tri(tind).tri((snum+2)%3)).info > -1) {
//                endpt = (1-tri(tind).sgn(snum))/2;
//            }
//            else if (pnt(tri(tind).tri((snum+1)%3)).info > -1  && pnt(tri(tind).tri((snum+2)%3)).info == -1) {
//                endpt = (1+tri(tind).sgn(snum))/2;
//            }
//            else {
			{
				/* KEEP POINT WHICH IS CLOSEST TO CENTER OF AREA */
				/* THIS WAY WORKS BEST BUT COSTS MORE */
				/* FIXME :: NEED TO ELIMINATE DUPLICATION OF WORK BETWEEN THIS & COLLAPSE */
				x = 0.0;
				y = 0.0;
				asum = 0.0;
				for(endpt=0;endpt<2;++endpt) {
					pnear = seg(sind).pnt(endpt);
					tind = pnt(pnear).tri;
					if (tind != seg(sind).tri(0) && tind != seg(sind).tri(1))
						prev = 0;
					else
						prev = 1;
					stoptri = tind;
					dir = 1;
					ntsrnd = 0;
					nssrnd = 0;
					nperim = 0;
					do {
						for(vn=0;vn<3;++vn)
							if (tri(tind).pnt(vn) == pnear) break;

						tind1 = tri(tind).tri((vn +dir)%3);
						if (tind1 < 0) {
							if (dir > 1) break;
							/* REVERSE DIRECTION AND GO BACK TO START */
							++dir;
							tind1 = pnt(pnear).tri;
							prev = 1;
							stoptri = -1;
						}

						if (tind1 != seg(sind).tri(0) && tind1 != seg(sind).tri(1)) {
							gbl->i2wk_lst1(ntsrnd++) = tind1;
							if (!prev) {
								gbl->i2wk_lst2(nssrnd++) = tri(tind).seg((vn +dir)%3);
							}
							prev = 0;
						}
						else {
							prev = 1;
						}

						tind = tind1;

					} while(tind != stoptri);

					tind = seg(sind).tri(endpt);
					if (tind > -1) {
						a = area(tind);
						asum += a;
						for(vn=0;vn<3;++vn) {
							x += a*pnts(tri(tind).pnt(vn))(0);
							y += a*pnts(tri(tind).pnt(vn))(1);
						}
					}
					for(j=0;j<ntsrnd;++j) {
						tind = gbl->i2wk_lst1(j);
						a = area(tind);
						asum += a;
						for(vn=0;vn<3;++vn) {
							x += a*pnts(tri(tind).pnt(vn))(0);
							y += a*pnts(tri(tind).pnt(vn))(1);
						}
					}
				}

				asum = 1./(3.*asum);
				x = x*asum;
				y = y*asum;
				p0 = seg(sind).pnt(0);
				p1 = seg(sind).pnt(1);
				dx = pnts(p0)(0) -x;
				dy = pnts(p0)(1) -y;
				l0 = dx*dx +dy*dy;
				dx = pnts(p1)(0) -x;
				dy = pnts(p1)(1) -y;
				l1 = dx*dx +dy*dy;

				endpt = (l0 > l1 ? 0 : 1);
			}
		}

#ifdef DEBUG_ADAPT
		std::cout << "collapsing interior" << cnt << ' ' << adapt_count << ' ' << sind << " endpt " << endpt << std::endl;
#endif

		/* COLLAPSE EDGE */
		collapse(sind,endpt);
		++cnt;


		/* RECLASSIFY AFFECTED TRIANGLES */
		for(i=0;i<gbl->i2wk_lst1(-1);++i) {
			tind = gbl->i2wk_lst1(i);
			if (pnt(tind).info > -1) tkoutlst(tind);
			minvl = lngth(tri(tind).pnt(0));
			minvl = MIN(minvl,lngth(tri(tind).pnt(1)));
			minvl = MIN(minvl,lngth(tri(tind).pnt(2)));
			gbl->fltwk(tind) = minvl/inscribedradius(tind);
			if (gbl->fltwk(tind) > tolsize) putinlst(tind);
		}

#ifdef DEBUG_ADAPT
		std::ostringstream nstr;
		nstr << adapt_count++ << std::flush;
		adapt_file = "adapt" +nstr.str() + "_" +gbl->idprefix;
		nstr.str("");
		output(adapt_file.c_str(),grid);
#endif
	}

	*gbl->log << "#Yaber finished: " << cnt << " sides coarsened" << std::endl;

	return;
}

void tri_mesh::checkintegrity() {
	int i,j,sind,dir;

	for(i=0;i<maxpst;++i) {
		if (gbl->intwk(i) > -1) {
			*gbl->log << "gbl->intwk check failed" << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
	}
	
	for(i=0;i<npnt;++i) {
		int tind = pnt(i).tri;
		for(j=0;j<3;++j)
			if (tri(tind).pnt(j) == i) goto next;
			
		*gbl->log << "tri.pnt is out of whack" <<  i << tind << std::endl;
		output("error");
		sim::abort(__LINE__,__FILE__,gbl->log);
		
		next: continue;
	}

	for(i=0;i<ntri;++i) {
		if (tri(i).info < 0) continue;

		if (area(i) < 0.0) *gbl->log << "negative area" << i << std::endl;

		for(j=0;j<3;++j) {
			sind = tri(i).seg(j);
			dir = -(tri(i).sgn(j) -1)/2;

//			if (seg(sind).info == -3) {
//				*gbl->log << "references deleted segment" <<  i << sind << std::endl;
//				for(i=0;i<nseg;++i)
//					seg(i).info += 2;
//				output("error");
//				sim::abort(__LINE__,__FILE__,gbl->log);
//			}

			if (seg(sind).pnt(dir) != tri(i).pnt((j+1)%3) && seg(sind).pnt(1-dir) != tri(i).pnt((j+2)%3)) {
				*gbl->log << "failed pnt check tind" << i << "sind" << sind << std::endl;
				for(i=0;i<nseg;++i)
					seg(i).info += 2;
				output("error");
				sim::abort(__LINE__,__FILE__,gbl->log);
			}

			if (seg(sind).tri(dir) != i) {
				*gbl->log << "failed segment check tind" << i << "sind" << sind << std::endl;
				for(i=0;i<nseg;++i)
					seg(i).info += 2;
				output("error");
				sim::abort(__LINE__,__FILE__,gbl->log);
			}

			if (tri(i).tri(j) != seg(sind).tri(1-dir)) {
				*gbl->log << "failed ttri check tind" << i << "sind" << sind << std::endl;
				for(i=0;i<nseg;++i)
					seg(i).info += 2;
				output("error");
				sim::abort(__LINE__,__FILE__,gbl->log);
			}

//			if (tri(i).tri(j) > 0) {
//				if(tri(tri(i).tri(j)).info < 0 && tri(tri(i).tri(j)).info != -2) {
//					*gbl->log << "triangle " << i << " references deleted tri " << tri(i).tri(j) << std::endl;
//					for(i=0;i<nseg;++i)
//						seg(i).info += 2;
//					output("error");
//					sim::abort(__LINE__,__FILE__,gbl->log);
//				}
//			}
		}
	}

	return;
}

void tri_mesh::bdry_yaber(FLT tolsize) {
	int sind,endpt,p0,p1,count;
	int bseg, pel, nel, sindprev, sindnext, saffect;

	/* COARSEN FIRST BOUNDARIES */
	for(int bnum=0;bnum<nebd;++bnum) {
		count = 0;

		if (!ebdry(bnum)->is_frst()) {
			ebdry(bnum)->sndtype() = boundary::int_msg;
			ebdry(bnum)->comm_prepare(boundary::all,0,boundary::master_slave);
			continue;
		}

		gbl->nlst = 0;
		for(int indx=0;indx<ebdry(bnum)->nseg;++indx) {
			sind = ebdry(bnum)->seg(indx);
			if (tri(sind).info&SDLTE) continue;
			gbl->fltwk(sind) = MIN(lngth(seg(sind).pnt(0)),lngth(seg(sind).pnt(1)))/distance(seg(sind).pnt(0),seg(sind).pnt(1));
			if (gbl->fltwk(sind) > tolsize) {
				putinlst(sind);
			}
		}

		/* SKIP FIRST SPOT SO CAN SEND LENGTH FIRST */
		ebdry(bnum)->sndsize() = 1;
		while (gbl->nlst > 0) {
			// START WITH LARGEST SIDE LENGTH RATIO
			sind = seg(gbl->nlst-1).info;
			bseg = getbdryseg(seg(sind).tri(1));

			/* ADJACENT SIDES */
			nel = ebdry(bnum)->next(bseg);
			pel = ebdry(bnum)->prev(bseg);
			sindnext = ebdry(bnum)->seg(nel);
			sindprev = ebdry(bnum)->seg(pel);

			/* PICK ENDPT */
			p0 = seg(sind).pnt(0);
			p1 = seg(sind).pnt(1);

			if (tri(p0).info&PSPEC && tri(p1).info&PSPEC) {
				tkoutlst(sind);
				continue;
			}
			else if (tri(p0).info&PSPEC) {
				endpt = 1;
				saffect = sindnext;
			}
			else if (tri(p1).info&PSPEC) {
				endpt = 0;
				saffect = sindprev;
			}
			else {
				/* PICK MORE ARPROPRIATE VERTEX TO DELETE */
				if (gbl->fltwk(sindprev) > gbl->fltwk(sindnext)) {
					endpt = 0;
					saffect = sindprev;
				}
				else {
					endpt = 1;
					saffect = sindnext;
				}
			}
#ifdef DEBUG_ADAPT
			std::cout << "collapsing boundary"  << adapt_count << ' ' << sind << " endpt " << endpt << std::endl;
#endif
			collapse(sind,endpt);
			tkoutlst(sind);

			/* UPDATE PREV/NEXT BOUNDARY POINTERS */
			ebdry(bnum)->prev(nel) = pel;
			ebdry(bnum)->next(pel) = nel;

			/* STORE ADAPTATION INFO FOR COMMUNICATION BOUNDARIES */
			ebdry(bnum)->isndbuf(ebdry(bnum)->sndsize()++) = bseg;
			ebdry(bnum)->isndbuf(ebdry(bnum)->sndsize()++) = endpt;

			/* UPDATE AFFECTED SIDE */
			if (pnt(saffect).info > -1) tkoutlst(saffect);
			gbl->fltwk(saffect) = MIN(lngth(seg(saffect).pnt(0)),lngth(seg(saffect).pnt(1)))/distance(seg(saffect).pnt(0),seg(saffect).pnt(1));
			if (gbl->fltwk(saffect) > tolsize) putinlst(saffect);
			
			/* Check for inverted interior triangles because of changing boundary geometry */
			int ntsrnd = gbl->i2wk_lst1(-1);
			std::vector<int> badpnt;
			int nbadpnt = 0;
			for(int i=0;i<ntsrnd;++i) {
				int tind = gbl->i2wk_lst1(i);
				if (area(tind) > 0.0) continue;
				
				for(int vn=0;vn<3;++vn) {
					if (tri(tind).pnt(vn) == seg(saffect).pnt(0) || tri(tind).pnt(vn) == seg(saffect).pnt(1)) continue;
					if (area(saffect,tri(tind).pnt(vn)) < FLT_EPSILON) {
						badpnt[nbadpnt++] = tri(tind).pnt(vn);
					}
				}
			}
			
			for(int i=0;i<nbadpnt;++i) {
				*gbl->log << "Warning deleting " << badpnt[i] << " because of boundary change to " << saffect << '\n';
				int tind = pnt(badpnt[i]).tri;
				int vn;
				for(vn=0;vn<3;++vn)  {
					if (tri(tind).pnt(vn) == badpnt[i])
						break;
				}
				collapse(tri(tind).seg((vn+1)%3),1);
			}
			
			++count;
#ifdef DEBUG_ADAPT
			std::ostringstream nstr;
			nstr << adapt_count++ << std::flush;
			adapt_file = "adapt" +nstr.str() + "_" +gbl->idprefix;
			nstr.str("");
			output(adapt_file.c_str(),grid);
#endif
		}
		ebdry(bnum)->isndbuf(0) = ebdry(bnum)->sndsize();
		ebdry(bnum)->sndtype() = boundary::int_msg;
		ebdry(bnum)->comm_prepare(boundary::all,0,boundary::master_slave);
		*gbl->log << "#Boundary coarsening finished, " << ebdry(bnum)->idnum << ' ' << count << " sides coarsened" << std::endl;
	}

	return;
}

void tri_mesh::bdry_yaber1() {
	int i,lseg,endpt,sind,sndsize;

	for(int bnum=0;bnum<nebd;++bnum) {

		ebdry(bnum)->comm_wait(boundary::all,0,boundary::master_slave);

		if (ebdry(bnum)->is_frst() || !ebdry(bnum)->is_comm()) continue;

		sndsize = ebdry(bnum)->ircvbuf(0,0);

		for(i=1;i<sndsize;i+=2) {
			lseg = ebdry(bnum)->nseg -1 -ebdry(bnum)->ircvbuf(0,i);
			endpt = 1 -ebdry(bnum)->ircvbuf(0,i+1);
			sind = ebdry(bnum)->seg(lseg);
#ifdef DEBUG_ADAPT
			std::cout << "collapsing boundary" << adapt_count << ' ' << sind << " endpt " << endpt << std::endl;
#endif
			collapse(sind,endpt);

			/* UPDATE PREV/NEXT BOUNDARY POINTERS */
			int nel = ebdry(bnum)->next(lseg);
			int pel = ebdry(bnum)->prev(lseg);
			ebdry(bnum)->prev(nel) = pel;
			ebdry(bnum)->next(pel) = nel;
			
			int saffect;
			if (endpt == 0) 
				saffect = ebdry(bnum)->seg(pel);
			else
				saffect = ebdry(bnum)->seg(nel);
				
			/* Check for inverted interior triangles because of changing boundary geometry */
			int ntsrnd = gbl->i2wk_lst1(-1);
			std::vector<int> badpnt;
			int nbadpnt = 0;
			for(int i=0;i<ntsrnd;++i) {
				int tind = gbl->i2wk_lst1(i);
				if (area(tind) > 0.0) continue;
				for(int vn=0;vn<3;++vn) {
					if (tri(tind).pnt(vn) == seg(saffect).pnt(0) || tri(tind).pnt(vn) == seg(saffect).pnt(1)) continue;
					if (area(saffect,tri(tind).pnt(vn)) < FLT_EPSILON) {
						badpnt[nbadpnt++] = tri(tind).pnt(vn);
					}
				}
			}
			
			for(int i=0;i<nbadpnt;++i) {
				*gbl->log << "Warning deleting " << badpnt[i] << " because of boundary change to " << saffect << '\n';
				int tind = pnt(badpnt[i]).tri;
				int vn;
				for(vn=0;vn<3;++vn)  {
					if (tri(tind).pnt(vn) == badpnt[i])
						break;
				}
				collapse(tri(tind).seg((vn+1)%3),1);
			}

#ifdef DEBUG_ADAPT
			std::ostringstream nstr;
			nstr << adapt_count++ << std::flush;
			adapt_file = gbl->idprefix +"_adapt" +nstr.str();
			nstr.str("");
			output(adapt_file.c_str(),grid);
#endif
		}
		*gbl->log << "#Slave Boundary coarsening finished, " << ebdry(bnum)->idnum << ' ' << (sndsize-1)/2 << " sides coarsened" << std::endl;
	}
	return;
}


void tri_mesh::checkintwk() const {
	int i;

	for(i=0;i<gbl->intwk.extent(firstDim)-1;++i)
		if (gbl->intwk(i) != -1) *gbl->log << "failed gbl->intwk check " << i << ' ' << gbl->intwk(i) << std::endl;

	return;
}
