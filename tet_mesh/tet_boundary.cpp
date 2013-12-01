/*
 *  boundaries.cpp
 *  mesh3d
 *
 *  Created by Michael Brazell on 7/31/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_mesh.h"
#include "tet_boundary.h"
#include <assert.h>
#include <float.h>
//#include <blitz/tinyvec-et.h>

/********************/
/* VERTEX FUNCTIONS */
/********************/
	
/* GENERIC VERTEX COMMUNICATIONS */
void vcomm::ploadbuff(boundary::groups grp,FLT *base,int bgn,int end, int stride) {
	int i,offset;

	if (!in_group(grp)) return;

	sndsize()=end-bgn+1;
	sndtype()=flt_msg;
	
	/* LOAD SEND BUFFER */    
	offset = pnt*stride +bgn;
	for (i=0;i<end-bgn+1;++i) 
		fsndbuf(i) = base[offset+i];
}

void vcomm::pfinalrcv(boundary::groups grp, int phi, comm_type type, operation op, FLT *base, int bgn, int end, int stride) {
	int i,offset;
		
	bool reload = comm_finish(grp,phi,type,op);
	if (!reload) return;
	
	offset = pnt*stride +bgn;
	for(i=0;i<end-bgn+1;++i) 
		base[offset++] = fsndbuf(i);
		
	return;
}

/**************************************/
/* GENERIC FUNCTIONS FOR EDGES        */
/**************************************/
void edge_bdry::alloc(int n) {
	maxseg = n;
	seg.resize(Range(-1,n));
}

void edge_bdry::copy(const edge_bdry& bin) {
	int i;
		
	if (!maxseg) alloc(bin.maxseg);
	else assert(bin.nseg <= maxseg);
	vbdry = bin.vbdry;
	
	nseg = bin.nseg;
	
	for(i=0;i<nseg;++i)
		seg(i) = bin.seg(i);
		
	return;
}

void edge_bdry::mvpttobdry(int indx, FLT psi, TinyVector<FLT,tet_mesh::ND> &pt) {
	/* FOR A LINEAR SIDE */
	int n;
	
	for (n=0;n<tet_mesh::ND;++n)
		pt(n) = 0.5*((1. -psi)*x.pnts(x.seg(seg(indx).gindx).pnt(0))(n) +(1.+psi)*x.pnts(x.seg(seg(indx).gindx).pnt(1))(n));
	
	return;
}

void edge_bdry::findbdrypt(const TinyVector<FLT,tet_mesh::ND> xpt, int &sidloc, FLT &psiloc) const {
	int k,sind,p0,p1,sidlocprev;
	FLT ol,psi,normdist;
	FLT psiprev,normdistprev;
	FLT mindist = 1.0e32;
	TinyVector<FLT,3> dx,normvect;
		
	if (x.seg(seg(0).gindx).pnt(0) == x.seg(seg(nseg-1).gindx).pnt(1)) {
		/* BOUNDARY IS A LOOP */
		sind = seg(nseg-1).gindx;
		p0 = x.seg(sind).pnt(0);
		p1 = x.seg(sind).pnt(1);
		dx = x.pnts(p1) -x.pnts(p0);
		ol = 2./dot(dx,dx);
		psi = ol*dot(xpt-x.pnts(p0),dx) -1.;
		normvect = xpt -x.pnts(p0) -dx*(1.+psi)*0.5;
		normdist = dot(normvect,normvect);
		normdist *= sqrt(ol*0.5);
		psiprev = psi;
		normdistprev = normdist;
		sidlocprev = nseg-1;
	} 
	else {
		psiprev = -1.0;
		sidlocprev = -1;  // FIXME: ENDPOINTS ARE NOT SEARCHED
		normdistprev = 1.0e10; // This shouldn't be used
	}
	
	for(k=0;k<nseg;++k) {
		sind = seg(k).gindx;
		p0 = x.seg(sind).pnt(0);
		p1 = x.seg(sind).pnt(1);
		dx = x.pnts(p1) -x.pnts(p0);
		ol = 2./dot(dx,dx);
		psi = ol*dot(xpt-x.pnts(p0),dx) -1.;
		normvect = xpt -x.pnts(p0) -dx*(1.+psi)*0.5;
		normdist = dot(normvect,normvect);
		normdist *= sqrt(ol*0.5);

		if (psi <= -1.0 && psiprev >= 1.0) {
			/* PREVIOUS & THIS SIDE ARE POTENTIAL MATCHES */
			if (fabs(normdist) < mindist) {
				mindist = fabs(normdist);
				sidloc = k;
				psiloc = -1.0;
			}
			if (fabs(normdistprev) < mindist) {
				mindist = fabs(normdistprev);
				sidloc = sidlocprev;
				psiloc = 1.0;
			}
		}
		else if (psi >= -1.0 && psi <= 1.0) {
			/* POTENTIAL SIDE */
			if (fabs(normdist) < mindist) {
				mindist = fabs(normdist);
				sidloc = k;
				psiloc = psi;
			}
		}
		psiprev = psi;
		normdistprev = normdist;
		sidlocprev = k;
	}
	
	return;
}

void edge_bdry::mgconnect(Array<tet_mesh::transfer,1> &cnnct,tet_mesh& tgt, int bnum) {
	int j,k,sind,ttind,p0,sidloc;
	FLT psiloc;
	
	TinyMatrix<int,6,2> vs;    
	vs(0,0)=2, vs(0,1)=3; 
	vs(1,0)=3, vs(1,1)=1; 
	vs(2,0)=2, vs(2,1)=1;
	vs(3,0)=1, vs(3,1)=0;
	vs(4,0)=2, vs(4,1)=0; 
	vs(5,0)=3, vs(5,1)=0;
		
	for(k=1;k<nseg;++k) {
		p0 = x.seg(seg(k).gindx).pnt(0);
		tgt.ebdry(bnum)->findbdrypt(x.pnts(p0), sidloc, psiloc);
		sind = tgt.ebdry(bnum)->seg(sidloc).gindx;
		ttind = tgt.seg(sind).tet;                            
		cnnct(p0).tet = ttind;
		for (j=0;j<6;++j) 
			if (tgt.tet(ttind).seg(j) == sind) break;
		assert(j < 6);
		
		psiloc = psiloc*tgt.tet(ttind).sgn(j);
		cnnct(p0).wt = 0.0;
		cnnct(p0).wt(vs(j,0)) = 0.5*(1.-psiloc);
		cnnct(p0).wt(vs(j,1)) = 0.5*(1.+psiloc);
	} 
	
	return;
}

/* SWAP ELEMENTS IN LIST */
void edge_bdry::swap(int s1, int s2) {
	segstruct ind;
		
	 if (s1 == s2) return;
	 
	/* FIXME NOT SURE HOW TO SWAP S VALUES */    
	ind = seg(s1);
	seg(s1) = seg(s2);
	seg(s2) = ind;
	
	/* PREV/NEXT IMPLEMENTATION */
	int prev1 = seg(s1).prev;
	int next1 = seg(s1).next;
	int prev2 = seg(s2).prev;
	int next2 = seg(s2).next;
	
	if (prev1 == s1) prev1 = s2;
	if (next1 == s1) next1 = s2;
	if (prev2 == s2) prev2 = s1;
	if (next2 == s2) next2 = s1;
	
	seg(prev1).next = s1;	 
	seg(next1).prev = s1;    
	 
	seg(prev2).next = s2;
	seg(next2).prev = s2;

	return;
}

void edge_bdry::setup_next_prev() {
	int i,next,prev,sind,sind2;

	/* DON'T ASSUME wk INITIALIZED TO -1 */
	for(i=0;i<nseg;++i) {
		sind = seg(i).gindx;
		seg(i).prev = -1;
		seg(i).next = -1;
		x.gbl->i1wk(x.seg(sind).pnt(0)) = -1;
		x.gbl->i1wk(x.seg(sind).pnt(1)) = -1;
	}
	
	/* FILL IN NEXT/PREV DATA */
	for(i=0; i < nseg; ++i) {
		sind = seg(i).gindx;
		int v0 = x.seg(sind).pnt(0);
		if (x.gbl->i1wk(v0) == -1) {
			x.gbl->i1wk(v0) = i;
		}
		else {
			/* found match */
			prev = x.gbl->i1wk(v0);
			seg(i).prev = prev;
			sind2 = seg(prev).gindx;
			if (x.seg(sind2).pnt(0) == v0)
				seg(prev).prev = i;
			else
				seg(prev).next = i;
		}

		int v1 = x.seg(sind).pnt(1);
		if (x.gbl->i1wk(v1) == -1) {
		   x.gbl->i1wk(v1) = i;
		}
		else {
			/* found match */
			next = x.gbl->i1wk(v1);
			seg(i).next = next;
			sind2 = seg(next).gindx;
			if (x.seg(sind2).pnt(1) == v1)
				seg(next).next = i;
			else
				seg(next).prev = i;
		}
	 }
	 
	 /* RESET gbl->i1wk TO -1 */
	for(i=0; i <nseg; ++i) {
		sind = seg(i).gindx;
		x.gbl->i1wk(x.seg(sind).pnt(1)) = -1;
		x.gbl->i1wk(x.seg(sind).pnt(0)) = -1;
	}
	
	return;
}


/* REORDERS BOUNDARIES TO BE SEQUENTIAL & REORIENTS EDGES TO ALL BE ALIGNED IN THE SAME DIRECTION */
/* USES gbl->i1wk & gbl->i2wk AS WORK ARRAYS */
void edge_bdry::reorder() {
	int i,count,next,sind,first;
	bool loop = false;
	
	/* FIND FIRST SIDE */    
	first = -1;
	for(i=0;i<nseg;++i) {
		if (seg(i).prev < 0 || seg(i).next < 0) {
			first = i;
			break;
		}
	}
	 
	 if (first < 0) {
		 /* EDGE LOOP */
		 first = 0;
		 loop = true;
	 }
	 else {
	    if (seg(first).prev != -1) {
		    /* Reverse orientation of first side */
			 seg(first).next = seg(first).prev;
			 seg(first).prev = -1;
			 sind = seg(first).gindx;
			 int v0 = x.seg(sind).pnt(0);
			 x.seg(sind).pnt(0) = x.seg(sind).pnt(1);
			 x.seg(sind).pnt(1) = v0;
			 //x.switch_edge_sign(sind);
		 }
	 }
	 
	/* First swap directions, then reorder */
	count = 0;
	int indx = first;	
	while(count < nseg && (next = seg(indx).next) > -1) {
		if (seg(next).prev != indx) {
			/* Reverse orientation side */
			seg(next).next = seg(next).prev;
			seg(next).prev = indx;
			sind = seg(next).gindx;
			int v0 = x.seg(sind).pnt(0);
			x.seg(sind).pnt(0) = x.seg(sind).pnt(1);
			x.seg(sind).pnt(1) = v0;
			//x.switch_edge_sign(sind);
		}
		++count;
		indx = next;
	}
	 
	/* insert -1 so while loop exits for an edge loop */
	if (loop == true) seg(seg(first).prev).next = -1;
	
	/* now reorder */
	count = 0;
	indx = first;	
	do {
		swap(count++,indx);
	} while (count < nseg && (indx = seg(count-1).next) > -1);    
	
	/* separate disconnected edge boundaries */
	 if (count < nseg) {
		++x.nebd;
		x.ebdry.resizeAndPreserve(x.nebd);
		x.ebdry(x.nebd-1) = create(x);
		x.ebdry(x.nebd-1)->copy(*this);
		for(i=0;i<nseg-count;++i)
			x.ebdry(x.nebd-1)->swap(i,i+count);
		x.ebdry(x.nebd-1)->nseg = nseg -count;
		*x.gbl->log << "#creating new " << mytype << " edge boundary: " << idnum << " num: " << x.ebdry(x.nebd-1)->nseg << std::endl;
		nseg = count;
	 }
	
	return;
}

void ecomm::match_numbering(int step) {
	int sind,sind2;
	TinyVector<FLT,tet_mesh::ND> mpnt1,mpnt2;

	switch(step) {
			
		case(0): {
			if (is_frst()) {
				return;
			}
			else {
				/* Slave receives master boundaries vertex positions */
				int count = 0;
				bool printstuff = false;	
				
				Array<int,1> gindx(nseg);
				
				for (int i=0;i<nseg;++i)
					gindx(i) = seg(i).gindx;

				for (int n=0;n<tet_mesh::ND;++n)
					mpnt2(n) = fsndbuf(count++);
				
				for (int i=0;i<nseg;++i) {
					mpnt1 = mpnt2;
					for (int n=0;n<tet_mesh::ND;++n)
						mpnt2(n) = fsndbuf(count++);

					/* brute force */
					FLT dist = 0.0;
					seg(i).gindx = -1;
					for (int j=0; j<nseg; ++j) {
						TinyVector<FLT,tet_mesh::ND> lclpnts1,lclpnts2;
						lclpnts1 = x.pnts(x.seg(gindx(j)).pnt(0));
						lclpnts2 = x.pnts(x.seg(gindx(j)).pnt(1));
						
						dist = 0.0;
						for (int k=0; k<tet_mesh::ND; ++k) {
							dist += pow(lclpnts1(k)-mpnt1(k), 2.0);
							dist += pow(lclpnts2(k)-mpnt2(k), 2.0);

						}

						if (dist <  10.*EPSILON) {
							seg(i).gindx = gindx(j);
							break;
						}
						
						dist = 0.0;
						for (int k=0; k<tet_mesh::ND; ++k) {
							dist += pow(lclpnts2(k)-mpnt1(k), 2.0);
							dist += pow(lclpnts1(k)-mpnt2(k), 2.0);
						}

						if (dist < 10.*EPSILON) {
							seg(i).gindx = gindx(j);
							break;
						}
						
					}
					
					if (dist > 10.*EPSILON) {
						printstuff = true;
						*x.gbl->log << "Matching edge numbering error, edge: " << idprefix << endl;						
						//*x.gbl->log << "Matching face numbering error: " << dist << ' ' << mpnt << ' ' << x.pnts(pnt(i).gindx) << '\n';
						//*x.gbl->log << "idnum " << idnum << " lcl point " << i << " gbl point " << pnt(i).gindx << '\n';
					}					
				}
				setup_next_prev();
				reorder();
				
				/* check to make sure pnts locations match perfectly */
				if (printstuff) {
					count = 0;
					for (int i=0; i<nseg; ++i) {
						for (int n=0; n<3; ++n) {
							*x.gbl->log << fsndbuf(count++)-x.pnts(x.seg(seg(i).gindx).pnt(0))(n) << endl;
						}
					}
					for (int n=0; n<3; ++n) 
						*x.gbl->log << fsndbuf(count++)-x.pnts(x.seg(seg(nseg-1).gindx).pnt(1))(n) << endl;
					
					exit(4);
				}
				
				
			}
			return;
		}
		case 1: {
			if (is_frst()) {
				sind = seg(0).gindx;
				isndbuf(0) = x.pnt(x.seg(sind).pnt(0)).info;
				isndbuf(1) = x.pnt(x.seg(sind).pnt(1)).info;
			}
			sndsize() = 2;
			sndtype() = int_msg;
			break;
		}
		case 2: {

			if (is_frst()) return;
			
			*x.gbl->log << "match numbering for edge: " << idprefix << "\n";

			sind = seg(0).gindx;
			if (isndbuf(0) == x.pnt(x.seg(sind).pnt(0)).info && isndbuf(1) == x.pnt(x.seg(sind).pnt(1)).info) {
				*x.gbl->log << "code returned 0 \n";
				return;
			}
				
			
			/* EDGE IS ORIENTED BACKWARDS FROM MASTER OR IS MISALIGNED LOOP OR BOTH */
			sind = seg(nseg-1).gindx;
			if (isndbuf(0) == x.pnt(x.seg(sind).pnt(1)).info && isndbuf(1) == x.pnt(x.seg(sind).pnt(0)).info) {
				/* just backwards, swap first & last, then reorder */
				swap(0,nseg-1);
				reorder();
				*x.gbl->log << "code returned 1 \n";

				return;
			}
			
			/* correct order but endpoints flipped */
			sind = seg(0).gindx;
			sind2 = seg(nseg-1).gindx;
			if (isndbuf(0) == x.pnt(x.seg(sind2).pnt(1)).info && isndbuf(1) == x.pnt(x.seg(sind).pnt(1)).info) {
				swap(0,nseg-1);
				setup_next_prev();
				reorder();
				*x.gbl->log << "code returned 2 \n";
				return;
			}

			/* reverse order and endpoints flipped */
			if (isndbuf(0) == x.pnt(x.seg(sind2).pnt(0)).info && isndbuf(1) == x.pnt(x.seg(sind).pnt(0)).info) {
				
//				swap(0,nseg-1);
//				
//				/* invert side direction */
//				int temp = seg(0).next;
//				seg(0).next = seg(0).prev;
//				seg(0).prev = temp;
//				sind = seg(0).gindx;
//				int v0 = x.seg(sind).pnt(0);
//				x.seg(sind).pnt(0) = x.seg(sind).pnt(1);
//				x.seg(sind).pnt(1) = v0;
				setup_next_prev();
				reorder();
				*x.gbl->log << "code returned 3 \n";
				
				return;
			}

			/* Find Start of Loop */
			for (int i=0;i<nseg;++i) {
				sind = seg(i).gindx;
				if (isndbuf(0) == x.pnt(x.seg(sind).pnt(0)).info && isndbuf(1) == x.pnt(x.seg(sind).pnt(1)).info) {
					swap(0,i);
					reorder();
					*x.gbl->log << "code returned 4 \n";

					return;
				}
				
				if (isndbuf(0) == x.pnt(x.seg(sind).pnt(1)).info && isndbuf(1) == x.pnt(x.seg(sind).pnt(0)).info) {
					swap(0,i);
					/* invert side direction */
					int temp = seg(0).next;
					seg(0).next = seg(0).prev;
					seg(0).prev = temp;
					sind = seg(0).gindx;
					int v0 = x.seg(sind).pnt(0);
					x.seg(sind).pnt(0) = x.seg(sind).pnt(1);
					x.seg(sind).pnt(1) = v0;
					
					reorder();
					*x.gbl->log << "code returned 5 \n";

					return;
				}
			}
		}
			
		*x.gbl->log << " Didn't find match for edge? " << idprefix;
		*x.gbl->log << " Trying to find " << isndbuf(0) << ' ' << isndbuf(1) << std::endl;
		*x.gbl->log << " These are my options " << std::endl;
		for (int i=0;i<nseg;++i) {
			sind = seg(i).gindx;
			*x.gbl->log << "global edge " << sind << "  pnts " << x.pnt(x.seg(sind).pnt(0)).info  << ' '<< x.pnt(x.seg(sind).pnt(1)).info <<  std::endl;
		}
		for (int i=0;i<nseg;++i) {
			sind = seg(i).gindx;
			for (int n=0; n<3; ++n) 
				*x.gbl->log << x.pnts(x.seg(sind).pnt(0))(n) << ", ";

			*x.gbl->log << "\n";
		}
		for (int n=0; n<3; ++n) 
			*x.gbl->log << x.pnts(x.seg(seg(nseg-1).gindx).pnt(1))(n) <<  ", ";
		*x.gbl->log << "\n";
			
		exit(9);

			
	 }
	

	return;
}
	
void ecomm::ploadbuff(boundary::groups grp,FLT *base,int bgn,int end, int stride) {
	int j,k,count,offset;
	int sind = 0;  // To avoid may be used uninitialized warnings
	
	if (!in_group(grp)) return;

	count = 0;
	for(j=0;j<nseg;++j) {
		sind = seg(j).gindx;
		offset = x.seg(sind).pnt(0)*stride;
		for (k=bgn;k<=end;++k) {
			 fsndbuf(count++) = base[offset+k];
		}
	}
	offset = x.seg(sind).pnt(1)*stride;
	for (k=bgn;k<=end;++k) 
		fsndbuf(count++) = base[offset+k]; 

		
	sndsize() = count;
	sndtype() = boundary::flt_msg;
}

void ecomm::pfinalrcv(boundary::groups grp, int phi, comm_type type, operation op, FLT *base,int bgn,int end, int stride) {
	int j,k,count,offset;
	int sind = 0;  // To avoid may be used uninitialized warnings

	
	bool reload = comm_finish(grp,phi,type,op);
	if (!reload) return;

	int ebp1 = end-bgn+1;
	count = 0;
	for(j=0;j<nseg;++j) {
		sind = seg(j).gindx;
		offset = x.seg(sind).pnt(0)*stride +bgn;
		for(k=0;k<ebp1;++k) {
			 base[offset+k] = fsndbuf(count++);
		}
	}
	offset = x.seg(sind).pnt(1)*stride +bgn;
	for(k=0;k<ebp1;++k) {
		base[offset+k] = fsndbuf(count++);
	}
		
	return;
}
	
void ecomm::sloadbuff(boundary::groups grp,FLT *base,int bgn,int end, int stride) {
	int j,k,count,sind,offset;
	
	if (!in_group(grp)) return;

	  count = 0;
	  for(j=0;j<nseg;++j) {
			sind = seg(j).gindx;
			offset = sind*stride;
			for (k=bgn;k<=end;++k) {
				 fsndbuf(count++) = base[offset+k];
			}
	  }

		
	sndsize() = count;
	sndtype() = boundary::flt_msg;
}

void ecomm::sfinalrcv(boundary::groups grp, int phi, comm_type type, operation op, FLT *base,int bgn,int end, int stride) {
	int j,k,count,offset,sind;
	
	/* ASSUMES REVERSE ORDERING OF SIDES */
	bool reload = comm_finish(grp,phi,type,op);
	if (!reload) return;

	  count = 0;
	  for(j=0;j<nseg;++j) {
			sind = seg(j).gindx;
			offset = sind*stride;
			for (k=bgn;k<=end;++k) {
				 base[offset+k] = fsndbuf(count++);
			}
	  }
 }

void epartition::mgconnect(Array<tet_mesh::transfer,1> &cnnct,tet_mesh& tgt, int bnum) {
	int i,j,k,p0;
	

	/* BOUNDARY IS AN INTERNAL PARTITION BOUNDARY */
	/* MAKE SURE ENDPOINTS ARE OK */
	i = x.seg(seg(0).gindx).pnt(0);
	if (cnnct(i).tet < 0) {
		tgt.otree.nearpt(x.pnts(i).data(),p0);
		cnnct(i).tet=tgt.pnt(p0).tet;
		for(j=0;j<4;++j) {
			cnnct(i).wt(j) = 0.0;
			if (tgt.tet(cnnct(i).tet).pnt(j) == p0) cnnct(i).wt(j) = 1.0;
		}
	}
	i = x.seg(seg(nseg-1).gindx).pnt(1);
	if (cnnct(i).tet < 0) {
		tgt.otree.nearpt(x.pnts(i).data(),p0);
		cnnct(i).tet=tgt.pnt(p0).tet;
		for(j=0;j<3;++j) {
			cnnct(i).wt(j) = 0.0;
			if (tgt.tet(cnnct(i).tet).pnt(j) == p0) cnnct(i).wt(j) = 1.0;
		}
	}
	
	if (first) {
		sndsize() = 0;
		sndtype() = int_msg;
		for(k=1;k<nseg;++k) {
			p0 = x.seg(seg(k).gindx).pnt(0);
			if (cnnct(p0).tet > 0) {
				isndbuf(sndsize()++) = -1;
			}
			else {
				isndbuf(sndsize()++) = +1;
				cnnct(p0).tet = 0;
				for(j=0;j<3;++j)
					cnnct(p0).wt(j) = 0.0;
			}
		}
	}
	
	comm_prepare(boundary::all,0,slave_master); 
	comm_exchange(boundary::all,0,slave_master);
	comm_wait(boundary::all,0,slave_master);
	
	if (!first) {
		i = 0;
		for(k=nseg-1;k>0;--k) {
			p0 = x.seg(seg(k).gindx).pnt(1);
			if (ircvbuf(0,i) < 0) {
				cnnct(p0).tet = 0;
				for(j=0;j<3;++j)
					cnnct(p0).wt(j) = 0.0;
			}
		}
	}                    
}




/********************/
/* FACE FUNCTIONS */
/********************/
void face_bdry::alloc(int mxsize) {

	/* SIDE INFO */
	maxpst = mxsize;
	seg.resize(Range(-1,maxpst));

	/* VERTEX INFO */                                
	pnt.resize(Range(-1,maxpst));

	/* TRI INFO */ 
	tri.resize(Range(-1,maxpst));
	
	otree.allocate(x.pnts, mxsize);

	return;
}

void face_bdry::copy(const face_bdry& tgt) {
	int i;
		
	if (!maxpst) {
		alloc(tgt.maxpst);
	}
	else {
		/* CHECK IF BIG ENOUGH */
		if (tgt.nseg > maxpst) {
			*x.gbl->log << "face_bdry is too big to copy" << std::endl;
			exit(1);
		}
	}
	
	/* COPY VERTEX INFO OVER */
	npnt = tgt.npnt;
	for(i=0;i<npnt;++i)
		pnt(i) = tgt.pnt(i);
			
	/* COPY SIDE INFORMATION */
	nseg = tgt.nseg;
	for(i=0;i<nseg;++i) {
		seg(i) = tgt.seg(i);
	}

	nebd = tgt.nebd;
	ebdry = tgt.ebdry;
	
	/* COPY ELEMENT DATA */
	ntri = tgt.ntri;
	for(i=0;i<ntri;++i) {
		tri(i) = tgt.tri(i);
	}
		
	otree.copy(tgt.otree);
	otree.change_vptr(tgt.x.pnts);
	
	return;  
}

void face_bdry::mgconnect(Array<tet_mesh::transfer,1> &cnnct,tet_mesh& tgt, int bnum) {
	int j,k,tind,ttind,p0,facloc;
	FLT r, s;
	
	// vertices of each face on a tet
	TinyMatrix<int,4,3> vf1;
	vf1(0,0)=1, vf1(0,1)=2, vf1(0,2)=3;
	vf1(1,0)=0, vf1(1,1)=3, vf1(1,2)=2;
	vf1(2,0)=0, vf1(2,1)=1, vf1(2,2)=3;
	vf1(3,0)=0, vf1(3,1)=2, vf1(3,2)=1;
		
	for(j=0;j<npnt;++j) {
		p0 = pnt(j).gindx;
		
		/* THIS IS THE RIGHT WAY, BUT FINDBDRYPT ISN'T WORKING YET */
		tgt.fbdry(bnum)->findbdrypt(x.pnts(p0), facloc, r, s);
		tind = tgt.fbdry(bnum)->tri(facloc).gindx;
		ttind = tgt.tri(tind).tet(0);                             
		cnnct(p0).tet = ttind;
		for (k=0;k<4;++k) 
			if (tgt.tet(ttind).tri(k) == tind) break;
		assert(k < 4);
		cnnct(p0).wt(k) = 0.0;

		// Assumes Definition of face is compatible with definition on tet, which should always be true :-)
		cnnct(p0).wt(vf1(k,0)) = 0.5*(1.+s);
		cnnct(p0).wt(vf1(k,2)) = 0.5*(-r -s);
		cnnct(p0).wt(vf1(k,1)) = 0.5*(1.+r);
	}
	
	return;
}

void face_bdry::mvpttobdry(int indx, FLT r, FLT s, TinyVector<FLT,tet_mesh::ND> &pt) {
	/* FOR A LINEAR SIDE */
	int n;
	TinyVector<int,3> v = tri(indx).pnt;
	for (int i=0; i<3;++i) 
		v(i) = pnt(v(i)).gindx;
		
	for (n=0;n<tet_mesh::ND;++n)
		pt(n) = 0.5*((1.+s)*x.pnts(v(0))(n) +(-r-s)*x.pnts(v(1))(n) +(1.+r)*x.pnts(v(2))(n));
	
	return;
}

void face_bdry::mvpttobdry(int indx, FLT psi, TinyVector<FLT,tet_mesh::ND> &pt) {
	/* FOR A LINEAR SIDE */
	int n;
	
	for (n=0;n<tet_mesh::ND;++n)
		pt(n) = 0.5*((1. -psi)*x.pnts(x.seg(seg(indx).gindx).pnt(0))(n) +(1.+psi)*x.pnts(x.seg(seg(indx).gindx).pnt(1))(n));
	
	return;
}

void fcomm::ploadbuff(boundary::groups grp,FLT *base,int bgn,int end, int stride) {
	int j,k,count,offset;
	
	if (!in_group(grp)) return;

	count = 0;
	for(j=0;j<npnt;++j) {
		offset = pnt(j).gindx*stride;
		for (k=bgn;k<=end;++k) {
			fsndbuf(count++) = base[offset+k];
		}
	}
		
	sndsize() = count;
	sndtype() = boundary::flt_msg;
}

void fcomm::pfinalrcv(boundary::groups grp, int phi, comm_type type, operation op, FLT *base,int bgn,int end, int stride) {
	int j,k,count,offset;
	/* ASSUMES REVERSE ORDERING OF SIDES */
	/* WON'T WORK IN 3D */
		
	bool reload = comm_finish(grp,phi,type,op);
	if (!reload) return;
	
//	*x.gbl->log << "I am here " << idprefix << ' ' << grp << ' ' << phi << ' ' << bgn << ' ' << end << ' ' << stride << std::endl;

	count = 0;
	for(j=0;j<npnt;++j) {
		offset = pnt(j).gindx*stride;
		for (k=bgn;k<=end;++k) {
			base[offset+k] = fsndbuf(count++);
		}
	}  
	return;
}

void fcomm::sloadbuff(boundary::groups grp,FLT *base,int bgn,int end, int stride) {
	int j,k,count,sind,offset;
	
	if (!in_group(grp)) return;

	count = 0;
	for(j=0;j<nseg;++j) {
		sind = seg(j).gindx;
		offset = sind*stride;
		for (k=bgn;k<=end;++k) {
			fsndbuf(count++) = base[offset+k];
		}
	}
		
	sndsize() = count;
	sndtype() = boundary::flt_msg;
}

void fcomm::sfinalrcv(boundary::groups grp, int phi, comm_type type, operation op, FLT *base,int bgn,int end, int stride) {
	int j,k,count,offset,sind;
	
	bool reload = comm_finish(grp,phi,type,op);
	if (!reload) return;
	
	count = 0;
	for(j=0;j<nseg;++j) {
		sind = seg(j).gindx;
		offset = sind*stride;
		for (k=bgn;k<=end;++k) {
			base[offset+k] = fsndbuf(count++);
		}

	}
	return;
}

void fcomm::tloadbuff(boundary::groups grp,FLT *base,int bgn,int end, int stride) {
	int j,k,count,tind,offset;
	
	if (!in_group(grp)) return;

	count = 0;
	for(j=0;j<ntri;++j) {
		tind = tri(j).gindx;
		offset = tind*stride;
		for (k=bgn;k<=end;++k) {
			fsndbuf(count++) = base[offset+k];
		}
	}
		
	sndsize() = count;
	sndtype() = boundary::flt_msg;
}

void fcomm::tfinalrcv(boundary::groups grp, int phi, comm_type type, operation op, FLT *base,int bgn,int end, int stride) {
	int j,k,count,offset,tind;
	/* ASSUMES REVERSE ORDERING OF SIDES */
	/* WON'T WORK IN 3D */
	
	bool reload = comm_finish(grp,phi,type,op);
	if (!reload) return;

	count = 0;
	for(j=0;j<ntri;++j) {
		tind = tri(j).gindx;
		offset = tind*stride;
		for (k=bgn;k<=end;++k) {
			base[offset+k] = fsndbuf(count++);
		}
	}
	return;
}

void fcomm::match_numbering(int step) {
	switch (step) {
		case(1): {
			if (is_frst()) {
				return;
			}
			else {
				/* Slave receives master boundaries vertex positions */
				int count = 0;
				bool printstuff = false;
				for (int i=0;i<npnt;++i) {
					TinyVector<FLT,tet_mesh::ND> mpnt;
					for (int n=0;n<tet_mesh::ND;++n)
						mpnt(n) = fsndbuf(count++);

					/* brute force */
//					FLT dist;
//					pnt(i).gindx = 0;
//					for (int j=0; j<x.npnt; ++j) {
//						TinyVector<FLT,tet_mesh::ND> lclpnts;
//						lclpnts = x.pnts(j);
//						dist = 0.0;
//						for (int k=0; k<tet_mesh::ND; ++k) {
//							dist += pow(lclpnts(k)-mpnt(k), 2.0);
//						}
//
//						if (dist < 1.0e-8) {
//							pnt(i).gindx = j;
//							break;
//						}
//					}
					
					FLT dist = x.otree.nearpt(mpnt.data(),pnt(i).gindx);

					if (dist > 10.*EPSILON) {
						printstuff = true;
						*x.gbl->log << "Matching face numbering error: " << dist << ' ' << mpnt << ' ' << x.pnts(pnt(i).gindx) << '\n';
						*x.gbl->log << "idnum " << idnum << " lcl point " << i << " gbl point " << pnt(i).gindx << '\n';
					}

				}
				if (printstuff) {
					otree.output("octree_error");
					exit(3);
				}
			}
			return;
		}
		
		case(2): {
			if (is_frst()) {
				/* Master loads integer data */
				int count = 0;
				for (int i=0;i<nseg;++i) {
					isndbuf(count++) = seg(i).pnt(0);
					isndbuf(count++) = seg(i).pnt(1);
				}
				
				for (int i=0;i<ntri;++i) {
					for(int j=0;j<3;++j) {
						isndbuf(count++) = tri(i).pnt(j);
					}
				}
				sndsize() = count;
				sndtype() = boundary::int_msg;
			}
			else {
				sndsize() = 2*nseg +3*ntri;
				sndtype() = boundary::int_msg;
			}
			return;
		}
		
		case(3): {
			if (is_frst()) return;
			
			int count = 0;
			for (int i=0;i<nseg;++i) {
				seg(i).pnt(0) = isndbuf(count++);
				seg(i).pnt(1) = isndbuf(count++);
			}
			create_seg_gindx();
			
			/* Triangles are reverse oriented */
			for (int i=0;i<ntri;++i) {
				tri(i).pnt(0) = isndbuf(count++);
				tri(i).pnt(2) = isndbuf(count++);
				tri(i).pnt(1) = isndbuf(count++);
			}
			create_tri_gindx();
		}
	}
	
	return;
}


void face_bdry::pull_apart_face_boundaries() {
	Array<int,2> listoftris(ntri,ntri);
	Array<int,1> numoftris(ntri);
	
	listoftris = -1;
	numoftris = 0;	
	
	for (int i=0; i < ntri; ++i) 
		tri(i).info = -1;
	
	int chunk = 0;
	
	for (int n=0; n<ntri; ++n) {
		if (tri(n).info < 0) {	
			
			tri(n).info = chunk;
			numoftris(chunk) = 1;
			listoftris(0,chunk) = n;
			
			for (int i=0; i < numoftris(chunk); ++i) {
				for (int j=0; j<3; ++j) {
					
					int tind = tri(listoftris(i,chunk)).tri(j);
					
					if (tri(tind).info < 0 && tind > -1) {
						tri(tind).info = chunk;
						listoftris(numoftris(chunk)++,chunk) = tind;
					}
				}
			}
			
			++chunk;
		}
	}
	
	if (chunk>1) {
		*x.gbl->log << "#creating " << chunk-1 << " new " << mytype << " face boundaries: " << idnum  << std::endl;
		
		Array<int,2> listofgindx(ntri,chunk);
		listofgindx = -1;
		
		for (int i=0; i<chunk; ++i) {
			for (int j=0; j<numoftris(i); ++j) {
				listofgindx(j,i) = tri(listoftris(j,i)).gindx;
			}
		}
		
		x.nfbd += chunk-1;
		x.fbdry.resizeAndPreserve(x.nfbd);
		ntri = numoftris(0);
		alloc(static_cast<int>(ntri*3));
		
		for (int i=0; i<ntri; ++i) {
			tri(i).gindx = listofgindx(i,0);
			for (int j=0; j<3; ++j) {
				tri(i).pnt(j) = x.tri(tri(i).gindx).pnt(j);
			}
		}
		
		create_from_tri();
		
		for (int n=1; n<chunk; ++n) {
			int fbd = x.nfbd-n;
			x.fbdry(fbd) = create(x);
			//x.fbdry(fbd)->copy(*this);
			x.fbdry(fbd)->alloc(static_cast<int>(numoftris(n)*3));
			x.fbdry(fbd)->ntri = numoftris(n);
			
			for (int i=0; i<numoftris(n); ++i) {
				x.fbdry(fbd)->tri(i).gindx = listofgindx(i,n);
				for (int j=0; j<3; ++j) {
					x.fbdry(fbd)->tri(i).pnt(j) = x.tri(x.fbdry(fbd)->tri(i).gindx).pnt(j);
				}
			}
			
			x.fbdry(fbd)->create_from_tri();
			
		}	
	}
	
	return;
}

void fpartition::mgconnect(Array<tet_mesh::transfer,1> &cnnct,tet_mesh& tgt, int bnum) {
	cout << "fpartition::mgconnect being called and doesnt work" << endl;
	
	//int i,j,k,p0;

	//	/* BOUNDARY IS AN INTERNAL PARTITION BOUNDARY */
	//	/* MAKE SURE ENDPOINTS ARE OK */
	//	i = x.seg(seg(0).gindx).pnt(0);
	//	if (cnnct(i).tet < 0) {
	//		tgt.otree.nearpt(x.pnts(i).data(),p0);
	//		cnnct(i).tet=tgt.pnt(p0).tet;
	//		for(j=0;j<4;++j) {
	//			cnnct(i).wt(j) = 0.0;
	//			if (tgt.tet(cnnct(i).tet).pnt(j) == p0) cnnct(i).wt(j) = 1.0;
	//		}
	//	}
	//	i = x.seg(seg(nseg-1).gindx).pnt(1);
	//	if (cnnct(i).tet < 0) {
	//		tgt.otree.nearpt(x.pnts(i).data(),p0);
	//		cnnct(i).tet=tgt.pnt(p0).tet;
	//		for(j=0;j<3;++j) {
	//			cnnct(i).wt(j) = 0.0;
	//			if (tgt.tet(cnnct(i).tet).pnt(j) == p0) cnnct(i).wt(j) = 1.0;
	//		}
	//	}
	//	
	//	if (first) {
	//		sndsize() = 0;
	//		sndtype() = int_msg;
	//		for(k=1;k<nseg;++k) {
	//			p0 = x.seg(seg(k).gindx).pnt(0);
	//			if (cnnct(p0).tet > 0) {
	//				isndbuf(sndsize()++) = -1;
	//			}
	//			else {
	//				isndbuf(sndsize()++) = +1;
	//				cnnct(p0).tet = 0;
	//				for(j=0;j<3;++j)
	//					cnnct(p0).wt(j) = 0.0;
	//			}
	//		}
	//	}
	//	
	//	comm_prepare(boundary::all,0,slave_master); 
	//	comm_exchange(boundary::all,0,slave_master);
	//	comm_wait(boundary::all,0,slave_master);
	//	
	//	if (!first) {
	//		i = 0;
	//		for(k=nseg-1;k>0;--k) {
	//			p0 = x.seg(seg(k).gindx).pnt(1);
	//			if (ircvbuf(0,i) < 0) {
	//				cnnct(p0).tet = 0;
	//				for(j=0;j<3;++j)
	//					cnnct(p0).wt(j) = 0.0;
	//			}
	//		}
	//	}                    
}

			




