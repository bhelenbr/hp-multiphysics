#include "tri_boundary.h"
#include <assert.h>
#include <float.h>

/********************/
/* VERTEX FUNCTIONS */
/********************/

/* GENERIC VERTEX COMMUNICATIONS */
void vcomm::vloadbuff(boundary::groups grp,FLT *base,int bgn,int end, int stride) {
	int i,offset;

	if (!in_group(grp)) return;

	sndsize()=end-bgn+1;
	sndtype()=flt_msg;

	/* LOAD SEND BUFFER */
	offset = pnt*stride +bgn;
	for (i=0;i<end-bgn+1;++i)
		fsndbuf(i) = base[offset+i];
}

void vcomm::vfinalrcv(boundary::groups grp, int phi, comm_type type, operation op, FLT *base, int bgn, int end, int stride) {
	int i,offset;
	bool reload = comm_finish(grp,phi,type,op);
	if (!reload) return;

	offset = pnt*stride +bgn;
	for (i=0;i<end-bgn+1;++i)
		base[offset+i] = fsndbuf(i);
}

/**************************************/
/* GENERIC FUNCTIONS FOR SIDES          */
/**************************************/
void edge_bdry::alloc(int n) {
	maxseg = n;
	seg.resize(Range(-1,n));
	prev.resize(Range(-1,n));
	next.resize(Range(-1,n));
	prev = 0;
	next = 0;
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

void edge_bdry::mvpttobdry(int indx, FLT psi, TinyVector<FLT,tri_mesh::ND> &pt) {
	/* FOR A LINEAR SIDE */
	int n;

	for (n=0;n<tri_mesh::ND;++n)
		pt(n) = 0.5*((1. -psi)*x.pnts(x.seg(seg(indx)).pnt(0))(n) +(1.+psi)*x.pnts(x.seg(seg(indx)).pnt(1))(n));

	return;
}

void edge_bdry::bdry_normal(int indx, FLT psi, TinyVector<FLT,tri_mesh::ND> &norm) {
	/* FOR A LINEAR SIDE */	
	for (int n=0;n<tri_mesh::ND;++n)
		norm(1-n) = x.pnts(x.seg(seg(indx)).pnt(1))(n) -x.pnts(x.seg(seg(indx)).pnt(0))(n);
	norm(1) *= -1.0;
	norm /= sqrt(norm(0)*norm(0) +norm(1)*norm(1));
	
	return;
}



void edge_bdry::findbdrypt(const TinyVector<FLT,tri_mesh::ND> xpt, int &sidloc, FLT &psiloc) const {
	int k,sind,p0,p1,sidlocprev;
	FLT dx,dy,ol,psi,normdist;
	FLT psiprev,normdistprev;
	FLT mindist = 1.0e32;

	if (x.seg(seg(0)).pnt(0) == x.seg(seg(nseg-1)).pnt(1)) {
		/* BOUNDARY IS A LOOP */
		sind = seg(nseg-1);
		p0 = x.seg(sind).pnt(0);
		p1 = x.seg(sind).pnt(1);
		dx = x.pnts(p1)(0) - x.pnts(p0)(0);
		dy = x.pnts(p1)(1) - x.pnts(p0)(1);
		ol = 2./(dx*dx +dy*dy);
		psi = ol*((xpt(0) -x.pnts(p0)(0))*dx +(xpt(1) -x.pnts(p0)(1))*dy) -1.;
		normdist = dx*(xpt(1)-x.pnts(p0)(1))-dy*(xpt(0)-x.pnts(p1)(0));
		normdist *= sqrt(ol/2.);
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
		sind = seg(k);
		p0 = x.seg(sind).pnt(0);
		p1 = x.seg(sind).pnt(1);
		dx = x.pnts(p1)(0) - x.pnts(p0)(0);
		dy = x.pnts(p1)(1) - x.pnts(p0)(1);
		ol = 2./(dx*dx +dy*dy);
		psi = ol*((xpt(0) -x.pnts(p0)(0))*dx +(xpt(1) -x.pnts(p0)(1))*dy) -1.;
		normdist = dx*(xpt(1)-x.pnts(p0)(1))-dy*(xpt(0)-x.pnts(p1)(0));
		normdist *= sqrt(ol/2.);

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



void edge_bdry::mgconnect(Array<tri_mesh::transfer,1> &cnnct, tri_mesh& tgt, int bnum) {
	int j,k,sind,tind,p0,sidloc;
	FLT psiloc;


	/* END POINTS SHOULD ALWAYS BE COINCIDENT */
	p0 = x.seg(seg(0)).pnt(0);
	sind = tgt.ebdry(bnum)->seg(0);
	tind = tgt.seg(sind).tri(0);
	cnnct(p0).tri = tind;
	for(j=0;j<3;++j) {
		cnnct(p0).wt(j) = 0.0;
		if (tgt.tri(tind).pnt(j) == tgt.seg(sind).pnt(0)) cnnct(p0).wt(j) = 1.0;
	}
	
	p0 = x.seg(seg(nseg-1)).pnt(1);
	sind = tgt.ebdry(bnum)->seg(tgt.ebdry(bnum)->nseg-1);
	tind = tgt.seg(sind).tri(0);
	cnnct(p0).tri = tind;
	for(j=0;j<3;++j) {
		cnnct(p0).wt(j) = 0.0;
		if (tgt.tri(tind).pnt(j) == tgt.seg(sind).pnt(1)) cnnct(p0).wt(j) = 1.0;
	}

	for(k=1;k<nseg;++k) {
		p0 = x.seg(seg(k)).pnt(0);
		tgt.ebdry(bnum)->findbdrypt(x.pnts(p0), sidloc, psiloc);
		sind = tgt.ebdry(bnum)->seg(sidloc);
		tind = tgt.seg(sind).tri(0);
		cnnct(p0).tri = tind;
		for (j=0;j<3;++j)
			if (tgt.tri(tind).seg(j) == sind) break;
		assert(j < 3);
		cnnct(p0).wt(j) = 0.0;
		cnnct(p0).wt((j+1)%3) = 0.5*(1.-psiloc);
		cnnct(p0).wt((j+2)%3) = 0.5*(1.+psiloc);
	}

	return;
}

/* SWAP ELEMENTS IN LIST */
void edge_bdry::swap(int s1, int s2) {
	int ind;

	ind = seg(s1);
	seg(s1) = seg(s2);
	seg(s2) = ind;

	/* PREV/NEXT IMPLEMENTATION */
	ind = prev(s1);
	prev(s1) = prev(s2);
	prev(s2) = ind;

	ind = next(s1);
	next(s1) = next(s2);
	next(s2) = ind;

	int prev1 = prev(s1);
	int next1 = next(s1);
	int prev2 = prev(s2);
	int next2 = next(s2);

	if (prev1 == s1) prev1 = s2;
	if (next1 == s1) next1 = s2;
	if (prev2 == s2) prev2 = s1;
	if (next2 == s2) next2 = s1;

	next(prev1) = s1;
	prev(next1) = s1;

	next(prev2) = s2;
	prev(next2) = s2;

	return;
}


/* REORDERS BOUNDARIES TO BE SEQUENTIAL */
/* USES gbl->intwk & gbl->i2wk AS WORK ARRAYS */
void edge_bdry::reorder() {
	int i,count,total,sind,minp,first;

	total = nseg;

	/* DON'T ASSUME wk INITIALIZED TO -1 */
	for(i=0;i<nseg;++i) {
		sind = seg(i);
		x.gbl->intwk(x.seg(sind).pnt(0)) = -1;
		x.gbl->i2wk(x.seg(sind).pnt(1)) = -1;
	}

	/* STORE SIDE INDICES BY VERTEX NUMBER */
	for(i=0; i < nseg; ++i) {
		sind = seg(i);
		x.gbl->intwk(x.seg(sind).pnt(1)) = i;
		x.gbl->i2wk(x.seg(sind).pnt(0)) = i;
	}
	
	/* Implement a break at vertex boundaries */
	/* For weird case of two boundaries with same number that must be split */
	for(i=0;i<x.nvbd;++i)
		x.gbl->i2wk(x.vbdry(i)->pnt) = -1;

	/* FIND FIRST SIDE */
	first = -1;
	for(i=0;i<nseg;++i) {
		sind = seg(i);
		if (x.gbl->intwk(x.seg(sind).pnt(0)) == -1) {
			first = i;
			break;
		}
	}

	/* SPECIAL CONSTRAINT IF LOOP */
	/* THIS IS TO ELIMINATE ANY INDEFINITENESS ABOUT SIDE ORDERING FOR LOOP */
	if (first < 0) {
		minp = x.npnt;
		for(i=0;i<nseg;++i) {
			sind = seg(i);
			if (x.seg(sind).pnt(0) < minp) {
				first = i;
				minp = x.seg(sind).pnt(0);
			}
		}
	}

	/* SWAP FIRST SIDE */
	count = 0;
	swap(count,first);
	x.gbl->intwk(x.seg(seg(first)).pnt(1)) = first;
	x.gbl->i2wk(x.seg(seg(first)).pnt(0)) = first;
	x.gbl->intwk(x.seg(seg(count)).pnt(1)) = count;
	x.gbl->i2wk(x.seg(seg(count)).pnt(0)) = -1;  // TO MAKE SURE LOOP STOPS

	/* REORDER LIST */
	while ((first = x.gbl->i2wk(x.seg(seg(count++)).pnt(1))) >= 0) {
		swap(count,first);
		x.gbl->intwk(x.seg(seg(first)).pnt(1)) = first;
		x.gbl->i2wk(x.seg(seg(first)).pnt(0)) = first;
		x.gbl->intwk(x.seg(seg(count)).pnt(1)) = count;
		x.gbl->i2wk(x.seg(seg(count)).pnt(0)) = count;
	}

	/* RESET gbl->intwk TO -1 */
	for(i=0; i <total; ++i) {
		sind = seg(i);
		x.gbl->intwk(x.seg(sind).pnt(1)) = -1;
	}

	if (count < total) {
		++x.nebd;
		x.ebdry.resizeAndPreserve(x.nebd);
		x.ebdry(x.nebd-1) = create(x);
		x.ebdry(x.nebd-1)->copy(*this);
		nseg = count;

		for(i=0;i<total-nseg;++i)
			x.ebdry(x.nebd-1)->swap(i,i+nseg);
		x.ebdry(x.nebd-1)->nseg = total-nseg;
		*x.gbl->log << "#creating new boundary: " << idnum << " num: " << x.ebdry(x.nebd-1)->nseg << std::endl;
	}

	for (i=0;i<nseg;++i) {
		prev(i) = i-1;
		next(i) = i+1;
	}
	next(nseg-1) = -1;

	return;
}

void ecomm::vloadbuff(boundary::groups grp,FLT *base,int bgn,int end, int stride) {
	int j,k,count,offset;
	int sind = 0;  // To avoid may be used uninitialized warnings

	if (!in_group(grp)) return;

	if (first) {
		count = 0;
		for(j=0;j<nseg;++j) {
			sind = seg(j);
			offset = x.seg(sind).pnt(0)*stride;
			for (k=bgn;k<=end;++k) {
				fsndbuf(count++) = base[offset+k];
			}
		}
		offset = x.seg(sind).pnt(1)*stride;
		for (k=bgn;k<=end;++k)
			fsndbuf(count++) = base[offset+k];
	}
	else {
		count = 0;
		for(j=nseg-1;j>=0;--j) {
			sind = seg(j);
			offset = x.seg(sind).pnt(1)*stride;
			for (k=bgn;k<=end;++k) {
				fsndbuf(count++) = base[offset+k];
			}
		}
		offset = x.seg(sind).pnt(0)*stride;
		for (k=bgn;k<=end;++k)
			fsndbuf(count++) = base[offset+k];
	}

	sndsize() = count;
	sndtype() = boundary::flt_msg;
}

void ecomm::vfinalrcv(boundary::groups grp, int phi, comm_type type, operation op, FLT *base,int bgn,int end, int stride) {
	int j,k,count,offset;
	int sind = 0;  // To avoid may be used uninitialized warnings


	bool reload = comm_finish(grp,phi,type,op);
	if (!reload) return;

	int ebp1 = end-bgn+1;
	if (first) {
		count = 0;
		for(j=0;j<nseg;++j) {
			sind = seg(j);
			offset = x.seg(sind).pnt(0)*stride +bgn;
			for(k=0;k<ebp1;++k) {
				base[offset+k] = fsndbuf(count++);
			}
		}
		offset = x.seg(sind).pnt(1)*stride +bgn;
		for(k=0;k<ebp1;++k) {
			base[offset+k] = fsndbuf(count++);
		}
	}
	else {
		count = 0;
		for(j=nseg-1;j>=0;--j) {
			sind = seg(j);
			offset = x.seg(sind).pnt(1)*stride;
			for (k=bgn;k<=end;++k) {
				base[offset+k] = fsndbuf(count++);
			}
		}
		offset = x.seg(sind).pnt(0)*stride;
		for (k=bgn;k<=end;++k)
			base[offset+k] = fsndbuf(count++);
	}

	return;
}

void ecomm::sloadbuff(boundary::groups grp,FLT *base,int bgn,int end, int stride) {
	int j,k,count,sind,offset;

	if (!in_group(grp)) return;

	if (first) {
		count = 0;
		for(j=0;j<nseg;++j) {
			sind = seg(j);
			offset = sind*stride;
			for (k=bgn;k<=end;++k) {
				fsndbuf(count++) = base[offset+k];
			}
		}
	}
	else {
		count = 0;
		for(j=nseg-1;j>=0;--j) {
			sind = seg(j);
			offset = sind*stride;
			for (k=bgn;k<=end;++k) {
				fsndbuf(count++) = base[offset+k];
			}
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

	if (first) {
		count = 0;
		for(j=0;j<nseg;++j) {
			sind = seg(j);
			offset = sind*stride;
			for (k=bgn;k<=end;++k) {
				base[offset+k] = fsndbuf(count++);
			}
		}
	}
	else {
		count = 0;
		for(j=nseg-1;j>=0;++j) {
			sind = seg(j);
			offset = sind*stride;
			for (k=bgn;k<=end;++k) {
				base[offset+k] = fsndbuf(count++);
			}
		}
	}
}

void epartition::alloc(int size) {
	ecomm::alloc(size);
	resize_buffers(4*size*(2+2+3)); // Four rows of pnts, seg pnts, tri pnts
	remote_halo.gbl = x.gbl;
	remote_halo.allocate(8*size);
	/* Make 2 Boundaries first is the matching boundary */
	remote_halo.nebd = 2;
	remote_halo.nvbd = 0;
	remote_halo.ebdry.resize(remote_halo.nebd);
	remote_halo.ebdry(0) = new ecomm(idnum,x);
	remote_halo.ebdry(0)->alloc(size);
	remote_halo.ebdry(1) = new edge_bdry(222,x);
	remote_halo.ebdry(1)->alloc(size);
	
	pnt_h.resize(2*size);
	seg_h.resize(2*size);
	sgn_h.resize(2*size);
	tri_h.resize(2*size);
}

void epartition::copy(const edge_bdry& bin) {
	ecomm::copy(bin);
	const epartition& tgt(dynamic_cast<const epartition&>(bin));
	ntri_h = tgt.ntri_h;
	nseg_h = tgt.nseg_h;
	npnt_h = tgt.npnt_h;
	if (ntri_h) {
		tri_h(Range(0,ntri_h-1)) = tgt.tri_h(Range(0,ntri_h-1)); 
		seg_h(Range(0,nseg_h-1)) = tgt.seg_h(Range(0,nseg_h-1));
		sgn_h(Range(0,nseg_h-1)) = tgt.sgn_h(Range(0,nseg_h-1));
		pnt_h(Range(0,npnt_h-1)) = tgt.pnt_h(Range(0,npnt_h-1));
		remote_halo.copy(tgt.remote_halo);
	}
	
	return;
}

void epartition::mgconnect(Array<tri_mesh::transfer,1> &cnnct, tri_mesh& tgt, int bnum) {
	int j,k,p0,sind,tind;

	if (tgt.npnt > x.npnt) {
		/* coarse to fine connection will be done differently because all points will be coincident */
		ecomm::mgconnect(cnnct,tgt,bnum);
		return;
	}
	
	/* END POINTS SHOULD ALWAYS BE COINCIDENT */
	p0 = x.seg(seg(0)).pnt(0);
	sind = tgt.ebdry(bnum)->seg(0);
	tind = tgt.seg(sind).tri(0);
	cnnct(p0).tri = tind;
	for(j=0;j<3;++j) {
		cnnct(p0).wt(j) = 0.0;
		if (tgt.tri(tind).pnt(j) == tgt.seg(sind).pnt(0)) cnnct(p0).wt(j) = 1.0;
	}
	
	p0 = x.seg(seg(nseg-1)).pnt(1);
	sind = tgt.ebdry(bnum)->seg(tgt.ebdry(bnum)->nseg-1);
	tind = tgt.seg(sind).tri(0);
	cnnct(p0).tri = tind;
	for(j=0;j<3;++j) {
		cnnct(p0).wt(j) = 0.0;
		if (tgt.tri(tind).pnt(j) == tgt.seg(sind).pnt(1)) cnnct(p0).wt(j) = 1.0;
	}

	if (first) {
		sndsize() = 0;
		sndtype() = int_msg;
		/* First and last point is always found */
		/* Only need to do interior points */
		for(k=1;k<nseg;++k) {
			p0 = x.seg(seg(k)).pnt(0);
			if (cnnct(p0).tri >= 0) {
				isndbuf(sndsize()++) = -1;
			}
			else {
				isndbuf(sndsize()++) = +1;
				cnnct(p0).tri = 0;
				for(j=0;j<3;++j)
					cnnct(p0).wt(j) = 0.0;
			}
		}
	}
	comm_prepare(boundary::partitions,0,master_slave);
}


void epartition::mgconnect1(Array<tri_mesh::transfer,1> &cnnct, tri_mesh& tgt, int bnum) {

	if (tgt.npnt > x.npnt) {
		/* coarse to fine connection will be done differently because all points will be coincident */
		ecomm::mgconnect1(cnnct,tgt,bnum);
		return;
	}
	
	comm_wait(boundary::partitions,0,master_slave);
	
	if (!first) {
		/* Check if endpoints are communication boundaries */
		if (!x.vbdry(vbdry(0))->is_comm()) {
			int p0 = x.seg(seg(0)).pnt(0);
			cnnct(p0).tri = 0;
			for(int j=0;j<3;++j)
				cnnct(p0).wt(j) = 0.0;
		}
		
		if (!x.vbdry(vbdry(1))->is_comm()) {
			int p0 = x.seg(seg(nseg-1)).pnt(1);
			cnnct(p0).tri = 0;
			for(int j=0;j<3;++j)
				cnnct(p0).wt(j) = 0.0;
		}
		
	
		int i = 0;
		for(int k=nseg-1;k>0;--k) {
			int p0 = x.seg(seg(k)).pnt(0);
			if (ircvbuf(0,i++) < 0) {
				cnnct(p0).tri = 0;
				for(int j=0;j<3;++j)
					cnnct(p0).wt(j) = 0.0;
			}
		}
	}
}

const int TSRCH = 0x100*0x4;
const int PSPEC = 0x4;
const int SSRCH = 0x10*0x4;
#if ((-1)&(0x100*0x4))
#define ISTSRCH(A) (!((A)&TSRCH))
#define SETTSRCH(A) A&=(~TSRCH)
#define CLRTSRCH(A) A|=(TSRCH)
#define ISSSRCH(A) (!((A)&SSRCH))
#define SETSSRCH(A) A&=(~SSRCH)
#define CLRSSRCH(A) A|=(SSRCH)
#define ISSPEC(A) (!((A)&PSPEC))
#define SETSPEC(A) A&=(~PSPEC)
#define CLRSPEC(A) A|=(PSPEC)
#else
#define ISTSRCH(A) (((A)&TSRCH))
#define SETTSRCH(A) A|=(TSRCH)
#define CLRTSRCH(A) A&=(~TSRCH)
#define ISSSRCH(A) (((A)&SSRCH))
#define SETSSRCH(A) A|=(SSRCH)
#define CLRSSRCH(A) A&=(~SSRCH)
#define ISSPEC(A) (((A)&PSPEC))
#define SETSPEC(A) A|=(PSPEC)
#define CLRSPEC(A) A&=(~PSPEC)
#endif


void epartition::calculate_halo() {
	int tind,vn;
	
	nseg_h = 0;
	/* MARK BOUNDARY SO DON'T GET INSERTED */
	/* FUNNY WAY OF MARKING SO CAN LEAVE gbl->intwk initialized to -1 */
	for(int j=0;j<nseg;++j) {
		int sind = seg(j);
		seg_h(nseg_h++) = sind;
		SETSSRCH(x.gbl->intwk(sind));
	}
	
	/* Find first side by going counter-clockwise */
	int sind = seg(0);
	int ppivot = x.seg(sind).pnt(0);
	int tindnext = x.seg(sind).tri(0);
	do {
		tind = tindnext;
		for(vn=0;vn<3;++vn)
			if (x.tri(tind).pnt(vn) == ppivot) break;
		
		tindnext = x.tri(tind).tri((vn +1)%3);
	} while(tindnext > -1);
				
	/* Found first side */
	seg_h(nseg_h++) = x.tri(tind).seg((vn +1)%3);
	SETSSRCH(x.gbl->intwk(seg_h(nseg_h-1)));
	
	ntri_h = 0;
	for(int scnt=0;scnt<nseg+1;++scnt) {
		tindnext = tind;
		/* Now go clockwise for the rest */
		do {
			tind = tindnext;
			for(vn=0;vn<3;++vn)
				if (x.tri(tind).pnt(vn) == ppivot) break;
			assert(vn < 3);
			
			if (!ISTSRCH(x.gbl->intwk(tind))) {
				tri_h(ntri_h++) = tind;
				SETTSRCH(x.gbl->intwk(tind));
				
				sind = x.tri(tind).seg(vn);
				if (!ISSSRCH(x.gbl->intwk(sind))) {
					seg_h(nseg_h++) = sind;
					SETSSRCH(x.gbl->intwk(sind));
				}

				sind = x.tri(tind).seg((vn+2)%3);
				if (!ISSSRCH(x.gbl->intwk(sind))) {
					seg_h(nseg_h++) = sind;
					SETSSRCH(x.gbl->intwk(sind));
				}
			}

			tindnext = x.tri(tind).tri((vn+2)%3);
		} while (tindnext > -1);
		
		/* Next ppivot */
		ppivot = x.tri(tind).pnt((vn+1)%3);
	}
	
	for(int i=0;i<ntri_h;++i)
		CLRTSRCH(x.gbl->intwk(tri_h(i)));
	
	for (int j=0;j<nseg_h;++j) {
		sind = seg_h(j);
		CLRSSRCH(x.gbl->intwk(sind));
	}
	
	/* Extract Halo Mesh for sending to remote partition */
	int pnt;
	/* load points */
	remote_halo.npnt = 0;
	for (int j=0;j<nseg_h;++j) {
		int pnt = x.seg(seg_h(j)).pnt(0);
		if (x.gbl->intwk(pnt) == -1) {
			x.gbl->intwk(pnt) = remote_halo.npnt;
			pnt_h(remote_halo.npnt) = pnt;
			remote_halo.pnts(remote_halo.npnt++) = x.pnts(pnt);
		}
		
		pnt = x.seg(seg_h(j)).pnt(1);
		if (x.gbl->intwk(pnt) == -1) {
			x.gbl->intwk(pnt) = remote_halo.npnt;
			pnt_h(remote_halo.npnt) = pnt;
			remote_halo.pnts(remote_halo.npnt++) = x.pnts(pnt);
		}
	}
	npnt_h = remote_halo.npnt;
	
	/* Create Side List */
	remote_halo.nseg = 0;
	for (int j=0;j<nseg_h;++j) {
		remote_halo.seg(remote_halo.nseg).pnt(0) = x.gbl->intwk(x.seg(seg_h(j)).pnt(0));
		remote_halo.seg(remote_halo.nseg).pnt(1) = x.gbl->intwk(x.seg(seg_h(j)).pnt(1));
		++remote_halo.nseg;
	}
	
	/* Create Tri List */
	remote_halo.ntri = ntri_h;
	for (int tind=0;tind<ntri_h;++tind) {
		for(int vn=0;vn<3;++vn) {
			remote_halo.tri(tind).pnt(vn) = x.gbl->intwk(x.tri(tri_h(tind)).pnt(vn));
		}
	}

	for(int j=0;j<npnt_h;++j)
		x.gbl->intwk(pnt_h(j)) = -1;
	
	/* Fill in connections between seg and tri so can figure out what segs are on boundary */
	remote_halo.createsegtri();
	
	for(int i=0;i<remote_halo.nseg;++i) {
		if (remote_halo.seg(i).tri(0) < 0) {
			/* Side is on boundary and direction must be reversed */
			int temp = remote_halo.seg(i).pnt(0);
			remote_halo.seg(i).pnt(0) = remote_halo.seg(i).pnt(1);
			remote_halo.seg(i).pnt(1) = temp;
			sgn_h(i) = -1;
		}
		else {
			sgn_h(i) = 1;
		}
	}
	
	/* pack up to send messages */
	sndtype() = boundary::flt_msg;
	sndsize() = 0;
	fsndbuf(sndsize()++) = remote_halo.npnt+100*FLT_EPSILON;
	fsndbuf(sndsize()++) = remote_halo.nseg+100*FLT_EPSILON;
	fsndbuf(sndsize()++) = remote_halo.ntri+100*FLT_EPSILON;
	for(int j=0;j<remote_halo.npnt;++j) {
		fsndbuf(sndsize()++) = remote_halo.pnts(j)(0);
		fsndbuf(sndsize()++) = remote_halo.pnts(j)(1);
	}
	for(int j=0;j<remote_halo.nseg;++j) {
		fsndbuf(sndsize()++) = remote_halo.seg(j).pnt(0)+100*FLT_EPSILON;
		fsndbuf(sndsize()++) = remote_halo.seg(j).pnt(1)+100*FLT_EPSILON;
	}
	
	for(int j=0;j<remote_halo.ntri;++j) {
		fsndbuf(sndsize()++) = remote_halo.tri(j).pnt(0)+100*FLT_EPSILON;
		fsndbuf(sndsize()++) = remote_halo.tri(j).pnt(1)+100*FLT_EPSILON;
		fsndbuf(sndsize()++) = remote_halo.tri(j).pnt(2)+100*FLT_EPSILON;
	}
}

void epartition::receive_halo() {
	
	comm_wait(boundary::partitions,0,boundary::symmetric);
	
	/* Now unpack what we have received */
	int cnt = 0;
	remote_halo.npnt = static_cast<int>(frcvbuf(0,cnt++));
	remote_halo.nseg = static_cast<int>(frcvbuf(0,cnt++));
	remote_halo.ntri = static_cast<int>(frcvbuf(0,cnt++));
	for(int j=0;j<remote_halo.npnt;++j) {
		remote_halo.pnts(j)(0) = frcvbuf(0,cnt++);
		remote_halo.pnts(j)(1) = frcvbuf(0,cnt++);
	}
	for(int j=0;j<remote_halo.nseg;++j) {
		remote_halo.seg(j).pnt(0) = static_cast<int>(frcvbuf(0,cnt++));
		remote_halo.seg(j).pnt(1) = static_cast<int>(frcvbuf(0,cnt++));
	}
	
	for(int j=0;j<remote_halo.ntri;++j) {
		remote_halo.tri(j).pnt(0) = static_cast<int>(frcvbuf(0,cnt++));
		remote_halo.tri(j).pnt(1) = static_cast<int>(frcvbuf(0,cnt++));
		remote_halo.tri(j).pnt(2) = static_cast<int>(frcvbuf(0,cnt++));
	}
	remote_halo.createsegtri();
	
	remote_halo.ebdry(0)->nseg = nseg;
	cnt = 0;
	for(int i=0;i<nseg;++i) {
		remote_halo.ebdry(0)->seg(i) = cnt++;
	}
	
	remote_halo.ebdry(1)->nseg = 0;
	for(int i=nseg;i<remote_halo.nseg;++i) {
		if (remote_halo.seg(i).tri(1) < 0) {
			remote_halo.ebdry(1)->seg(remote_halo.ebdry(1)->nseg++) = i;
		}
	}
	
	remote_halo.bdrylabel();
	remote_halo.createtritri();
	remote_halo.createpnttri();
	remote_halo.cnt_nbor();
	remote_halo.treeinit();
	
	
	if (x.gbl->adapt_output) {	
		std::string adapt_file;
		std::ostringstream nstr;
		nstr << x.gbl->tstep << std::flush;
		adapt_file = "halo" +nstr.str() +"_" +idprefix;
		remote_halo.output(adapt_file);
	}

}
