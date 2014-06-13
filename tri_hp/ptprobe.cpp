#include "tri_hp.h"
#include "hp_boundary.h"
#include <assert.h>

bool tri_hp::ptprobe(TinyVector<FLT,2> xp, Array<FLT,1> uout, int& tind, int tlvl) {
	FLT r,s;
	bool found;

	found = findinteriorpt(xp,tind,r,s);
	ugtouht(tind,tlvl);  
	basis::tri(log2p)->ptprobe(NV,uout.data(),&uht(0)(0),MXTM);
	return(found);
	
}

void tri_hp::findandmvptincurved(TinyVector<FLT,2>& xp, int &tind, FLT &r, FLT &s) {
	TinyVector<FLT,3> wgt;
	int v0;
	bool found;

	qtree.nearpt(xp,v0);
	found = findtri(xp,v0,tind);
	if (!found) {
		*gbl->log << "#Warning: couldn't find tri " << xp << " nearpt " << v0 << " neartri " << tind << std::endl;
	}

	getwgts(wgt);
	/* TRIANGLE COORDINATES */    
	s = wgt(0)*2 -1.0;
	r = wgt(2)*2 -1.0;

	if (tri(tind).info < 0) {
		basis::tri(log2p)->ptvalues_rs(r,s);
		return;
	}

	/* MOVE POINT WITH SIDE CURVATURE */
	crdtocht(tind);
	basis::tri(log2p)->ptprobe_bdry(ND,xp.data(),r,s,&cht(0,0),MXTM);

	/* need to do this because ptprobe_bdry only calculates boundary function */
	basis::tri(log2p)->ptvalues_rs(r,s);

	return;
}

int mistake_counter = 0;

bool tri_hp::findinteriorpt(TinyVector<FLT,ND> xp, int &tind, FLT &r, FLT &s) {
	FLT dr,ds,dx,dy,det,roundoff;
	TinyVector<FLT,3> wgt;
	TinyVector<FLT,ND> x,dxmax,ddr,dds;
	int n,iter,v0;
	bool found = true;
	const FLT ubound = 1.0+10.0*FLT_EPSILON;

	if (tind < 0) {
		qtree.nearpt(xp,v0);
		found = findtri(xp,v0,tind);
	}
	else {
		found = findtri(xp,tind);
	}
	getwgts(wgt);

	/* TRIANGLE COORDINATES */    
	s = wgt(0)*2 -1.0;
	r = wgt(2)*2 -1.0;

	if (tri(tind).info >= 0) {
		/* DEAL WITH CURVED SIDES */
		crdtocht(tind);

		for(n=0;n<ND;++n)
			dxmax(n) = fabs(cht(n,0)-cht(n,1)) +fabs(cht(n,1)-cht(n,2));
		roundoff = 10.0*EPSILON*(1.0 +(fabs(xp(0))*dxmax(1) +fabs(xp(1))*dxmax(0))/area(tind));

		iter = 0;
		do {
			basis::tri(log2p)->ptprobe_bdry(ND,x.data(),ddr.data(),dds.data(),r,s,&cht(0,0),MXTM);
			det = 1.0/(fabs(ddr(0)*dds(1) - ddr(1)*dds(0)) +10.0*EPSILON);
			dx = xp(0)-x(0);
			dy = xp(1)-x(1);
			dr =  (dds(1)*dx -dds(0)*dy)*det;
			ds = -(ddr(1)*dx -ddr(0)*dy)*det;

			r += dr;
			s += ds;
			if (iter++ > 100) {
				*gbl->log << "#Warning: max iterations for curved triangle " << tind << "find tri?" << found << " from near pt " << v0 << " loc: " << xp << " x: " << x << " r: " << r << " s: " << s << " dr: " << dr << " ds: " << ds <<std::endl;
				/* TRIANGLE COORDINATES */    
				s = wgt(0)*2 -1.0;
				r = wgt(2)*2 -1.0;
				*gbl->log  << "#Warning: this was the first guess " << r << ' ' << s << ' ' << '\n';
				basis::tri(log2p)->ptvalues_rs(r,s);
				return(false);
			}
		} while (fabs(dr) +fabs(ds) > roundoff);
			
		
		if (fabs(r+1.0) +fabs(s+1.0) +fabs(-r-s) -2.0 < roundoff ) {
			found = true;
		}
		else {
			*gbl->log << "#Warning: point outside triangle " << tind << "find tri?" << found << " loc: " << xp << " x: " << x << " r: " << r << " s: " << s << " dr: " << dr << " ds: " << ds <<std::endl;
			found = false;
		}

		/* need to do this because ptprobe_bdry only calculates boundary function */
		basis::tri(log2p)->ptvalues_rs(r,s);

		return(found);
	}
	
	/* CHECK FOR STRAIGHT EDGED TRIANGLES */
	if (!found) {
		*gbl->log << "#Warning point outside of straight edged triangle " << tind << " loc: " << xp << " r: " << r << " s: " << s << std::endl;
	}
	/* need to do this because ptprobe_bdry only calculates boundary function */
	basis::tri(log2p)->ptvalues_rs(r,s);
	return(found);
}
