#include "tet_hp.h"
#include "boundary.h"
#include "hp_boundary.h"
#include <assert.h>

void tet_hp::ptprobe(TinyVector<FLT,3> xp, Array<FLT,1> uout, int tlvl) {
	FLT r,s,t;
	int tind;
	
	findinteriorpt(xp,tind,r,s,t);
	ugtouht(tind,tlvl);  
	basis::tet(log2p).ptprobe(NV,uout.data(),&uht(0)(0),MXTM);
}

//void tet_hp::ptprobe_bdry(int bnum, TinyVector<FLT,2> xp, Array<FLT,1> uout,int tlvl) {
//	FLT psi;
//	int sind;
//	
//	hp_ebdry(bnum)->findandmovebdrypt(xp,sind,psi);
//	sind = ebdry(bnum)->seg(sind);
//	ugtouht1d(sind,tlvl);  
//	basis::tri(log2p).ptprobe1d(NV,uout.data(),&uht(0)(0),MXTM);
//}

//void tet_hp::findandmvptincurved(TinyVector<FLT,2>& xp, int &tind, FLT &r, FLT &s) {
//	TinyVector<FLT,3> wgt;
//	int v0;
//	bool found;
//	
//	qtree.nearpt(xp.data(),v0);
//	found = findtri(xp,v0,tind);
//	if (!found) {
//		*gbl->log << "#Warning: couldn't find tri " << xp << " nearpt " << v0 << " neartri " << tind << std::endl;
//	}
//	
//	 getwgts(wgt);
//	/* TRIANGLE COORDINATES */    
//	s = wgt(0)*2 -1.0;
//	r = wgt(2)*2 -1.0;
//	
//	if (tri(tind).info < 0) {
//		basis::tri(log2p).ptvalues_rs(r,s);
//		return;
//	}
//	
//	/* MOVE POINT WITH SIDE CURVATURE */
//	crdtocht(tind);
//	basis::tri(log2p).ptprobe_bdry(ND,xp.data(),r,s,&cht(0,0),MXTM);
//		
//	/* need to do this because ptprobe_bdry only calculates boundary function */
//	basis::tri(log2p).ptvalues_rs(r,s);
//
//	return;
//}

int mistake_counter = 0;

int tet_hp::findinteriorpt(TinyVector<FLT,ND> xp, int &tind, FLT &r, FLT &s, FLT &t) {
//	FLT dr,ds,dx,dy,det,roundoff;
//	TinyVector<FLT,ND> x,dxmax,ddr,dds;
//	int n,iter,tind1;

	int ierr = 0;
	TinyVector<FLT,4> wgt;
	int v0;
	bool found;

	otree.nearpt(xp.data(),v0);
	found = findtet(xp,v0,tind);
	getwgts(wgt);

	/* TRIANGLE COORDINATES */    
	t = 2.0*wgt(0)-1.0;
	s = 2.0*wgt(1)-1.0;
	r = 2.0*wgt(3)-1.0;
	
//	if (tri(tind).info >= 0) {
//		/* DEAL WITH CURVED SIDES */
//		crdtocht(tind);
//		
//		for(n=0;n<ND;++n)
//			dxmax(n) = fabs(cht(n,0)-cht(n,1)) +fabs(cht(n,1)-cht(n,2));
//		roundoff = 10.0*EPSILON*(1.0 +(fabs(xp(0))*dxmax(1) +fabs(xp(1))*dxmax(0))/area(tind));
//		
//		iter = 0;
//		do {
//			basis::tri(log2p).ptprobe_bdry(ND,x.data(),ddr.data(),dds.data(),r,s,&cht(0,0),MXTM);
//			det = 1.0/(fabs(ddr(0)*dds(1) - ddr(1)*dds(0)) +10.0*EPSILON);
//			dx = xp(0)-x(0);
//			dy = xp(1)-x(1);
//			dr =  (dds(1)*dx -dds(0)*dy)*det;
//			ds = -(ddr(1)*dx -ddr(0)*dy)*det;
//
//			r += dr;
//			s += ds;
//			if (iter++ > 100) {
//				*gbl->log << "#Warning: max iterations for curved triangle " << tind << "find tri?" << found << " from near pt " << v0 << " loc: " << xp << " x: " << x << " r: " << r << " s: " << s << " dr: " << dr << " ds: " << ds <<std::endl;
//				std::ostringstream fname;
//				fname << "target_solution" << gbl->tstep << '_' << gbl->idprefix;
//				tet_mesh::output(fname.str().c_str(),tet_mesh::grid);
//				tet_hp::output(fname.str().c_str(),tet_hp::tecplot);
//				/* TRIANGLE COORDINATES */    
//				s = wgt(0)*2 -1.0;
//				r = wgt(2)*2 -1.0;
//				*gbl->log  << "#Warning: this was the first guess " << r << ' ' << s << ' ' << '\n';
//				ierr = 1;
//				break;
//			}
//		} while (fabs(dr) +fabs(ds) > roundoff);
//		
//		if (r < -(1.0+10.0*FLT_EPSILON) || r > (1.0+10.0*FLT_EPSILON) || s < -(1.0+10.0*FLT_EPSILON) || s > (1.0+10.0*FLT_EPSILON)) {
//			*gbl->log << "#Warning: point outside triangle " << tind << "find tri?" << found << " loc: " << xp << " x: " << x << " r: " << r << " s: " << s << " dr: " << dr << " ds: " << ds <<std::endl;
//			std::ostringstream fname;
//			fname << "target_solution" << gbl->tstep << '_' << gbl->idprefix;
//			tet_mesh::output(fname.str().c_str(),tet_mesh::grid);
//			tet_hp::output(fname.str().c_str(),tet_hp::tecplot);
//			ierr = 1;
//		}
//		/* need to do this because ptprobe_bdry only calculates boundary function */
//		basis::tri(log2p).ptvalues_rs(r,s);
//
//		return(ierr);
//	}
//	else if (tind1 < 0) {
//		*gbl->log << "#Warning point outside of straight edged triangle " << tind << " loc: " << xp << " x: " << x << " r: " << r << " s: " << s << std::endl;
//		ierr = 1;
//	}
		
	/* need to do this because ptprobe_bdry only calculates boundary function */
	basis::tet(log2p).ptvalues_rst(r,s,t);
	return(ierr);
}
