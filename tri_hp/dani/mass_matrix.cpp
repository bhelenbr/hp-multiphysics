//
//  mass_matrix.cpp
//  tri_hp
//
//  Created by Brian Helenbrook on 3/30/17.
//
//

#include <stdio.h>
#include "tri_hp_dani.h"

void tri_hp_dani::create_mass_matrix(sparse_row_major& mass) {
	
	/* Calculate Storage for Sparse Mass Matrix */
	const int sm=basis::tri(log2p)->sm();
	const int im=basis::tri(log2p)->im();
	const int tm=basis::tri(log2p)->tm();
	const int lgpx = basis::tri(log2p)->gpx(), lgpn = basis::tri(log2p)->gpn();
	
	/* rank of mass matrix */
	const int size = npnt +nseg*sm +ntri*im;
	
	/* find number of non-zeros for each row */
	Array<int,1> nnzero(size);
	const int begin_seg = npnt;
	const int begin_tri = begin_seg+nseg*sm;
	const int ndofs = begin_tri +ntri*im;
	
	/* SELF CONNECTIONS */
	nnzero(Range(0,begin_seg-1)) = 1;
	if (sm) {
		nnzero(Range(begin_seg,begin_tri-1)) = (2 +sm);
	}
	if (im) nnzero(Range(begin_tri,ndofs-1)) = tm;
	
	/* edges and vertices connected to a vertex */
	for(int i=0; i<npnt; ++i)
		nnzero(i) += pnt(i).nnbor*(sm+1);
	
	for(int i=0; i<ntri; ++i) {
		/* interior and opposing side mode for each vertex */
		for(int j=0;j<3;++j)
			nnzero(tri(i).pnt(j)) += (im+sm);
		
		
		/* interior modes,opposing side modes, and opposing vertex to each side */
		for(int j=0;j<3;++j)
			for(int m=0;m<sm;++m)
				nnzero(begin_seg+tri(i).seg(j)*sm+m) += im +2*sm +1;
	}
	mass.resize(size,nnzero);
	
	/* Calculate Sparse Mass Matrix */
	uht(0) = 0.0;
	
	for(int tind=0;tind<ntri;++tind) {
		
		/* Create vector of global indices */
		Array<int,1> gindx(tm),gsign(tm);
		gsign=1;
		
		/* VERTEX MODES */
		int indx = 0;
		for (int m=0; m<3; ++m) {
			gindx(indx++) = tri(tind).pnt(m);
		}
		
		if (sm > 0) {
			/* SIDE MODES */
			for(int i=0;i<3;++i) {
				int sind = npnt +tri(tind).seg(i)*sm;
				int sgn = tri(tind).sgn(i);
				int msgn = 1;
				for (int m = 0; m < sm; ++m) {
					gindx(indx) = sind;
					gsign(indx++) = msgn;
					msgn *= sgn;
					++sind;
				}
			}
			
			int iind = npnt +nseg*sm +tind*im;
			for(int m=0;m<im;++m) {
				gindx(indx++) = iind++;
			}
		}
		
		/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
		crdtocht(tind);
		
		/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
		for(int n=0;n<ND;++n)
			basis::tri(log2p)->proj_bdry(&cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), &dcrd(n,1)(0,0),MXGP);
		
		for(int i=0;i<lgpx;++i) {
			for(int j=0;j<lgpn;++j) {
				cjcb(i,j) = sin(crd(0)(i,j))*(dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j));
			}
		}
		
		
		for(int m=0;m<tm;++m) {
			uht(0)(m) = 1.0*gsign(m);
			basis::tri(log2p)->proj(&uht(0)(0),&u(0)(0,0),MXGP);
			
			
			for(int i=0;i<lgpx;++i) {
				for(int j=0;j<lgpn;++j) {
					u(0)(i,j) *= cjcb(i,j);
				}
			}
			basis::tri(log2p)->intgrt(&lf(0)(0),&u(0)(0,0),MXGP);
			uht(0)(m) = 0.0;
			
			/* store in mass matrix */
			for(int k=0;k<tm;++k)
				mass.add_values(gindx(m),gindx(k),gsign(k)*lf(0)(k));
		}
	}
	mass.check_for_unused_entries(*gbl->log);
	
	return;
}
