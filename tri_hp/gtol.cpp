/*
 *  gtol.cpp
 *  planar++
 *
 *  Created by helenbrk on Sun Oct 14 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp.h"
#include "hp_boundary.h"
#include <assert.h>

void tri_hp::ugtouht(int tind) {
	int i,k,m,n,indx,sind,cnt;
	int sign, msgn;

	/* THIS IS FOR FLOW VARIABLES ON ANY MESH */
	/* VERTICES */    
	for (i=0; i<3; ++i) {
		indx = tri(tind).pnt(i);
		for(n=0; n<NV; ++n)
			uht(n)(i) = ug.v(indx,n);
	}


	/* SIDES */
	cnt = 3;
	for(i=0;i<3;++i) {
		sind = tri(tind).seg(i);
		sign = tri(tind).sgn(i);
		msgn = 1;
		for (m = 0; m < basis::tri(log2p)->sm(); ++m) {
			for(n=0; n<NV; ++n)
				uht(n)(cnt) = msgn*ug.s(sind,m,n);
			msgn *= sign;
			++cnt;
		}
	}

	/* INTERIORS */    
	if (basis::tri(log2p)->im() > 0) {
		indx = 0;
		for(m = 1; m < basis::tri(log2p)->sm(); ++m) {
			for(k = 0; k < basis::tri(log2p)->sm()-m; ++k) {
				for(n=0; n<NV; ++n)
					uht(n)(cnt) = ug.i(tind,indx,n);
				++cnt; ++indx;
			}
			indx += sm0 -basis::tri(log2p)->sm();
		}
	}

	return;
}

void tri_hp::ugtouht(int tind, int tlvl) {
	int i,k,m,n,sind,indx,cnt;
	int sign, msgn;
	vsi &ug = ugbd(tlvl);

	/* THIS IS FOR FLOW VARIABLES ON ANY MESH */
	/* VERTICES */    
	for (i=0; i<3; ++i) {
		indx = tri(tind).pnt(i);
		for(n=0; n<NV; ++n)
			uht(n)(i) = ug.v(indx,n);
	}

	/* SIDES */
	cnt = 3;
	for(i=0;i<3;++i) {
		sind = tri(tind).seg(i);
		sign = tri(tind).sgn(i);
		msgn = 1;
		for (m = 0; m < basis::tri(log2p)->sm(); ++m) {
			for(n=0; n<NV; ++n)
				uht(n)(cnt) = msgn*ug.s(sind,m,n);
			msgn *= sign;
			++cnt;
		}
	}

	/* INTERIORS */    
	if (basis::tri(log2p)->im() > 0) {    
		indx = 0;
		for(m = 1; m < basis::tri(log2p)->sm(); ++m) {
			for(k = 0; k < basis::tri(log2p)->sm()-m; ++k) {
				for(n=0; n<NV; ++n)
					uht(n)(cnt) = ug.i(tind,indx,n);
				++cnt; ++indx;
			}
			indx += sm0 -basis::tri(log2p)->sm();
		}
	}

	return;
}

 void tri_hp::ugtouht_bdry(int tind) {
	int i,m,n,indx,sind,cnt;
	int sign, msgn;

	for (i=0; i<3; ++i) {
		indx = tri(tind).pnt(i);
		for(n=0; n<NV; ++n)
			uht(n)(i) = ug.v(indx,n);
	}

	/* SIDES */
	cnt = 3;
	for(i=0;i<3;++i) {
		sind = tri(tind).seg(i);
		sign = tri(tind).sgn(i);
		msgn = 1;
		for (m = 0; m < basis::tri(log2p)->sm(); ++m) {
			for(n=0; n<NV; ++n)
				uht(n)(cnt) = msgn*ug.s(sind,m,n);
			msgn *= sign;
			++cnt;
		}
	}

	return;
}

 void tri_hp::ugtouht_bdry(int tind, int tlvl) {
	int i,m,n,indx,sind,cnt;
	int sign, msgn;
	vsi &ug = ugbd(tlvl);

	for (i=0; i<3; ++i) {
		indx = tri(tind).pnt(i);
		for(n=0; n<NV; ++n)
			uht(n)(i) = ug.v(indx,n);
	}

	/* SIDES */
	cnt = 3;
	for(i=0;i<3;++i) {
		sind = tri(tind).seg(i);
		sign = tri(tind).sgn(i);
		msgn = 1;
		for (m = 0; m < basis::tri(log2p)->sm(); ++m) {
			for(n=0; n<NV; ++n)
				uht(n)(cnt) = msgn*ug.s(sind,m,n);
			msgn *= sign;
			++cnt;
		}
	}

	return;
}


void tri_hp::ugtouht1d(int sind) {
	int m,n,v0,v1;

	v0 = seg(sind).pnt(0);
	v1 = seg(sind).pnt(1);
	for(n=0;n<NV;++n) {
		uht(n)(0) = ug.v(v0,n);
		uht(n)(1) = ug.v(v1,n);
	}

	for(m=0;m<basis::tri(log2p)->sm();++m)
		for(n=0;n<NV;++n) 
			uht(n)(m+2) = ug.s(sind,m,n);

	return;
}

void tri_hp::ugtouht1d(int sind, int tlvl) {
	int m,n,v0,v1;
	vsi &ug = ugbd(tlvl);

	v0 = seg(sind).pnt(0);
	v1 = seg(sind).pnt(1);
	for(n=0;n<NV;++n) {
		uht(n)(0) = ug.v(v0,n);
		uht(n)(1) = ug.v(v1,n);
	}

	for(m=0;m<basis::tri(log2p)->sm();++m)
		for(n=0;n<NV;++n) 
			uht(n)(m+2) = ug.s(sind,m,n);

	return;
}

void tri_hp::crdtocht(int tind) {
	int i,m,n,cnt,bnum,sind,indx;

	/* VERTICES */    
	for (i=0; i<3; ++i) {
		indx = tri(tind).pnt(i);
		for(n=0; n<ND; ++n)
			cht(n,i) = pnts(indx)(n);
	}

	if (basis::tri(log2p)->sm() == 0) return;

	/* SIDES */
	cnt = 3;
	for (i=0; i<3;++i) {    
		sind = tri(tind).seg(i);
		if (seg(sind).tri(1) >= 0) {
			for(m=0;m<basis::tri(log2p)->sm();++m) {
				for(n=0;n<ND;++n)
					cht(n,cnt) = 0.0;
				++cnt;
			}
		}
		else {
			bnum = getbdrynum(seg(sind).tri(1));
			indx = getbdryseg(seg(sind).tri(1));
			if (hp_ebdry(bnum)->is_curved()) {
				for (m = 0; m < basis::tri(log2p)->sm(); ++m) {
					for(n=0; n<ND; ++n)
						cht(n,cnt) = hp_ebdry(bnum)->crds(indx,m,n);
					++cnt;
				}
			}
			else {
				for(m=0;m<basis::tri(log2p)->sm();++m) {
					for(n=0;n<ND;++n)
						cht(n,cnt) = 0.0;
					++cnt;
				}
			}                
		}
	}

	return;
}

void tri_hp::crdtocht(int tind, int tlvl) {
	int i,m,n,cnt,bnum,sind,indx;

	/* VERTICES */    
	for (i=0; i<3; ++i) {
		indx = tri(tind).pnt(i);
		for(n=0; n<ND; ++n)
			cht(n,i) = vrtxbd(tlvl)(indx)(n);
	}

	if (basis::tri(log2p)->sm() == 0) return;

	/* SIDES */
	cnt = 3;
	for (i=0; i<3;++i) {    
		sind = tri(tind).seg(i);
		if (seg(sind).tri(1) >= 0) {
			for(m=0;m<basis::tri(log2p)->sm();++m) {
				for(n=0;n<ND;++n)
					cht(n,cnt) = 0.0;
				++cnt;
			}
		}
		else {
			bnum = getbdrynum(seg(sind).tri(1));
			indx = getbdryseg(seg(sind).tri(1));
			if (hp_ebdry(bnum)->is_curved()) {
				for (m = 0; m < basis::tri(log2p)->sm(); ++m) {
					for(n=0; n<ND; ++n)
						cht(n,cnt) = hp_ebdry(bnum)->crdsbd(tlvl,indx,m,n);
					++cnt;
				}
			}
			else {
				for(m=0;m<basis::tri(log2p)->sm();++m) {
					for(n=0;n<ND;++n)
						cht(n,cnt) = 0.0;
					++cnt;
				}
			} 
		}
	}

	return;
}

void tri_hp::crdtocht1d(int sind) {
	int m,n,bnum,indx,v0,v1;

	v0 = seg(sind).pnt(0);
	v1 = seg(sind).pnt(1);
	for(n=0;n<ND;++n) {
		cht(n,0) = pnts(v0)(n);
		cht(n,1) = pnts(v1)(n);
	}

	if (seg(sind).tri(1) >= 0) {
		for(m=0;m<basis::tri(log2p)->sm();++m)
			for(n=0;n<ND;++n) 
				cht(n,m+2) = 0.0;
	}
	else {
		bnum = getbdrynum(seg(sind).tri(1));
		indx = getbdryseg(seg(sind).tri(1));
		if (hp_ebdry(bnum)->is_curved()) {
			for(m=0;m<basis::tri(log2p)->sm();++m)
				for(n=0;n<ND;++n) 
					cht(n,m+2) = hp_ebdry(bnum)->crds(indx,m,n);
		}
		else {
			for(m=0;m<basis::tri(log2p)->sm();++m)
				for(n=0;n<ND;++n) 
					cht(n,m+2) = 0.0;
		}  
	}

	return;
}

void tri_hp::crdtocht1d(int sind,int tlvl) {
	int m,n,bnum,indx,v0,v1;

	v0 = seg(sind).pnt(0);
	v1 = seg(sind).pnt(1);
	for(n=0;n<ND;++n) {
		cht(n,0) = vrtxbd(tlvl)(v0)(n);
		cht(n,1) = vrtxbd(tlvl)(v1)(n);
	}

	if (seg(sind).tri(1) >= 0) {
		for(m=0;m<basis::tri(log2p)->sm();++m)
			for(n=0;n<ND;++n) 
				cht(n,m+2) = 0.0;
	}
	else {
		bnum = getbdrynum(seg(sind).tri(1));
		indx = getbdryseg(seg(sind).tri(1));
		if (hp_ebdry(bnum)->is_curved()) {
			for(m=0;m<basis::tri(log2p)->sm();++m)
				for(n=0;n<ND;++n) {
					cht(n,m+2) = hp_ebdry(bnum)->crdsbd(tlvl,indx,m,n);                    
				}
		}
		else {
			for(m=0;m<basis::tri(log2p)->sm();++m)
				for(n=0;n<ND;++n) 
					cht(n,m+2) = 0.0;
		}  
	}

	return;
}

void tri_hp::lftog(int tind, struct vsi g) {
	int i,m,n,indx,gindx,sind,sgn,msgn;

	/* VERTEX MODES */
	for (m=0; m<3; ++m) {
		gindx = tri(tind).pnt(m);
		for(n=0;n<NV;++n)
			g.v(gindx,n) += lf(n)(m);
	}


	if (basis::tri(log2p)->p() > 1) {
		/* SIDE MODES */
		indx = 3;
		for(i=0;i<3;++i) {
			sind = tri(tind).seg(i);
			sgn = tri(tind).sgn(i);
			msgn = 1;
			for (m = 0; m < basis::tri(log2p)->sm(); ++m) {
				for(n=0;n<NV;++n)
					g.s(sind,m,n) += msgn*lf(n)(indx);
				msgn *= sgn;
				++gindx;
				++indx;
			}
		}

		gindx = 0;
		indx = basis::tri(log2p)->bm();
		for(m=0;m<basis::tri(log2p)->im();++m) {
			for(n=0;n<NV;++n)
				g.i(tind,gindx,n) += lf(n)(indx);
			++gindx;
			++indx;
		}
	}

	return;
}
