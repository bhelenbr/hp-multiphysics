//
//  partition.cpp
//  tri_hp
//
//  Created by Brian Helenbrook on 4/7/16.
//
//

#include <stdio.h>
#include "tri_hp.h"
#include "hp_boundary.h"


void tri_hp::partition(multigrid_interface& xmesh, int npart, int maxenum, int maxvnum) {
	tri_hp& xin(dynamic_cast<tri_hp&>(xmesh));
	
	// Delete existing boundaries
	for(int i=0;i<nebd;++i) {
		delete hp_ebdry(i);
	}
	
	// Delete existing boundaries
	for(int i=0;i<nvbd;++i) {
		delete hp_vbdry(i);
	}
	
	tri_mesh::partition(xmesh,npart,maxenum,maxvnum);

	for(int t=0;t<gbl->nadapt;++t) {
		for(int i=0;i<npnt;++i) {
			ugbd(t).v(i,Range::all()) = xin.ugbd(t).v(pnt(i).info,Range::all());
			
			for(int n=0;n<ND;++n)
				vrtxbd(t)(i)(n) = xin.vrtxbd(t)(pnt(i).info)(n);
		}
		
		if (sm0) {
			for(int i=0;i<nseg;++i) {
				if (pnt(seg(i).pnt(0)).info == xin.seg(seg(i).info).pnt(0)) {
					ugbd(t).s(i,Range::all(),Range::all()) = xin.ugbd(t).s(seg(i).info,Range::all(),Range::all());
				}
				else {
					int msgn = 1;
					for(int m=0;m<sm0;++m) {
						ugbd(t).s(i,m,Range::all()) = msgn*xin.ugbd(t).s(seg(i).info,m,Range::all());
						msgn *= -1;
					}
				}
			}
		
			if (im0) {
				for(int i=0;i<ntri;++i) {
					ugbd(t).i(i,Range::all(),Range::all()) = xin.ugbd(t).i(tri(i).info,Range::all(),Range::all());
				}
			}
		}
	}
	
	hp_ebdry.resize(nebd);
	for(int i=0;i<nebd;++i) {
		int sind = seg(ebdry(i)->seg(0)).info;
		if (xin.seg(sind).tri(1) < 0) {
			/* physical boundary */
			int bnum = getbdrynum(xin.seg(sind).tri(1));
			hp_ebdry(i) = xin.hp_ebdry(bnum)->create(*this,*ebdry(i));
			for (int k=0;k<ebdry(i)->nseg;++k) {
				int sind = ebdry(i)->seg(k);
				int tgt = seg(sind).info;
				tgt = xin.getbdryseg(xin.seg(tgt).tri(1));
				hp_ebdry(i)->movesdata_bdry(k, xin.hp_ebdry(bnum), tgt);
			}
		}
		else {
			hp_ebdry(i) = new hp_edge_bdry(*this,*ebdry(i));
		}
	}
	
	hp_vbdry.resize(nvbd);
	for(int i=0;i<nvbd;++i) {
		int p0 = pnt(vbdry(i)->pnt).info;
		for(int j=0;j<xin.nvbd;++j) {
			if (xin.vbdry(j)->pnt == p0) {
				hp_vbdry(i) = xin.hp_vbdry(j)->create(*this,*vbdry(i));
				goto next1;
			}
		}
		/* New vertex boundary must be communication type */
		/* Try to match physics of side??? */
		hp_vbdry(i) = new hp_vrtx_bdry(*this,*vbdry(i));
		
	next1: continue;
	}
}
