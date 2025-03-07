/*
 *  copy.cpp
 *  planar++
 *
 *  Created by helenbrk on Thu Oct 25 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp.h"
#include "hp_boundary.h"

void tri_hp::copy(const tri_hp& tgt) {
	int i,n,t;

	/* COPY MESH INFORMATION */
	tri_mesh::copy(tgt);

	for(t=0;t<gbl->nadapt;++t) {
		ugbd(t).v(Range(0,npnt-1),Range::all()) = tgt.ugbd(t).v(Range(0,npnt-1),Range::all());
		if (sm0) ugbd(t).s(Range(0,nseg-1),Range::all(),Range::all()) = tgt.ugbd(t).s(Range(0,nseg-1),Range::all(),Range::all());
		if (im0) ugbd(t).i(Range(0,ntri-1),Range::all(),Range::all()) = tgt.ugbd(t).i(Range(0,ntri-1),Range::all(),Range::all());

		for(i=0;i<npnt;++i)
			for(n=0;n<ND;++n)
				vrtxbd(t)(i)(n) = tgt.vrtxbd(t)(i)(n);
        
#ifdef ALLCURVED
        if (allcurved) {
            if (sm0) crvbd(t).s(Range(0,nseg-1),Range::all(),Range::all()) = tgt.crvbd(t).s(Range(0,nseg-1),Range::all(),Range::all());
            if (im0) crvbd(t).i(Range(0,ntri-1),Range::all(),Range::all()) = tgt.crvbd(t).i(Range(0,ntri-1),Range::all(),Range::all());
        }
#endif
	}
    pmetric.reset(tgt.pmetric->create(*this).release());
    
	for(i=0;i<nebd;++i)
		hp_ebdry(i)->copy(*tgt.hp_ebdry(i));

	for(i=0;i<nvbd;++i)
		hp_vbdry(i)->copy(*tgt.hp_vbdry(i));

	return;
}

void tri_hp::transfer_halo_solutions() {
	
	for(int i=0;i<nebd;++i) {
		if (hp_partition *tgt = dynamic_cast<hp_partition *>(hp_ebdry(i))) {
			tgt->snd_solution();
		}
	}
	
	smsgpass(boundary::partitions, 0, boundary::symmetric);
	
	for(int i=0;i<nebd;++i) {
		if (hp_partition *tgt = dynamic_cast<hp_partition *>(hp_ebdry(i))) {
			tgt->rcv_solution();
		}
	}
}

void tri_hp::append_halos() {
	
	/* nebd can change during this loop because of appending process */
	int i = 0;
	while(i < nebd) {
		if (epartition *tgt = dynamic_cast<epartition *>(ebdry(i))) {
			/* Add solution values to mesh */
			/* If values are appended first then tri_mesh::append will automatically move them */
			hp_partition *hp_tgt = dynamic_cast<hp_partition *>(hp_ebdry(i));
			for(int j=0;j<gbl->nadapt;++j) {
				for(int i=0;i<tgt->remote_halo.npnt;++i) {
					for(int n=0;n<NV;++n) {
						ugbd(j).v(npnt+i,n) = hp_tgt->ugbd(j).v(i,n);
					}
					for(int n=0;n<ND;++n) {
						vrtxbd(j)(npnt+i)(n) = hp_tgt->vrtxbd(j)(i)(n);
					}
				}
			
				for(int i=0;i<tgt->remote_halo.nseg;++i) {
					for(int m=0;m<sm0;++m) {
						for(int n=0;n<NV;++n) {
							ugbd(j).s(nseg+i,m,n) = hp_tgt->ugbd(j).s(i,m,n);
						}
					}
				}
				
				for(int i=0;i<tgt->remote_halo.ntri;++i) {
					for(int m=0;m<im0;++m) {
						for(int n=0;n<NV;++n) {
							ugbd(j).i(ntri+i,m,n) = hp_tgt->ugbd(j).i(i,m,n);
						}
					}
				}
			}
			
			/* Append mesh */
			append(tgt->remote_halo);
			
			/* Make hp B.C.'s the same */
			delete hp_ebdry(i);
			for(int j=i;j<nebd-1;++j) {
				hp_ebdry(j) = hp_ebdry(j+1);
			}
			hp_ebdry(nebd-1) = getnewedgeobject(nebd-1,"plain");
			
			cleanup_after_adapt();

			continue;
		}
		++i;
	}
}





