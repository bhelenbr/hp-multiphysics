/*
 *  copy.cpp
 *  planar++
 *
 *  Created by helenbrk on Thu Oct 25 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp.h"
#include "hp_boundary.h"

void tet_hp::copy(const tet_hp& tgt) {
	int i,n,t;
	
	/* COPY MESH INFORMATION */
	tet_mesh::copy(tgt);

	for(t=0;t<gbl->nadapt;++t) {
		ugbd(t).v(Range(0,npnt-1),Range::all()) = tgt.ugbd(t).v(Range(0,npnt-1),Range::all());
		if (em0) ugbd(t).e(Range(0,nseg-1),Range::all(),Range::all()) = tgt.ugbd(t).e(Range(0,nseg-1),Range::all(),Range::all());
		if (fm0) ugbd(t).f(Range(0,ntri-1),Range::all(),Range::all()) = tgt.ugbd(t).f(Range(0,ntri-1),Range::all(),Range::all());
		if (im0) ugbd(t).i(Range(0,ntet-1),Range::all(),Range::all()) = tgt.ugbd(t).i(Range(0,ntet-1),Range::all(),Range::all());
		
		for(i=0;i<npnt;++i)
			for(n=0;n<ND;++n)
				vrtxbd(t)(i)(n) = tgt.vrtxbd(t)(i)(n);
	}

	for(i=0;i<nfbd;++i)
		hp_fbdry(i)->copy(*tgt.hp_fbdry(i));
		
	for(i=0;i<nebd;++i)
		hp_ebdry(i)->copy(*tgt.hp_ebdry(i));
		
	for(i=0;i<nvbd;++i)
		hp_vbdry(i)->copy(*tgt.hp_vbdry(i));
		
	return;
}

	


