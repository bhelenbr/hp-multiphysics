#include "tri_hp.h"
#include "hp_boundary.h"


//#define DEBUG

void tri_hp::mg_prolongate() {
	int i,j,ind,tind;

	if(!coarse_level) {
		++log2p;
		if (log2p == log2pmax) coarse_flag = false;
		return;
	}

	if (mmovement == coupled_deformable) {
		r_tri_mesh::mg_prolongate();
	}

	for(i=0;i<nebd;++i)
		hp_ebdry(i)->mg_prolongate();

    /* CALCULATE CORRECTIONS */
    vug_frst(Range(0,npnt-1),Range::all()) -= ug.v(Range(0,npnt-1),Range::all());

#ifdef DEBUG
    *gbl->log << vug_frst(Range(0,npnt-1),Range::all());
#endif

    /* LOOP THROUGH FINE VERTICES    */
    /* TO DETERMINE CHANGE IN SOLUTION */   
    tri_hp *fmesh = dynamic_cast<tri_hp *>(fine);
    int fnvrtx = fmesh->npnt;
    for(i=0;i<fnvrtx;++i) {
		tind = fmesh->ccnnct(i).tri;

		gbl->res.v(i,Range::all()) = 0.0;

		for(j=0;j<3;++j) {
			ind = tri(tind).pnt(j);
			gbl->res.v(i,Range::all()) -= fmesh->ccnnct(i).wt(j)*vug_frst(ind,Range::all());
		}
	}

#ifdef DEBUG
	*gbl->log << gbl->res.v(Range(0,fnvrtx-1),Range::all());
#endif


	/* Need to communicate for partition boundaries */
	for(int last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		for(i=0;i<nvbd;++i)
			vbdry(i)->vloadbuff(boundary::partitions,gbl->res0.v.data(),0,NV-1,NV);
		for(i=0;i<nebd;++i)
			ebdry(i)->vloadbuff(boundary::partitions,gbl->res0.v.data(),0,NV-1,NV);
		
		for(i=0;i<nebd;++i) 
			ebdry(i)->comm_prepare(boundary::partitions,mp_phase,boundary::symmetric);
		for(i=0;i<nvbd;++i)
			vbdry(i)->comm_prepare(boundary::partitions,0,boundary::symmetric);
		
		for(i=0;i<nebd;++i) 
			ebdry(i)->comm_exchange(boundary::partitions,mp_phase,boundary::symmetric);
		for(i=0;i<nvbd;++i)
			vbdry(i)->comm_exchange(boundary::partitions,0,boundary::symmetric);
		
		last_phase = true;
		for(i=0;i<nebd;++i) {
			last_phase &= ebdry(i)->comm_wait(boundary::partitions,mp_phase,boundary::symmetric);
			ebdry(i)->vfinalrcv(boundary::partitions,mp_phase,boundary::symmetric,boundary::sum,gbl->res0.v.data(),0,NV-1,NV);
		}
		for(i=0;i<nvbd;++i) {
			vbdry(i)->comm_wait(boundary::partitions,0,boundary::symmetric);
			vbdry(i)->vfinalrcv(boundary::partitions,0,boundary::symmetric,boundary::average,gbl->res0.v.data(),0,NV-1,NV);
		}
	}

	/* ADD CORRECTION */
	fmesh->ug.v(Range(0,fnvrtx-1),Range::all()) += gbl->res.v(Range(0,fnvrtx-1),Range::all());      

	return;
}






