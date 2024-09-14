#include "tri_hp.h"
#include "hp_boundary.h"

void tri_hp::mg_restrict() {
	int i,j,k,m,n,tind,v0,indx,indx1;

	isfrst = true;
	coarse_flag = true;

	for(i=0;i<nebd;++i)
		hp_ebdry(i)->mg_restrict();

	helper->mg_restrict();

	if(!coarse_level) {
		--log2p;

		/* TRANSFER IS ON FINE MESH */
		hp_gbl->res0.v(Range(0,npnt-1),Range::all()) = hp_gbl->res.v(Range(0,npnt-1),Range::all());

		if (basis::tri(log2p)->p() > 1) {
			hp_gbl->res0.s(Range(0,nseg-1),Range(0,basis::tri(log2p)->sm()-1),Range::all()) = hp_gbl->res.s(Range(0,nseg-1),Range(0,basis::tri(log2p)->sm()-1),Range::all());

			if (basis::tri(log2p)->p() > 2) {

				for(tind=0;tind<ntri;++tind) {
					indx = 0;
					indx1 = 0;
					for(m=1;m<basis::tri(log2p)->sm();++m) {
						for(k=0;k<basis::tri(log2p)->sm()-m;++k) {
							hp_gbl->res0.i(tind,indx,Range::all()) = hp_gbl->res.i(tind,indx1,Range::all());
							++indx;
							++indx1;
						}
						indx1 += basis::tri(log2p)->p();
					}
				}
			}
		}
	}
	else {
		if (mmovement == coupled_deformable) { 
			r_tri_mesh::mg_restrict();
		}

		tri_hp *fmesh = dynamic_cast<tri_hp *>(fine);

		hp_gbl->res0.v(Range(0,npnt-1),Range::all()) = 0.0;

		/* LOOP THROUGH FINE VERTICES TO CALCULATE RESIDUAL  */
		for(i=0;i<fmesh->npnt;++i) {
			tind = fmesh->ccnnct(i).tri;
			for(j=0;j<3;++j) {
				v0 = tri(tind).pnt(j);
				for(n=0;n<NV;++n)
					hp_gbl->res0.v(v0,n) += fmesh->ccnnct(i).wt(j)*hp_gbl->res.v(i,n);
			}
		}
		
		/* Need to communicate for partition boundaries */
		for(int last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
			for(i=0;i<nvbd;++i)
				vbdry(i)->vloadbuff(boundary::partitions,hp_gbl->res0.v.data(),0,NV-1,NV);
			for(i=0;i<nebd;++i)
				ebdry(i)->vloadbuff(boundary::partitions,hp_gbl->res0.v.data(),0,NV-1,NV);
			
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
				ebdry(i)->vfinalrcv(boundary::partitions,mp_phase,boundary::symmetric,boundary::sum,hp_gbl->res0.v.data(),0,NV-1,NV);
			}
			for(i=0;i<nvbd;++i) {
				vbdry(i)->comm_wait(boundary::partitions,0,boundary::symmetric);
				vbdry(i)->vfinalrcv(boundary::partitions,0,boundary::symmetric,boundary::average,hp_gbl->res0.v.data(),0,NV-1,NV);
			}
		}	
		

		/* LOOP THROUGH COARSE VERTICES    */
		/* TO CALCULATE VUG ON COARSE MESH */
		for(i=0;i<npnt;++i) {
			tind = fcnnct(i).tri;

			ug.v(i,Range::all()) = 0.0;

			for(j=0;j<3;++j) {
				ug.v(i,Range::all()) += fcnnct(i).wt(j)*fmesh->ug.v(fmesh->tri(tind).pnt(j),Range::all());
			}

			vug_frst(i,Range::all()) = ug.v(i,Range::all());
		}
	}


	return;
}

