/*
 *  tadvance.cpp
 *  tri_hp
 *
 *  Created by helenbrk on Wed Dec 05 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp.h"
#include "hp_boundary.h"

/* DIRK SCHEMES */
#ifdef DIRK
void tri_hp::tadvance() {
	int i,j,n,s,tind,stage;

	/* DO STUFF FOR DEFORMABLE MESH FIRST */
	if (log2p == log2pmax && gbl->substep == 0 && (mmovement == coupled_deformable || mmovement == uncoupled_deformable)) {
		r_tri_mesh::tadvance();
	}

	stage = gbl->substep +gbl->esdirk;
	if (!coarse_level) {
		if (stage > 0) {
			/* BACK CALCULATE K TERM */
			ugbd(stage+1).v(Range(0,npnt-1),Range::all()) = (ug.v(Range(0,npnt-1),Range::all()) -ugbd(1).v(Range(0,npnt-1),Range::all()))*gbl->adirk(stage-1,stage-1);
			if (basis::tri(log2p)->sm()) {
				ugbd(stage+1).s(Range(0,nseg-1),Range::all(),Range::all()) = (ug.s(Range(0,nseg-1),Range::all(),Range::all()) -ugbd(1).s(Range(0,nseg-1),Range::all(),Range::all()))*gbl->adirk(stage-1,stage-1);
				if (basis::tri(log2p)->im()) {
					ugbd(stage+1).i(Range(0,ntri-1),Range::all(),Range::all()) = (ug.i(Range(0,ntri-1),Range::all(),Range::all()) -ugbd(1).i(Range(0,ntri-1),Range::all(),Range::all()))*gbl->adirk(stage-1,stage-1);
				}
			}
			for(i=0;i<npnt;++i)
				for(n=0;n<ND;++n)
					vrtxbd(stage+1)(i)(n) = (pnts(i)(n)-vrtxbd(1)(i)(n))*gbl->adirk(stage-1,stage-1);
		}

		if (gbl->substep == 0) {
			/* STORE TILDE W */
			ugbd(1).v(Range(0,npnt-1),Range::all()) = ug.v(Range(0,npnt-1),Range::all());
			if (basis::tri(log2p)->sm()) {
				ugbd(1).s(Range(0,nseg-1),Range::all(),Range::all()) = ug.s(Range(0,nseg-1),Range::all(),Range::all());
				if (basis::tri(log2p)->im()) {
					ugbd(1).i(Range(0,ntri-1),Range::all(),Range::all()) = ug.i(Range(0,ntri-1),Range::all(),Range::all());
				}
			}

			/* SAME FOR MESH INFORMATION */
			for(i=0;i<npnt;++i)
				for(n=0;n<ND;++n)
					vrtxbd(1)(i)(n) = pnts(i)(n);
		}

		/* UPDATE TILDE W */
		for (s=0;s<stage;++s) {
			ugbd(1).v(Range(0,npnt-1),Range::all()) += gbl->adirk(stage,s)*ugbd(s+2).v(Range(0,npnt-1),Range::all());
			if (basis::tri(log2p)->sm()) {
				ugbd(1).s(Range(0,nseg-1),Range::all(),Range::all()) += gbl->adirk(stage,s)*ugbd(s+2).s(Range(0,nseg-1),Range::all(),Range::all());
				if (basis::tri(log2p)->im()) {
					ugbd(1).i(Range(0,ntri-1),Range::all(),Range::all()) += gbl->adirk(stage,s)*ugbd(s+2).i(Range(0,ntri-1),Range::all(),Range::all());
				}
			}
			for(i=0;i<npnt;++i) 
				for(n=0;n<ND;++n)
					vrtxbd(1)(i)(n) += gbl->adirk(stage,s)*vrtxbd(s+2)(i)(n);
		}

		/* EXTRAPOLATE */
		if (stage  && gbl->dti > 0.0) {
			FLT constant = gbl->extrapolate*gbl->cdirk(gbl->substep);  // Extrapolation seems to cause more trouble than it is worth
			ugbd(0).v(Range(0,npnt-1),Range::all()) += constant*ugbd(stage+1).v(Range(0,npnt-1),Range::all());
			if (basis::tri(log2p)->sm()) {
				ugbd(0).s(Range(0,nseg-1),Range::all(),Range::all()) += constant*ugbd(stage+1).s(Range(0,nseg-1),Range::all(),Range::all());
				if (basis::tri(log2p)->im()) {
					ugbd(0).i(Range(0,ntri-1),Range::all(),Range::all()) += constant*ugbd(stage+1).i(Range(0,ntri-1),Range::all(),Range::all());
				}
			}
			if (((mmovement == coupled_deformable) || (mmovement == uncoupled_deformable))) 
				vrtxbd(0)(Range(0,npnt-1)) += constant*vrtxbd(stage+1)(Range(0,npnt-1));               
		}
	}
	else {
		tri_hp* fmesh = dynamic_cast<tri_hp *>(fine);

		/* CALCULATE UNSTEADY SOURCE TERMS ON COARSE MESHES */
		for(i=0;i<npnt;++i) {
			tind = fcnnct(i).tri;

			ugbd(1).v(i,Range::all()) = 0.0;

			for(n=0;n<ND;++n)
				vrtxbd(1)(i)(n) = 0.0;


			for(j=0;j<3;++j) {
				ugbd(1).v(i,Range::all()) += fcnnct(i).wt(j)*fmesh->ugbd(1).v(fmesh->tri(tind).pnt(j),Range::all());
				for(n=0;n<ND;++n)
					vrtxbd(1)(i)(n) += fcnnct(i).wt(j)*fmesh->vrtxbd(1)(fmesh->tri(tind).pnt(j))(n);
			}
		}
	}

	for(i=0;i<nvbd;++i)
		hp_vbdry(i)->tadvance();

	for(i=0;i<nebd;++i) 
		hp_ebdry(i)->tadvance();

	helper->tadvance();

	calculate_unsteady_sources();
	
	
	
	return;

}

/* A GENERIC CALCULATION OF SOURCES FOR AN AUTONOMOUS SYSTEM IN STANDARD FORM */
/* WILL NEED TO BE OVERRIDDEN FOR SPECIAL CASES */
void tri_hp::calculate_unsteady_sources() {
	int i,j,n,tind;
#ifdef petsc
	int start = log2pmax;
#else
	int start = 0;
#endif
	
	for (log2p=start;log2p<=log2pmax;++log2p) {
		for(tind=0;tind<ntri;++tind) {
            pmetric->calc_metrics(tind, crd, dcrd, 1);
            
			ugtouht(tind,1);
			for(n=0;n<NV;++n)
				basis::tri(log2p)->proj(&uht(n)(0),&u(n)(0,0),MXGP);


			for(i=0;i<basis::tri(log2p)->gpx();++i) {
				for(j=0;j<basis::tri(log2p)->gpn();++j) {    
					cjcb(i,j) = -gbl->bd(0)*RAD(crd(0)(i,j))*(dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j));
					for(n=0;n<NV;++n)
						dugdt(log2p)(tind,n,i,j) = u(n)(i,j)*cjcb(i,j);
					for(n=0;n<ND;++n)
						dxdt(log2p)(tind,n,i,j) = crd(n)(i,j);
				}				
			}
		}
	}
	log2p = log2pmax;

	return;
}

void tri_hp::reset_timestep() {
    /* Rewind ugbd(1) to value from previous time step */
    for (int substep = gbl->substep; substep >= 0; --substep) {
        const int stage = substep +gbl->esdirk;
        for (int s=0;s<stage;++s) {
            ugbd(1).v(Range(0,npnt-1),Range::all()) -= gbl->adirk(stage,s)*ugbd(s+2).v(Range(0,npnt-1),Range::all());
            if (basis::tri(log2p)->sm()) {
                ugbd(1).s(Range(0,nseg-1),Range::all(),Range::all()) -= gbl->adirk(stage,s)*ugbd(s+2).s(Range(0,nseg-1),Range::all(),Range::all());
                if (basis::tri(log2p)->im()) {
                    ugbd(1).i(Range(0,ntri-1),Range::all(),Range::all()) -= gbl->adirk(stage,s)*ugbd(s+2).i(Range(0,ntri-1),Range::all(),Range::all());
                }
            }
            for(int i=0;i<npnt;++i)
                for(int n=0;n<ND;++n)
                    vrtxbd(1)(i)(n) -= gbl->adirk(stage,s)*vrtxbd(s+2)(i)(n);
        }
    }
    
    /* update ugbd(0) to be solution from previous time step */
    ug.v(Range(0,npnt-1),Range::all()) = ugbd(1).v(Range(0,npnt-1),Range::all());
    if (basis::tri(log2p)->sm()) {
        ug.s(Range(0,nseg-1),Range::all(),Range::all()) = ugbd(1).s(Range(0,nseg-1),Range::all(),Range::all());
        if (basis::tri(log2p)->im()) {
            ug.i(Range(0,ntri-1),Range::all(),Range::all()) = ugbd(1).i(Range(0,ntri-1),Range::all(),Range::all());
        }
    }

    /* SAME FOR MESH INFORMATION */
    for(int i=0;i<npnt;++i)
        for(int n=0;n<ND;++n)
            pnts(i)(n) = vrtxbd(1)(i)(n);
    
    
    const int stage = gbl->esdirk;
    if (stage) {
        /* Rewind ugbd(1) */
        ugbd(1).v(Range(0,npnt-1),Range::all()) = ug.v(Range(0,npnt-1),Range::all()) -ugbd(stage+1).v(Range(0,npnt-1),Range::all())/gbl->adirk(stage-1,stage-1);
        if (basis::tri(log2p)->sm()) {
            ugbd(1).s(Range(0,nseg-1),Range::all(),Range::all()) = ug.s(Range(0,nseg-1),Range::all(),Range::all()) -ugbd(stage+1).s(Range(0,nseg-1),Range::all(),Range::all())/gbl->adirk(stage-1,stage-1);
            if (basis::tri(log2p)->im()) {
                ugbd(1).i(Range(0,ntri-1),Range::all(),Range::all()) = ug.i(Range(0,ntri-1),Range::all(),Range::all()) -ugbd(stage+1).i(Range(0,ntri-1),Range::all(),Range::all())/gbl->adirk(stage-1,stage-1);
            }
        }
        for(int i=0;i<npnt;++i)
            for(int n=0;n<ND;++n)
                vrtxbd(1)(i)(n) = pnts(i)(n)-vrtxbd(stage+1)(i)(n)/gbl->adirk(stage-1,stage-1);
    }
    
    for(int i=0;i<nebd;++i)
        hp_ebdry(i)->reset_timestep();
}
#else

void tri_hp::tadvance() {
	int i,j,n,tind;

	/* DO STUFF FOR DEFORMABLE MESH FIRST */    
	if (log2p == log2pmax && gbl->substep == 0 && (mmovement == coupled_deformable || mmovement == uncoupled_deformable)) {
		r_tri_mesh::tadvance();
	}

	if (!coarse_level) {
		/* SHIFT DATA */
		for (int lvl=gbl->nhist-2;lvl>=0; --lvl) {            
			ugbd(lvl+1).v(Range(0,npnt-1),Range::all()) = ugbd(lvl).v(Range(0,npnt-1),Range::all());
			if (basis::tri(log2p)->sm()) {
				ugbd(lvl+1).s(Range(0,nseg-1),Range::all(),Range::all()) = ugbd(lvl).s(Range(0,nseg-1),Range::all(),Range::all());
				if (basis::tri(log2p)->im()) {
					ugbd(lvl+1).i(Range(0,ntri-1),Range::all(),Range::all()) = ugbd(lvl).i(Range(0,ntri-1),Range::all(),Range::all());
				}
			}

			/* SAME FOR MESH INFORMATION */
			for(i=0;i<npnt;++i)
				for(n=0;n<ND;++n)
					vrtxbd(lvl+1)(i)(n) = vrtxbd(lvl)(i)(n);
		}

		/* EXTRAPOLATE */
		if (gbl->dti > 0.0 && gbl->tstep > 1) {
			ugbd(0).v(Range(0,npnt-1),Range::all()) += ugbd(1).v(Range(0,npnt-1),Range::all()) -ugbd(2).v(Range(0,npnt-1),Range::all());
			if (basis::tri(log2p)->sm()) {
				ugbd(0).s(Range(0,nseg-1),Range::all(),Range::all()) += ugbd(1).s(Range(0,nseg-1),Range::all(),Range::all()) -ugbd(2).s(Range(0,nseg-1),Range::all(),Range::all());
				if (basis::tri(log2p)->im()) {
					ugbd(0).i(Range(0,ntri-1),Range::all(),Range::all()) += ugbd(1).i(Range(0,ntri-1),Range::all(),Range::all()) -ugbd(2).i(Range(0,ntri-1),Range::all(),Range::all());
				}
			}
			if (((mmovement == coupled_deformable) || (mmovement == uncoupled_deformable))) 
				vrtxbd(0)(Range(0,npnt-1)) += vrtxbd(1)(Range(0,npnt-1)) -vrtxbd(2)(Range(0,npnt-1));               
		}
	}
	else {
		tri_hp* fmesh = dynamic_cast<tri_hp *>(fine);

		/* CALCULATE UNSTEADY SOURCE TERMS ON COARSE MESHES */
		for(int lvl=0;lvl<gbl->nhist;++lvl) {
			for(i=0;i<npnt;++i) {
				tind = fcnnct(i).tri;

				ugbd(lvl).v(i,Range::all()) = 0.0;

				for(n=0;n<ND;++n)
					vrtxbd(lvl)(i)(n) = 0.0;


				for(j=0;j<3;++j) {
					ugbd(lvl).v(i,Range::all()) += fcnnct(i).wt(j)*fmesh->ugbd(lvl).v(fmesh->tri(tind).pnt(j),Range::all());
					for(n=0;n<ND;++n)
						vrtxbd(lvl)(i)(n) += fcnnct(i).wt(j)*fmesh->vrtxbd(lvl)(fmesh->tri(tind).pnt(j))(n);
				}
			}
		}
	}

	for(i=0;i<nvbd;++i) 
		hp_vbdry(i)->tadvance();

	for(i=0;i<nebd;++i) 
		hp_ebdry(i)->tadvance();

	helper->tadvance();

	calculate_unsteady_sources();

	return;

}

/* A GENERIC CALCULATION OF SOURCES FOR AN AUTONOMOUS SYSTEM IN STANDARD FORM */
/* WILL NEED TO BE OVERRIDDEN FOR SPECIAL CASES */
void tri_hp::calculate_unsteady_sources() {
	int i,j,n,tind;
#ifdef petsc
	int start = log2pmax;
#else
	int start = 0;
#endif
	
	for (log2p=start;log2p<=log2pmax;++log2p) {
		for (int level=1;level<min(gbl->nhist,gbl->tstep+1);++level) {
			for(tind=0;tind<ntri;++tind) {
                pmetric->calc_metrics(tind, crd, dcrd, level);

				ugtouht(tind,level);
				for(n=0;n<NV;++n)
					basis::tri(log2p)->proj(&uht(n)(0),&u(n)(0,0),MXGP);


				for(i=0;i<basis::tri(log2p)->gpx();++i) {
					for(j=0;j<basis::tri(log2p)->gpn();++j) {    
						cjcb(i,j) = -gbl->bd(level)*RAD(crd(0)(i,j))*(dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j));
						for(n=0;n<NV;++n)
							dugdt(log2p)(tind,n,i,j) += u(n)(i,j)*cjcb(i,j);
						for(n=0;n<ND;++n)
							dxdt(log2p)(tind,n,i,j) += gbl->bd(level)/gbl->bd(0)*crd(n)(i,j);
					}                
				}
			}
		}
	}

	log2p = log2pmax;

	return;
}
#endif


