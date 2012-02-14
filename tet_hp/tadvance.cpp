/*
 *  tadvance.cpp
 *  tet_hp
 *
 *  Created by helenbrk on Wed Dec 05 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp.h"
#include "hp_boundary.h"
#include <blitz/tinyvec-et.h>

/* DIRK SCHEMES */
#ifdef DIRK
void tet_hp::tadvance() {
	int i,j,n,s,tind,stage;
	
	/* DO STUFF FOR DEFORMABLE MESH FIRST */    
//    if (log2p == log2pmax && gbl->substep == 0 && (mmovement == coupled_deformable || mmovement == uncoupled_deformable)) {
//        r_tet_mesh::tadvance();
//    }

	stage = gbl->substep +gbl->esdirk;
	if (!coarse_level) {
		if (stage > 0) {
			/* BACK CALCULATE K TERM */
			ugbd(stage+1).v(Range(0,npnt-1),Range::all()) = (ug.v(Range(0,npnt-1),Range::all()) -ugbd(1).v(Range(0,npnt-1),Range::all()))*gbl->adirk(stage-1,stage-1);
			if (basis::tet(log2p).em) {
				ugbd(stage+1).e(Range(0,nseg-1),Range::all(),Range::all()) = (ug.e(Range(0,nseg-1),Range::all(),Range::all()) -ugbd(1).e(Range(0,nseg-1),Range::all(),Range::all()))*gbl->adirk(stage-1,stage-1);
				if (basis::tet(log2p).fm) {
					ugbd(stage+1).f(Range(0,ntri-1),Range::all(),Range::all()) = (ug.f(Range(0,ntri-1),Range::all(),Range::all()) -ugbd(1).f(Range(0,ntri-1),Range::all(),Range::all()))*gbl->adirk(stage-1,stage-1);
					if (basis::tet(log2p).im) {
						ugbd(stage+1).i(Range(0,ntet-1),Range::all(),Range::all()) = (ug.i(Range(0,ntet-1),Range::all(),Range::all()) -ugbd(1).i(Range(0,ntet-1),Range::all(),Range::all()))*gbl->adirk(stage-1,stage-1);
					}
				}
			}
			for(i=0;i<npnt;++i)
				for(n=0;n<ND;++n)
					vrtxbd(stage+1)(i)(n) = (pnts(i)(n)-vrtxbd(1)(i)(n))*gbl->adirk(stage-1,stage-1);
		}
		
		if (gbl->substep == 0) {
			/* STORE TILDE W */
			ugbd(1).v(Range(0,npnt-1),Range::all()) = ug.v(Range(0,npnt-1),Range::all());
			if (basis::tet(log2p).em) {
				ugbd(1).e(Range(0,nseg-1),Range::all(),Range::all()) = ug.e(Range(0,nseg-1),Range::all(),Range::all());
				if (basis::tet(log2p).fm) {
					ugbd(1).f(Range(0,ntri-1),Range::all(),Range::all()) = ug.f(Range(0,ntri-1),Range::all(),Range::all());
					if (basis::tet(log2p).im) {
						ugbd(1).i(Range(0,ntet-1),Range::all(),Range::all()) = ug.i(Range(0,ntet-1),Range::all(),Range::all());
					}
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
			if (basis::tet(log2p).em) {
				ugbd(1).e(Range(0,nseg-1),Range::all(),Range::all()) += gbl->adirk(stage,s)*ugbd(s+2).e(Range(0,nseg-1),Range::all(),Range::all());
				if (basis::tet(log2p).fm) {
					ugbd(1).f(Range(0,ntri-1),Range::all(),Range::all()) += gbl->adirk(stage,s)*ugbd(s+2).f(Range(0,ntri-1),Range::all(),Range::all());
					if (basis::tet(log2p).im) {
						ugbd(1).i(Range(0,ntet-1),Range::all(),Range::all()) += gbl->adirk(stage,s)*ugbd(s+2).i(Range(0,ntet-1),Range::all(),Range::all());
					}
				}
			}
			for(i=0;i<npnt;++i) 
				for(n=0;n<ND;++n)
					vrtxbd(1)(i)(n) += gbl->adirk(stage,s)*vrtxbd(s+2)(i)(n);
		}
		
		/* EXTRAPOLATE */
		if (stage  && gbl->dti > 0.0) {
			FLT constant =  0*gbl->cdirk(gbl->substep);  // Extrapolation seems to cause more trouble than it is worth
			ugbd(0).v(Range(0,npnt-1),Range::all()) += constant*ugbd(stage+1).v(Range(0,npnt-1),Range::all());
			if (basis::tet(log2p).em) {
				ugbd(0).e(Range(0,nseg-1),Range::all(),Range::all()) += constant*ugbd(stage+1).e(Range(0,nseg-1),Range::all(),Range::all());
				if (basis::tet(log2p).fm) {
					ugbd(0).f(Range(0,ntri-1),Range::all(),Range::all()) += constant*ugbd(stage+1).f(Range(0,ntri-1),Range::all(),Range::all());
					if (basis::tet(log2p).im) {
						ugbd(0).i(Range(0,ntet-1),Range::all(),Range::all()) += constant*ugbd(stage+1).i(Range(0,ntet-1),Range::all(),Range::all());
					}
				}
			}
			if (((mmovement == coupled_deformable) || (mmovement == uncoupled_deformable))) 
				vrtxbd(0)(Range(0,npnt-1)) += constant*vrtxbd(stage+1)(Range(0,npnt-1));               
		}
	}
	else {
		tet_hp* fmesh = dynamic_cast<tet_hp *>(fine);

		/* CALCULATE UNSTEADY SOURCE TERMS ON COARSE MESHES */
		for(i=0;i<npnt;++i) {
			tind = fcnnct(i).tet;

			ugbd(1).v(i,Range::all()) = 0.0;

			for(n=0;n<ND;++n)
				vrtxbd(1)(i)(n) = 0.0;

			for(j=0;j<4;++j) {
				ugbd(1).v(i,Range::all()) += fcnnct(i).wt(j)*fmesh->ugbd(1).v(fmesh->tet(tind).pnt(j),Range::all());
				for(n=0;n<ND;++n)
					vrtxbd(1)(i)(n) += fcnnct(i).wt(j)*fmesh->vrtxbd(1)(fmesh->tet(tind).pnt(j))(n);
			}
		}
	}


		

		
	for(i=0;i<nfbd;++i) 
		hp_fbdry(i)->tadvance();
	
	for(i=0;i<nebd;++i) 
		hp_ebdry(i)->tadvance();
	
	for(i=0;i<nvbd;++i) 
		hp_vbdry(i)->tadvance();
	
	helper->tadvance();

	calculate_unsteady_sources();

	return;
}

/* A GENERIC CALCULATION OF SOURCES FOR AN AUTONOMOUS SYSTEM IN STANDARD FORM */
/* WILL NEED TO BE OVERRIDDEN FOR SPECIAL CASES */
void tet_hp::calculate_unsteady_sources() {
	int lgpx = basis::tet(log2p).gpx, lgpy = basis::tet(log2p).gpy, lgpz = basis::tet(log2p).gpz;
	int i,j,k,n,tind;
	int stridey = MXGP;
	int stridex = MXGP*MXGP;
	TinyVector<int,4> v;

	for (log2p=0;log2p<=log2pmax;++log2p) {
		for(tind=0;tind<ntet;++tind) {
			v = tet(tind).pnt;

			if (tet(tind).info > -1) {
				crdtocht(tind,1);
				for(n=0;n<ND;++n)
					basis::tet(log2p).proj_bdry(&cht(n)(0), &crd(n)(0)(0)(0), &dcrd(n)(0)(0)(0)(0), &dcrd(n)(1)(0)(0)(0),&dcrd(n)(2)(0)(0)(0),stridex,stridey);
			}
			else {
				for(n=0;n<ND;++n)
					basis::tet(log2p).proj(vrtxbd(1)(v(0))(n),vrtxbd(1)(v(1))(n),vrtxbd(1)(v(2))(n),vrtxbd(1)(v(3))(n),&crd(n)(0)(0)(0),stridex,stridey);

				for(i=0;i<lgpx;++i) {
					for(j=0;j<lgpy;++j) {
						for(k=0;k<lgpz;++k) {
							for(n=0;n<ND;++n) {
								dcrd(n)(0)(i)(j)(k) = 0.5*(vrtxbd(1)(v(3))(n) -vrtxbd(1)(v(2))(n));
								dcrd(n)(1)(i)(j)(k) = 0.5*(vrtxbd(1)(v(1))(n) -vrtxbd(1)(v(2))(n));
								dcrd(n)(2)(i)(j)(k) = 0.5*(vrtxbd(1)(v(0))(n) -vrtxbd(1)(v(2))(n));
							}
						}
					}
				}
			}
		
			ugtouht(tind,1);
			for(n=0;n<NV;++n)
				basis::tet(log2p).proj(&uht(n)(0),&u(n)(0)(0)(0),stridex, stridey);

			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpy;++j) {
					for(k=0;k<lgpz;++k) {
						cjcb(i)(j)(k) = -gbl->bd(0)*(dcrd(0)(0)(i)(j)(k)*(dcrd(1)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(1)(i)(j)(k))-dcrd(0)(1)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(0)(i)(j)(k))+dcrd(0)(2)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)-dcrd(1)(1)(i)(j)(k)*dcrd(2)(0)(i)(j)(k)));
						for(n=0;n<NV;++n)
							dugdt(log2p,tind,n)(i)(j)(k) = u(n)(i)(j)(k)*cjcb(i)(j)(k);
						for(n=0;n<ND;++n)
							dxdt(log2p,tind,n)(i)(j)(k) = crd(n)(i)(j)(k);
					}
				}	
			}
		}
		
	}
	log2p = log2pmax;
	
	return;
}

#else
void tet_hp::tadvance() {
	int i,j,n,s,ttind,stage;
	
	/* DO STUFF FOR DEFORMABLE MESH FIRST */    
	//    if (log2p == log2pmax && gbl->substep == 0 && (mmovement == coupled_deformable || mmovement == uncoupled_deformable)) {
	//        r_tet_mesh::tadvance();
	//    }
	
	if (!coarse_level) {
		/* SHIFT DATA */
		for (int lvl=gbl->nhist-2;lvl>=0; --lvl) {            
			ugbd(lvl+1).v(Range(0,npnt-1),Range::all()) = ugbd(lvl).v(Range(0,npnt-1),Range::all());
			if (basis::tet(log2p).em) {
				ugbd(lvl+1).e(Range(0,nseg-1),Range::all(),Range::all()) = ugbd(lvl).e(Range(0,nseg-1),Range::all(),Range::all());
				if (basis::tet(log2p).fm) {
					ugbd(lvl+1).f(Range(0,ntri-1),Range::all(),Range::all()) = ugbd(lvl).f(Range(0,ntri-1),Range::all(),Range::all());
					if (basis::tet(log2p).im) {
						ugbd(lvl+1).i(Range(0,ntet-1),Range::all(),Range::all()) = ugbd(lvl).i(Range(0,ntet-1),Range::all(),Range::all());
					}
				}
			}
			/* SAME FOR MESH INFORMATION */
			for(i=0;i<npnt;++i)
				for(n=0;n<ND;++n)
					vrtxbd(lvl+1)(i)(n) = vrtxbd(lvl)(i)(n);
		}
		

		
		/* EXTRAPOLATE */
		if (gbl->dti > 0.0 && gbl->tstep > 1) {
			ugbd(0).v(Range(0,npnt-1),Range::all()) += ugbd(1).v(Range(0,npnt-1),Range::all())-ugbd(2).v(Range(0,npnt-1),Range::all());
			if (basis::tet(log2p).em) {
				ugbd(0).e(Range(0,nseg-1),Range::all(),Range::all()) += ugbd(1).e(Range(0,nseg-1),Range::all(),Range::all())-ugbd(2).e(Range(0,nseg-1),Range::all(),Range::all());
				if (basis::tet(log2p).fm) {
					ugbd(0).f(Range(0,ntri-1),Range::all(),Range::all()) += ugbd(1).f(Range(0,ntri-1),Range::all(),Range::all())-ugbd(2).f(Range(0,ntri-1),Range::all(),Range::all());
					if (basis::tet(log2p).im) {
						ugbd(0).i(Range(0,ntet-1),Range::all(),Range::all()) += ugbd(1).i(Range(0,ntet-1),Range::all(),Range::all())-ugbd(2).i(Range(0,ntet-1),Range::all(),Range::all());
					}
				}
			}
			if (((mmovement == coupled_deformable) || (mmovement == uncoupled_deformable))) 
				vrtxbd(0)(Range(0,npnt-1)) += vrtxbd(1)(Range(0,npnt-1))-vrtxbd(2)(Range(0,npnt-1));               
		}
	}
	else {
		tet_hp* fmesh = dynamic_cast<tet_hp *>(fine);
		
		/* CALCULATE UNSTEADY SOURCE TERMS ON COARSE MESHES */
		for(int lvl=0;lvl<gbl->nhist;++lvl) {
			for(i=0;i<npnt;++i) {
				ttind = fcnnct(i).tet;
				
				ugbd(lvl).v(i,Range::all()) = 0.0;
				
				for(n=0;n<ND;++n)
					vrtxbd(lvl)(i)(n) = 0.0;
				
				for(j=0;j<4;++j) {
					ugbd(lvl).v(i,Range::all()) += fcnnct(i).wt(j)*fmesh->ugbd(lvl).v(fmesh->tet(ttind).pnt(j),Range::all());
					for(n=0;n<ND;++n)
						vrtxbd(lvl)(i)(n) += fcnnct(i).wt(j)*fmesh->vrtxbd(lvl)(fmesh->tet(ttind).pnt(j))(n);
				}
			}
		}
	}
	

	

	
	for(i=0;i<nfbd;++i) 
		hp_fbdry(i)->tadvance();
	
	for(i=0;i<nebd;++i) 
		hp_ebdry(i)->tadvance();
	
	for(i=0;i<nvbd;++i) 
		hp_vbdry(i)->tadvance();
	
	helper->tadvance();
	
	calculate_unsteady_sources();
	
	return;
}

/* A GENERIC CALCULATION OF SOURCES FOR AN AUTONOMOUS SYSTEM IN STANDARD FORM */
/* WILL NEED TO BE OVERRIDDEN FOR SPECIAL CASES */
void tet_hp::calculate_unsteady_sources() {
	int lgpx = basis::tet(log2p).gpx, lgpy = basis::tet(log2p).gpy, lgpz = basis::tet(log2p).gpz;
	int i,j,k,n,tind;
	int stridey = MXGP;
	int stridex = MXGP*MXGP;
	TinyVector<int,4> v;
	
	for (log2p=0;log2p<=log2pmax;++log2p) {
		for (int level=1;level<min(gbl->nhist,gbl->tstep+1);++level) {
			for(tind=0;tind<ntet;++tind) {
				v = tet(tind).pnt;
				
				if (tet(tind).info > -1) {
					crdtocht(tind,level);
					for(n=0;n<ND;++n)
						basis::tet(log2p).proj_bdry(&cht(n)(0), &crd(n)(0)(0)(0), &dcrd(n)(0)(0)(0)(0), &dcrd(n)(1)(0)(0)(0),&dcrd(n)(2)(0)(0)(0),stridex,stridey);
				}
				else {
					for(n=0;n<ND;++n)
						basis::tet(log2p).proj(vrtxbd(level)(v(0))(n),vrtxbd(level)(v(1))(n),vrtxbd(level)(v(2))(n),vrtxbd(level)(v(3))(n),&crd(n)(0)(0)(0),stridex,stridey);
					
					for(i=0;i<lgpx;++i) {
						for(j=0;j<lgpy;++j) {
							for(k=0;k<lgpz;++k) {
								for(n=0;n<ND;++n) {
									dcrd(n)(0)(i)(j)(k) = 0.5*(pnts(tet(tind).pnt(3))(n) -pnts(tet(tind).pnt(2))(n));
									dcrd(n)(1)(i)(j)(k) = 0.5*(pnts(tet(tind).pnt(1))(n) -pnts(tet(tind).pnt(2))(n));
									dcrd(n)(2)(i)(j)(k) = 0.5*(pnts(tet(tind).pnt(0))(n) -pnts(tet(tind).pnt(2))(n));
								}
							}
						}
					}
				}
				
				ugtouht(tind,level);
				for(n=0;n<NV;++n)
					basis::tet(log2p).proj(&uht(n)(0),&u(n)(0)(0)(0),stridex, stridey);
				
				for(i=0;i<lgpx;++i) {
					for(j=0;j<lgpy;++j) {
						for(k=0;k<lgpz;++k) {
							cjcb(i)(j)(k) = -gbl->bd(level)*(dcrd(0)(0)(i)(j)(k)*(dcrd(1)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(1)(i)(j)(k))-dcrd(0)(1)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(0)(i)(j)(k))+dcrd(0)(2)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)-dcrd(1)(1)(i)(j)(k)*dcrd(2)(0)(i)(j)(k)));
							for(n=0;n<NV;++n)
								dugdt(log2p,tind,n)(i)(j)(k) += u(n)(i)(j)(k)*cjcb(i)(j)(k);
							for(n=0;n<ND;++n)
								dxdt(log2p,tind,n)(i)(j)(k) += gbl->bd(level)/gbl->bd(0)*crd(n)(i)(j)(k);
						}
					}	
				}
			}
		}
	}
	log2p = log2pmax;
	
	return;
}

#endif
