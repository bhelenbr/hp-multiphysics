/*
 *  adapt.cpp
 *  planar++
 *
 *  Created by helenbrk on Tue Oct 23 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp.h"
#include "hp_boundary.h"
#include <myblas.h>

void tri_hp::adapt() {
    treeinit();  // FIXME??
    gbl->pstr->copy(*this);
    tri_mesh::adapt();
    setinfo();
    return;
}

void tri_hp::updatepdata(int v0) {
    int n,tind,step; 
    FLT r,s;      
        
    gbl->pstr->findandmvptincurved(pnts(v0),tind,r,s);
        
    for(step=0;step<gbl->nadapt;++step) {
        gbl->pstr->ugtouht(tind,step);
        basis::tri(log2p).ptprobe(NV,&ugbd(step).v(v0,0),&gbl->pstr->uht(0)(0),MXTM);
    }
    
    if (gbl->pstr->tri(tind).info > -1) {
        for(step=1;step<gbl->nadapt;++step) {
            gbl->pstr->crdtocht(tind,step);
            basis::tri(log2p).ptprobe_bdry(ND,&vrtxbd(step)(v0)(0),&gbl->pstr->cht(0,0),MXTM);
        }
    }
    else {
        for(step=1;step<gbl->nadapt;++step) {
            for(n=0;n<ND;++n) 
                vrtxbd(step)(v0)(n) = gbl->pstr->vrtxbd(step)(gbl->pstr->tri(tind).pnt(0))(n)*(s +1.)/2.
                                            +gbl->pstr->vrtxbd(step)(gbl->pstr->tri(tind).pnt(1))(n)*(-r -s)/2.
                                            +gbl->pstr->vrtxbd(step)(gbl->pstr->tri(tind).pnt(2))(n)*(r +1.)/2.;
        }
    }

    return;
}

void tri_hp::updatepdata_bdry(int bnum, int bel, int endpt) {
    int n,sind,sidloc,v0,step;
    FLT psi;
    
    v0 = seg(ebdry(bnum)->seg(bel)).pnt(endpt);
    gbl->pstr->hp_ebdry(bnum)->findandmovebdrypt(pnts(v0),sidloc,psi);
    sind = gbl->pstr->ebdry(bnum)->seg(sidloc);
    
    for(step=0;step<gbl->nadapt;++step) {
        gbl->pstr->ugtouht1d(sind,step);
        basis::tri(log2p).ptprobe1d(NV,&ugbd(step).v(v0,0),&gbl->pstr->uht(0)(0),MXTM);
    }

    if (hp_ebdry(bnum)->is_curved()) {
        for(step=1;step<gbl->nadapt;++step) {
            gbl->pstr->crdtocht1d(sind,step);
            basis::tri(log2p).ptprobe1d(ND,&vrtxbd(step)(v0)(0),&gbl->pstr->cht(0,0),MXTM);
        }
    }
    else {
        for(step=1;step<gbl->nadapt;++step) {
            for(n=0;n<ND;++n) 
                vrtxbd(step)(v0)(n) = gbl->pstr->vrtxbd(step)(gbl->pstr->seg(sind).pnt(0))(n)*(1. -psi)/2.
                                                 +gbl->pstr->vrtxbd(step)(gbl->pstr->seg(sind).pnt(1))(n)*(1. +psi)/2.;
        }
    }
    
    /* FOR INTERNALLY STORED DATA */
    hp_ebdry(bnum)->updatepdata_bdry(bel,endpt,gbl->pstr->hp_ebdry(bnum));
    
    return;
}

void tri_hp::movepdata(int from, int to) {
    int n,step;
            
    for(step=0;step<gbl->nadapt;++step) {
        for(n=0;n<NV;++n)
            ugbd(step).v(to,n) = ugbd(step).v(from,n);
            
        for(n=0;n<ND;++n)
            vrtxbd(step)(to)(n) = vrtxbd(step)(from)(n);
    }
    
    return;
}

void tri_hp::movepdata_bdry(int bnum,int bel,int endpt) {
    /* This is just for internal data (if any) */
    hp_ebdry(bnum)->movepdata_bdry(bel,endpt,gbl->pstr->hp_ebdry(bnum));
}

static int error_count = 0;

void tri_hp::updatesdata(int sind) {
    int i,m,n,v0,v1,tind,step,info;
    FLT r,s,upt[NV];
    char uplo[] = "U";
    TinyVector<FLT,2> pt;
	int ierr;
    
    if (!sm0) return;
    
    v0 = seg(sind).pnt(0);
    v1 = seg(sind).pnt(1);
    
    for(n=0;n<ND;++n)
        basis::tri(log2p).proj1d(pnts(v0)(n),pnts(v1)(n),&crd(n)(0,0));

    for(step=0;step<gbl->nadapt;++step)
        for(n=0;n<NV;++n)
            basis::tri(log2p).proj1d(ugbd(step).v(v0,n),ugbd(step).v(v1,n),&bdwk(step,n)(0,0));
        
    for(i=0;i<basis::tri(log2p).gpx;++i) {
        pt(0) = crd(0)(0,i);
        pt(1) = crd(1)(0,i);
        ierr = gbl->pstr->findinteriorpt(pt,tind,r,s);
		if (ierr) {
			*gbl->log << "Warning #" << error_count << ": didn't find interior point in updatesdata for " << sind << ' ' << pt << std::endl;
			std::ostringstream fname;
			fname << "current_solution" << error_count++ << '_' << gbl->idprefix;
			tri_mesh::output(fname.str().c_str(),tri_mesh::grid);
			tri_hp::output(fname.str().c_str(),tri_hp::tecplot);
		}
        
        for(step=0;step<gbl->nadapt;++step) {
            gbl->pstr->ugtouht(tind,step);
            basis::tri(log2p).ptprobe(NV,upt,&gbl->pstr->uht(0)(0),MXTM);
            for(n=0;n<NV;++n)    {
                bdwk(step,n)(0,i) -= upt[n];
            }
        }
    }

    for(step=0;step<gbl->nadapt;++step) {
        for(n=0;n<NV;++n)
            basis::tri(log2p).intgrt1d(&lf(n)(0),&bdwk(step,n)(0,0));

        for(n=0;n<NV;++n) {
            PBTRS(uplo,basis::tri(log2p).sm,basis::tri(log2p).sbwth,1,&basis::tri(log2p).sdiag1d(0,0),basis::tri(log2p).sbwth+1,&lf(n)(2),basis::tri(log2p).sm,info);
            for(m=0;m<basis::tri(log2p).sm;++m) 
                ugbd(step).s(sind,m,n) = -lf(n)(2+m);
        }
    }
    return;
}

void tri_hp::updatesdata_bdry(int bnum,int bel) {
    int m,n,sind,v0,v1,step,stgt,info;
    TinyVector<FLT,2> pt;
    FLT psi;
    FLT upt[NV];
    char uplo[] = "U";
    
    if (!sm0) return;
    
    sind = ebdry(bnum)->seg(bel);
    v0 = seg(sind).pnt(0);
    v1 = seg(sind).pnt(1);

    for(step=0;step<gbl->nadapt;++step)
        for(n=0;n<NV;++n)
            basis::tri(log2p).proj1d(ugbd(step).v(v0,n),ugbd(step).v(v1,n),&bdwk(step,n)(0,0));
            
    for(step=0;step<gbl->nadapt;++step)
        for(n=0;n<ND;++n)
            basis::tri(log2p).proj1d(vrtxbd(step)(v0)(n),vrtxbd(step)(v1)(n),&bdwk(step,n)(1,0));
            
    
    if (hp_ebdry(bnum)->is_curved()) {

        for(m=0;m<basis::tri(log2p).gpx;++m) {
            pt(0) = bdwk(0,0)(1,m);
            pt(1) = bdwk(0,1)(1,m);
            gbl->pstr->hp_ebdry(bnum)->findandmovebdrypt(pt,stgt,psi);
            stgt = gbl->pstr->ebdry(bnum)->seg(stgt);
              
            for(step=0;step<gbl->nadapt;++step) {
                gbl->pstr->ugtouht1d(stgt,step);
                basis::tri(log2p).ptprobe1d(NV,upt,&gbl->pstr->uht(0)(0),MXTM);
                for(n=0;n<NV;++n)    
                    bdwk(step,n)(0,m) -= upt[n];
          
                gbl->pstr->crdtocht1d(stgt,step);
                basis::tri(log2p).ptprobe1d(ND,upt,&gbl->pstr->cht(0,0),MXTM);
                for(n=0;n<ND;++n)    
                    bdwk(step,n)(1,m) -= upt[n];
            }                          
        }      

        for(step=0;step<gbl->nadapt;++step) {
            for(n=0;n<ND;++n) {
                basis::tri(log2p).intgrt1d(&lf(n)(0),&bdwk(step,n)(1,0));
                PBTRS(uplo,basis::tri(log2p).sm,basis::tri(log2p).sbwth,1,&basis::tri(log2p).sdiag1d(0,0),basis::tri(log2p).sbwth+1,&lf(n)(2),basis::tri(log2p).sm,info);
            
                for(m=0;m<basis::tri(log2p).sm;++m)
                    hp_ebdry(bnum)->crdsbd(step,bel,m,n) = -lf(n)(m+2);
            }
        }
    }
    else {
        for(m=0;m<basis::tri(log2p).gpx;++m) {
            pt(0) = bdwk(0,0)(1,m);
            pt(1) = bdwk(0,1)(1,m);
            
            /* FIND PSI */                
            gbl->pstr->hp_ebdry(bnum)->findandmovebdrypt(pt,stgt,psi);
            stgt = gbl->pstr->ebdry(bnum)->seg(stgt);
            
            /* CALCULATE VALUE OF SOLUTION AT POINT */
            for(step=0;step<gbl->nadapt;++step) {
                gbl->pstr->ugtouht1d(stgt,step);
                basis::tri(log2p).ptprobe1d(NV,upt,&gbl->pstr->uht(0)(0),MXTM);
                for(n=0;n<NV;++n)    
                    bdwk(step,n)(0,m) -= upt[n];
            }
        }
    }

    for(step=0;step<gbl->nadapt;++step) {
        for(n=0;n<NV;++n)
            basis::tri(log2p).intgrt1d(&lf(n)(0),&bdwk(step,n)(0,0));

        for(n=0;n<NV;++n) {
            PBTRS(uplo,basis::tri(log2p).sm,basis::tri(log2p).sbwth,1,&basis::tri(log2p).sdiag1d(0,0),basis::tri(log2p).sbwth+1,&lf(n)(2),basis::tri(log2p).sm,info);
            for(m=0;m<basis::tri(log2p).sm;++m) {
                ugbd(step).s(sind,m,n) = -lf(n)(2+m);
            }
        }
    }
    
    /* UPDATE INTERNAL INFORMATION */
    hp_ebdry(bnum)->updatesdata_bdry(bel,gbl->pstr->hp_ebdry(bnum));
    
    return;
}

void tri_hp::movesdata(int from, int to) {
    int step;
    
    if (!sm0) return;
    
    for(step=0;step<gbl->nadapt;++step)
        ugbd(step).s(to,Range::all(),Range::all()) = ugbd(step).s(from,Range::all(),Range::all());

    return;
}

void tri_hp::movesdata_bdry(int bnum,int bel) {
    
    hp_ebdry(bnum)->movesdata_bdry(bel,gbl->pstr->hp_ebdry(bnum));

    return;
}

void tri_hp::updatetdata(int tind) {
    int i,j,n,ttgt,step,info;
    FLT r,s;
    FLT upt[NV];
    char uplo[] = "U";
    TinyVector<FLT,2> pt;
	int ierr;
    
    if (!im0) return;  /* FIXME NEED TO FIX THIS IN MESH SO CAN TURN OFF ENTIRE LOOP */
        
    for(step=0;step<gbl->nadapt;++step) {
        ugtouht_bdry(tind,step);
        for(n=0;n<NV;++n)
            basis::tri(log2p).proj_bdry(&uht(n)(0),&bdwk(step,n)(0,0),MXGP);
    }
    
    crdtocht(tind);
    for(n=0;n<ND;++n)
        basis::tri(log2p).proj_bdry(&cht(n,0),&crd(n)(0,0),MXGP);
    
    for (i=0; i < basis::tri(log2p).gpx; ++i ) {
        for (j=0; j < basis::tri(log2p).gpn; ++j ) {
            pt(0) = crd(0)(i,j);
            pt(1) = crd(1)(i,j);
            ierr = gbl->pstr->findinteriorpt(pt,ttgt,r,s);
			if (ierr) {
				*gbl->log << "Warning #" << error_count << ": didn't find interior point in updatetdata for " << tind << ' ' << pt << std::endl;
				std::ostringstream fname;
				fname << "current_solution" << error_count++ << '_' << gbl->idprefix;
				tri_mesh::output(fname.str().c_str(),tri_mesh::grid);
				tri_hp::output(fname.str().c_str(),tri_hp::tecplot);
			}            
            for(step=0;step<gbl->nadapt;++step) {
                gbl->pstr->ugtouht(ttgt,step);
                basis::tri(log2p).ptprobe(NV,upt,&gbl->pstr->uht(0)(0),MXTM);
                for(n=0;n<NV;++n)
                    bdwk(step,n)(i,j) -= upt[n];
            }
        }
    }
                        
    for(step=0;step<gbl->nadapt;++step) {
        for(n=0;n<NV;++n) {
            basis::tri(log2p).intgrt(lf(n).data(),&bdwk(step,n)(0,0),MXGP);
            PBTRS(uplo,basis::tri(log2p).im,basis::tri(log2p).ibwth,1,&basis::tri(log2p).idiag(0,0),basis::tri(log2p).ibwth+1,&lf(n)(basis::tri(log2p).bm),basis::tri(log2p).im,info);
            for(i=0;i<basis::tri(log2p).im;++i)
                ugbd(step).i(tind,i,n) = -lf(n)(basis::tri(log2p).bm+i);
        }
    }
        
    return;
}

void tri_hp::movetdata(int from, int to) {
    int step;
    
    if (!im0) return;  /* FIXME NEED TO FIX THIS IN MESH SO CAN TURN OFF ENTIRE LOOP */

    for(step=0;step<gbl->nadapt;++step) {
        ugbd(step).i(to,Range::all(),Range::all()) = ugbd(step).i(from,Range::all(),Range::all());
    }
                
    return;
}
        

    

