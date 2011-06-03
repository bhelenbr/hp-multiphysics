#include "tet_mesh.h"
#include <utilities.h>
#include<assert.h>

void tet_mesh::copy(const tet_mesh& tgt) {
	int i,n;
		
	if (!initialized) {
		init(tgt);
	}
	else {
		/* CHECK IF BIG ENOUGH */
		if (tgt.ntri > maxvst) {
			*gbl->log << "mesh is too big to copy: " << tgt.ntri << ' ' << maxvst << std::endl;
			exit(1);
		}
	}
	
	/* COPY VERTEX INFO OVER */
	npnt = tgt.npnt;
	for(i=0;i<npnt;++i) {
		for(n=0;n<ND;++n)
			pnts(i)(n) = tgt.pnts(i)(n);
		lngth(i) = tgt.lngth(i);
		pnt(i) = tgt.pnt(i);
	}
		
	/* COPY VERTEX BOUNDARY INFO */
	for(i=0;i<nvbd;++i)
		vbdry(i)->copy(*tgt.vbdry(i));
			
	/* COPY SIDE INFORMATION */
	nseg = tgt.nseg;
	for(i=0;i<nseg;++i)
		seg(i) = tgt.seg(i);
		
	/* COPY SIDE BOUNDARY INFO */
	for(i=0;i<nebd;++i)
		ebdry(i)->copy(*tgt.ebdry(i));
	
	/* COPY TRI DATA */
	ntri = tgt.ntri;
	for(i=0;i<ntri;++i)
		tri(i) = tgt.tri(i);
		
	/* COPY FACE BOUNDARY INFO */
	for(i=0;i<nfbd;++i)
		fbdry(i)->copy(*tgt.fbdry(i));
		
	/* COPY ELEMENT DATA */
	ntet = tgt.ntet;
	for(i=0;i<ntet;++i)
		tet(i) = tgt.tet(i);
		
	otree.copy(tgt.otree);
	otree.change_vptr(pnts);
	
	return;  
}

//void tri_mesh::append(const tri_mesh &z) {
//	int i,j,k,n,vrt,sind,flip;
//    int nvrtxold, nsideold,ntriold;
//    int sind1,tind1,v1a,v1b;
//    int sind2,tind2,v2a,v2b;
//            
//	for(i=0;i<z.npnt;++i) {
//		for(n=0;n<ND;++n) 
//			pnts(i+npnt)(n) = z.pnts(i)(n);
//        qtree.addpt(i+npnt);
//    }
//        
//    for(i=0;i<z.npnt;++i)
//        lngth(i+npnt) = z.lngth(i);
//
//    for(i=0;i<z.npnt;++i) {
//        pnt(i+npnt).nnbor = z.pnt(i).nnbor;
//        pnt(i+npnt).tri = z.pnt(i).tri +ntri;
//    }
//        
//    /* MOVE BOUNDARY INFO */
//    vbdry.resizeAndPreserve(nvbd+z.nvbd);
//    for(i=0;i<z.nvbd;++i) {
//        vbdry(nvbd) = z.vbdry(i)->create(*this);
//        vbdry(nvbd)->alloc(4);
//        vbdry(nvbd)->pnt = z.vbdry(i)->pnt +npnt;
//        ++nvbd;
//    }
//	
//	for(i=0;i<z.nseg;++i) {
//		for(n=0;n<2;++n) {
//			seg(i+nseg).pnt(n) = z.seg(i).pnt(n) +npnt;
//            seg(i+nseg).tri(n) = z.seg(i).tri(n) +ntri;
//        }
//    }
//    
//    /* MOVE BOUNDARY INFO */
//    ebdry.resizeAndPreserve(nebd+z.nebd);
//    for(i=0;i<z.nebd;++i) {
//		ebdry(nebd++) = z.ebdry(i)->create(*this);
//        ebdry(nebd-1)->alloc(z.ebdry(i)->maxseg);
//        ebdry(nebd-1)->nseg = z.ebdry(i)->nseg;
//		for(j=0;j<z.ebdry(i)->nseg;++j)
//			ebdry(nebd-1)->seg(j) = z.ebdry(i)->seg(j) +nseg;
//	}
//
//    /* MOVE TRI INFO */
//    for(i=0;i<z.ntri;++i) {
//		for(n=0;n<3;++n) {
//			tri(i+ntri).pnt(n) = z.tri(i).pnt(n) +npnt;
//            tri(i+ntri).seg(n) = z.tri(i).seg(n)+nseg;
//            tri(i+ntri).sgn(n) = z.tri(i).sgn(n);
//            tri(i+ntri).tri(n) = z.tri(i).tri(n)+ntri;
//        }
//    }
//	
//    nvrtxold = npnt;
//    nsideold = nseg;
//    ntriold = ntri;
//    
//	npnt += z.npnt;
//	nseg += z.nseg;
//	ntri += z.ntri;		
//    
//    bdrylabel();  // CHANGES STRI / TTRI ON BOUNDARIES TO POINT TO GROUP/ELEMENT    
//    
//    /* KEEPS TRACK OF DELETED SIDES = -3, TOUCHED SIDES =-2, UNTOUCHED SIDES =-1 */
//    for(i=nsideold;i<nseg;++i)
//        seg(i).info = -1;
//
//    /* VINFO TO KEEP TRACK OF DELETED POINTS (-3) */
//    for(i=nvrtxold;i<npnt;++i)
//        pnt(i).info = -1;
//    
//    /* FIND MATCHING COMMUNICATION BOUNDARIES */
//    for(i=0;i<nebd -z.nebd;++i) {
//        if (!ebdry(i)->is_comm()) continue;
//        
//        for(j=0;j<z.nebd;++j) {
//            if (ebdry(i)->idnum == z.ebdry(j)->idnum) {
//                nseg = ebdry(i)->nseg;
//                
//                /* FIRST POINT DONE REVERSE */
//                sind1 = ebdry(i)->seg(nseg-1);
//                tind1 = seg(sind1).tri(0);                
//                v1b = seg(sind1).pnt(1);
//                sind2 = z.ebdry(j)->seg(0);
//                tind2 = z.seg(sind2).tri(0) +ntriold;
//                v2a = z.seg(sind2).pnt(0) +nvrtxold;
//                pnt(v1b).nnbor += z.pnt(v2a-nvrtxold).nnbor-1;
//                qtree.dltpt(v2a);
//                pnt(v2a).info = -3;
//                do {
//                    for(vrt=0;vrt<3;++vrt) 
//                        if (tri(tind2).pnt(vrt) == v2a) break; 
//                    assert(vrt != 3);
//                    
//                    tri(tind2).pnt(vrt) = v1b;
//                    sind = tri(tind2).seg((vrt +1)%3);
//                    flip = (1 +tri(tind2).sgn((vrt +1)%3))/2;
//                    assert(seg(sind).pnt(flip) == v2a);
//                    seg(sind).pnt(flip) = v1b;
//                    tind2 = tri(tind2).tri((vrt +1)%3);
//                    
//                } while(tind2 > 0); 
//                
//                for(k=0;k<nseg;++k) {
//                    sind1 = ebdry(i)->seg(nseg-k-1);
//                    tind1 = seg(sind1).tri(0);
//                    v1a = seg(sind1).pnt(0);
//                    sind2 = z.ebdry(j)->seg(k);
//                    tind2 = z.seg(sind2).tri(0) +ntriold;
//                    v2b = z.seg(sind2).pnt(1) +nvrtxold;
//                    sind2 += nsideold;                    
//                    pnt(v1a).nnbor += z.pnt(v2b-nvrtxold).nnbor-2;
//
//                    for(vrt=0;vrt<3;++vrt) 
//                        if (tri(tind1).seg(vrt) == sind1) break;
//                    assert(vrt < 3);
//                    tri(tind1).tri(vrt) = tind2;
//                    seg(sind1).tri(1) = tind2;
//                    
//                    for(vrt=0;vrt<3;++vrt) 
//                        if (tri(tind2).pnt(vrt) == v2b) break;
//                    assert(vrt < 3);
//                    tri(tind2).seg((vrt+1)%3) = sind1;
//                    tri(tind2).sgn((vrt+1)%3) = -1;
//                    tri(tind2).tri((vrt+1)%3) = tind1;
//                    seg(sind2).info = -3;
//                    pnt(v2b).info = -3;
//                    qtree.dltpt(v2b);
//                    
//                    for(;;) {
//                        tri(tind2).pnt(vrt) = v1a;
//                        sind = tri(tind2).seg((vrt +2)%3);
//                        flip = (1 -tri(tind2).sgn((vrt +2)%3))/2;
//                        assert(seg(sind).pnt(flip) == v2b);
//                        seg(sind).pnt(flip) = v1a;
//                        tind2 = tri(tind2).tri((vrt +2)%3);
//                        
//                        if (tind2 < 0) break;
//                        
//                        for(vrt=0;vrt<3;++vrt) 
//                            if (tri(tind2).pnt(vrt) == v2b) break; 
//                        assert(vrt != 3);
//                    }
//                }
//                /* LAST POINT ONLY SHARES ONE EDGE */
//                pnt(v1a).nnbor += 1;
//                
//                delete ebdry(i);
//                for(k=i;k<nebd-1;++k)
//                    ebdry(k) = ebdry(k+1);
//                delete ebdry(nebd-z.nebd+j-1);
//                for(k=nebd-z.nebd+j-1;k<nebd-2;++k)
//                    ebdry(k) = ebdry(k+1);
//                nebd -= 2;
//                break;
//            }
//        }
//    }
//    
//        /* DELETE LEFTOVER POINTS */
//    /* VINFO > NPOINT STORES POINT MOVEMENT HISTORY */
//    for(i=0;i<npnt;++i) 
//        if (pnt(i).info == -3) 
//            dltpnt(i);
//    
//    /* FIX BOUNDARY CONDITION POINTERS */
//    for(i=nvbd-z.nvbd;i<nvbd;++i)
//        if (vbdry(i)->pnt >= npnt) 
//            vbdry(i)->pnt = pnt(vbdry(i)->pnt).info;  
//                                
//    /* CLEAN UP SIDES */
//    /* SINFO WILL END UP STORING -1 UNTOUCHED, -2 TOUCHED, or INITIAL INDEX OF UNTOUCHED SIDE */
//    /* SINFO > NSIDE WILL STORE MOVEMENT HISTORY */
//    for(i=nsideold;i<nseg;++i) 
//        if (seg(i).info == -3) dltseg(i);
//    
//    /* FIX BOUNDARY CONDITION POINTERS */
//    for(i=nebd-z.nebd+1;i<nebd;++i)
//        for(j=0;j<ebdry(i)->nseg;++j) 
//            if (ebdry(i)->seg(j) >= nseg) 
//                ebdry(i)->seg(j) = seg(ebdry(i)->seg(j)).info; 
//
//    for (i=0;i<nebd;++i) {
//        ebdry(i)->reorder();
//    }
//    
//  	bdrylabel();  // CHANGES STRI / TTRI ON BOUNDARIES TO POINT TO GROUP/ELEMENT                    
//        
//    return;
//}
//
//
