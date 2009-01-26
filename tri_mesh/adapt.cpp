/*
 *  adapt.cpp
 *  mesh
 *
 *  Created by Brian Helenbrook on 12/2/04.
 *  Copyright 2004 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_mesh.h"


/* Memory Usage */
/* tri(i).info: (Should be able to merge this with intwk, except for reorder?) */
/* For vertices 0 = untouched, 1 = touched, 2 = deleted, 3 = special */
/* For sides 0 = untouced, 1 = touched, 2 = deleted */
/* For triangles 0 = untouched, 1 = touched, 2 = deleted, 3 = searched */

/* List of entities needing refinement */
/* seg(i).info */

/* Back reference into list of entities needing refinement */
/* pnt(i).info */

/* ratio of triangle size to target */
/* fscr1(i) */

/* in insert & bdry_insert */
/* gbl->i2wk_lst1(ntdel++) = tnum; */
/* gbl->i2wk_lst2(nskeep++) = sind; */
/* gbl->i2wk_lst3(nsdel++) = sind; */
		
/* in findtri */
/* intwk mark for searches */
/* i2wk list of triangles needing search */

/* in yaber & collapse */
/* gbl->i2wk_lst1, gbl->i2wk_lst2, gbl->i2wk_lst3 */
		
/* in reorder */
/* intwk + i2wk */

void tri_mesh::adapt() {
    int i;
    
    /* CALCULATE TARGET LENGTH */
    length();
    for(bool last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
        pmsgload(boundary::all_phased,mp_phase,boundary::symmetric,lngth.data(),0,0,1);
        pmsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
        last_phase = true;
        last_phase &= pmsgwait_rcv(boundary::all_phased,mp_phase, boundary::symmetric, boundary::average,lngth.data(),0,0,1);
    }    

    /* SET FLAGS ETC... */
    setup_for_adapt();
    
    /* SWAP EGES TO MAKE DELAUNAY */
    swap(1.0e-10); 
    
    /* COARSEN FIRST EDGES & SEND MESSAGES */
    bdry_yaber(gbl->tolerance);
    
    /* TRANSFER MESSAGES */
    for(i=0;i<nebd;++i) 
        ebdry(i)->comm_exchange(boundary::all,0,boundary::master_slave);
    
    /* COARSEN MATCHING BOUNDARIES */
    bdry_yaber1();
    
    /* COARSEN INTERIOR */
    yaber(gbl->tolerance);
            
    /* REFINE FIRST EDGES */
    bdry_rebay(gbl->tolerance);
            
    /* TRANSFER MESSAGES */
    for(i=0;i<nebd;++i) 
        ebdry(i)->comm_exchange(boundary::all,0,boundary::master_slave);
            
    /* REFINE MATCHING EDGES */
    bdry_rebay1();
            
    /* REFINE INTERIOR */
    rebay(gbl->tolerance);
            
    /* REMOVE DELETED ENTITIES */
    cleanup_after_adapt();
}


void tri_mesh::setup_for_adapt() {
    int i, p0, p1;

    /* SET-UP ADAPTION TRACKING STUFF */
    /* For vertices 0 = untouched, 1 = touched, 2 = deleted, 3 = special */
    /* For sides 0 = untouced, 1 = touched, 2 = deleted */
    /* For triangles 0 = untouched, 1 = touched, 2 = deleted, 3 = searched */
    for(i=0;i<maxpst;++i)
        tri(i).info = 0;
        
    /* Back reference into list of entities needing refinement */
    for(i=0;i<maxpst;++i)
        pnt(i).info = -1;
        
    /* MARK BEGINNING/END OF SIDE GROUPS & SPECIAL VERTEX POINTS */
    /* THESE SHOULD NOT BE DELETED */
    for(i=0;i<nebd;++i) {
        p0 = seg(ebdry(i)->seg(0)).pnt(0);
        p1 = seg(ebdry(i)->seg(ebdry(i)->nseg-1)).pnt(1);
        tri(p0).info |=  PSPEC;
        tri(p1).info |=  PSPEC;
    }
    
    for(i=0;i<nvbd;++i)
        tri(vbdry(i)->pnt).info |= PSPEC;
        
    return;
}
        
        
void tri_mesh::cleanup_after_adapt() {
    int i,j,sind,p0;
    
    if (gbl->adapt_output) {
        std::string adapt_file;
        std::ostringstream nstr;
        nstr << gbl->tstep << std::flush;
        adapt_file = gbl->idprefix +"_adapt" +nstr.str();
        nstr.str("");
        output(adapt_file.c_str(),debug_adapt);
    }
    
    /* DELETE SIDES FROM BOUNDARY CONDITIONS */
    for(i=0;i<nebd;++i) {		
        for(j=ebdry(i)->nseg-1;j>=0;--j) {
            if (tri(ebdry(i)->seg(j)).info&SDLTE) {
                ebdry(i)->seg(j) = ebdry(i)->seg(--ebdry(i)->nseg);
				ebdry(i)->prev(j) = ebdry(i)->prev(ebdry(i)->nseg);
				ebdry(i)->next(j) = ebdry(i)->next(ebdry(i)->nseg);
				ebdry(i)->prev(ebdry(i)->next(j)) = j;
				ebdry(i)->next(ebdry(i)->prev(j)) = j;
			}
        }
        ebdry(i)->reorder();
    }
    bdrylabel();

    
    /* UPDATE BOUNDARY DATA */
	/* AT THIS POINT NO SIDE OR PNTS HAVE BEEN MOVED */
	/* SO CAN USE SIND TO LOOK UP BOUNDARY ELEMENT NUMBER IN OLD MESH FOR UNTOUCHED SIDES */
	/* THIS PROBABLY WON'T WORK FOR VERTICES.  WILL HAVE TO CHANGE A LITTLE */
	/* DON'T HAVE ANY DIRECT WAY OF SAYING HOW BOUNDARY ENTITIES HAVE CHANGE POSITION */
    for (i=0;i<nebd;++i) {
        /* THIS IS PROBABLY UNNECESSARY SINCE FIRST */
        /* & LAST VERTEX SHOULD NEVER BE CHANGED */
        sind = ebdry(i)->seg(0);
        p0 = seg(sind).pnt(0);
        if (tri(p0).info&PTOUC) {
            updatepdata_bdry(i,0,0);
            tri(p0).info &= ~PTOUC;
        }
        else movepdata_bdry(i,0,0);
        
        for(j=0;j<ebdry(i)->nseg;++j) {
            sind = ebdry(i)->seg(j);
            p0 = seg(sind).pnt(1);
            if (tri(p0).info&PTOUC) {
                updatepdata_bdry(i,j,1);
                tri(p0).info &= ~PTOUC;
            }
            else movepdata_bdry(i,j,1);
            if (tri(sind).info&STOUC) {
                updatesdata_bdry(i,j);
                tri(sind).info &= ~STOUC;
            }
            else movesdata_bdry(i,j);
        }        
    }
    
    /* DELETE LEFTOVER POINTS */
    /* VINFO > NPOINT STORES POINT MOVEMENT HISTORY */            
    for(i=0;i<npnt;++i) {
        if (tri(i).info&PDLTE) dltpnt(i);
        else if (tri(i).info&PTOUC) {
            pnt(i).info = -2;
            updatepdata(i);
        }
    }
    
    /* FIX BOUNDARY CONDITION POINTERS */
    for(i=0;i<nvbd;++i)
        if (pnt(vbdry(i)->pnt).info > -1) 
            vbdry(i)->pnt = pnt(vbdry(i)->pnt).info;  
                                
    /* CLEAN UP SIDES */
    /* SINFO WILL END UP STORING -1 UNTOUCHED, -2 TOUCHED, or INITIAL INDEX OF UNTOUCHED SIDE */
    /* SINFO > NSIDE WILL STORE MOVEMENT LOCATION */  /* FIXME: HAVEN"T TESTED THIS */
    for(i=0;i<nseg;++i) {
        seg(i).info = -1;
        if (tri(i).info&SDLTE) dltseg(i);
        else if (tri(i).info&STOUC) {
            seg(i).info = -2;
            updatesdata(i);
        }
    }
            
    /* FIX BOUNDARY CONDITION POINTERS */
    for(i=0;i<nebd;++i)
        for(j=0;j<ebdry(i)->nseg;++j) 
            if (ebdry(i)->seg(j) >= nseg) 
                ebdry(i)->seg(j) = seg(ebdry(i)->seg(j)).info; 
        
    /* CLEAN UP DELETED TRIS */
    /* TINFO < NTRI STORES INDEX OF ORIGINAL TRI ( > 0), TINFO = 0 -> UNMOVED */
    /* TINFO > NTRI STORES TRI MOVEMENT HISTORY */
    for(i=0;i<ntri;++i) {
        if (tri(i).info&TDLTE) dlttri(i);
        else if (tri(i).info&TTOUC) {
            tri(i).info = -2;
            updatetdata(i);
        }
    }
    
    cnt_nbor();
        
    return;
}


void tri_mesh::putinlst(int sind) {
    int i, temp, top, bot, mid;
        
    /* CREATE ORDERED LIST OF SIDES SMALLEST gbl->fltwk TO LARGEST */
    bot = 0;
    if (gbl->nlst > 0) {
        top = 0;
        bot = gbl->nlst-1;
        if (gbl->fltwk(sind) < gbl->fltwk(seg(top).info)) {
            bot = 0;
        }
        else if (gbl->fltwk(sind) > gbl->fltwk(seg(bot).info)) {
            bot = gbl->nlst;
        }
        else {
            while(top < bot-1) {
                mid = top + (bot -top)/2;
                if (gbl->fltwk(sind) > gbl->fltwk(seg(mid).info))
                    top = mid;
                else
                    bot = mid;
            }
        }
        for(i=gbl->nlst-1;i>=bot;--i) {
            temp = seg(i).info;
            seg(i+1).info = temp;
            pnt(temp).info = i+1;
        }
    }
    seg(bot).info= sind;
    pnt(sind).info = bot;
    ++gbl->nlst;

    assert(gbl->nlst < maxpst -1);
    
    return;
}

void tri_mesh::tkoutlst(int sind) {
    int bgn,temp,i;
    
    bgn = pnt(sind).info;
    for(i=bgn+1;i<gbl->nlst;++i) {
        temp = seg(i).info;
        seg(i-1).info = temp;
        pnt(temp).info = i-1;
    }
    pnt(sind).info = -1;
    --gbl->nlst;
    seg(gbl->nlst).info = -1;
    
    return;
}


        

    
    


        

