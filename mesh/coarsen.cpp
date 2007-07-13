#include "mesh.h"
#include <utilities.h>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <assert.h>
#include <typeinfo>


int mesh::coarsen(FLT factor, const class mesh& inmesh) {
    int i,j,n,sind;
    int v0, v1, odd;
    FLT mindist;
    
    if (!initialized) {
        /* VERTEX STORAGE ALLOCATION */
        allocate_duplicate(1.9,inmesh);
    }
    
    for(i=0;i<maxvst;++i)
        td(i).info = 0;

    /* PREPARE FOR COARSENING */
    for(i=0;i<inmesh.nvrtx;++i)
        fscr1(i) = 1.0e8;
    
    for(i=0;i<inmesh.nside;++i) {
        v0 = inmesh.sd(i).vrtx(0);
        v1 = inmesh.sd(i).vrtx(1);
        fscr1(v0) = MIN(inmesh.distance(v0,v1),fscr1(v0));
        fscr1(v1) = MIN(inmesh.distance(v0,v1),fscr1(v1));
    }

    for(i=0;i<inmesh.nvrtx;++i)
        fscr1(i) *= factor; /* SHOULD BE BETWEEN 1.5 and 2.0 */
        
    nvrtx = 0;
    nside = 0;
    ntri  = 0;

    /* USE I1WK TO KEEP TRACK OF INDICES */
    
    int inmeshbdrysides = 2*inmesh.nside -3*inmesh.ntri;
    if (maxvst < inmeshbdrysides/2) {
        *sim::log << "coarse mesh is not big enough: " << inmeshbdrysides/2 << ' ' << maxvst << std::endl;
    }
    
    /* COARSEN SIDES    */
    for(i=0;i<nsbd;++i) {
        sbdry(i)->nel = 0;
        if (typeid(sbdry(i)) != typeid(inmesh.sbdry(i))) {
            *sim::log << "can't coarsen into object with different boundaries" << std::endl;
            exit(1);
        }

        /* CHECK IF FIRST POINT INSERTED*/
        v0 = inmesh.sd(inmesh.sbdry(i)->el(0)).vrtx(0);
        if (i1wk(v0) < 0) {
            for(n=0;n<ND;++n)
                vrtx(nvrtx)(n) = inmesh.vrtx(v0)(n);
            i1wk(v0) = nvrtx;
            sd(nside).vrtx(0) = nvrtx;
            ++nvrtx;
        }
        else { 
            sd(nside).vrtx(0) = i1wk(v0);
        }

        odd = inmesh.sbdry(i)->nel%2;
        if (odd) {
            for(j=2;j<inmesh.sbdry(i)->nel/2;j+=2) {
                v0 = inmesh.sd(inmesh.sbdry(i)->el(j)).vrtx(0);
                for(n=0;n<ND;++n)
                    vrtx(nvrtx)(n) = inmesh.vrtx(v0)(n);
                i1wk(v0) = nvrtx;
                sd(nside).vrtx(1) = nvrtx;
                sbdry(i)->el(sbdry(i)->nel) = nside;
                sd(nside).tri(1) = -1;
                ++nside; 
                ++sbdry(i)->nel;
                sd(nside).vrtx(0) = nvrtx;
                ++nvrtx;
            } 

            /* MIDDLE POINT OF ODD NUMBERED SIDE */
            if (inmesh.sbdry(i)->nel > 1) {
                j = inmesh.sbdry(i)->nel/2;
                v0 = inmesh.sd(inmesh.sbdry(i)->el(j)).vrtx(0);
                v1 = inmesh.sd(inmesh.sbdry(i)->el(j)).vrtx(1);
                for(n=0;n<ND;++n)
                    vrtx(nvrtx)(n) = 0.5*(inmesh.vrtx(v0)(n) +inmesh.vrtx(v1)(n));
                i1wk(v0) = nvrtx;
                i1wk(v1)= nvrtx;
                sd(nside).vrtx(1) = nvrtx;
                sbdry(i)->el(sbdry(i)->nel) = nside;
                sd(nside).tri(1) = -1;
                ++nside; 
                ++sbdry(i)->nel;
                sd(nside).vrtx(0) = nvrtx;
                ++nvrtx;
            }
            
            for(j = inmesh.sbdry(i)->nel -((inmesh.sbdry(i)->nel-2)/4)*2;j<inmesh.sbdry(i)->nel;j+=2) {
                v0 = inmesh.sd(inmesh.sbdry(i)->el(j)).vrtx(0);
                for(n=0;n<ND;++n)
                    vrtx(nvrtx)(n) = inmesh.vrtx(v0)(n);
                i1wk(v0) = nvrtx;
                sd(nside).vrtx(1) = nvrtx;
                sbdry(i)->el(sbdry(i)->nel) = nside;
                sd(nside).tri(1) = -1;
                ++nside; 
                ++sbdry(i)->nel;
                sd(nside).vrtx(0) = nvrtx;
                ++nvrtx;
            }
        }
        else {
            for(j=2;j<inmesh.sbdry(i)->nel;j+=2) {
                v0 = inmesh.sd(inmesh.sbdry(i)->el(j)).vrtx(0);
                for(n=0;n<ND;++n)
                    vrtx(nvrtx)(n) = inmesh.vrtx(v0)(n);
                i1wk(v0) = nvrtx;
                sd(nside).vrtx(1) = nvrtx;
                sbdry(i)->el(sbdry(i)->nel) = nside;
                sd(nside).tri(1) = -1;
                ++nside; 
                ++sbdry(i)->nel;
                sd(nside).vrtx(0) = nvrtx;
                ++nvrtx;
            }
        }
        
        /* INSERT LAST POINT */
        v0 = inmesh.sd(inmesh.sbdry(i)->el(inmesh.sbdry(i)->nel-1)).vrtx(1);
        if (i1wk(v0) < 0) {
            for(n=0;n<ND;++n)
                vrtx(nvrtx)(n) = inmesh.vrtx(v0)(n);
            i1wk(v0) = nvrtx;
            sd(nside).vrtx(1) = nvrtx;
            sbdry(i)->el(sbdry(i)->nel) = nside;
            sd(nside).tri(1) = -1;
            ++nside;
            ++sbdry(i)->nel;
            ++nvrtx;
        }
        else {
            sd(nside).vrtx(1) = i1wk(v0);
            sbdry(i)->el(sbdry(i)->nel) = nside;
            sd(nside).tri(1) = -1;
            ++nside;
            ++sbdry(i)->nel;
        }
    }
    
    /* MOVE VERTEX BDRY INFORMATION */
    for(i=0;i<inmesh.nvbd;++i) {
        vbdry(i)->copy(*inmesh.vbdry(i));
        vbdry(i)->v0 = i1wk(inmesh.vbdry(i)->v0);
    }
    
    if (maxvst < nside/2+nside) {
        *sim::log << "coarse mesh is not large enough: " << nside/2 +nside << ' ' << maxvst << std::endl;
    }

    
    treeinit();
    
    for(i=0;i<nside;++i)
        i2wk_lst1(i) = i+1;
            
    triangulate(nside);
    
/* FUNNY WAY OF MARKING BECAUSE OF -1 AS UNMARKED */
/* VSPEC IS 0x4 */
#if ((-1)&0x4)
#define ISSPEC(A) (!((A)&VSPEC))
#define SETSPEC ((-1)&(~VSPEC))
#else
#define ISSPEC(A) (((A)&VSPEC))
#define SETSPEC ((-1)|(VSPEC))
#endif

    

    /****************************************************/            
    /* Boyer-Watson Algorithm to insert interior points */
    /****************************************************/
    /* MARK BOUNDARY SO DON'T GET INSERTED */
    /* FUNNY WAY OF MARKING SO CAN LEAVE i1wk initialized to -1 */
    for(i=0;i<inmesh.nsbd;++i) {
        for(j=0;j<inmesh.sbdry(i)->nel;++j) {
            sind = inmesh.sbdry(i)->el(j);
            i1wk(inmesh.sd(sind).vrtx(0)) = SETSPEC;
            i1wk(inmesh.sd(sind).vrtx(1)) = SETSPEC;
        }
    }
    
    /* maxsrch must be high for triangulated domain with no interior points */
    maxsrch = MIN(maxvst,inmesh.ntri);
    
    for(i=0;i<inmesh.nvrtx;++i) {
        if (ISSPEC(i1wk(i))) continue;
        
        mindist = qtree.nearpt(inmesh.vrtx(i).data(),j);
        if (sqrt(mindist) < fscr1(i)) continue;
                
        insert(inmesh.vrtx(i));
    }
    cnt_nbor();
    
    /* reset maxsrch */
    maxsrch = 100;
    
    /* RESET i1wk */
    for(i=0;i<inmesh.nsbd;++i) {
        for(j=0;j<inmesh.sbdry(i)->nel;++j) {
            sind = inmesh.sbdry(i)->el(j);
            i1wk(inmesh.sd(sind).vrtx(0)) = -1;
            i1wk(inmesh.sd(sind).vrtx(1)) = -1;
        }
    }
    bdrylabel();
    initvlngth();
    
    /* PRINT SOME GENERAL DEBUGGING INFO */
    *sim::log << "#" << std::endl << "#COARSE MESH " << std::endl;
    *sim::log << "#MAXVST:" << maxvst << " VERTICES:" << nvrtx << " SIDES:" << nside << " ELEMENTS:" << ntri << std::endl;    
    /* PRINT BOUNDARY INFO */
    for(i=0;i<nsbd;++i)
        *sim::log << "#" << sbdry(i)->idprefix << " " << sbdry(i)->mytype << " " << sbdry(i)->nel << std::endl;

    return(1);
}


block::ctrl mesh::coarsen2(block::ctrl ctrl_message, FLT factor, const class mesh &inmesh, FLT size_reduce) {
    int i;
    
    if (ctrl_message == block::begin) excpt = 0;
    else ++excpt;
        
    switch (excpt) {
        case(0):     
            if (!initialized) {
                allocate_duplicate(size_reduce,inmesh);
            }        
            copy(inmesh);
            initvlngth();
            for(i=0;i<nvrtx;++i)
                vlngth(i) = factor*vlngth(i);      
            setup_for_adapt();
            return(block::advance);

        case(1):            
            bdry_yaber(1.414);
            return(block::advance);

        case(2):
            for(i=0;i<nsbd;++i) 
                sbdry(i)->comm_exchange(boundary::all,0,boundary::master_slave);
            return(block::advance);
            
        case(3):
            /* COARSEN MATCHING BOUNDARIES */
            bdry_yaber1();
            return(block::advance);
            
        case(4):
            /* INTERIOR SIDE COARSENING */
            yaber(1.414);
            return(block::advance);
            
        case(5):
            cleanup_after_adapt();
            qtree.reinit();  // REMOVES UNUSED QUADS
            return(block::stop);
    }
 
    return(block::stop);
}

void mesh::coarsen3() {
    int i,j,v0,tind,sind,node,cnt=0;
    
    /* SET-UP ADAPTION TRACKING STUFF */
    /* For vertices 0 = untouched, 1 = touched, 2 = deleted, 3 = special */
    /* For sides 0 = untouced, 1 = touched, 2 = deleted */
    /* For triangles 0 = untouched, 1 = touched, 2 = deleted, 3 = searched */
    for(i=0;i<maxvst;++i)
        td(i).info = 0;
    
    /* COARSEN SIDE EDGES FIRST */
    for(i=0;i<nsbd;++i) {
        *sim::log << "coarsening boundary " << i << ": type " << sbdry(i)->mytype << " sides " << sbdry(i)->nel << std::endl;
        for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(0);
            if (vd(v0).info > 0) {
                /* COARSEN SIDE */
                collapse(sind,0);
                vd(v0).info = -1;
                ++cnt;
            }
        }
    }
    
    
    for(i=0;i<nvrtx;++i) {
        if (vd(i).info > 0) {
            tind = vd(i).tri;
            for(j=0;j<3;++j) {
                if (td(tind).vrtx(j) == i) {
                    break;
                }
            }
            j = (j+1)%3;
            sind = td(tind).side(j);
            node = (1+td(tind).sign(j))/2;
            collapse(sind,node);
            vd(v0).info = -1;
            ++cnt;
        }
    }
    
    cleanup_after_adapt();
    
    *sim::log << "#Coarsen finished: " << cnt << " sides coarsened" << std::endl;
    for(i=0;i<nsbd;++i)
        *sim::log << "boundary " << i << ": type " << sbdry(i)->mytype << " sides " << sbdry(i)->nel << std::endl;
    
    return;
}
    
    

