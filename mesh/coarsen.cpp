#include "mesh.h"
#include <utilities.h>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <assert.h>
#include <typeinfo>


void mesh::coarsen(FLT factor, const class mesh& inmesh) {
    int i,j,n,sind;
    int v0, v1, odd;
    FLT mindist;
    
    if (!initialized) {
        /* VERTEX STORAGE ALLOCATION */
        init(inmesh,1.9);
    }
    
    for(i=0;i<maxvst;++i)
        td(i).info = 0;

    /* PREPARE FOR COARSENING */
    for(i=0;i<inmesh.nvrtx;++i)
        gbl->fltwk(i) = 1.0e8;
    
    for(i=0;i<inmesh.nside;++i) {
        v0 = inmesh.sd(i).vrtx(0);
        v1 = inmesh.sd(i).vrtx(1);
        gbl->fltwk(v0) = MIN(inmesh.distance(v0,v1),gbl->fltwk(v0));
        gbl->fltwk(v1) = MIN(inmesh.distance(v0,v1),gbl->fltwk(v1));
    }

    for(i=0;i<inmesh.nvrtx;++i)
        gbl->fltwk(i) *= factor; /* SHOULD BE BETWEEN 1.5 and 2.0 */
        
    nvrtx = 0;
    nside = 0;
    ntri  = 0;

    /* USE gbl->intwk TO KEEP TRACK OF INDICES */
    
    int inmeshbdrysides = 2*inmesh.nside -3*inmesh.ntri;
    if (maxvst < inmeshbdrysides/2) {
        *gbl->log << "coarse mesh is not big enough: " << inmeshbdrysides/2 << ' ' << maxvst << std::endl;
    }
    
    /* COARSEN SIDES    */
    for(i=0;i<nsbd;++i) {
        sbdry(i)->nel = 0;
        if (typeid(sbdry(i)) != typeid(inmesh.sbdry(i))) {
            *gbl->log << "can't coarsen into object with different boundaries" << std::endl;
            exit(1);
        }

        /* CHECK IF FIRST POINT INSERTED*/
        v0 = inmesh.sd(inmesh.sbdry(i)->el(0)).vrtx(0);
        if (gbl->intwk(v0) < 0) {
            for(n=0;n<ND;++n)
                vrtx(nvrtx)(n) = inmesh.vrtx(v0)(n);
            gbl->intwk(v0) = nvrtx;
            sd(nside).vrtx(0) = nvrtx;
            ++nvrtx;
        }
        else { 
            sd(nside).vrtx(0) = gbl->intwk(v0);
        }

        odd = inmesh.sbdry(i)->nel%2;
        if (odd) {
            for(j=2;j<inmesh.sbdry(i)->nel/2;j+=2) {
                v0 = inmesh.sd(inmesh.sbdry(i)->el(j)).vrtx(0);
                for(n=0;n<ND;++n)
                    vrtx(nvrtx)(n) = inmesh.vrtx(v0)(n);
                gbl->intwk(v0) = nvrtx;
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
                gbl->intwk(v0) = nvrtx;
                gbl->intwk(v1)= nvrtx;
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
                gbl->intwk(v0) = nvrtx;
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
                gbl->intwk(v0) = nvrtx;
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
        if (gbl->intwk(v0) < 0) {
            for(n=0;n<ND;++n)
                vrtx(nvrtx)(n) = inmesh.vrtx(v0)(n);
            gbl->intwk(v0) = nvrtx;
            sd(nside).vrtx(1) = nvrtx;
            sbdry(i)->el(sbdry(i)->nel) = nside;
            sd(nside).tri(1) = -1;
            ++nside;
            ++sbdry(i)->nel;
            ++nvrtx;
        }
        else {
            sd(nside).vrtx(1) = gbl->intwk(v0);
            sbdry(i)->el(sbdry(i)->nel) = nside;
            sd(nside).tri(1) = -1;
            ++nside;
            ++sbdry(i)->nel;
        }
    }
    
    /* MOVE VERTEX BDRY INFORMATION */
    for(i=0;i<inmesh.nvbd;++i) {
        vbdry(i)->copy(*inmesh.vbdry(i));
        vbdry(i)->v0 = gbl->intwk(inmesh.vbdry(i)->v0);
    }
    
    if (maxvst < nside/2+nside) {
        *gbl->log << "coarse mesh is not large enough: " << nside/2 +nside << ' ' << maxvst << std::endl;
    }

    
    treeinit();
    
    for(i=0;i<nside;++i)
        gbl->i2wk_lst1(i) = i+1;
            
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
    /* FUNNY WAY OF MARKING SO CAN LEAVE gbl->intwk initialized to -1 */
    for(i=0;i<inmesh.nsbd;++i) {
        for(j=0;j<inmesh.sbdry(i)->nel;++j) {
            sind = inmesh.sbdry(i)->el(j);
            gbl->intwk(inmesh.sd(sind).vrtx(0)) = SETSPEC;
            gbl->intwk(inmesh.sd(sind).vrtx(1)) = SETSPEC;
        }
    }
    
    /* maxsrch must be high for triangulated domain with no interior points */
    gbl->maxsrch = MIN(maxvst,inmesh.ntri);
    
    for(i=0;i<inmesh.nvrtx;++i) {
        if (ISSPEC(gbl->intwk(i))) continue;
        
        mindist = qtree.nearpt(inmesh.vrtx(i).data(),j);
        if (sqrt(mindist) < gbl->fltwk(i)) continue;
                
        insert(inmesh.vrtx(i));
    }
    cnt_nbor();
    
    /* reset maxsrch */
    gbl->maxsrch = 100;
    
    /* RESET gbl->intwk */
    for(i=0;i<inmesh.nsbd;++i) {
        for(j=0;j<inmesh.sbdry(i)->nel;++j) {
            sind = inmesh.sbdry(i)->el(j);
            gbl->intwk(inmesh.sd(sind).vrtx(0)) = -1;
            gbl->intwk(inmesh.sd(sind).vrtx(1)) = -1;
        }
    }
    bdrylabel();
    initvlngth();
    
    /* PRINT SOME GENERAL DEBUGGING INFO */
    *gbl->log << "#" << std::endl << "#COARSE MESH " << std::endl;
    *gbl->log << "#MAXVST:" << maxvst << " VERTICES:" << nvrtx << " SIDES:" << nside << " ELEMENTS:" << ntri << std::endl;    
    /* PRINT BOUNDARY INFO */
    for(i=0;i<nsbd;++i)
        *gbl->log << "#" << sbdry(i)->idprefix << " " << sbdry(i)->mytype << " " << sbdry(i)->nel << std::endl;

    return;
}


void mesh::coarsen2(FLT factor, const class mesh &inmesh, FLT size_reduce) {
    int i;
  
    if (!initialized) {
        init(inmesh,size_reduce);
    }        
    copy(inmesh);
    initvlngth();
    for(i=0;i<nvrtx;++i)
        vlngth(i) = factor*vlngth(i);      
    setup_for_adapt();
         
    bdry_yaber(1.414);

    for(i=0;i<nsbd;++i) 
        sbdry(i)->comm_exchange(boundary::all,0,boundary::master_slave);

    /* COARSEN MATCHING BOUNDARIES */
    bdry_yaber1();

    /* INTERIOR SIDE COARSENING */
    yaber(1.414);

    cleanup_after_adapt();
    qtree.reinit();  // REMOVES UNUSED QUADS
    
    return;
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
        *gbl->log << "coarsening boundary " << i << ": type " << sbdry(i)->mytype << " sides " << sbdry(i)->nel << std::endl;
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
    
    *gbl->log << "#Coarsen finished: " << cnt << " sides coarsened" << std::endl;
    for(i=0;i<nsbd;++i)
        *gbl->log << "boundary " << i << ": type " << sbdry(i)->mytype << " sides " << sbdry(i)->nel << std::endl;
    
    return;
}
    
    

