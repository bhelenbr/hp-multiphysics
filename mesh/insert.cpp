/*
 *  insert.cpp
 *  mblock
 *
 *  Created by helenbrk on Fri Aug 31 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#include "mesh.h"
#include <math.h>
#include <float.h>    
#include <assert.h>
#include <stdlib.h>

int tri_mesh::insert(const TinyVector<FLT,ND> &x) {
    int n,tind,vnear,err;
    
    for(n=0;n<ND;++n)
        vrtx(nvrtx)(n) = x(n);
    
    qtree.addpt(nvrtx);
    qtree.nearpt(nvrtx,vnear);

    /* FIND TRIANGLE CONTAINING POINT */        
    tind = findtri(x,vnear);
    if (tind < 0) {
        std::cerr << "couldn't find triangle for point: " << x(0) << ' ' << x(1) << " vnear: " << vnear << std::endl;
        std::cerr << "maxsrch: " << gbl->maxsrch << "vtri: " << vd(vnear).tri << std::endl;
        output("error");
        exit(1);
    }     
    if (nvrtx >= maxvst) {
        *gbl->log << "need to use larger growth factor: too many vertices" << std::endl;
        output("error");
        exit(1);
    }
    err = insert(nvrtx,tind);
    nvrtx += 1 -err;
    
    return(err);
}

int tri_mesh::insert(int vnum, int tnum) {
    int ntdel, nskeep, nsdel;
    int i,j,tind,tin,tnxt,v0,sstart,dir,rstrt;
    int sind,sind1,snum;
    
    /* SEARCH SURROUNDING FOR NONDELAUNEY TRIANGLES */
    /* SIDES ON BOUNDARY OF HOLE (SKEEP) */
    /* SIDES IN INTERIOR OF HOLE (SDEL) (STARTS AT HALFWAY THROUGH LIST) */
    ntdel = 0;
    nskeep = 0;
    nsdel = 0;

    gbl->i2wk_lst1(ntdel++) = tnum;
    td(tnum).info |= TSRCH;
    
    tind = tnum;
    v0 = td(tnum).vrtx(0);
    /* FIRST SIDE DONE SEPARATELY SO KNOW WHEN TO STOP */
    snum = 2;
    tnxt = td(tind).tri(snum);
    sind = td(tind).side(snum);
    dir = (1+td(tind).sign(snum))/2;
    sstart = sind;
    rstrt = 1;
    
    if (tnxt < 0) {
        gbl->i2wk_lst2(nskeep++) = sind;
        gbl->i2wk_lst2(nskeep++) = dir;
        v0 = sd(sind).vrtx(dir);
        snum = (snum+1)%3;
    }
    else if (incircle(tnxt,vrtx(vnum)) > 0.0) {
        if (dir > 0) gbl->i2wk_lst3(nsdel++) = sind;
        rstrt = 2;
        gbl->i2wk_lst1(ntdel++) = tnxt;
        td(tnxt).info |= TSRCH;
        tind = tnxt;
        for(snum=0;snum<3;++snum) {
            if (td(tind).vrtx(snum) == v0) break;
        }
        snum = (snum+2)%3;
    }
    else {
        gbl->i2wk_lst2(nskeep++) = sind;
        gbl->i2wk_lst2(nskeep++) = dir;
        v0 = sd(sind).vrtx(dir);
        snum = (snum+1)%3;
    }
    tnxt = td(tind).tri(snum);
    sind = td(tind).side(snum);
    dir = (1+td(tind).sign(snum))/2;    
    
    /* GO COUNTER-CLOCKWISE AROUND VERTICES OF HOLE */
    /* IF START SIDE IS IN THE INTERIOR MUST HIT IT TWICE */
    for(j=0;j<rstrt;++j) {
        do  {
            if (tnxt < 0) {
                gbl->i2wk_lst2(nskeep++) = sind;
                gbl->i2wk_lst2(nskeep++) = dir;
                v0 = sd(sind).vrtx(dir);
                snum = (snum+1)%3;
            }
            else if (td(tnxt).info&TSRCH) {
                if (dir > 0) gbl->i2wk_lst3(nsdel++) = sind;
                tind = tnxt;
                for(snum=0;snum<3;++snum) {
                    if (td(tind).vrtx(snum) == v0) break;
                }
                snum = (snum+2)%3;
            }
            else if (incircle(tnxt,vrtx(vnum)) > 0.0) {
                if (dir > 0) gbl->i2wk_lst3(nsdel++) = sind;
                gbl->i2wk_lst1(ntdel++) = tnxt;
                td(tnxt).info |= TSRCH;
                tind = tnxt;
                for(snum=0;snum<3;++snum) {
                    if (td(tind).vrtx(snum) == v0) break;
                }
                snum = (snum+2)%3;
            }
            else {
                gbl->i2wk_lst2(nskeep++) = sind;
                gbl->i2wk_lst2(nskeep++) = dir;
                v0 = sd(sind).vrtx(dir);
                snum = (snum+1)%3;
            }
            tnxt = td(tind).tri(snum);
            sind = td(tind).side(snum);
            dir = (1+td(tind).sign(snum))/2;    
        } while (sind != sstart);
    }
        
    /* RESET TSRCH FLAGS */
    for(i=0;i<ntdel;++i)
        td(gbl->i2wk_lst1(i)).info &= ~TSRCH;
        
    nskeep = nskeep >> 1;
            
    /*	CHECK THAT WE AREN'T INSERTING POINT VERY CLOSE TO BOUNDARY */
//    for(i=0;i<nskeep;++i) {
//        sind = gbl->i2wk_lst2(i);
//        if(fabs(minangle(vnum, sd(sind).vrtx(0) , sd(sind).vrtx(1))) < 6.0*M_PI/180.0) {
//            *gbl->log << "#Warning: inserting close to boundary" << std::endl;
//        }
//    }
    
    /* APPEND NEW INDICES TO LIST OF INDICES TO USE FOR NEW SIDES & TRIS*/
    for(i=nsdel;i<nskeep;++i)
        gbl->i2wk_lst3(i) = nside +(i-nsdel);
    
    for(i=ntdel;i<nskeep;++i) 
        gbl->i2wk_lst1(i) = ntri +(i-ntdel);
    
    /* PERIODIC TRIANGLE */
    gbl->i2wk_lst1(nskeep) = gbl->i2wk_lst1(0);

    ntri += 2;
    nside += 3;
    
    if (ntri > maxvst || nside > maxvst) {
        *gbl->log << "need to use bigger growth factor: too many sides/tris" << nside << ntri << std::endl;
        output("error");
        exit(1);
    }
    
    for(i=0;i<nskeep;++i) {
        tind = gbl->i2wk_lst1(i);
        td(tind).info |= TTOUC;

        tnxt = gbl->i2wk_lst1(i+1);
        sind = gbl->i2wk_lst2(i<<1);
        dir = gbl->i2wk_lst2((i<<1) +1);
            
        /* CREATE NEW INFO */
        v0 = sd(sind).vrtx(1-dir);
        td(tind).vrtx(0) = v0;
        vd(v0).tri = tind;
        v0 = sd(sind).vrtx(dir);
        td(tind).vrtx(1) = v0;
        vd(v0).tri = tind;
        td(tind).vrtx(2) = vnum;
        vd(vnum).tri = tind;

        /* SIDE 2 INFO */        
        td(tind).side(2) = sind;
        td(tind).sign(2) = -1 +2*dir;
        
        sd(sind).tri(1-dir) = tind;
        tin = sd(sind).tri(dir);
        td(tind).tri(2) = tin;
        if (tin > -1) {
            for(j=0;j<3;++j) {
                if (td(tin).side(j) == sind) {
                    td(tin).tri(j) = tind;
                    break;
                }
            }
        }
            
        /* CREATE SIDE 0 */
        v0 = sd(sind).vrtx(dir);
        sind1 = gbl->i2wk_lst3(i);
        sd(sind1).tri(0) = tind;
        sd(sind1).vrtx(0) = v0;
        sd(sind1).vrtx(1) = vnum;
        td(sind1).info |= STOUC;


        td(tind).side(0) = sind1;
        td(tind).sign(0) = 1;
        td(tind).tri(0) = tnxt;

        sd(sind1).tri(1) = tnxt;
        td(tnxt).side(1) = sind1;
        td(tnxt).sign(1) = -1;
        td(tnxt).tri(1) = tind;
    }
    
    gbl->i2wk_lst1(-1) = ntdel;
    gbl->i2wk_lst2(-1) = nskeep;
    gbl->i2wk_lst3(-1) = nsdel;
        
    return(0);
}

void tri_mesh::bdry_insert(int vnum, int sind, int endpt) {
    int ntdel, nsdel, nskeep;
    int sstart;
    int i,j,tin,tind,tnxt,v0,dir;
    int snum, sind1;
    
    /* ADD POINT TO QUADTREE */
    qtree.addpt(nvrtx);
    td(vnum).info |= VTOUC;
    
    /* SEARCH SURROUNDING FOR NONDELAUNEY TRIANGLES */
    /* SIDES ON BOUNDARY OF HOLE (SKEEP) */
    /* SIDES IN INTERIOR OF HOLE (SDEL) (STARTS AT HALFWAY THROUGH LIST) */
    ntdel = 0;
    nskeep = 0;
    nsdel = 0;

    tind = sd(sind).tri(0);
    gbl->i2wk_lst1(ntdel++) = tind;
    td(tind).info |= TSRCH;
    for(snum=0;snum<3;++snum)
        if (td(tind).side(snum) == sind) break;

    sstart = sind;
    v0 = sd(sind).vrtx(1);
    snum = (snum+1)%3;
    tnxt = td(tind).tri(snum);
    sind = td(tind).side(snum);
    dir = (1+td(tind).sign(snum))/2;
    
    /* GO COUNTER-CLOCKWISE AROUND VERTICES OF HOLE */
    do  {
        if (tnxt < 0) {
            gbl->i2wk_lst2(nskeep++) = sind;
            gbl->i2wk_lst2(nskeep++) = dir;
            v0 = sd(sind).vrtx(dir);
            snum = (snum+1)%3;
        }
        else if (td(tnxt).info&TSRCH) {
            if (dir > 0) gbl->i2wk_lst3(nsdel++) = sind;
            tind = tnxt;
            for(snum=0;snum<3;++snum) {
                if (td(tind).vrtx(snum) == v0) break;
            }
            snum = (snum+2)%3;
        }
        else if (incircle(tnxt,vrtx(vnum)) > 0.0) {
            if (dir > 0) gbl->i2wk_lst3(nsdel++) = sind;
            gbl->i2wk_lst1(ntdel++) = tnxt;
            td(tnxt).info |= TSRCH;
            tind = tnxt;
            for(snum=0;snum<3;++snum) {
                if (td(tind).vrtx(snum) == v0) break;
            }
            snum = (snum+2)%3;
        }
        else {
            gbl->i2wk_lst2(nskeep++) = sind;
            gbl->i2wk_lst2(nskeep++) = dir;
            v0 = sd(sind).vrtx(dir);
            snum = (snum+1)%3;
        }
        tnxt = td(tind).tri(snum);
        sind = td(tind).side(snum);
        dir = (1+td(tind).sign(snum))/2;    
    } while (sind != sstart);

    /* RESET TSRCH FLAGS */
    for(i=0;i<ntdel;++i)
        td(gbl->i2wk_lst1(i)).info &= ~TSRCH;

    /* ALTER OLD BOUNDARY SIDE & CREATE NEW SIDE */
    sd(nside).vrtx(endpt) = vnum;
    sd(nside).vrtx(1-endpt) = sd(sind).vrtx(1-endpt);
    sd(sind).vrtx(1-endpt) = vnum;
    td(sind).info |= STOUC;
    td(nside).info |= STOUC;
    
    /* ADD NEW SIDE TO BOUNDARY GROUP */
    /* NEED TO REORDER WHEN FINISHED */
    i = getbdrynum(sd(sind).tri(1));
    sd(nside).tri(1) = trinumatbdry(i,sbdry(i)->nel);
    if (sbdry(i)->nel >= sbdry(i)->maxel) {
        output("error");
        *gbl->log << "need to use bigger growth factor (too many boundary sides)" << std::endl;
        exit(1);
    }
    sbdry(i)->el(sbdry(i)->nel++) = nside;
    ++nside;
    
    nskeep = nskeep>>1;
        
    /* APPEND NEW INDICES TO LIST OF INDICES TO USE FOR NEW SIDES & TRIS*/
    for(i=nsdel;i<nskeep;++i)
        gbl->i2wk_lst3(i) = nside +(i-nsdel);
    
    for(i=ntdel;i<nskeep;++i) 
        gbl->i2wk_lst1(i) = ntri +(i-ntdel);
            
    ++ntri;
    ++nside;  
    
    if (ntri > maxvst || nside > maxvst) {
        *gbl->log << "need to use bigger growth factor: too many sides/tris:" << nside << ntri << std::endl;
        output("error");
        exit(1);
    }
    
    /* FIX FIRST AND LAST TRIANGLE BOUNDARY SIDE */
    if (endpt) {
        tind = gbl->i2wk_lst1(0);
        td(tind).side(1) = sind;
        td(tind).sign(1) = 1;
        td(tind).tri(1) = sd(sind).tri(1);
        sd(sind).tri(0) = tind;
        
        tind = gbl->i2wk_lst1(nskeep-1);
        td(tind).side(0) = nside-2;
        td(tind).sign(0) = 1;
        td(tind).tri(0) = sd(nside-2).tri(1);
        sd(nside-2).tri(0) = tind;
 
    }
    else {
        tind = gbl->i2wk_lst1(0);
        td(tind).side(1) = nside-2;
        td(tind).sign(1) = 1;
        td(tind).tri(1) = sd(nside-2).tri(1);  
        sd(nside-2).tri(0) = tind;
        
        tind = gbl->i2wk_lst1(nskeep-1);
        td(tind).side(0) = sind;
        td(tind).sign(0) = 1;
        td(tind).tri(0) = sd(sind).tri(1);
        sd(sind).tri(0) = tind;
    }
    
    
    for(i=0;i<nskeep-1;++i) {
        tind = gbl->i2wk_lst1(i);
        td(tind).info |= TTOUC;

        tnxt = gbl->i2wk_lst1(i+1);
        sind = gbl->i2wk_lst2(i<<1);
        dir = gbl->i2wk_lst2((i<<1) +1);
            
        /* CREATE NEW INFO */
        v0 = sd(sind).vrtx(1-dir);
        td(tind).vrtx(0) = v0;
        vd(v0).tri = tind;
        v0 = sd(sind).vrtx(dir);
        td(tind).vrtx(1) = v0;
        vd(v0).tri = tind;
        td(tind).vrtx(2) = vnum;
        vd(vnum).tri = tind;

        /* SIDE 2 INFO */        
        td(tind).side(2) = sind;
        td(tind).sign(2) = -1 +2*dir;
        
        sd(sind).tri(1-dir) = tind;
        tin = sd(sind).tri(dir);
        td(tind).tri(2) = tin;
        if (tin > -1) {
            for(j=0;j<3;++j) {
                if (td(tin).side(j) == sind) {
                    td(tin).tri(j) = tind;
                    break;
                }
            }
        }
            
        /* CREATE SIDE 0 */
        v0 = sd(sind).vrtx(dir);
        sind1 = gbl->i2wk_lst3(i);
        sd(sind1).tri(0) = tind;
        sd(sind1).vrtx(0) = v0;
        sd(sind1).vrtx(1) = vnum;
        td(sind1).info |= STOUC;


        td(tind).side(0) = sind1;
        td(tind).sign(0) = 1;
        td(tind).tri(0) = tnxt;

        sd(sind1).tri(1) = tnxt;
        td(tnxt).side(1) = sind1;
        td(tnxt).sign(1) = -1;
        td(tnxt).tri(1) = tind;
    }
    
    /* LAST TRIANGLE */
    i = nskeep-1;
    tind = gbl->i2wk_lst1(i);
    td(tind).info |= TTOUC;

    sind = gbl->i2wk_lst2(i<<1);
    dir = gbl->i2wk_lst2((i<<1) +1);
        
    /* CREATE NEW INFO */
    v0 = sd(sind).vrtx(1-dir);
    td(tind).vrtx(0) = v0;
    vd(v0).tri = tind;
    v0 = sd(sind).vrtx(dir);
    td(tind).vrtx(1) = v0;
    vd(v0).tri = tind;
    td(tind).vrtx(2) = vnum;
    vd(vnum).tri = tind;

    /* SIDE 2 INFO */        
    td(tind).side(2) = sind;
    td(tind).sign(2) = -1 +2*dir;
    
    sd(sind).tri(1-dir) = tind;
    tin = sd(sind).tri(dir);
    td(tind).tri(2) = tin;
    if (tin > -1) {
        for(j=0;j<3;++j) {
            if (td(tin).side(j) == sind) {
                td(tin).tri(j) = tind;
                break;
            }
        }
    }
    
    return;
}

int tri_mesh::findtri(const TinyVector<FLT,ND> x, int vnear) {
    int i,j,vn,dir,stoptri,tin,tind;
    int ntdel;
    int tclose,nsurround;
    FLT minclosest,closest;
    int v0,v1,v2;
    FLT dx0,dy0,dx1,dy1,dx2,dy2;
    TinyVector<FLT,3> a;
    
/* TSRCH = 0x100*0x4 */
#if ((-1)&(0x100*0x4))
#define ISSRCH(A) (!((A)&TSRCH))
#define SETSRCH(A) A&=(~TSRCH)
#define CLRSRCH(A) A|=(TSRCH)
#else
#define ISSRCH(A) (((A)&TSRCH))
#define SETSRCH(A) A|=(TSRCH)
#define CLRSRCH(A) A&=(~TSRCH)
#endif
    
    /* HERE WE USE gbl->intwk & gbl->i2wk THIS MUST BE -1 BEFORE USING */
    tind = vd(vnear).tri;
    stoptri = tind;
    dir = 1;
    ntdel = 0;
    do {
        if (intri(tind,x) < area(tind)*10.*EPSILON) goto FOUND;
        SETSRCH(gbl->intwk(tind));
        gbl->i2wk(ntdel++) = tind;
    
        for(vn=0;vn<3;++vn) 
            if (td(tind).vrtx(vn) == vnear) break;
        
        tind = td(tind).tri((vn +dir)%3);
        if (tind < 0) {
            if (dir > 1) break;
            /* REVERSE DIRECTION AND GO BACK TO START */
            ++dir;
            tind = vd(vnear).tri;
            for(vn=0;vn<3;++vn) 
                if (td(tind).vrtx(vn) == vnear) break;

            tind = td(tind).tri((vn +dir)%3);
            if (tind < 0) break;
            stoptri = -1;
        }
    } while(tind != stoptri); 
    
    nsurround = ntdel;
        
    /* DIDN'T FIND TRIANGLE */
    /* NEED TO SEARCH SURROUNDING TRIANGLES */
    for(i=0;i<ntdel;++i) {
        tin = gbl->i2wk(i);
        for(j=0;j<3;++j) {
            tind = td(tin).tri(j);
            if (tind < 0) continue;
            if (ISSRCH(gbl->intwk(tind))) continue;
            SETSRCH(gbl->intwk(tind));
            gbl->i2wk(ntdel++) = tind;            
            if (intri(tind,x) < area(tind)*10.*EPSILON) goto FOUND;
        }
        if (ntdel >= gbl->maxsrch-4) break;
    }
//    std::cerr << "couldn't find tri for point " << x[0] << ' ' << x[1] << ' ' << vnear << std::endl;
    minclosest = -1.0e16;
    tclose = -1;
    for (i=0;i<ntdel;++i) {
        tind = gbl->i2wk(i);
        
        v0 = td(tind).vrtx(0);
        v1 = td(tind).vrtx(1);
        v2 = td(tind).vrtx(2);
        
        dx0 =  (x(0) -vrtx(v0)(0));
        dy0 =  (x(1) -vrtx(v0)(1)); 
        dx1 =  (x(0) -vrtx(v1)(0));
        dy1 =  (x(1) -vrtx(v1)(1));
        dx2 =  (x(0) -vrtx(v2)(0));
        dy2 =  (x(1) -vrtx(v2)(1));
        
        a(0) = (dy2*dx1 -dx2*dy1);
        a(1) = (dy0*dx2 -dx0*dy2);
        a(2) = (dy1*dx0 -dx1*dy0);
        
        /* FIND NEGATIVE SIDE */
        /* CHECK IF 2 SIDES POSITIVE & 1 NEGATIVE */
        if (a(0)*a(1)*a(2) < 0) {
            a(0) /= distance(v2,v1);
            a(1) /= distance(v0,v2);
            a(2) /= distance(v1,v0);
        
            for (int j=0;j<3;++j) {
                if (a(j) < 0.0 && a(j) > minclosest) {
                    minclosest = a(j);
                    tclose = tind;
                    break;
                }
            }
        }
    }
    if (tclose < 0) {
        *gbl->log << "Major Trouble in Findtri " << x << ' ' << vnear << '\n';
        exit(1);
    }

    intri(tclose,x);
    tind = -tclose;
        
FOUND:
    /* RESET gbl->intwkW1 TO -1 */
    for(i=0;i<ntdel;++i) {
        CLRSRCH(gbl->intwk(gbl->i2wk(i)));
    }
 
    return(tind);
}

