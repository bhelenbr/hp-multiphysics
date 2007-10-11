/*
 *  srebay.cpp
 *  mblock
 *
 *  Created by helenbrk on Wed Sep 12 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */
#include "mesh.h"
#include <utilities.h>
#include <assert.h>
#include <math.h>
#include <blitz/tinyvec-et.h>

#define REBAY
#define NO_DEBUG_ADAPT
#define VERBOSE

#ifdef DEBUG_ADAPT
int adapt_count = 0;
static std::string adapt_file;
#endif

/*  fscr1 = ratio of triangle size to target
    sd(i).info = triangle refinement queue
    vd(tind).info = back reference from triangle into queue
*/
    
void tri_mesh::rebay(FLT tolsize) {
    int i,j,n,tind,tfind,v0,v1,v2,vnear,nsnew,ntnew,snum,intrcnt,err;
    TinyVector<FLT,ND> xpt,dx,xmid,xdif,rn;
    TinyVector<FLT,3> wt;
    FLT maxvl;
    FLT norm,p,q,s1sq,s2sq,rad1,rad2,rs;
    FLT densty,cirrad,arg,dist,rsign;
    
    /* TO ADJUST FOR CIRCUMSCRIBED RADIUS */
    tolsize /= sqrt(3.);
         
    /* COUNTER TO SEE HOW MANY VERTICES INSERTED */
    intrcnt = 0;
    
    /* CLASSIFY TRIANGLES AS ACCEPTED OR UNACCEPTED */
    nlst = 0;
    for(i=0;i<ntri;++i) {
        if (td(i).info&TDLTE) continue;
        maxvl = vlngth(td(i).vrtx(0));
        maxvl = MAX(maxvl,vlngth(td(i).vrtx(1)));
        maxvl = MAX(maxvl,vlngth(td(i).vrtx(2)));
        gbl->fltwk(i) = circumradius(i)/maxvl;        
        if (gbl->fltwk(i) > tolsize) putinlst(i);
    }

    /* BEGIN REFINEMENT ALGORITHM */
    while (nlst > 0) {
        for(i=nlst-1;i>=0;--i) {
             /* FIND TRIANGLE FACES ON BOUNDARY FIRST */
            for(j=0;j<3;++j) {
                tind = td(sd(i).info).tri(j);
                if (tind < 0)  {
                    snum = j;
                    tind = sd(i).info;
                    goto TFOUND;
                }
            }
            /* FIND TRIANGLE FACES ON BOUNDARY OF ACCEPTED REGIONS */
            for(j=0;j<3;++j) {
                tind = td(sd(i).info).tri(j);
                if (vd(tind).info == -1)  {
                    snum = j;
                    tind = sd(i).info;
                    goto TFOUND;
                }
            }
        }
        *gbl->log << "Didn't find triangle???" << std::endl;
                
        TFOUND:
        
        /* CHECK THAT NOT ALL SIDES ARE ACCEPTED */
        tfind = 0;
        for (j=0;j<3;++j) 
            if (td(tind).tri(j) < 0 || vd(td(tind).tri(j)).info == -1) ++tfind;

        if (tfind == 3) {
            tkoutlst(tind);
            continue;
        }
        
        if (nvrtx > maxvst -2) {
            *gbl->log << "too many vertices" << std::endl;
            exit(1);
        }
        
        /* THIS IS TRIANGLE BASED REFINEMENT */
        /*	FIND ACCEPTED EDGE ON TRIANGLE & BUILD OFF OF IT */
        v0 = td(tind).vrtx((snum+1)%3);
        v1 = td(tind).vrtx((snum+2)%3);
        v2 = td(tind).vrtx(snum);
        
#ifdef REBAY
        /* USE REBAY'S ALGORITHM FOR INSERT POINT */
        p = 0.0;
        for(n=0;n<ND;++n) {
            dx(n) = vrtx(v0)(n) -vrtx(v1)(n);
            p += pow(dx(n),2);
        }
        p = 0.5*sqrt(p);
        
        s1sq = 0.0;
        for(n=0;n<ND;++n) {
            dx(n) = vrtx(v2)(n) -vrtx(v1)(n);
            s1sq += pow(dx(n),2);
        }
        s1sq *= 0.25;

        s2sq = 0.0;
        for(n=0;n<ND;++n) {
            dx(n) = vrtx(v2)(n) -vrtx(v0)(n);
            s2sq += pow(dx(n),2);
        }
        s2sq *= 0.25;            
        circumcenter(tind,xpt);
        if (p*p > s1sq +s2sq) goto INSRT;
        
        q = 0.0;
        for(n=0;n<ND;++n) {
            xmid(n) = .5*(vrtx(v0)(n) +vrtx(v1)(n));
            dx(n) = xpt(n) -xmid(n);
            q += pow(dx(n),2);
        }
        q = sqrt(q);
        if (q < p) goto INSRT;

        densty = 0.5*(vlngth(v0)  +vlngth(v1))/sqrt(3.0);
        rad1 = MAX(densty,p);
        rad2 = .5*(p*p  +q*q)/q;
        cirrad = MIN(rad1,rad2);
        arg = fabs(cirrad*cirrad  -p*p);
        dist = cirrad  +sqrt(arg);
        rs = 0.0;
        for(n=0;n<ND;++n) {
            dx(n) = vrtx(v1)(n) -vrtx(v0)(n);
            rs += pow(dx(n),2);
        }
        rs = 1./sqrt(rs);
        rn(0) =  dx(1)*rs;
        rn(1) = -dx(0)*rs;
        for(n=0;n<ND;++n)
            xdif(n) = vrtx(v2)(n)  -xmid(n);
        rsign = 1.;
        if (xdif(0)*rn(0) +xdif(1)*rn(1) < 0.) rsign = -1.;
        for(n=0;n<ND;++n)
            xpt(n) = xmid(n) +rsign*dist*rn(n);
#elif defined(MIDPOINT)
            /* MIDPOINT RULE (VERY SIMPLE) */
        for(n=0;n<ND;++n)
            xpt(n) = 0.5*(vrtx(v1)(n) +vrtx(v2)(n));
#else
        densty = 0.5*(vlngth(v0)  +vlngth(v1))/sqrt(3.0);
        rs = 0.0;
        for(n=0;n<ND;++n) {
            xmid(n) = .5*(vrtx(v0)(n) +vrtx(v1)(n));
            dx(n) = vrtx(v2)(n) -xmid(n);
            rs += pow(dx(n),2);
        }
        rs = sqrt(rs);
        if (densty > rs) {
            *gbl->log << "just checking\n";
            densty = 0.0;
            for(n=0;n<ND;++n) {
                densty += pow(vrtx(v1)(n) -vrtx(v0)(n),2);
            }
            densty = sqrt(densty);
        }
        rs = densty/rs;

        xpt = xmid +rs*dx;
#endif

INSRT:
        /* INSERT POINT */
        for(n=0;n<ND;++n)
            vrtx(nvrtx)(n) = xpt(n);

#ifdef DEBUG_ADAPT
        *gbl->log << "Inserting interior side " << intrcnt << ' ' << adapt_count << ' ';
        for(n=0;n<ND;++n)
            *gbl->log << vrtx(nvrtx)(n) << ' ';
        *gbl->log << std::endl;
#endif
        
        dist = qtree.nearpt(vrtx(nvrtx).data(),vnear);
        norm = 0.0;
        for (n=0;n<ND;++n)
            norm += fabs(vrtx(nvrtx)(n));
        if (dist < 100.0*EPSILON*norm) {
#ifdef VERBOSE
            *gbl->log << "#Point to close to insert " << dist << std::endl;
            *gbl->log << vrtx(v0) << std::endl; 
            *gbl->log << vrtx(v1) << std::endl;
            *gbl->log << vrtx(v2) << std::endl;
            *gbl->log << vrtx(nvrtx) << std::endl;
            *gbl->log << td(tind).vrtx <<  snum << std::endl;
            for(j=0;j<3;++j) {
                if (td(tind).tri(j) < 0 || vd(td(tind).tri(j)).info == -1) {
                    *gbl->log << "side " << j << " is accepted\n";
                }
            }            
#endif
            tkoutlst(tind);
            continue;
        }
            
        tfind = findtri(xpt,vnear);
        if (tfind < 0) {
#ifdef VERBOSE
            *gbl->log << "#Warning: Trying to insert outside domain " << std::endl;
            *gbl->log << vrtx(v0) << std::endl; 
            *gbl->log << vrtx(v1) << std::endl;
            *gbl->log << vrtx(v2) << std::endl;
            *gbl->log << vrtx(nvrtx) << std::endl;
            *gbl->log << td(tind).vrtx << snum << std::endl;
            for(j=0;j<3;++j) {
                if (td(tind).tri(j) < 0 || vd(td(tind).tri(j)).info == -1) {
                    *gbl->log << "side " << j << " is accepted\n";
                }
            }            
#endif
            tkoutlst(tind);
            continue;
        }

        getwgts(wt);
        vlngth(nvrtx) = 0.0;
        for(i=0;i<3;++i)
            vlngth(nvrtx) += wt(i)*vlngth(td(tfind).vrtx(i));

        err = insert(nvrtx,tfind);
        if (!err) {
            /* ADD POINT TO QUADTREE */
            td(nvrtx).info |= VTOUC;
            qtree.addpt(nvrtx);
            nsnew = gbl->i2wk_lst3(-1) +3;
            ntnew = gbl->i2wk_lst1(-1) +2;
            ++intrcnt;
        }
        else {
#ifdef VERBOSE
            *gbl->log << "#Warning: Makes Bad Triangle " << std::endl;
            *gbl->log << vrtx(v0) << std::endl; 
            *gbl->log << vrtx(v1) << std::endl;
            *gbl->log << vrtx(v2) << std::endl;
            *gbl->log << vrtx(nvrtx) << std::endl;
            *gbl->log << td(tind).vrtx << snum << std::endl;
            for(j=0;j<3;++j) {
                if (td(tind).tri(j) < 0 || vd(td(tind).tri(j)).info == -1) {
                    *gbl->log << "side " << j << " is accepted\n";
                }
            }            
#endif
            tkoutlst(tind);
            continue;
        }
        ++nvrtx;
            
        for(i=0;i<gbl->i2wk_lst1(-1);++i) 
            if (vd(gbl->i2wk_lst1(i)).info > -1) tkoutlst(gbl->i2wk_lst1(i));
        
        if (vd(tind).info > -1) tkoutlst(tind);
            
            
        for(i=0;i<ntnew;++i) {
            tind = gbl->i2wk_lst1(i);
            maxvl = vlngth(td(tind).vrtx(0));
            maxvl = MAX(maxvl,vlngth(td(tind).vrtx(1)));
            maxvl = MAX(maxvl,vlngth(td(tind).vrtx(2)));
            gbl->fltwk(tind) = circumradius(tind)/maxvl;
            if (gbl->fltwk(tind) > tolsize) putinlst(tind);
        }
#ifdef DEBUG_ADAPT
        std::ostringstream nstr;
        nstr << adapt_count++ << std::flush;
        adapt_file = "adapt" +nstr.str() + "_" +gbl->idprefix;
        nstr.str("");
        output(adapt_file.c_str(),grid);
#endif
        
    }
        
    *gbl->log << "#Rebay finished: new interior points " << intrcnt << std::endl;

    return;
}

void tri_mesh::bdry_rebay(FLT tolsize) {
    int sind,v0,count,el,psifxpt;
    FLT psi;
    TinyVector<FLT,tri_mesh::ND> endpt;
    
    /* REFINE BOUNDARY SIDES */
    for(int bnum=0;bnum<nsbd;++bnum) {
        count = 0;

        if (!sbdry(bnum)->is_frst()) {
            sbdry(bnum)->sndtype() = boundary::int_msg;
            sbdry(bnum)->sndsize() = 0;
            sbdry(bnum)->comm_prepare(boundary::all,0,boundary::master_slave);
            continue;
        }
        
        /* CHECK THAT ENDPOINTS ARE ON CURVE JUST TO BE SURE */
        sind = sbdry(bnum)->el(0);
        v0 = sd(sind).vrtx(0);
        endpt = vrtx(v0);
        sbdry(bnum)->mvpttobdry(0,-1.0,endpt);
        endpt -= vrtx(v0);
        if (fabs(endpt(0)) +fabs(endpt(1)) > FLT_EPSILON*(fabs(qtree.xmax(0)-qtree.xmin(0)) +fabs(qtree.xmax(1)-qtree.xmin(1)))) {
            *gbl->log << "first endpoint of boundary " << sbdry(bnum)->idprefix << " does not seem to be on curve\n";
        }
        
        sind = sbdry(bnum)->el(sbdry(bnum)->nel-1);
        v0 = sd(sind).vrtx(1);
        endpt = vrtx(v0);
        sbdry(bnum)->mvpttobdry(sbdry(bnum)->nel-1,1.0,endpt);
        endpt -= vrtx(v0);
        if (fabs(endpt(0)) +fabs(endpt(1)) > FLT_EPSILON*(fabs(qtree.xmax(0)-qtree.xmin(0)) +fabs(qtree.xmax(1)-qtree.xmin(1)))) {
            *gbl->log << "last endpoint of boundary " << sbdry(bnum)->idprefix << " does not seem to be on curve\n";
        }
        
        nlst = 0;
        for(int indx=0;indx<sbdry(bnum)->nel;++indx) {
            sind = sbdry(bnum)->el(indx);
            if (td(sind).info&SDLTE) continue;
            gbl->fltwk(sind) = distance(sd(sind).vrtx(0),sd(sind).vrtx(1))/MAX(vlngth(sd(sind).vrtx(0)),vlngth(sd(sind).vrtx(1)));
            if (gbl->fltwk(sind) > tolsize) putinlst(sind);
        }
        
        /* SKIP FIRST SPOT SO CAN SEND LENGTH FIRST */
        sbdry(bnum)->sndsize() = 1;
        while (nlst > 0) {
            // START WITH LARGEST SIDE LENGTH RATIO
            sind = sd(nlst-1).info;
            el = getbdryel(sd(sind).tri(1));
            
            /* FOR NOW INSERTION POINT IN MIDDLE */
            /* FIXED POINT ARITHMETIC SO I CAN PASS AN INTEGER */
            psi = 0.0;
            psifxpt = static_cast<int>(256*psi);
            psi = psifxpt/256.0;
            sbdry(bnum)->mvpttobdry(el,psi,vrtx(nvrtx));
            vlngth(nvrtx) = 0.5*((1.-psi)*vlngth(sd(sind).vrtx(0)) +(1.+psi)*vlngth(sd(sind).vrtx(1)));

#ifdef DEBUG_ADAPT
            *gbl->log << "Inserting boundary side " << count << ' ' << adapt_count << ' ';
            for(int n=0;n<ND;++n)
                *gbl->log << vrtx(nvrtx)(n) << ' ';
            *gbl->log << " el " << el << " psi " << psi << std::endl;
#endif
            /* INSERT POINT */
            bdry_insert(nvrtx,sind);
            ++nvrtx;
            ++count;
            
            /* STORE ADAPTATION INFO FOR COMMUNICATION BOUNDARIES */
            sbdry(bnum)->isndbuf(sbdry(bnum)->sndsize()++) = el;
            sbdry(bnum)->isndbuf(sbdry(bnum)->sndsize()++) = psifxpt;

            /* UPDATE MODIFIED SIDE */
            tkoutlst(sind);
            gbl->fltwk(sind) = distance(sd(sind).vrtx(0),sd(sind).vrtx(1))/MAX(vlngth(sd(sind).vrtx(0)),vlngth(sd(sind).vrtx(1)));
            if (gbl->fltwk(sind) > tolsize) putinlst(sind);
            
            /* UPDATE NEW BOUNDARY SIDE */
            sind = sbdry(bnum)->el(sbdry(bnum)->nel -1);
            gbl->fltwk(sind) = distance(sd(sind).vrtx(0),sd(sind).vrtx(1))/MAX(vlngth(sd(sind).vrtx(0)),vlngth(sd(sind).vrtx(1)));
            if (gbl->fltwk(sind) > tolsize) putinlst(sind);
#ifdef DEBUG_ADAPT
            std::ostringstream nstr;
            nstr << adapt_count++ << std::flush;
            adapt_file = "adapt" +nstr.str() + "_" +gbl->idprefix;
            nstr.str("");
            output(adapt_file.c_str(),grid);
#endif
        }
        sbdry(bnum)->isndbuf(0) = sbdry(bnum)->sndsize();
        sbdry(bnum)->sndtype() = boundary::int_msg;
        sbdry(bnum)->comm_prepare(boundary::all,0,boundary::master_slave);
        *gbl->log << "#Boundary refinement finished, " << sbdry(bnum)->idnum << ' ' << count << " sides added" << std::endl;

    }
    
    return;
}

void tri_mesh::bdry_rebay1() {
    int i,sind;
    int el, sndsize;
    int nel_bgn, psifxpt;
    FLT psi;
    
    /* REFINE MATCHING BOUNDARIES */
    for(int bnum=0;bnum<nsbd;++bnum) {
        
        sbdry(bnum)->comm_wait(boundary::all,0,boundary::master_slave);

        if (sbdry(bnum)->is_frst() || !sbdry(bnum)->is_comm()) continue;        
        
        nel_bgn = sbdry(bnum)->nel;
        
        sndsize = sbdry(bnum)->ircvbuf(0,0);
        for(i=1;i<sndsize;i+=2) {
            el = sbdry(bnum)->ircvbuf(0,i);
            if (el < nel_bgn) el = nel_bgn -1 -el;

            sind = sbdry(bnum)->el(el);
            psifxpt = -sbdry(bnum)->ircvbuf(0,i+1);
            psi = psifxpt/256.0;
            sbdry(bnum)->mvpttobdry(el,psi,vrtx(nvrtx));
            vlngth(nvrtx) = 0.5*((1.-psi)*vlngth(sd(sind).vrtx(0)) +(1.+psi)*vlngth(sd(sind).vrtx(1)));
#ifdef DEBUG_ADAPT
            *gbl->log << "Inserting boundary side " << i << ' ' << adapt_count << ' ';
            for(int n=0;n<ND;++n)
                *gbl->log << vrtx(nvrtx)(n) << ' ';
            *gbl->log << " el " << el << " psi " << psi << std::endl;
#endif
            bdry_insert(nvrtx,sind,1);
            ++nvrtx;
            
#ifdef DEBUG_ADAPT
            std::ostringstream nstr;
            nstr << adapt_count++ << std::flush;
            adapt_file = "adapt" +nstr.str() + "_" +gbl->idprefix;
            nstr.str("");
            output(adapt_file.c_str(),grid);
#endif
        }
        *gbl->log << "#Slave boundary refinement finished, " << sbdry(bnum)->idnum << ' ' << (sndsize-1)/2 << " sides added" << std::endl;

    }
    
    return;
}    

    
