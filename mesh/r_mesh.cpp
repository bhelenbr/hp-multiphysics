#include "r_mesh.h"
#include "r_boundary.h"
#include "block.h"
#include <utilities.h>
#include <iostream>
#include <cmath>
#include <input_map.h>
#include <fstream>


void r_mesh::init(input_map& input, r_mesh::gbl *gin) {
    std::string keyword;
    std::istringstream data;
    std::string filename;
    int coarse,ival;
    
    mesh::init(input,gin);
        
    keyword = idprefix + "_coarse";
    if (!input.get(keyword,coarse)) {
        *sim::log << "#Error: not sure whether to initialize for r_mesh for multigrid" << std::endl;
    }
    
    keyword = idprefix + "_r_fadd";
    if (!input.get(keyword,fadd)) {
        input.getwdefault("r_fadd",fadd,1.0);
    }

    keyword = idprefix + "_r_cfl";
    if (!input.get(keyword,r_cfl)) {
        input.getwdefault("r_cfl",r_cfl,0.5);
    }
    
    keyword = idprefix + "_r_output_type";
    if (input.get(keyword,ival)) {
        output_type = static_cast<mesh::filetype>(ival);
    }
    else {
        if (input.get("r_output_type",ival)) {
            output_type = static_cast<mesh::filetype>(ival);
        }
        else {
            output_type = mesh::grid;
        }
    }
    
    /* local storage */    
    ksprg.resize(maxvst);
    kvol.resize(maxvst);
    src.resize(maxvst);
    if (coarse) vrtx_frst.resize(maxvst);
    isfrst = false;
    
    /* BLOCK SHARED INFORMATION */
    gbl_ptr = gin;
    if (!coarse) {
        gbl_ptr->diag.resize(maxvst);
        gbl_ptr->res.resize(maxvst);
        gbl_ptr->res1.resize(maxvst);
    }

    r_sbdry.resize(nsbd);
    for(int i=0;i<nsbd;++i)
        r_sbdry(i) = getnewsideobject(i,input);
    
    return;
}

r_mesh::~r_mesh() {
    for(int i=0;i<nsbd;++i)
        delete r_sbdry(i);
}


void r_mesh::rklaplace() {
    int sind,tind,v0,v1,k;
    FLT dx,dy,l;
    
    for(sind=0;sind<nside;++sind)
        ksprg(sind) = 0.0;

    /* COEFFICIENTS FOR LAPLACE EQUATION */
    /* THIS REQUIRES 2 EVALUATIONS OF SIDE LENGTH FOR EACH SIDE */
    /* BUT IS LOGISTICALLY SIMPLE          */            
    for(tind=0;tind<ntri;++tind) {
        for(k=0;k<3;++k) {
            sind = td(tind).side(k);
            v0 = sd(sind).vrtx(0);
            v1 = sd(sind).vrtx(1);
            dx = vrtx(v1)(0) -vrtx(v0)(0);
            dy = vrtx(v1)(1) -vrtx(v0)(1);        
            l  = (dx*dx +dy*dy)/area(tind);

            ksprg(sind) -= l;
            sind = td(tind).side((k+1)%3);
            ksprg(sind) += l;
            sind = td(tind).side((k+2)%3);
            ksprg(sind) += l;
        }
    }
    
    return;
}

#ifdef FOURTH
void r_mesh::calc_kvol() {
    int last_phase, mp_phase;
    
    for(int i=0;i<nvrtx;++i) 
        kvol(i) = 0.0;

    for(int tind=0;tind<ntri;++tind) 
        for(int i=0;i<3;++i) 
            kvol(td(tind).vrtx(i)) += area(tind);

    for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
        vmsgload(boundary::all_phased,mp_phase,boundary::symmetric,kvol.data(),0,0,1);
        vmsgpass(boundary::all_phased,mp_phase, boundary::symmetric);
        last_phase = true;
        last_phase &= vmsgwait_rcv(boundary::all_phased,mp_phase, boundary::symmetric, boundary::average, kvol.data(),0,0,1);
    }

    for(int i=0;i<nvrtx;++i)
        kvol(i) = 1./kvol(i);
}
#endif

void r_mesh::rksprg() {
    int sind,v0,v1;
    double dx,dy;

    /* 2D SPRING CONSTANTS FINE MESH*/
    for(sind=0;sind<nside;++sind) {
        v0 = sd(sind).vrtx(0);
        v1 = sd(sind).vrtx(1);
        dx = vrtx(v1)(0) -vrtx(v0)(0);
        dy = vrtx(v1)(1) -vrtx(v0)(1);
        ksprg(sind) = 1.0/(dx*dx +dy*dy);
    }

    return;
}

void r_mesh::rkmgrid(Array<mesh::transfer,1> &fv_to_ct, r_mesh *fmesh) {
    int i,j,sind,tind,tind0,tind1,v0,v1;    
    
    /* Load Diag Locally */
    Array<FLT,1> diag;
    diag.reference(gbl_ptr->diag);
    Array<TinyVector<FLT,ND>,1> res1;
    res1.reference(gbl_ptr->res1);
    
    /* TEMPORARILY USE DIAG TO STORE DIAGONAL SUM */
    for(i=0;i<fmesh->nvrtx;++i)
        diag(i) = 0.0;

    /* FORM KIJ SUM AT VERTICES */
    for(sind=0;sind<fmesh->nside;++sind) {
        v0 = fmesh->sd(sind).vrtx(0);
        v1 = fmesh->sd(sind).vrtx(1);
        diag(v0) += fmesh->ksprg(sind);
        diag(v1) += fmesh->ksprg(sind);
    }

    for(i=0;i<nside;++i)
        ksprg(i) = 0.0;

    /* LOOP THROUGH FINE VERTICES    */
    /* TO CALCULATE KSPRG ON COARSE MESH */    
    for(i=0;i<fmesh->nvrtx;++i) {
        tind = fv_to_ct(i).tri;
        for(j=0;j<3;++j) {
            sind = td(tind).side(j);
            ksprg(sind) -= fv_to_ct(i).wt(j)*fv_to_ct(i).wt((j+1)%3)*diag(i);
        }
    }

    /* LOOP THROUGH FINE SIDES */
    for(i=0;i<fmesh->nside;++i) {
        v0 = fmesh->sd(i).vrtx(0);
        v1 = fmesh->sd(i).vrtx(1);
        tind0 = fv_to_ct(v0).tri;
        tind1 = fv_to_ct(v1).tri;
                        
        /* TEMPORARILY STORE WEIGHTS FOR FINE POINTS (0,1) */
        /* FOR EACH COARSE VERTEX */
        for(j=0;j<3;++j)  {
            res1(td(tind1).vrtx(j))(0) = 0.0;
            res1(td(tind0).vrtx(j))(1) = 0.0;
        }

        for(j=0;j<3;++j)  {
            res1(td(tind0).vrtx(j))(0) = fv_to_ct(v0).wt(j);
            res1(td(tind1).vrtx(j))(1) = fv_to_ct(v1).wt(j);
        }
        
        /* LOOP THROUGH COARSE TRIANGLE 0 SIDES */
        for(j=0;j<3;++j) {
            sind = td(tind0).side(j);
            ksprg(sind) += fmesh->ksprg(i)*
                (res1(sd(sind).vrtx(0))(0)*res1(sd(sind).vrtx(1))(1)
                +res1(sd(sind).vrtx(1))(0)*res1(sd(sind).vrtx(0))(1));
        }

        if (tind0 != tind1) {
            for(j=0;j<3;++j) {
                sind = td(tind1).side(j);
                if (sd(sind).tri(0) +sd(sind).tri(1) != tind0 +tind1) {
                    ksprg(sind) += fmesh->ksprg(i)*
                        (res1(sd(sind).vrtx(0))(0)*res1(sd(sind).vrtx(1))(1)
                        +res1(sd(sind).vrtx(1))(0)*res1(sd(sind).vrtx(0))(1));

                }
            }
        }
    }

    return;
}


void r_mesh::update() {
    int i,n;
    
    r_mesh::rsdl();
        
    Array<FLT,1> diag;
    diag.reference(gbl_ptr->diag);
    Array<TinyVector<FLT,ND>,1> res;
    res.reference(gbl_ptr->res);

    for(i=0;i<nvrtx;++i)
        for(n=0;n<ND;++n)
            vrtx(i)(n) -= diag(i)*res(i)(n);
}

void r_mesh::zero_source() {
    int i,n;
    
    for(i=0;i<nvrtx;++i) 
        for(n=0;n<ND;++n) 
            src(i)(n) = 0.0;
                
    return;
}

void r_mesh::sumsrc() {
    int i,n;
    
    Array<TinyVector<FLT,ND>,1> res;
    res.reference(gbl_ptr->res);

    for(i=0;i<nvrtx;++i) 
        for(n=0;n<ND;++n)
            src(i)(n) = -1.0*res(i)(n);
    
    return;
}

void r_mesh::mg_getfres(Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, r_mesh *fmesh) {
    int i,j,n,tind,v0;
    
    Array<TinyVector<FLT,ND>,1> fres;
    fres.reference(gbl_ptr->res);
        
    for(i=0;i<nvrtx;++i)
        for(n=0;n<ND;++n)
            src(i)(n) = 0.0;
            
    /* LOOP THROUGH FINE VERTICES TO CALCULATE RESIDUAL  */
    for(i=0;i<fmesh->nvrtx;++i) {
        tind = fv_to_ct(i).tri;
        for(j=0;j<3;++j) {
            v0 = td(tind).vrtx(j);
            for(n=0;n<ND;++n)
                src(v0)(n) += fadd*fv_to_ct(i).wt(j)*fres(i)(n);
        }
    }
    
    /* LOOP THROUGH fv_to_ct VERTICES    */
    /* TO CALCULATE VRTX ON fv_to_ct MESH */
    for(i=0;i<nvrtx;++i) {
        tind = cv_to_ft(i).tri;

        for(n=0;n<ND;++n)
            vrtx(i)(n) = 0.0;
            
        for(j=0;j<3;++j) {
            for(n=0;n<ND;++n)
                vrtx(i)(n) += cv_to_ft(i).wt(j)*fmesh->vrtx(fmesh->td(tind).vrtx(j))(n);
        }
        
        for(n=0;n<ND;++n)
            vrtx_frst(i)(n) = vrtx(i)(n);
    }
    isfrst = true;
    
    return;
}

void r_mesh::mg_getcchng(Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, r_mesh *cmesh) {
    int i,j,n,ind,tind;
    int last_phase, mp_phase;  

    Array<TinyVector<FLT,ND>,1> res;
    res.reference(gbl_ptr->res);
    
    /* DETERMINE CORRECTIONS ON COARSE MESH    */    
    for(i=0;i<cmesh->nvrtx;++i)
        for(n=0;n<ND;++n) 
            cmesh->vrtx_frst(i)(n) -= cmesh->vrtx(i)(n);

    /* LOOP THROUGH FINE VERTICES    */
    /* TO DETERMINE CHANGE IN SOLUTION */    
    for(i=0;i<nvrtx;++i) {
        
        for(n=0;n<ND;++n)
            res(i)(n) = 0.0;
        
        tind = fv_to_ct(i).tri;
        
        for(j=0;j<3;++j) {
            ind = cmesh->td(tind).vrtx(j);
            for(n=0;n<ND;++n) 
                res(i)(n) -= fv_to_ct(i).wt(j)*cmesh->vrtx_frst(ind)(n);
        }
    }
    
    for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
        for(i=0;i<nsbd;++i)
            sbdry(i)->vloadbuff(boundary::partitions,(FLT *) res.data(),0,1,2);
        
        for(i=0;i<nsbd;++i) 
            sbdry(i)->comm_prepare(boundary::partitions,mp_phase, boundary::symmetric);

        for(i=0;i<nsbd;++i) 
            sbdry(i)->comm_exchange(boundary::partitions,mp_phase, boundary::symmetric);
                        
        last_phase = true;
        for(i=0;i<nsbd;++i) {
            last_phase &= sbdry(i)->comm_wait(boundary::partitions,mp_phase, boundary::symmetric);
            sbdry(i)->vfinalrcv(boundary::partitions,mp_phase, boundary::symmetric,boundary::average,(FLT *) res.data(),0,1,2);
        }
    }

    for(i=0;i<nvrtx;++i)
        for(n=0;n<ND;++n) 
            vrtx(i)(n) += res(i)(n);
   
    return;
}

FLT r_mesh::maxres() {
    int i,n;
    FLT mxr[ND];
    FLT sum;
    
    Array<TinyVector<FLT,ND>,1> res;
    res.reference(gbl_ptr->res);

    for(n=0;n<ND;++n)
        mxr[n] = 0.0;

    for(i=0;i<nvrtx;++i)
        for(n=0;n<ND;++n)
            mxr[n] = MAX(mxr[n],fabs(res(i)(n)));
            
    sum = 0.0;
    for(n=0;n<ND;++n) {
        *sim::log << ' ' << mxr[n] << ' ';
        sum += mxr[n];
    }
    
            
    return(sum);
}


void r_mesh::tadvance(bool coarse, Array<mesh::transfer,1> &fv_to_ct,Array<mesh::transfer,1> &cv_to_ft, r_mesh *fmesh) {

    if (!coarse) {
        rklaplace();
#ifdef FOURTH
        calc_kvol();
#endif         
        zero_source();
        rsdl();
        sumsrc();
        moveboundaries();
    }
    else {
#ifdef GEOMETRIC
        rklaplace();
#else
        /* USE MULTIGRID INTERPOLATION (ALGEBRAIC) */
        /* MUST BE DONE THIS WAY FOR SPRING METHOD */
        /* SETUP FIRST MESH */
        rkmgrid(fv_to_ct,fmesh);
#endif

#ifdef  FOURTH
        calc_kvol();
#endif
    }

    return;
}

void r_mesh::moveboundaries() {
    
    /* MOVE BOUNDARY POSITIONS */
    for(int i=0;i<nsbd;++i)
        r_sbdry(i)->tadvance();
    
    return;
}

void r_mesh::rsdl() {
    int last_phase, mp_phase;
    int i,n,sind,v0,v1;
    FLT dx,dy;
    Array<TinyVector<FLT,ND>,1> res;
    res.reference(gbl_ptr->res);
    
#ifdef FOURTH
    /*************************************/ 
    /* FOURTH ORDER MESH MOVEMENT SCHEME */
    /*************************************/
    Array<TinyVector<FLT,ND>,1> res1;
    res1.reference(gbl_ptr->res1);
    
    for(i=0;i<nvrtx;++i)
        for(n=0;n<ND;++n)
            res1(i)(n) = 0.0;

    for (sind = 0; sind < nside; ++sind) {
        v0 = sd(sind).vrtx(0);
        v1 = sd(sind).vrtx(1);

        dx = ksprg(sind)*(vrtx(v1)(0)-vrtx(v0)(0));
        dy = ksprg(sind)*(vrtx(v1)(1)-vrtx(v0)(1));

        res1(v0)(0) -= dx;
        res1(v0)(1) -= dy;

        res1(v1)(0) += dx;
        res1(v1)(1) += dy;        
    }
    
    /* APPLY DIRICHLET BOUNDARY CONDITIONS */
    for(i=0;i<nsbd;++i)
        r_sbdry(i)->fixdx2();

    for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
        vmsgload(boundary::all_phased, mp_phase, boundary::symmetric,(FLT *) gbl_ptr->res1.data(),0,1,2);
        vmsgpass(boundary::all_phased, mp_phase, boundary::symmetric);
        last_phase = true;
        last_phase &= vmsgwait_rcv(boundary::all_phased, mp_phase/3, boundary::symmetric,  boundary::average, (FLT *) gbl_ptr->res1.data(),0,1,2);
    }
  
    /* DIVIDE BY VOLUME FOR AN APPROXIMATION TO D^2/DX^2 */
    for(i=0;i<nvrtx;++i)
        for(n=0;n<ND;++n) 
            res1(i)(n) *= kvol(i);
    
    for(i=0;i<nvrtx;++i)
        for(n=0;n<ND;++n)
            res(i)(n) = 0.0;

    for (sind = 0; sind < nside; ++sind) {
        v0 = sd(sind).vrtx(0);
        v1 = sd(sind).vrtx(1);

        dx = ksprg(sind)*(res1(v1)(0)-res1(v0)(0));
        dy = ksprg(sind)*(res1(v1)(1)-res1(v0)(1));

        res(v0)(0) -= dx;
        res(v0)(1) -= dy;

        res(v1)(0) += dx;
        res(v1)(1) += dy;        
    }
#else   
    res(Range(0,nvrtx-1)) = 0.0;

    int lnside = nside;
    FLT lksprg;
    for(int sind=0;sind<lnside;++sind) {
        v0 = sd(sind).vrtx(0);
        v1 = sd(sind).vrtx(1);

        lksprg = ksprg(sind);
        dx = lksprg*(vrtx(v1)(0)-vrtx(v0)(0));
        dy = lksprg*(vrtx(v1)(1)-vrtx(v0)(1));

        res(v0)(0) -= dx;
        res(v0)(1) -= dy;

        res(v1)(0) += dx;
        res(v1)(1) += dy;
    }
#endif

    /* CALCULATE DRIVING TERM ON FIRST ENTRY TO COARSE MESH */
    if (isfrst) {
        for(i=0;i<nvrtx;++i) 
            for(n=0;n<ND;++n)
                src(i)(n) -= res(i)(n);

        isfrst = false;
    }

    /* ADD IN MULTIGRID SOURCE OR FMESH SOURCE */
    for(i=0;i<nvrtx;++i) 
        for(n=0;n<ND;++n)
            res(i)(n) += src(i)(n);

    /* APPLY DIRICHLET BOUNDARY CONDITIONS */
    for(i=0;i<nsbd;++i)
        r_sbdry(i)->dirichlet();

    for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
        vmsgload(boundary::all_phased,mp_phase, boundary::symmetric,(FLT *) gbl_ptr->res.data(),0,1,2);
        vmsgpass(boundary::all_phased,mp_phase, boundary::symmetric);
        last_phase = true;
        last_phase &= vmsgwait_rcv(boundary::all_phased,mp_phase, boundary::symmetric, boundary::average, (FLT *) gbl_ptr->res.data(),0,1,2);
    }    

    return;
}


void r_mesh::setup_preconditioner() {
    int last_phase, mp_phase;
    int i,v0,v1,sind;
    Array<FLT,1> diag;
    diag.reference(gbl_ptr->diag);
    
#ifdef FOURTH
    /**************************************************/
    /* DETERMINE MESH MOVEMENT TIME STEP              */
    /**************************************************/
    for(i=0;i<nvrtx;++i)
        diag(i) = 0.0;
        
    for(sind=0;sind<nside;++sind) {
        v0 = sd(sind).vrtx(0);
        v1 = sd(sind).vrtx(1);
        diag(v0) += fabs(ksprg(sind));
        diag(v1) += fabs(ksprg(sind));
    }

    
    for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
        vmsgload(boundary::all_phased, mp_phase, boundary::symmetric,gbl_ptr->diag.data(),0,0,1);
        vmsgpass(boundary::all_phased, mp_phase, boundary::symmetric);
        last_phase = true;
        last_phase &= vmsgwait_rcv(boundary::all_phased, mp_phase, boundary::symmetric,  boundary::average, gbl_ptr->diag.data(),0,0,1);
    }
        
    Array<TinyVector<FLT,ND>,1> res1;
    res1.reference(gbl_ptr->res1);

    for(i=0;i<nvrtx;++i)
        res1(i)(0) = diag(i)*kvol(i);
        
    for(i=0;i<nvrtx;++i)
        diag(i) = 0.0;
    
    for(sind=0;sind<nside;++sind) {
        v0 = sd(sind).vrtx(0);
        v1 = sd(sind).vrtx(1);
        diag(v0) += fabs(ksprg(sind))*(res1(v0)(0) +fabs(ksprg(sind))*kvol(v1));
        diag(v1) += fabs(ksprg(sind))*(res1(v1)(0) +fabs(ksprg(sind))*kvol(v0));
    }
#else
    /**************************************************/
    /* DETERMINE MESH MOVEMENT TIME STEP              */
    /**************************************************/
    for(i=0;i<nvrtx;++i)
        diag(i) = 0.0;

    /* FORM TIME STEP FOR MV_UPDATE */
    for(sind=0;sind<nside;++sind) {
        v0 = sd(sind).vrtx(0);
        v1 = sd(sind).vrtx(1);
        diag(v0) += fabs(ksprg(sind));
        diag(v1) += fabs(ksprg(sind));
    }
#endif

     for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
        vmsgload(boundary::all_phased,mp_phase, boundary::symmetric,gbl_ptr->diag.data(),0,0,1);
        vmsgpass(boundary::all_phased,mp_phase, boundary::symmetric);
        last_phase = true;
        last_phase &= vmsgwait_rcv(boundary::all_phased,mp_phase, boundary::symmetric, boundary::average, gbl_ptr->diag.data(),0,0,1);
    }   

    for(i=0;i<nvrtx;++i)
        diag(i) = r_cfl/diag(i);
    
    return;
}
