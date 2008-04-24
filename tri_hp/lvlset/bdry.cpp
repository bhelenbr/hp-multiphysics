#include "bdry_lvlset.h"
#include <myblas.h>

/*************************************************/
/* SET DIRICHLET BOUNDARY VALUES & FLUXES ********/
/* (THINGS THAT ARE INDEPENDENT OF THE SOLUTION) */
/* BUT NOT INDEPENDENT OF TIME/MESH POSITION     */
/*************************************************/
using namespace bdry_lvlset;

class tri_hp_lvlset_stype {
    public:
        static const int ntypes = 5;
        enum ids {unknown=-1,inflow,flow_inflow,outflow,characteristic,euler};
        static const char names[ntypes][40];
        static int getid(const char *nin) {
            for(int i=0;i<ntypes;++i)
                if (!strcmp(nin,names[i])) return(i);
            return(-1);
        }
};

const char tri_hp_lvlset_stype::names[ntypes][40] = {"inflow","flow_inflow","outflow","characteristic","euler"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_edge_bdry* tri_hp_lvlset::getnewsideobject(int bnum, input_map &bdrydata) {
    std::string keyword,val;
    std::istringstream data;
    int type;          
    hp_edge_bdry *temp;  
    
        if (bdrydata.get(ebdry(bnum)->idprefix + "_ins_type",val)) {
        type = tri_hp_lvlset_stype::getid(val.c_str());
        if (type == tri_hp_lvlset_stype::unknown)  {
            *gbl->log << "unknown side type:" << val << std::endl;
            exit(1);
        }
    }
    else {
        type = tri_hp_lvlset_stype::unknown;
    }

    switch(type) {
        case tri_hp_lvlset_stype::inflow: {
            temp = new inflow(*this,*ebdry(bnum));
            break;
        }
        case tri_hp_lvlset_stype::flow_inflow: {
#ifdef LOCALIZED_WITH_DISTANCE_FUNCTION
            temp = new flow_inflow(*this,*ebdry(bnum));
#else
            *gbl->log << "can't use flow_inflow type with convective level-set\n";
#endif
            break;
        }
        case tri_hp_lvlset_stype::outflow: {
            temp = new neumann(*this,*ebdry(bnum));
            break;
        }
        case tri_hp_lvlset_stype::characteristic: {
            temp = new characteristic(*this,*ebdry(bnum));
            break;
        }
        case tri_hp_lvlset_stype::euler: {
            temp = new euler(*this,*ebdry(bnum));
            break;
        }
        default: {
            temp = tri_hp_ins::getnewsideobject(bnum,bdrydata);
            break;
        }
    }    
    return(temp);
}

void characteristic::flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, Array<FLT,1>& flx) {
    FLT ul,vl,ur,vr,pl,pr,cl,cr,rho,rhoi;
    FLT s,um,v,c,den,lam0,lam1,lam2,mag,mu;
    FLT nu,gam,qmax;
    Array<FLT,1> ub(x.NV), uvp(x.NV);
    
    /* CHARACTERISTIC FAR-FIELD B.C. */  
    rho = x.gbl->rho + (x.gbl->rho2-x.gbl->rho)*x.heavyside_if(u(2)/x.gbl->width);
    mu = x.gbl->mu + (x.gbl->mu2-x.gbl->mu)*x.heavyside_if(u(2)/x.gbl->width);
    nu = mu/rho;
    rhoi = 1./rho;
    mag = sqrt(norm(0)*norm(0) + norm(1)*norm(1));
    qmax = pow(u(0)-0.5*mv(0),2.0) +pow(u(1)-0.5*mv(1),2.0);
    gam = 3.0*qmax +(0.5*mag*x.gbl->bd[0] +2.*nu/mag)*(0.5*mag*x.gbl->bd[0] +2.*nu/mag);

    norm(0) /= mag;
    norm(1) /= mag;
    
    ul =  u(0)*norm(0) +u(1)*norm(1);
    vl = -u(0)*norm(1) +u(1)*norm(0);
    pl =  u(x.NV-1);
    
    /* FREESTREAM CONDITIONS */
    for(int n=0;n<x.NV;++n)
        ub(n) = x.gbl->ibc->f(n,xpt,x.gbl->time);

    ur =  ub(0)*norm(0) +ub(1)*norm(1);
    vr = -ub(0)*norm(1) +ub(1)*norm(0);
    pr =  ub(x.NV-1);
        
    um = mv(0)*norm(0) +mv(1)*norm(1);
    
    cl = sqrt((ul-.5*um)*(ul-.5*um) +gam);
    cr = sqrt((ur-.5*um)*(ur-.5*um) +gam);
    c = 0.5*(cl+cr);
    s = 0.5*(ul+ur);
    v = 0.5*(vl+vr);
    
    den = 1./(2*c);
    lam0 = s -um;
    lam1 = s-.5*um +c; /* always positive */
    lam2 = s-.5*um -c; /* always negative */
        
    /* PERFORM CHARACTERISTIC SWAP */
    /* BASED ON LINEARIZATION AROUND UL,VL,PL */
    uvp(0) = ((pl-pr)*rhoi +(ul*lam1 -ur*lam2))*den;
    if (lam0 > 0.0) {
        uvp(1) = v*((pr-pl)*rhoi +lam2*(ur-ul))*den/(lam0-lam2) +vl;
        for(int n=tri_mesh::ND;n<x.NV-1;++n)
            uvp(n) = u(n);
    }
    else {
        uvp(1) = v*((pr-pl)*rhoi +lam1*(ur-ul))*den/(lam0-lam1) +vr;
        for(int n=tri_mesh::ND;n<x.NV-1;++n)
            uvp(n) = ub(n);
    }
    uvp(x.NV-1) = (rho*(ul -ur)*gam - lam2*pl +lam1*pr)*den;
    
    /* CHANGE BACK TO X,Y COORDINATES */
    ub(0) =  uvp(0)*norm(0) -uvp(1)*norm(1);
    ub(1) =  uvp(0)*norm(1) +uvp(1)*norm(0);
    
    for(int n=tri_mesh::ND;n<x.NV;++n)
        ub(n) =uvp(n);
            
    norm *= mag;
    
    flx(x.NV-1) = rho*((ub(0) -mv(0))*norm(0) +(ub(1) -mv(1))*norm(1));
    flx(x.NV-2) = flx(x.NV-1)*ub(x.NV-2);

    for(int n=0;n<tri_mesh::ND;++n)
        flx(n) = flx(x.NV-1)*ub(n) +ub(x.NV-1)*norm(n);
        
#ifndef CONSERVATIVE
    flx(x.NV-2) = 0.0;
#endif
        
    return;
}

#ifdef LOCALIZED_WITH_DISTANCE_FUNCTION
void flow_inflow::tadvance() {
    int j,k,m,n,v0,v1,sind,indx,info;
    TinyVector<FLT,tri_mesh::ND> pt;
    char uplo[] = "U";
    
    hp_edge_bdry::tadvance();
        
    /* UPDATE BOUNDARY CONDITION VALUES */
    for(j=0;j<base.nseg;++j) {
        sind = base.seg(j);
        v0 = x.seg(sind).pnt(0);
        for(n=0;n<x.NV-2;++n)
            x.ug.v(v0,n) = x.gbl->ibc->f(n,x.pnts(v0),x.gbl->time);
    }
    v0 = x.seg(sind).pnt(1);
    for(n=0;n<x.NV-2;++n)
        x.ug.v(v0,n) = x.gbl->ibc->f(n,x.pnts(v0),x.gbl->time);

    /*******************/    
    /* SET SIDE VALUES */
    /*******************/
    for(j=0;j<base.nseg;++j) {
        sind = base.seg(j);
        v0 = x.seg(sind).pnt(0);
        v1 = x.seg(sind).pnt(1);
        
        if (is_curved()) {
            x.crdtocht1d(sind);
            for(n=0;n<tri_mesh::ND;++n)
                basis::tri(x.log2p).proj1d(&x.cht(n,0),&x.crd(n)(0,0),&x.dcrd(n,0)(0,0));
        }
        else {
            for(n=0;n<tri_mesh::ND;++n) {
                basis::tri(x.log2p).proj1d(x.pnts(v0)(n),x.pnts(v1)(n),&x.crd(n)(0,0));
                
                for(k=0;k<basis::tri(x.log2p).gpx;++k)
                    x.dcrd(n,0)(0,k) = 0.5*(x.pnts(v1)(n)-x.pnts(v0)(n));
            }
        }

        if (basis::tri(x.log2p).sm) {
            for(n=0;n<x.NV-2;++n)
                basis::tri(x.log2p).proj1d(x.ug.v(v0,n),x.ug.v(v1,n),&x.res(n)(0,0));

            for(k=0;k<basis::tri(x.log2p).gpx; ++k) {
                pt(0) = x.crd(0)(0,k);
                pt(1) = x.crd(1)(0,k);
                for(n=0;n<x.NV-2;++n)
                    x.res(n)(0,k) -= x.gbl->ibc->f(n,pt,x.gbl->time);
            }
            for(n=0;n<x.NV-2;++n)
                basis::tri(x.log2p).intgrt1d(&x.lf(n)(0),&x.res(n)(0,0));

            indx = sind*x.sm0;
            for(n=0;n<x.NV-2;++n) {
                PBTRS(uplo,basis::tri(x.log2p).sm,basis::tri(x.log2p).sbwth,1,&basis::tri(x.log2p).sdiag1d(0,0),basis::tri(x.log2p).sbwth+1,&x.lf(n)(2),basis::tri(x.log2p).sm,info);
                for(m=0;m<basis::tri(x.log2p).sm;++m) 
                    x.ug.s(sind,m,n) = -x.lf(n)(2+m);
            }
        }
    }
    
    return;
}
#endif


