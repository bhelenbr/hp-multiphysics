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
hp_side_bdry* tri_hp_lvlset::getnewsideobject(int bnum, input_map &bdrydata) {
    std::string keyword,val;
    std::istringstream data;
    int type;          
    hp_side_bdry *temp;  
    
    
    keyword =  sbdry(bnum)->idprefix + "_ins_type";
    if (bdrydata.get(keyword,val)) {
        type = tri_hp_lvlset_stype::getid(val.c_str());
        if (type == tri_hp_lvlset_stype::unknown)  {
            *sim::log << "unknown side type:" << val << std::endl;
            exit(1);
        }
    }
    else {
        type = tri_hp_lvlset_stype::unknown;
    }

    switch(type) {
        case tri_hp_lvlset_stype::inflow: {
            temp = new inflow(*this,*sbdry(bnum));
            break;
        }
        case tri_hp_lvlset_stype::flow_inflow: {
            temp = new flow_inflow(*this,*sbdry(bnum));
            break;
        }
        case tri_hp_lvlset_stype::outflow: {
            temp = new neumann(*this,*sbdry(bnum));
            break;
        }
        case tri_hp_lvlset_stype::characteristic: {
            temp = new characteristic(*this,*sbdry(bnum));
            break;
        }
        case tri_hp_lvlset_stype::euler: {
            temp = new euler(*this,*sbdry(bnum));
            break;
        }
        default: {
            temp = tri_hp_ins::getnewsideobject(bnum,bdrydata);
            break;
        }
    }    
    return(temp);
}

void characteristic::flux(Array<FLT,1>& u, TinyVector<FLT,mesh::ND> xpt, TinyVector<FLT,mesh::ND> mv, TinyVector<FLT,mesh::ND> norm, Array<FLT,1>& flx) {
    FLT ul,vl,ur,vr,pl,pr,cl,cr,rho,rhoi;
    FLT s,um,v,c,den,lam0,lam1,lam2,mag,mu;
    FLT nu,gam,qmax;
    Array<FLT,1> ub(x.NV), uvp(x.NV);
    
    /* CHARACTERISTIC FAR-FIELD B.C. */  
    rho = x.gbl_ptr->rho + (x.gbl_ptr->rho2-x.gbl_ptr->rho)*x.heavyside_if(u(2)/x.gbl_ptr->width);
    mu = x.gbl_ptr->mu + (x.gbl_ptr->mu2-x.gbl_ptr->mu)*x.heavyside_if(u(2)/x.gbl_ptr->width);
    nu = mu/rho;
    rhoi = 1./rho;
    mag = sqrt(norm(0)*norm(0) + norm(1)*norm(1));
    qmax = pow(u(0)-0.5*mv(0),2.0) +pow(u(1)-0.5*mv(1),2.0);
    gam = 3.0*qmax +(0.5*mag*sim::bd[0] +2.*nu/mag)*(0.5*mag*sim::bd[0] +2.*nu/mag);

    norm(0) /= mag;
    norm(1) /= mag;
    
    ul =  u(0)*norm(0) +u(1)*norm(1);
    vl = -u(0)*norm(1) +u(1)*norm(0);
    pl =  u(x.NV-1);
    
    /* FREESTREAM CONDITIONS */
    for(int n=0;n<x.NV;++n)
        ub(n) = x.gbl_ptr->ibc->f(n,xpt);

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
        for(int n=mesh::ND;n<x.NV-1;++n)
            uvp(n) = u(n);
    }
    else {
        uvp(1) = v*((pr-pl)*rhoi +lam1*(ur-ul))*den/(lam0-lam1) +vr;
        for(int n=mesh::ND;n<x.NV-1;++n)
            uvp(n) = ub(n);
    }
    uvp(x.NV-1) = (rho*(ul -ur)*gam - lam2*pl +lam1*pr)*den;
    
    /* CHANGE BACK TO X,Y COORDINATES */
    ub(0) =  uvp(0)*norm(0) -uvp(1)*norm(1);
    ub(1) =  uvp(0)*norm(1) +uvp(1)*norm(0);
    
    for(int n=mesh::ND;n<x.NV-1;++n)
        ub(n) =uvp(n);
    
    ub(x.NV-1) =  uvp(x.NV-1);
    
    norm *= mag;
    
    flx(x.NV-1) = rho*((ub(0) -mv(0))*norm(0) +(ub(1) -mv(1))*norm(1));
    
    for(int n=0;n<mesh::ND;++n)
        flx(n) = flx(x.NV-1)*ub(n) +ub(x.NV-1)*norm(n);
    flx(2) = 0.0;
        
    return;
}

block::ctrl flow_inflow::tadvance(bool coarse, block::ctrl ctrl_message) {
    int j,k,m,n,v0,v1,sind,indx,info;
    TinyVector<FLT,mesh::ND> pt;
    char uplo[] = "U";
    block::ctrl state;
    
    if (ctrl_message == block::begin) excpt1 = 0;
    
    if (excpt1 == 0) {
        if ((state = hp_side_bdry::tadvance(coarse,ctrl_message)) != block::stop) return(state);
        ++excpt1;
        
        /* UPDATE BOUNDARY CONDITION VALUES */
        for(j=0;j<base.nel;++j) {
            sind = base.el(j);
            v0 = x.sd(sind).vrtx(0);
            for(n=0;n<x.NV-2;++n)
                x.ug.v(v0,n) = x.gbl_ptr->ibc->f(n,x.vrtx(v0));
        }
        v0 = x.sd(sind).vrtx(1);
        for(n=0;n<x.NV-2;++n)
            x.ug.v(v0,n) = x.gbl_ptr->ibc->f(n,x.vrtx(v0));

        /*******************/    
        /* SET SIDE VALUES */
        /*******************/
        for(j=0;j<base.nel;++j) {
            sind = base.el(j);
            v0 = x.sd(sind).vrtx(0);
            v1 = x.sd(sind).vrtx(1);
            
            if (is_curved()) {
                x.crdtocht1d(sind);
                for(n=0;n<mesh::ND;++n)
                    basis::tri(x.log2p).proj1d(&x.cht(n,0),&x.crd(n)(0,0),&x.dcrd(n,0)(0,0));
            }
            else {
                for(n=0;n<mesh::ND;++n) {
                    basis::tri(x.log2p).proj1d(x.vrtx(v0)(n),x.vrtx(v1)(n),&x.crd(n)(0,0));
                    
                    for(k=0;k<basis::tri(x.log2p).gpx;++k)
                        x.dcrd(n,0)(0,k) = 0.5*(x.vrtx(v1)(n)-x.vrtx(v0)(n));
                }
            }

            if (basis::tri(x.log2p).sm) {
                for(n=0;n<x.NV-2;++n)
                    basis::tri(x.log2p).proj1d(x.ug.v(v0,n),x.ug.v(v1,n),&x.res(n)(0,0));

                for(k=0;k<basis::tri(x.log2p).gpx; ++k) {
                    pt(0) = x.crd(0)(0,k);
                    pt(1) = x.crd(1)(0,k);
                    for(n=0;n<x.NV-2;++n)
                        x.res(n)(0,k) -= x.gbl_ptr->ibc->f(n,pt);
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
    }
    return(block::stop);
}



