#include "bdry_ins.h"
#include <myblas.h>

/*************************************************/
/* SET DIRICHLET BOUNDARY VALUES & FLUXES ********/
/* (THINGS THAT ARE INDEPENDENT OF THE SOLUTION) */
/* BUT NOT INDEPENDENT OF TIME/MESH POSITION     */
/*************************************************/
using namespace bdry_ins;

//void generic::output(std::ostream& fout, tri_hp::filetype typ,int tlvl) {
//	int i,m,n,ind,sind,tind,seg;
//	FLT visc[tri_mesh::ND+1][tri_mesh::ND+1][tri_mesh::ND][tri_mesh::ND];
//    TinyVector<FLT,tri_mesh::ND> norm, mvel;
//    FLT convect;
//    
//    switch(typ) {
//        case(tri_hp::text): case(tri_hp::binary): {
//            hp_edge_bdry::output(fout,typ,tlvl);
//            break;
//        }
//        case(tri_hp::tecplot): {
//            if (!report_flag) return;
//            
//            conv_flux = 0.0;
//            diff_flux = 0.0;
//            moment = 0.0;
//            circumference = 0.0;
//            circulation = 0.0;
//                        
//            for(ind=0; ind < base.nseg; ++ind) {
//                sind = base.seg(ind);
//                tind = x.seg(sind).tri(0);        
//                
//                for(seg=0;seg<3;++seg)
//                    if (x.tri(tind).seg(seg) == sind) break;
//                assert(seg != 3);
//                
//                x.crdtocht(tind);
//                for(m=basis::tri(x.log2p).bm;m<basis::tri(x.log2p).tm;++m)
//                    for(n=0;n<tri_mesh::ND;++n)
//                        x.cht(n,m) = 0.0;
//                        
//                for(n=0;n<tri_mesh::ND;++n)
//                    basis::tri(x.log2p).proj_side(seg,&x.cht(n,0), &x.crd(n)(0,0), &x.dcrd(n,0)(0,0), &x.dcrd(n,1)(0,0));
//
//                x.ugtouht(tind);
//                for(n=0;n<x.NV;++n)
//                    basis::tri(x.log2p).proj_side(seg,&x.uht(n)(0),&x.u(n)(0,0),&x.du(n,0)(0,0),&x.du(n,1)(0,0));
//
//                for (i=0;i<basis::tri(x.log2p).gpx;++i) {
//                    circumference += basis::tri(x.log2p).wtx(i)*sqrt(x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i) +x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i));
//                    x.cjcb(0,i) = x.gbl->mu*RAD(x.crd(0)(0,i))/(x.dcrd(0,0)(0,i)*x.dcrd(1,1)(0,i) -x.dcrd(1,0)(0,i)*x.dcrd(0,1)(0,i));
//                    
//                    /* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
//                    /* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
//                    visc[0][0][0][0] =  x.cjcb(0,i)*(2.*x.dcrd(1,1)(0,i)*x.dcrd(1,1)(0,i) +x.dcrd(0,1)(0,i)*x.dcrd(0,1)(0,i));
//                    visc[0][0][1][1] =  x.cjcb(0,i)*(2.*x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i));
//                    visc[0][0][0][1] = -x.cjcb(0,i)*(2.*x.dcrd(1,1)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,1)(0,i)*x.dcrd(0,0)(0,i));
//#define             viscI0II0II1II0I visc[0][0][0][1]
//
//                    visc[1][1][0][0] =  x.cjcb(0,i)*(x.dcrd(1,1)(0,i)*x.dcrd(1,1)(0,i) +2.*x.dcrd(0,1)(0,i)*x.dcrd(0,1)(0,i));
//                    visc[1][1][1][1] =  x.cjcb(0,i)*(x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i) +2.*x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i));
//                    visc[1][1][0][1] = -x.cjcb(0,i)*(x.dcrd(1,1)(0,i)*x.dcrd(1,0)(0,i) +2.*x.dcrd(0,1)(0,i)*x.dcrd(0,0)(0,i));
//#define             viscI1II1II1II0I visc[1][1][0][1]
//                    
//                    visc[0][1][0][0] = -x.cjcb(0,i)*x.dcrd(0,1)(0,i)*x.dcrd(1,1)(0,i);
//                    visc[0][1][1][1] = -x.cjcb(0,i)*x.dcrd(0,0)(0,i)*x.dcrd(1,0)(0,i);
//                    visc[0][1][0][1] =  x.cjcb(0,i)*x.dcrd(0,1)(0,i)*x.dcrd(1,0)(0,i);
//                    visc[0][1][1][0] =  x.cjcb(0,i)*x.dcrd(0,0)(0,i)*x.dcrd(1,1)(0,i);
//
//                    /* OTHER SYMMETRIES     */                
//#define             viscI1II0II0II0I visc[0][1][0][0]
//#define             viscI1II0II1II1I visc[0][1][1][1]
//#define             viscI1II0II0II1I visc[0][1][1][0]
//#define             viscI1II0II1II0I visc[0][1][0][1]
//
//                    /* DIFFUSIVE FLUXES ( FOR EXTRA VARIABLES) */
//                    visc[2][2][1][0] =  x.cjcb(0,i)*(x.dcrd(1,1)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,1)(0,i)*x.dcrd(0,0)(0,i));
//                    visc[2][2][1][1] = -x.cjcb(0,i)*(x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i));
//
//                    for (n=tri_mesh::ND;n<x.NV-1;++n) 
//                        diff_flux(n) -= x.gbl->D(n)/x.gbl->mu*basis::tri(x.log2p).wtx(i)*(-visc[2][2][1][0]*x.du(2,0)(0,i) -visc[2][2][1][1]*x.du(2,1)(0,i));
//
//                    diff_flux(0) -=    basis::tri(x.log2p).wtx(i)*(-x.u(2)(0,i)*RAD(x.crd(0)(0,i))*x.dcrd(1,0)(0,i) 
//                                    -viscI0II0II1II0I*x.du(0,0)(0,i) -visc[0][1][1][0]*x.du(1,0)(0,i)
//                                    -visc[0][0][1][1]*x.du(0,1)(0,i) -visc[0][1][1][1]*x.du(1,1)(0,i));															
//                    diff_flux(1) -=    basis::tri(x.log2p).wtx(i)*( x.u(2)(0,i)*RAD(x.crd(0)(0,i))*x.dcrd(0,0)(0,i)
//                                    -viscI1II0II1II0I*x.du(0,0)(0,i) -viscI1II1II1II0I*x.du(1,0)(0,i)
//                                    -viscI1II0II1II1I*x.du(0,1)(0,i) -visc[1][1][1][1]*x.du(1,1)(0,i));
//                                    
//                                    
//                    norm(0) = x.dcrd(1,0)(0,i);
//                    norm(1) = -x.dcrd(0,0)(0,i);                
//                    for(n=0;n<tri_mesh::ND;++n) {
//                        mvel(n) = x.gbl->bd(0)*(x.crd(n)(0,i) -dxdt(x.log2p,ind)(n,i));
//#ifdef DROP
//                        mvel(n) += tri_hp_ins::mesh_ref_vel(n);
//#endif
//                    }
//                    
//                    circulation += -norm(1)*(x.u(0)(0,i)-mvel(0)) +norm(0)*(x.u(1)(0,i)-mvel(1));
//                    
//                    convect = basis::tri(x.log2p).wtx(i)*RAD(x.crd(0)(0,i))*((x.u(0)(0,i)-mvel(0))*norm(0) +(x.u(1)(0,i)-mvel(1))*norm(1));
//                    conv_flux(2) -= convect;
//                    conv_flux(0) -= x.u(0)(0,i)*convect;
//                    conv_flux(1) -= x.u(1)(0,i)*convect;
//                }				
//            }
//            fout << base.idprefix << " circumference: " << circumference << std::endl;
//            fout << base.idprefix << " viscous/pressure flux: " << diff_flux << std::endl;
//            fout << base.idprefix << " convective flux: " << conv_flux << std::endl;
//            fout << base.idprefix << " circulation: " << circulation << std::endl;
//                
//            /* OUTPUT AUXILIARY FLUXES */
//            fout << base.idprefix << "total fluxes: " << total_flux << std::endl;
//            break;
//        }
//                
//        case(tri_hp::auxiliary): {
//            if (!report_flag) return;                
//            
//            /* AUXILIARY FLUX METHOD */
//            int v0;
//            total_flux = 0.0;
//            for(ind=0; ind < base.nseg; ++ind) {
//                sind = base.seg(ind);
//                v0 = x.seg(sind).pnt(0);
//                total_flux += x.gbl->res.v(v0,Range::all());
//            }
//            v0 = x.seg(sind).pnt(1);
//            total_flux += x.gbl->res.v(v0,Range::all());
//        }
//    }
//    
//	return;
//}




//void neumann::rsdl(int stage) {
//    int j,k,n,v0,v1,sind;
//    TinyVector<FLT,2> pt,mvel,nrm;
//    Array<FLT,1> u(x.NV),flx(x.NV);
//    
//    for(j=0;j<base.nseg;++j) {
//        sind = base.seg(j);
//        v0 = x.seg(sind).pnt(0);
//        v1 = x.seg(sind).pnt(1);
//        
//        x.crdtocht1d(sind);
//        for(n=0;n<tri_mesh::ND;++n)
//            basis::tri(x.log2p).proj1d(&x.cht(n,0),&x.crd(n)(0,0),&x.dcrd(n,0)(0,0));
//        
//        x.crdtocht1d(sind,1);
//        for(n=0;n<tri_mesh::ND;++n)
//            basis::tri(x.log2p).proj1d(&x.cht(n,0),&x.crd(n)(1,0));
//        
//        x.ugtouht1d(sind);
//        for(n=0;n<x.NV;++n)
//            basis::tri(x.log2p).proj1d(&x.uht(n)(0),&x.u(n)(0,0));
//
//        for(k=0;k<basis::tri(x.log2p).gpx;++k) {
//            nrm(0) = x.dcrd(1,0)(0,k);
//            nrm(1) = -x.dcrd(0,0)(0,k);                
//            for(n=0;n<tri_mesh::ND;++n) {
//                pt(n) = x.crd(n)(0,k);
//                mvel(n) = x.gbl->bd(0)*(x.crd(n)(0,k) -dxdt(x.log2p,j)(n,k));
//#ifdef DROP
//                mvel(n) += tri_hp_ins::mesh_ref_vel(n);
//#endif
//            }
//            
//            for(n=0;n<x.NV;++n)
//                u(n) = x.u(n)(0,k);
//            
//            flux(u,pt,mvel,nrm,flx);
//                        
//            for(n=0;n<x.NV;++n)
//                x.res(n)(0,k) = RAD(x.crd(0)(0,k))*flx(n);
//
//        }
//
//        for(n=0;n<x.NV;++n)
//            basis::tri(x.log2p).intgrt1d(&x.lf(n)(0),&x.res(n)(0,0));
//                    
//        for(n=0;n<x.NV;++n)
//            x.gbl->res.v(v0,n) += x.lf(n)(0);
//
//        for(n=0;n<x.NV;++n)
//            x.gbl->res.v(v1,n) += x.lf(n)(1);
//        
//        for(k=0;k<basis::tri(x.log2p).sm;++k) {
//            for(n=0;n<x.NV;++n)
//                x.gbl->res.s(sind,k,n) += x.lf(n)(k+2);
//        }
//    }
//    return;
//}

void applied_stress::init(input_map& inmap,void* gbl_in) {
    std::string keyword;
    std::ostringstream nstr;

    neumann::init(inmap,gbl_in);
    
    stress.resize(tet_mesh::ND);

    for(int n=0;n<tet_mesh::ND;++n) {
        nstr.str("");
        nstr << base.idprefix << "_stress" << n << std::flush;
        if (inmap.find(nstr.str()) != inmap.end()) {
            stress(n).init(inmap,nstr.str());
        }
        else {
            *x.gbl->log << "couldn't find stress function " << nstr.str() << '\n';
            exit(1);
        }
    }
    
    return;
}

//void inflow::setvalues(init_bdry_cndtn *ibc) {
//    int j,k,m,n,v0,v1,sind,indx,info;
//	TinyVector<FLT,tri_mesh::ND> pt;
//    char uplo[] = "U";
//
//    /* UPDATE BOUNDARY CONDITION VALUES */
//    for(j=0;j<base.nseg;++j) {
//        sind = base.seg(j);
//        v0 = x.seg(sind).pnt(0);
//        for(n=0;n<x.NV-1;++n)
//            x.ug.v(v0,n) = ibc->f(n,x.pnts(v0),x.gbl->time);
//    }
//    v0 = x.seg(sind).pnt(1);
//    for(n=0;n<x.NV-1;++n)
//        x.ug.v(v0,n) = ibc->f(n,x.pnts(v0),x.gbl->time);
//
//    /*******************/    
//    /* SET SIDE VALUES */
//    /*******************/
//    for(j=0;j<base.nseg;++j) {
//        sind = base.seg(j);
//        v0 = x.seg(sind).pnt(0);
//        v1 = x.seg(sind).pnt(1);
//        
//        if (is_curved()) {
//            x.crdtocht1d(sind);
//            for(n=0;n<tri_mesh::ND;++n)
//                basis::tri(x.log2p).proj1d(&x.cht(n,0),&x.crd(n)(0,0),&x.dcrd(n,0)(0,0));
//        }
//        else {
//            for(n=0;n<tri_mesh::ND;++n) {
//                basis::tri(x.log2p).proj1d(x.pnts(v0)(n),x.pnts(v1)(n),&x.crd(n)(0,0));
//                
//                for(k=0;k<basis::tri(x.log2p).gpx;++k)
//                    x.dcrd(n,0)(0,k) = 0.5*(x.pnts(v1)(n)-x.pnts(v0)(n));
//            }
//        }
//
//        if (basis::tri(x.log2p).sm) {
//            for(n=0;n<x.NV-1;++n)
//                basis::tri(x.log2p).proj1d(x.ug.v(v0,n),x.ug.v(v1,n),&x.res(n)(0,0));
//
//            for(k=0;k<basis::tri(x.log2p).gpx; ++k) {
//                pt(0) = x.crd(0)(0,k);
//                pt(1) = x.crd(1)(0,k);
//                for(n=0;n<x.NV-1;++n)
//                    x.res(n)(0,k) -= ibc->f(n,pt,x.gbl->time);
//            }
//            for(n=0;n<x.NV-1;++n)
//                basis::tri(x.log2p).intgrt1d(&x.lf(n)(0),&x.res(n)(0,0));
//
//            indx = sind*x.sm0;
//            for(n=0;n<x.NV-1;++n) {
//                PBTRS(uplo,basis::tri(x.log2p).sm,basis::tri(x.log2p).sbwth,1,&basis::tri(x.log2p).sdiag1d(0,0),basis::tri(x.log2p).sbwth+1,&x.lf(n)(2),basis::tri(x.log2p).sm,info);
//                for(m=0;m<basis::tri(x.log2p).sm;++m) 
//                    x.ug.s(sind,m,n) = -x.lf(n)(2+m);
//            }
//        }
//    }
//    return;
//}

//void symmetry::tadvance() {
//    int j,m,v0,sind;
//    TinyVector<FLT,tri_mesh::ND> pt;
//    
//    hp_edge_bdry::tadvance();
//    
//    /* UPDATE BOUNDARY CONDITION VALUES */
//    for(j=0;j<base.nseg;++j) {
//        sind = base.seg(j);
//        v0 = x.seg(sind).pnt(0);
//        x.ug.v(v0,dir) = 0.0;
//    }
//    v0 = x.seg(sind).pnt(1);
//    x.ug.v(v0,dir) = 0.0;
//
//    /*******************/    
//    /* SET SIDE VALUES */
//    /*******************/
//    for(j=0;j<base.nseg;++j) {
//        sind = base.seg(j);
//        for(m=0;m<basis::tri(x.log2p).sm;++m) {
//            x.ug.s(sind,m,dir) = 0.0;
//        }
//    }
//    
//    return;
//}


//void characteristic::flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, Array<FLT,1>& flx) {
//    FLT ul,vl,ur,vr,pl,pr,cl,cr,rho,rhoi;
//    FLT s,um,v,c,den,lam0,lam1,lam2,mag,hmax;
//    FLT nu,gam,qmax;
//    Array<FLT,1> ub(x.NV), uvp(x.NV);
//    
//    /* CHARACTERISTIC FAR-FIELD B.C. */   
//	
//    rho = x.gbl->rho;
//    nu = x.gbl->mu/x.gbl->rho;
//    rhoi = 1./rho;
//    mag = sqrt(norm(0)*norm(0) + norm(1)*norm(1));
//	hmax = mag*2.0/(0.25*(basis::tri(x.log2p).p +1)*(basis::tri(x.log2p).p+1));
//    qmax = pow(u(0)-0.5*mv(0),2.0) +pow(u(1)-0.5*mv(1),2.0);
//	gam = 3.0*qmax +(0.5*hmax*x.gbl->bd(0) +2.*nu/hmax)*(0.5*hmax*x.gbl->bd(0) +2.*nu/hmax);
//
//    norm(0) /= mag;
//    norm(1) /= mag;
//    
//    ul =  u(0)*norm(0) +u(1)*norm(1);
//    vl = -u(0)*norm(1) +u(1)*norm(0);
//    pl =  u(x.NV-1);
//    
//    /* FREESTREAM CONDITIONS */
//    for(int n=0;n<x.NV;++n)
//        ub(n) = ibc->f(n,xpt,x.gbl->time);
//        
//    ur =  ub(0)*norm(0) +ub(1)*norm(1);
//    vr = -ub(0)*norm(1) +ub(1)*norm(0);
//    pr =  ub(x.NV-1);
//        
//    um = mv(0)*norm(0) +mv(1)*norm(1);
//    
//    cl = sqrt((ul-.5*um)*(ul-.5*um) +gam);
//    cr = sqrt((ur-.5*um)*(ur-.5*um) +gam);
//    c = 0.5*(cl+cr);
//    s = 0.5*(ul+ur);
//    v = 0.5*(vl+vr);
//    
//    den = 1./(2*c);
//    lam0 = s -um;
//    lam1 = s-.5*um +c; /* always positive */
//    lam2 = s-.5*um -c; /* always negative */
//        
//    /* PERFORM CHARACTERISTIC SWAP */
//    /* BASED ON LINEARIZATION AROUND UL,VL,PL */
//    uvp(0) = ((pl-pr)*rhoi +(ul*lam1 -ur*lam2))*den;
//    if (lam0 > 0.0) {
//        uvp(1) = v*((pr-pl)*rhoi +lam2*(ur-ul))*den/(lam0-lam2) +vl;
//        for(int n=tri_mesh::ND;n<x.NV-1;++n)
//            uvp(n) = u(n);
//    }
//    else {
//        uvp(1) = v*((pr-pl)*rhoi +lam1*(ur-ul))*den/(lam0-lam1) +vr;
//        for(int n=tri_mesh::ND;n<x.NV-1;++n)
//            uvp(n) = ub(n);
//    }
//    uvp(x.NV-1) = (rho*(ul -ur)*gam - lam2*pl +lam1*pr)*den;
//    
//    /* CHANGE BACK TO X,Y COORDINATES */
//    ub(0) =  uvp(0)*norm(0) -uvp(1)*norm(1);
//    ub(1) =  uvp(0)*norm(1) +uvp(1)*norm(0);
//    
//    for(int n=tri_mesh::ND;n<x.NV;++n)
//        ub(n) =uvp(n);
//    
//    norm *= mag;
//    
//    flx(x.NV-1) = rho*((ub(0) -mv(0))*norm(0) +(ub(1) -mv(1))*norm(1));
//    
//    for(int n=0;n<tri_mesh::ND;++n)
//        flx(n) = flx(x.NV-1)*ub(n) +ub(x.NV-1)*norm(n);
//        
//    for(int n=tri_mesh::ND;n<x.NV-1;++n)
//        flx(n) = flx(x.NV-1)*ub(n);
//		
//	// *x.gbl->log << x.npnt << '\t' << u << '\t' << xpt << '\t' << mv << '\t' << norm << '\t' << flx << '\n';
//
//        
//    return;
//}