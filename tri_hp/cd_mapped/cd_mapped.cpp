//
//  cd_mapped.cpp
//  tri_hp
//
//  Created by Brian Helenbrook on 6/8/24.
//

#include "cd_mapped.h"

void cd_mapped::init(input_map& input, void *gin) {
    tri_hp_cd::init(input,gin);
    std::string mapval;
    if (!input.get(gbl->idprefix+"_mapping",mapval)) {
        *gbl->log << "Couldn't read mapping " << gbl->idprefix+"_mapping" << std::endl;
        sim::abort(__LINE__,__FILE__,gbl->log);
    }
    
    if (mapval == "polar") {
        map = make_shared<polar_mapping>();
    }
    else if (mapval == "polar_log") {
        map = make_shared<polar_log_mapping>();
    }
    else {
        *gbl->log << "Unrecognized mapping " << mapval << std::endl;
    }
    
    map->init(input,gbl->idprefix,gbl->log);
}

void cd_mapped::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
    tri_hp_cd::init(in,why,sizereduce1d);
    const cd_mapped& cdm = dynamic_cast<const cd_mapped&>(in);
    map = cdm.map;
}

void cd_mapped::element_rsdl(int tind, int stage, Array<TinyVector<FLT,MXTM>,1> &uht,Array<TinyVector<FLT,MXTM>,1> &lf_re,Array<TinyVector<FLT,MXTM>,1> &lf_im){
    FLT visc[ND][ND][ND][ND], tres[NV];
    FLT cv00[MXGP][MXGP],cv01[MXGP][MXGP];
    FLT e00[MXGP][MXGP],e01[MXGP][MXGP];
    const int lgpx = basis::tri(log2p)->gpx(), lgpn = basis::tri(log2p)->gpn();
    
    TinyVector<TinyMatrix<FLT,MXGP,MXGP>,ND> crd, mvel;
    TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,ND,ND> dcrd;
    calc_metrics(tind, crd, dcrd, mvel);
    
    basis::tri(log2p)->proj(&uht(0)(0),&u(0)(0,0),&du(0,0)(0,0),&du(0,1)(0,0),MXGP);
    
    /* lf IS WHERE I WILL STORE THE ELEMENT RESIDUAL */
    for(int n=0;n<NV;++n){
        for(int i=0;i<basis::tri(log2p)->tm();++i){
            lf_re(n)(i) = 0.0;
            lf_im(n)(i) = 0.0;
        }
    }
    
    /* CONVECTION */
    for(int i=0;i<lgpx;++i) {
        for(int j=0;j<lgpn;++j) {
            
#ifdef CONST_A
            FLT fluxx = gbl->rhocv*RAD(crd(0)(i,j))*(gbl->ax -mvel(0)(i,j))*u(0)(i,j);
            FLT fluxy = gbl->rhocv*RAD(crd(0)(i,j))*(gbl->ay -mvel(1)(i,j))*u(0)(i,j);
#else
            pt(0) = crd(0)(i,j);
            pt(1) = crd(1)(i,j);
            FLT fluxx = gbl->rhocv*RAD(crd(0)(i,j))*(gbl->a->f(0,pt,gbl->time) -mvel(0)(i,j))*u(0)(i,j);
            FLT fluxy = gbl->rhocv*RAD(crd(0)(i,j))*(gbl->a->f(1,pt,gbl->time) -mvel(1)(i,j))*u(0)(i,j);
#endif
            cv00[i][j] = +dcrd(1,1)(i,j)*fluxx -dcrd(0,1)(i,j)*fluxy;
            cv01[i][j] = -dcrd(1,0)(i,j)*fluxx +dcrd(0,0)(i,j)*fluxy;
        }
    }
    basis::tri(log2p)->intgrtrs(&lf_im(0)(0),&cv00[0][0],&cv01[0][0],MXGP);
    
    /* NEGATIVE REAL TERMS */
    if (gbl->beta(stage) > 0.0) {
        
        /* TIME DERIVATIVE TERMS */
        for(int i=0;i<lgpx;++i) {
            for(int j=0;j<lgpn;++j) {
                TinyVector<FLT,ND> pt(crd(0)(i,j),crd(1)(i,j));
                cjcb(i,j) = dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j);
                res(0)(i,j) = gbl->rhocv*RAD(crd(0)(i,j))*gbl->bd(0)*u(0)(i,j)*cjcb(i,j) +dugdt(log2p)(tind,0,i,j);
                res(0)(i,j) -= RAD(crd(0)(i,j))*cjcb(i,j)*gbl->src->f(0,pt,gbl->time);
            }
        }
        basis::tri(log2p)->intgrt(&lf_re(0)(0),&res(0)(0,0),MXGP);
        
        /* DIFFUSIVE TERMS  */
        for(int i=0;i<lgpx;++i) {
            for(int j=0;j<lgpn;++j) {
                
                cjcb(i,j) = gbl->kcond*RAD(crd(0)(i,j))/cjcb(i,j);
                
                /* DIFFUSION TENSOR (LOTS OF SYMMETRY THOUGH)*/
                /* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
                visc[0][0][0][0] = -cjcb(i,j)*(dcrd(1,1)(i,j)*dcrd(1,1)(i,j) +dcrd(0,1)(i,j)*dcrd(0,1)(i,j));
                visc[0][0][1][1] = -cjcb(i,j)*(dcrd(1,0)(i,j)*dcrd(1,0)(i,j) +dcrd(0,0)(i,j)*dcrd(0,0)(i,j));
                visc[0][0][0][1] =  cjcb(i,j)*(dcrd(1,1)(i,j)*dcrd(1,0)(i,j) +dcrd(0,1)(i,j)*dcrd(0,0)(i,j));
#define          viscI0II0II1II0I visc[0][0][0][1]
                
                e00[i][j] = +visc[0][0][0][0]*du(0,0)(i,j) +visc[0][0][0][1]*du(0,1)(i,j);
                e01[i][j] = +viscI0II0II1II0I*du(0,0)(i,j) +visc[0][0][1][1]*du(0,1)(i,j);
                
                cv00[i][j] += e00[i][j];
                cv01[i][j] += e01[i][j];
            }
        }
        basis::tri(log2p)->derivr(cv00[0],&res(0)(0,0),MXGP);
        basis::tri(log2p)->derivs(cv01[0],&res(0)(0,0),MXGP);
        
        /* THIS IS BASED ON CONSERVATIVE LINEARIZED MATRICES */
        for(int i=0;i<lgpx;++i) {
            for(int j=0;j<lgpn;++j) {
                TinyVector<FLT,ND> pt(crd(0)(i,j),crd(1)(i,j));
                tres[0] = gbl->tau(tind)*res(0)(i,j);
                
#ifdef CONST_A
                e00[i][j] -= (dcrd(1,1)(i,j)*(gbl->ax-mvel(0)(i,j))
                              -dcrd(0,1)(i,j)*(gbl->ay-mvel(1)(i,j)))*tres[0];
                e01[i][j] -= (-dcrd(1,0)(i,j)*(gbl->ax-mvel(0)(i,j))
                              +dcrd(0,0)(i,j)*(gbl->ay-mvel(1)(i,j)))*tres[0];
#else
                e00[i][j] -= (dcrd(1,1)(i,j)*(gbl->a->f(0,pt,gbl->time)-mvel(0)(i,j))
                              -dcrd(0,1)(i,j)*(gbl->a->f(1,pt,gbl->time)-mvel(1)(i,j)))*tres[0];
                e01[i][j] -= (-dcrd(1,0)(i,j)*(gbl->a->f(0,pt,gbl->time)-mvel(0)(i,j))
                              +dcrd(0,0)(i,j)*(gbl->a->f(1,pt,gbl->time)-mvel(1)(i,j)))*tres[0];
#endif
            }
        }
        basis::tri(log2p)->intgrtrs(&lf_re(0)(0),e00[0],e01[0],MXGP);
        
        for(int n=0;n<NV;++n)
            for(int i=0;i<basis::tri(log2p)->tm();++i)
                lf_re(n)(i) *= gbl->beta(stage);
        
    }
    
    
    return;
}

void cd_mapped::l2error(init_bdry_cndtn *comparison) {
    int i,j,n,tind;
    FLT err;
    Array<int,1> loc(NV);
    Array<FLT,1> mxr(NV),l2r(NV);
    TinyVector<FLT,2> pt;
    
    for(n=0;n<NV;++n) {
        mxr(n) = 0.0;
        l2r(n) = 0.0;
    }
    
    for(tind=0;tind<ntri;++tind) {
        
        TinyVector<TinyMatrix<FLT,MXGP,MXGP>,ND> crd, mvel;
        TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,ND,ND> dcrd;
        calc_metrics(tind, crd, dcrd, mvel);
        
        /* This assumes that the comparison function is written in polar coordinates */
        /* So I am just going to reproject the polar coordinates */
        /* LOAD INDICES OF VERTEX POINTS */
        TinyVector<int,3> v;
        v = tri(tind).pnt;
        
        for(int n=0;n<ND;++n)
            basis::tri(log2p)->proj(pnts(v(0))(n),pnts(v(1))(n),pnts(v(2))(n),&crd(n)(0,0),MXGP);
        
        ugtouht(tind);
        for(n=0;n<NV;++n)
            basis::tri(log2p)->proj(&uht(n)(0),&u(n)(0,0),MXGP);
        
        for (i=0;i<basis::tri(log2p)->gpx();++i) {
            for (j=0;j<basis::tri(log2p)->gpn();++j) {
                cjcb(i,j) = (dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j));
                pt(0) = crd(0)(i,j);
                pt(1) = crd(1)(i,j);
                for(n=0;n<NV;++n) {
                    err =  fabs(u(n)(i,j)-comparison->f(n,pt,gbl->time));
                    if (err >= mxr(n)) {
                        mxr(n) = err;
                        loc(n) = tind;
                    }
                    l2r(n) += err*err*basis::tri(log2p)->wtx(i)*basis::tri(log2p)->wtn(j)*cjcb(i,j);
                }
            }
        }
    }
    
    for(n=0;n<NV;++n) {
        l2r(n) = sqrt(l2r(n));
        *gbl->log << "#L_2: " << l2r(n) << " L_inf " << mxr(n) <<  ' ' << loc(n);
    }
    *gbl->log << '\n';
    
    return;
}

void cd_mapped::calc_metrics(int tind, TinyVector<TinyMatrix<FLT,MXGP,MXGP>,ND>& crd, TinyMatrix<TinyMatrix<FLT,MXGP,MXGP>,ND,ND>& dcrd, TinyVector<TinyMatrix<FLT,MXGP,MXGP>,ND>& mvel) const {
    const int lgpx = basis::tri(log2p)->gpx(), lgpn = basis::tri(log2p)->gpn();
    // for local mesh velocity info
    
    /* LOAD INDICES OF VERTEX POINTS */
    TinyVector<int,3> v;
    v = tri(tind).pnt;
    
    /* PROJECT VERTEX COORDINATES TO GAUSS POINTS */
    for(int n=0;n<ND;++n)
        basis::tri(log2p)->proj(pnts(v(0))(n),pnts(v(1))(n),pnts(v(2))(n),&crd(n)(0,0),MXGP);
    
    /* CALCULATE COORDINATE DERIVATIVES A SIMPLE WAY */
    for(int i=0;i<basis::tri(log2p)->gpx();++i) {
        for(int j=0;j<basis::tri(log2p)->gpn();++j) {
            const TinyVector<FLT,ND> pt(crd(0)(i,j),crd(1)(i,j));
            TinyVector<FLT,ND> xpt;
            
            map->to_physical_frame(pt, xpt);
            crd(0)(i,j) = xpt(0);
            crd(1)(i,j) = xpt(1);
            
            
            TinyMatrix<FLT,ND,ND> dxdrt, drtdxn, jac;
            map->calc_metrics(pt, dxdrt);
            
            for(int n=0;n<ND;++n) {
                drtdxn(n,0) = 0.5*(pnts(v(2))(n) -pnts(v(1))(n));
                drtdxn(n,1) = 0.5*(pnts(v(0))(n) -pnts(v(1))(n));
            }
            
            // dx/dxieta = dx/drtheta*drtheta/dxieta
            // dx/dxi = [dx/dtheta, dx/dr].*[dtheta/dxi, dr/dxi]
            for (int i1 = 0; i1 < ND; ++i1 ) {
                for (int j1 = 0; j1 < ND; ++j1 ) {
                    FLT sum = 0.0;
                    for (int k1 = 0; k1 < ND; ++k1 ) {
                        sum += dxdrt(i1,k1)*drtdxn(k1,j1);
                    }
                    dcrd(i1,j1)(i,j) = sum;
                }
            }
        }
    }
    
    /* CALCULATE MESH VELOCITY */
    for(int i=0;i<lgpx;++i) {
        for(int j=0;j<lgpn;++j) {
            mvel(0)(i,j) = gbl->bd(0)*(crd(0)(i,j) -dxdt(log2p)(tind,0,i,j));
            mvel(1)(i,j) = gbl->bd(0)*(crd(1)(i,j) -dxdt(log2p)(tind,1,i,j));
#ifdef MESH_REF_VEL
            mvel(0)(i,j) += gbl->mesh_ref_vel(0);
            mvel(1)(i,j) += gbl->mesh_ref_vel(1);
#endif
        }
    }
}

void cd_mapped::output(const std::string& fname, block::output_purpose why) {
    tri_hp::output(fname,why);

    if (why == block::display) {
        Array<TinyVector<FLT,ND>,1> mapped_pnts(npnt);
        for (int i = 0; i < npnt; ++i) {
            map->to_physical_frame(pnts(i), mapped_pnts(i));
        }
        Array<TinyVector<FLT,ND>,1> temp;
        temp.reference(vrtxbd(0));
        vrtxbd(0).reference(mapped_pnts);
        auto fname2 = fname +"physical";
        tri_hp::output(fname2,why);
        vrtxbd(0).reference(temp);
    }
}
