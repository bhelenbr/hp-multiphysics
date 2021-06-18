//
//  shock.cpp
//  tri_hp
//
//  Created by Brian Helenbrook on 12/30/14.
//
//

#include "shock.h"
#include "bdry_cns.h"
#include <myblas.h>
#include <iostream>

#define NEW_STABILIZATION
//#define WAY1

using namespace bdry_cns;

// extern FLT body[ND];


void shock::init(input_map& inmap,void* gin) {
	std::string keyword,matching_block,side_id,master_block,master_id;
	std::istringstream data;
	std::string filename;

	const int sm = basis::tri(x.log2p)->sm();
    gbl = static_cast<global *>(gin);

    /* Both sides have to do work for this Boundary Condition */
    inmap[base.idprefix+"_symmetric"] = "1";
    
	hp_coupled_bdry::init(inmap,gin);
#ifdef WAY1
    c0_indices_xy.clear();
#endif
    
	u_opp_v.resize(base.nseg+1,x.NV);
	u_opp_e.resize(base.nseg,sm,x.NV);
	
	return;
}


void shock::rsdl(int stage) {
    
	const int sm = basis::tri(x.log2p)->sm();
    
    base.vloadbuff(boundary::all,x.ug.v.data(),0,x.NV-1,x.NV);
	base.comm_prepare(boundary::all,0,boundary::symmetric);
	base.comm_exchange(boundary::all,0,boundary::symmetric);
	base.comm_wait(boundary::all,0,boundary::symmetric);
	base.comm_finish(boundary::all,0,boundary::symmetric,boundary::replace);
	
    if (base.is_frst()) {
        for (int i=0,count=0; i<base.nseg+1; i++){
            for(int n=0; n<x.NV; n++){
                u_opp_v(i,n) = base.fsndbuf(count++);
            }
        }
    }
    else {
        for (int i=base.nseg,count=0; i>=0; i--){
            for(int n=0; n<x.NV; n++){
                u_opp_v(i,n) = base.fsndbuf(count++);
            }
        }
    }
    
    
    int j,k,count,offset,sind;
    int bgn = 0;
    int end = sm-1;
    int stride = sm;
    FLT *sdata = x.ug.s.data();

    int ebp1 = end-bgn+1;
    count = 0;
    for(j=0;j<base.nseg;++j) {
        sind = base.seg(j);
        offset = (sind*stride +bgn)*x.NV;
        FLT sign = 1.0;
        for(k=0;k<ebp1;++k) {
            for(int n = 0; n<x.NV; n++){
                base.fsndbuf(count++) = sign*sdata[offset +n];
            }
            offset+=x.NV;
//            sign*=-1.0;
        }
    }

    base.sndsize() = count;
    base.sndtype() = boundary::flt_msg;


    base.comm_prepare(boundary::all,0,boundary::symmetric);
    base.comm_exchange(boundary::all,0,boundary::symmetric);
    base.comm_wait(boundary::all,0,boundary::symmetric);
    base.comm_finish(boundary::all,0,boundary::symmetric,boundary::replace);


    for (int j=base.nseg-1,count=0;j>=0; j--){
        for (int k=0;k<sm;k++){
            for (int n=0;n<x.NV;n++){
                u_opp_e(j,k,n) = base.fsndbuf(count++);
            }
        }
    }
    
	hp_coupled_bdry::rsdl(stage);

}

void shock::element_rsdl(int indx, Array<TinyVector<FLT,MXTM>,1> lf) {

	const int sm = basis::tri(x.log2p)->sm();

	int i,n,sind;
	TinyVector<FLT,tri_mesh::ND> norm, rp, norm_meshc;
	Array<FLT,1> ubar(x.NV);
	FLT jcb; //jacobian 1d on edge
	Array<TinyVector<FLT,MXGP>,1> u(x.NV),u_opp(x.NV);
	TinyMatrix<FLT,tri_mesh::ND,MXGP> crd, dcrd, mvel_u, mvel_d;
	TinyMatrix<FLT,7,MXGP> res;
	Array<TinyVector<FLT,MXTM>,1> uht_opp(x.NV);

	for (int n=0; n<x.NV; n++){
		uht_opp(n)(0) = u_opp_v(indx,n);
		uht_opp(n)(1) = u_opp_v(indx+1,n);
        FLT sign = 1.0;
		for (int j=0; j<sm; j++){
			uht_opp(n)(j+2) = sign*u_opp_e(indx,j,n);
            sign *= -1.0;
		}
	}

	sind = base.seg(indx);
	x.crdtocht1d(sind);
	for(n=0;n<tri_mesh::ND;++n)
		basis::tri(x.log2p)->proj1d(&x.cht(n,0),&crd(n,0),&dcrd(n,0));
	
	for(n=0;n<x.NV;++n)
		basis::tri(x.log2p)->proj1d(&x.uht(n)(0),&u(n)(0));
    
	for(n=0;n<x.NV;++n)
		basis::tri(x.log2p)->proj1d(&uht_opp(n)(0),&u_opp(n)(0));
    
#ifdef NEW_STABILIZATION
    /* Calculate stabilization constant based on analysis of linear elements and constant tau */
    int v0 = x.seg(sind).pnt(0);
    int v1 = x.seg(sind).pnt(1);
    norm(0) =  (x.pnts(v1)(1) -x.pnts(v0)(1));
    norm(1) = -(x.pnts(v1)(0) -x.pnts(v0)(0));
    FLT h = sqrt(norm(0)*norm(0) +norm(1)*norm(1));
    FLT hsm = h/(.25*(basis::tri(x.log2p)->p()+1)*(basis::tri(x.log2p)->p()+1));
    TinyVector<FLT,tri_mesh::ND> mvel;
    mvel(0) = x.uht(0)(0)-(x.gbl->bd(0)*(x.pnts(v0)(0) -x.vrtxbd(1)(v0)(0)));
    mvel(1) = x.uht(1)(0)-(x.gbl->bd(0)*(x.pnts(v0)(1) -x.vrtxbd(1)(v0)(1)));
#ifdef MESH_REF_VEL
    mvel(0) -= x.gbl->mesh_ref_vel(0);
    mvel(1) -= x.gbl->mesh_ref_vel(1);
#endif  
    FLT vslp0 = (-mvel(0)*norm(1) +mvel(1)*norm(0))/h;

    mvel(0) = x.uht(0)(1)-(x.gbl->bd(0)*(x.pnts(v1)(0) -x.vrtxbd(1)(v1)(0)));
    mvel(1) = x.uht(1)(1)-(x.gbl->bd(0)*(x.pnts(v1)(1) -x.vrtxbd(1)(v1)(1)));
#ifdef MESH_REF_VEL
    mvel(0) -= x.gbl->mesh_ref_vel(0);
    mvel(1) -= x.gbl->mesh_ref_vel(1);
#endif
    FLT vslp1 = (-mvel(0)*norm(1) +mvel(1)*norm(0))/h;
    
    FLT meshc = gbl->adis*hsm*(3*(abs(vslp0)+abs(vslp1)) +vslp0-vslp1)/(4*(vslp0*vslp0+vslp0*vslp1+vslp1*vslp1)+100*FLT_EPSILON)*2/h;
    FLT meshc1 = gbl->adis*hsm*(3*(abs(vslp0)+abs(vslp1)) -vslp0+vslp1)/(4*(vslp0*vslp0+vslp0*vslp1+vslp1*vslp1)+100*FLT_EPSILON)*2/h;
#endif
    
	for(i=0;i<basis::tri(x.log2p)->gpx();++i) {
		norm(0) =  dcrd(1,i);
		norm(1) = -dcrd(0,i);
		jcb = sqrt(norm(0)*norm(0) +norm(1)*norm(1));
        
        TinyVector<FLT,2> urel, urel_opp, mvel;
        for (int n=0;n<x.ND;++n) {
            mvel(n) = x.gbl->bd(0)*(crd(n,i) -dxdt(x.log2p,indx)(n,i));
#ifdef MESH_REF_VEL
            mvel(n) -= x.gbl->mesh_ref_vel(n);
#endif
            urel(n) = u(n+1)(i) -mvel(n);
            urel_opp(n) = u_opp(n+1)(i)-mvel(n);
        }
        
        FLT mvel_norm = (mvel(0)*norm(0) +mvel(1)*norm(1))/jcb;
        FLT unorm = (u(1)(i)*norm(0) +u(2)(i)*norm(1))/jcb;
        FLT unorm_rel = unorm-mvel_norm;

        FLT unorm_opp, pu, uu, vu, RTu, cd, vslpu, sign;
        if (unorm_rel > 0.0) {
            /* This is upstream side */
            sign = 1.0;
            unorm_opp = -(u_opp(1)(i)*norm(0) +u(2)(i)*norm(1))/jcb;
            pu = u(0)(i);
            uu = u(1)(i);
            vu = u(2)(i);
            RTu = u(3)(i);
            cd = sqrt(x.gbl->gamma*u_opp(3)(i));
            vslpu = (-urel(0)*norm(1) +urel(1)*norm(0))/jcb;
        } else {
            /* This is downstream side */
            sign = -1;
            unorm_opp = unorm;
            unorm = -(u_opp(1)(i)*norm(0) +u(2)(i)*norm(1))/jcb;
            mvel_norm *= -1.0;
            unorm_rel = unorm-mvel_norm;
            pu = u_opp(0)(i);
            uu = u_opp(1)(i);
            vu = u_opp(2)(i);
            RTu = u_opp(3)(i);
            cd = sqrt(x.gbl->gamma*u(3)(i));
            vslpu = -(-urel_opp(0)*norm(1) +urel_opp(1)*norm(0))/jcb;
#ifdef NEW_STABILIZATION
            meshc = meshc1;
#endif
        }
        FLT rhou = pu/RTu;
        FLT Eu = RTu/(x.gbl->gamma-1.0)+0.5*(uu*uu +vu*vu);
        FLT cu = sqrt(x.gbl->gamma*RTu);
        FLT Mu = unorm_rel/cu;
        if (int err = shock_mach(Mu, cu, cd, unorm+unorm_opp)) {
            *x.gbl->log << "#Trouble finding Mach number\n";
        }
        
        /* TANGENTIAL SPACING */
        res(0,i) = -sign*ksprg(indx)*jcb;
        
        FLT shock_vel = unorm-cu*Mu;
        res(1,i) = (mvel_norm-shock_vel)*jcb;
        
//        (*x.gbl->log).precision(16);
//        *x.gbl->log << Mu << ' ' << pu << ' ' << RTu << ' ' << rhou << ' ' << uu << ' ' << vu << ' ' << cu << ' ' << cd << ' ' << mvel_norm << ' ' << unorm << ' ' << shock_vel << std::endl;

        
        /* UPWINDING BASED ON TANGENTIAL VELOCITY */
#ifndef NEW_STABILIZATION
        // gbl->meshc(indx) = gbl->adis/((basis::tri(x.log2p)->p()+1)*(basis::tri(x.log2p)->p()+1)*fabs(vslp));
        FLT meshc = gbl->adis/(fabs(vslpu) +100.0*FLT_EPSILON);
#endif
        res(2,i) = -sign*res(1,i)*vslpu*meshc;
        
        /* CONTINUITY FLUX */
        FLT mdot = sign*rhou*unorm_rel*jcb;
        res(3,i) = mdot;
        /* U-MOMENTUM FLUX */
        res(4,i) = mdot*uu +pu*norm(0);
        /* V-MOMENTUM FLUX */
        res(5,i) = mdot*vu +pu*norm(1);
        /* ENERGY FLUX */
        res(6,i) = mdot*Eu +pu*sign*unorm*jcb;
        
//        *x.gbl->log << indx << ' ' << i << ' ' << res(0,i) << ' ' << res(1,i) << ' ' << res(2,i) << ' ' << res(3,i) << ' ' << res(4,i) << ' ' << res(5,i) << ' ' << res(6,i) << std::endl;
        
	}
    
	lf = 0.0;
	/* INTEGRATE & STORE shock SOURCE TERM */
	basis::tri(x.log2p)->intgrt1d(&lf(0)(0),&res(3,0));
	basis::tri(x.log2p)->intgrt1d(&lf(1)(0),&res(4,0));
	basis::tri(x.log2p)->intgrt1d(&lf(2)(0),&res(5,0));
	basis::tri(x.log2p)->intgrt1d(&lf(3)(0),&res(6,0));
    
	/* INTEGRATE & STORE MESH MOVEMENT RESIDUALS */
    basis::tri(x.log2p)->intgrtx1d(&lf(x.NV)(0),&res(0,0));
	basis::tri(x.log2p)->intgrt1d(&lf(x.NV+1)(0),&res(1,0));
	basis::tri(x.log2p)->intgrtx1d(&lf(x.NV+1)(0),&res(2,0));
    
	return;
}

int shock::shock_mach(FLT &M, FLT cu, FLT cd, FLT vdiff) {
    const FLT gam = x.gbl->gamma;
    const FLT RHS = (gam+1)/(2*cu)*(2*cd/(gam-1)+vdiff);
    const int maxit = 100;
    const FLT tol = 1.0e-13;

    int it = 0;
    FLT dM = 1.0;
    while(fabs(dM)>tol && it<maxit) {
        FLT M2 = pow(M,2.0);
        FLT f = 1/M*(1/(gam -1)*sqrt((2*gam*M2 -(gam -1))*((gam -1)*M2 +2)) +M2 -1) -RHS;

        FLT fp = ((2*(gam -1)*(2*gam*M2 -gam +1))/M + (4*gam*((gam -1)*M2 +2))/M - (2*((gam -1)*M2 +2)*(2*gam*M2 -gam +1))/(M2*M))/(2*(gam -1)*sqrt(((M2*(gam -1) +2)*(2*M2*gam -gam +1))/M2)) - (M2 -1)/M2 +2;
        
        dM = -f/fp;
        M = M +dM;
        it = it+1;

    }
    if (it >= maxit) {
        return(1);
    }
    return(0);
}

void shock::element_jacobian_opp(int indx, Array<FLT,2>& K) {
	
	const int sm = basis::tri(x.log2p)->sm();
#ifdef WAY1
    const int nrows = x.NV+NV;
#else
    const int nrows = x.NV;
#endif

	Array<TinyVector<FLT,MXTM>,1> Rbar(x.NV+tri_mesh::ND),lf(x.NV+tri_mesh::ND);
	
	/* Calculate and store initial residual */
	int sind = base.seg(indx);
	x.ugtouht1d(sind);
	element_rsdl(indx,Rbar);
	
	Array<FLT,1> dw(x.NV);
	dw = 0.0;
    for(int i=0;i<2;++i){
        for(int n=0;n<x.NV;++n){
            dw(n) = dw(n) + fabs(u_opp_v(i+indx,n));
        }
    }
	
	dw *= eps_r;
	dw += eps_a;
    
	/* Numerically create Jacobian */
	int kcol = 0;
	for(int mode = 0; mode < 2; ++mode) {
		for(int var = 0; var < x.NV; ++var) {
			u_opp_v(mode+indx,var) += dw(var);

			element_rsdl(indx,lf);
			
			int krow = 0;
			for(int k=0;k<sm+2;++k)
                for(int n=0;n<nrows;++n)
					K(krow++,kcol) = (lf(n)(k)-Rbar(n)(k))/dw(var);
			++kcol;
			u_opp_v(mode+indx,var) -= dw(var);
		}
	}
    
    dw = 0.0;
    for(int i=0;i<sm;i++){
        for(int n=0;n<x.NV;n++){
            dw(n) = dw(n)+fabs(u_opp_e(indx,i,n));
        }
    }

    dw *= eps_r;
    dw += eps_a;

	for(int mode = 0; mode < sm; ++mode) {
		for(int var = 0; var < x.NV; ++var) {
			u_opp_e(indx,mode,var) += dw(var);
			
			element_rsdl(indx,lf);
			
			int krow = 0;
			for(int k=0;k<sm+2;++k)
                for(int n=0;n<nrows;++n)
                    K(krow++,kcol) = (lf(n)(k)-Rbar(n)(k))/dw(var);
			++kcol;
			u_opp_e(indx,mode,var) -= dw(var);
		}
	}
    
	return;
}

void shock::element_jacobian(int indx, Array<FLT,2>& K) {
    hp_coupled_bdry::element_jacobian(indx, K);
    
    const int sm = basis::tri(x.log2p)->sm();

    /* Multiply flow Jacobian terms for shock velocity equation by 2 because they are missing the remote parts */
    /* Don't do mesh jacobian terms because these do not have cross terms */
    int kcol = 0;
    for(int mode = 0; mode < 2; ++mode) {
        for(int var = 0; var < x.NV; ++var) {
            int krow = x.NV+NV-1;
            for(int k=0;k<sm+2;++k) {
                K(krow,kcol) *= 2.0;
                krow += x.NV+NV;
            }
            ++kcol;
        }
        kcol += tri_mesh::ND;
    }
    
    
    for(int mode = 2; mode < sm+2; ++mode) {
        for(int var = 0; var < x.NV; ++var) {
            int krow = x.NV+NV-1;;
            for(int k=0;k<sm+2;++k) {
                K(krow,kcol) *= 2.0;
                krow += x.NV+NV;
            }
            ++kcol;
        }
        kcol += tri_mesh::ND;
    }
}


#ifdef petsc
void shock::petsc_jacobian() {
    
    hp_coupled_bdry::petsc_jacobian();
    
    const int sm = basis::tri(x.log2p)->sm();
    const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
#ifdef WAY1
    const int nrows = vdofs;
#else
    const int nrows = x.NV;
#endif
    Array<FLT,2> K(nrows*(sm+2),x.NV*(sm+2));
    Array<int,1> loc_to_glo_row(nrows*(sm+2));
    Array<int,1> loc_to_glo_col(x.NV*(sm+2));


    /* Send flow indices to opposite boundary */
    int cnt = 0;
    int sind;
    for (int j=base.nseg-1;j>=0;--j) {
        sind = base.seg(j);
        base.isndbuf(cnt++) = x.seg(sind).pnt(1)*vdofs+x.jacobian_start;
        base.isndbuf(cnt++) = x.npnt*vdofs+x.NV*sm*sind+x.jacobian_start;
    }
    base.isndbuf(cnt++) = x.seg(sind).pnt(0)*vdofs+x.jacobian_start;
    base.sndsize() = cnt;
    base.sndtype() = boundary::int_msg;

    base.comm_prepare(boundary::all,0,boundary::symmetric);
    base.comm_exchange(boundary::all,0,boundary::symmetric);
    base.comm_wait(boundary::all,0,boundary::symmetric);
    base.comm_finish(boundary::all,0,boundary::symmetric,boundary::replace);

    if (gbl->symmetric || is_master) {
        /* This is effect of variables u,v,p,x,y on */
        /* source terms added to flow residuals */
        /* and x,y mesh movement equations */
        for (int j=0;j<base.nseg;++j) {
            int sind = base.seg(j);

            /* CREATE COLUMN NUMBERING LIST */
            int ind = 0;
            int gindx = base.isndbuf(2*j);
            for(int var = 0; var < x.NV; ++var)
                loc_to_glo_col(ind++) = gindx++;

            gindx = base.isndbuf(2*(j+1));
            for(int var = 0; var < x.NV; ++var)
                loc_to_glo_col(ind++) = gindx++;

            int count = 0;
            for(int mode = 0; mode < sm; ++mode) {
                for(int var = 0; var < x.NV; ++var)
                    loc_to_glo_col(ind++) = base.isndbuf(2*j+1)+count++;
            }

	    /* CREATE ROW NUMBERING LIST */
            ind = 0;
            for(int mode = 0; mode < 2; ++mode) {
                int gindx = vdofs*x.seg(sind).pnt(mode);
                for(int var = 0; var < nrows; ++var)
                    loc_to_glo_row(ind++) = gindx++;
            }

            int gindxNV = x.npnt*vdofs +x.NV*sind*sm;
#ifdef WAY1
            int gindxND = jacobian_start +j*tri_mesh::ND*sm;
#endif
            for(int mode = 0; mode < sm; ++mode) {
                for(int var = 0; var < x.NV; ++var)
                    loc_to_glo_row(ind++) = gindxNV++;
#ifdef WAY1
                for(int var = 0; var < tri_mesh::ND; ++var)
                    loc_to_glo_row(ind++) = gindxND++;
#endif
            }
            element_jacobian_opp(j,K);

#ifdef MY_SPARSE
            x.J_mpi.add_values(nrows*(sm+2),loc_to_glo_row,x.NV*(sm+2),loc_to_glo_col,K);
#else
            MatSetValuesLocal(x.petsc_J,nrows*(sm+2),loc_to_glo_row.data(),x.NV*(sm+2),loc_to_glo_col.data(),K.data(),ADD_VALUES);
#endif
        }
    }

}

#ifdef WAY1
void shock::non_sparse_snd(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi) {
    base.sndsize() = 0;
    base.sndtype() = boundary::int_msg;
    return;
}
#endif


int shock::non_sparse_rcv(int phase, Array<int,1> &nnzero, Array<int,1> &nnzero_mpi) {
    
    if (!base.is_comm() || base.matchphase(boundary::all_phased,0) != phase) return(0);
    
    int count = hp_coupled_bdry::non_sparse_rcv(phase, nnzero, nnzero_mpi);
    
    const int sm=basis::tri(x.log2p)->sm();
    const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
    const int begin_seg = x.npnt*vdofs;

    /* if local match then entries will go in to local jacobian */
    Array<int,1> target;
    if (base.is_local(0))  // only 1 matching boundary
        target.reference(nnzero);
    else
        target.reference(nnzero_mpi);

    int p0,p1;
    for (int i=0;i<base.nseg;++i) {
        int sind = base.seg(i);
        p0 = x.seg(sind).pnt(0)*vdofs;
        p1 = x.seg(sind).pnt(1)*vdofs;
#ifdef WAY1
        for(int n=0;n<vdofs;++n) {
#else
        for(int n=0;n<x.NV;++n) {
#endif
            target(p0 +n) += (2+sm)*x.NV;
            target(p1 +n) += (1+sm)*x.NV;
        }
    }
#ifdef WAY1
    for(int n=0;n<vdofs;++n) {
#else
    for(int n=0;n<x.NV;++n) {
#endif
        target(p1 +n) += x.NV;
    }

    // Now add to side degrees of freedom
    if (sm) {
        for (int i=0;i<base.nseg;++i) {
            int sind = base.seg(i);
            for (int mode=0;mode<sm;++mode) {
                for(int n=0;n<x.NV;++n) {
                    target(begin_seg +sind*sm*x.NV +mode*x.NV +n) += (2+sm)*x.NV;
                }
            }
        }
    }

    // Last thing to receive is nnzero for edge equations
#ifdef WAY1
    if (sm) {
        for (int i=0;i<base.nseg;++i) {
            for (int m=0;m<sm;++m) {
                for(int n=0;n<NV;++n) {
                    target(jacobian_start+i*sm*NV +m*NV +n) = (2+sm)*x.NV;
//                    target(jacobian_start+i*sm*NV +m*NV +n) = (2+sm)*vdofs;
                }
            }
        }
    }
#endif
    
    return(0);
}
#endif
