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
//#include "../../tri_mesh/tri_boundary.h"

#include <iostream>

//#define MPDEBUG

//#define DEBUG

//#define BODYFORCE


using namespace bdry_cns;

// extern FLT body[ND];


void shock::init(input_map& inmap,void* gin) {
	std::string keyword,matching_block,side_id,master_block,master_id;
	std::istringstream data;
	std::string filename;

	const int sm = basis::tri(x.log2p)->sm();
    gbl = static_cast<global *>(gin);

    inmap[base.idprefix+"_symmetric"] = "1";
    
	hp_coupled_bdry::init(inmap,gin);
    c0_indices_xy.clear();
    
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
        int sign = 1;
        for(k=0;k<ebp1;++k) {
            for(int n = 0; n<x.NV; n++){
                base.fsndbuf(count++) = sign*sdata[offset +n];
            }
            offset+=x.NV;
            sign*=-1;
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
	TinyMatrix<FLT,tri_mesh::ND,MXGP> crd, dcrd, mvel_u, mvel, mvel_d;
	TinyMatrix<FLT,7,MXGP> res;
	Array<TinyVector<FLT,MXTM>,1> uht_opp(x.NV);
    Array<FLT,1> mvel_u_new(MXGP), shock_vel(MXGP),  mvel_u_norm(MXGP), mvel_norm(MXGP);

	for (int n=0; n<x.NV; n++){
		uht_opp(n)(0) = u_opp_v(indx,n);
		uht_opp(n)(1) = u_opp_v(indx+1,n);
		for (int j=0; j<sm; j++){
			uht_opp(n)(j+2) = u_opp_e(indx,j,n);
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
    
	for(i=0;i<basis::tri(x.log2p)->gpx();++i) {
		norm(0) =  dcrd(1,i);
		norm(1) = -dcrd(0,i);
        
		jcb = sqrt(norm(0)*norm(0) +norm(1)*norm(1));
        
        FLT pu, uu, vu, RTu, rhou, Eu, cu, pd, ud, vd, RTd, rhod, cd, shock_sign;
		if (u(0)(i)<u_opp(0)(i)){
			pu = u(0)(i);
			uu = u(1)(i);
			vu = u(2)(i);
			RTu = u(3)(i);

			pd = u_opp(0)(i);
			ud = u_opp(1)(i);
			vd = u_opp(2)(i);
			RTd = u_opp(3)(i);
            
            shock_sign = 1.0;
            
		}
		else{
			pd = u(0)(i);
			ud = u(1)(i);
			vd = u(2)(i);
			RTd = u(3)(i);

			pu = u_opp(0)(i);
			uu = u_opp(1)(i);
			vu = u_opp(2)(i);
			RTu = u_opp(3)(i);
            
            shock_sign = -1.0;
            
		}
        
        
        rhou = pu/RTu;
		rhod = pd/RTd;
		cu = sqrt(x.gbl->gamma*RTu);
		cd = sqrt(x.gbl->gamma*RTd);
        Eu = RTu/(x.gbl->gamma-1.0)+0.5*(uu*uu+vu*vu);
        
        //Is flow in + or - direction?
        int flag = 0;
//        if(uu > 0){
//            flag = 0;
//        }
//        else{
//            flag = 1;
//        }
        
		/* Call NewtZou to find Mu, use in shock velocity evaluation */
        FLT normal_uu = (uu*norm(0)+vu*norm(1))/jcb;
        FLT normal_ud = (ud*norm(0)+vd*norm(1))/jcb;
        
        /* RELATIVE VELOCITY STORED IN MVEL(N)*/
        if (flag == 0){
            mvel_u(0,i) = uu -(x.gbl->bd(0)*(crd(0,i) -dxdt(x.log2p,indx)(0,i)));
#ifdef MESH_REF_VEL
            mvel_u(0,i) -= x.gbl->mesh_ref_vel(0);
#endif
            mvel_u(1,i) = vu -(x.gbl->bd(0)*(crd(1,i) -dxdt(x.log2p,indx)(1,i)));
#ifdef MESH_REF_VEL
            mvel_u(1,i) -= x.gbl->mesh_ref_vel(1);
#endif
            
            
            
            mvel_d(0,i) = ud -(x.gbl->bd(0)*(crd(0,i) -dxdt(x.log2p,indx)(0,i)));
#ifdef MESH_REF_VEL
            mvel_d(0,i) -= x.gbl->mesh_ref_vel(0);
#endif
            mvel_d(1,i) = vd -(x.gbl->bd(0)*(crd(1,i) -dxdt(x.log2p,indx)(1,i)));
#ifdef MESH_REF_VEL
            mvel_d(1,i) -= x.gbl->mesh_ref_vel(1);
#endif
        }
        else{
            /* RELATIVE VELOCITY STORED IN MVEL(N)*/
            mvel_u(0,i) = (x.gbl->bd(0)*(crd(0,i) -dxdt(x.log2p,indx)(0,i)))-uu;
#ifdef MESH_REF_VEL
            mvel_u(0,i) -= x.gbl->mesh_ref_vel(0);
#endif
            mvel_u(1,i) = (x.gbl->bd(0)*(crd(1,i) -dxdt(x.log2p,indx)(1,i)))-vu;
#ifdef MESH_REF_VEL
            mvel_u(1,i) -= x.gbl->mesh_ref_vel(1);
#endif
            
            
            mvel_d(0,i) = (x.gbl->bd(0)*(crd(0,i) -dxdt(x.log2p,indx)(0,i)))-ud;
#ifdef MESH_REF_VEL
            mvel_d(0,i) -= x.gbl->mesh_ref_vel(0);
#endif
            mvel_d(1,i) = (x.gbl->bd(0)*(crd(1,i) -dxdt(x.log2p,indx)(1,i)))-vd;
#ifdef MESH_REF_VEL
            mvel_d(1,i) -= x.gbl->mesh_ref_vel(1);
#endif
        }
        
        norm(0) = shock_sign*norm(0);
        norm(1) = shock_sign*norm(1);
        
        normal_uu = (uu*norm(0)+vu*norm(1))/jcb;
        normal_ud = (ud*norm(0)+vd*norm(1))/jcb;

        mvel_norm(i) = ((x.gbl->bd(0)*(crd(0,i) -dxdt(x.log2p,indx)(0,i)))*norm(0)+(x.gbl->bd(0)*(crd(1,i) -dxdt(x.log2p,indx)(1,i)))*norm(1))/jcb;
        mvel_u_norm(i) = (mvel_u(0,i)*norm(0)+mvel_u(1,i)*norm(1))/jcb;
        
        
        FLT Mu = mvel_u_norm(i)/cu;
        NewtZou(Mu, normal_ud, cd, normal_uu, cu, flag, crd(0,i), crd(1,i));


        if (flag == 0){
            shock_vel(i) = normal_uu-cu*Mu;
        }
        else{
            shock_vel(i) = normal_uu+cu*Mu;
        }
        

        /* TANGENTIAL SPACING */
        res(0,i) = -shock_sign*ksprg(indx)*jcb;
        /* NORMAL FLUX */
        res(1,i) = (mvel_norm(i)-shock_vel(i))*jcb;
        /* UPWINDING BASED ON TANGENTIAL VELOCITY */
        norm(0) = shock_sign*norm(0);
        norm(1) = shock_sign*norm(1);
        
        /* Calculate stabilization constant based on analysis of linear elements and constant tau */
        FLT dMda = -(x.gbl->gamma + 1.0)/(2.0*cu*(((2*(x.gbl->gamma - 1.0)*(2.0*x.gbl->gamma*Mu*Mu - x.gbl->gamma + 1.0))/Mu + (4.0*x.gbl->gamma*((x.gbl->gamma - 1.0)*Mu*Mu + 2.0))/Mu - (2.0*((x.gbl->gamma - 1.0)*Mu*Mu + 2.0)*(2.0*x.gbl->gamma*Mu*Mu - x.gbl->gamma + 1.0))/(Mu*Mu*Mu))/(2.0*(x.gbl->gamma - 1.0)*sqrt((((Mu*Mu*(x.gbl->gamma - 1.0) + 2.0)*(2.0*Mu*Mu*x.gbl->gamma - x.gbl->gamma + 1.0))/(Mu*Mu)))) - (Mu*Mu - 1.0)/(Mu*Mu) + 2.0));
        FLT dMdb = -dMda;
        FLT tan_rel = -mvel_u(0,i)*norm(1) +mvel_u(1,i)*norm(0);
        FLT tan_u = -norm(1)*uu +norm(0)*vu;
        FLT tan_d = -norm(1)*ud +norm(0)*vd;
        FLT vslp = (tan_rel-cu*tan_u*dMda+cu*tan_d*dMdb)/jcb;
        gbl->meshc(indx) = gbl->adis/((basis::tri(x.log2p)->p()+1)*(basis::tri(x.log2p)->p()+1)*fabs(vslp));
        
        res(2,i) = -res(1,i)*(tan_rel-cu*tan_u*dMda+cu*tan_d*dMdb)/jcb*gbl->meshc(indx);


        normal_uu = (uu*norm(0)+vu*norm(1))/jcb;
        normal_ud = (ud*norm(0)+vd*norm(1))/jcb;

        mvel_norm(i) = ((x.gbl->bd(0)*(crd(0,i) -dxdt(x.log2p,indx)(0,i)))*norm(0)+(x.gbl->bd(0)*(crd(1,i) -dxdt(x.log2p,indx)(1,i)))*norm(1))/jcb;
        if (flag == 0){
            mvel_u_new(i) = (normal_uu-mvel_norm(i))*jcb;
        }
        else{
            mvel_u_new(i) = (mvel_norm(i)-normal_uu)*jcb;
        }
        
        /* CONTINUITY FLUX */
        res(3,i) = rhou*mvel_u_new(i);
        /* U-MOMENTUM FLUX */
        res(4,i) = rhou*uu*mvel_u_new(i)+pu*norm(0);
        /* V-MOMENTUM FLUX */
        res(5,i) = rhou*vu*mvel_u_new(i)+pu*norm(1);
        /* ENERGY FLUX */
        res(6,i) = rhou*Eu*mvel_u_new(i)+pu*normal_uu*jcb;
        
        
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


void shock::NewtZou(FLT &Mu, FLT vd, FLT cd, FLT vu, FLT cu, int flag, FLT crdx, FLT crdy){
    
    FLT Jr, f, fp;
    if(flag == 0){
        Jr = 2.0*cd/(x.gbl->gamma-1.0)-vd;
    f=(1.0/(x.gbl->gamma-1.0))*sqrt(((2.0*x.gbl->gamma*pow(Mu,2.0))-(x.gbl->gamma-1.0))*(((x.gbl->gamma-1.0)*pow(Mu,2.0))+2.0)/pow(Mu,2.0))+(((pow(Mu,2.0))-1.0)/Mu)-((x.gbl->gamma+1.0)/(2.0*cu))*(Jr+vu);
    }
    else{
        Jr = 2.0*cd/(x.gbl->gamma-1.0)+vd;
    f=(1.0/(x.gbl->gamma-1.0))*sqrt(((2.0*x.gbl->gamma*pow(Mu,2.0))-(x.gbl->gamma-1.0))*(((x.gbl->gamma-1.0)*pow(Mu,2.0))+2.0)/pow(Mu,2.0))+(((pow(Mu,2.0))-1.0)/Mu)-((x.gbl->gamma+1.0)/(2.0*cu))*(Jr-vu);
    }
    fp = ((2.0*(x.gbl->gamma - 1.0)*(2.0*x.gbl->gamma*pow(Mu,2.0) - x.gbl->gamma + 1.0))/Mu + (4.0*x.gbl->gamma*((x.gbl->gamma - 1.0)*pow(Mu,2.0) + 2.0))/Mu - (2.0*((x.gbl->gamma - 1.0)*pow(Mu,2.0) + 2.0)*(2.0*x.gbl->gamma*pow(Mu,2.0) - x.gbl->gamma + 1.0))/pow(Mu,3.0))/(2.0*(x.gbl->gamma - 1.0)*sqrt(((pow(Mu,2.0)*(x.gbl->gamma - 1.0) + 2.0)*(2.0*pow(Mu,2.0)*x.gbl->gamma - x.gbl->gamma + 1.0))/pow(Mu,2.0))) - (pow(Mu,2.0) - 1.0)/pow(Mu,2.0) + 2.0;
    
    int it = 0;
    int maxit = 100;
    FLT tol = 1.0e-13;
    while(fabs(f)>tol && it<maxit){
        Mu = Mu-(f/fp);
        
        if(flag == 0){
        f=(1.0/(x.gbl->gamma-1.0))*sqrt(((2.0*x.gbl->gamma*pow(Mu,2.0))-(x.gbl->gamma-1.0))*(((x.gbl->gamma-1.0)*pow(Mu,2.0))+2.0)/pow(Mu,2.0))+(((pow(Mu,2.0))-1.0)/Mu)-((x.gbl->gamma+1.0)/(2.0*cu))*(Jr+vu);
        }
        else{
        f=(1.0/(x.gbl->gamma-1.0))*sqrt(((2.0*x.gbl->gamma*pow(Mu,2.0))-(x.gbl->gamma-1.0))*(((x.gbl->gamma-1.0)*pow(Mu,2.0))+2.0)/pow(Mu,2.0))+(((pow(Mu,2.0))-1.0)/Mu)-((x.gbl->gamma+1.0)/(2.0*cu))*(Jr-vu);
        }
    
        fp = ((2.0*(x.gbl->gamma - 1.0)*(2.0*x.gbl->gamma*pow(Mu,2.0) - x.gbl->gamma + 1.0))/Mu + (4.0*x.gbl->gamma*((x.gbl->gamma - 1.0)*pow(Mu,2.0) + 2.0))/Mu - (2.0*((x.gbl->gamma - 1.0)*pow(Mu,2.0) + 2.0)*(2.0*x.gbl->gamma*pow(Mu,2.0) - x.gbl->gamma + 1.0))/pow(Mu,3.0))/(2.0*(x.gbl->gamma - 1.0)*sqrt(((pow(Mu,2.0)*(x.gbl->gamma - 1.0) + 2.0)*(2.0*pow(Mu,2.0)*x.gbl->gamma - x.gbl->gamma + 1.0))/pow(Mu,2.0))) - (pow(Mu,2.0) - 1.0)/pow(Mu,2.0) + 2.0;
        it = it+1;
    }
    if (it == maxit){
        std::cout << "Newton Solver in Zou exceeded iteration limit at (" << crdx << "," << crdy << ")" << std::endl;
    }
}



void shock::element_jacobian_opp(int indx, Array<FLT,2>& K) {
	
	int sm = basis::tri(x.log2p)->sm();

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
                for(int n=0;n<x.NV+NV;++n)
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
                for(int n=0;n<x.NV+NV;++n)
                    K(krow++,kcol) = (lf(n)(k)-Rbar(n)(k))/dw(var);
			++kcol;
			u_opp_e(indx,mode,var) -= dw(var);
		}
	}
    
	return;
}


#ifdef petsc
void shock::petsc_jacobian() {
    
    hp_coupled_bdry::petsc_jacobian();
    
    const int sm = basis::tri(x.log2p)->sm();
    const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;

    Array<FLT,2> K(vdofs*(sm+2),x.NV*(sm+2));
    Array<int,1> loc_to_glo_row(vdofs*(sm+2));
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
                for(int var = 0; var < vdofs; ++var)
                    loc_to_glo_row(ind++) = gindx++;
            }

            int gindxNV = x.npnt*vdofs +x.NV*sind*sm;
            int gindxND = jacobian_start +j*tri_mesh::ND*sm;
            for(int mode = 0; mode < sm; ++mode) {
                for(int var = 0; var < x.NV; ++var)
                    loc_to_glo_row(ind++) = gindxNV++;

                for(int var = 0; var < tri_mesh::ND; ++var)
                    loc_to_glo_row(ind++) = gindxND++;
            }
            element_jacobian_opp(j,K);

#ifdef MY_SPARSE
            x.J_mpi.add_values(vdofs*(sm+2),loc_to_glo_row,x.NV*(sm+2),loc_to_glo_col,K);
#else
            MatSetValuesLocal(x.petsc_J,vdofs*(sm+2),loc_to_glo_row.data(),vdofs*(sm+2),loc_to_glo_col.data(),K.data(),ADD_VALUES);
#endif
        }
    }
}

void shock::non_sparse_snd(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi) {
    base.sndsize() = 0;
    base.sndtype() = boundary::int_msg;
    return;
}


int shock::non_sparse_rcv(int phase, Array<int,1> &nnzero, Array<int,1> &nnzero_mpi) {
    
    if (!base.is_comm() || base.matchphase(boundary::all_phased,0) != phase) return(0);
    
//    int count = hp_coupled_bdry::non_sparse_rcv(phase, nnzero, nnzero_mpi);
    
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
        for(int n=0;n<vdofs;++n) {
//        for(int n=0;n<x.NV;++n) {
            target(p0 +n) += (2+sm)*x.NV;
            target(p1 +n) += (1+sm)*x.NV;
        }
    }
    for(int n=0;n<vdofs;++n) {
//    for(int n=0;n<x.NV;++n) {
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
    if (sm) {
        for (int i=0;i<base.nseg;++i) {
            for (int m=0;m<sm;++m) {
                for(int n=0;n<NV;++n) {
                    //target(jacobian_start+i*sm*NV +m*NV +n) = (2+sm)*x.NV;
                    target(jacobian_start+i*sm*NV +m*NV +n) = (2+sm)*vdofs;
                }
            }
        }
    }
    
    return(0);
}
#endif
