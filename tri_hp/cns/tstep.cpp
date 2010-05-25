//#include <utilities.h>
//#include <myblas.h>

#include "tri_hp_cns.h"
#include "../hp_boundary.h"

// #define TIMEACCURATE
//#define REFINED_WAY

void tri_hp_cns::setup_preconditioner() {
	/* SET-UP PRECONDITIONER */
	int tind,i,j,side,v0;
	FLT jcb,h,hmax,q,qmax,lam1,gam,pmax,rtmax;
	TinyVector<int,3> v;
	Array<double,1> umax(NV);
	
	//FLT nu = gbl->mu/gbl->rho;

	/***************************************/
	/** DETERMINE FLOW PSEUDO-TIME STEP ****/
	/***************************************/
	gbl->vprcn_ut(Range(0,npnt-1),Range::all(),Range::all()) = 0.0;
	if (basis::tri(log2p)->sm() > 0) {
		gbl->sprcn_ut(Range(0,nseg-1),Range::all(),Range::all()) = 0.0;
	}
	gbl->tprcn_ut(Range(0,ntri-1),Range::all(),Range::all()) = 0.0;

//#ifdef TIMEACCURATE
//	FLT dtstari = 0.0;
//#endif

	for(tind = 0; tind < ntri; ++tind) {
		jcb = 0.25*area(tind);  // area is 2 x triangle area
		v = tri(tind).pnt;
		
		/* LOAD SOLUTION COEFFICIENTS FOR THIS ELEMENT */
		/* PROJECT SOLUTION TO GAUSS POINTS WITH DERIVATIVES IF NEEDED FOR VISCOUS TERMS */
		ugtouht(tind);
		for(int n=0;n<ND;++n)
			basis::tri(log2p)->proj(&uht(n)(0),&u(n)(0,0),MXGP);
			
		/* IF TINFO > -1 IT IS CURVED ELEMENT */
		if (tri(tind).info > -1) {
			/* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
			crdtocht(tind);

			/* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
			for(int n=0;n<ND;++n)
				basis::tri(log2p)->proj_bdry(&cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), &dcrd(n,1)(0,0),MXGP);
		
			TinyVector<FLT,ND> mvel;
			qmax = 0.0;
			hmax = 0.0;
			pmax = 0.0;
			rtmax = 0.0;
			umax = 0.0;
			FLT jcbmin = jcb;
			int lgpx = basis::tri(log2p)->gpx(), lgpn = basis::tri(log2p)->gpn();
			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpn;++j) {
					
					mvel(0) = gbl->bd(0)*(crd(0)(i,j) -dxdt(log2p,tind,0)(i,j));
					mvel(1) = gbl->bd(0)*(crd(1)(i,j) -dxdt(log2p,tind,1)(i,j));                  
					jcbmin = MIN(jcbmin,dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j));
					
					/* CALCULATE CURVED SIDE LENGTHS */
					h = 0.0;
					for (int n=0;n<ND;++n)
						h += dcrd(n,0)(i,j)*dcrd(n,0)(i,j);
					hmax = MAX(h,hmax);
						
					h = 0.0;
					for (int n=0;n<ND;++n)
						h += dcrd(n,1)(i,j)*dcrd(n,1)(i,j);
					hmax = MAX(h,hmax);

					h = 0.0;
					for (int n=0;n<ND;++n)
						h += (dcrd(n,1)(i,j) -dcrd(n,0)(i,j))*(dcrd(n,1)(i,j) -dcrd(n,0)(i,j));
					hmax = MAX(h,hmax);
				
					q = pow(u(1)(i,j)-0.5*mvel(0),2.0)  +pow(u(2)(i,j)-0.5*mvel(1),2.0);
					qmax = MAX(qmax,q);
					
					pmax = MAX(pmax,fabs(u(0)(i,j)));
					rtmax = MAX(rtmax,fabs(u(3)(i,j)));
					
					umax(0) = MAX(umax(0),fabs(u(0)(i,j)));
					umax(1) = MAX(umax(1),fabs(u(1)(i,j)-0.5*mvel(0)));
					umax(2) = MAX(umax(2),fabs(u(2)(i,j)-0.5*mvel(1)));
					umax(3) = MAX(umax(3),fabs(u(3)(i,j)));
						
					
				}
			}	
			hmax = 2.*sqrt(hmax);
			h = 4.*jcbmin/(0.25*(basis::tri(log2p)->p() +1)*(basis::tri(log2p)->p()+1)*hmax);
			hmax = hmax/(0.25*(basis::tri(log2p)->p() +1)*(basis::tri(log2p)->p()+1));
		}
		else {
			/* PROJECT VERTEX COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
			for(int n=0;n<ND;++n)
				basis::tri(log2p)->proj(pnts(v(0))(n),pnts(v(1))(n),pnts(v(2))(n),&crd(n)(0,0),MXGP);
		
			TinyVector<FLT,ND> mvel;
			qmax = 0.0;
			hmax = 0.0;
			pmax = 0.0;
			rtmax = 0.0;
			umax = 0.0;
			int lgpx = basis::tri(log2p)->gpx(), lgpn = basis::tri(log2p)->gpn();
			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpn;++j) {
					
					mvel(0) = gbl->bd(0)*(crd(0)(i,j) -dxdt(log2p,tind,0)(i,j));
					mvel(1) = gbl->bd(0)*(crd(1)(i,j) -dxdt(log2p,tind,1)(i,j));                       
					q = pow(u(1)(i,j)-0.5*mvel(0),2.0)  +pow(u(2)(i,j)-0.5*mvel(1),2.0);
					qmax = MAX(qmax,q);
					pmax = MAX(pmax,fabs(u(0)(i,j)));
					rtmax = MAX(rtmax,fabs(u(3)(i,j)));
					
					umax(0) = MAX(umax(0),fabs(u(0)(i,j)));
					umax(1) = MAX(umax(1),fabs(u(1)(i,j)-0.5*mvel(0)));
					umax(2) = MAX(umax(2),fabs(u(2)(i,j)-0.5*mvel(1)));
					umax(3) = MAX(umax(3),fabs(u(3)(i,j)));
				
				}
			}
			for(j=0;j<3;++j) {
				h = pow(pnts(v(j))(0) -pnts(v((j+1)%3))(0),2.0) +pow(pnts(v(j))(1) -pnts(v((j+1)%3))(1),2.0);
				hmax = (h > hmax ? h : hmax);
			}
			hmax = sqrt(hmax);
			h = 4.*jcb/(0.25*(basis::tri(log2p)->p() +1)*(basis::tri(log2p)->p()+1)*hmax);
			hmax = hmax/(0.25*(basis::tri(log2p)->p() +1)*(basis::tri(log2p)->p()+1));
		}			

		if (!(h > 0.0)) { 
			*gbl->log << "negative triangle area caught in tstep. Problem triangle is : " << tind << std::endl;
			*gbl->log << "approximate location: " << pnts(v(0))(0) << ' ' << pnts(v(0))(1) << std::endl;
			tri_mesh::output("negative",grid);
			sim::abort(__LINE__,__FILE__,gbl->log);
		}

		if  (std::isnan(qmax)) { 
			*gbl->log << "flow solution has nan's" << std::endl;
			output("nan",tecplot);
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
		q = sqrt(qmax);

		FLT tstep;
		Array<double,2> tprcn(NV,NV),tau(NV,NV);		
		
		pennsylvania_peanut_butter(umax,hmax,tprcn,tau,tstep);
		
		gbl->tprcn_ut(tind,Range::all(),Range::all())=tprcn;
		gbl->tau(tind,Range::all(),Range::all())=adis*tau/jcb;
		
	
		for(i=0;i<3;++i) {
			gbl->vprcn_ut(v(i),Range::all(),Range::all())  += gbl->tprcn_ut(tind,Range::all(),Range::all());
			if (basis::tri(log2p)->sm() > 0) {
				side = tri(tind).seg(i);
				gbl->sprcn_ut(side,Range::all(),Range::all()) += gbl->tprcn_ut(tind,Range::all(),Range::all());
			}
		}
	}
	
	tri_hp::setup_preconditioner();
}


void tri_hp_cns::pennsylvania_peanut_butter(Array<double,1> u, FLT hmax, Array<FLT,2> &Pinv, Array<FLT,2> &Tau, FLT &timestep) {
	
	Array<FLT,2> P(NV,NV),A(NV,NV),B(NV,NV),S(NV,NV),Tinv(NV,NV),temp(NV,NV);
	
	FLT gm1 = gbl->gamma-1;
	FLT gogm1 = gbl->gamma/gm1;
	FLT pr = u(0);
	FLT uv = u(1);
	FLT vv = u(2);
	FLT rt = u(3);	
	FLT ke = 0.5*(uv*uv+vv*vv);
	
	/* Preconditioner */
	P = ke*gm1,            -uv*gm1,       -vv*gm1,       gm1,
		-uv/pr*rt,         rt/pr,         0.0,           0.0,
		-vv*rt/pr,         0.0,           1.0/pr*rt,     0.0,
		rt*(gm1*ke-rt)/pr, -rt*uv*gm1/pr, -rt*vv*gm1/pr, rt*gm1/pr;
	
	/* Inverse of Preconditioner */
	Pinv = 1.0/rt,        0.0,      0.0,      -pr/rt/rt,
		uv/rt,         pr/rt,    0.0,      -uv*pr/rt/rt,
		vv/rt,         0.0,      pr/rt,    -vv*pr/rt/rt,
		1.0/gm1+ke/rt, uv*pr/rt, vv*pr/rt, -pr/rt/rt*ke;	
	
	/* df/dw */
	A = uv/rt,               pr/rt,                           0.0,         -pr/rt/rt*uv,
		uv*uv/rt+1.0,        2.0*pr/rt*uv,                    0.0,         -pr/rt/rt*uv*uv,
		uv*vv/rt,            pr/rt*vv,                        pr/rt*uv,    -pr/rt/rt*uv*vv,
		uv*(gogm1*rt+ke)/rt, pr/rt*(gogm1*rt+ke)+pr/rt*uv*uv, pr/rt*uv*vv, -pr/rt/rt*uv*(gogm1*rt+ke)+pr/rt*uv*gogm1;
	
	temp = 0.0;
	for(int i=0; i<NV; ++i)
		for(int j=0; j<NV; ++j)
			for(int k=0; k<NV; ++k)
				temp(i,j)+=P(i,k)*A(k,j);
	A = temp;
	
	matrix_absolute_value(A);
	
	/* dg/dw */
	B = vv/rt,               0.0,         pr/rt,                           -pr/rt/rt*vv,
		uv*vv/rt,            pr/rt*vv,    pr/rt*uv,                        -pr/rt/rt*uv*vv,
		vv*vv/rt+1.0,        0.0,         2.0*pr/rt*vv,                    -pr/rt/rt*vv*vv,
		vv*(gogm1*rt+ke)/rt, pr/rt*uv*vv, pr/rt*(gogm1*rt+ke)+pr/rt*vv*vv, -pr/rt/rt*vv*(gogm1*rt+ke)+pr/rt*vv*gogm1;
	
	temp = 0.0;
	for(int i=0; i<NV; ++i)
		for(int j=0; j<NV; ++j)
			for(int k=0; k<NV; ++k)
				temp(i,j)+=P(i,k)*B(k,j);
	
	B = temp;
	
	matrix_absolute_value(B);
	
	S = 0.0, 0.0,     0.0,     0.0,
		0.0, gbl->mu, 0.0,     0.0,
		0.0, 0.0,     gbl->mu, 0.0,
		0.0, 0.0,     0.0,     gbl->kcond;
	
	S = Pinv*gbl->dti+1.0/hmax/hmax*S;
	
	temp = 0.0;
	for(int i=0; i<NV; ++i)
		for(int j=0; j<NV; ++j)
			for(int k=0; k<NV; ++k)
				temp(i,j)+=P(i,k)*S(k,j);
	
	S = temp;
	
	Tinv = 2.0/hmax*(A+B+hmax*S);
	
	/* smallest eigenvalue of Tau tilde */
	timestep = 1.0/spectral_radius(Tinv);
	
	/*  LU factorization  */
	int info,ipiv[NV];
	GETRF(NV, NV, Tinv.data(), NV, ipiv, info);
	
	for (int i = 0; i < NV; ++i)
		for (int j = 0; j < NV; ++j)
			temp(i,j)=P(j,i);
	
	/* Solve transposed system temp' = inv(Tinv')*temp' */
	char trans[] = "T";
	GETRS(trans,NV,NV,Tinv.data(),NV,ipiv,temp.data(),NV,info);
	
	for (int i = 0; i < NV; ++i)
		for (int j = 0; j < NV; ++j)
			Tau(i,j)=temp(j,i);
	
	
	return;
}
