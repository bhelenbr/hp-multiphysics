//#include <utilities.h>
//#include <myblas.h>

#include "tri_hp_cns.h"
#include "../hp_boundary.h"

void tri_hp_cns::setup_preconditioner() {
	/* SET-UP PRECONDITIONER */
	int tind,i,j,side,v0;
	FLT jcb,h,hmax,q,qmax,lam1,gam,pmax,rtmax;
	TinyVector<int,3> v;
	Array<double,1> umax(NV);
	
	/***************************************/
	/** DETERMINE FLOW PSEUDO-TIME STEP ****/
	/***************************************/
	gbl->vprcn_ut(Range(0,npnt-1),Range::all(),Range::all()) = 0.0;
	if (basis::tri(log2p)->sm() > 0) {
		gbl->sprcn_ut(Range(0,nseg-1),Range::all(),Range::all()) = 0.0;
	}
	gbl->tprcn_ut(Range(0,ntri-1),Range::all(),Range::all()) = 0.0;

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
			hmax = 0.0;
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
			hmax = 0.0;
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

		FLT tstep;
		Array<double,2> tprcn(NV,NV),tau(NV,NV);		
		
		pennsylvania_peanut_butter(umax,hmax,tprcn,tau,tstep);
		
		gbl->tprcn_ut(tind,Range::all(),Range::all()) = jcb*tprcn/tstep;

		gbl->tau(tind,Range::all(),Range::all()) = adis*tau/jcb;
		
		for(i=0;i<3;++i) {
			gbl->vprcn_ut(v(i),Range::all(),Range::all())  += basis::tri(log2p)->vdiag()*gbl->tprcn_ut(tind,Range::all(),Range::all());
			if (basis::tri(log2p)->sm() > 0) {
				side = tri(tind).seg(i);
				gbl->sprcn_ut(side,Range::all(),Range::all()) += gbl->tprcn_ut(tind,Range::all(),Range::all());
			}
		}
	}
	
	tri_hp::setup_preconditioner();
}


void tri_hp_cns::pennsylvania_peanut_butter(Array<double,1> pvu, FLT hmax, Array<FLT,2> &Pinv, Array<FLT,2> &Tau, FLT &timestep) {
	
	Array<FLT,2> P(NV,NV), A(NV,NV), V(NV,NV), VINV(NV,NV), B(NV,NV), S(NV,NV), Tinv(NV,NV), temp(NV,NV);
	Array<FLT,1> Aeigs(NV),Beigs(NV);
	
	FLT gam = gbl->gamma;
	FLT gm1 = gam-1.0;
	FLT gogm1 = gam/gm1;
	FLT pr = pvu(0);
	FLT u = pvu(1);
	FLT v = pvu(2);
	FLT rt = pvu(NV-1);
	FLT rho = pr/rt;
	FLT ke = 0.5*(u*u+v*v);
	FLT E = rt/gm1+ke;
	FLT c2 = gam*rt;
	FLT c = sqrt(c2);
	
//	/* Preconditioner */
//	P = ke*gm1,          -u*gm1,     -v*gm1,      gm1,
//	    -u/rho,          1.0/rho,    0.0,         0.0,
//	    -v/rho,          0.0,        1.0/rho,     0.0,
//		(gm1*ke-rt)/rho, -u*gm1/rho, -v*gm1/rho, gm1/rho;	
//
//	/* Inverse of Preconditioner */
//	Pinv = 1.0/rt,               0.0,   0.0,   -rho/rt,
//		   u/rt,                 rho,   0.0,   -rho*u/rt,
//		   v/rt,                 0.0,   rho,   -rho*v/rt,
//		   (rt+gm1*ke)/(gm1*rt), rho*u, rho*v, -rho*ke/rt;		
	
	
	/* Preconditioner */
	P = 1,0,0,0,
		0,1,0,0,
		0,0,1,0,
		0,0,0,1;
	
	/* Inverse of Preconditioner */
	Pinv = 1,0,0,0,
		   0,1,0,0,
		   0,0,1,0,
		   0,0,0,1;
	
//	/* eigenvectors of P*df/dw */
//	V = 0.0, 0.0, 1.0, 1.0,
//	    0.0, 0.0, c/(gam*pr), -c/(gam*pr),
//	    1.0, 0.0, 0.0, 0.0,
//	    0.0, 1.0, gm1/(gam*rho), gm1/(gam*rho);
//	
//	/* take absolute value, u and c are already positive */
//	if (u > c) {
//		Aeigs = u,u,u+c,u-c;
//	}
//	else {
//		Aeigs = u,u,u+c,c-u;	
//	}
//
//	/* inverse of eigenvectors (P*df/dw)^-1 */		
//	VINV = 0.0,            0.0,           1.0, 0.0,
//		   -gm1/(gam*rho), 0.0,           0.0, 1.0,
//	       0.5,            0.5*gam*pr/c,  0.0, 0.0,
//		   0.5,            -0.5*gam*pr/c, 0.0, 0.0;
//	
//	for(int i=0; i < NV; ++i)
//		for(int j=0; j < NV; ++j)
//			VINV(i,j) = Aeigs(i)*VINV(i,j);
//	
//	A = 0.0;
//	for(int i=0; i<NV; ++i)
//		for(int j=0; j<NV; ++j)
//			for(int k=0; k<NV; ++k)
//				A(i,j)+=V(i,k)*VINV(k,j);
	
	
	
	
	/* df/dw derivative of fluxes wrt primitive variables */
	A = u/rt,               rho,                       0.0,     -rho*u/rt,
		u*u/rt+1.0,         2.0*rho*u,                 0.0,     -rho*u*u/rt,
		u*v/rt,             rho*v,                     rho*u,   -rho*u*v/rt,
		u*(gogm1*rt+ke)/rt, rho*(gogm1*rt+ke)+rho*u*u, rho*u*v, -rho*u*(gogm1*rt+ke)/rt+rho*u*gogm1;
	
	temp = 0.0;
	for(int i=0; i<NV; ++i)
		for(int j=0; j<NV; ++j)
			for(int k=0; k<NV; ++k)
				temp(i,j)+=P(i,k)*A(k,j);
	A = temp;
	matrix_absolute_value(A);
	
	
	

//	/* df/dw derivative of fluxes wrt conservative variables*/
//	A =  0.0, 1.0, 0.0, 0.0,
//		 -u*u+gm1*ke, (2.0-gm1)*u, -v*gm1, gm1,
//	     -u*v, v, u, 0.0,
//	     u*(-gam*E+2.0*gm1*ke), gam*E-gm1*u*u-gm1*ke, -u*v*gm1, u*gam;
//
//	temp = 0.0;
//	for(int i=0; i<NV; ++i)
//		for(int j=0; j<NV; ++j)
//			for(int k=0; k<NV; ++k)
//				temp(i,j)+=P(i,k)*B(k,j);
//	
//	matrix_absolute_value(A);
	
	
	
	
	
//	/* eigenvectors of P*dg/dw */
//	V = 0.0, 0.0, 1.0,           1.0,
//		1.0, 0.0, 0.0,           0.0,
//		0.0, 0.0, c/(gam*pr),    -c/(gam*pr),
//		0.0, 1.0, gm1/(gam*rho), gm1/(gam*rho);
//	
//	if (v > c) {
//		Beigs = v,v,v+c,v-c;
//	}
//	else {
//		Beigs = v,v,v+c,c-v;	
//	}
//	
//	VINV = 0.0,            1.0, 0.0,           0.0,
//	       -gm1/(gam*rho), 0.0, 0.0,           1.0,
//		   0.5,            0.0, 0.5*gam*pr/c,  0.0,
//		   0.5,            0.0, -0.5*gam*pr/c, 0.0;	
//	
//	for(int i=0; i < NV; ++i)
//		for(int j=0; j < NV; ++j)
//			VINV(i,j) = Beigs(i)*VINV(i,j);
//	
//	B = 0.0;
//	for(int i=0; i<NV; ++i)
//		for(int j=0; j<NV; ++j)
//			for(int k=0; k<NV; ++k)
//				B(i,j)+=V(i,k)*VINV(k,j);	
	
	
	
	/* dg/dw derivative of fluxes wrt primitive variables*/
	B = v/rt,               0.0,       rho,                         -rho*v/rt,
		u*v/rt,            rho*v,    rho*u,                      -rho*u*v/rt,
		v*v/rt+1.0,        0.0,       2.0*rho*v,                  -rho*v*v/rt,
		v*(gogm1*rt+ke)/rt, rho*u*v, rho*(gogm1*rt+ke)+rho*v*v, -rho*v*(gogm1*rt+ke)/rt+rho*v*gogm1;
	
	temp = 0.0;
	for(int i=0; i<NV; ++i)
		for(int j=0; j<NV; ++j)
			for(int k=0; k<NV; ++k)
				temp(i,j)+=P(i,k)*B(k,j);
	B = temp;
	matrix_absolute_value(B);
	
	
	
//	/* dg/dw derivative of fluxes wrt conservative variables*/
//	B =   0.0, 0.0, 1.0, 0.0,
//	      -u*v, v, u, 0.0,
//		  -v*v+gm1*ke, -u*gm1, (2.0-gm1)*v, gm1,
//		  v*(-gam*E+2.0*gm1*ke),  -u*v*gm1, gam*E-gm1*v*v-gm1*ke, v*gam;
//	
//	temp = 0.0;
//	for(int i=0; i<NV; ++i)
//		for(int j=0; j<NV; ++j)
//			for(int k=0; k<NV; ++k)
//				temp(i,j)+=P(i,k)*B(k,j);
//	matrix_absolute_value(B);
	
	
	
	FLT nu = gbl->mu/rho;
	
	FLT cp = gogm1*gbl->R;
	FLT alpha = gbl->kcond/(rho*cp);
	
	S = 0.0, 0.0, 0.0, 0.0,
		0.0, nu,  0.0, 0.0,
		0.0, 0.0, nu,  0.0,
		0.0, 0.0, 0.0, alpha;
	
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
	
	if (info != 0) {
		*gbl->log << "DGETRF FAILED FOR CNS EXPLICIT TSTEP" << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
	
	for (int i = 0; i < NV; ++i)
		for (int j = 0; j < NV; ++j)
			temp(i,j)=P(j,i);
	
	/* Solve transposed system temp' = inv(Tinv')*temp' */
	char trans[] = "T";
	GETRS(trans,NV,NV,Tinv.data(),NV,ipiv,temp.data(),NV,info);
	
	if (info != 0) {
		*gbl->log << "DGETRS FAILED FOR CNS EXPLICIT TSTEP" << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
	
	for (int i = 0; i < NV; ++i)
		for (int j = 0; j < NV; ++j)
			Tau(i,j)=temp(j,i);
	
	return;
}

