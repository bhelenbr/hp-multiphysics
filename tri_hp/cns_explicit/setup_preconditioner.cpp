//#include <utilities.h>
//#include <myblas.h>

#include "tri_hp_cns_explicit.h"
#include "../hp_boundary.h"

//#define TIMEACCURATE

void tri_hp_cns_explicit::setup_preconditioner() {
	/* SET-UP PRECONDITIONER */
	int tind,i,j,side;
	FLT jcb,h,hmax;
	TinyVector<int,3> v;
	Array<double,1> umax(NV),cvu(NV);
    int err = 0;
	
	/***************************************/
	/** DETERMINE FLOW PSEUDO-TIME STEP ****/
	/***************************************/
	gbl->vprcn(Range(0,npnt-1),Range::all()) = 0.0;
	if (basis::tri(log2p)->sm() > 0) {
		gbl->sprcn(Range(0,nseg-1),Range::all()) = 0.0;
	}
	gbl->tprcn(Range(0,ntri-1),Range::all()) = 0.0;

#ifdef TIMEACCURATE
	FLT dtstari = 0.0;
#endif

	for(tind = 0; tind < ntri; ++tind) {
		jcb = 0.25*area(tind);  // area is 2 x triangle area
		v = tri(tind).pnt;
	
		/* LOAD SOLUTION COEFFICIENTS FOR THIS ELEMENT */
		/* PROJECT SOLUTION TO GAUSS POINTS WITH DERIVATIVES IF NEEDED FOR VISCOUS TERMS */
		ugtouht(tind);
		for(int n=0;n<ND;++n)
			basis::tri(log2p)->proj(&uht(n)(0),&u(n)(0,0),MXGP);
	
		int lgpx = basis::tri(log2p)->gpx(), lgpn = basis::tri(log2p)->gpn();

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
			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpn;++j) {
					
					mvel(0) = gbl->bd(0)*(crd(0)(i,j) -dxdt(log2p)(tind,0,i,j));
					mvel(1) = gbl->bd(0)*(crd(1)(i,j) -dxdt(log2p)(tind,1,i,j)); 
#ifdef MESH_REF_VEL
					mvel(0) += gbl->mesh_ref_vel(0);
					mvel(1) += gbl->mesh_ref_vel(1);
#endif
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
				
					umax(0) = MAX(umax(0),fabs(u(0)(i,j)));
					umax(1) = MAX(umax(1),fabs(u(1)(i,j)-0.5*mvel(0)));// temp fix because u(1) is rhou not u
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
					
					mvel(0) = gbl->bd(0)*(crd(0)(i,j) -dxdt(log2p)(tind,0,i,j));
					mvel(1) = gbl->bd(0)*(crd(1)(i,j) -dxdt(log2p)(tind,1,i,j));
#ifdef MESH_REF_VEL
					mvel(0) += gbl->mesh_ref_vel(0);
					mvel(1) += gbl->mesh_ref_vel(1);
#endif
				
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
            err = 1;
            break;
		}

		if  (std::isnan(umax(0)) && std::isnan(umax(1)) && std::isnan(umax(2)) && std::isnan(umax(3))) { 
			*gbl->log << "flow solution has nan's" << std::endl;
            err = 1;
            break;
		}
		
		FLT tstep;
		Array<double,2> tprcn(NV,NV),tau(NV,NV);		
		
		calculate_tau(umax,hmax,tprcn,tau,tstep);

		/* SET UP DISSIPATIVE COEFFICIENTS */
		gbl->tau(tind,Range::all(),Range::all())=adis*tau/jcb;
		
		/* SET UP DIAGONAL PRECONDITIONER */
		jcb /= tstep;  // temp fix me

#ifdef TIMEACCURATE
		dtstari = MAX(1./tstep,dtstari);
		
	}
	
	/* find max dtstari for all blocks and use on every block  */
	FLT dtstari_recv;
	sim::blks.allreduce(&dtstari,&dtstari_recv,1,blocks::flt_msg,blocks::max);
	dtstari = dtstari_recv;
	
	*gbl->log << "#iterative to physical time step ratio: " << gbl->bd(0)/dtstari << ' ' << gbl->bd(0) << ' ' << dtstari << '\n';

	
	for(tind=0;tind<ntri;++tind) {
		v = tri(tind).pnt;
		jcb = 0.25*area(tind)*dtstari; 
#endif
		
		jcb *= RAD((pnts(v(0))(0) +pnts(v(1))(0) +pnts(v(2))(0))/3.);

		gbl->tprcn(tind,0) = jcb;
		gbl->tprcn(tind,1) = jcb;
		gbl->tprcn(tind,2) = jcb;
		gbl->tprcn(tind,3) = jcb;

		for(i=0;i<3;++i) {
			gbl->vprcn(v(i),Range::all())  += gbl->tprcn(tind,Range::all());
			if (basis::tri(log2p)->sm() > 0) {
				side = tri(tind).seg(i);
				gbl->sprcn(side,Range::all()) += gbl->tprcn(tind,Range::all());
			}
		}
	}
	
	return(tri_hp::setup_preconditioner()+err);
}

void tri_hp_cns_explicit::calculate_tau(Array<double,1> cvu, FLT h, Array<FLT,2> &Pinv, Array<FLT,2> &Tau, FLT &timestep) {
	
	Array<FLT,2> P(NV,NV), A(NV,NV), V(NV,NV), VINV(NV,NV), B(NV,NV), S(NV,NV), Tinv(NV,NV), temp(NV,NV);
	Array<FLT,1> Aeigs(NV),Beigs(NV);
	
	FLT gam = gbl->gamma;
	FLT gm1 = gam-1.0;
	FLT gogm1 = gam/gm1;
	FLT rho = cvu(0);
	FLT rhou = cvu(1);
	FLT u = rhou/rho;
	FLT rhov = cvu(2);
	FLT v = rhov/rho;
	FLT rhoE = cvu(3);
	FLT E = rhoE/rho;
	FLT ke = 0.5*(u*u+v*v);
	FLT c2 = gam*gm1*(E-ke);
	FLT c = sqrt(c2);
	
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
	
	/* eigenvectors of df/dw */
	V = 0.0, 1.0,           1.0,              1.0,
		0.0, u,             u+c,              u-c,
		1.0, 0.0,           v,                v,
		v,   0.5*(u*u-v*v), gam*E-gm1*ke+u*c, gam*E-gm1*ke-u*c;
	
	/* take absolute value, u and c are already positive */
	if (u > c) {
		Aeigs = u,u,u+c,u-c;
	}
	else {
		Aeigs = u,u,u+c,c-u;	
	}
		
	/* inverse of eigenvectors of df/dw */
	VINV= -v*ke*gm1,         v*u*gm1,         c2+gm1*v*v, -v*gm1,
	      -gm1*ke+c2,        u*gm1,           v*gm1,      -gm1,
	      0.5*(-u*c+gm1*ke), -0.5*(-c+u*gm1), -0.5*v*gm1, 0.5*gm1,
	      0.5*(u*c+gm1*ke),  -0.5*(c+u*gm1),  -0.5*v*gm1, 0.5*gm1;
			
	VINV /= c2;
	
	for(int i=0; i < NV; ++i)
		for(int j=0; j < NV; ++j)
			VINV(i,j) = Aeigs(i)*VINV(i,j);

	A = 0.0;
	for(int i=0; i<NV; ++i)
		for(int j=0; j<NV; ++j)
			for(int k=0; k<NV; ++k)
				A(i,j)+=V(i,k)*VINV(k,j);

	V = 0.0, 1.0, 1.0, 1.0,
		1.0, 0.0, u, u,
		0.0, v, v+c, v-c,
		u,-0.5*(u*u-v*v), gam*E-gm1*ke+v*c, gam*E-gm1*ke-v*c;
	
	if (v > c) {
		Beigs = v,v,v+c,v-c;
	}
	else {
		Beigs = v,v,v+c,c-v;	
	}
	
	VINV = -u*ke*gm1,c2+gm1*u*u, v*u*gm1, -u*gm1,
	       -gm1*ke+c2, u*gm1, v*gm1, -gm1,
		   0.5*(-v*c+gm1*ke),-0.5*u*gm1, -0.5*(-c+v*gm1), 0.5*gm1,
	       0.5*(v*c+gm1*ke), -0.5*u*gm1, -0.5*(c+v*gm1), 0.5*gm1;
		
	VINV /= c2;

	for(int i=0; i < NV; ++i)
		for(int j=0; j < NV; ++j)
			VINV(i,j) = Beigs(i)*VINV(i,j);
	
	B = 0.0;
	for(int i=0; i<NV; ++i)
		for(int j=0; j<NV; ++j)
			for(int k=0; k<NV; ++k)
				B(i,j)+=V(i,k)*VINV(k,j);	
	
	FLT nu = gbl->mu/rho;
	
	FLT cp = gogm1*gbl->R;
	FLT alpha = gbl->kcond/(rho*cp);
	
	S = 0.0, 0.0,      0.0,      0.0,
		0.0, nu/(h*h), 0.0,      0.0,
		0.0, 0.0,      nu/(h*h), 0.0,
		0.0, 0.0,      0.0,      alpha/(h*h);

	for(int i=0; i<NV; ++i)
		S(i,i) += gbl->bd(0);
	
	Tinv = 2.0/h*(A+B+h*S);

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
