#include <math.h>

#include "tri_hp_cns.h"
#include "../hp_boundary.h"

// #define TIMEACCURATE
#define REFINED_WAY

void tri_hp_cns::setup_preconditioner() {
	/* SET-UP PRECONDITIONER */
	int tind,i,j,side,v0;
	FLT jcb,h,hmax,q,qmax,lam1,gam,pmax,rtmax;
	TinyVector<int,3> v;
	TinyMatrix<FLT,2> P(4,4),Pinv(4,4),A(4,4),B(4,4),S(4,4),Tinv(4,4),temp(4,4);
	int info,ipiv[NV];
	char trans[] = "T";
	
	FLT nu = gbl->mu/gbl->rho;

	/***************************************/
	/** DETERMINE FLOW PSEUDO-TIME STEP ****/
	/***************************************/
	gbl->vprcn_ut(Range(0,npnt-1),Range::all(),Range::all()) = 0.0;
	if (basis::tri(log2p)->sm() > 0) {
		gbl->sprcn_ut(Range(0,nseg-1),Range::all(),Range::all()) = 0.0;
	}
	gbl->tprcn_ut(Range(0,ntri-1),Range::all(),Range::all()) = 0.0;

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
			int lgpx = basis::tri(log2p)->gpx(), lgpn = basis::tri(log2p)->gpn();
			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpn;++j) {
					
					mvel(0) = gbl->bd(0)*(crd(0)(i,j) -dxdt(log2p,tind,0)(i,j));
					mvel(1) = gbl->bd(0)*(crd(1)(i,j) -dxdt(log2p,tind,1)(i,j));                       
					q = pow(u(1)(i,j)-0.5*mvel(0),2.0)  +pow(u(2)(i,j)-0.5*mvel(1),2.0);
					qmax = MAX(qmax,q);
					pmax = MAX(pmax,fabs(u(0)(i,j)));
					rtmax = MAX(rtmax,fabs(u(3)(i,j)));
				
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
			exit(1);
		}

		if  (std::isnan(qmax)) { 
			*gbl->log << "flow solution has nan's" << std::endl;
			output("nan",tecplot);
			exit(1);
		}
		q = sqrt(qmax);
		FLT op = 1.0/pmax;
		FLT rt = rtmax;
		FLT gm1 = gam-1;
		
		P =       qmax*gm1,        -q*gm1,       -q*gm1,     gm1,
			      -q*op*rt,         op*rt,        0.0,       0.0,
			      -q*op*rt           0.0,        op*rt,      0.0,
			op*rt*(gm1*qmax-rt),-op*rt*q*gm1,-op*rt*q*gm1,op*rt*gm1;

		
		/* Calculate and store tprcn (4x4) matrix */	
		
		FLT ort = 1.0/rtmax;
		
		Pinv =         ort,        0.0,          0.0       -pmax*ort*ort,
				      ort*q,     pmax*ort,       0.0,     -pmax*ort*ort*q,
				      ort*q,       0.0,       pmax*ort,   -pmax*ort*ort*q,
				1.0/gm1+rt*qmax, pmax*ort*q, pmax*ort*q,-pmax*ort*ort*qmax;


		gbl->tprcn_ut(tind,Range::all(),Range::all())=Pinv;//temp replace Pinv with this

		A =         ort*q,                 pmax*ort,                    0.0,                        -pmax*ort*ort*q,
				ort*qmax+1.0,          2.0*pmax*ort*q,                  0.0,                     -pmax*ort*ort*qmax,
			     ort*qmax,               pmax*ort*q,                pmax*ort*q,                    -pmax*ort*ort*qmax,
		   q*(gam/gm1+qmax*ort),pmax*((gam/gm1+qmax*ort)+ort*qmax),pmax*ort*qmax,-pmax*ort*ort*q*(gam/gm1*rtmax+qmax)+pmax*ort*q*gam/gm1;
		
		temp = 0.0;
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				for (int k = 0; k < n; ++k)
					temp(i,j)+=P(i,k)*A(k,j);	
		A=temp;
		matrix_absolute_value(A, NV);
	

		B =         ort*q,            0.0,						pmax*ort,                -pmax*ort*ort*q,
				  ort*qmax,         pmax*ort*q,                pmax*ort*q,                 -pmax*ort*ort*qmax,
				ort*qmax+1.0,         0.0,                     2.0*pmax*ort*q,               -pmax*ort*ort*qmax,
			q*(gam/gm1+qmax*ort),pmax*ort*qmax,pmax*((gam/gm1+qmax*ort)+ort*qmax),-pmax*ort*ort*q*(gam/gm1*rtmax+qmax)+pmax*ort*q*gam/gm1;
		
		temp = 0.0;
		for (int i = 0; i < NV; ++i)
			for (int j = 0; j < NV; ++j)
				for (int k = 0; k < NV; ++k)
					temp(i,j)+=P(i,k)*B(k,j);	
		B=temp;
		matrix_absolute_value(B, NV);
		
		S=0,0, 0, 0,
		  0,nu,0, 0,
		  0,0, nu,0,
		  0,0, 0, nu;
			
		S=pmax/hmax*S;//temp fix me
		
		Tinv=2.0/hmax*(A+B+hmax*S);// temp fix me, cancel out hmax?
		
		/* Write a routine given tprcn, a, b, s, returns dt & tau */
		dt=spectral_radius(Tinv,NV);
		
		/*  LU factorization  */
		GETRF(NV, NV, Tinv.data(), NV, ipiv, info);
		
		for (int i = 0; i < NV; ++i)
			for (int j = 0; j < NV; ++j)
				temp(i,j)=P(j,i);
		
		/* Solve transposed system temp = inv(Tinv)*temp */
		GETRS(trans,NV,NV,Tinv.data(),NV,ipiv,temp.data(),n,info);
		
		for (int i = 0; i < NV; ++i)
			for (int j = 0; j < NV; ++j)
				gbl->tau(tind,i,j)=temp(j,i);
					  
		/* FROM HERE BELOW WILL CHANGE */
		lam1 = q + sqrt(qmax +gam);

		/* SET UP DISSIPATIVE COEFFICIENTS */
		gbl->tau(tind,0) = adis*h/(jcb*sqrt(gam));
		gbl->tau(tind,NV-1) = qmax*gbl->tau(tind,0);

		/* SET UP DIAGONAL PRECONDITIONER */
		// jcb *= 8.*nu*(1./(hmax*hmax) +1./(h*h)) +2*lam1/h +2*sqrt(gam)/hmax +gbl->bd(0);
		jcb *= 2.*nu*(1./(hmax*hmax) +1./(h*h)) +3*lam1/h;  // heuristically tuned

#ifdef TIMEACCURATE
		dtstari = MAX((nu/(h*h) +lam1/h +gbl->bd(0)),dtstari);

	}
	printf("#iterative to physical time step ratio: %f\n",gbl->bd(0)/dtstari);

	for(tind=0;tind<ntri;++tind) {
		v = tri(tind).pnt;
		jcb = 0.25*area(tind)*dtstari;
#endif

		jcb *= RAD((pnts(v(0))(0) +pnts(v(1))(0) +pnts(v(2))(0))/3.);

		gbl->tprcn_ut(tind,0,0) = gbl->rho*jcb;    
		gbl->tprcn_ut(tind,1,1) = gbl->rho*jcb;      
		gbl->tprcn_ut(tind,2,2) = jcb/gam;
		gbl->tprcn_ut(tind,0,2) = jcb*ubar/gam;
		gbl->tprcn_ut(tind,1,2) = jcb*vbar/gam;
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


