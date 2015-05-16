//#include <utilities.h>
//#include <myblas.h>

#include "tri_hp_cns.h"
#include "../hp_boundary.h"

//#define TIMEACCURATE

void tri_hp_cns::setup_preconditioner() {
	/* SET-UP PRECONDITIONER */
	int tind,i,j,side;
	FLT jcb,h,hmax,q,qmax,tstep,dtstari;
	TinyVector<int,3> v;
	Array<double,1> umax(NV),ubar(NV);
	Array<double,2> tprcn(NV,NV),tau(NV,NV),dpdc(NV,NV);

	/***************************************/
	/** DETERMINE FLOW PSEUDO-TIME STEP ****/
	/***************************************/
	if(gbl->diagonal_preconditioner){
		gbl->vprcn(Range(0,npnt-1),Range::all()) = 0.0;
		gbl->vpreconditioner(Range(0,npnt-1),Range::all(),Range::all()) = 0.0;
		if (basis::tri(log2p)->sm() > 0) {
			gbl->sprcn(Range(0,nseg-1),Range::all()) = 0.0;
			gbl->spreconditioner(Range(0,nseg-1),Range::all(),Range::all()) = 0.0;
		}
		gbl->tprcn(Range(0,ntri-1),Range::all()) = 0.0;
		gbl->tpreconditioner(Range(0,ntri-1),Range::all(),Range::all()) = 0.0;

		
	} else {

		gbl->vprcn_ut(Range(0,npnt-1),Range::all(),Range::all()) = 0.0;
		if (basis::tri(log2p)->sm() > 0) {
			gbl->sprcn_ut(Range(0,nseg-1),Range::all(),Range::all()) = 0.0;
		}
		gbl->tprcn_ut(Range(0,ntri-1),Range::all(),Range::all()) = 0.0;
	}

#ifdef TIMEACCURATE
	FLT dtstarimax = 0.0;
#endif
	
	for(tind = 0; tind < ntri; ++tind) {
		jcb = 0.25*area(tind);  // area is 2 x triangle area
		v = tri(tind).pnt;
		
		/* LOAD SOLUTION COEFFICIENTS FOR THIS ELEMENT */
		/* PROJECT SOLUTION TO GAUSS POINTS WITH DERIVATIVES IF NEEDED FOR VISCOUS TERMS */
		ugtouht(tind);
		for(int n=0;n<NV;++n)
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
			ubar = 0.0;
			FLT jcbmin = jcb;
			int lgpx = basis::tri(log2p)->gpx(), lgpn = basis::tri(log2p)->gpn();
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
				
					q = pow(u(1)(i,j)-0.5*mvel(0),2.0)  +pow(u(2)(i,j)-0.5*mvel(1),2.0);
					qmax = MAX(qmax,q);
					
					umax(0) = MAX(umax(0),fabs(u(0)(i,j)));
					umax(1) = MAX(umax(1),fabs(u(1)(i,j)-0.5*mvel(0)));
					umax(2) = MAX(umax(2),fabs(u(2)(i,j)-0.5*mvel(1)));
					umax(3) = MAX(umax(3),fabs(u(3)(i,j)));
					
					for(int n=0;n<NV;++n)
						ubar(n) += u(n)(i,j);
						
					
				}
			}	
			hmax = 2.*sqrt(hmax);
			h = 4.*jcbmin/(0.25*(basis::tri(log2p)->p() +1)*(basis::tri(log2p)->p()+1)*hmax);
			hmax = hmax/(0.25*(basis::tri(log2p)->p() +1)*(basis::tri(log2p)->p()+1));
			ubar /= lgpx*lgpn;
		}
		else {
			/* PROJECT VERTEX COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
			for(int n=0;n<ND;++n)
				basis::tri(log2p)->proj(pnts(v(0))(n),pnts(v(1))(n),pnts(v(2))(n),&crd(n)(0,0),MXGP);
		
			TinyVector<FLT,ND> mvel;
			hmax = 0.0;
			umax = 0.0;
			ubar = 0.0;
			int lgpx = basis::tri(log2p)->gpx(), lgpn = basis::tri(log2p)->gpn();
			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpn;++j) {
					
					mvel(0) = gbl->bd(0)*(crd(0)(i,j) -dxdt(log2p)(tind,0,i,j));
					mvel(1) = gbl->bd(0)*(crd(1)(i,j) -dxdt(log2p)(tind,1,i,j));
#ifdef MESH_REF_VEL
					mvel(0) += gbl->mesh_ref_vel(0);
					mvel(1) += gbl->mesh_ref_vel(1);
#endif
					q = pow(u(1)(i,j)-0.5*mvel(0),2.0)  +pow(u(2)(i,j)-0.5*mvel(1),2.0);
					qmax = MAX(qmax,q);
					
					umax(0) = MAX(umax(0),fabs(u(0)(i,j)));
					umax(1) = MAX(umax(1),fabs(u(1)(i,j)-0.5*mvel(0)));
					umax(2) = MAX(umax(2),fabs(u(2)(i,j)-0.5*mvel(1)));
					umax(3) = MAX(umax(3),fabs(u(3)(i,j)));
					
					for(int n=0;n<NV;++n)
						ubar(n) += u(n)(i,j);
				
				}
			}
			for(j=0;j<3;++j) {
				h = pow(pnts(v(j))(0) -pnts(v((j+1)%3))(0),2.0) +pow(pnts(v(j))(1) -pnts(v((j+1)%3))(1),2.0);
				hmax = (h > hmax ? h : hmax);
			}
			hmax = sqrt(hmax);
			h = 4.*jcb/(0.25*(basis::tri(log2p)->p() +1)*(basis::tri(log2p)->p()+1)*hmax);
			hmax = hmax/(0.25*(basis::tri(log2p)->p() +1)*(basis::tri(log2p)->p()+1));
			ubar /= lgpx*lgpn;

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

		
		pennsylvania_peanut_butter(umax,hmax,tprcn,tau,tstep);

		dtstari = 1.0/tstep;
		
		if(gbl->diagonal_preconditioner) {
			gbl->tpreconditioner(tind,Range::all(),Range::all()) = tprcn;
		}
		else {
			gbl->tprcn_ut(tind,Range::all(),Range::all()) = tprcn;
		}
		
		gbl->tau(tind,Range::all(),Range::all()) = adis*tau/jcb;

#ifdef TIMEACCURATE
		dtstarimax = MAX(dtstari,dtstarimax);
		
	}
	
	/* find max dtstari for all blocks and use on every block  */
	FLT dtstari_recv;
	sim::blks.allreduce(&dtstarimax,&dtstari_recv,1,blocks::flt_msg,blocks::max);
	dtstarimax = dtstari_recv;
	
	*gbl->log << "#iterative to physical time step ratio: " << gbl->bd(0)/dtstarimax << ' ' << gbl->bd(0) << ' ' << dtstarimax << '\n';
	
	
	for(tind=0;tind<ntri;++tind) {
		v = tri(tind).pnt;
		jcb = 0.25*area(tind); 
		dtstari = dtstarimax;
#endif
		
		if(gbl->diagonal_preconditioner){
			gbl->tprcn(tind,Range::all()) = jcb;
			gbl->tpreconditioner(tind,Range::all(),Range::all()) *= dtstari;
			
			for(i=0;i<3;++i) {
			
				gbl->vprcn(v(i),Range::all())  += gbl->tprcn(tind,Range::all());
				gbl->vpreconditioner(v(i),Range::all(),Range::all())  += gbl->tpreconditioner(tind,Range::all(),Range::all())*jcb;

				if (basis::tri(log2p)->sm() > 0) {
					side = tri(tind).seg(i);
					
					gbl->sprcn(side,Range::all()) += gbl->tprcn(tind,Range::all());
					gbl->spreconditioner(side,Range::all(),Range::all()) += gbl->tpreconditioner(tind,Range::all(),Range::all());
				}
			}
		}
		else {
			gbl->tprcn_ut(tind,Range::all(),Range::all()) *= jcb*dtstari;
			
			for(i=0;i<3;++i) {
				
				gbl->vprcn_ut(v(i),Range::all(),Range::all())  += gbl->tprcn_ut(tind,Range::all(),Range::all());
				
				if (basis::tri(log2p)->sm() > 0) {
					side = tri(tind).seg(i);
					gbl->sprcn_ut(side,Range::all(),Range::all()) += gbl->tprcn_ut(tind,Range::all(),Range::all());
				}
			}

		}
	}
	
	if(gbl->diagonal_preconditioner){
		for(i=0;i<npnt;++i) {
			gbl->vpreconditioner(i,Range::all(),Range::all()) /=  gbl->vprcn(i,0);
		}
	}
	
	// remember to do parallel communication
	
	tri_hp::setup_preconditioner();
	
	int last_phase,mp_phase;
	
	if(gbl->diagonal_preconditioner){
		for(int stage = 0; stage<NV; ++stage) {
			for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
				vc0load(mp_phase,gbl->vpreconditioner.data() +stage*NV,NV);
				pmsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
				last_phase = true;
				last_phase &= vc0wait_rcv(mp_phase,gbl->vpreconditioner.data()+stage*NV,NV);
			}
			if (log2p) {
				sc0load(gbl->spreconditioner.data()+stage*NV,0,0,NV);
				smsgpass(boundary::all,0,boundary::symmetric);
				sc0wait_rcv(gbl->spreconditioner.data()+stage*NV,0,0,NV);
			}
		}
	}
	
}

void tri_hp_cns::pennsylvania_peanut_butter(Array<double,1> pvu, FLT h, Array<FLT,2> &Pinv, Array<FLT,2> &Tau, FLT &timestep) {
	
	Array<double,2> P(NV,NV), V(NV,NV), VINV(NV,NV), dpdc(NV,NV), dcdp(NV,NV), A(NV,NV), B(NV,NV), S(NV,NV), Tinv(NV,NV), temp(NV,NV);
	Array<FLT,1> Aeigs(NV),Beigs(NV);
	
	FLT gam = gbl->gamma;
	FLT gm1 = gam-1.0;
	FLT gogm1 = gam/gm1;
	FLT pr = pvu(0);
	FLT u = pvu(1);
	FLT v = pvu(2);
	FLT rt = pvu(3);
	FLT rho = pr/rt;
	FLT ke = 0.5*(u*u+v*v);
	FLT c2 = gam*rt;
	FLT c = sqrt(c2);

	FLT nu = gbl->mu/rho;
	FLT cp = gogm1*gbl->R;
	FLT alpha = gbl->kcond/(rho*cp);
	
	/* need to tune better */
//	FLT hdt = 0.25*pow(h*gbl->bd(0),2.0);
//	FLT vel = 1.0*(u*u+v*v);
//	FLT nuh = 4.0*(pow(nu/h,2.0)+pow(alpha/h,2.0));
//
//	FLT umag = sqrt(hdt+vel+nuh); // with reynolds and prandtl dependence
//	
//	FLT M = MAX(1.0e-5,umag/c);
//	
//	FLT b2,alph;
//
//	if(M > .6) { // turn off preconditioner
//		b2 = 1.0;
//		alph = 0.0;
//	} else {
//		b2 = M*M/(1.0-M*M);
//		alph = 1.0+b2;
//	}
	
	
	

//	FLT hdt = .25*h*gbl->bd(0)/c;
//	FLT umag = sqrt(u*u+v*v);
//	FLT M = MAX(umag/c,1.0e-5);
//	FLT nuh = 4.0*(nu/h+alpha/h)/c;
//	
//	FLT beta,b2,alph;
//	
//	if(M > .6) { // turn off preconditioner
//		b2 = 1.0;
//		alph = 0.0;
//	} else {
//		beta = sqrt(M*M/(1.0-M*M))+hdt/(hdt+1.0)+nuh;
//		b2 = beta*beta;
//		alph = 1.0+b2;
//	}
	
	
	
	
	
	FLT hdt = 0.5*h*gbl->bd(0)/c;
	FLT umag = sqrt(u*u+v*v);
	FLT M = MAX(umag/c,1.0e-5);
	FLT nuh = 4.0*nu/(h*c);
	FLT alh = 2.0*alpha/(h*c);//maybe it should be smaller?
	
	FLT b2 = MIN(M*M/(1.0-M*M) + hdt*hdt + nuh*nuh + alh*alh,1.0);
	FLT alph = 1.0+b2;
	//cout  << b2 << ' ' <<  M*M << ' ' << M*M/(1.0-M*M) << ' ' << hdt*hdt << ' ' << nuh*nuh << ' ' << alh*alh << endl;

#ifdef petsc
	b2 = 1.0;
	alph = 0.0;
#endif

	alph = 0.0; // prevents wiggles when residual gets small, not sure why
	
	/* Preconditioner */
	P = b2,                   0.0, 0.0, 0.0,
	    -alph*u/(pr*gam),     1.0, 0.0, 0.0,
	    -alph*v/(pr*gam),     0.0, 1.0, 0.0,
	    (b2-1.0)/(gogm1*rho), 0.0, 0.0, 1.0;
	
	/* Inverse of Preconditioner */
	Pinv = 1.0/b2,					 0.0, 0.0, 0.0,
		   alph*u/(pr*gam*b2),		 1.0, 0.0, 0.0,
		   alph*v/(pr*gam*b2),		 0.0, 1.0, 0.0,
		   -(b2-1.0)/(gogm1*rho*b2), 0.0, 0.0, 1.0;

	/* jacobian of primitive wrt conservative */
	dpdc = ke*gm1,          -u*gm1,     -v*gm1,      gm1,
		    -u/rho,          1.0/rho,    0.0,        0.0,
		    -v/rho,          0.0,        1.0/rho,    0.0,
			(gm1*ke-rt)/rho, -u*gm1/rho, -v*gm1/rho, gm1/rho;	

	/* jacobian of primitive wrt conservative */
	dcdp = 1.0/rt,               0.0,   0.0,   -rho/rt,
		   u/rt,                 rho,   0.0,   -rho*u/rt,
	       v/rt,                 0.0,   rho,   -rho*v/rt,
	       (rt+gm1*ke)/(gm1*rt), rho*u, rho*v, -rho*ke/rt;	
	
	temp = 0.0;
	for(int i=0; i<NV; ++i)
		for(int j=0; j<NV; ++j)
			for(int k=0; k<NV; ++k)
				temp(i,j)+=P(i,k)*dpdc(k,j);
	P = temp;

//	/* df/dw derivative of fluxes wrt primitive variables */
//	A = u/rt,               rho,                       0.0,     -rho*u/rt,
//		u*u/rt+1.0,         2.0*rho*u,                 0.0,     -rho*u*u/rt,
//		u*v/rt,             rho*v,                     rho*u,   -rho*u*v/rt,
//		u*(gogm1*rt+ke)/rt, rho*(gogm1*rt+ke)+rho*u*u, rho*u*v, -rho*u*(gogm1*rt+ke)/rt+rho*u*gogm1;
//
//	temp = 0.0;
//	for(int i=0; i<NV; ++i)
//		for(int j=0; j<NV; ++j)
//			for(int k=0; k<NV; ++k)
//				temp(i,j)+=P(i,k)*A(k,j);
//	A = temp;	
//	matrix_absolute_value(A);
	
	FLT temp1 = sqrt(u*u*(1.0-2.0*b2+b2*b2)+4.0*b2*c2);
	
	V = 0.5*(u*(b2-1.0)+temp1)*rho, 0.5*(u*(b2-1.0)-temp1)*rho, 0.0, 0.0,
	    1.0,1.0,0.0,0.0,
	    0.0,0.0,1.0,0.0,
	    0.5*(u*(b2-1.0)*gm1+gm1*temp1)/gam, 0.5*(u*(b2-1.0)*gm1-gm1*temp1)/gam, 0.0, 1.0;
	
	Aeigs = 0.5*(u+u*b2+temp1), 0.5*(u+u*b2-temp1),u,u;

	for(int i=0; i<NV; ++i)
		Aeigs(i) = abs(Aeigs(i));

	VINV =  1.0/(temp1*rho), -0.5*(u*(b2-1.0)-temp1)/temp1,0.0,0.0,
		    -1.0/(temp1*rho), 0.5*(u*(b2-1.0)+temp1)/temp1,0.0,0.0,
			0.0, 0.0, 1.0, 0.0,
			-gm1/(gam*rho), 0.0, 0.0, 1.0;
	
	for(int i=0; i < NV; ++i)
		for(int j=0; j < NV; ++j)
			VINV(i,j) = Aeigs(i)*VINV(i,j);
	
	A = 0.0;
	for(int i=0; i<NV; ++i)
		for(int j=0; j<NV; ++j)
			for(int k=0; k<NV; ++k)
				A(i,j)+=V(i,k)*VINV(k,j);
	
//	/* dg/dw derivative of fluxes wrt primitive variables*/
//	B = v/rt,               0.0,     rho,                       -rho*v/rt,
//		u*v/rt,             rho*v,   rho*u,                     -rho*u*v/rt,
//		v*v/rt+1.0,         0.0,     2.0*rho*v,                 -rho*v*v/rt,
//		v*(gogm1*rt+ke)/rt, rho*u*v, rho*(gogm1*rt+ke)+rho*v*v, -rho*v*(gogm1*rt+ke)/rt+rho*v*gogm1;
//	
//	temp = 0.0;
//	for(int i=0; i<NV; ++i)
//		for(int j=0; j<NV; ++j)
//			for(int k=0; k<NV; ++k)
//				temp(i,j)+=P(i,k)*B(k,j);
//	B = temp;
//	matrix_absolute_value(B);
	
	
	FLT temp2 = sqrt(v*v*(1.0-2.0*b2+b2*b2)+4.0*b2*c2);
	
	V = 0.0, 0.0, 0.5*(v*(b2-1.0)+temp2)*rho, 0.5*(v*(b2-1.0)-temp2)*rho,
		0.0, 1.0, 0.0, 0.0,
		0.0, 0.0, 1.0, 1.0,
		1.0, 0.0, 0.5*(gm1*v*(b2-1.0)+gm1*temp2)/gam, 0.5*(gm1*v*(b2-1.0)-gm1*temp2)/gam;
	
	Beigs = v, v, 0.5*(v+v*b2+temp2), 0.5*(v+v*b2-temp2);
	
	for(int i=0; i<NV; ++i)
		Beigs(i) = abs(Beigs(i));

	VINV = -gm1/(gam*rho), 0.0, 0.0, 1.0,
		   0.0, 1.0, 0.0, 0.0,
		   1.0/(temp2*rho), 0.0, 0.5*(-v*(b2-1.0)+temp2)/temp2, 0.0,
		   -1.0/(rho*temp2), 0.0, 0.5*(v*(b2-1.0)+temp2)/temp2, 0.0;

	for(int i=0; i < NV; ++i)
		for(int j=0; j < NV; ++j)
			VINV(i,j) = Beigs(i)*VINV(i,j);
	
	B = 0.0;
	for(int i=0; i<NV; ++i)
		for(int j=0; j<NV; ++j)
			for(int k=0; k<NV; ++k)
				B(i,j)+=V(i,k)*VINV(k,j);
	
	S = 0.0, 0.0,      0.0,      0.0,
		0.0, nu/(h*h), 0.0,      0.0,
		0.0, 0.0,      nu/(h*h), 0.0,
		0.0, 0.0,      0.0,      alpha/(h*h);
	
	temp = 0.0;
	for(int i=0; i<NV; ++i)
		for(int j=0; j<NV; ++j)
			for(int k=0; k<NV; ++k)
				temp(i,j)+=P(i,k)*S(k,j);
	S = temp;
	
	for(int i=0; i<NV; ++i)
		S(i,i) += gbl->bd(0);
	
	Tinv = 2.0/h*(A+B+h*S);
		
	/* smallest eigenvalue of Tau tilde */
	timestep = 1.0/spectral_radius(Tinv);

	/*  LU factorization  */
	int info,ipiv[NV];
	GETRF(NV, NV, Tinv.data(), NV, ipiv, info);
	
	if (info != 0) {
		cout << "P " << P << endl;
		*gbl->log << "DGETRF FAILED FOR CNS TSTEP" << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
	
	for (int i = 0; i < NV; ++i)
		for (int j = 0; j < NV; ++j)
			temp(i,j)=P(j,i);
	
	/* Solve transposed system temp' = inv(Tinv')*temp' */
	char trans[] = "T";
	GETRS(trans,NV,NV,Tinv.data(),NV,ipiv,temp.data(),NV,info);
	
	if (info != 0) {
		*gbl->log << "DGETRS FAILED FOR CNS TSTEP" << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
	
	for (int i = 0; i < NV; ++i)
		for (int j = 0; j < NV; ++j)
			Tau(i,j)=temp(j,i);

//	/* jacobi smoothing need to modify P too? */
//	temp = 0.0;
//	for (int i = 0; i < NV; ++i)
//		for (int j = 0; j < NV; ++j)
//			for (int k = 0; k < NV; ++k)
//				temp(i,j) += Pinv(i,k)*(A(k,j)/h+B(k,j)/h);
//				
//	Pinv = temp;
	
	return;
}



