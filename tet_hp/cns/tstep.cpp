//#include <utilities.h>
//#include <myblas.h>

#include "tet_hp_cns.h"
#include "../hp_boundary.h"
#include<blitz/tinyvec-et.h>

//#define TIMEACCURATE

void tet_hp_cns::setup_preconditioner() {
	/* SET-UP PRECONDITIONER */
	int gpx = basis::tet(log2p).gpx;
	int gpy = basis::tet(log2p).gpy;
	int gpz = basis::tet(log2p).gpz;
	int stridey = MXGP;
	int stridex = MXGP*MXGP;
	FLT jcb,hmin,hmax,qmax,tstep,dtstari;
	TinyVector<FLT,ND> mvel;
	TinyVector<int,4> v;
	TinyVector<int,3> vtri;
	TinyVector<double,3>vec1,vec2,vec3;

	Array<double,1> umax(NV),ubar(NV);
	Array<double,2> tprcn(NV,NV),tau(NV,NV);
	
	/***************************************/
	/** DETERMINE FLOW PSEUDO-TIME STEP ****/
	/***************************************/
	Array<double,3> tpreconditioner(ntet,NV,NV);
	
	gbl->vpreconditioner(Range(0,npnt-1),Range::all(),Range::all()) = 0.0;
	gbl->epreconditioner(Range(0,nseg-1),Range::all(),Range::all()) = 0.0;

	gbl->vprcn(Range(0,npnt-1),Range::all()) = 0.0;
	if (basis::tet(log2p).em > 0) {
		gbl->eprcn(Range(0,nseg-1),Range::all()) = 0.0;	
		if (basis::tet(log2p).fm > 0) {
			gbl->fprcn(Range(0,ntri-1),Range::all()) = 0.0;
		}
	}

	
#ifdef TIMEACCURATE
	FLT dtstarimax = 0.0;
#endif
	
	for(int tind = 0; tind < ntet; ++tind) {
		jcb = 0.125*tet(tind).vol;   /* volume is 8 x tet volume */
		v = tet(tind).pnt;
		
		ugtouht(tind);
		for(int n=0;n<NV;++n)
			basis::tet(log2p).proj(&uht(n)(0),&u(n)(0)(0)(0),stridex,stridey);

		umax = 0.0;
		ubar = 0.0;
		qmax = 0.0;

		for(int i=0;i<gpx;++i) {
			for(int j=0;j<gpy;++j) {
				for(int k=0;k<gpz;++k) {
				
					mvel(0) = 0.0;//gbl->bd(0)*(crd(0)(i,j) -dxdt(log2p,tind,0)(i)(j)(k));
					mvel(1) = 0.0;//gbl->bd(0)*(crd(1)(i,j) -dxdt(log2p,tind,1)(i)(j)(k));                       
					mvel(2) = 0.0;//gbl->bd(0)*(crd(2)(i,j) -dxdt(log2p,tind,2)(i)(j)(k)); 
						
					FLT q = pow(u(1)(i)(j)(k)-0.5*mvel(0),2.0)  +pow(u(2)(i)(j)(k)-0.5*mvel(1),2.0) +pow(u(3)(i)(j)(k)-0.5*mvel(2),2.0);
					qmax = MAX(qmax,q);
					
					umax(0) = MAX(umax(0),fabs(u(0)(i)(j)(k)));
					umax(1) = MAX(umax(1),fabs(u(1)(i)(j)(k)-0.5*mvel(0)));
					umax(2) = MAX(umax(2),fabs(u(2)(i)(j)(k)-0.5*mvel(1)));
					umax(3) = MAX(umax(3),fabs(u(3)(i)(j)(k)-0.5*mvel(2)));
					umax(4) = MAX(umax(4),fabs(u(4)(i)(j)(k)));
					
					for(int n=0;n<NV;++n)
						ubar(n) += u(n)(i)(j)(k);
					}				
			}
		}

		ubar /= gpx*gpy*gpz;					

		FLT amin = 1.0e99;
		FLT amax = 0.0;
		for(int j=0;j<4;++j) { /* FIND MAX FACE AREA AND THEN DIVIDE VOLUME BY IT */
			vtri = tri(tet(tind).tri(j)).pnt;
			vec1 = pnts(vtri(0))-pnts(vtri(1));
			vec2 = pnts(vtri(0))-pnts(vtri(2));
			vec3 = cross(vec1,vec2);
			FLT a = 0.5*sqrt(vec3(0)*vec3(0)+vec3(1)*vec3(1)+vec3(2)*vec3(2));
			amax = (a > amax ? a : amax);
			amin = (a < amin ? a : amin);
		}
		
		hmin = 4.0*jcb/(0.25*(basis::tet(log2p).p+1)*(basis::tet(log2p).p+1)*amax); /* 3*8/6=4 */
		hmax = 4.0*jcb/(0.25*(basis::tet(log2p).p+1)*(basis::tet(log2p).p+1)*amin); /* 3*8/6=4 */
		
		if (!(hmin > 0.0)) { 
			*gbl->log << "negative tetrahedral area caught in tstep. Problem tet is : " << tind << std::endl;
			*gbl->log << "approximate location: " << pnts(v(0))(0) << ' ' << pnts(v(0))(1) << ' ' << pnts(v(0))(2) << std::endl;
			tet_mesh::output("negative",grid);
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
		
		if  (std::isnan(qmax)) { 
			*gbl->log << "flow solution has nan's" << std::endl;
			output("nan",tecplot);
			sim::abort(__LINE__,__FILE__,gbl->log);
		}		
		
		pennsylvania_peanut_butter(umax,hmin,tprcn,tau,tstep);

		dtstari = 1.0/tstep;

		tpreconditioner(tind,Range::all(),Range::all()) = tprcn;

		gbl->tau(tind,Range::all(),Range::all()) = adis*tau/jcb;
		
#ifdef TIMEACCURATE
		dtstarimax = MAX(dtstari,dtstarimax);
		
	}
	
	/* find max dtstari for all blocks and use on every block  */
	FLT dtstari_recv;
	sim::blks.allreduce(&dtstarimax,&dtstari_recv,1,blocks::flt_msg,blocks::max);
	dtstarimax = dtstari_recv;
	
	*gbl->log << "#iterative to physical time step ratio: " << gbl->bd(0)/dtstarimax << ' ' << gbl->bd(0) << ' ' << dtstarimax << '\n';
	
	
	for(int tind=0;tind<ntet;++tind) {
		v = tet(tind).pnt;
		jcb = 0.125*tet(tind).vol;   /* volume is 8 x tet volume */
		dtstari = dtstarimax;
#endif
		
		gbl->iprcn(tind,Range::all()) = jcb;
		tpreconditioner(tind,Range::all(),Range::all()) *= dtstari;
		
		for(int i=0;i<4;++i) {
			
			gbl->vprcn(v(i),Range::all())  += gbl->iprcn(tind,Range::all());
			gbl->vpreconditioner(v(i),Range::all(),Range::all())  += tpreconditioner(tind,Range::all(),Range::all())*jcb;
			
			if (basis::tet(log2p).fm > 0) {
				gbl->fprcn(tet(tind).tri(i),Range::all()) += gbl->iprcn(tind,Range::all());
			}
		}

		if (basis::tet(log2p).em > 0) {
			for(int i=0;i<6;++i) {
				gbl->eprcn(tet(tind).seg(i),Range::all()) += gbl->iprcn(tind,Range::all());
				gbl->epreconditioner(tet(tind).seg(i),Range::all(),Range::all())  += tpreconditioner(tind,Range::all(),Range::all())*jcb;
			}
		}

	}
	
	for(int i=0;i<npnt;++i) {
		gbl->vpreconditioner(i,Range::all(),Range::all()) /=  gbl->vprcn(i,0);
	}
	for(int i=0;i<nseg;++i) {
		gbl->epreconditioner(i,Range::all(),Range::all()) /=  gbl->eprcn(i,0);
	}
	
	// note temp: remember to do parallel communication
	
	tet_hp::setup_preconditioner();
	
	int last_phase,mp_phase;
	
	if(gbl->diagonal_preconditioner){
		for(int stage = 0; stage<NV; ++stage) {
			for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
				pc0load(mp_phase,gbl->vpreconditioner.data() +stage*NV,NV);
				pmsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
				last_phase = true;
				last_phase &= pc0wait_rcv(mp_phase,gbl->vpreconditioner.data()+stage*NV,NV);
			}
		}
	}
	
}

void tet_hp_cns::pennsylvania_peanut_butter(Array<double,1> pvu, FLT h, Array<FLT,2> &Pinv, Array<FLT,2> &Tau, FLT &timestep) {
	
	Array<double,2> P(NV,NV), V(NV,NV), VINV(NV,NV), dpdc(NV,NV), dcdp(NV,NV), A(NV,NV), B(NV,NV), C(NV,NV), S(NV,NV), Tinv(NV,NV), temp(NV,NV);
	Array<FLT,1> Aeigs(NV),Beigs(NV),Ceigs(NV);

	FLT gam = gbl->gamma;
	FLT gm1 = gam-1.0;
	FLT gogm1 = gam/gm1;
	FLT pr = pvu(0);
	FLT u = pvu(1);
	FLT v = pvu(2);
	FLT w = pvu(3);
	FLT rt = pvu(4);
	FLT rho = pr/rt;
	FLT ke = 0.5*(u*u+v*v+w*w);
	FLT c2 = gam*rt;
	FLT c = sqrt(c2);
	
	FLT nu = gbl->mu/rho;
	FLT cp = gogm1*gbl->R;
	FLT alpha = gbl->kcond/(rho*cp);
	
	FLT hdt = 0.5*h*gbl->bd(0)/c;
	FLT umag = sqrt(u*u+v*v+w*w);
	FLT M = MAX(umag/c,1.0e-5);
	FLT nuh = 4.0*nu/(h*c);
	FLT alh = 2.0*alpha/(h*c);//maybe it should be smaller?
	
	FLT b2 = MIN(M*M/(1.0-M*M) + hdt*hdt + nuh*nuh + alh*alh,1.0);
	FLT alph = 1.0+b2;
	//cout  << b2 << ' ' <<  M*M << ' ' << M*M/(1.0-M*M) << ' ' << hdt*hdt << ' ' << nuh*nuh << ' ' << alh*alh << endl;

	alph = 0.0; // prevents wiggles when residual gets small, not sure why
	//b2 = 1.0; // turn off preconditioner for now
	
	/* Preconditioner */
	P = b2,					  0.0, 0.0, 0.0, 0.0,
		-alph*u/(pr*gam),     1.0, 0.0, 0.0, 0.0,
		-alph*v/(pr*gam),     0.0, 1.0, 0.0, 0.0,
		-alph*w/(pr*gam),     0.0, 0.0, 1.0, 0.0,
		(b2-1.0)/(gogm1*rho), 0.0, 0.0, 0.0, 1.0;
	
	/* Inverse of Preconditioner */
	Pinv = 1.0/b2,					 0.0, 0.0, 0.0, 0.0,
		   alph*u/(pr*gam*b2),		 1.0, 0.0, 0.0, 0.0,
		   alph*v/(pr*gam*b2),		 0.0, 1.0, 0.0, 0.0,
		   alph*w/(pr*gam*b2),		 0.0, 0.0, 1.0, 0.0,
		   -(b2-1.0)/(gogm1*rho*b2), 0.0, 0.0, 0.0, 1.0;
	
	/* jacobian of primitive wrt conservative */
	dpdc = ke*gm1,          -u*gm1,     -v*gm1,    -w*gm1,      gm1,
		   -u/rho,          1.0/rho,    0.0,        0.0,        0.0,
		   -v/rho,          0.0,        1.0/rho,    0.0,        0.0,
		   -w/rho,          0.0,        0.0,		1.0/rho,    0.0,
		   (gm1*ke-rt)/rho, -u*gm1/rho, -v*gm1/rho, -w*gm1/rho, gm1/rho;	
	
	/* jacobian of primitive wrt conservative */
	dcdp = 1.0/rt,               0.0,   0.0,   0.0,   -rho/rt,
		   u/rt,                 rho,   0.0,   0.0,   -rho*u/rt,
		   v/rt,                 0.0,   rho,   0.0,   -rho*v/rt,
		   w/rt,                 0.0,   0.0,   rho,   -rho*w/rt,
	       (rt+gm1*ke)/(gm1*rt), rho*u, rho*v, rho*w, -rho*ke/rt;	
	
	temp = 0.0;
	for(int i=0; i<NV; ++i)
		for(int j=0; j<NV; ++j)
			for(int k=0; k<NV; ++k)
				temp(i,j)+=P(i,k)*dpdc(k,j);
	P = temp;
	
	
	FLT temp1 = sqrt(u*u*(1.0-2.0*b2+b2*b2)+4.0*b2*c2);

	V = 0.5*(u*(b2-1.0)+temp1)*rho,			0.5*(u*(b2-1.0)-temp1)*rho,			0.0, 0.0, 0.0,
	    1.0,								1.0,								0.0, 0.0, 0.0,
	    0.0,								0.0,								1.0, 0.0, 0.0,
	    0.0,								0.0,								0.0, 1.0, 0.0,
		0.5*(u*(b2-1.0)*gm1+gm1*temp1)/gam, 0.5*(u*(b2-1.0)*gm1-gm1*temp1)/gam, 0.0, 0.0, 1.0;
	
	Aeigs = 0.5*(u+u*b2+temp1), 0.5*(u+u*b2-temp1),u,u,u;

	for(int i=0; i<NV; ++i)
		Aeigs(i) = abs(Aeigs(i));
	
	VINV = 1.0/(rho*temp1),	 0.5*(temp1-u*(b2-1.0))/temp1, 0.0, 0.0, 0.0,
		   -1.0/(rho*temp1), 0.5*(u*(b2-1.0)+temp1)/temp1, 0.0, 0.0, 0.0,
		   0.0,				 0.0,						   1.0, 0.0, 0.0,
		   0.0,				 0.0,						   0.0, 1.0, 0.0,
		   -gm1/(gam*rho),	 0.0,						   0.0, 0.0, 1.0;
	
	for(int i=0; i < NV; ++i)
		for(int j=0; j < NV; ++j)
			VINV(i,j) = Aeigs(i)*VINV(i,j);
	
	A = 0.0;
	for(int i=0; i<NV; ++i)
		for(int j=0; j<NV; ++j)
			for(int k=0; k<NV; ++k)
				A(i,j)+=V(i,k)*VINV(k,j);
	
	FLT temp2 = sqrt(v*v*(1.0-2.0*b2+b2*b2)+4.0*b2*c2);
	
	V = 0.5*(v*(b2-1.0)+temp2)*rho,		    0.5*(v*(b2-1.0)-temp2)*rho,			0.0, 0.0, 0.0, 
		0.0,							    0.0,								1.0, 0.0, 0.0,
		1.0,							    1.0,								0.0, 0.0, 0.0,
		0.0,							    0.0,								0.0, 1.0, 0.0,
		0.5*(gm1*v*(b2-1.0)+gm1*temp2)/gam, 0.5*(gm1*v*(b2-1.0)-gm1*temp2)/gam, 0.0, 0.0, 1.0;
	
	Beigs = 0.5*(v+v*b2+temp2), 0.5*(v+v*b2-temp2), v, v, v;

	for(int i=0; i<NV; ++i)
		Beigs(i) = abs(Beigs(i));
	
	VINV = 	1.0/(rho*temp2),  0.0, 0.5*(temp2-v*(b2-1.0))/temp2, 0.0, 0.0,
			-1.0/(rho*temp2), 0.0, 0.5*(v*(b2-1.0)+temp2)/temp2, 0.0, 0.0,
			0.0,			  1.0, 0.0,						     0.0, 0.0,
			0.0,			  0.0, 0.0,							 1.0, 0.0,
			-gm1/(gam*rho),	  0.0, 0.0,							 0.0, 1.0;

	for(int i=0; i < NV; ++i)
		for(int j=0; j < NV; ++j)
			VINV(i,j) = Beigs(i)*VINV(i,j);
	
	B = 0.0;
	for(int i=0; i<NV; ++i)
		for(int j=0; j<NV; ++j)
			for(int k=0; k<NV; ++k)
				B(i,j)+=V(i,k)*VINV(k,j);

	FLT temp3 = sqrt(w*w*(1.0-2.0*b2+b2*b2)+4.0*b2*c2);
	
	V = 0.5*(w*(b2-1.0)+temp3)*rho,			0.5*(w*(b2-1.0)-temp3)*rho,			0.0, 0.0, 0.0,
		0.0,								0.0,								1.0, 0.0, 0.0,
		0.0,								0.0,								0.0, 1.0, 0.0,
		1.0,								1.0,								0.0, 0.0, 0.0,
		0.5*(gm1*w*(b2-1.0)+gm1*temp3)/gam, 0.5*(gm1*w*(b2-1.0)-gm1*temp2)/gam, 0.0, 0.0, 1.0;
	
	Ceigs = 0.5*(w+w*b2+temp3), 0.5*(w+w*b2-temp3), w, w, w;

	for(int i=0; i<NV; ++i)
		Ceigs(i) = abs(Ceigs(i));

	VINV = 1.0/(rho*temp3),  0.0, 0.0, 0.5*(temp3-w*(b2-1.0))/temp3, 0.0,
		   -1.0/(rho*temp3), 0.0, 0.0, 0.5*(w*(b2-1.0)+temp3)/temp3, 0.0,
		   0.0,				 1.0, 0.0, 0.0,							 0.0,
		   0.0,				 0.0, 1.0, 0.0,							 0.0,
		   -gm1/(rho*gam),	 0.0, 0.0, 0.0,							 1.0;
	
	
	for(int i=0; i < NV; ++i)
		for(int j=0; j < NV; ++j)
			VINV(i,j) = Ceigs(i)*VINV(i,j);
	
	C = 0.0;
	for(int i=0; i<NV; ++i)
		for(int j=0; j<NV; ++j)
			for(int k=0; k<NV; ++k)
				C(i,j)+=V(i,k)*VINV(k,j);

	S = 0.0, 0.0,      0.0,      0.0,      0.0,
		0.0, nu/(h*h), 0.0,      0.0,      0.0,
		0.0, 0.0,      nu/(h*h), 0.0,      0.0,
		0.0, 0.0,      0.0,      nu/(h*h), 0.0,
		0.0, 0.0,      0.0,      0.0,      alpha/(h*h);
	
	temp = 0.0;
	for(int i=0; i<NV; ++i)
		for(int j=0; j<NV; ++j)
			for(int k=0; k<NV; ++k)
				temp(i,j)+=P(i,k)*S(k,j);
	S = temp;
	
	for(int i=0; i<NV; ++i)
		S(i,i) += gbl->bd(0);

	Tinv = 2.0/h*(A+B+C+h*S);

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

	return;
}



