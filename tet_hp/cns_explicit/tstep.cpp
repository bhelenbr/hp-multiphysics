//#include <utilities.h>
//#include <myblas.h>

#include "tet_hp_cns_explicit.h"
#include "../hp_boundary.h"
#include<blitz/tinyvec-et.h>

//#define TIMEACCURATE

void tet_hp_cns_explicit::setup_preconditioner() {
	/* SET-UP PRECONDITIONER */
	int gpx = basis::tet(log2p).gpx;
	int gpy = basis::tet(log2p).gpy;
	int gpz = basis::tet(log2p).gpz;
	int stridey = MXGP;
	int stridex = MXGP*MXGP;
	FLT jcb,hmin,hmax,tstep;
	TinyVector<FLT,ND> mvel;
	TinyVector<int,4> v;
	TinyVector<int,3> vtri;
	TinyVector<double,3> vec1,vec2,vec3;

	Array<double,1> umax(NV),ubar(NV);
	Array<double,2> tprcn(NV,NV),tau(NV,NV);
	
	/***************************************/
	/** DETERMINE FLOW PSEUDO-TIME STEP ****/
	/***************************************/
	gbl->vprcn(Range(0,npnt-1),Range::all()) = 0.0;
	if (basis::tet(log2p).em > 0) {
		gbl->eprcn(Range(0,nseg-1),Range::all()) = 0.0;	
		if (basis::tet(log2p).fm > 0) {
			gbl->fprcn(Range(0,ntri-1),Range::all()) = 0.0;
		}
	}

	
#ifdef TIMEACCURATE
	FLT dtstari = 0.0;
#endif
	
	for(int tind = 0; tind < ntet; ++tind) {
		jcb = 0.125*tet(tind).vol;   /* volume is 8 x tet volume */
		v = tet(tind).pnt;
		
		ugtouht(tind);
		for(int n=0;n<NV;++n)
			basis::tet(log2p).proj(&uht(n)(0),&u(n)(0)(0)(0),stridex,stridey);

		umax = 0.0;
		ubar = 0.0;

		for(int i=0;i<gpx;++i) {
			for(int j=0;j<gpy;++j) {
				for(int k=0;k<gpz;++k) {
					for(int n=0;n<NV;++n){
						umax(n) = MAX(umax(n),fabs(u(n)(i)(j)(k)));
						ubar(n) += u(n)(i)(j)(k);
					}
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
		
		if  (std::isnan(umax(1)+umax(2)+umax(3))) { 
			*gbl->log << "flow solution has nan's" << std::endl;
			output("nan",tecplot);
			sim::abort(__LINE__,__FILE__,gbl->log);
		}		
		//cout << hmin << ' ' << hmax << endl;
		pennsylvania_peanut_butter(umax,hmin,tprcn,tau,tstep);
		tstep*=.1;
		gbl->tau(tind,Range::all(),Range::all()) = adis*tau/jcb;
			 
	    jcb /= tstep;
			 
#ifdef TIMEACCURATE
		dtstari = MAX(1.0/tstep,dtstari);
		
	}
	
	/* find max dtstari for all blocks and use on every block  */
	FLT dtstari_recv;
	sim::blks.allreduce(&dtstari,&dtstari_recv,1,blocks::flt_msg,blocks::max);
	dtstari = dtstari_recv;
	
	*gbl->log << "#iterative to physical time step ratio: " << gbl->bd(0)/dtstari << ' ' << gbl->bd(0) << ' ' << dtstari << '\n';
	
	
	for(int tind=0;tind<ntet;++tind) {
		v = tet(tind).pnt;
		jcb = 0.125*tet(tind).vol*dtstari;   /* volume is 8 x tet volume */
#endif
		
		gbl->iprcn(tind,Range::all()) = jcb;
		
		for(int i=0;i<4;++i) {
			
			gbl->vprcn(v(i),Range::all())  += gbl->iprcn(tind,Range::all());
			
			if (basis::tet(log2p).fm > 0) {
				gbl->fprcn(tet(tind).tri(i),Range::all()) += gbl->iprcn(tind,Range::all());
			}
		}

		if (basis::tet(log2p).em > 0) {
			for(int i=0;i<6;++i) {
				gbl->eprcn(tet(tind).seg(i),Range::all()) += gbl->iprcn(tind,Range::all());
			}
		}

	}
	
	tet_hp::setup_preconditioner();

}

void tet_hp_cns_explicit::pennsylvania_peanut_butter(Array<double,1> cvu, FLT h, Array<FLT,2> &Pinv, Array<FLT,2> &Tau, FLT &timestep) {
	
	Array<double,2> P(NV,NV), V(NV,NV), VINV(NV,NV), A(NV,NV), B(NV,NV), C(NV,NV), S(NV,NV), Tinv(NV,NV), temp(NV,NV);
	Array<FLT,1> eigs(NV);

	FLT gam = gbl->gamma;
	FLT gm1 = gam-1.0;
	FLT gogm1 = gam/gm1;
	
	FLT rho = cvu(0);
	FLT rhou = cvu(1);
	FLT rhov = cvu(2);
	FLT rhow = cvu(3);
	FLT rhoE = cvu(4);
	
	FLT u = rhou/rho;
	FLT v = rhov/rho;
	FLT w = rhow/rho;
	FLT E = rhoE/rho;	
	FLT ke = 0.5*(u*u+v*v+w*w);
	FLT c2 = gam*gm1*(E-ke);
	FLT c = sqrt(c2);
	
	FLT nu = gbl->mu/rho;
	FLT cp = gogm1*gbl->R;
	FLT alpha = gbl->kcond/(rho*cp);
	
	/* Preconditioner */
	P = 1.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 1.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 1.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 1.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 1.0;
	
	/* Inverse of Preconditioner */
	Pinv = 1.0, 0.0, 0.0, 0.0, 0.0,
		   0.0, 1.0, 0.0, 0.0, 0.0,
		   0.0, 0.0, 1.0, 0.0, 0.0,
		   0.0, 0.0, 0.0, 1.0, 0.0,
		   0.0, 0.0, 0.0, 0.0, 1.0;

	/* eigenvectors of df_x/dw */
	V = 0.0, 1.0,    0.0, 1.0,                     1.0,
		0.0, u,      0.0, u+c,                     u-c,
		0.0, 0.0,    1.0, v,                       v,
		1.0, 0.0,    0.0, w,                       w,
		w,   u*u-ke, v,   (gm1*u*c+gm1*ke+c2)/gm1, (-gm1*u*c+gm1*ke+c2)/gm1;
		
	eigs = u,u,u,u+c,abs(u-c);

	VINV = -w*ke*gm1,        u*w*gm1,        v*w*gm1,    gm1*w*w+c2, -w*gm1,
	       c2-gm1*ke,        u*gm1,          v*gm1,      w*gm1,      -gm1,
	       -v*ke*gm1,        u*v*gm1,        gm1*v*v+c2, v*w*gm1,    -v*gm1,
	       0.5*(gm1*ke-u*c), -0.5*(u*gm1-c), -0.5*v*gm1, -0.5*w*gm1, 0.5*gm1,
		   0.5*(gm1*ke+u*c), -0.5*(u*gm1+c), -0.5*v*gm1, -0.5*w*gm1, 0.5*gm1;
	
	VINV /= c2;
	
	for(int i=0; i < NV; ++i)
		for(int j=0; j < NV; ++j)
			VINV(i,j) = eigs(i)*VINV(i,j);
	
	A = 0.0;
	for(int i=0; i<NV; ++i)
		for(int j=0; j<NV; ++j)
			for(int k=0; k<NV; ++k)
				A(i,j)+=V(i,k)*VINV(k,j);

	V =  0.0, 1.0,    0.0, 1.0,                     1.0,
		 1.0, 0.0,    0.0, u,                       u,
		 0.0, v,      0.0, v+c,                     v-c,
         0.0, 0.0,    1.0, w,                       w,
		 u,   v*v-ke, w,   (gm1*v*c+gm1*ke+c2)/gm1, (-gm1*v*c+gm1*ke+c2)/gm1;
	
	eigs = v,v,v,v+c,abs(v-c);
	
	VINV = -u*ke*gm1,        c2+gm1*u*u, u*v*gm1,        gm1*w*u,    -gm1*u,
	       c2-gm1*ke,        gm1*u,      gm1*v,          w*gm1,      -gm1,
		   -w*ke*gm1,        gm1*w*u,    gm1*v*w,        c2+gm1*w*w, -w*gm1,
	       0.5*(gm1*ke-v*c), -0.5*gm1*u, -0.5*(v*gm1-c), -0.5*w*gm1, 0.5*gm1,
		   0.5*(gm1*ke+v*c), -0.5*gm1*u, -0.5*(v*gm1+c), -0.5*w*gm1, 0.5*gm1;
	
	VINV /= c2;

	for(int i=0; i < NV; ++i)
		for(int j=0; j < NV; ++j)
			VINV(i,j) = eigs(i)*VINV(i,j);
	
	B = 0.0;
	for(int i=0; i<NV; ++i)
		for(int j=0; j<NV; ++j)
			for(int k=0; k<NV; ++k)
				B(i,j)+=V(i,k)*VINV(k,j);

	V = 0.0, 0.0, 1.0,    1.0,                     1.0,
	    1.0, 0.0, 0.0,    u,                       u,
        0.0, 1.0, 0.0,    v,                       v,
        0.0, 0.0, w,      w+c,                     w-c,
	    u,   v,   w*w-ke, (gm1*w*c+gm1*ke+c2)/gm1, (-gm1*w*c+gm1*ke+c2)/gm1;
	
	eigs = w,w,w,w+c,abs(w-c);
	
	VINV = -u*ke*gm1,        c2+gm1*u*u, gm1*u*v,    gm1*u*w,        -(gam-1)*u,
	       -v*ke*gm1,        gm1*u*v,    c2+gm1*v*v, gm1*v*w,        -gm1*v,
	       c2-gm1*ke,        gm1*u,      gm1*v,      gm1*w,          -gm1,
	       0.5*(gm1*ke-w*c), -0.5*gm1*u, -0.5*gm1*v, -0.5*(w*gm1-c), 0.5*gm1,
		   0.5*(gm1*ke+w*c), -0.5*gm1*u, -0.5*gm1*v, -0.5*(w*gm1+c), 0.5*gm1;
	
	VINV /= c2;

	for(int i=0; i < NV; ++i)
		for(int j=0; j < NV; ++j)
			VINV(i,j) = eigs(i)*VINV(i,j);
	
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
	
	for(int i=0; i<NV; ++i)
		S(i,i) += gbl->bd(0);
	
	Tinv = 2.0/h*(A+B+C+h*S);
	
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



