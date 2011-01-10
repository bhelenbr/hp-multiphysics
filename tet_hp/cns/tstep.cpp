#include <math.h>

#include "tet_hp_cns.h"
#include "../hp_boundary.h"
#include<blitz/tinyvec-et.h>
// temp fix entire routine
void tet_hp_cns::setup_preconditioner() {
	if (gbl->diagonal_preconditioner) {
		/* SET-UP DIAGONAL PRECONDITIONER */
		int tind,i,j,side,find;
		FLT jcb,h,hmax,q,qmax,lam1,gam,a,amax,amin;
		TinyVector<int,4> v;
		TinyVector<int,3> vtri;

		TinyVector<double,3>vec1,vec2,vec3;
		FLT havg = 0.0;
		FLT nu = gbl->mu/gbl->rho;

		/***************************************/
		/** DETERMINE FLOW PSEUDO-TIME STEP ****/
		/***************************************/
		gbl->vprcn(Range::all(),Range::all())=0.0;
		if (basis::tet(log2p).em > 0) {
			gbl->eprcn(Range::all(),Range::all())=0.0;
			if (basis::tet(log2p).fm > 0) {
				gbl->fprcn(Range::all(),Range::all())=0.0;
			}
		}
		
#ifdef TIMEACCURATE
		gam = 10.0;
		FLT dtstari = 0.0;
#endif

		for(tind = 0; tind < ntet; ++tind) {
			jcb = 0.125*tet(tind).vol;   /* volume is 8 x tet volume */
			v = tet(tind).pnt;	
			amin=1.0e99;
			amax = 0.0;
			for(j=0;j<4;++j) { /* FIND MAX FACE AREA AND THEN DIVIDE VOLUME BY IT */
				vtri = tri(tet(tind).tri(j)).pnt;
				vec1=pnts(vtri(0))-pnts(vtri(1));
				vec2=pnts(vtri(0))-pnts(vtri(2));
				vec3=cross(vec1,vec2);
				a=.5*sqrt(vec3(0)*vec3(0)+vec3(1)*vec3(1)+vec3(2)*vec3(2));
				amax = (a > amax ? a : amax);
				amin = (a < amin ? a : amin);
			}

			float c = 1.0;
			h = c*4.0*jcb/(0.25*(basis::tet(log2p).p+1)*(basis::tet(log2p).p+1)*amax); /* 3*8/6=4 */
			hmax = c*4.0*jcb/(0.25*(basis::tet(log2p).p+1)*(basis::tet(log2p).p+1)*amin); /* 3*8/6=4 */
//			h = 4.0*jcb/(0.25*(p0+1)*(p0+1)*amax); /* 3*8/6=4 temp FIXME */
//			hmax = 4.0*jcb/(0.25*(p0+1)*(p0+1)*amin); /* 3*8/6=4 temp FIXME */
			
			havg+=hmax;
			qmax = 0.0;
			for(j=0;j<4;++j) {
				q = pow(ug.v(v(j),0)-0.5*(gbl->bd(0)*(pnts(v(j))(0) -vrtxbd(1)(v(j))(0))),2.0) 
					+pow(ug.v(v(j),1)-0.5*(gbl->bd(0)*(pnts(v(j))(1) -vrtxbd(1)(v(j))(1))),2.0)  
					+pow(ug.v(v(j),2)-0.5*(gbl->bd(0)*(pnts(v(j))(2) -vrtxbd(1)(v(j))(2))),2.0);
				qmax = MAX(qmax,q);
			}
			if (!(jcb > 0.0)) { 
				*gbl->log << "negative tetrahedral volume caught in tstep. Problem tet is : " << tind << std::endl;
				*gbl->log << "approximate location: " << pnts(v(0))(0) << ' ' << pnts(v(0))(1)<< ' ' << pnts(v(0))(2) << std::endl;
				tet_mesh::output("negative",grid);
				exit(1);
			}

			if  (!(qmax >= 0.0)) {  // THIS CATCHES NAN'S TOO
				*gbl->log << "flow solution has nan's" << std::endl;
				output("nan",tecplot);
				exit(1);
			}

#ifndef INERTIALESS

#ifndef TIMEACCURATE
			gam = 3.0*qmax +(0.5*hmax*gbl->bd(0) +2.*nu/hmax)*(0.5*hmax*gbl->bd(0) +2.*nu/hmax);
			if (gbl->mu + gbl->bd(0) == 0.0) gam = MAX(gam,0.1);
#endif
			q = sqrt(qmax);
			lam1 = q + sqrt(qmax +gam);

			/* SET UP DISSIPATIVE COEFFICIENTS */
			gbl->tau(tind,0) = adis*h/(jcb*sqrt(gam));
			gbl->tau(tind,NV-1) = qmax*gbl->tau(tind,0);

			/* SET UP DIAGONAL PRECONDITIONER */
			// jcb *= 8.*nu*(1./(hmax*hmax) +1./(h*h)) +2*lam1/h +2*sqrt(gam)/hmax +gbl->bd(0);
			jcb *= 2.*nu*(1./(h*h) +1./(h*h)+1./(h*h)) +3*lam1/h;  // heuristically tuned
#else
			gam = pow(2.*nu/hmax,2); 
			lam1 = sqrt(gam);
			/* SET UP DISSIPATIVE COEFFICIENTS */
			gbl->tau(tind,0)  = adis*h/(jcb*sqrt(gam));
			gbl->tau(tind,NV-1) = 0.0;

			jcb *= 8.*nu*(1./(h*h) +1./(h*h)+1./(h*h)) +2*lam1/h +2*sqrt(gam)/hmax;
#endif
#ifdef TIMEACCURATE
			dtstari = MAX((nu/(h*h) +lam1/h +gbl->bd(0)),dtstari);

		}
		printf("#iterative to physical time step ratio: %f\n",gbl->bd(0)/dtstari);

		for(tind=0;tind<ntet;++tind) {
			v = tet(tind).pnt;
			jcb = 0.125*tet(tind).vol*dtstari;
#endif

			gbl->iprcn(tind,0) = gbl->rho*jcb;    
			gbl->iprcn(tind,1) = gbl->rho*jcb;
			gbl->iprcn(tind,2) = gbl->rho*jcb;      
			gbl->iprcn(tind,3) =  jcb/gam;
			for(i=0;i<4;++i) 
				gbl->vprcn(v(i),Range::all())  += gbl->iprcn(tind,Range::all());            
			if (basis::tet(log2p).em > 0) {
				for(i=0;i<6;++i){
					side = tet(tind).seg(i);
					gbl->eprcn(side,Range::all()) += gbl->iprcn(tind,Range::all());
				}
			}
			if (basis::tet(log2p).fm > 0) {
				for(i=0;i<4;++i){
					side = tet(tind).tri(i);
					gbl->fprcn(side,Range::all()) += gbl->iprcn(tind,Range::all());
				}
			}
		}
		//cout <<"havg = " <<  havg/ntet << endl;
	}
	else {
		cout << "matrix preconditioner being called and doesn't work " << endl;
//        /* SET-UP MATRIX PRECONDITIONER */
//        int tind,i,j,side,v0;
//        FLT jcb,h,hmax,q,qmax,lam1,gam,ubar,vbar;
//        TinyVector<int,4> v;
//        
//        FLT nu = gbl->mu/gbl->rho;
//              
//        /***************************************/
//        /** DETERMINE FLOW PSEUDO-TIME STEP ****/
//        /***************************************/
//        gbl->vprcn_ut(Range(0,npnt-1),Range::all(),Range::all()) = 0.0;
//        if (basis::tet(log2p).em > 0) {
//            gbl->eprcn_ut(Range(0,nseg-1),Range::all(),Range::all()) = 0.0;
//        }
//		if (basis::tet(log2p).fm > 0) {
//            gbl->fprcn_ut(Range(0,ntri-1),Range::all(),Range::all()) = 0.0;
//        }
//        gbl->iprcn_ut(Range(0,ntet-1),Range::all(),Range::all()) = 0.0;
//
//        for(tind = 0; tind < ntet; ++tind) {
//            jcb = tet(tind).vol/8.0; 
//            hmax = 0.0;
//            for(j=0;j<3;++j) {
//                h = pow(pnts(v(j))(0) -pnts(v((j+1)%3))(0),2.0) + 
//                pow(pnts(v(j))(1) -pnts(v((j+1)%3))(1),2.0);
//                hmax = (h > hmax ? h : hmax);
//            }
//            hmax = sqrt(hmax);
//            
//            if (!(jcb > 0.0)) {  // THIS CATCHES NAN'S TOO
//                *gbl->log << "negative triangle area caught in tstep. Problem triangle is : " << tind << std::endl;
//                *gbl->log << "approximate location: " << pnts(v(0))(0) << ' ' << pnts(v(0))(1) << std::endl;
//                tri_mesh::output("negative",grid);
//                exit(1);
//            }
//            h = 4.*jcb/(0.25*(basis::tri(log2p).p +1)*(basis::tri(log2p).p+1)*hmax);
//            hmax = hmax/(0.25*(basis::tri(log2p).p +1)*(basis::tri(log2p).p+1));
//        
//            qmax = 0.0;
//            ubar = 0.0;
//            vbar = 0.0;
//            for(j=0;j<3;++j) {
//                v0 = v(j);
//
//                q = pow(ug.v(v0,0)-(gbl->bd(0)*(pnts(v0)(0) -vrtxbd(1)(v0)(0))),2.0) 
//                    +pow(ug.v(v0,1)-(gbl->bd(0)*(pnts(v0)(1) -vrtxbd(1)(v0)(1))),2.0);  
//                ubar += ug.v(v0,0);
//                vbar += ug.v(v0,1);
//                qmax = MAX(qmax,q);
//            }
//            
//            ubar /= 3.0;
//            vbar /= 3.0;
//
//            gam = 3.0*qmax +(0.5*hmax*gbl->bd(0) +2.*nu/hmax)*(0.5*hmax*gbl->bd(0) +2.*nu/hmax);
//            if (gbl->mu + gbl->bd(0) == 0.0) gam = MAX(gam,0.1);
//
//            q = sqrt(qmax);
//            lam1 = q + sqrt(qmax +gam);
//            
//            /* SET UP DISSIPATIVE COEFFICIENTS */
//            gbl->tau(tind,1) = adis*h/(jcb*sqrt(gam));
//            gbl->tau(tind,0) = qmax*gbl->tau(tind,1);
//            
//            /* SET UP DIAGONAL PRECONDITIONER */
//            // jcb *= 8.*nu*(1./(hmax*hmax) +1./(h*h)) +2*lam1/h +2*sqrt(gam)/hmax +gbl->bd(0);
//            jcb *= 2.*nu*(1./(hmax*hmax) +1./(h*h)) +3*lam1/h;  // heuristically tuned
//            jcb *= RAD((pnts(v(0))(0) +pnts(v(1))(0) +pnts(v(2))(0))/3.);
//
//            gbl->tprcn_ut(tind,0,0) = gbl->rho*jcb;    
//            gbl->tprcn_ut(tind,1,1) = gbl->rho*jcb;      
//            gbl->tprcn_ut(tind,2,2) = jcb/gam;
//            gbl->tprcn_ut(tind,0,2) = jcb*ubar/gam;
//            gbl->tprcn_ut(tind,1,2) = jcb*vbar/gam;
//            for(i=0;i<3;++i) {
//                gbl->vprcn_ut(v(i),Range::all(),Range::all())  += basis::tri(log2p).vdiag*gbl->tprcn_ut(tind,Range::all(),Range::all());
//                if (basis::tri(log2p).sm > 0) {
//                    side = tri(tind).seg(i);
//                    gbl->sprcn_ut(side,Range::all(),Range::all()) += gbl->tprcn_ut(tind,Range::all(),Range::all());
//                }
//            }
//        }
	}


	
	tet_hp::setup_preconditioner();
}

void tet_hp_cns::pennsylvania_peanut_butter(Array<double,1> u, FLT hmax, Array<FLT,2> &Pinv, Array<FLT,2> &Tau, FLT &timestep) {
	
	Array<FLT,2> P(NV,NV),A(NV,NV),B(NV,NV),S(NV,NV),Tinv(NV,NV),temp(NV,NV);
	
	FLT gm1 = gbl->gamma-1;
	FLT gogm1 = gbl->gamma/gm1;
	FLT pr = u(0);
	FLT uv = u(1);
	FLT vv = u(2);
	FLT wv = u(3);
	FLT rt = u(4);	
	FLT ke = 0.5*(uv*uv+vv*vv+wv*wv);
	
	/* Preconditioner */
	P = ke*gm1,            -uv*gm1,           -vv*gm1,           gm1,
	-uv/pr*rt,         1.0/pr*rt,         0.0,               0.0,
	-1.0/pr*vv*rt,     0.0,               1.0/pr*rt,         0.0,
	pr*(gm1*ke-rt)*rt, -1.0/pr*rt*uv*gm1, -1.0/pr*rt*vv*gm1, 1.0/pr*rt*gm1;
	
	
	/* Inverse of Preconditioner */
	Pinv = 1.0/rt,        0.0,      0.0,      -pr/rt/rt,
	uv/rt,         pr/rt,    0.0,      -uv*pr/rt/rt,
	vv/rt,         0.0,      pr/rt,    -vv*pr/rt/rt,
	1.0/gm1+ke/rt, uv*pr/rt, vv*pr/rt, -pr/rt/rt*ke;	
	
	/* df/dw */
	A = 1.0/rt*uv,               pr/rt,                           0.0,         -pr/rt/rt*uv,
	1.0/rt*uv*uv+1.0,        2.0*pr/rt*uv,                    0.0,         -pr/rt/rt*uv*uv,
	1.0/rt*uv*vv,            pr/rt*vv,                        pr/rt*uv,    -pr/rt/rt*uv*vv,
	1.0/rt*uv*(gogm1*rt+ke), pr/rt*(gogm1*rt+ke)+pr/rt*uv*uv, pr/rt*uv*vv, -pr/rt/rt*uv*(gogm1*rt+ke)+pr/rt*uv*gogm1;
	
	temp = 0.0;
	for(int i=0; i<NV; ++i)
		for(int j=0; j<NV; ++j)
			for(int k=0; k<NV; ++k)
				temp(i,j)+=P(i,k)*A(k,j);
	A = temp;
	
	matrix_absolute_value(A);
	
	/* dg/dw */
	B = 1.0/rt*vv,               0.0,         pr/rt,                           -pr/rt/rt*vv,
	1.0/rt*uv*vv,            pr/rt*vv,    pr/rt*uv,                        -pr/rt/rt*uv*vv,
	1.0/rt*vv*vv+1.0,        0.0,         2.0*pr/rt*vv,                    -pr/rt/rt*vv*vv,
	1.0/rt*vv*(gogm1*rt+ke), pr/rt*uv*vv, pr/rt*(gogm1*rt+ke)+pr/rt*vv*vv, -pr/rt/rt*vv*(gogm1*rt+ke)+pr/rt*vv*gogm1;
	
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
