/*
 *  particle.cpp
 *  tet_hp
 *
 *  Created by michael brazell on 2/26/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp.h"
#include "hp_boundary.h"
#include <myblas.h>
#include <iostream>
#include <fstream>
#include <math.h>

void tet_hp::getu(const TinyVector<FLT,3> &xp, int &tind,TinyVector<FLT,4> &uout) {
	int tetsearch,tlvl;
	double rnx,rny,rnz,minvol,volinv,r,s,t;
	TinyVector<int,4> v;	
	TinyVector<FLT,4> vol;	
	TinyMatrix<FLT,5,3> pnt;
	
	for(int i=0;i < 3; ++i){
		pnt(4,i)=xp(i);
	}

	tetsearch=1;
	while(tetsearch > -1){
		tetsearch=-1;		
		v=tet(tind).pnt;// vertices
		for(int i=0;i < 4; ++i){
			for(int j=0;j < 3; ++j){
				pnt(i,j)=pnts(v(i))(j); // load vertex locations
			}
		}

		minvol=0.0;
		for(int i=0;i < 4; ++i){// find volume of tetrahedral created from a face of tetrahedral and point
			v(0)=0;
			v(1)=1;
			v(2)=2;
			v(3)=3;
			v(i)=4;
			rnx = (pnt(v(1),1)-pnt(v(0),1))*(pnt(v(2),2)-pnt(v(0),2))-(pnt(v(2),1)-pnt(v(0),1))*(pnt(v(1),2)-pnt(v(0),2));
			rny = (pnt(v(1),2)-pnt(v(0),2))*(pnt(v(2),0)-pnt(v(0),0))-(pnt(v(2),2)-pnt(v(0),2))*(pnt(v(1),0)-pnt(v(0),0));
			rnz = (pnt(v(1),0)-pnt(v(0),0))*(pnt(v(2),1)-pnt(v(0),1))-(pnt(v(2),0)-pnt(v(0),0))*(pnt(v(1),1)-pnt(v(0),1));
			vol(i)=-rnx*(pnt(v(3),0)-pnt(v(0),0))-rny*(pnt(v(3),1)-pnt(v(0),1))-rnz*(pnt(v(3),2)-pnt(v(0),2)); // calculate volume
			if(vol(i) < minvol){
				minvol=vol(i); // keep track of smallest volume 
				tetsearch=i;
			}
		}
		if(tetsearch > -1){
			tind=tet(tind).tet(tetsearch);// if tet is not found then search in direction of smallest volume
		}
		if(tind < 0) return;
		//cout << tind+1 << ", " ;
	}

	volinv=1.0/tet(tind).vol; // volume of tetrahedral

	/* TETRAHEDRAL COORDINATES */  
	t = 2.0*vol(0)*volinv -1.0;
	s = 2.0*vol(1)*volinv -1.0;
	r = 2.0*vol(3)*volinv -1.0;

	basis::tet(log2p).ptvalues_rst(r,s,t);// value of basis at r,s,t
	tlvl = 0;
	ugtouht(tind,tlvl); // take global coefficients and load them to local coefficients

	basis::tet(log2p).ptprobe(NV,uout.data(),&uht(0)(0),MXTM);	// use basis and local coefficients to interpolate value of function at point	

	return;
}

//void findtet(const TinyVector<FLT,3> &xp, int &seedvrtx, int &tind){
//	
//	int tetsearch;
//	double rnx,rny,rnz,minvol,volinv,r,s,t;
//	TinyVector<int,4> v;	
//	TinyVector<FLT,4> vol;	
//	TinyMatrix<FLT,5,3> pnt;
//	
//	for(int i=0;i < 3; ++i){
//		pnt(4,i)=xp(i);
//	}
//		
//	tind = pnt(seedvertex).tet;
//
//	tetsearch=1;
//	while(tetsearch > -1){
//		tetsearch=-1;		
//		v=tet(tind).pnt;// vertices
//		for(int i=0;i < 4; ++i){
//			for(int j=0;j < 3; ++j){
//				pnt(i,j)=pnts(v(i))(j); // load vertex locations
//			}
//		}
//		
//		minvol=0.0;
//		for(int i = 0; i < 4; ++i){// find volume of tetrahedral created from a face of tetrahedral and point
//			v(0)=0;
//			v(1)=1;
//			v(2)=2;
//			v(3)=3;
//			v(i)=4;
//			rnx = (pnt(v(1),1)-pnt(v(0),1))*(pnt(v(2),2)-pnt(v(0),2))-(pnt(v(2),1)-pnt(v(0),1))*(pnt(v(1),2)-pnt(v(0),2));
//			rny = (pnt(v(1),2)-pnt(v(0),2))*(pnt(v(2),0)-pnt(v(0),0))-(pnt(v(2),2)-pnt(v(0),2))*(pnt(v(1),0)-pnt(v(0),0));
//			rnz = (pnt(v(1),0)-pnt(v(0),0))*(pnt(v(2),1)-pnt(v(0),1))-(pnt(v(2),0)-pnt(v(0),0))*(pnt(v(1),1)-pnt(v(0),1));
//			vol(i)=-rnx*(pnt(v(3),0)-pnt(v(0),0))-rny*(pnt(v(3),1)-pnt(v(0),1))-rnz*(pnt(v(3),2)-pnt(v(0),2)); // calculate volume
//			if(vol(i) < minvol){
//				minvol=vol(i); // keep track of smallest volume 
//				tetsearch=i;
//			}
//		}
//		if(tetsearch > -1){
//			tind=tet(tind).tet(tetsearch);// if tet is not found then search in direction of smallest volume
//		}
//		if(tind < 0){
//			*gbl->log << "Error pnt outside domain: seed vertex = " << seedvrtx " tet = " << tind << " location = " << xp << std::endl;   
//			sim::abort(__LINE__,__FILE__,gbl->log);
//		}
//
//	}
//	
//	return;
//}


void tet_hp::getudx(const TinyVector<FLT,3> &xp, int &tind,TinyVector<FLT,4> &uout,TinyVector<FLT,4> &dudx,TinyVector<FLT,4> &dudy,TinyVector<FLT,4> &dudz){
	
	int tetsearch,tlvl;
	double rnx,rny,rnz,minvol,volinv,r,s,t;
	TinyVector<int,4> v;	
	TinyVector<FLT,4> vol;	
	TinyMatrix<FLT,5,3> pnt;
	
	for(int i=0;i < 3; ++i){
		pnt(4,i)=xp(i);
	}

		
		tetsearch=1;
		while(tetsearch > -1){
			tetsearch=-1;		
			v=tet(tind).pnt;
			for(int i=0;i < 4; ++i){
				for(int j=0;j < 3; ++j){
					pnt(i,j)=pnts(v(i))(j);
				}
			}

			minvol=0.0;
			for(int i=0;i < 4; ++i){
				v(0)=0;
				v(1)=1;
				v(2)=2;
				v(3)=3;
				v(i)=4;
				rnx = (pnt(v(1),1)-pnt(v(0),1))*(pnt(v(2),2)-pnt(v(0),2))-(pnt(v(2),1)-pnt(v(0),1))*(pnt(v(1),2)-pnt(v(0),2));
				rny = (pnt(v(1),2)-pnt(v(0),2))*(pnt(v(2),0)-pnt(v(0),0))-(pnt(v(2),2)-pnt(v(0),2))*(pnt(v(1),0)-pnt(v(0),0));
				rnz = (pnt(v(1),0)-pnt(v(0),0))*(pnt(v(2),1)-pnt(v(0),1))-(pnt(v(2),0)-pnt(v(0),0))*(pnt(v(1),1)-pnt(v(0),1));
				vol(i)=-rnx*(pnt(v(3),0)-pnt(v(0),0))-rny*(pnt(v(3),1)-pnt(v(0),1))-rnz*(pnt(v(3),2)-pnt(v(0),2));
				if(vol(i) < minvol){
					minvol=vol(i);
					tetsearch=i;
				}
			}
			if(tetsearch > -1){
				tind=tet(tind).tet(tetsearch);
			}
			if(tind < 0) return;
			//cout << tind+1 << ", " ;
		}

		volinv=1.0/tet(tind).vol;

		/* TETRAHEDRAL COORDINATES */  
		t = 2.0*vol(0)*volinv -1.0;
		s = 2.0*vol(1)*volinv -1.0;
		r = 2.0*vol(3)*volinv -1.0;

		//basis::tet(log2p).ptvalues_rst(r,s,t);
		tlvl = 0;
		ugtouht(tind,tlvl);  
//		basis::tet(log2p).ptprobe(NV, uout.data(), dudx.data(), dudy.data(), dudz.data(), r, s, t ,&uht(0)(0),MXTM);

	return;
}


void tet_hp::particle() {
	TinyVector<FLT,4> uout,up;
	double dudy,dwdy,lcl;
	TinyVector<FLT,3> x0,xp,g;
	TinyVector<FLT,3> temp,temp2;
	TinyVector<FLT,4> k1,k2,k3;
	TinyVector<FLT,4> k11,k12,k21,k22,k31,k32,k41,k42;
	int tind=30;
	int t;
	
	/* fixed constants */
	double avou=0.2699; // alpha nu over u star
	double mu = 1.84e-5; // viscosity
	double nu = 1.502e-5; // kinematic viscosity
	double q=0.0688; // constant for u velocity
	double Dh=0.03; // Diameter of channel 
	
	/* variable inputs */
	double tau=1.e-2;//relaxation time
	double dt=0.05;//2.5*tau;// time step
	double S=965.0;//density ratio
	double Re=9588.0;//Reynolds number	
	
	x0(0)=6.0,x0(1)=12.0,x0(2)=.000000001;//initial position

	/* dependent variables */
	double C=1.0/tau;// 1/tau
	double U=Re*nu/Dh;//mean flow velocity
	lcl=6.9/Re;
	lcl=-1.8*log10(lcl);
	double ustar=sqrt(U*U/8./pow(lcl,2.));//friction velocity
	double taud=tau*nu/ustar/ustar;// dimensional relaxation time
	double lamb = .07e-6; // mean free path
	double di = .0001; // initial guess at non dimensional diameter
	for(int i = 0; i < 30; ++i){// newtons method to solve for non dimensional diameter
		lcl=-1.1*di/lamb/2;
		di=di-(di*di+1.257*lamb*2.0*di+0.8*di*lamb*exp(lcl)-18.0*nu*nu*tau/S/ustar/ustar)/(2.0*di+lamb*2.514+0.8*lamb*exp(lcl)-0.44*di*exp(lcl));
		//cout << di << endl;
	}
	double Kn = 2.0*lamb/di; // Knudson Number
	lcl = -1.1/Kn;
	double Cc = 1.0+1.257*Kn+0.4*Kn*exp(lcl); // Cunningham correction factor
	double dp=sqrt(18.0*tau/S/Cc);// nondimensional diameter
	double L1=3.08/S/dp;// Lift coefficient for dudy
	double L2=3.08/S/dp; // Lift coefficient for dwdy 
	double gp=-9.807*nu/ustar/ustar/ustar;
	//double di=dp*nu/ustar;// particle diameter
	double D=1.38e-23*(273.15+15.0)*Cc/3.0/3.14/mu/di;// Brownian Diffusivity
	double Sc=nu/D;// Schmidt number
	lcl = 2.0/3.0;
	double ud=2*x0(2)*0.7/100+.084/pow(Sc,lcl); // inertia-interception deposition velocity
	
	cout << "tau          ud"<< endl;
	cout << tau << ' ' << ud << endl;
	cout << "di = " << di << " dp = " << dp << " taud = " << taud << " Kn = " << Kn << " Cc = " << Cc<< " g = " << gp << " ustar = " << ustar << " U = " << U << endl;
	ofstream myfile;
		
	cout.precision(15);
	myfile.open("x.dat");
	myfile.precision(15);	

	g(0)=gp,g(1)=0.0,g(2)=0.0;// gravity
	xp=x0;
	getu(xp, tind, uout);
	up=uout; // initialize particle velocity with fluid velocity	
	tind = 0;
	myfile << xp(0) << ' ' << xp(1) << ' ' << xp(2) << ' ' << tind << ' ' << sqrt(uout(0)*uout(0)+uout(1)*uout(1)+uout(2)*uout(2)) << endl;
	for(t=0; t < 1000000; ++t){
		getu(xp, tind, uout); // calculate fluid velocity at location xp
		dudy=1.0-tanh(q*xp(1))*tanh(q*xp(1));//4/pow((exp(q*xp(1))+exp(-q*xp(1))),2);
		dwdy=avou*avou*avou*xp(2)*(-4.43967969967749e-06*pow(avou*xp(1),10.)
							  +8.65534353334536e-05*pow(avou*xp(1),9.)
							  -6.62974554556394e-04*pow(avou*xp(1),8.)
							  +2.47307460136993e-03*pow(avou*xp(1),7.)
							  -5.42318009049220e-03*pow(avou*xp(1),6.)
							  +1.95264285001060e-02*pow(avou*xp(1),5.)
							  -1.01372917497318e-01*pow(avou*xp(1),4.)
							  +2.46838310150933e-01*pow(avou*xp(1),3.)
							  -2.51943937552724e-04*pow(avou*xp(1),2.)
							  -9.89977577245599e-01*avou*xp(1)
							  +1.23258728815415e+00);
		if(tind < 0) break;
		for(int i = 0; i < 3; ++i){
			k11(i)=up(i);
			temp(i)=xp(i)+0.5*dt*k11(i);
		}
		
		
		k12(0)=C*(uout(0)-up(0))+g(0);
		k12(1)=C*(uout(1)-up(1))+L1*dudy/sqrt(fabs(dudy))*(uout(0)-up(0))+L2*dwdy/sqrt(fabs(dwdy))*(uout(2)-up(2));
		k12(2)=C*(uout(2)-up(2));
		
		for(int i = 0; i < 3; ++i)
			temp2(i)=up(i)+0.5*dt*k12(i);
		
		getu(temp, tind, uout);

		dudy=1.0-tanh(q*temp(1))*tanh(q*temp(1));//4/pow((exp(q*temp(1))+exp(-q*temp(1))),2);
		dwdy=avou*avou*avou*temp(2)*(-4.43967969967749e-06*pow(avou*temp(1),10.)
								+8.65534353334536e-05*pow(avou*temp(1),9.)
								-6.62974554556394e-04*pow(avou*temp(1),8.)
								+2.47307460136993e-03*pow(avou*temp(1),7.)
								-5.42318009049220e-03*pow(avou*temp(1),6.)
								+1.95264285001060e-02*pow(avou*temp(1),5.)
								-1.01372917497318e-01*pow(avou*temp(1),4.)
								+2.46838310150933e-01*pow(avou*temp(1),3.)
								-2.51943937552724e-04*pow(avou*temp(1),2.)
								-9.89977577245599e-01*avou*temp(1)
								+1.23258728815415e+00);
		if(tind < 0) break;
		for(int i = 0; i < 3; ++i){
			k21(i)=temp2(i);
			temp(i)=xp(i)+0.5*dt*k21(i);
		}
		k22(0)=C*(uout(0)-temp2(0))+g(0);
		k22(1)=C*(uout(1)-temp2(1))+L1*dudy/sqrt(fabs(dudy))*(uout(0)-temp2(0))+L2*dwdy/sqrt(fabs(dwdy))*(uout(2)-temp2(2));
		k22(2)=C*(uout(2)-temp2(2));
		
		for(int i = 0; i < 3; ++i)
			temp2(i)=up(i)+0.5*dt*k22(i);		
		
		getu(temp, tind, uout);

		dudy=1.0-tanh(q*temp(1))*tanh(q*temp(1));//4/pow((exp(q*temp(1))+exp(-q*temp(1))),2);
		dwdy=avou*avou*avou*temp(2)*(-4.43967969967749e-06*pow(avou*temp(1),10.)
								+8.65534353334536e-05*pow(avou*temp(1),9.)
								-6.62974554556394e-04*pow(avou*temp(1),8.)
								+2.47307460136993e-03*pow(avou*temp(1),7.)
								-5.42318009049220e-03*pow(avou*temp(1),6.)
								+1.95264285001060e-02*pow(avou*temp(1),5.)
								-1.01372917497318e-01*pow(avou*temp(1),4.)
								+2.46838310150933e-01*pow(avou*temp(1),3.)
								-2.51943937552724e-04*pow(avou*temp(1),2.)
								-9.89977577245599e-01*avou*temp(1)
								+1.23258728815415e+00);
		if(tind < 0) break;
		for(int i = 0; i < 3; ++i){
			k31(i)=temp2(i);
			temp(i)=xp(i)+dt*k31(i);
		}
		k32(0)=C*(uout(0)-temp2(0))+g(0);
		k32(1)=C*(uout(1)-temp2(1))+L1*dudy/sqrt(fabs(dudy))*(uout(0)-temp2(0))+L2*dwdy/sqrt(fabs(dwdy))*(uout(2)-temp2(2));
		k32(2)=C*(uout(2)-temp2(2));
		
		for(int i = 0; i < 3; ++i)
			temp2(i)=up(i)+dt*k32(i);
		
		getu(temp, tind, uout);

		dudy=1.0-tanh(q*temp(1))*tanh(q*temp(1));//4/pow((exp(q*temp(1))+exp(-q*temp(1))),2);
		dwdy=avou*avou*avou*temp(2)*(-4.43967969967749e-06*pow(avou*temp(1),10.)
								+8.65534353334536e-05*pow(avou*temp(1),9.)
								-6.62974554556394e-04*pow(avou*temp(1),8.)
								+2.47307460136993e-03*pow(avou*temp(1),7.)
								-5.42318009049220e-03*pow(avou*temp(1),6.)
								+1.95264285001060e-02*pow(avou*temp(1),5.)
								-1.01372917497318e-01*pow(avou*temp(1),4.)
								+2.46838310150933e-01*pow(avou*temp(1),3.)
								-2.51943937552724e-04*pow(avou*temp(1),2.)
								-9.89977577245599e-01*avou*temp(1)
								+1.23258728815415e+00);
		if(tind < 0) break;
		for(int i = 0; i < 3; ++i)
			k41(i)=temp2(i);
		k42(0)=C*(uout(0)-temp2(0))+g(0);
		k42(1)=C*(uout(1)-temp2(1))+L1*dudy/sqrt(fabs(dudy))*(uout(0)-temp2(0))+L2*dwdy/sqrt(fabs(dwdy))*(uout(2)-temp2(2));
		k42(2)=C*(uout(2)-temp2(2));
		
		/* Runge Kutta 4 */
		for(int i = 0; i < 3; ++i){
			xp(i)=xp(i)+dt/6.0*(k11(i)+2.0*k21(i)+2.0*k31(i)+k41(i));
			up(i)=up(i)+dt/6.0*(k12(i)+2.0*k22(i)+2.0*k32(i)+k42(i));
		}
		xp(0) = x0(0);
		
		myfile << xp(0) << ' ' << xp(1) << ' ' << xp(2) << ' ' << tind << ' ' << sqrt(uout(0)*uout(0)+uout(1)*uout(1)+uout(2)*uout(2)) << endl;
		if(xp(1) < dp/2){ 
			++t;
			break;
		}

		//cout << up(0)<<' ' << up(1) << ' ' << up(2) << endl;
		//cout << tind+1 << ", " ;
	}
	
	//cout << endl;
	if(tind < 0){
		for(int i = 1; i < 3; ++i){
			xp(i)=xp(i)+dt*up(i);
		}
		myfile << xp(0) << ' ' << xp(1) << ' ' << xp(2) << ' ' << tind << ' ' << sqrt(uout(0)*uout(0)+uout(1)*uout(1)+uout(2)*uout(2)) << endl,++t;
	}
	myfile << "0 0 " << t+1 <<  " 0 0" << endl;
	
	if(xp(1) < dp/2 && xp(2) < 25.0) cout << "particle deposited on bottom wall"<< endl;
	
	/* Runge Kutta Streamlines */
	dt=.2;
	xp=x0;	
	tind = 0;
	myfile << xp(0) << ' ' << xp(1) << ' ' << xp(2) << ' ' << tind << ' ' << sqrt(uout(0)*uout(0)+uout(1)*uout(1)+uout(2)*uout(2)) << endl;
	for(t=0; t < 10; ++t){
			getu(xp, tind, uout);
			if(tind < 0) break;
			k1=uout;
			for(int i=0; i < 3; ++i)
				temp(i)=xp(i)+0.5*dt*k1(i);
			getu(temp, tind, uout);
			if(tind < 0) break;

			k2=uout;
			for(int i=0; i < 3; ++i)
				temp(i)=xp(i)+0.5*dt*k2(i);
			getu(temp, tind, uout);
			if(tind < 0) break;

			k3=uout;
			for(int i=0; i < 3; ++i)
				temp(i)=xp(i)+dt*k3(i);
			getu(temp, tind, uout);
			if(tind < 0) break;

			for(int i=0; i < 3; ++i)
				xp(i)=xp(i)+dt/6.0*(k1(i)+2.0*k2(i)+2.0*k3(i)+uout(i));				
			xp(0) = x0(0);

			myfile << xp(0) << ' ' << xp(1) << ' ' << xp(2) << ' ' << tind << ' ' << sqrt(uout(0)*uout(0)+uout(1)*uout(1)+uout(2)*uout(2)) << endl;
		//cout << tind+1 << ", " ;
	}
	//cout << endl;
	if(tind < 0) myfile << xp(0) << ' ' << xp(1) << ' ' << xp(2) << ' ' << tind << ' ' << sqrt(uout(0)*uout(0)+uout(1)*uout(1)+uout(2)*uout(2)) << endl,++t;
	myfile << "0 0 " << t+1 <<  " 0 0" << endl;


	
	myfile.close();
	
	return;
	
}


//void tet_hp::particle() {
//	TinyVector<FLT,4> uout,up;
//	double dudy,dwdy;
//	TinyVector<FLT,3> xp,g;
//	TinyVector<FLT,3> temp,temp2;
//	TinyVector<FLT,4> k1,k2,k3;
//	TinyVector<FLT,4> k11,k12,k21,k22,k31,k32,k41,k42;
//	int tind=30;
//	int t;
//	double tau=0.005;//relaxation time
//	double C=1.0/tau;//should be 1/tau
//	double dt=tau;// time step
//
//	ofstream myfile;
//
//	g(0)=0.0,g(1)=0.0,g(2)=-9.81;//gravity
//	
//	xp(0)=-0.2,xp(1)=0.0,xp(2)=0.2;
////	getu(xp, tind, uout);
////	up=uout;
//	up(0)=0.0,up(1)=0.0,up(2)=0.0;
//
//	
//	cout.precision(15);
//	myfile.open("x.dat");
//	myfile.precision(15);
//	
//	dt = .01;
///* Backwards Euler streamlines */
////	for(t=0; t < 1000; ++t){
////		getu(xp, tind, uout);
////		if(tind < 0) break;
////		for(int i=0; i < 3; ++i)
////			xp(i)=xp(i)+dt*uout(i);			
////		myfile << xp(0) << ' ' << xp(1) << ' ' << xp(2) << ' ' << tind << endl;
////	}
////	//if(tind < 0) myfile << xp(0) << ' ' << xp(1) << ' ' << xp(2) << ' ' << tind << endl,++t;
////	myfile << "0 0 " << t <<  " 0" << endl;
//	
//	/* Runge Kutta Streamlines */
//	dt=.01;
//   for(int part=0; part < 9; ++part){
//	   xp(0)=-.4+part*.1,xp(1)=0.2,xp(2)=0.2;
//		tind = 0;
//		myfile << xp(0) << ' ' << xp(1) << ' ' << xp(2) << ' ' << tind << endl;
//		for(t=0; t < 5000; ++t){
//			for(int time2=0; time2 < 2; ++time2){
//				getu(xp, tind, uout);
//				if(tind < 0) break;
//				k1=uout;
//				for(int i=0; i < 3; ++i)
//					temp(i)=xp(i)+0.5*dt*k1(i);
//				getu(temp, tind, uout);
//				k2=uout;
//				for(int i=0; i < 3; ++i)
//					temp(i)=xp(i)+0.5*dt*k2(i);
//				getu(temp, tind, uout);
//				k3=uout;
//				for(int i=0; i < 3; ++i)
//					temp(i)=xp(i)+dt*k3(i);
//				getu(temp, tind, uout);
//				for(int i=0; i < 3; ++i)
//					xp(i)=xp(i)+dt/6.0*(k1(i)+2.0*k2(i)+2.0*k3(i)+uout(i));				
//
//			}
//			myfile << xp(0) << ' ' << xp(1) << ' ' << xp(2) << ' ' << tind << endl;
//			//cout << tind+1 << ", " ;
//		}
//	    //cout << endl;
//		if(tind < 0) myfile << xp(0) << ' ' << xp(1) << ' ' << xp(2) << ' ' << tind << endl,++t;
//		myfile << "0 0 " << t+1 <<  " 0" << endl;
//	}
//
//	
////	/* Aerosol Kinetics including drag and lift with Runge Kutta */
////	dt=0.1*tau;//2.5*tau
////	for(int part=0; part < 9; ++part){
////		xp(0)=-.4+part*.1,xp(1)=0.2,xp(2)=0.2;
////		up(0)=0.0,up(1)=0.0,up(2)=0.0;
////		tind = 0;
////		myfile << xp(0) << ' ' << xp(1) << ' ' << xp(2) << ' ' << tind << endl;
////		for(t=0; t < 10000; ++t){
////			getudx(xp, tind, uout,dudx,dudy,dudz);
////			cout << dudx << dudy << dudz << endl;
////			if(tind < 0) break;
////			for(int i = 0; i < 3; ++i){
////				k11(i)=up(i);
////				temp(i)=xp(i)+0.5*dt*k11(i);
////			}
////			k12(0)=(C+L*(dudx(1)/sqrt(fabs(dudx(1)))+dudx(2)/sqrt(fabs(dudx(2)))))*(uout(0)-up(0))+g(0);
////			k12(1)=(C+L*(dudy(0)/sqrt(fabs(dudy(0)))+dudy(2)/sqrt(fabs(dudy(2)))))*(uout(1)-up(1))+g(1);
////			k12(2)=(C+L*(dudz(0)/sqrt(fabs(dudz(0)))+dudz(1)/sqrt(fabs(dudz(1)))))*(uout(2)-up(2))+g(2);
////			
////			for(int i = 0; i < 3; ++i){
////				temp2(i)=up(i)+0.5*dt*k12(i);
////			}
////			
////			getudx(temp, tind, uout,dudx,dudy,dudz);
////			for(int i = 0; i < 3; ++i){
////				k21(i)=temp2(i);
////				temp(i)=xp(i)+0.5*dt*k21(i);
////			}
////			k22(0)=(C+L*(dudx(1)/sqrt(fabs(dudx(1)))+dudx(2)/sqrt(fabs(dudx(2)))))*(uout(0)-temp2(0))+g(0);
////			k22(1)=(C+L*(dudy(0)/sqrt(fabs(dudy(0)))+dudy(2)/sqrt(fabs(dudy(2)))))*(uout(1)-temp2(1))+g(1);
////			k22(2)=(C+L*(dudz(0)/sqrt(fabs(dudz(0)))+dudz(1)/sqrt(fabs(dudz(1)))))*(uout(2)-temp2(2))+g(2);
////			
////			for(int i = 0; i < 3; ++i){
////				temp2(i)=up(i)+0.5*dt*k22(i);
////			}
////			
////			getudx(temp, tind, uout,dudx,dudy,dudz);
////			if(tind < 0) break;
////			for(int i = 0; i < 3; ++i){
////				k31(i)=temp2(i);
////				temp(i)=xp(i)+dt*k31(i);
////			}
////			k32(0)=(C+L*(dudx(1)/sqrt(fabs(dudx(1)))+dudx(2)/sqrt(fabs(dudx(2)))))*(uout(0)-temp2(0))+g(0);
////			k32(1)=(C+L*(dudy(0)/sqrt(fabs(dudy(0)))+dudy(2)/sqrt(fabs(dudy(2)))))*(uout(1)-temp2(1))+g(1);
////			k32(2)=(C+L*(dudz(0)/sqrt(fabs(dudz(0)))+dudz(1)/sqrt(fabs(dudz(1)))))*(uout(2)-temp2(2))+g(2);
////			
////			for(int i = 0; i < 3; ++i){
////				temp2(i)=up(i)+dt*k32(i);
////			}
////			
////			getudx(temp, tind, uout,dudx,dudy,dudz);
////			if(tind < 0) break;
////			for(int i = 0; i < 3; ++i){
////				k41(i)=temp2(i);
////			}
////			k42(0)=(C+L*(dudx(1)/sqrt(fabs(dudx(1)))+dudx(2)/sqrt(fabs(dudx(2)))))*(uout(0)-temp2(0))+g(0);
////			k42(1)=(C+L*(dudy(0)/sqrt(fabs(dudy(0)))+dudy(2)/sqrt(fabs(dudy(2)))))*(uout(1)-temp2(1))+g(1);
////			k42(2)=(C+L*(dudz(0)/sqrt(fabs(dudz(0)))+dudz(1)/sqrt(fabs(dudz(1)))))*(uout(2)-temp2(2))+g(2);
////			
////			
////			for(int i = 0; i < 3; ++i){
////				xp(i)=xp(i)+dt/6.0*(k11(i)+2.0*k21(i)+2.0*k31(i)+k41(i));
////				up(i)=up(i)+dt/6.0*(k12(i)+2.0*k22(i)+2.0*k32(i)+k42(i));
////			}
////
////			myfile << xp(0) << ' ' << xp(1) << ' ' << xp(2) << ' ' << tind << endl;
////			//cout << up(0)<<' ' << up(1) << ' ' << up(2) << endl;
////			//cout << tind+1 << ", " ;
////		}
////		//cout << endl;
////		if(tind < 0) myfile << xp(0) << ' ' << xp(1) << ' ' << xp(2) << ' ' << tind << endl,++t;
////		myfile << "0 0 " << t+1 <<  " 0" << endl;
////
////	}
//	
//		/* Aerosol Kinetics including drag with Runge Kutta */
//	dt=2.5*tau;
//	for(int part=0; part < 9; ++part){
//		xp(0)=-.4+part*.1,xp(1)=0.2,xp(2)=0.2;
//		up(0)=0.0,up(1)=0.0,up(2)=0.0;
//		tind = 0;
//		myfile << xp(0) << ' ' << xp(1) << ' ' << xp(2) << ' ' << tind << endl;
//		for(t=0; t < 20000; ++t){
//			getu(xp, tind, uout);
//			if(tind < 0) break;
//			for(int i = 0; i < 3; ++i){
//				k11(i)=up(i);
//				k12(i)=C*(uout(i)-up(i))+g(i);
//				temp(i)=xp(i)+0.5*dt*k11(i);
//				temp2(i)=up(i)+0.5*dt*k12(i);
//			}
//			getu(temp, tind, uout);
//			if(tind < 0) break;
//			for(int i = 0; i < 3; ++i){
//				k21(i)=temp2(i);
//				k22(i)=C*(uout(i)-temp2(i))+g(i);
//				temp(i)=xp(i)+0.5*dt*k21(i);
//				temp2(i)=up(i)+0.5*dt*k22(i);
//			}
//			getu(temp, tind, uout);
//			if(tind < 0) break;
//			for(int i = 0; i < 3; ++i){
//				k31(i)=temp2(i);
//				k32(i)=C*(uout(i)-temp2(i))+g(i);
//				temp(i)=xp(i)+dt*k31(i);
//				temp2(i)=up(i)+dt*k32(i);
//			}
//			getu(temp, tind, uout);
//			if(tind < 0) break;
//			for(int i = 0; i < 3; ++i){
//				k41(i)=temp2(i);
//				k42(i)=C*(uout(i)-temp2(i))+g(i);
//			}
//			for(int i = 0; i < 3; ++i){
//				xp(i)=xp(i)+dt/6.0*(k11(i)+2.0*k21(i)+2.0*k31(i)+k41(i));
//				up(i)=up(i)+dt/6.0*(k12(i)+2.0*k22(i)+2.0*k32(i)+k42(i));
//			}
//
//			myfile << xp(0) << ' ' << xp(1) << ' ' << xp(2) << ' ' << tind << endl;
//			//cout << up(0)<<' ' << up(1) << ' ' << up(2) << endl;
//			//cout << tind+1 << ", " ;
//		}
//		//cout << endl;
//		if(tind < 0) myfile << xp(0) << ' ' << xp(1) << ' ' << xp(2) << ' ' << tind << endl,++t;
//		myfile << "0 0 " << t+1 <<  " 0" << endl;
//
//	}
//	
//	
//	myfile.close();
//
//	return;
//
//}

