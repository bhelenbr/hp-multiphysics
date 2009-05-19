/*
 *  outputstuff.cpp
 *  tet_basis
 *
 *  Created by Michael Brazell on 3/9/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#include "tet_basis.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <utilities.h>
#include <myblas.h>

void tet_basis::outputstuff(int check) {

	/*	case 1:	 Basis function at Gauss points (plotting)
		case 2:  Basis function -1.0 -> 1.0 (plotting)
		case 3:  Basis function at 1 point 
		case 4:  Build Mass Matrix 
		case 5:  check proj.cpp using ptprobe
		case 6:  Basis function at Gauss points using proj.cpp (plotting)
		case 7:  Build Mass Matrix using proj.cpp (faster)
		case 8:  check lagrangian derivatives using proj(dx,dy,dz)
		case 9:  check lagrangian polynomial deriv with known function's derivative
		case 10: Check intgrtrst using proj and intgrt 
		case 11: Build Mass Matrix and check it using brute force method
		case 12: search for Orthogonal Mass Matrix  
		case 13: check to see if the basis can approximate pascals tetrahedral
		case 14: Build Stiffness Matrix using proj.cpp */
	
//	check = 13;

	char trans[] = "T";
	int info,n;	
	double fsave, g1,g2;
	FLT f,r,s,t,x,y,z;
	FLT N,sum;
	ofstream myfile2,fout;
	
	int stridey = gpz;
	int stridex = gpz*gpy;
	Array<int,1> ipiv(2*tm);
	Array<FLT,1> lin(tm);
	Array<FLT,1> lin2(tm);
	Array<double,1> rslt1(tm);
	Array<double,1> rslt2(tm);
	Array<double,1> rslt3(tm);
	Array<double,1> rslt4(tm);
	Array<double,1> rslt5(tm);
	Array<double,1> rslt6(tm);
	Array<double,2> wk0(tm,tm);
	Array<double,2> wk1(tm,tm);
	Array<double,2> wk2(tm,tm);
	Array<double,2> rsltm(tm,tm);
	Array<double,2> rsltm2(tm,tm);
	Array<double,3> diff(gpx,gpy,gpz);
	Array<double,3> f1(gpx,gpy,gpz);
	Array<double,3> f2(gpx,gpy,gpz);
	Array<double,3> dx(gpx,gpy,gpz);
	Array<double,3> dx1(gpx,gpy,gpz);
	Array<double,3> dy(gpx,gpy,gpz);
	Array<double,3> dy1(gpx,gpy,gpz);
	Array<double,3> dz(gpx,gpy,gpz);
	Array<double,3> dz1(gpx,gpy,gpz);
	Array<double,3> zeros(gpx,gpy,gpz);
	
	Array<double,2> f12d(gpx,gpy);


	switch (check){
	
	case 1:	
	
		/* Basis function at Gauss points  */
		for (n = 0; n < tm; ++n) {
			std::cout << "VARIABLES = \"x\" \"y\" \"z\" \"f\"" << std::endl;
			std::cout << "ZONE K=" << gpx << ", J=" <<gpx << ", I=" << gpx<< ", F=POINT" << std::endl;

			lin = 0.0;
			lin(n) = 1.0;		
		
			for(int i = 0; i < gpx; ++i) {
				for(int j = 0; j < gpy; ++j) {
					for(int k = 0; k < gpz; ++k) {
				  
						t = zp(k);
						s = 0.5*(yp(j)+1.0)*(1.0-t)-1.0;
						r = 0.5*(xp(i)+1.0)*(-s -t)-1.0;				
					
						ptprobe(1, &f,r,s,t,lin.data(),tm);
					
						std::cout << r << ' ' << s << ' ' << t << ' ' << f << std::endl;
					}
				}
			}
			std::cout << std::endl << std::endl;		
		}
		break;
	
	
//----------------------------------------------------------------------
	
	case 2:
	
	
		myfile2.open("basis.dat");
		myfile2.precision(15);
		
		
		/* Basis function -1.0 -> 1.0  */
		N = 9;
		n = 4;
		//for (n = 0; n < tm; ++n) {
			myfile2 << "VARIABLES = \"x\" \"y\" \"z\" \"f\"" << std::endl;
			myfile2 << "ZONE K=" << N+1 << ", J=" <<N+1 << ", I=" << N+1<< ", F=POINT" << std::endl;

			lin = 0.0;
			lin(n) = 1.0;
			x=-1.0;	
			for(int i = 0; i <= N; ++i) {
				y=-1.0;
				for(int j = 0; j <= N; ++j) {
					z=-1.0;
					for(int k = 0; k <= N; ++k) {
				  
						t = z;
						s = 0.5*(y+1.0)*(1.0-t)-1.0;
						r = 0.5*(x+1.0)*(-s-t)-1.0;				
					
						ptprobe(1, &f,r,s,t,lin.data(),tm);
					
						myfile2 << r << ' ' << s << ' ' << t << ' ' << f << std::endl;
						z += 2/N;
					}
					y += 2/N;
				}
				x += 2/N;
			}
			myfile2 << std::endl << std::endl;
		//}
		
		myfile2.close();
		break;

//----------------------------------------------------------------------
	
	case 3:
		/* Basis function at 1 point */	
		for (n = 0; n < tm; ++n) {
			lin = 0.0;
			lin(n) = 1.0;
			
			x=-.2;
			y=-.2;
			z=-.2;
			
			t = z;
			s = 0.5*(y+1.0)*(1.0-t)-1.0;
			r = 0.5*(x+1.0)*(-s -t)-1.0;				
		
			ptprobe(1, &f,r,s,t,lin.data(),tm);
					
			std::cout << n+1 << ' ' << f << std::endl;
			std::cout << std::endl ;
		}
		break;

//---------------------------------------------------------------------------------	

	case 4:
		/* Build Mass Matrix */
		for (n = 0; n < tm; ++n) {
			lin = 0.0;
			lin(n) = 1.0;		
			//std::cout << "VARIABLES = \"x\" \"y\" \"z\" \"f\"" << std::endl;
			//std::cout << "ZONE K=" << gpx << ", J=" << gpx << ", I=" << gpx << ", F=POINT" << std::endl;
			for(int i = 0; i < gpx; ++i) {
				for(int j = 0; j < gpy; ++j) {
					for(int k = 0; k < gpz; ++k){
				  
						t = zp(k);
						s = 0.5*(yp(j)+1.0)*(1.0-t)-1.0;
						r = 0.5*(xp(i)+1.0)*(-s-t)-1.0;
										
						ptprobe(1, &f1(i,j,k),r,s,t,lin.data(),tm);
						//f1(i,j,k)=1.0;				
						//std::cout << r << ' ' << s << ' ' << t << ' ' << f1(i,j,k) << std::endl;

					}
				}
			}
			
			//std::cout << n << std::endl;
			intgrt(rslt1.data(), f1.data(), stridex, stridey);
			
			for(int m = 0; m < tm; ++m){
				rsltm(m,n) = rslt1(m);
			}
				
		}
		
		//std::cout << rsltm<<endl;
		
		myfile2.open("mass.dat");
		myfile2.precision(15);
		
		for(int i = 0; i <tm; ++i){
			for(int j = 0; j <tm; ++j){

				myfile2 << rsltm(i,j) << ' ' ;	
				
			}
			myfile2 << endl;
			
		}
		
		myfile2.close();
		break;

	
//---------------------------------------------------------------------------------	


	case 5:
		/* check proj.cpp using ptprobe */
			
		for (n = 0; n < tm; ++n) {
			fsave = 0.0;
			lin = 0.0;
			lin(n) = 1.0;		
			proj(lin.data(), &f1(0,0,0), stridex, stridey);

			for(int i = 0; i < gpx; ++i) {
				for(int j = 0; j < gpy; ++j) {
					for(int k = 0; k < gpz; ++k) {
				  
						t = zp(k);
						s = 0.5*(yp(j)+1.0)*(1.0-t)-1.0;
						r = 0.5*(xp(i)+1.0)*(-s-t)-1.0;				
					
						ptprobe(1, &f,r,s,t,lin.data(),tm);
					
	//					std::cout << f-f1(i,j,k) << std::endl;
						if (fabs(f-f1(i,j,k)) > fsave){
							fsave = fabs(f-f1(i,j,k));
						}
					}
				}
			}
	 
			std::cout << fsave << endl;
	//		std::cout << std::endl << std::endl;		
		}
		break;
	
	
	
//---------------------------------------------------------------------------------	

	case 6:
		/* Basis function at Gauss points using proj.cpp  */
		
		myfile2.open("basisproj.dat");
		myfile2.precision(15);
		n=15;
		//for (n = 0; n < tm; ++n) {
			myfile2 << "VARIABLES = \"x\" \"y\" \"z\" \"f\"" << std::endl;
			myfile2 << "ZONE K=" << gpx << ", J=" <<gpx << ", I=" << gpx<< ", F=POINT" << std::endl;

			lin = 0.0;
			lin(n) = 1.0;		
			proj(lin.data(), &f1(0,0,0), stridex, stridey);
		
			for(int i = 0; i < gpx; ++i) {
				for(int j = 0; j < gpy; ++j) {
					for(int k = 0; k < gpz; ++k) {
				  
						t = zp(k);
						s = 0.5*(yp(j)+1.0)*(1.0-t)-1.0;
						r = 0.5*(xp(i)+1.0)*(-s -t)-1.0;		
					
						myfile2 << r << ' ' << s << ' ' << t << ' ' << f1(i,j,k) << std::endl;
					}
				}
			}
			myfile2 << std::endl << std::endl;		
		//}
		break;
		
		myfile2.close();



//---------------------------------------------------------------------------------	

	case 7:

		/* Build Mass Matrix using proj.cpp */

		for (n = 0; n < tm; ++n) {
			lin = 0.0;
			lin(n) = 1.0;	
			proj(lin.data(), &f1(0,0,0), stridex, stridey);

			intgrt(rslt1.data(), f1.data(), stridex, stridey);
			
			for(int m = 0; m < tm; ++m){
				rsltm(m,n) = rslt1(m);
			}
				
		}
		
		
		
	//	ofstream myfile2;
		myfile2.open("p16mass.dat");
		myfile2.precision(15);
		
		for(int i = 0; i <tm; ++i){
			for(int j = 0; j <tm; ++j){

				myfile2 << rsltm(i,j) << ' ' ;	
				
			}
			myfile2 << endl;
			
		}
		
		myfile2.close();
		break;


//---------------------------------------------------------------------------------	
	
	case 8:
		/* check lagrangian derivatives using proj(dx,dy,dz) */

	//	lin = 0.0;
	//	lin(2) = 1.0;
		lin=1.0;
		dx = 0.0;
		dy = 0.0;	
		dz = 0.0;
		dx1 = 0.0;
		dy1 = 0.0;
		dz1 = 0.0;
		proj(lin.data(), &f1(0,0,0), &dx(0,0,0), &dy(0,0,0), &dz(0,0,0),  stridex,  stridey);
	//	std::cout << f1 << std::endl;
	//	std::cout << dx << std::endl;
		derivr(f1.data(), dx1.data(), stridex, stridey);	
		derivs(f1.data(), dy1.data(), stridex, stridey);		
		derivt(f1.data(), dz1.data(), stridex, stridey);
		
		diff = dx-dx1+dy-dy1+dz-dz1;
		
		std::cout << max(fabs(diff)) << std::endl;

		break;
		
//---------------------------------------------------------------------------------	

	case 9:
		/*  check lagrangian polynomial deriv with known function's derivative  */

		for(int i = 0; i < gpx; ++i) {
			for(int j = 0; j < gpy; ++j) {
				for(int k = 0; k < gpz; ++k) {
					t = zp(k);
					s = 0.5*(yp(j)+1.0)*(1.0-t)-1.0;
					r = 0.5*(xp(i)+1.0)*(-s-t)-1.0;
					f1(i,j,k) = pow(r,p); 
					f2(i,j,k) = p*pow(r,p-1);
				}
			}
		}
		
		dx = 0;	
		dy = 0;
		dz = 0;
		derivr(f1.data(), dx.data(), stridex, stridey);
		derivs(f1.data(), dy.data(), stridex, stridey);
		derivt(f1.data(), dz.data(), stridex, stridey);

	   
	//   for(int i = 0; i < gpx; ++i) {
	//		for(int j = 0; j < gpy; ++j) {
	//			for(int k = 0; k < gpz; ++k) {
	//				t = zp(k);
	//				s = 0.5*(yp(j)+1.0)*(1.0-t)-1.0;
	//				r = 0.5*(xp(i)+1.0)*(-s-t)-1.0;
	//            std::cout << dx(i,j,k)*(1+xp(i))/2. - dy(i,j,k) << '\n';
	//			}
	//		}
	//	}
	   
		diff= dx-f2;
//		diff= dy-f2;
//		diff= dz-f2;

	   
		std::cout << max(fabs(diff)) << std::endl;
		break;


//---------------------------------------------------------------------------------	

	case 10:
		/* Check intgrtrst using proj and intgrt   */
		
		dx = 0.0;
		dy = 0.0;	
		dz = 0.0;
		zeros = 0.0;
		
		for(int i = 0; i < tm; ++i){
			lin = 0;
			lin(i) = 1;
			rslt4 = 0.0;
			rslt5 = 0.0;
			rslt6 = 0.0;

			proj(lin.data(), &f1(0,0,0), &dx(0,0,0), &dy(0,0,0), &dz(0,0,0),  stridex,  stridey);
			intgrt(rslt1.data(),&dx(0,0,0), stridex, stridey);
			intgrt(rslt2.data(),&dy(0,0,0), stridex, stridey);
			intgrt(rslt3.data(),&dz(0,0,0), stridex, stridey);

			intgrtrst(rslt4.data(),&f1(0,0,0),&zeros(0,0,0),&zeros(0,0,0), stridex, stridey);
			intgrtrst(rslt5.data(),&zeros(0,0,0),&f1(0,0,0),&zeros(0,0,0), stridex, stridey);
			intgrtrst(rslt6.data(),&zeros(0,0,0),&zeros(0,0,0),&f1(0,0,0), stridex, stridey);

			
			for(int j = 0; j < tm; ++j){
				wk0(i,j) = rslt1(j)+rslt2(j)+rslt3(j);
				wk1(j,i) = rslt4(j)+rslt5(j)+rslt6(j);		
			}
		}
		
		for(int i = 0; i < tm; ++i){
			for(int j = 0; j < tm; ++j){
				sum=wk0(i,j)+wk1(i,j);
				if(fabs(sum) >= 10e-15)
					std::cout << i << ' ' << j << ' ' << sum << '\n';
			}
		}
		
		sum=max(wk0+wk1);
		std::cout << sum << endl;
		break;	


//---------------------------------------------------------------------------------	


	case 11:

		/* Build Mass Matrix and check it using brute force method */


		for (n = 0; n < tm; ++n) {
			lin = 0.0;
			lin(n) = 1.0;		
			//std::cout << "VARIABLES = \"x\" \"y\" \"z\" \"f\"" << std::endl;
			//std::cout << "ZONE K=" << gpx << ", J=" << gpx << ", I=" << gpx << ", F=POINT" << std::endl;
			for(int i = 0; i < gpx; ++i) {
				for(int j = 0; j < gpy; ++j) {
					for(int k = 0; k < gpz; ++k){
				  
						t = zp(k);
						s = 0.5*(yp(j)+1.0)*(1.0-t)-1.0;
						r = 0.5*(xp(i)+1.0)*(-s-t)-1.0;
										
						ptprobe(1, &f1(i,j,k),r,s,t,lin.data(),tm);
						//f1(i,j,k)=1.0;				
						//std::cout << r << ' ' << s << ' ' << t << ' ' << f1(i,j,k) << std::endl;

					}
				}
			}
			
			//std::cout << n << std::endl;
			intgrt(rslt1.data(), f1.data(), stridex, stridey);
			
			for(int m = 0; m < tm; ++m){
				rsltm(m,n) = rslt1(m);
			}
				
		}
		
		//std::cout << rsltm<<endl;
		
	   rsltm2 = 0.0;
	   
	   for (int m = 0; m < tm; ++m) {
		  for (n = 0; n < tm; ++n) {
			 lin = 0.0;
			 lin(n) = 1.0;
			 lin2 = 0.0;
			 lin2(m) = 1.0;
			 rsltm2(m,n) = 0.0;
			 
			 for(int i = 0; i < gpx; ++i) {
				for(int j = 0; j < gpy; ++j) {
				   for(int k = 0; k < gpz; ++k){
				  
					  t = zp(k);
					  s = 0.5*(yp(j)+1.0)*(1.0-t)-1.0;
					  r = 0.5*(xp(i)+1.0)*(-s-t)-1.0;
								  
					  ptprobe(1, &g1,r,s,t,lin.data(),tm);
					  ptprobe(1, &g2,r,s,t,lin2.data(),tm);
					  rsltm2(m,n) += wtx(i)*wty(j)*wtz(k)*g1*g2;
				   }
				}
			 }
		  }
		}
				
		for(int i = 0; i <tm; ++i){
			for(int j = 0; j <tm; ++j){
				if(fabs(rsltm2(i,j) -rsltm(i,j)) > 10e-14)
					std::cout << i << ' ' << j << ' ' << rsltm(i,j) << ' ' << rsltm2(i,j) -rsltm(i,j) << '\n' ;	
				
			}
		
			
		}
		break;
	


//---------------------------------------------------------------------------------	

	case 12:

		/* Orthogonal Mass Matrix  */

		lin = 0.0;

		//vertex 0 p4
		lin=0;
		lin(0)=1.0;
		lin(13)=-2.314285714285524;
		lin(14)=-2.399999999999402;
		lin(15)=-1.500000000001130;
		lin(16)=-2.314285714285482;
		lin(17)=-2.399999999999338;
		lin(18)=-1.500000000001286;
		lin(19)=-2.314285714285487;
		lin(20)=-2.399999999999480;
		lin(21)=-1.500000000000970;
		
//		//edge 1 p4
//		lin=0;
//		lin(4)=1.0;
//		lin(5)=0.0;
//		lin(6)=2.1;
//		lin(22)=-4.0;
//		lin(23)=1.2;
//		lin(24)=0.0;
//		lin(25)=-4.0;
//		lin(26)=1.2;
//		lin(27)=0.0;
//		lin(34)=16.6;
//		
//		// new face mode
//		lin = 0;
//		lin(22) = 1;
//		lin(34) = -11;
		
//		// face 0 vertex ball
//		lin = 0;
//		lin(4)=-0.02681;//edge 1
//		lin(5)=-0.00000;
//		lin(6)=1.10135;
//		lin(7)=-0.08130;//edge 2
//		lin(8)=0.16345;
//		lin(9)=-0.12452;
//		lin(10)=-0.08130;//edge 3
//		lin(11)=0.16345;
//		lin(12)=-0.12452;
//		lin(22)=1.0000;//face 0
//		lin(23)=-0.97008;
//		lin(24)=0.0;
//		lin(34)=2.16795;//interior
	
		proj(lin.data(), &f1(0,0,0), stridex, stridey);
		intgrt(rslt1.data(), f1.data(), stridex, stridey);

		
		std::cout << rslt1 << std::endl;
		for(int i = 0; i<tm; ++i){
			cout << i+1 << ' ' << rslt1(i) << endl;
		}
		
		N = 20;
		
		fout.open("basis.dat");
		fout.precision(15);		
		
		
		fout << "VARIABLES = \"x\" \"y\" \"z\" \"f\"" << std::endl;
		fout << "ZONE K=" << N+1 << ", J=" <<N+1 << ", I=" << N+1<< ", F=POINT" << std::endl;

		x=-1.0;	
		for(int i = 0; i <= N; ++i) {
			y=-1.0;
			for(int j = 0; j <= N; ++j) {
				z=-1.0;
				for(int k = 0; k <= N; ++k) {
			  
					t = z;
					s = 0.5*(y+1.0)*(1.0-t)-1.0;
					r = 0.5*(x+1.0)*(-s -t)-1.0;				
				
					ptprobe(1, &f,r,s,t,lin.data(),tm);
				
					fout << r << ' ' << s << ' ' << t << ' ' << f << std::endl;
					z += 2/N;
				}
				y += 2/N;
			}
			x += 2/N;
		}
		fout << std::endl << std::endl;
		
	//	// 1d plot of basis
	//	N = 100;
	//
	//	x=-1.0;	
	//	y=-1.0;
	//	z=-1.0;
	//
	//	for(int k = 0; k <= N; ++k) {
	//
	//		t = z;
	//		s = 0.5*(y+1.0)*(1.0-t)-1.0;
	//		r = 0.5*(x+1.0)*(-s-t)-1.0;				
	//
	//		ptprobe(1, &f,r,s,t,lin.data(),tm);
	//
	//		std::cout << f << std::endl;
	//		
	//		z += 2/N;
	//	}
	//
	//	std::cout << std::endl << std::endl;
	
	
		fout.close();
		
		//cout << lgrnge3d(0,1,1,1) << endl;
		break;



//---------------------------------------------------------------------------------	

	case 13:
	 // check to see if the basis can approximate pascals tetrahedral

		for (n = 0; n < tm; ++n) {
			lin = 0.0;
			lin(n) = 1.0;	
			proj(lin.data(), &f1(0,0,0), stridex, stridey);

			intgrt(rslt1.data(), f1.data(), stridex, stridey);
			
			for(int m = 0; m < tm; ++m){
				rsltm(m,n) = rslt1(m);
			}
				
		}
		
		for(int i = 0; i < gpx; ++i) {
			for(int j = 0; j < gpy; ++j) {
				for(int k = 0; k < gpz; ++k){
			  
					t = zp(k);
					s = 0.5*(yp(j)+1.0)*(1.0-t)-1.0;
					r = 0.5*(xp(i)+1.0)*(-s-t)-1.0;
									
					f1(i,j,k)=r*s*t*r;

				}
			}
		}
		
		intgrt(rslt1.data(), f1.data(), stridex, stridey);
		
	 //	cout << rslt1 << endl;
		
		GETRF(tm,tm,&rsltm(0,0),tm,&ipiv(0),info);

		GETRS(trans,tm,1,&rsltm(0,0),tm,&ipiv(0),&rslt1(0),tm,info);
		
	 //	cout << rslt1 << endl;
		
		proj(rslt1.data(), &f2(0,0,0), stridex, stridey);
		
	//	for(int i = 0; i < gpx; ++i) {
	//		for(int j = 0; j < gpy; ++j) {
	//			for(int k = 0; k < gpz; ++k){
	//		  
	//				cout << f1(i,j,k)-f2(i,j,k) << endl;
	//
	//			}
	//		}
	//	}
		
		cout << max(fabs(f1-f2)) << endl;
		break;
	
//---------------------------------------------------------------------------------	
	
	case 14:
	
		/* Build Stiffness Matrix using proj.cpp */

		for (n = 0; n < tm; ++n) {
			lin = 0.0;
			lin(n) = 1.0;	
			rslt1=0.0;

			proj(lin.data(), &f1(0,0,0), &dx(0,0,0),&dy(0,0,0),&dz(0,0,0),stridex, stridey);
			intgrtrst(rslt1.data(), dx.data(),dy.data(),dz.data(), stridex, stridey);

			for(int m = 0; m < tm; ++m){
				rsltm(m,n) = rslt1(m);
			}
				
		}
		
			
		myfile2.open("stiff.dat");
		myfile2.precision(15);
		
		for(int i = 0; i <tm; ++i){
			for(int j = 0; j <tm; ++j){

				myfile2 << rsltm(i,j) << ' ' ;	
				
			}
			myfile2 << endl;
			
		}
		
		myfile2.close();
		break;	
	
	

//---------------------------------------------------------------------------------	


	case 15:
		/* check proj2d.cpp using ptprobe2d */
			
		for (n = 0; n < 3+3*em+fm; ++n) {
			fsave = 0.0;
			lin = 0.0;
			lin(n) = 1.0;		
			proj2d(lin.data(), &f12d(0,0), stridey);

			for(int i = 0; i < gpx; ++i) {
				for(int j = 0; j < gpy; ++j) {				  
					
						s = yp(j);
						r = 0.5*(xp(i)+1.0)*(1-s)-1.0;									
					
						ptprobe2d(1, &f,r,s,lin.data(),3+3*em+fm);

	//					std::cout << f-f1(i,j,k) << std::endl;
						if (fabs(f-f1(i,j)) > fsave){
							fsave = fabs(f-f12d(i,j));
						}	
				}
			}
	 
			std::cout << fsave << endl;
	//		std::cout << std::endl << std::endl;		
		}
		break;
		
		
		//---------------------------------------------------------------------------------	


	case 16:

		/* Build Mass Matrix2d and check it using brute force method */


		for (n = 0; n < 3+3*em+fm; ++n) {
			lin = 0.0;
			lin(n) = 1.0;		
			//std::cout << "VARIABLES = \"x\" \"y\" \"z\" \"f\"" << std::endl;
			//std::cout << "ZONE K=" << gpx << ", J=" << gpx << ", I=" << gpx << ", F=POINT" << std::endl;
			for(int i = 0; i < gpx; ++i) {
				for(int j = 0; j < gpy; ++j) {				  
						s = yp(j);
						r = 0.5*(xp(i)+1.0)*(1-s)-1.0;
										
						ptprobe2d(1, &f12d(i,j),r,s,lin.data(),3+3*em+fm);
						//f1(i,j,k)=1.0;				
						//std::cout << r << ' ' << s << ' ' << t << ' ' << f1(i,j,k) << std::endl;
				}
			}
			
			//std::cout << n << std::endl;
			intgrt2d(rslt1.data(), f12d.data(), stridey);
			
			for(int m = 0; m < 3+3*em+fm; ++m){
				rsltm(m,n) = rslt1(m);
			}
				
		}
		
		//std::cout << rsltm<<endl;
		
	   rsltm2 = 0.0;
	   
	   for (int m = 0; m < 3+3*em+fm; ++m) {
		  for (n = 0; n < 3+3*em+fm; ++n) {
			 lin = 0.0;
			 lin(n) = 1.0;
			 lin2 = 0.0;
			 lin2(m) = 1.0;
			 rsltm2(m,n) = 0.0;
			 
			 for(int i = 0; i < gpx; ++i) {
				for(int j = 0; j < gpy; ++j) {
				  
						s = yp(j);
						r = 0.5*(xp(i)+1.0)*(1-s)-1.0;
								  
					  ptprobe2d(1, &g1,r,s,lin.data(),3+3*em+fm);
					  ptprobe2d(1, &g2,r,s,lin2.data(),3+3*em+fm);
					  rsltm2(m,n) += wtx(i)*wty(j)*g1*g2;

				}
			 }
		  }
		}
				
		for(int i = 0; i <3+3*em+fm; ++i){
			for(int j = 0; j <3+3*em+fm; ++j){
				if (fabs(rsltm2(i,j) -rsltm(i,j)) > 10e-15){
						std::cout << i << ' ' << j << ' ' << rsltm(i,j) << ' ' << rsltm2(i,j) -rsltm(i,j) << '\n' ;
				}	
				
			}
		
			
		}
		break;


		//---------------------------------------------------------------------------------	


	case 17:
		/* Check basis one point at a time p = 4 */
		r = 0.1;
		s = 0.1;
		t = 0.1;
		
		x = 2.0*(1+r)/(-s-t) -1.0;
		y = 2.0*(1+s)/(1-t)-1.0;
		z = t;
				
		lin = 0;

//		// vertex 0 works
//		lin(0) = 1.0;
//		g2 = (1+t)/2;

//		// vertex 1 works
//		lin(1) = 1.0;
//		g2 = (1+s)/2;

//		// vertex 2 works
//		lin(2) = 1.0;
//		g2 = (-r-s-t-1)/2;

//		// vertex 3 works
//		lin(3) = 1.0;
//		g2 = (1+r)/2;

//		// edge 1 mode 1 works
//		lin(4) = 1.0;
//		g2 = (-r-s-t-1)/2*(1+r)/2;

//		// edge 2 mode 1 works
//		lin(7) = 1.0;
//		g2 = (1+s)/2*(1+r)/2;

//		// edge 3 mode 1 works
//		lin(10) = 1.0;
//		g2 = (-r-s-t-1)/2*(1+s)/2;

//		// edge 4 mode 1 works
//		lin(13) = 1.0;
//		g2 = (1+t)/2*(1+s)/2;
		
//		// edge 5 mode 1 works
//		lin(16) = 1.0;
//		g2 = (-r-s-t-1)/2*(1+t)/2;

//		// edge 6 mode 1 works
//		lin(19) = 1.0;
//		g2 = (1+r)/2*(1+t)/2;

//		// face 0 mode 1 works
//		lin(22) = 1.0;
//		g2 = (-r-s-t-1)/2*(1+r)/2*(1+s)/2;

//		// face 1 mode 1 works
//		lin(25) = 1.0;
//		g2 = (1-x)*(1+x)/4*(1-y)*(1-y)/4*(1-z)*(1-z)/4*(1+z)/2;

//		// face 2 mode 1 works
//		lin(28) = 1.0;
//		g2 = (1+x)/2*(1-y)*(1+y)/4*(1-z)*(1-z)/4*(1+z)/2;

		// face 3 mode 1 works
		lin(31) = 1.0;
		g2 = (1-x)/2*(1-y)*(1+y)/4*(1-z)*(1-z)/4*(1+z)/2;

		ptprobe(1, &g1,r,s,t,lin.data(),tm) ;
		
		cout << "function = " << g2 << " probe = " << g1 << endl;
		//cout << x << ' ' << y << ' ' << z << endl;
	
		


		break;

	
	
}


}
