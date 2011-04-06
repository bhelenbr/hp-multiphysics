/*
 *  tet_nodal_basis.cpp
 *  tet_basis
 *
 *  Created by michael brazell on 4/6/11.
 *  Copyright 2011 Clarkson University. All rights reserved.
 *
 */

#include "tet_nodal_basis.h"

//Array<tet_basis,1> basis::tet;

void tet_nodal_basis::initialize(int pdegree, int gpoints) {
	
	tm = (pdegree+1)*(pdegree+2)*(pdegree+3)/6;
	
	ipiv.resize(tm);
	int info;
	P.resize(tm,tm);
	P = 0;
	for(int i = 0; i < tm; ++i)
		P(i,i) = 1.0;
	
	GETRF(tm,tm,&P(0,0),tm,&ipiv(0),info);
	
	tet_basis::initialize(pdegree, gpoints);	
	
	FLT lcl;
	Array<FLT,1> lin(tm);
	Array<FLT,2> M(tm,tm);
	Array<FLT,1> u(tm),l(tm);   
	int stridey = gpx;
	int stridex = gpx*gpx;
	Array<FLT,3> wk(gpx,gpx,gpx);
	Array<FLT,1> nodes(p-1);
	int gl = 1;
	
	if(p == 2){
		nodes(0)=0.0;
	}
	if(p == 3){
		if(gl == 1){
			nodes(0)=-0.447213595499958;
			nodes(1)=0.447213595499958;
		}
		else{
			nodes(0)=-1.0/3.0;
			nodes(1)=1.0/3.0;
		}
	}
	/* VERTEX 0 */
	for(int m = 0; m < tm; ++m){
		lin = 0.0;
		lin(m) = 1.0;
		ptprobe(1, &P(0,m),-1.0,-1.0,1.0,lin.data(),tm);			
	}
	
	/* VERTEX 1 */
	for(int m = 0; m < tm; ++m){
		lin = 0.0;
		lin(m) = 1.0;
		ptprobe(1, &P(1,m),-1.0,1.0,-1.0,lin.data(),tm);			
	}
	
	/* VERTEX 2 */
	for(int m = 0; m < tm; ++m){
		lin = 0.0;
		lin(m) = 1.0;
		ptprobe(1, &P(2,m),-1.0,-1.0,-1.0,lin.data(),tm);			
	}
	
	/* VERTEX 3 */
	for(int m = 0; m < tm; ++m){
		lin = 0.0;
		lin(m) = 1.0;
		ptprobe(1, &P(3,m),1.0,-1.0,-1.0,lin.data(),tm);			
	}
	
	if(p > 1){
		/* EDGE 1 */
		for(int i = 0; i < em; ++i){
			for(int m = 0; m < tm; ++m){
				lin = 0.0;
				lin(m) = 1.0;
				ptprobe(1, &P(vm+i,m),nodes(i),-1.0,-1.0,lin.data(),tm);			
			}
		}
		
		/* EDGE 2 */
		for(int i = 0; i < em; ++i){
			for(int m = 0; m < tm; ++m){
				lin = 0.0;
				lin(m) = 1.0;
				ptprobe(1, &P(vm+em+i,m),nodes(em-i-1),nodes(i),-1.0,lin.data(),tm);			
			}
		}
		
		/* EDGE 3 */
		for(int i = 0; i < em; ++i){
			for(int m = 0; m < tm; ++m){
				lin = 0.0;
				lin(m) = 1.0;
				ptprobe(1, &P(vm+2*em+i,m),-1.0,nodes(i),-1.0,lin.data(),tm);			
			}
		}
		
		/* EDGE 4 */
		for(int i = 0; i < em; ++i){
			for(int m = 0; m < tm; ++m){
				lin = 0.0;
				lin(m) = 1.0;
				ptprobe(1, &P(vm+4*em-i-1,m),-1.0,nodes(i),nodes(em-i-1),lin.data(),tm);			
			}
		}
		
		/* EDGE 5 */
		for(int i = 0; i < em; ++i){
			for(int m = 0; m < tm; ++m){
				lin = 0.0;
				lin(m) = 1.0;
				ptprobe(1, &P(vm+4*em+i,m),-1.0,-1.0,nodes(i),lin.data(),tm);			
			}
		}
		
		/* EDGE 6 */
		for(int i = 0; i < em; ++i){
			for(int m = 0; m < tm; ++m){
				lin = 0.0;
				lin(m) = 1.0;
				ptprobe(1, &P(vm+6*em-i-1,m),nodes(i),-1.0,nodes(em-i-1),lin.data(),tm);			
			}
		}
		if(p == 3){
			/* harcoded in needs more thought for higher order */
			/* FACE 0*/				
			for(int m = 0; m < tm; ++m){
				lin = 0.0;
				lin(m) = 1.0;
				ptprobe(1, &P(vm+6*em,m),-1.0/3.0,-1.0/3.0,-1.0,lin.data(),tm);			
			}		
			/* FACE 1 */				
			for(int m = 0; m < tm; ++m){
				lin = 0.0;
				lin(m) = 1.0;
				ptprobe(1, &P(vm+6*em+fm,m),-1.0/3.0,-1.0,-1.0/3.0,lin.data(),tm);			
			}
			/* FACE 2 */				
			for(int m = 0; m < tm; ++m){
				lin = 0.0;
				lin(m) = 1.0;
				ptprobe(1, &P(vm+6*em+2*fm,m),-1.0/3.0,-1.0/3.0,-1.0/3.0,lin.data(),tm);			
			}
			/* FACE 3 */				
			for(int m = 0; m < tm; ++m){
				lin = 0.0;
				lin(m) = 1.0;
				ptprobe(1, &P(vm+6*em+3*fm,m),-1.0,-1.0/3.0,-1.0/3.0,lin.data(),tm);			
			}		
		}			
	}
	if(p > 3){
		cout << "p too high for nodal basis due to orientation issues" << endl;
		//exit(3); 
	}
	
	GETRF(tm,tm,&P(0,0),tm,&ipiv(0),info);
	
	/* MAKE MASS MATRIX */
	for(int m=0;m<tm;++m) {
		u=0;
		u(m) = 1.0;
		proj(&u(0),&wk(0,0,0),stridex,stridey); 
		intgrt(l.data(),wk.data(),stridex,stridey); 
		
		for(int i=0;i<tm;++i) 
			M(m,i) = l(i);      
	}	
	
	/* LUMPED MASS MATRIX*/
	lcl=0.0;
	for(int i = 0; i < tm; ++i)
		lcl+=fabs(M(0,i));
	vdiag = 1.0/lcl;
	
	if(p > 1){
		lcl = 0.0;
		for(int i = 0; i < tm; ++i)
			lcl+=fabs(M(4,i));
		ediag(0)=1.0/lcl;
	}
	if(p > 2){
		lcl = 0.0;
		for(int i = 0; i < tm; ++i)
			lcl+=fabs(M(5,i));
		ediag(1)=1.0/lcl;		
		
		lcl = 0.0;
		for(int i = 0; i < tm; ++i)
			lcl+=fabs(M(16,i));
		fdiag(0)=1.0/lcl;
	}
	
	/* DIAGONAL MASS MATRIX */	
	
	vdiag = 1.0/M(0,0);
	
	if(p > 1){
		ediag(0)=1.0/M(4,4);
	}
	
	if(p > 2){
		ediag(1)=1.0/M(5,5);
		fdiag(0)=1.0/M(16,16);
	}
	
	//		cout << vdiag << endl << ediag <<endl<< fdiag << endl;
	
	
}

