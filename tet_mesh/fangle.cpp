/*
 *  fangle.cpp
 *  mesh
 *
 *  Created by Michael Brazell on 6/13/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

 #include <math.h>
 #include "tet_mesh.h"

void tet_mesh::fangle(int i, double Q) {

	int N1, N2, N3;
	double RAD, PROD, D1, D2, D3, C, CHALF, T, ANGL1, ANGL2, ANGL3, S1, S2, S3 

	RAD = 45.0/atan(1.0);
	//RAD = 57.295779513082323;
	N1 = tri(i).pnt(0);
	N2 = tri(i).pnt(1);
	N3 = tri(i).pnt(2);
	PROD = (pnts(N2)(0)-pnts(N1)(0))*(pnts(N3)(0)-pnts(N1)(0))+(pnts(N2)(1)-pnts(N1)(1))*(pnts(N3)(1)-pnts(N1)(1))+(pnts(N2)(2)-pnts(N1)(2))*(pnts(N3)(2)-pnts(N1)(2));
	D1 = sqrt(pow(pnts(N2)(0)-pnts(N1)(0),2)+pow(pnts(N2)(1)-pnts(N1)(1),2)+pow(pnts(N2)(2)-pnts(N1)(2),2)); 
	D2 = sqrt(pow(pnts(N3)(0)-pnts(N1)(0),2)+pow(pnts(N3)(1)-pnts(N1)(1),2)+pow(pnts(N3)(2)-pnts(N1)(2),2)); 
	C = PROD1/D1/D2;
	CHALF = 0.5*(1.0+C);
	if (CHALF <= 1.0e-6) {
		cout << "face " << i << " has an angle of 180 degrees" << endl;
		cout << "routine fangle has stopped" << endl;
		exit(1);
	}
		
	T = sqrt(1.0/CHALF -1.0);
	ANGL1 = 2.0*atan(T);
	S1 = sin(ANGL1);
	ANGL1 = RAD*ANGL1;
	PROD = (pnts(N3)(0)-pnts(N2)(0))*(pnts(N1)(0)-pnts(N2)(0))+(pnts(N3)(1)-pnts(N2)(1))*(pnts(N1)(1)-pnts(N2)(1))+(pnts(N3)(2)-pnts(N2)(2))*(pnts(N1)(2)-pnts(N2)(2));
	D3 = sqrt(pow(pnts(N3)(0)-pnts(N2)(0),2)+pow(pnts(N3)(1)-pnts(N2)(1),2)+pow(pnts(N3)(2)-pnts(N2)(2),2)); 
	C = PROD/D3/D1;
	CHALF = 0.5*(1.0 +C);
	if (CHALF <= 1.0e-6) {
		cout << "face " << i << " has an angle of 180 degrees" << endl;
		cout << "routine fangle has stopped" << endl;
		exit(1);
	}
	T = sqrt(1.0/CHALF-1.0);
	ANGL2 = 2.0*atan(T);
	S2 = sin(ANGL2);
	ANGL2 = RAD*ANGL2;
	PROD = (pnts(N1)(0)-pnts(N3)(0))*(pnts(N2)(0)-pnts(N3)(0))+(pnts(N1)(1)-pnts(N3)(1))*(pnts(N2)(1)-pnts(N3)(1))+(pnts(N1)(2)-pnts(N3)(2))*(pnts(N2)(2)-pnts(N3)(2));
	C = PROD/D2/D3;
	CHALF = 0.5*(1.0+C);
	if (CHALF <= 1.0e-6) {
		cout << "face " << i << " has an angle of 180 degrees" << endl;
		cout << "routine fangle has stopped" << endl;
		exit(1);
	}
	
	T = sqrt(1.0/CHALF-1.0);
	ANGL3 = 2.0*atan(T)
	S3  = sin(ANGL3);
	ANGL3 = RAD*ANGL3;
	Q = 0.5*(S1  +S2  +S3)/S1/S2/S3;
			
	return;
}

