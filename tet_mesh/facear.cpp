/*
 *  fangle.cpp
 *  mesh
 *
 *  Created by Michael Brazell on 6/13/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

 #include <math.h>
 #include "mesh.h"

void tet_mesh::facear(int N1, int N2, int N3) {
	double AX, AY, AZ, AREA
	
	AX = (pnts(N2)(1)-pnts(N1)(1))*(pnts(N3)(2)-pnts(N1)(2))-(pnts(N2)(2)-pnts(N1)(2))*(pnts(N3)(1)-pnts(N1)(1));

	AY = (pnts(N2)(2)-pnts(N1)(2))*(pnts(N3)(0)-pnts(N1)(0))-(pnts(N2)(0)-pnts(N1)(0))*(pnts(N3)(2)-pnts(N1)(2));

	AZ = (pnts(N2)(0)-pnts(N1)(0))*(pnts(N3)(1)-pnts(N1)(1))-(pnts(N2)(1)-pnts(N1)(1))*(pnts(N3)(0)-pnts(N1)(0));
	
	return(AREA);
}

