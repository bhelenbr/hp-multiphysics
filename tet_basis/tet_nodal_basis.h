/*
 *  tet_nodal_basis.h
 *  tet_basis
 *
 *  Created by michael brazell on 4/6/11.
 *  Copyright 2011 Clarkson University. All rights reserved.
 *
 */

#ifndef _tet_nodal_basis_h_
#define _tet_nodal_basis_h_

#include "tet_basis.h"
#include <myblas.h>

class tet_nodal_basis : public tet_basis {
public:
	Array<FLT,2> P;
	Array<int,1> ipiv;
	
	void initialize(int pdegree, int gpoints);
	inline void initialize(int pdegree) { initialize(pdegree, pdegree+1);}
	
	
	inline void proj(FLT *lin1, FLT *f1, FLT *dx, FLT *dy, FLT *dz, int stridex, int stridey) {tet_basis::proj(lin1,f1,dx,dy,dz,stridex,stridey);}

	inline void proj(FLT u1, FLT u2, FLT u3, FLT u4, FLT *f, int stridex, int stridey) {tet_basis::proj(u1,u2,u3,u4,f,stridex,stridey);}
	
	
	void proj(FLT *lin, FLT *f1, int stridex, int stridey) {
		int info;
		char trans[] = "T";

		GETRS(trans,tm,1,P.data(),tm,&ipiv(0),lin,tm,info);
		tet_basis::proj(lin,f1,stridex,stridey);
	}
	
	void intgrt(FLT *lin, FLT *f1, int stridex, int stridey){
		int info;
		char normal[] = "N";
		tet_basis::intgrt(lin, f1, stridex, stridey);
		GETRS(normal,tm,1,P.data(),tm,&ipiv(0),lin,tm,info);
		
		
	}
	
	
	
};

/*uncomment this for nodal 2/4 */
//namespace basis {
//   extern Array<tet_nodal_basis,1> tet;
//}


#endif

