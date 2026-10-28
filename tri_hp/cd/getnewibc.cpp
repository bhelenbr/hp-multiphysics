/*
 *  getnewibc_ins.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 2/6/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_cd.h"
#include "ibc_cd.h"
#include <cstring>
#include <math.h>
namespace ibc_cd {
	
#ifdef CONVDIFF
	
#if (CASE == 0)
	FLT forcing(FLT x,FLT y) { return(-cos(2.*M_PI*x));}
#else
	FLT forcing(FLT x,FLT y) {return(0.0);}
#endif
	FLT blayer = 0.0;
	FLT axext = AMP*cos(M_PI*LAM/180.0), ayext = AMP*sin(M_PI*LAM/180.0);
	FLT nuext = (1.0 - AMP);
	
	FLT f1(int n, FLT x, FLT y) {
		FLT nux = nuext*4.*M_PI*M_PI;
		FLT axx = axext*2.*M_PI;
		FLT xx = x*2.*M_PI;
		FLT yx = y*2.*M_PI;
		double eps = 0.05;
		
		switch(n) {
			case((0+CASE)%3):
				return(axx*sin(xx)/(axx*axx +nux*nux)
					   +nux*cos(xx)/(axx*axx +nux*nux)
					   +(nux > 0.0 ? blayer*exp(axx/nux*xx) : 0.0)
					   +startup*((sin(xx) +sin(100*xx))*(sin(yx)+sin(100*yx))));
			case((1+CASE)%3):
				return(startup*((sin(xx) +sin(100*xx))*(sin(yx)+sin(100*yx))));
			case((2+CASE)%3):
				if (x > 2.*eps)
					return(0.0);
				return(1+cos(M_PI*(x-eps)/eps));
		}
		return(0.0);
	}
	
	FLT df1d(int n, FLT x, FLT y) {
		FLT nux = nuext*4.*M_PI*M_PI;
		FLT axx = axext*2.*M_PI;
		FLT xx = x*2.*M_PI;
		switch(n) {
			case((0+CASE)%3):
				return(axx*cos(xx)/(axx*axx +nux*nux) 
					   -nuext*sin(x)/(axx*axx +nux*nux)
					   +(nux > 0.0 ? blayer*axx/nux*exp(axx/nux*xx) : 0.0));
			case(!CASE):
				return(0.0);
		}
		return(0.0);
	}
	
#endif
	
	class ibc_type {
	public:
		const static int ntypes = 2;
		enum ids {power,soi};
		const static char names[ntypes][40];
		static int getid(const char *nin) {
			int i;
			for(i=0;i<ntypes;++i) 
				if (!strcmp(nin,names[i])) return(i);
			return(-1);
		}
	};
	const char ibc_type::names[ntypes][40] = {"power","soi"};
	
}


init_bdry_cndtn *tri_hp_cd::getnewibc(std::string name) {
	std::string keyword,ibcname;
	init_bdry_cndtn *temp;
	int type;
	

	type = ibc_cd::ibc_type::getid(name.c_str());
	switch(type) {
		case(ibc_cd::ibc_type::power):
			temp = new ibc_cd::power_src;
			break;
		case(ibc_cd::ibc_type::soi):
			temp = new ibc_cd::soi_src(*this);
			break;
		default: {
			return(tri_hp::getnewibc(name));
		}
	}
	return(temp);
}





