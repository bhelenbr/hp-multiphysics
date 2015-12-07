#include "tet_basis.h"
#include <stdio.h>
#include <myblas.h>
#include <time.h>

#define P 4
#define GP P+1

int main (int argc, const char * argv[]) {
	clock_t cpu_time;
//	Array<double,2> f(5,5);
//	Array<double,1> lin(3+3*3+3);
	
	/* CREATE 1 BASIS FOR GENERAL USE */
	basis::tet.resize(1);
	basis::tet(0).initialize(P,GP);

	basis::tet(0).outputstuff(12);
	// case 7: mass matrix
	
	
	
	clock();
//	double uht[MXTM];
//	double u[GP][GP][GP];
//	double dx[GP][GP][GP];
//	double dy[GP][GP][GP];
//	double dz[GP][GP][GP];
//
//	for (int i = 0; i < 100; ++i) {
//		basis::tet(0).proj(uht,u[0],dx[0],dy[0],dz[0],GP*GP,GP);
//		basis::tet(0).intgrtrst(uht,dx[0],dy[0],dz[0],GP*GP,GP); 
//		basis::tet(0).intgrt(uht,dx[0],GP*GP,GP); 
//	}

	cpu_time = clock();
	printf("#that took %ld cpu time\n",cpu_time);      

	return 0;
}
