#include "tri_basis.h"
#include <stdio.h>
#include <myblas.h>
#include <utilities.h>
#include <time.h>

// #define OP_COUNT
#ifdef OP_COUNT
#include <CHUD/CHUD.h>
#endif

#define P 1
#define GP P+1
#define EP 0

int main (int argc, const char * argv[]) {
	clock_t bgn_time, end_time;
	
	/* CREATE 1 BASIS FOR GENERAL USE */
	tri_basis<P,EP> tri;	
	tri.initialize();
	
	bgn_time = clock();
	double uht[tri.tm()];
	double u[GP][GP];
	double dx[GP][GP];
	double dy[GP][GP];
	
#ifdef OP_COUNT
	chudInitialize(); 
	chudAcquireRemoteAccess(); 
	chudStartRemotePerfMonitor("fpucount"); 
#endif
	
	for (int i = 0; i < 10000; ++i) {
		// basis::tri(0).proj(uht,u[0],dx[0],dy[0],GP);
		// basis::tri(0).intgrtrs(uht,dx[0],dy[0],GP); 
		tri.intgrt(uht,dx[0],GP); 
	}
	
	
#ifdef OP_COUNT
	chudStopRemotePerfMonitor(); 
	chudReleaseRemoteAccess();
#endif
	
	end_time = clock();
	printf("#that took %lf seconds\n",static_cast<double>(end_time - bgn_time)/ CLOCKS_PER_SEC);\
	
	return 0;
}
