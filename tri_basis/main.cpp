#include "tri_basis.h"
#include <stdio.h>
#include <myblas.h>
#include <time.h>

// #define OP_COUNT
#ifdef OP_COUNT
#include <CHUD/CHUD.h>
#endif

#define P 4
#define EP 1
#define GP P+1+EP


int main (int argc, const char * argv[]) {
	clock_t bgn_time, end_time;
	
//#define SKIP
#ifdef SKIP
	tri_basis_array<1> bases;
	Array<FLT,1> f;
	Array<FLT,1> rslt;
	f.resize(bases(2)->gpx());
	rslt.resize(bases(2)->sm()+2);
	f = 0.0;
	bases(2)->intgrt1d(rslt.data(),f.data());
	std::cout << rslt << std::endl;
#endif
	
#ifdef OP_COUNT
	chudInitialize(); 
	chudAcquireRemoteAccess(); 
	chudStartRemotePerfMonitor("fpucount"); 
#endif

#define TEST	
#ifdef TEST
	/* CREATE 1 BASIS FOR GENERAL USE */
	tri_basis<P,EP> tri;	
	tri.initialize();
	
	bgn_time = clock();
	double uht[tri.tm()];
	double u[GP][GP];
	double dx[GP][GP];
	double dy[GP][GP];
	
	for (int i = 0; i < 10000; ++i) {
		basis::tri(0).proj(uht,u[0],dx[0],dy[0],GP);
		basis::tri(0).intgrtrs(uht,dx[0],dy[0],GP); 
		tri.intgrt(uht,dx[0],GP); 
	}
#endif
	
#ifdef OP_COUNT
	chudStopRemotePerfMonitor(); 
	chudReleaseRemoteAccess();
#endif
	
	end_time = clock();
	printf("#that took %lf seconds\n",static_cast<double>(end_time - bgn_time)/ CLOCKS_PER_SEC);\
	
	return 0;
}
