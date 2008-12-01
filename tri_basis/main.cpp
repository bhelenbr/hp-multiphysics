#include "tri_basis.h"
#include <stdio.h>
#include <myblas.h>
#include <utilities.h>
#include <time.h>

#define P 4
#define GP P+2

#define OP_COUNT
#include <CHUD/CHUD.h>

int main (int argc, const char * argv[]) {
    clock_t cpu_time;
    
    /* CREATE 1 BASIS FOR GENERAL USE */
    basis::tri.resize(1);
    basis::tri(0).initialize(P,GP);
    
    clock();
    double uht[MXTM];
    double u[GP][GP];
    double dx[GP][GP];
    double dy[GP][GP];
	
#ifdef OP_COUNT
	chudInitialize(); 
	chudAcquireRemoteAccess(); 
	chudStartRemotePerfMonitor("fpucount"); 
#endif
	
    for (int i = 0; i < 10000; ++i) {
        basis::tri(0).proj(uht,u[0],dx[0],dy[0],GP);
      	basis::tri(0).intgrtrs(uht,dx[0],dy[0],GP); 
        basis::tri(0).intgrt(uht,dx[0],GP); 
    }
	
	
#ifdef OP_COUNT
	chudStopRemotePerfMonitor(); 
	chudReleaseRemoteAccess();
#endif
		
    cpu_time = clock();
    printf("#that took %ld cpu time\n",cpu_time);        

    return 0;
}
