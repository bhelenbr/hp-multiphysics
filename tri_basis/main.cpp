#include "hpbasis.h"
#include <stdio.h>
#include <myblas.h>
#include <utilities.h>
#include <time.h>

#define P 4
#define GP P+2

int main (int argc, const char * argv[]) {
   clock_t cpu_time;
   hpbasis b;

   b.initialize(P,GP);
   clock();
   double uht[MXTM];
   double u[GP][GP];
   double dx[GP][GP];
   double dy[GP][GP];
   
   for (int i = 0; i < 4000; ++i) {
      b.proj(uht,u[0],dx[0],dy[0],GP);
     	b.intgrtrs(uht,dx[0],dy[0],GP); 
      b.intgrt(uht,dx[0],GP); 
   }
   cpu_time = clock();
   printf("#that took %ld cpu time\n",cpu_time);      

   return 0;
}
