#include "utilities.h"
#include <stdio.h>
#include <new>
#include <fstream>

int main (int argc, const char * argv[]) {
   char out[100];
     
   number_str(out,"name",5,3);

   printf("%s\n",out);

   sharedmem mymem1(10*sizeof(int));
   sharedmem mymem2(10*sizeof(int));

   mymem1.checkout();
   
   int *myints = new (mymem1.data()) int;

   for(int i=0; i<10; ++i)
      myints[i] = i;
      
   for(int i=0;i<10;++i)
      printf("%d %d\n",i,myints[i]);
      
   mymem1.checkin();

   mymem2.reference(mymem1);
   mymem2.checkout();
   myints = new (mymem2.data()) int;

   for(int i=0;i<10;++i)
      printf("%d %d\n",i,myints[i]);
      
   mymem2.checkin();

   
   return(0);

}
