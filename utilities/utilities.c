#include<stdlib.h>
#include<string.h>
#include<stdio.h>

void *xmalloc(size_t a) {
	void *temp;

   if (a <= 0) return(NULL);
   
	temp = malloc(a);
	if (temp == NULL) {
		printf("couldn't allocate %ld\n",a);
		exit(1);
	}
	
/*	static int i =0;
   printf("xmalloc %d:  %d bytes\n",i++,a); */
   
	return(temp);
}

void number_str(char *out, const char *in, int n, int dgts) {
   int i,denom;
   char temp[100];
   
   strcpy(out,in);
   
   denom = 1;
   for(i=0;i<dgts-1;++i) denom*=10;
   
   for(i=0;i<dgts;++i) {
      if (n/denom == 0)
         strcat(out,"0");
      else break;
      
      denom /= 10;
   }
   
   if (n != 0) {
      sprintf(temp,"%d",n);
      strcat(out,temp);
   }
   
   return;
}

