#include "utilities.h"
#include<stdio.h>

int main (int argc, const char * argv[]) {
    char out[100];
        
    number_str(out,"name",5,3);
    
    printf("%s\n",out);
    
    return(0);
    
}
