/*
 *  main.cpp
 *  planar++
 *
 *  Created by helenbrk on Tue Oct 09 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */


#include "blocks.h"
#include"osdepend.h"
#include<stdio.h>
#include<utilities.h>
#include<time.h>

#ifdef PV3
extern "C" int MAIN_(int argc, char **argv);

int MAIN_(int argc, char**argv) {
#else
int main(int argc, char **argv) {
#endif
   class blocks myblock;
   clock_t cpu_time;

   /* START CLOCK TIMER */
   clock();

   myblock.init(argv[1]);
   myblock.go();
   
   cpu_time = clock();
   printf("that took %ld cpu time\n",cpu_time);
   
   return(0);


   /* FOR TEST ADAPTATION */
   myblock.init(argv[1]);
   myblock.output(10,tecplot);
   
   return(0);

}



