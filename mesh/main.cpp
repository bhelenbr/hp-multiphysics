#include<cstdio>
#include"blocks.h"
#include<iostream>

/*	DEFINES/THINGS TO TIDDLE WITH */
/*	r_mesh.h: FOURTH: SWITCH FROM LAPLACIAN TO BIHARMONIC */
/*	blocks.cpp: GEOMETRIC: SWITCH FOR DETERMING K_IJ IN MULTIGRID (USE GEOMETRIC)*/
/* mesh.h: USEOLDBTYPE: remaps boundary identifies from my old types to new */
/*	r_mesh.h: FIX?_MASK: DETERMINES DIRICHLET MESH MOVMENT BC'S */
/* r_mesh.h: FIX2?_MASK: DETERMINES DIRICHLET B.C.'S on Nabla^2 for Biharmonic */
/* mesh.h: ?DIR_MP: mask for communication boundaries */
/*	r_mesh.cpp: different perturb functions for deformed mesh */
/*	r_mesh.cpp: vnn/fadd - jacobi/multigrid parameters */

#include<cstring>
#include<math.h>
#include<utilities.h>
#include<time.h>

FLT center;
FLT amp;

int main(int argc, char *argv[]) {
   int i,step = 0;
   clock_t cpu_time;
   class blocks z;
   char outname[100];
   char *inname[2];
   char *out2[2];

/*	START CLOCK TIMER */
   clock();
   inname[0] = "/Volumes/work/helenbrk/Codes/spectral_hp/build/mesh068.1";
   inname[1] = "/Volumes/work/helenbrk/Codes/spectral_hp/build/mesh068.0";  
   z.init(2,1,inname,easymesh,10.0);
   z.out_mesh("mesh",tecplot);
   return(0);  
   
   
/*	THIS DEFORMS A MESH */
//	inname[0] = "../../grids/TIM/tim";
//	inname[0] = "/Network/Servers/shelob.camp.clarkson.edu/home/helenbrk/codes/grids/WAVE/PRDC/wave5";

   inname[0] = "/Volumes/work/helenbrk/Codes/grids/BOAT/boat2";
   
   z.init(1,3,inname,easymesh,10.0);
   
   z.out_mesh("begin",tecplot);

   for(step = 1; step<=1;++step) {
      center += step/FLT(10);  // FOR THE MOVING CYLINDER PROBLEM
      amp = 0.25*step/10.;  // FOR THE DEFORMING SURFACE 
      z.ksrc();
      z.perturb();
      

      for(i=0;i<200;++i) {
         z.cycle(2);
         printf("%d ",i);
         z.maxres();
         printf("\n");
      }
      
      z.out_mesh("deformed",tecplot);
      
      z.restructure(0.66);
      number_str(outname, "test", step, 2);
      z.out_mesh(outname);
   }
   
   out2[0] = "bamp0.375";
   out2[1] = "tamp0.375";
   z.out_mesh(out2,easymesh);
   
   cpu_time = clock();
   printf("that took %ld cpu time\n",cpu_time);

   return(0);
}
