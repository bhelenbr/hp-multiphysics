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

#define WAVE

FLT center;
/* THIS NEEDS TO BE MODIFED DEPENDING ON DESIRED RESULT */
/* SIDE BASED REFINEMENT CRITERION */
void mesh::density() {
   int i;
   
//   for(i=0;i<nvrtx;++i)
//      vlngth[i] = 1.05*3.1415/20.0*(1.0  - 0.875*exp(-((vrtx[i][0] - center)*(vrtx[i][0] -center) + vrtx[i][1]*vrtx[i][1]) +0.5*0.5));

//   for(i=0;i<nvrtx;++i)
//      vlngth[i] = 1./32.;
      
   for(i=0;i<nside;++i)
      fltwk[i] = (vlngth[svrtx[i][0]] +vlngth[svrtx[i][1]])/
      (2.*distance(svrtx[i][0],svrtx[i][1]));
   
   return;
}

#ifdef CYLINDER
void r_mesh::perturb(int step) {
   int i,j,sind,v0;

   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & (EULR_MASK)) {
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            vrtx[v0][0]  += step/FLT(10);
         }
      }
   }
   center += step/FLT(10);
   
   return;
}
#endif

#ifdef WAVE
FLT amp;

void r_mesh::perturb(int step) {
   int i,j,sind,v0;

   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & (FSRF_MASK +IFCE_MASK)) {
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            vrtx[v0][1]  = +0.25*step/20.*sin(2.*M_PI*(vrtx[v0][0]-0.0));
         }
         v0 = svrtx[sind][1];
         vrtx[v0][1]  = +0.25*step/20.*sin(2.*M_PI*(vrtx[v0][0]-0.0));
      }
   }
   amp = 0.25*step/10.;

   return;
}
#endif

int main(int argc, char *argv[]) {
   int i,step = 0;
   char outname[100];
   clock_t cpu_time;
   class mesh x, y;
   class blocks z;
   char *inname[2];
   char *out2[2];

/*	START CLOCK TIMER */
   clock();
   
/*   
   x.in_mesh("/Network/Servers/shelob.camp.clarkson.edu/home/helenbrk/codes/grids/WAVE/PRDC/wave5",easymesh,5.0);//   y = x;
*/
/*
   x.in_mesh("cstart",easymesh,8);
   x.density();
   x.swap();
   x.yaber(1/0.66);
   x.rebay(0.66);
   x.setbcinfo();
   x.out_mesh("start");
   exit(10);
*/

//   x.in_mesh("start",easymesh);
//   x.out_mesh("start",tecplot);

//   y.out_mesh("test");
//   exit(0);
//   x.rebay(0.8);
//   x.setbcinfo();
//   x.out_mesh("test");
//   exit(0);
/*
   x.in_mesh("../../data264",gambit,2.0);
   x.out_mesh("test");
   x.swap();
   x.out_mesh("test1");
   exit(0);
*/
#ifdef SKIP
      y.in_mesh("../../grids/TIM/tim",easymesh,100.0);
      y.density();
      y.yaber(1.5);
      y.out_mesh("test1");
      y.rebay(0.66);
      y.out_mesh("test");
      exit(0);
#endif

/*	THIS DEFORMS A MESH */

   
//	inname[0] = "../../grids/TIM/tim";
//	inname[0] = "/Network/Servers/shelob.camp.clarkson.edu/home/helenbrk/codes/grids/WAVE/PRDC/wave5";
   
   inname[0] = "/Volumes/work/helenbrk/Codes/grids/WIND/PRDC/bot1";
   inname[1] = "/Volumes/work/helenbrk/Codes/grids/WIND/PRDC/top1";
   
   z.init(2,3,inname,easymesh,50.0);

#define NOT_ONE

#ifdef ONE
/*	All in one shot */
   z.out_mesh("start");
   z.ksrc();
   z.perturb(5);
   for(i=0;i<100;++i)
      z.cycle(2);
   z.out_mesh("deformed");
   z.restructure(0.66);
#else
/* For Incremental Changes */
   for(step = 1; step<=30;++step) {
      z.ksrc();
      z.perturb(step);
#ifdef SKIP
      for(i=0;i<40;++i) {
         z.cycle(2);
         printf("%d ",i);
         z.maxres();
         printf("\n");
      }
#endif
      z.restructure(0.66);
      number_str(outname, "test", step, 2);
      z.out_mesh(outname);
   }
   
   out2[0] = "bamp0.375";
   out2[1] = "tamp0.375";
   z.out_mesh(out2,easymesh);
   
#endif
   cpu_time = clock();
   printf("that took %ld cpu time\n",cpu_time);

   return(0);
}
