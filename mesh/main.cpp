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

FLT center;
/* THIS NEEDS TO BE MODIFED DEPENDING ON DESIRED RESULT */
/* SIDE BASED REFINEMENT CRITERION */
void mesh::density() {
   int i;
   
//   for(i=0;i<nvrtx;++i)
//      vlngth[i] = 1.05*3.1415/20.0*(1.0  - 0.875*exp(-((vrtx[i][0] - center)*(vrtx[i][0] -center) + vrtx[i][1]*vrtx[i][1]) +0.5*0.5));

//   for(i=0;i<nvrtx;++i)
//      vlngth[i] = 0.5;

      
   for(i=0;i<nside;++i)
      fltwk[i] = (vlngth[svrtx[i][0]] +vlngth[svrtx[i][1]])/
      (2.*distance(svrtx[i][0],svrtx[i][1]));
   
   return;
}

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

#ifdef SKIP
void r_mesh::perturb(int step) {
   int i,j,sind,v0;

   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & (FSRF_MASK +IFCE_MASK)) {
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            vrtx[v0][1]  = +0.25*step/10.*cos(2.*M_PI*vrtx[v0][0]);
         }
         v0 = svrtx[sind][1];
         vrtx[v0][1]  = +0.25*step/10.*cos(2.*M_PI*vrtx[v0][0]);
      }
   }

   return;
}
#endif

#ifdef SKIP
void r_mesh::perturb(int step) {
   int i,j,sind,v0;

   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & INFL_MASK) {
         for(j=0;j<sbdry[i].num;++j) {
            sind = sbdry[i].el[j];
            v0 = svrtx[sind][0];
            vrtx[v0][1]  *= 2.0;
         }
         v0 = svrtx[sind][1];
         vrtx[v0][1]  *= 2.0;
      }
   }

   return;
}
#endif

int main(int argc, char *argv[]) {
   int i,step;
   char name[100];

/* THIS READS IN A MESH THEN COARSENS IT & OUTPUTS */
   class mesh x, y;
   
//   x.in_mesh("../../grids/TIM/tim",easymesh,10.0);
//   y = x;
//   y.out_mesh("test");
//   exit(0);
//   x.rebay(0.8);
//   x.bcinfo();
//   x.out_mesh("test");
//   exit(0);
   x.in_mesh("../../data264",gambit,2.0);
   x.out_mesh("test");
   x.swap();
   x.out_mesh("test1");
   exit(0);
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
   class blocks z;

  z.init(1,3,"../../grids/TIM/tim",easymesh,100.0);
//  z.init(1,3,"test7",easymesh,5.0);

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
   for(step = 1; step<=3;++step) {
      z.ksrc();
      z.perturb(step);
      for(i=0;i<100;++i)
         z.cycle(2);
      z.restructure(0.66);
      
      strcpy(name,"test");
      sprintf(&name[4],"%d",step);
      z.out_mesh(name);
   }
#endif

   return(0);
}
