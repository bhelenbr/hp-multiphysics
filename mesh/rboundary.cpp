/*
 *  mvpttobdry.cpp
 *  mblock
 *
 *  Created by helenbrk on Wed Oct 24 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "r_mesh.h"
#include<float.h>
#include<math.h>
#include<stdlib.h>

#define R 0.5

class rboat1 : public rfixd {
   public:
      rboat1(r_mesh &x,int type): rfixd(x,type) {}
      void tadvance() {
         int v0;
         for(int j=0;j<nsd();++j) {
            v0 = b().svrtx[sd(j)][0];
            b().vrtx[v0][0] -= 0.0;
            b().vrtx[v0][1] -= 0.5;
         }
      }
};

class rboat2 : public rfixd {
   public:
      rboat2(r_mesh &x,int type): rfixd(x,type) {}
       void tadvance() {
         int v0;
         for(int j=0;j<nsd();++j) {
            v0 = b().svrtx[sd(j)][0];
            b().vrtx[v0][0] -= 0.9;
         }
      }
};

/* THIS IS SET UP FOR RMESH PROBLEMS */
side_boundary* r_mesh::getnewsideobject(int type) {
   if (type & PRDX_MASK) 
      return(new rprdx(*this,type));
   else if (type & PRDY_MASK) 
      return(new rprdy(*this,type));
   else if (type & COMX_MASK) 
      return(new rcomx(*this,type));
   else if (type & COMY_MASK) 
      return(new rcomy(*this,type));
   else if (type & FSRF_MASK)
      return(new rboat1(*this,type));
   else if (type & EULR_MASK)
      return(new rboat2(*this,type)); 
   else return(new rfixd(*this,type));
}
