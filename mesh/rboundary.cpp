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
void r_mesh::getnewsideobject(int i, int type) {
   
   do {
      if (type & PRDX_MASK) {
         sbdry[i] = new rprdx(*this,type);
         break;
      }
      else if (type & PRDY_MASK) {
         sbdry[i] = new rprdy(*this,type);
         break;
      }
      else if (type & COMX_MASK) {
         rcomm *temp = new rcomm(*this,type);
         temp->setphase(0);
         sbdry[i] = temp;
         break;
      }
      else if (type & COMY_MASK) {
         rcomm *temp = new rcomm(*this,type);
         temp->setphase(0);
         sbdry[i] = temp;
         break;
      }
      else if (type & FSRF_MASK) {
         sbdry[i] = new rboat1(*this,type);
         break;
      }
      else if (type & EULR_MASK) {
         sbdry[i] = new rboat2(*this,type); 
         break;
      }
      else {
         sbdry[i] = new rfixd(*this,type);
         break;
      }
   } while(0);
   
   rbdry[i] = dynamic_cast<rbdry_interface *>(sbdry[i]);
   
   return;
}
