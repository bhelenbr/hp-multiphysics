/*
 *  mvpttobdry.cpp
 *  mblock
 *
 *  Created by helenbrk on Wed Oct 24 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "r_mesh.h"
#include "rblock.h"
#include "blocks.h"
#include<float.h>
#include<math.h>
#include<stdlib.h>

block * blocks::getnewblock(int type) {
   return(new rblock<r_mesh>);
}

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
void mesh::getnewvrtxobject(int i, int type) {
/*
   if (type & (PRDX_MASK +PRDY_MASK +COMX_MASK +COMY_MASK) ) {
      vbdry[i] = new vcom_boundary(*this,type);
      return;
   }
*/

   if (type & (PRDX_MASK)) { 
      vbdry[i] = new prdc_template<vcomm<mesh>,mesh>(type,*this);
      return;
   }
   vbdry[i] = new vrtx_template<mesh>(type,*this);

   return;
}


/* FUNCTION TO CREATE BOUNDARY OBJECTS */
void mesh::getnewsideobject(int i, int type) {
/*
   if (type & PRDX_MASK) {
      sbdry[i] = new prdx_boundary(*this,type);
      return;
   }
   if (type & PRDY_MASK) {
      sbdry[i] = new prdy_boundary(*this,type);
      return;
   }
   if (type & (COMX_MASK+COMY_MASK)) {
      sbdry[i] = new scomm_boundary(*this,type); 
      return;
   }
*/
   sbdry[i] = new side_template<mesh>(type,*this);

   return;
}

class rboat1 : public rfixd {
   public:
      rboat1(int type,r_mesh &x): rfixd(type,x) {}
      void tadvance() {
         int v0;
         for(int j=0;j<nsd();++j) {
            v0 = b().svrtx[sd(j)][0];
            b().vrtx[v0][0] -= 0.0;
            b().vrtx[v0][1] -= 0.05;
         }
         v0 = b().svrtx[sd(nsd()-1)][1];
         b().vrtx[v0][0] -= 0.0;
         b().vrtx[v0][1] -= 0.05;
      }
};

class rboat2 : public rfixd {
   public:
       rboat2(int type,r_mesh &x): rfixd(type,x) {}
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
         rprdx *temp = new rprdx(type,*this);
         temp->setdir() = 0;
         sbdry[i] = temp;
         break;
      }
      else if (type & PRDY_MASK) {
         rprdy *temp = new rprdy(type,*this);
         temp->setdir() = 1;  
         sbdry[i] = temp;
         break;
      }
      else if (type & IFCE_MASK) {
         rcomm *temp = new rcomm(type,*this);
         sbdry[i] = temp;
         break;
      }
/*
      else if (type & COMX_MASK) {
         rcomm *temp = new rcomm(type,*this);
         sbdry[i] = temp;
         break;
      }
      else if (type & COMY_MASK) {
         rcomm *temp = new rcomm(type,*this);
         sbdry[i] = temp;
         break;
      }
*/
      else if (type & OUTF_MASK) {
         sbdry[i] = new rboat1(type,*this);
         break;
      }
/*
      else if (type & EULR_MASK) {
         sbdry[i] = new rboat2(type,*this); 
         break;
      }
*/
      else {
         sbdry[i] = new rfixd(type,*this);
         break;
      }
   } while(0);
   
   rbdry[i] = dynamic_cast<rbdry_interface *>(sbdry[i]);
   
   return;
}
