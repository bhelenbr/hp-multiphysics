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
#include "sblock.h"

/* OLD BOUNDARY TYPES */
#define FSRF_MASK (1<<0)
#define IFCE_MASK (1<<1)
#define INFL_MASK (1<<2)
#define OUTF_MASK (1<<3)
#define SYMM_MASK (1<<4)
#define EULR_MASK (1<<5)
#define PRDX_MASK (1<<6)   
#define PRDY_MASK (1<<7)
#define COMX_MASK (1<<8)
#define COMY_MASK (1<<9)
#define CURV_MASK (1<<10)

enum {vcommx = 1,vcommy,vprdx,vprdy};
enum {splain=1, scommx, scommy, sprdx, sprdy, sb23, stranslating1, stranslating2};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
/* INSTANTIATED IMPLICITLY BECAUSE OF REFERENCE ABOVE TO 2,3?? */
template<int ND> void mesh<ND>::getnewvrtxobject(int bnum, int idnum) {
   /* STRIP OF IDNTY INFO TO JUST GET TYPE */
   int type = idnum&0xFFFF;
   
#ifdef MAPPING
   /* SET UP MAPPING IF NEEDED */
   int mapping[4][2] = {{COMX_MASK,vcommx},{COMY_MASK,vcommy},{PRDX_MASK,vprdx},{PRDY_MASK,vprdy}};   
   for(int i=0;i<4;++i) {
      if (mapping[i][0] == type) {
         type = mapping[i][1];
         break;
      }
   }
#endif
   
   switch(type) {
      case vcommx: case vcommy: {
         vbdry[bnum] = new vcomm<mesh<ND> >(idnum,*this);
      	break;
      }
      case vprdx: {
         prdc_template<vcomm<mesh<ND> >,mesh<ND> > *temp = new prdc_template<vcomm<mesh<ND> >,mesh<ND> >(idnum,*this);
         (*temp).setdir() = 1;
         vbdry[bnum] = temp;
      	break;
      }
      case vprdy: {
         vbdry[bnum] = new prdc_template<vcomm<mesh<ND> >,mesh>(idnum,*this);
         break;
      }
      default: {
        *log << "Don't know this vertex type\n" << std::endl;
        vbdry[bnum] = new vrtx_template<mesh<ND> >(idnum,*this);
        break;
      }
   }
   return;
}


/* FUNCTION TO CREATE BOUNDARY OBJECTS */
template<int ND> void mesh<ND>::getnewsideobject(int bnum, int idnum) {

/* STRIP OF IDNTY INFO TO JUST GET TYPE */
   int type = idnum&0xFFFF;
   
#ifdef MAPPING
   /* SET UP MAPPING IF NEEDED */
   int mapping[4][2] = {{COMX_MASK,scommx},{COMY_MASK,scommy},{PRDX_MASK,sprdx},{PRDY_MASK,sprdy}};   
   for(int i=0;i<4;++i) {
      if (mapping[i][0] == type) {
         type = mapping[i][1];
         break;
      }
   }
#endif
   
   switch(type) {
      case splain: {
         sbdry[bnum] = new side_template<mesh<ND> >(idnum,*this);
         break;
      }
      case scommx: case scommy: {
         sbdry[bnum] = new scomm<mesh<ND> >(idnum,*this); 
         break;
      }
      case sprdx: {
         sbdry[bnum] = new prdc_template<scomm<mesh<ND> >,mesh<ND> >(idnum,*this);
         break;
      }
      case sprdy: {
         prdc_template<scomm<mesh<ND> >,mesh<ND> > *temp = new prdc_template<scomm<mesh<ND> >,mesh<ND> >(idnum,*this);
         (*temp).setdir() = 1;
         sbdry[bnum] = temp;
         break;
      }
      case sb23: {
         if (ND == 3)
            sbdry[bnum] = new threetotwo<mesh<ND> >(idnum,*this);
         else
            sbdry[bnum] = new twotothree<mesh<ND> >(idnum,*this);
         break;
      }
      default: {
         *log << "not sure what kind of sboundary this is for mesh: " << idnum << std::endl;
         sbdry[bnum] = new side_template<mesh<ND> >(idnum,*this);
      }
   }

   return;
}


class rtranslating : public rfixd {
   private:
      FLT dx[2];
   public:
      rtranslating(int type,r_mesh &x, FLT dxin[2]): rfixd(type,x) { 
         for(int n=0;n<2;++n)
            dx[n] = dxin[n];
      }
      void tadvance() {
         int n,v0;
         for(int j=0;j<nsd();++j) {
            v0 = b().svrtx[sd(j)][0];
            for(n=0;n<2;++n)
               b().vrtx[v0][n] += dx[n];
         }
         v0 = b().svrtx[sd(nsd()-1)][1];
         for(n=0;n<2;++n)
         b().vrtx[v0][n] += dx[n];
      }
};

/* THIS IS SET UP FOR RMESH PROBLEMS */
void r_mesh::getnewsideobject(int bnum, int idnum) {
   /* STRIP OF IDNTY INFO TO JUST GET TYPE */
   int type = idnum&0xFFFF;
   
   /* SET UP MAPPING IF NEEDED */
#ifdef MAPPING
   int mapping[7][2] = {{COMX_MASK,scommx},{COMY_MASK,scommy},{PRDX_MASK,sprdx},{PRDY_MASK,sprdy},
                        {OUTF_MASK,stranslating1},{EULR_MASK,stranslating2},{FSRF_MASK,sb23}};   
   for(int i=0;i<7;++i) {
      if (mapping[i][0] == type) {
         type = mapping[i][1];
         break;
      }
   }
#endif

   switch(type) {
      case splain: {
         sbdry[bnum] = new rfixd(idnum,*this);
         break;
      }
      case scommx: case scommy: {
         rcomm *temp = new rcomm(idnum,*this);
         sbdry[bnum] = temp;
         break;
      }
      case sprdx: {
         rprdx *temp = new rprdx(idnum,*this);
         temp->setdir() = 0;
         sbdry[bnum] = temp;
         break;
      }
      case sprdy: {
         rprdy *temp = new rprdy(idnum,*this);
         temp->setdir() = 1;  
         sbdry[bnum] = temp;
         break;
      }
      case sb23: {
         sbdry[bnum] = new rgeneric<twotothree<mesh<2> >,rfixx,rfixy>(idnum,*this);
         break;
      }
      case stranslating1: {
         FLT dx[2] = {0.0,-0.05};
         sbdry[bnum] = new rtranslating(idnum,*this,dx);
         break;
      }
      case stranslating2: {
         FLT dx[2] = {-0.9,0.0};
         sbdry[bnum] = new rtranslating(idnum,*this,dx);
         break;
      }
      default: {
         *log << "Don't know this rboundary type:" << idnum << std::endl;
         sbdry[bnum] = new rfixd(idnum,*this);
         break;
      }
   };
   
   rbdry[bnum] = dynamic_cast<rbdry_interface *>(sbdry[bnum]);
   
   return;
}


block * blocks::getnewblock(int type) {
   enum {twoD=1, threeD, rblock2D, capri};
   switch(type) {
      case (twoD):
         return(new mgrid<mesh<2> >);
      case(threeD):
         return(new mgrid<mesh<3> >);
      case(rblock2D):
         return(new rblock<r_mesh>);
      case(capri):
         return(new capriblock);
   }
   
   return(0);
}