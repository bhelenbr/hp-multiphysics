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

#define MAPPING

enum {vcommx = 1,vcommy,vprdx,vprdy};
enum {splain=1, scommx, scommy, sprdx, sprdy, sb23, stranslating1, stranslating2, capri_edge};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
/* INSTANTIATED IMPLICITLY BECAUSE OF REFERENCE ABOVE TO 2,3?? */
template<int ND> void mesh<ND>::getnewvrtxobject(int bnum, int idnum) {
   /* STRIP OF IDNTY INFO TO JUST GET TYPE */
   int type = idnum&0xFFFF;
   
#ifdef MAPPING
   /* SET UP MAPPING IF NEEDED */
   int mapping[4][2] = {{COMX_MASK,vcommx},{COMY_MASK,vcommy},{PRDX_MASK+COMX_MASK,vprdx},{PRDY_MASK+COMY_MASK,vprdy}};   
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
         vbdry[bnum] = new prdc_template<vcomm<mesh<ND> >,mesh>(idnum,*this);
      	break;
      }
      case vprdy: {
         prdc_template<vcomm<mesh<ND> >,mesh<ND> > *temp = new prdc_template<vcomm<mesh<ND> >,mesh<ND> >(idnum,*this);
         (*temp).setdir() = 1;
         vbdry[bnum] = temp;
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
   int mapping[5][2] = {{COMX_MASK,scommx},{COMY_MASK,scommy},{PRDX_MASK+COMX_MASK,sprdx},{PRDY_MASK+COMY_MASK,sprdy},{IFCE_MASK+CURV_MASK,scommy}};   
   for(int i=0;i<5;++i) {
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
#ifdef CAPRI
      case capri_edge: {
         sbdry[bnum] = new capri_edge<mesh<ND> >(idnum,*this);
         break;
      }
#endif
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
      }
};

class rtransandrot : public rfixd {
   private:
      FLT v0, amp, omega;
      FLT time;
   public:
      rtransandrot(int type,r_mesh &x, FLT v0in, FLT ampin, FLT omegain) : 
         rfixd(type,x), v0(v0in), amp(ampin), omega(omegain), time(0.0) {}
         
      void tadvance() {
         int n,vrt;
         FLT center[2], center1[2], dx[2],xp[2];
         FLT theta, theta1,dtheta;
         
         center[0] = v0*time;
         center[1] = amp*(1-cos(omega*time));
         theta = atan(omega*amp*sin(omega*time)/v0);
         
         

         time = time +1.0;
         
         center1[0] = v0*time;
         center1[1] = amp*(1-cos(omega*time));
         theta1 = atan(omega*amp*sin(omega*time)/v0);
         dtheta = theta1-theta;
                  
         for(int j=0;j<nsd();++j) {
            vrt = b().svrtx[sd(j)][0];
            xp[0] = b().vrtx[vrt][0]-center[0];
            xp[1] = b().vrtx[vrt][1]-center[1];         
            dx[0] = center1[0]-center[0] +xp[0]*cos(dtheta)-xp[1]*sin(dtheta) -xp[0];
            dx[1] = center1[1]-center[1] +xp[0]*sin(dtheta)+xp[1]*cos(dtheta) -xp[1];
            for(n=0;n<2;++n)
               b().vrtx[vrt][n] += dx[n];
         }
      }
};

/* THIS IS SET UP FOR RMESH PROBLEMS */
void r_mesh::getnewsideobject(int bnum, int idnum) {
   /* STRIP OF IDNTY INFO TO JUST GET TYPE */
   int type = idnum&0xFFFF;
   
   /* SET UP MAPPING IF NEEDED */
#ifdef MAPPING
   int mapping[10][2] = {{COMX_MASK,scommx},{COMY_MASK,scommy},{IFCE_MASK+CURV_MASK,scommy},{PRDX_MASK+COMX_MASK,sprdx},{PRDY_MASK+COMY_MASK,sprdy},
                        {OUTF_MASK,splain},{EULR_MASK,stranslating2},{FSRF_MASK,sb23},{SYMM_MASK,stranslating1},{INFL_MASK,stranslating2}};   
   for(int i=0;i<10;++i) {
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
#ifndef FOURTH
         sbdry[bnum] = new rgeneric<twotothree<mesh<2> >,rfixx,rfixy>(idnum,*this);
#else
         sbdry[bnum] = new rgeneric<twotothree<mesh<2> >,rfixx,rfixy,no_rfix,no_rfix>(idnum,*this);
#endif
         break;
      }
      case stranslating1: {
         sbdry[bnum] = new rtransandrot(idnum,*this,0.02,2.5,1./300.);
         break;
      }
      case stranslating2: {
         FLT dx[2] = {0.0,0.1};
         sbdry[bnum] = new rtranslating(idnum,*this,dx);
         break;
      }
      default: {
         *log << "Don't know this rboundary type:" << idnum << std::endl;
#ifndef FOURTH
         sbdry[bnum] = new rfixd(idnum,*this);
#else
         sbdry[bnum] = new rfarfield(idnum,*this);
#endif
         break;
      }
   }
   
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


   
   
