/*
 *  hp_boundary.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/2/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp.h"
#include "hp_boundary.h"

#define NO_MPDEBUG

FLT hp_vrtx_bdry::dummy;

void hp_side_bdry::smatchsolution_rcv(int phi, FLT *sdata, int bgn, int end, int stride) {
   
   if (!base.is_comm()) return;
      
   /* CAN'T USE sfinalrcv BECAUSE OF CHANGING SIGNS */
   int j,k,m,n,count,countdn,countup,offset,sind,sign;
   FLT mtchinv;
   
   /* ASSUMES REVERSE ORDERING OF SIDES */
   /* WON'T WORK IN 3D */
   
   int matches = 1;
   
   int bgnsign = (bgn % 2 ? -1 : 1);
   
   /* RELOAD FROM BUFFER */
   /* ELIMINATES V/S/F COUPLING IN ONE PHASE */
   /* FINALRCV SHOULD BE CALLED F,S,V ORDER (V HAS FINAL AUTHORITY) */   
   for(m=0;m<base.matches();++m) {   
      if (base.matchphase(m) != phi) continue;
      
      ++matches;
      
      int ebp1 = end-bgn+1;
      countdn = (base.nel-1)*ebp1*x.NV;
      countup = 0;
      for(j=0;j<base.nel;++j) {
         sign = bgnsign;
         for(k=0;k<ebp1;++k) {
            for(n=0;n<x.NV;++n)
               base.fsndbuf(countup++) += sign*base.frcvbuf(m,countdn +k);
            sign *= -1;
         }
         countdn -= ebp1*x.NV;
      }
   }
   
   if (matches > 1) {
      mtchinv = 1./matches;

#ifdef MPDEBUG
      *sim::log << "finalrcv"  << base.idnum << " " << base.is_frst() << std::endl;
#endif
      count = 0;
      for(j=0;j<base.nel;++j) {
         sind = base.el(j);
         offset = (sind*stride +bgn)*x.NV;
         for (k=bgn;k<=end;++k) {
            for(n=0;n<x.NV;++n) {
               sdata[offset] = base.fsndbuf(count++)*mtchinv;
#ifdef MPDEBUG
               *sim::log << "\t" << sdata[offset] << std::endl;
#endif
               ++offset;
            }
         }

      }
   }
   return;
}

void hp_side_bdry::copy_data(const hp_side_bdry &bin) {
   
   if (!curved || !x.sm0) return;
   
   for(int i=0;i<sim::nadapt; ++i)
      crvbd(i)(Range(0,base.nel-1),Range::all()) = bin.crvbd(i)(Range(0,base.nel-1),Range::all());
}

void hp_side_bdry::init(input_map& inmap) {
   int i,coarse;
   std::string keyword;
   std::istringstream data;
   std::string filename;
   
   dxdt.resize(x.log2pmax,base.maxel);
   
   keyword = base.idprefix + ".curved";
   inmap.getwdefault(keyword,curved,false);

   keyword = base.idprefix + ".coarse";
   inmap.getwdefault(keyword,coarse,0);
   
   keyword = base.idprefix + ".coupled";
   inmap.getwdefault(keyword,coupled,false);
   
   if (curved && !coarse) {
      crv.resize(base.maxel,x.sm0);
      for(i=1;i<sim::nhist+1;++i)
         crvbd(i).resize(base.maxel,x.sm0);
      crvbd(0).reference(crv);
   }

   return;
}

void hp_side_bdry::output(std::ostream& fout, tri_hp::filetype typ,int tlvl) {
   int j,m,n;
   
   switch(typ) {
      case(tri_hp::text):
         fout << base.idprefix << " " << mytype << std::endl;
         if (curved) {
            fout << "p0: " << x.p0 << std::endl;
            
            for(j=0;j<base.nel;++j) {
               for(m=0;m<x.sm0;++m) {
                  for(n=0;n<mesh::ND;++n)
                     fout << crvbd(tlvl)(j,m)(n) << ' ';
                  fout << std::endl;
               }
            }
         }
         
      default:
         break;
   }
   return;
}
   
void hp_side_bdry::input(ifstream& fin,tri_hp::filetype typ,int tlvl) {
   int j,m,n,pmin;
   std::string idin, mytypein;

   switch(typ) {
      case(text):
         fin >> idin >> mytypein;
         if (curved) {
            fin.ignore(80,':');
            fin >>  pmin;
            for(j=0;j<base.nel;++j) {
               for(m=0;m<pmin -1;++m) {
                  for(n=0;n<mesh::ND;++n)
                     fin >> crvbd(tlvl)(j,m)(n);
               }
               for(m=pmin-1;m<x.sm0;++m) {
                  for(n=0;n<mesh::ND;++n)
                     crvbd(tlvl)(j,m)(n) = 0.0;
               }
            }
         }
         break;
         
      default:
         break;
   }
}

void hp_side_bdry::curv_init(int tlvl) {
   int i,j,m,n,v0,v1,sind,info;
   TinyVector<FLT,2> pt;
   char uplo[] = "U";

   if (!curved) return;
   
//   for(bind=0;bind<nvbd;++bind) {
//      if (!(vbdry[bind].type&MVPT_MASK)) continue;
//      
//      v0 = vbdry[bind].el[0];
//      x = vrtx(v0)(0);
//      y = vrtx(v0)(1);
//      mvpttobdry(vbdry[bind].type,x,y);
//      vrtx(v0)(0) = x;
//      vrtx(v0)(1) = y;
//   }

   /* MESS WITH END VERTICES? */
   for(j=0;j<base.nel;++j) {
      sind = base.el(j);
      v0 = x.sd(sind).vrtx(0);
      base.mvpttobdry(j,0.0, x.vrtx(v0));
   }
   v0 = x.sd(sind).vrtx(1);
   base.mvpttobdry(base.nel-1,1.0, x.vrtx(v0));

   if (basis::tri(x.log2p).p == 1) return;
      
   /*****************************/
   /* SET UP HIGHER ORDER MODES */
   /*****************************/
   for(j=0;j<base.nel;++j) {
      sind = base.el(j);

      v0 = x.sd(sind).vrtx(0);
      v1 = x.sd(sind).vrtx(1);
      
      for(n=0;n<mesh::ND;++n) 
         basis::tri(x.log2p).proj1d(x.vrtx(v0)(n),x.vrtx(v1)(n),&x.crd(n)(0,0));
   
      for(i=0;i<basis::tri(x.log2p).gpx;++i) {
         pt(0) = x.crd(0)(0,i);
         pt(1) = x.crd(1)(0,i);
         base.mvpttobdry(j,basis::tri(x.log2p).gx(i,2),pt);
         x.crd(0)(0,i) -= pt(0);
         x.crd(1)(0,i) -= pt(1);
      }
      
      for(n=0;n<mesh::ND;++n) {
         basis::tri(x.log2p).intgrt1d(&x.cf(n,0),&x.crd(n)(0,0));
         DPBTRS(uplo,basis::tri(x.log2p).sm,basis::tri(x.log2p).sbwth,1,&basis::tri(x.log2p).sdiag1d(0,0),basis::tri(x.log2p).sbwth+1,&x.cf(n,2),basis::tri(x.log2p).sm,info);
      
         for(m=0;m<basis::tri(x.log2p).sm;++m)
            crvbd(tlvl)(j,m)(n) = -x.cf(n,m+2);
      }
      
      /* TEST FOR A CIRCLE 
      basis::tri(x.log2p).ptvalues1d(0.0);
      basis::tri(x.log2p).ptprobe1d(mesh::ND,&pt(0),&x.cht(0,0),MXTM);
      *sim::log << pt << ' ' << pt(0)*pt(0) +pt(1)*pt(1) << std::endl;
      */

   }
   return;
}

void hp_side_bdry::findbdrypt(const TinyVector<FLT,2> xp,int &bel,FLT &psi) {
   int sind,v0,v1,iter;
   FLT dx,dy,ol,roundoff,dpsi;
   TinyVector<FLT,2> pt;
      
   base.findbdrypt(xp,bel,psi);
   if (!curved) {
      basis::tri(x.log2p).ptvalues1d(psi);
      return;
   }
   
   sind = base.el(bel);
   v0 = x.sd(sind).vrtx(0);
   v1 = x.sd(sind).vrtx(1);
   dx = x.vrtx(v1)(0) - x.vrtx(v0)(0);
   dy = x.vrtx(v1)(1) - x.vrtx(v0)(1);
   ol = 2./(dx*dx +dy*dy);
   dx *= ol;
   dy *= ol;
   
   /* FIND PSI SUCH THAT TANGENTIAL POSITION ALONG LINEAR SIDE STAYS THE SAME */
   /* THIS WAY, MULTIPLE CALLS WILL NOT GIVE DIFFERENT RESULTS */ 
   x.crdtocht1d(sind);
   
   iter = 0;
   roundoff = 10.0*EPSILON*(1.0 +(fabs(xp(0)*dx) +fabs(xp(1)*dy)));
   do {
      basis::tri(x.log2p).ptprobe1d(x.ND,pt.data(),psi,&x.cht(0,0),MXTM);
      dpsi = (pt(0) -xp(0))*dx +(pt(1) -xp(1))*dy;
      psi -= dpsi;
      if (iter++ > 100) {
         *sim::log << "#Warning: max iterations for curved side in bdry_locate type: " << base.idnum << " el: " << bel << " sind: " << sind << " loc: " << xp << " dpsi: " << dpsi << std::endl;
         break;
      }  
   } while (fabs(dpsi) > roundoff);
}

block::ctrl hp_side_bdry::tadvance(int excpt) {
   int stage = sim::substep +sim::esdirk;  
   
   if (x.p0 == 1 || !curved) return(block::stop);
   
   switch(excpt) {
      case(0): {
         if (stage) {
            /* BACK CALCULATE K TERM */
            for(int j=0;j<base.nel;++j) {
               for(int m=0;m<basis::tri(x.log2p).sm;++m) {
                  for(int n=0;n<mesh::ND;++n)
                     crvbd(stage+1)(j,m)(n) = (crvbd(0)(j,m)(n)-crvbd(1)(j,m)(n))*sim::adirk[sim::substep][sim::substep];
               }
            }
         }
         
         if (sim::substep == 0) {
            /* STORE TILDE W */
            for(int j=0;j<base.nel;++j) {
               for(int m=0;m<basis::tri(x.log2p).sm;++m) {
                  for(int n=0;n<mesh::ND;++n)
                     crvbd(1)(j,m)(n) = crvbd(0)(j,m)(n);
               }
            }
         }
         
         /* UPDATE TILDE W */
         for (int s=0;s<stage;++s) {         
            for(int j=0;j<base.nel;++j) {
               for(int m=0;m<basis::tri(x.log2p).sm;++m) {
                  for(int n=0;n<mesh::ND;++n)
                     crvbd(1)(j,m)(n) += sim::adirk[stage][s]*crvbd(s+2)(j,m)(n);
               }
            }
         }
         return(block::advance);
      }
   
      case(1): {
         /* CALCULATE MESH VELOCITY SOURCE TERMS ON BOUNDARY? */
         return(block::advance);
      }
      case(2): {
         /* EXTRAPOLATE SOLUTION OR MOVE TO NEXT TIME POSITION */
         if (!coupled) curv_init();
         return(block::advance);
      }
   }
   
   return(block::stop);
}





void tri_hp::vc0load(int phase, FLT *vdata) {
   int i;
         
   /* SEND COMMUNICATIONS TO ADJACENT MESHES */\
   for(i=0;i<nsbd;++i) 
      hp_sbdry(i)->vmatchsolution_snd(phase,vdata);
   for(i=0;i<nvbd;++i)
      hp_vbdry(i)->vmatchsolution_snd(phase,vdata);
   
   for(i=0;i<nsbd;++i)
      sbdry(i)->comm_prepare(phase);
   for(i=0;i<nvbd;++i)
      vbdry(i)->comm_prepare(phase);
   
   return;
}
int tri_hp::vc0wait_rcv(int phase, FLT *vdata) {
   int stop = 1;
   int i;
   
   for(i=0;i<nsbd;++i)
      stop &= sbdry(i)->comm_wait(phase);
   for(i=0;i<nvbd;++i)
      stop &= vbdry(i)->comm_wait(phase);
      

   for(i=0;i<nsbd;++i) 
      hp_sbdry(i)->vmatchsolution_rcv(phase,vdata);
   for(i=0;i<nvbd;++i)
      hp_vbdry(i)->vmatchsolution_rcv(phase,vdata);
      
   return(stop);
}

int tri_hp::vc0rcv(int phase, FLT *vdata) {
   int stop = 1,i;
   
   for(i=0;i<nsbd;++i)
      stop &= sbdry(i)->comm_nowait(phase);
   for(i=0;i<nvbd;++i)
      stop &= vbdry(i)->comm_nowait(phase);
      
   for(i=0;i<nsbd;++i) 
      hp_sbdry(i)->vmatchsolution_rcv(phase,vdata);
   for(i=0;i<nvbd;++i)
      hp_vbdry(i)->vmatchsolution_rcv(phase,vdata);
      
   return(stop);
}

void tri_hp::sc0load(int phase, FLT *sdata, int bgnmode, int endmode, int modestride) {
   int i;
      
   /* SEND COMMUNICATIONS TO ADJACENT MESHES */\
   for(i=0;i<nsbd;++i) 
      hp_sbdry(i)->smatchsolution_snd(phase,sdata,bgnmode,endmode,modestride);
   
   for(i=0;i<nsbd;++i)
      sbdry(i)->comm_prepare(phase);
   
   return;
}

int tri_hp::sc0wait_rcv(int phase,FLT *sdata, int bgnmode, int endmode, int modestride) {
   int stop = 1;
   int i;
   
   for(i=0;i<nsbd;++i)
      stop &= sbdry(i)->comm_wait(phase);

   for(i=0;i<nsbd;++i) 
      hp_sbdry(i)->smatchsolution_rcv(phase,sdata,bgnmode,endmode,modestride);
      
   return(stop);
}

int tri_hp::sc0rcv(int phase, FLT *sdata, int bgnmode, int endmode, int modestride) {
   int stop = 1,i;
   
   for(i=0;i<nsbd;++i)
      stop &= sbdry(i)->comm_nowait(phase);
      
   for(i=0;i<nsbd;++i) 
      hp_sbdry(i)->smatchsolution_rcv(phase,sdata,bgnmode,endmode,modestride);
      
   return(stop);
}
