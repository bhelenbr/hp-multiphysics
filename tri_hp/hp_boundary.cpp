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

hp_vrtx_bdry* tri_hp::getnewvrtxobject(int bnum, input_map &bdrydata) {
   hp_vrtx_bdry *temp = new hp_vrtx_bdry(*this,*vbdry(bnum));   
   return(temp);
}

hp_side_bdry* tri_hp::getnewsideobject(int bnum, input_map &bdrydata) {
   hp_side_bdry *temp = new hp_side_bdry(*this,*sbdry(bnum));
   return(temp);
}

FLT hp_vrtx_bdry::dummy;

void hp_side_bdry::smatchsolution_rcv(FLT *sdata, int bgn, int end, int stride) {
   /* CAN'T USE sfinalrcv BECAUSE OF CHANGING SIGNS */
   int j,k,m,n,count,countdn,countup,offset,sind,sign;
   FLT mtchinv;
   
   if (!base.is_comm()) return;

   int matches = 1;
   int bgnsign = (bgn % 2 ? -1 : 1);
   
   /* ASSUMES REVERSE ORDERING OF SIDES */
   for(m=0;m<base.matches();++m) {   
         
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

void hp_side_bdry::init(input_map& inmap,void* &gbl_in) {
   int i;
   bool coarse;
   std::string keyword;
   std::istringstream data;
   std::string filename;
   
   dxdt.resize(x.log2pmax+1,base.maxel);
   
   keyword = base.idprefix + ".curved";
   inmap.getwdefault(keyword,curved,false);

   keyword = x.idprefix + ".coarse";
   inmap.getwdefault(keyword,coarse,false);
   
   keyword = base.idprefix + ".coupled";
   inmap.getwdefault(keyword,coupled,false);
   
   if (curved && !coarse) {
      crv.resize(base.maxel,x.sm0);
      for(i=1;i<sim::nhist+1;++i)
         crvbd(i).resize(base.maxel,x.sm0);
      crvbd(0).reference(crv);
   }
   
   base.resize_buffers(base.maxel*(x.sm0+2)*x.NV);

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
      case(tri_hp::text):
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

void hp_side_bdry::mvpttobdry(int bel,FLT psi,TinyVector<FLT,2> &xp) {
   int sind;
   sind = base.el(bel);
   basis::tri(x.log2p).ptvalues1d(psi);
   x.crdtocht1d(sind);
   basis::tri(x.log2p).ptprobe1d(x.ND,xp.data(),psi,&x.cht(0,0),MXTM);
   return;
}


block::ctrl hp_side_bdry::tadvance(bool coarse, block::ctrl ctrl_message) {
   int stage = sim::substep +sim::esdirk;  
      
   if (ctrl_message == block::begin) {
      if (x.p0 > 1 && curved) {
         if (stage) {
            /* BACK CALCULATE K TERM */
            for(int j=0;j<base.nel;++j) {
               for(int m=0;m<basis::tri(x.log2p).sm;++m) {
                  for(int n=0;n<mesh::ND;++n)
                     crvbd(stage+1)(j,m)(n) = (crvbd(0)(j,m)(n)-crvbd(1)(j,m)(n))*sim::adirk[stage-1][stage-1];
               }
            }
         }
         
         if (sim::substep == 0) {
            /* STORE TILDE W */
            for(int j=0;j<base.nel;++j) {
               for(int m=0;m<basis::tri(x.log2p).sm;++m) {
                  for(int n=0;n<mesh::ND;++n)
                     crvbd(1)(j,m)(n) = crv(j,m)(n);
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
      }
      
      calculate_unsteady_sources(coarse);
      
      /* EXTRAPOLATE SOLUTION OR MOVE TO NEXT TIME POSITION */
      if (!coupled && curved) curv_init();
   }
   return(block::stop);
}

void hp_side_bdry::calculate_unsteady_sources(bool coarse) {
   int i,j,n,sind;
   
   for(i=0;i<=x.log2pmax;++i) {
      for(j=0;j<base.nel;++j) {
         sind = base.el(j);
         x.crdtocht1d(sind,1);
         for(n=0;n<mesh::ND;++n)
            basis::tri(i).proj1d(&x.cht(n,0),&dxdt(i,j)(n,0));
      }
   }
   
   return;
}




void tri_hp::vc0load(int phase, FLT *vdata) {
   int i;
         
   /* SEND COMMUNICATIONS TO ADJACENT MESHES */\
   for(i=0;i<nsbd;++i) 
      hp_sbdry(i)->vmatchsolution_snd(phase,vdata);
   for(i=0;i<nvbd;++i)
      hp_vbdry(i)->vmatchsolution_snd(phase,vdata);
   
   for(i=0;i<nsbd;++i)
      sbdry(i)->comm_prepare(boundary::all_phased,phase,boundary::symmetric);
   for(i=0;i<nvbd;++i)
      vbdry(i)->comm_prepare(boundary::all_phased,phase,boundary::symmetric);
   
   return;
}
int tri_hp::vc0wait_rcv(int phase, FLT *vdata) {
   int stop = 1;
   int i;
   
   for(i=0;i<nsbd;++i)
      stop &= sbdry(i)->comm_wait(boundary::all_phased,phase,boundary::symmetric);
   for(i=0;i<nvbd;++i)
      stop &= vbdry(i)->comm_wait(boundary::all_phased,phase,boundary::symmetric);
      

   for(i=0;i<nsbd;++i) 
      hp_sbdry(i)->vmatchsolution_rcv(phase,vdata);
   for(i=0;i<nvbd;++i)
      hp_vbdry(i)->vmatchsolution_rcv(phase,vdata);
      
   return(stop);
}

int tri_hp::vc0rcv(int phase, FLT *vdata) {
   int stop = 1,i;
   
   for(i=0;i<nsbd;++i)
      stop &= sbdry(i)->comm_nowait(boundary::all_phased,phase,boundary::symmetric);
   for(i=0;i<nvbd;++i)
      stop &= vbdry(i)->comm_nowait(boundary::all_phased,phase,boundary::symmetric);
      
   for(i=0;i<nsbd;++i) 
      hp_sbdry(i)->vmatchsolution_rcv(phase,vdata);
   for(i=0;i<nvbd;++i)
      hp_vbdry(i)->vmatchsolution_rcv(phase,vdata);
      
   return(stop);
}

void tri_hp::sc0load(FLT *sdata, int bgnmode, int endmode, int modestride) {
   int i;
      
   /* SEND COMMUNICATIONS TO ADJACENT MESHES */\
   for(i=0;i<nsbd;++i) 
      hp_sbdry(i)->smatchsolution_snd(sdata,bgnmode,endmode,modestride);
   
   for(i=0;i<nsbd;++i)
      sbdry(i)->comm_prepare(boundary::all,0,boundary::symmetric);
   
   return;
}

int tri_hp::sc0wait_rcv(FLT *sdata, int bgnmode, int endmode, int modestride) {
   int stop = 1;
   int i;
   
   for(i=0;i<nsbd;++i)
      stop &= sbdry(i)->comm_wait(boundary::all,0,boundary::symmetric);

   for(i=0;i<nsbd;++i) 
      hp_sbdry(i)->smatchsolution_rcv(sdata,bgnmode,endmode,modestride);
      
   return(stop);
}

int tri_hp::sc0rcv(FLT *sdata, int bgnmode, int endmode, int modestride) {
   int stop = 1,i;
   
   for(i=0;i<nsbd;++i)
      stop &= sbdry(i)->comm_nowait(boundary::all,0,boundary::symmetric);
      
   for(i=0;i<nsbd;++i) 
      hp_sbdry(i)->smatchsolution_rcv(sdata,bgnmode,endmode,modestride);
      
   return(stop);
}
   

block::ctrl tri_hp::matchboundaries(block::ctrl ctrl_message) {
   int i, m, n, msgn, bnum, count;
   block::ctrl state;
   
   if (ctrl_message == block::begin) excpt = 0;
   
   switch(excpt) {
      case 0: {
         /* Match boundary vertices */
         state = mesh::matchboundaries(ctrl_message);
         if (state != block::stop) return(state);
         else {
            ++excpt;
            mp_phase = -1;
            ctrl_message = block::stay;
         }
      }
      case 1: {
         if (ctrl_message == block::stay) {
            
            if (!sm0) {
               excpt += 2;
            }
            else {
               /* Match curved sides */
               for(bnum=0;bnum<nsbd;++bnum) {
                  if (sbdry(bnum)->is_comm() && hp_sbdry(bnum)->is_curved()) {            
                     count = 0;
                     for(i=0;i<sbdry(bnum)->nel;++i) {
                        for(m=0;m<basis::tri(log2p).sm;++m) {
                           for(n=0;n<ND;++n)
                              sbdry(bnum)->fsndbuf(count++) = hp_sbdry(bnum)->crds(i,m,n);
                        }
                     }
                  }
                  sbdry(bnum)->sndsize() = count;
                  sbdry(bnum)->sndtype() = boundary::flt_msg;
                  sbdry(bnum)->comm_prepare(boundary::all,0,boundary::master_slave);
               }
            }
         }
         ++excpt;
         return(block::advance);
      }
      
      case 2: {
         for(bnum=0;bnum<nsbd;++bnum) {
            if (sbdry(bnum)->is_comm() && hp_sbdry(bnum)->is_curved()) {            
               sbdry(bnum)->comm_exchange(boundary::all,0,boundary::master_slave);
            }
         }
         ++excpt;
         return(block::advance);
      }
      case 3: {
         for(bnum=0;bnum<nsbd;++bnum) {
            if (!sbdry(bnum)->is_frst() && hp_sbdry(bnum)->is_curved()) {            
               sbdry(bnum)->comm_wait(boundary::all,0,boundary::master_slave);
               
               count = 0;
               for(i=sbdry(bnum)->nel-1;i>=0;--i) {
                  msgn = 1;
                  for(m=0;m<basis::tri(log2p).sm;++m) {
                     for(n=0;n<ND;++n)
                        hp_sbdry(bnum)->crds(i,m,n) = msgn*sbdry(bnum)->frcvbuf(0,count++);
                     msgn *= -1;
                  }
               }
            }
         }
      }
      
      case 4: {
         mp_phase = -1;
         ++excpt;
         ctrl_message = block::stay;
      }
            
      case 5: {
         ++mp_phase;
         excpt += ctrl_message;
         switch(mp_phase%3) {
            case(0):
               vc0load(mp_phase/3,ug.v.data());
               return(block::stay);
            case(1):
               vmsgpass(boundary::all_phased,mp_phase/3,boundary::symmetric);
               return(block::stay);
            case(2):
               return(static_cast<block::ctrl>(vc0wait_rcv(mp_phase/3,ug.v.data())));
         }
      }
      case 6: {
         mp_phase = -1;
         ++excpt;
         ctrl_message = block::stay;
      }
      case 7: {
         if (ctrl_message == block::stay) {
               
            if (!sm0) return(block::advance);
            
            ++mp_phase;
            switch(mp_phase%3) {
               case(0):
                  sc0load(ug.s.data(),0,sm0-1,ug.s.extent(secondDim));
                  return(block::stay);
               case(1):
                  smsgpass(boundary::all,0,boundary::symmetric);
                  return(block::stay);
               case(2):
                  sc0wait_rcv(ug.s.data(),0,sm0-1,ug.s.extent(secondDim));
                  return(block::advance);
            }
         }
         else
            ++excpt;
      }
   }
   
   return(block::stop);
}