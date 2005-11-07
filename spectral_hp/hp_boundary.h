/*
 *  hp_boundary.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 9/3/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */
 
#include <boundary.h>
#include <hpbasis.h>
#include <myblas.h>

class hp_vrtx_bdry {
   protected:
      std::string mytype;
      tri_hp& x;
      vrtx_bdry& base;
      static FLT dummy;

   public:
      hp_vrtx_bdry(tri_hp& xin, vrtx_bdry &bin) : x(xin), base(bin) {mytype = "plain";}
      virtual void init(std::map <std::string,std::string>& input, std::string prefix) {}
      virtual void create(const hp_vrtx_bdry &tgt);
      virtual void copy_data(const hp_vrtx_bdry *tgt);
      virtual ~hp_vrtx_bdry() {};
            
      /* BOUNDARY CONDITION FUNCTIONS */
      virtual void vdirichlet() {}
      virtual void vmatchsolution_snd(int phase, FLT *vdata) {base.vloadbuff(vdata,0,x.NV-1,x.NV);}
      virtual void vmatchsolution_rcv(int phase, FLT *vdata) {base.vfinalrcv(phase,vdata,0,x.NV-1,x.NV);}
            
      /* FOR COUPLED DYNAMIC BOUNDARIES */
      virtual block::ctrl setup_preconditioner(int excpt) {return(block::stop);}
      virtual block::ctrl tadvance(int excpt) {return(block::stop);}
      virtual block::ctrl rsdl(int excpt) {return(block::stop);}
      virtual block::ctrl update(int excpt) {return(block::stop);}
      virtual block::ctrl minvrt(int excpt) {return(block::stop);}
      virtual void mg_getfres(int excpt, Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, hp_side_bdry *fbdry) {} 
      virtual void mg_getcchng(int excpt, Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, hp_side_bdry *fbdry) {}   
      
      /* input output functions */
      virtual void output(std::ostream& fout, tri_hp::filetype typ,int tlvl = 0) {
         switch(typ) {
            case(text):
               fout << base.idprefix << " " << mytype << std::endl;
               break;
            default:
               break;
         }
         return;
      }
         
      virtual void input(ifstream& fin,tri_hp::filetype typ,int tlvl = 0) {
         std::string keyword,value;
         switch(typ) {
            case(text):
               fin >> keyword >> value;
               if (keyword != base.idprefix) {
                  *sim::log << "types don't match " << base.idprefix << ' ' << keyword << std::endl;
                  exit(1);
               }
               if (value != mytype) {
                  *sim::log << "types don't match " << value << ' ' << mytype << std::endl;
                  exit(1);
               }
               break;
            default:
               break;
         }

         return;
      }
};

class hp_side_bdry {
   protected:
      std::string mytype;
      tri_hp& x;
      side_bdry &base;
      Array<TinyMatrix<FLT,mesh::ND,MXGP>,2> dxdt;
      static FLT dummy;

   public:
      hp_side_bdry(tri_hp& xin, side_bdry &bin) : x(xin), base(bin) {mytype = "plain";}
      virtual void init(std::map <std::string,std::string>& input, std::string prefix) {
         dxdt.resize(x.log2pmax,base.maxel);
      }
      virtual void create(const hp_side_bdry &tgt);
      virtual void copy_data(const hp_side_bdry *tgt);
      virtual ~hp_side_bdry() {};
      
      /* CURVATURE FUNCTIONS */
      virtual bool is_curved() {return(false);}
      virtual FLT& crds(int ind, int mode, int dir) {return(dummy=0);}
      virtual FLT& crdsbd(int ind, int mode, int dir, int tlvl) {return(dummy=0);}
      virtual void curvinit(int tlvl) {}
      
      /* BOUNDARY CONDITION FUNCTIONS */
      virtual void addbflux(bool mgrid) {}
      virtual void vdirichlet() {}
      virtual void sdirichlet(int mode) {}
      void vmatchsolution_snd(int phase, FLT *vdata) {base.vloadbuff(vdata,0,x.NV-1,x.NV);}
      void vmatchsolution_rcv(int phase, FLT *vdata) {base.vfinalrcv(phase,vdata,0,x.NV-1,x.NV);}
      void smatchsolution_snd(int phase, FLT *sdata, int bgnmode, int endmode, int modestride) {
         base.sloadbuff(sdata,bgnmode*x.NV,(endmode+1)*x.NV-1,x.NV*modestride);
         return;
      }
      void smatchsolution_rcv(int phi, FLT *sdata, int bgn, int end, int stride) {
         
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
            std::cout << "finalrcv"  << idnum << " " << is_frst() << std::endl;
   #endif
            count = 0;
            for(j=0;j<base.nel;++j) {
               sind = base.el(j);
               offset = sind*stride*x.NV;
               for (k=bgn;k<=end;++k) {
                  for(n=0;n<x.NV;++n) {
                     sdata[offset] = base.fsndbuf(count++)*mtchinv;
   #ifdef MPDEBUG
                     std::cout << "\t" << sdata[offset+k] << std::endl;
   #endif
                     ++offset;
                  }
               }

            }
         }
         return;
      }

   
            
      /* FOR COUPLED DYNAMIC BOUNDARIES */
      virtual block::ctrl setup_preconditioner(int excpt) {return(block::stop);}
      virtual block::ctrl tadvance(int excpt) {return(block::stop);}
      virtual block::ctrl rsdl(int excpt) {return(block::stop);}
      virtual block::ctrl update(int excpt) {return(block::stop);}
      virtual block::ctrl minvrt(int excpt) {return(block::stop);}
      virtual void mg_getfres(int excpt, Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, hp_side_bdry *fbdry) {} 
      virtual void mg_getcchng(int excpt, Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, hp_side_bdry *fbdry) {}   
      
      /* ADAPTATION FUNCTIONS */
      virtual void updatevdata_bdry(int bel,int endpt,hp_side_bdry *bin) {}
      virtual void movevdata_bdry(int bel,int endpt,hp_side_bdry *bin) {}
      virtual void updatesdata_bdry(int bel,hp_side_bdry *bin) {}
      virtual void movesdata_bdry(int bel,hp_side_bdry *bin) {}
      
      /* SEARCH FUNCTIONS */
      virtual void findbdrypt(TinyVector<FLT,2> xp,int &bel,FLT &psi) {
         base.findbdrypt(xp,bel,psi);
         basis::tri(x.log2p).ptvalues1d(psi);
      }
      
      /* input output functions */
      virtual void output(std::ostream& fout, tri_hp::filetype typ,int tlvl = 0) {
         switch(typ) {
            case(text):
               fout << base.idprefix << " " << mytype << std::endl;
               break;
            default:
               break;
         }
         return;
      }
         
      virtual void input(ifstream& fin,tri_hp::filetype typ,int tlvl = 0) {
         std::string keyword,value;
         switch(typ) {
            case(text):
               fin >> keyword >> value;
               if (keyword != base.idprefix) {
                  *sim::log << "types don't match " << base.idprefix << ' ' << keyword << std::endl;
                  exit(1);
               }
               if (value != mytype) {
                  *sim::log << "types don't match " << value << ' ' << mytype << std::endl;
                  exit(1);
               }
               break;
            default:
               break;
         }

         return;
      }
      
      /* ENPOINT COMMUNICATION ROUTINES */
      virtual void findendpts(class mesh& x) {}
      virtual void msgload_vrtx(int endpt,FLT *base,int bgn,int end,int stride) {}
      virtual void msgpass_vrtx(int endpt) {}
      virtual void msgwait_rcv_vrtx(int endpt,FLT *base,int bgn,int end,int stride) {}
};


class hp_curved : public hp_side_bdry {
   bool coupled;
   Array<TinyVector<FLT,mesh::ND>,2> crv;
   TinyVector<Array<TinyVector<FLT,mesh::ND>,2>,sim::nhist+1> crvbd;
   
   hp_curved(tri_hp& hpin, side_bdry &bin) : hp_side_bdry(hpin,bin) { 
      mytype = "curved";
   }
   
   void init(std::map <std::string,std::string>& input, std::string prefix) {
      int i,coarse;
      std::string keyword;
      std::istringstream data;
      std::string filename;
      
      hp_side_bdry::init(input,prefix);
   
      keyword = prefix + ".coarse";
      data.str(input[keyword]);
      if (!(data >> coarse)) {
         coarse = 0;
      }
      data.clear();
      
      keyword = prefix + ".coupled";
      data.str(input[keyword]);
      if (!(data >> coupled)) {
         coupled = false;
      }
      data.clear();
      
      if (!coarse) {
         crv.resize(base.maxel,x.sm0);
         for(i=1;i<sim::nhist+1;++i)
            crvbd(i).resize(base.maxel,x.sm0);
         crvbd(0).reference(crv);
      }
         
      return;
   }
   
   void output(std::ostream& fout, tri_hp::filetype typ,int tlvl = 0) {
      int j,m,n;
      
      hp_side_bdry::output(fout,typ,tlvl);
      
      switch(typ) {
         case(text):
            fout << "p0: " << x.p0;
            for(j=0;j<base.nel;++j) {
               for(m=0;m<x.sm0;++m) {
                  for(n=0;n<mesh::ND;++n)
                     fout << crvbd(tlvl)(j,m)(n) << ' ';
                  fout << std::endl;
               }
            }
            
         default:
            break;
      }
      return;
   }
         
   virtual void input(ifstream& fin,tri_hp::filetype typ,int tlvl = 0) {
      int j,m,n,pmin;
      hp_side_bdry::input(fin,typ,tlvl);
      
      switch(typ) {
         case(text):
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
            
         default:
            break;
      }
      return;
   }
   
   void curv_init(int tlvl = 0) {
      int i,j,m,n,v0,v1,sind,info;
      TinyVector<FLT,2> pt;
      char uplo[] = "U";

      
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
         base.mvpttobdry(j,-1.0, x.vrtx(v0));
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
            base.mvpttobdry(j,(1.-2.*basis::tri(x.log2p).gx(1,i)),pt);
            x.crd(0)(0,i) -= pt(0);
            x.crd(1)(0,i) -= pt(1);
         }
         
         for(n=0;n<mesh::ND;++n) {
            basis::tri(x.log2p).intgrt1d(&x.cf(n,0),&x.crd(n)(0,0));
            DPBTRS(uplo,basis::tri(x.log2p).sm,basis::tri(x.log2p).sbwth,1,&basis::tri(x.log2p).sdiag1d(0,0),basis::tri(x.log2p).sbwth+1,&x.cf(n,2),basis::tri(x.log2p).sm,info);
         
            for(m=0;m<basis::tri(x.log2p).sm;++m)
               crvbd(tlvl)(j,m)(n) = -x.cf(n,m+2);
         }
      }
      return;
   }
   
   void findbdrypt(TinyVector<FLT,2> xp,int &bel,FLT &psi) {
      int sind,v0,v1,iter;
      FLT dx,dy,ol,roundoff,dpsi;
      TinyVector<FLT,2> pt;
      
      base.findbdrypt(xp,bel,psi);
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
         psi += dpsi;
         if (iter++ > 100) {
            *sim::log << "#Warning: max iterations for curved side in bdry_locate type: " << base.idnum << " loc: " << xp << std::endl;
            break;
         }  
      } while (fabs(dpsi) > roundoff);
   }
   
   block::ctrl tadvance(int excpt) {
      int stage = sim::dirkstage +sim::esdirk;
   
      switch(excpt) {
         case(0): {
            if (stage) {
               /* BACK CALCULATE K TERM */
               for(int j=0;j<base.nel;++j) {
                  for(int m=0;m<basis::tri(x.log2p).sm;++m) {
                     for(int n=0;n<mesh::ND;++n)
                        crvbd(stage)(j,m)(n) = (crvbd(0)(j,m)(n)-crvbd(1)(j,m)(n))*sim::adirk[sim::dirkstage][sim::dirkstage];
                  }
               }
            }
            
            if (sim::dirkstage == 0) {
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

};

