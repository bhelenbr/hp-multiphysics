/*
 *  hp_boundary.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 9/3/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */
 
#include <boundary.h>

class hp_vrtx_bdry {
   protected:
      std::string mytype;
      vrtx_bdry &base;
      static FLT dummy;

   public:
      hp_vrtx_bdry(vrtx_bdry &bin) : base(bin) {mytype = "plain";}
      virtual void init(std::map <std::string,std::string>& input, std::string prefix) {}
      virtual void create(const hp_vrtx_bdry &tgt);
      virtual void copy_data(const hp_vrtx_bdry *tgt);
      virtual ~hp_vrtx_bdry() {};
            
      /* BOUNDARY CONDITION FUNCTIONS */
      virtual void vdirichlet() {}
      virtual void vmatchsolution_snd(int phase, FLT *vdata) = 0;
      virtual void vmatchsolution_rcv(int phase, FLT *vdata) = 0;
            
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
      side_bdry &base;
      static FLT dummy;

   public:
      hp_side_bdry(side_bdry &bin) : base(bin) {mytype = "plain";}
      virtual void init(std::map <std::string,std::string>& input, std::string prefix) {}
      virtual void create(const hp_side_bdry &tgt);
      virtual void copy_data(const hp_side_bdry *tgt);
      virtual ~hp_side_bdry() {};
      
      /* CURVATURE FUNCTIONS */
      virtual bool is_curved() {return(false);}
      virtual FLT crds(int ind, int mode, int dir) {return(0.0);}
      virtual FLT crdsbd(int ind, int mode, int dir, int tlvl) {return(0.0);}
      virtual FLT& crdsin(int ind, int mode, int dir) {return(dummy);}
      virtual FLT& crdsbdin(int ind, int mode, int dir, int tlvl) {return(dummy);}
      virtual void curvinit(int tlvl) {}
      
      /* BOUNDARY CONDITION FUNCTIONS */
      virtual void addbflux(bool mgrid) {}
      virtual void vdirichlet() {}
      virtual void sdirichlet(int mode) {}
      virtual void vmatchsolution_snd(int phase, FLT *vdata) = 0;
      virtual void vmatchsolution_rcv(int phase, FLT *vdata)= 0;
      virtual void smatchsolution_snd(int phase, FLT *sdata, int bgnmode, int endmode, int modestride) = 0;
      virtual void smatchsolution_rcv(int phase, FLT *sdata, int bgnmode, int endmode, int modestride) = 0;
            
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
      virtual void findbdrypt(TinyVector<FLT,2> xp,int &bel,FLT &psi) {}
      
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

/* TEMPLATE CLASS TO MAKE BOUNDARY CONDITIONS FOR AN <HPBASE> CLASS */
template<class HPBASE> class hp_vgeneric : public hp_vrtx_bdry { 
   HPBASE &x;
   
   hp_vgeneric(HPBASE& hpin, vrtx_bdry &bin) : hp_vrtx_bdry(bin), x(hpin) {mytype = "generic";}
  
   void vmatchsolution_snd(int phase, Array<FLT,2> &vdata) {base.vloadbuff(&vdata(0,0),0,x.nv-1,x.nv);}
   void vmatchsolution_rcv(int phase, Array<FLT,2> &vdata) {base.vfinalrcv(phase,&vdata(0,0),0,0,x.nv-1,x.nv);}
};


/* TEMPLATE CLASS TO MAKE BOUNDARY CONDITIONS FOR AN <HPBASE> CLASS */
template<class HPBASE> class hp_sgeneric : public hp_side_bdry { 
   HPBASE &x;
   
   hp_sgeneric(HPBASE& hpin, side_bdry &bin) : hp_side_bdry(bin), x(hpin) {mytype = "generic";}
   void findbdrypt(TinyVector<FLT,2> xp,int &bel,FLT &psi) {
      base.findbdrypt(xp,bel,psi);
      basis::tri(x.log2p).ptvalues1d(psi);
   }
   void vmatchsolution_snd(int phase, FLT *vdata) {base.vloadbuff(vdata,0,x.nv-1,x.nv);}
   void vmatchsolution_rcv(int phase, FLT *vdata) {base.vfinalrcv(phase,vdata,0,0,x.nv-1,x.nv);}
   void smatchsolution_snd(int phase, FLT *sdata, int bgnmode, int endmode, int modestride) {
      sloadbuff(sdata,bgnmode*x.nv,(endmode+1)*x.nv-1,x.nv*modestride);
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
};

template<class HPBASE> class hp_curved : public hp_sgeneric<HPBASE> {
   Array<TinyVector<FLT,HPBASE::ND>,2> crv;
   TinyVector<Array<TinyVector<FLT,HPBASE::ND>,2>,sim::nhist+1> crvbd;
   
   hp_curved(HPBASE& hpin, side_bdry &bin) : hp_sgeneric<HPBASE>(hpin,bin) { 
      hp_sgeneric<HPBASE>::mytype = "curved";
   }
   
   void init(std::map <std::string,std::string>& input, std::string prefix) {
      int i,p,coarse,adapt_storage;
      std::string keyword;
      std::istringstream data;
      std::string filename;
   
      keyword = prefix + ".coarse";
      data.str(input[keyword]);
      if (!(data >> coarse)) {
         coarse = 0;
      }
      data.clear();
      
      if (!coarse) {
         crv.resize(hp_sgeneric<HPBASE>::base.maxel,HPBASE::sm0);
         for(i=1;i<sim::nhist+1;++i)
            crvbd(i).resize(hp_sgeneric<HPBASE>::base.maxel,HPBASE::sm0);
         crvbd(0).reference(crv);
      }
      return;
   }
   
   void output(std::ostream& fout, tri_hp::filetype typ,int tlvl = 0) {
      int j,m,n;
      
      hp_sgeneric<HPBASE>::output(fout,typ,tlvl);
      
      switch(typ) {
         case(text):
            fout << "p0: " << HPBASE::p0;
            for(j=0;j<hp_sgeneric<HPBASE>::base.nel;++j) {
               for(m=0;m<HPBASE::sm0;++m) {
                  for(n=0;n<HPBASE::ND;++n)
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
      hp_sgeneric<HPBASE>::input(fin,typ,tlvl);
      
      switch(typ) {
         case(text):
            fin.ignore(80,':');
            fin >>  pmin;
            for(j=0;j<hp_sgeneric<HPBASE>::base.nel;++j) {
               for(m=0;m<pmin -1;++m) {
                  for(n=0;n<HPBASE::ND;++n)
                     fin >> crvbd(tlvl)(j,m)(n);
               }
               for(m=pmin-1;m<HPBASE::sm0;++m) {
                  for(n=0;n<HPBASE::ND;++n)
                     crvbd(tlvl)(j,m)(n) = 0.0;
               }
            }
            
         default:
            break;
      }
      return;
   }
   
   void curv_init(int tlvl) {
      int i,j,m,n,v0,v1,sind,info;
      TinyVector<FLT,2> pt;
      char uplo = 'U';
      
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
      for(j=0;j<hp_sgeneric<HPBASE>::base.nel;++j) {
         sind = hp_sgeneric<HPBASE>::base.el(j);
         v0 = HPBASE::sd(sind).vrtx(0);
         hp_sgeneric<HPBASE>::base.mvpttobdry(j,-1.0, HPBASE::vrtx(v0));
      }
      v0 = HPBASE::sd(sind).vrtx(1);
      hp_sgeneric<HPBASE>::base.mvpttobdry(hp_sgeneric<HPBASE>::base.nel-1,1.0, HPBASE::vrtx(v0));

      if (basis::tri(HPBASE::log2p).p == 1) return;
         
      /*****************************/
      /* SET UP HIGHER ORDER MODES */
      /*****************************/
      for(j=0;j<hp_sgeneric<HPBASE>::base.nel;++j) {
         sind = hp_sgeneric<HPBASE>::base.el(j);

         v0 = HPBASE::sd(sind).vrtx(0);
         v1 = HPBASE::sd(sind).vrtx(1);
         
         for(n=0;n<HPBASE::ND;++n) 
            basis::tri(HPBASE::log2p).proj1d(HPBASE::vrtx(v0)(n),HPBASE::vrtx(v1)(n),&HPBASE::crd(n)(0,0));
      
         for(i=0;i<basis::tri(HPBASE::log2p).gpx;++i) {
            pt(0) = HPBASE::crd(0)(0,i);
            pt(1) = HPBASE::crd(1)(0,i);
            hp_sgeneric<HPBASE>::base.mvpttobdry(j,(1.-2.*basis::tri(HPBASE::log2p).gx(1,i)),pt);
            HPBASE::crd(0)(0,i) -= pt(0);
            HPBASE::crd(1)(0,i) -= pt(1);
         }
         
         for(n=0;n<HPBASE::ND;++n) {
            basis::tri(HPBASE::log2p).intgrt1d(&HPBASE::cf(n,0),&HPBASE::crd(n)(0,0));
            PBTRS(uplo,basis::tri(HPBASE::log2p).sm,basis::tri(HPBASE::log2p).sbwth,1,&basis::tri(HPBASE::log2p).sdiag1d(0,0),basis::tri(HPBASE::log2p).sbwth+1,&HPBASE::cf(n,2),basis::tri(HPBASE::log2p).sm,info);
         
            for(m=0;m<basis::tri(HPBASE::log2p).sm;++m)
               crvbd(tlvl)(j,m)(n) = -HPBASE::cf(n,m+2);
         }
      }
      return;
   }
   
   void findbdrypt(TinyVector<FLT,2> xp,int &bel,FLT &psi) {
      int sind,v0,v1,iter;
      FLT dx,dy,ol,roundoff,dpsi;
      TinyVector<FLT,2> pt;
      
      hp_sgeneric<HPBASE>::base.findbdrypt(xp,bel,psi);
      sind = hp_sgeneric<HPBASE>::base.el(bel);
      v0 = hp_sgeneric<HPBASE>::x.sd(sind).vrtx(0);
      v1 = hp_sgeneric<HPBASE>::x.sd(sind).vrtx(1);
      dx = hp_sgeneric<HPBASE>::x.vrtx(v1)(0) - hp_sgeneric<HPBASE>::x.vrtx(v0)(0);
      dy = hp_sgeneric<HPBASE>::x.vrtx(v1)(1) - hp_sgeneric<HPBASE>::x.vrtx(v0)(1);
      ol = 2./(dx*dx +dy*dy);
      dx *= ol;
      dy *= ol;
      
      /* FIND PSI SUCH THAT TANGENTIAL POSITION ALONG LINEAR SIDE STAYS THE SAME */
      /* THIS WAY, MULTIPLE CALLS WILL NOT GIVE DIFFERENT RESULTS */ 
      hp_sgeneric<HPBASE>::x.crdtocht1d(sind);

      iter = 0;
      roundoff = 10.0*EPSILON*(1.0 +(fabs(xp(0)*dx) +fabs(xp(1)*dy)));
      do {
         basis::tri(HPBASE::log2p).ptprobe1d(HPBASE::ND,pt.data(),psi,&hp_sgeneric<HPBASE>::x.cht(0,0),MXTM);
         dpsi = (pt(0) -xp(0))*dx +(pt(1) -xp(1))*dy;
         psi += dpsi;
         if (iter++ > 100) {
            *sim::log << "#Warning: max iterations for curved side in bdry_locate type: " << hp_sgeneric<HPBASE>::base.idnum << " loc: " << xp << std::endl;
            break;
         }  
      } while (fabs(dpsi) > roundoff);
   }

};

