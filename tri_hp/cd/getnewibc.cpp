      
/*
 *  getnewibc_ins.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 2/6/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_cd.h"

namespace ibc_cd {

#ifdef CONVDIFF

#if (CASE == 0)
FLT forcing(FLT x,FLT y) { return(-cos(2.*M_PI*x));}
#else
FLT forcing(FLT x,FLT y) {return(0.0);}
#endif
FLT blayer = 0.0;
FLT axext = AMP*cos(M_PI*LAM/180.0), ayext = AMP*sin(M_PI*LAM/180.0);
FLT nuext = (1.0 - AMP);

FLT f1(int n, FLT x, FLT y) {
   FLT nux = nuext*4.*M_PI*M_PI;
   FLT axx = axext*2.*M_PI;
   FLT xx = x*2.*M_PI;
   FLT yx = y*2.*M_PI;
   double eps = 0.05;
   
   switch(n) {
      case((0+CASE)%3):
         return(axx*sin(xx)/(axx*axx +nux*nux)
           +nux*cos(xx)/(axx*axx +nux*nux)
           +(nux > 0.0 ? blayer*exp(axx/nux*xx) : 0.0)
           +startup*((sin(xx) +sin(100*xx))*(sin(yx)+sin(100*yx))));
      case((1+CASE)%3):
         return(startup*((sin(xx) +sin(100*xx))*(sin(yx)+sin(100*yx))));
      case((2+CASE)%3):
         if (x > 2.*eps)
            return(0.0);
         return(1+cos(M_PI*(x-eps)/eps));
   }
   return(0.0);
}

FLT df1d(int n, FLT x, FLT y) {
   FLT nux = nuext*4.*M_PI*M_PI;
   FLT axx = axext*2.*M_PI;
   FLT xx = x*2.*M_PI;
   switch(n) {
      case((0+CASE)%3):
         return(axx*cos(xx)/(axx*axx +nux*nux) 
           -nuext*sin(x)/(axx*axx +nux*nux)
           +(nux > 0.0 ? blayer*axx/nux*exp(axx/nux*xx) : 0.0));
      case(!CASE):
         return(0.0);
   }
   return(0.0);
}

#endif



   class ibc_type {
      public:
         const static int ntypes = 0;
         enum ids {};
         const static char names[ntypes][40];
         static int getid(const char *nin) {
            int i;
            for(i=0;i<ntypes;++i) 
               if (!strcmp(nin,names[i])) return(i);
            return(-1);
      }
   };
   const char ibc_type::names[ntypes][40] = {};

   class zero_src : public init_bdry_cndtn {
      private:
      public:
         FLT f(int n, TinyVector<FLT,mesh::ND> x) {
            return(0.0);
         }
         void input(input_map &inmap,std::string idnty) {}
   };

   
   class power_src : public init_bdry_cndtn {
      private:
         Array<FLT,1> c;
         Array<FLT,2> a;
         Array<FLT,2> spd;
      public:
         FLT f(int n, TinyVector<FLT,mesh::ND> x) {
            return(c(n)*pow(x(0)-spd(n,0)*sim::time,a(n,0))*pow(x(1)-spd(n,1)*sim::time,a(n,1)));
         }
         void input(input_map &inmap,std::string idnty) {
            std::string keyword,val;
            std::istringstream data;
            std::ostringstream nstr;
            int nvar;
            
            keyword = idnty +"_nvariable";

            if (!inmap.get(keyword,nvar))
               inmap.getwdefault("nvariable",nvar,1);
            
            c.resize(nvar);
            a.resize(nvar,2);
            spd.resize(nvar,2);

            for(int n=0;n<nvar;++n) {
               nstr.str("");
               nstr << idnty << "_src_coeff" << n << std::flush;
               if (!inmap.get(nstr.str(),c(n))) {
                  inmap.getwdefault("src_coeff",c(n),0.0);
               }

               nstr.str("");
               nstr << idnty << "_src_powers" << n << std::flush;
               if (!inmap.getline(nstr.str(),val)) {
                  nstr.str("");
                  nstr << "src_powers" << n << std::flush;
                  inmap.getlinewdefault(nstr.str(),val,"0.0 0.0");
               }            
               data.str(val);
               data >> a(n,0) >> a(n,1);  
               data.clear(); 
               
               nstr.str("");
               nstr << idnty << "_src_speeds" << n << std::flush;
               if (!inmap.getline(nstr.str(),val)) {
                  nstr.str("");
                  nstr << "src_speeds" << n << std::flush;
                  inmap.getlinewdefault(nstr.str(),val,"0.0 0.0");
               }            
               data.str(val);
               data >> spd(n,0) >> spd(n,1);  
               data.clear(); 
            }
         }
   };

   class soi_src : public init_bdry_cndtn {

       private:
         double xl,xr,yt,yb,c;

       public:
          FLT f(int n, TinyVector<FLT,mesh::ND> x){

            
            if( x(0) >= xl && x(0) <= xr && x(1) <= yt && x(1) >= yb )		//xdl < xdr < xgl < xj < xgr < xsl < xsr
               return c;
            else
               return 0;

          }


          void input(input_map &blockdata,std::string idnty) {
             std::string keyword,val;
             std::istringstream data;
            
            keyword = idnty+"_xl";
            blockdata.getwdefault(keyword,xl,0.68);
            
            keyword = idnty+"_xr";
            blockdata.getwdefault(keyword,xr,0.7);

            keyword = idnty+"_yt";
            blockdata.getwdefault(keyword,yt,0.0);

            keyword = idnty+"_yb";
            blockdata.getwdefault(keyword,yb,-0.07);
            
            keyword = idnty+"_constant";
            blockdata.getwdefault(keyword,c,0.41);

          }

   };
   
   class src_type {
      public:
         const static int ntypes = 3;
         enum ids {zero,power,soi};
         const static char names[ntypes][40];
         static int getid(const char *nin) {
            int i;
            for(i=0;i<ntypes;++i) 
               if (!strcmp(nin,names[i])) return(i);
            return(-1);
      }
   };
   const char src_type::names[ntypes][40] = {"zero","power","soi"};

   
}

init_bdry_cndtn *tri_hp_cd::getnewibc(input_map& inmap) {
   std::string keyword,ibcname;
   int type;

   /* FIND INITIAL CONDITION TYPE */
   keyword = std::string(idprefix) + "_ibc";
   if (!inmap.get(keyword,ibcname)) {
      if (!inmap.get("ibc",ibcname)) {
         *sim::log << "couldn't find initial condition type" << std::endl;
      }
   }
   type = ibc_cd::ibc_type::getid(ibcname.c_str());

      
   switch(type) {
      default: {
         return(tri_hp::getnewibc(inmap));
      }
   }
}

init_bdry_cndtn *tri_hp_cd::getnewsrc(input_map& inmap) {
   std::string keyword,ibcname;
   int type;

   /* FIND INITIAL CONDITION TYPE */
   keyword = std::string(idprefix) + "_src";
   if (!inmap.get(keyword,ibcname)) {
      if (!inmap.get("src",ibcname)) {
         *sim::log << "couldn't find source type" << std::endl;
      }
   }
   type = ibc_cd::src_type::getid(ibcname.c_str());
   
   

      
   switch(type) {
      case(ibc_cd::src_type::zero):
         return(new ibc_cd::zero_src);
      case(ibc_cd::src_type::power):
         return(new ibc_cd::power_src);
      case(ibc_cd::src_type::soi):
         return(new ibc_cd::soi_src);
      default: {
         return(new ibc_cd::zero_src);
      }
   }
}

