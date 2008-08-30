        
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

    class zero_src : public init_bdry_cndtn {
        private:
        public:
            FLT f(int n, TinyVector<FLT,tri_mesh::ND> x,FLT time) {
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
            FLT f(int n, TinyVector<FLT,tri_mesh::ND> x, FLT time) {
                return(c(n)*pow(x(0)-spd(n,0)*time,a(n,0))*pow(x(1)-spd(n,1)*time,a(n,1)));
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
        tri_hp_cd &x;
        
         public:
            double xl,xr,yt,yb,c;
            FLT f(int n, TinyVector<FLT,tri_mesh::ND> x, FLT time){

                
                if( x(0) >= xl && x(0) <= xr && x(1) <= yt && x(1) >= yb )		//xdl < xdr < xgl < xj < xgr < xsl < xsr
                    return c;
                else
                    return 0;

             }
             
            soi_src(tri_hp_cd& xin) : x(xin) {}

            void input(input_map &blockdata,std::string idnty); 

    };
    
    class source_changer : public tri_hp_helper {
        protected:
            tri_hp_cd &x;
            FLT delta_c;
            int interval;
        public:
            soi_src *src;

        public:
            source_changer(tri_hp_cd& xin) : tri_hp_helper(xin), x(xin) {}
            
            void init(input_map& input, std::string idnty) {
                std::string keyword, val;
                keyword = idnty + "_delta_c";
                if (!input.get(keyword,delta_c)) {
                    input.getwdefault("delta_c",delta_c,0.0);
                }
                input.getwdefault("parameter_interval",interval,1);
            }
            
            tri_hp_helper* create(tri_hp& xin) { return new source_changer(dynamic_cast<tri_hp_cd&>(xin));}            

            void tadvance() {
                if (x.coarse_level) return;
            
                if ( (x.gbl->tstep % interval) +x.gbl->substep == 0) {
                    src->c += delta_c;
                    *x.gbl->log << "new source strength is " << src->c << std::endl;
                }            
                return;
            }
    };
    
    void soi_src::input(input_map &blockdata,std::string idnty) {
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
        
        if (source_changer *sptr = dynamic_cast<source_changer *>(x.helper)) {
            sptr->src = this;
        }
     }
    
    class helper_type {
        public:
            const static int ntypes = 1;
            enum ids {source_changer};
            const static char names[ntypes][40];
            static int getid(const char *nin) {
                int i;
                for(i=0;i<ntypes;++i) 
                    if (!strcmp(nin,names[i])) return(i);
                return(-1);
            }
    };
    const char helper_type::names[ntypes][40] = {"source_changer"};

    
    class ibc_type {
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
    const char ibc_type::names[ntypes][40] = {"zero","power","soi"};

    
}


init_bdry_cndtn *tri_hp_cd::getnewibc(std::string suffix, input_map& inmap) {
    std::string keyword,ibcname;
	init_bdry_cndtn *temp;
    int type;

	/* FIND INITIAL CONDITION TYPE */
	keyword = gbl->idprefix + "_" +suffix;
	if (!inmap.get(keyword,ibcname)) {
		keyword = suffix;
		if (!inmap.get(keyword,ibcname)) {
			*gbl->log << "couldn't find cd initial condition type" << std::endl;
		}
	}
    type = ibc_cd::ibc_type::getid(ibcname.c_str());
 
    switch(type) {
        case(ibc_cd::ibc_type::zero):
            temp = new ibc_cd::zero_src;
			break;
        case(ibc_cd::ibc_type::power):
            temp = new ibc_cd::power_src;
			break;
        case(ibc_cd::ibc_type::soi):
            temp = new ibc_cd::soi_src(*this);
			break;
        default: {
            return(tri_hp::getnewibc(suffix,inmap));
        }
    }
	temp->input(inmap,keyword);
	return(temp);
}

tri_hp_helper *tri_hp_cd::getnewhelper(input_map& inmap) {
    std::string keyword,helpername;
    int type;
    
    /* FIND INITIAL CONDITION TYPE */
    keyword = std::string(gbl->idprefix) + "_helper";
    if (!inmap.get(keyword,helpername)) {
        if (!inmap.get("helper",helpername)) {
            type = -1;
        }
    }
    
    type = ibc_cd::helper_type::getid(helpername.c_str());
        
    switch(type) {
        case ibc_cd::helper_type::source_changer: {
            tri_hp_helper *temp = new ibc_cd::source_changer(*this);
            return(temp);
        }
        default: {
            return(tri_hp::getnewhelper(inmap));
        }
    }
}

