/*
 *  getnewibc_ins.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 2/6/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp_cd.h"
#include <cstring>
#include <math.h>
#include <fstream>
#include <iostream>
using std::ifstream;
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
	public:
		FLT f(int n, TinyVector<FLT,tri_mesh::ND> x,FLT time) {
			return(0.0);
		}
		void init(input_map &inmap,std::string idnty) {}
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
		void init(input_map &inmap,std::string idnty) {
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
		ifstream infile;
		
	public:
		double fq,wth,c,S; //	(frequency (fq), width (wth), power input (c))
		int cycle,source_switch,icount;// R is the random key;g and incount are increment counters;L is the option of data loading
		string binary_filename;
		
		FLT f(int n, TinyVector<FLT,tri_mesh::ND> xpt, FLT time) {
			std::string keyword,val,inf,blockidnum;
			
			
			//*x.gbl->log << fq << ' ' << time << ' ' << source_switch << std::endl;
			int current_cycle = static_cast<int>(fq*(time+1.0e-8));
			if (current_cycle > cycle) {
				// I'm on a new time step!
				
				
				/*Option to load data from a file*/
				if (S==0.0){ 
					
					//Operation to be performed at the end of the file
					if (infile.eof()) {
						infile.close();
						*x.gbl->log << "REACHED THE END OF THE FILE" << std::endl;
						sim::abort(__LINE__,__FILE__,&std::cerr);
					}
					//making sure the file is open
					if (!infile) {
						
						*x.gbl->log << "Not in file "<< std::endl;
						sim::abort(__LINE__,__FILE__,&std::cerr);
						
					}
					else {
						*x.gbl->log << "line in file =  " << icount << std::endl; 
						infile >> source_switch;
						
						
						
						*x.gbl->log <<" Source_switch is = " << source_switch << std::endl;
						if (source_switch==0){
							*x.gbl->log << "#Switch is off this cycle\n"; 
							
						}
						if (source_switch==1){
							*x.gbl->log << "#Switch is on this cycle\n";
							
						}
						
					}
					
				}
				
				/*Switch between Random and Continuous*/
				/*Continuous*/
				if (S==1.0) 	{
					source_switch=1;
					*x.gbl->log << "#Switch is on this cycle\n";
				}
				/*Random*/
				if (S==2.0){
					
	                FLT random_number=rand() % 100;
					if (random_number < 50) {
						source_switch = 0;
						*x.gbl->log << "#Switch is off this cycle\n" ;
						
					} else {
						source_switch = 1;
						*x.gbl->log << "#Switch is on this cycle\n" ;
						
					}
				}   
				
				icount++;
				cycle = current_cycle;
			}
			
			FLT omegat = (2.*M_PI*fq*time+(1/2)*M_PI);
			return(source_switch*c*(exp(cos(omegat)/wth))/(exp(1/wth))); 
			
		}
		soi_src(tri_hp_cd& xin) : x(xin), cycle(-1),icount(0) {}
		void init(input_map &blockdata,std::string idnty) {
			std::string keyword,val;
			std::istringstream data;
			
			
			keyword = idnty+"_simopt";
			blockdata.getwdefault(keyword,S,1.0);
			
			keyword = idnty+"_frequency";
			blockdata.getwdefault(keyword,fq,0.595);
			
			keyword = idnty+"_width";
			blockdata.getwdefault(keyword,wth,0.1);
			
			keyword = idnty+"_constant";
			blockdata.getwdefault(keyword,c,0.89);
			
			
			keyword = x.gbl->idprefix + "_binary_data";
			blockdata.get(keyword,binary_filename);
			infile.open(binary_filename.c_str());
			
		}	
		
		
	};
	
	
	
	
	
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


init_bdry_cndtn *tri_hp_cd::getnewibc(std::string name) {
	std::string keyword,ibcname;
	init_bdry_cndtn *temp;
	int type;
	

	type = ibc_cd::ibc_type::getid(name.c_str());
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
			return(tri_hp::getnewibc(name));
		}
	}
	return(temp);
}





