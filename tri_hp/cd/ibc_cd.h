/*
 *  ibc_cd.h
 *  tri_hp
 *
 *  Source term functions for the temperature equation.  Shared by
 *  the cd and thermoelastic modules.
 *
 */

#ifndef _ibc_cd_h_
#define _ibc_cd_h_

#include "../tri_hp.h"
#include <math.h>
#include <fstream>
#include <iostream>
using std::ifstream;

namespace ibc_cd {
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
		tri_hp &x;
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
			
			FLT omegat = (2.*M_PI*fq*time+(1./2.)*M_PI);
			return(source_switch*c*(exp(cos(omegat)/wth))/(exp(1/wth)));
			
		}
		soi_src(tri_hp& xin) : x(xin), cycle(-1),icount(0) {}
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
}
#endif
