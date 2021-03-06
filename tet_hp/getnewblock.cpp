/*
 *  initfunc.cpp
 *  planar++
 *
 *  Created by helenbrk on Wed Oct 24 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp.h"

//#define MASS
#define CD
#define CD_MULTI
#define INS
#define CNS
//#define CNS_EXPLICIT


#define POD
#ifdef MASS
#include "mass/tet_hp_mass.h"
#endif

#ifdef CD
#include "cd/tet_hp_cd.h"
#endif

#ifdef CD_MULTI
#include "cd_multi/tet_hp_cd_multi.h"
#endif

#ifdef INS
#include "ins/tet_hp_ins.h"
#endif

#ifdef CNS
#include "cns/tet_hp_cns.h"
#endif

#ifdef CNS_EXPLICIT
#include "cns_explicit/tet_hp_cns_explicit.h"
#endif

#ifdef POD
#include "pod/pod_simulate.h"
#include "pod/pod_generate.h"
#include "hp_boundary.h"
#endif

class btype {
	public:
		const static int ntypes = 11;
		enum ids {plain,mass,cd,cd_multi,ins,cns,cns_explicit,pod_ins_gen,pod_cd_gen,pod_ins_sim,pod_cd_sim};
		const static char names[ntypes][40];
		static int getid(const char *nin) {
			int i;
			for(i=0;i<ntypes;++i) 
				if (!strcmp(nin,names[i])) return(i);
			return(-1);
		}
};
const char btype::names[ntypes][40] = {"plain","mass","cd","cd_multi","ins","cns","cns_explicit","pod_ins_gen","pod_cd_gen","pod_ins_sim","pod_cd_sim"};

multigrid_interface* block::getnewlevel(input_map& inmap) {
	std::string keyword,val,ibcname,srcname;
	std::istringstream data;
	int type;          
	
	/* FIND BLOCK TYPE */
	if (inmap.get(idprefix+"_type",val)) {
		type = btype::getid(val.c_str());
	}
	else {
		if (!inmap.get("blocktype",val)) {
			std::cerr << "couldn't find block type" << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
		type = btype::getid(val.c_str());
	}

	switch(type) {
		case btype::plain: {
			tet_hp *temp = new tet_hp();
			return(temp);
		}
			
#ifdef MASS
		case btype::mass: {
			tet_hp_mass *temp = new tet_hp_mass();
			return(temp);
		}
#endif
			
#ifdef CD
		case btype::cd: {
			tet_hp_cd *temp = new tet_hp_cd();
			return(temp);
		}
#endif
			
#ifdef CD_MULTI
		case btype::cd_multi: {
			tet_hp_cd_multi *temp = new tet_hp_cd_multi();
			return(temp);
		}
#endif

#ifdef INS
		case btype::ins: {
			tet_hp_ins *temp = new tet_hp_ins();
			return(temp);
		}
#endif   
	
#ifdef CNS
		case btype::cns: {
			tet_hp_cns *temp = new tet_hp_cns();
			return(temp);
		}
#endif 
		
#ifdef CNS_EXPLICIT
		case btype::cns_explicit: {
			tet_hp_cns_explicit *temp = new tet_hp_cns_explicit();
			return(temp);
		}
#endif 
			
#if (defined(POD) && defined(CD))
		case btype::pod_cd_gen: {
			pod_generate<tet_hp_cd> *temp = new pod_generate<tet_hp_cd>();
			return(temp);
		}
			
		case btype::pod_cd_sim: {
			pod_simulate<tet_hp_cd> *temp = new pod_simulate<tet_hp_cd>();
			return(temp);
		}
#endif
			
#if (defined(POD) && defined(INS))
		case btype::pod_ins_gen: {
			pod_generate<tet_hp_ins> *temp = new pod_generate<tet_hp_ins>();
			return(temp);
		}
			
		case btype::pod_ins_sim: {
			pod_simulate<tet_hp_ins> *temp = new pod_simulate<tet_hp_ins>();
			return(temp);
		}
#endif
		default: {
			std::cerr << "unrecognizable block type: " <<  type << std::endl;
			tet_hp *temp = new tet_hp();
			return(temp);
		}
	} 
		
	return(0);
}

