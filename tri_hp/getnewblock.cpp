/*
 *  initfunc.cpp
 *  planar++
 *
 *  Created by helenbrk on Wed Oct 24 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include <blocks.h>
#include <block.h>
#include <r_tri_mesh.h>

//#define CD
#define INS
//#define PS
//#define SWIRL
#define BUOYANCY
//#define SWE
//#define LVLSET
//#define EXPLICIT
//#define CNS

#define POD

#ifdef CD
#include "cd/tri_hp_cd.h"
#endif

#ifdef INS
#include "ins/tri_hp_ins.h"
#endif

#ifdef PS
#include "planestrain/tri_hp_ps.h"
#endif

#ifdef SWIRL
#include "swirl/tri_hp_swirl.h"
#endif

#ifdef BUOYANCY
#include "buoyancy/tri_hp_buoyancy.h"
#endif

#ifdef SWE
#include "swe/tri_hp_swe.h"
#endif

#ifdef LVLSET
#include "lvlset/tri_hp_lvlset.h"
#endif

#ifdef POD
#include "pod/pod_simulate.h"
#include "pod/pod_generate.h"
#include "hp_boundary.h"
#endif

#ifdef EXPLICIT
#include "explicit/tri_hp_explicit.h"
#endif

#ifdef CNS
#include "cns/tri_hp_cns.h"
#endif

class btype {
	public:
		const static int ntypes = 16;
		enum ids {r_tri_mesh,cd,ins,ps,swirl,buoyancy,pod_ins_gen,pod_cd_gen,pod_cns_gen,pod_ins_sim,pod_cns_sim,pod_cd_sim,swe,lvlset,explct,cns};
		const static char names[ntypes][40];
		static int getid(const char *nin) {
			int i;
			for(i=0;i<ntypes;++i) 
				if (!strcmp(nin,names[i])) return(i);
			return(-1);
		}
};
const char btype::names[ntypes][40] = {"r_tri_mesh","cd","ins","ps","swirl","buoyancy",
    "pod_ins_gen","pod_cd_gen","pod_cns_gen","pod_ins_sim","pod_cns_sim","pod_cd_sim","swe","lvlset","explicit","cns"};

multigrid_interface* block::getnewlevel(input_map& input) {
	std::string keyword,val,ibcname;
	std::istringstream data;
	int type;          

	/* FIND BLOCK TYPE */
	if (input.get(idprefix+"_type",val)) {
		type = btype::getid(val.c_str());
	}
	else {
		if (!input.get("blocktype",val)) {
			std::cerr << "couldn't find block type" << std::endl;
			exit(1);
		}
		type = btype::getid(val.c_str());
	}

	switch(type) {
		case btype::r_tri_mesh: {
			r_tri_mesh *temp = new r_tri_mesh();
			return(temp);
		}
#ifdef CD
		case btype::cd: {
			tri_hp_cd *temp = new tri_hp_cd();
			return(temp);
		}
#endif

#ifdef INS
		case btype::ins: {
			tri_hp_ins *temp = new tri_hp_ins();
			return(temp);
		}
#endif

#ifdef PS
		case btype::ps: {
			tri_hp_ps *temp = new tri_hp_ps();
			return(temp);
		}
#endif

#ifdef SWIRL
		case btype::swirl: {
			tri_hp_swirl *temp = new tri_hp_swirl();
			return(temp);
		}
#endif

#ifdef BUOYANCY
		case btype::buoyancy: {
			tri_hp_buoyancy *temp = new tri_hp_buoyancy();
			return(temp);
		}
#endif

#ifdef SWE
		case btype::swe: {
			tri_hp_swe *temp = new tri_hp_swe();
			return(temp);
		}
#endif

#ifdef LVLSET
		case btype::lvlset: {
			tri_hp_lvlset *temp = new tri_hp_lvlset();
			return(temp);
		}
#endif

#if (defined(POD) && defined(CD))
		case btype::pod_cd_gen: {
			pod_generate<tri_hp_cd> *temp = new pod_generate<tri_hp_cd>();
			return(temp);
		}

		case btype::pod_cd_sim: {
			pod_simulate<tri_hp_cd> *temp = new pod_simulate<tri_hp_cd>();
			return(temp);
		}
#endif

#if (defined(POD) && defined(INS))
		case btype::pod_ins_gen: {
			pod_generate<tri_hp_ins> *temp = new pod_generate<tri_hp_ins>();
			return(temp);
		}

		case btype::pod_ins_sim: {
			pod_simulate<tri_hp_ins> *temp = new pod_simulate<tri_hp_ins>();
			return(temp);
		}
#endif

#if (defined(POD) && defined(CNS))
		case btype::pod_cns_gen: {
			pod_generate<tri_hp_cns> *temp = new pod_generate<tri_hp_cns>();
			return(temp);
		}
			
		case btype::pod_cns_sim: {
			pod_simulate<tri_hp_cns> *temp = new pod_simulate<tri_hp_cns>();
			return(temp);
		}
#endif
			
#ifdef EXPLICIT
		case btype::explct: {
			tri_hp_explicit *temp = new tri_hp_explicit();
			return(temp);
		}
#endif

#ifdef CNS
		case btype::cns: {
			tri_hp_cns *temp = new tri_hp_cns();
			return(temp);
		}
#endif
		default: {
			std::cerr << "unrecognizable block type: " <<  type << std::endl;
			r_tri_mesh *temp = new r_tri_mesh();
			return(temp);
		}
	} 

	return(0);
}

