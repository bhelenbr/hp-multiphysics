/*
 *  allocate.cpp
 *  planar++
 *
 *  Created by helenbrk on Thu Oct 11 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp.h"
#include <input_map.h>
#include "hp_boundary.h"

const int nmovetypes = 5;
const char movetypes[nmovetypes][80] = {"fixed","uncoupled_rigid","coupled_rigid","uncoupled_deformable","coupled_deformable"};

#ifdef AXISYMMETRIC
tri_basis_array<1> basis::tri;
#else
tri_basis_array<0> basis::tri;
#endif

void tri_hp::init(input_map& inmap, shared_ptr<block_global> gin) {
	int i,ival;
	std::string keyword, line;
	std::istringstream data;
	std::string filename;

	gbl = gin;
    hp_gbl = make_shared<hp_global>();

	coarse_flag = false;
	isfrst = true;

	keyword = gbl->idprefix + "_mesh_movement";
	if (!inmap.get(keyword,line)) {
		keyword = "mesh_movement";
		inmap.getwdefault(keyword,line,std::string("fixed"));
	}
	for (i=0;i<nmovetypes;++i)
		if (line == movetypes[i]) break;
	if (i == nmovetypes) 
		*gbl->log << "unrecognized mesh movement type" << std::endl;
	mmovement = static_cast<movementtype>(i);

	/* Initialize stuff for r_tri_mesh */
	if ((mmovement == coupled_deformable) || (mmovement == uncoupled_deformable)) 
		r_tri_mesh::init(inmap,gin);
	else
		tri_mesh::init(inmap,gin);

	keyword = gbl->idprefix + "_nvariable";

	inmap.getwdefault(keyword,NV,1);

	keyword = gbl->idprefix + "_log2p";
	if (!inmap.get(keyword,log2p)) {
		inmap.getwdefault("log2p",log2p,0);
	}
	if (coarse_level) log2p = 0;

#ifdef AXISYMMETRIC
	*gbl->log << "#AXISYMMETRIC is defined" << std::endl;
#else
	*gbl->log << "#AXISYMMETRIC is NOT defined" << std::endl;
#endif
	
#ifdef MESH_REF_VEL
	*gbl->log << "#MESH_REF_VEL is defined" << std::endl;
#else
	*gbl->log << "#MESH_REF_VEL is NOT defined" << std::endl;
#endif

	p0 = basis::tri(log2p)->p();
	sm0 = basis::tri(log2p)->sm();
	im0 = basis::tri(log2p)->im();
	log2pmax = log2p;

	TinyVector<std::string,3> output_purposes;
	TinyVector<int,3> defaults;
	output_purposes(0) = "display_type";
	defaults(0) = tri_hp::tecplot;
	output_purposes(1) = "restart_type";
	defaults(1) = tri_hp::netcdf;
	output_purposes(2) = "debug_type";
	defaults(2) = tri_hp::tecplot;
	for(int i=0;i<3;++i) {
		if (!inmap.get(gbl->idprefix + "_" + output_purposes(i),ival)) inmap.getwdefault(output_purposes(i),ival,defaults(i));
		output_type(i) = static_cast<filetype>(ival);
	}
	if (!inmap.get(gbl->idprefix + "_reload_type",ival)) inmap.getwdefault("reload_type",ival,static_cast<int>(tri_hp::netcdf));
	reload_type = static_cast<filetype>(ival); 

	/* Check that static work arrays are big enough */
	if (u.extent(firstDim) < NV) {
		u.resize(NV);
		res.resize(NV);
		du.resize(NV,ND);
		uht.resize(NV);
		lf.resize(MAX(NV,ND));
		bdwk.resize(gbl->nhist+1,MAX(NV,ND));
	}

	/* Allocate solution vector storage */
	ug.v.resize(maxpst,NV);
	ug.s.resize(maxpst,sm0,NV);
	ug.i.resize(maxpst,im0,NV);

	/* For ease of access have level 0 in time history reference ug */
	ugbd.resize(gbl->nhist+1);
	vrtxbd.resize(gbl->nhist+1);
	ugbd(0).v.reference(ug.v);
	ugbd(0).s.reference(ug.s);
	ugbd(0).i.reference(ug.i);
	vrtxbd(0).reference(pnts); 

	for(i=1;i<gbl->nhist+1;++i) {
		ugbd(i).v.resize(maxpst,NV);
		ugbd(i).s.resize(maxpst,sm0,NV);
		ugbd(i).i.resize(maxpst,im0,NV);
		vrtxbd(i).resize(maxpst);
	}
    
#ifdef ALLCURVED
    inmap.getwdefault(gbl->idprefix+"_allcurved",allcurved,false);
    if (allcurved) {
        /* Allocate solution vector storage */
        crv.s.resize(maxpst,sm0,ND);
        crv.i.resize(maxpst,im0,ND);

        /* For ease of access have level 0 in time history reference ug */
        crvbd.resize(gbl->nhist+1);
        crvbd(0).s.reference(crv.s);
        crvbd(0).i.reference(crv.i);

        for(i=1;i<gbl->nhist+1;++i) {
            crvbd(i).s.resize(maxpst,sm0,ND);
            crvbd(i).i.resize(maxpst,im0,ND);
        }
        pcrdtocht = crdtocht_allcurved;
        pcrdtocht_nhist = crdtocht_nhist_allcurved;
        pcrdtocht1d = crdtocht1d_allcurved;
        pcrdtocht1d_nhist = crdtocht1d_nhist_allcurved;
        
        /* Load curvature coefficients */
        load_curvatures("curvatures",reload_type,crv);
        for(i=1;i<gbl->nhist+1;++i) {
            crvbd(i).s = crv.s;
            crvbd(i).i = crv.i;
        }
    } else {
        pcrdtocht = crdtocht_standard;
        pcrdtocht_nhist = crdtocht_nhist_standard;
        pcrdtocht1d = crdtocht1d_standard;
        pcrdtocht1d_nhist = crdtocht1d_nhist_standard;
    }
#endif

	/* GET INITIAL CONDITION FUNCTION */
	std::string ibcname;
	keyword = gbl->idprefix + "_ibc";
	if (!inmap.get(keyword,ibcname)) {
		keyword = "ibc";
		if (!inmap.get(keyword,ibcname)) {
			*gbl->log << "couldn't find ibc" << std::endl;
		}
	}
	hp_gbl->ibc = getnewibc(ibcname);
	hp_gbl->ibc->init(inmap,keyword);
    
    pmetric = getnewmetric(inmap);
    pmetric->init(inmap);
    
	/* SET MESH VELOCITY TO ZERO */
#ifdef MESH_REF_VEL
	hp_gbl->mesh_ref_vel = 0;
#endif

	/* ALLOCATE BOUNDARY CONDITION STUFF */
	hp_ebdry.resize(nebd);
	hp_vbdry.resize(nvbd);
	for(i=0;i<nebd;++i) {
		keyword =  ebdry(i)->idprefix + "_hp_type";
		std::string val;
		inmap.getwdefault(keyword,val,string("plain"));
		hp_ebdry(i) = getnewedgeobject(i,val);
	}
	for(i=0;i<nvbd;++i) {
		keyword =  vbdry(i)->idprefix + "_hp_type";
		std::string val;
		inmap.getwdefault(keyword,val,string("plain"));
		hp_vbdry(i) = getnewvrtxobject(i,val);
	}
	for(i=0;i<nebd;++i) hp_ebdry(i)->init(inmap);
	for(i=0;i<nvbd;++i) hp_vbdry(i)->init(inmap);
	pmetric->setinfo();

	inmap.getwdefault("hp_fadd",fadd,1.0);

	/* GET HELPER FUNCTION */
	std::string helpername;
	keyword = gbl->idprefix + "_helper";
	if (!inmap.get(keyword,helpername)) {
		keyword = "helper";
		inmap.getwdefault(keyword,helpername,std::string("plain"));
	}
	helper = getnewhelper(helpername);
	helper->init(inmap,keyword);

	/* UNSTEADY SOURCE TERMS */
	dugdt.resize(log2pmax+1);
	dxdt.resize(log2pmax+1);
#ifdef petsc
	int start = log2pmax;
#else
	int start = 0;
#endif	
	for(int i=start;i<=log2pmax;++i) {
		dugdt(i).resize(maxpst,NV,basis::tri(i)->gpx(),basis::tri(i)->gpn());
		dxdt(i).resize(maxpst,ND,basis::tri(i)->gpx(),basis::tri(i)->gpn());
	}
	
	/* To make a vector of residuals that is contiguous in memory */
	//	hp_gbl->res1d.resize(maxpst*(1 +sm0 +im0)*NV);
	//	Array<FLT,2> vtoreference(hp_gbl->res1d.data(),shape(npnt,NV),neverDeleteData);
	//	hp_gbl->res.v.reference(vtoreference);
	//	Array<FLT,3> storeference(hp_gbl->res1d.data() +npnt*NV,shape(nseg,sm0,NV),neverDeleteData);
	//	hp_gbl->res.s.reference(storeference);
	//	Array<FLT,3> itoreference(hp_gbl->res1d.data() +npnt*NV +nseg*sm0*NV,shape(ntri,im0,NV),neverDeleteData);
	//	hp_gbl->res.i.reference(itoreference);
	
	hp_gbl->res.v.resize(maxpst,NV);
	hp_gbl->res.s.resize(maxpst,sm0,NV);
	hp_gbl->res.i.resize(maxpst,im0,NV);
	
	hp_gbl->res_r.v.resize(maxpst,NV);
	hp_gbl->res_r.s.resize(maxpst,sm0,NV);
	hp_gbl->res_r.i.resize(maxpst,im0,NV);
	hp_gbl->res_r.v = 0.;
	hp_gbl->res_r.s = 0.;
	hp_gbl->res_r.i = 0.;

	
#ifndef petsc
	/* Multigrid Storage all except highest order (log2p+1)*/
	dres.resize(log2p);
	
	for(i=0;i<log2p;++i) {
		dres(i).v.resize(maxpst,NV);
		dres(i).s.resize(maxpst,basis::tri(i)->sm(),NV);
		dres(i).i.resize(maxpst,basis::tri(i)->im(),NV);
	}
	
	/* Allocate block stuff */
	hp_gbl->ug0.v.resize(maxpst,NV);
	hp_gbl->ug0.s.resize(maxpst,sm0,NV);
	hp_gbl->ug0.i.resize(maxpst,im0,NV);

	hp_gbl->res0.v.resize(maxpst,NV);
	hp_gbl->res0.s.resize(maxpst,basis::tri(log2p)->sm(),NV);
	hp_gbl->res0.i.resize(maxpst,basis::tri(log2p)->im(),NV); 
#endif
	
	inmap.getwdefault("diagonal_preconditioner",hp_gbl->diagonal_preconditioner,true);
	if (hp_gbl->diagonal_preconditioner) {
		hp_gbl->vprcn.resize(maxpst,NV);
		hp_gbl->sprcn.resize(maxpst,NV);
		hp_gbl->tprcn.resize(maxpst,NV);
	} else {
		hp_gbl->vprcn_ut.resize(maxpst,NV,NV);
		hp_gbl->sprcn_ut.resize(maxpst,NV,NV);
		hp_gbl->tprcn_ut.resize(maxpst,NV,NV);
	}

	double CFLdflt[4] = {2.5, 1.5, 1.0, 0.5};
	if (!inmap.get(gbl->idprefix +"_cfl",hp_gbl->cfl.data(),log2pmax+1)) inmap.getwdefault("cfl",hp_gbl->cfl.data(),log2pmax+1,CFLdflt); 

	/***************************************************/
	/* RESTART SEQUENCE OR INITIAL CONDITION SEQUENCE */
	/**************************************************/
	int restartfile;
	if (inmap.get("restart",restartfile)) {
		std::ostringstream nstr;
		std::string fname;
		nstr << restartfile << std::flush;
		fname = "rstrt" +nstr.str();
		input(fname);
	} 
	else {
		for(i=0;i<nebd;++i)
			hp_ebdry(i)->curv_init();  /* FIXME WILL NEED TO CHANGE THIS TO "tobasis" */

		/* USE TOBASIS TO INITALIZE SOLUTION */
		tobasis(hp_gbl->ibc);
	}

	/*********************************/
	/* ALLOCATE ADAPTATION STORAGE    */
	/*********************************/
	if (gbl->adapt_interval) {
		std::string estring;
		if (!inmap.get(gbl->idprefix + "_error_estimator",estring)) inmap.getwdefault("error_estimator",estring,std::string("none"));
		if (estring == "none") 
			hp_gbl->error_estimator = hp_global::none;
		else if (estring == "energy_norm")
			hp_gbl->error_estimator = hp_global::energy_norm;
		else if (estring == "scale_independent")
			hp_gbl->error_estimator = hp_global::scale_independent;
		else {
			*gbl->log << "Error estimator not recognized" << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
		
		inmap.getwdefault("curvature_sensitivity",hp_gbl->curvature_sensitivity,0.0);
	}
	
#ifdef petsc
	inmap.getwdefault("under_relaxation",under_relax,1.0);
	string petsc_options;
	if (inmap.getline("petsc",petsc_options)) {
		//PetscErrorCode err = PetscOptionsInsertString(petsc_options.c_str());
		PetscErrorCode err = PetscOptionsInsertString(NULL,petsc_options.c_str());
		CHKERRABORT(MPI_COMM_WORLD,err);
	}
	petsc_initialize();
#endif
    
    *gbl->log << "# Initial mesh with DOF: " << npnt +nseg*sm0 +ntri*im0 << std::endl;
	
	return;
}

void tri_hp::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
	std::string keyword;
	std::istringstream data;
	std::string filename;

	const tri_hp& inmesh = dynamic_cast<const tri_hp &>(in);
	gbl = inmesh.gbl;
    hp_gbl = inmesh.hp_gbl;
    pmetric.reset(inmesh.pmetric->create(*this).release());

	/* Initialize stuff for r_tri_mesh */
	mmovement = inmesh.mmovement;
	if (((mmovement == coupled_deformable) || (mmovement == uncoupled_deformable))) 
		r_tri_mesh::init(in,why,sizereduce1d); 
	else
		tri_mesh::init(in,why,sizereduce1d);

	NV = inmesh.NV;

	if (why == multigrid) {
		log2p = 0; 
		p0 = 1;
		sm0 = 0;
		im0 = 0;
		log2pmax = 0;
		coarse_flag = true;
		isfrst = true;
	}
	else {
		log2p = inmesh.log2p;
		p0 = inmesh.p0;
		sm0 = inmesh.sm0;
		im0 = inmesh.im0;
		log2pmax = inmesh.log2pmax;
		coarse_flag = false;
		isfrst = true;
	}
    
#ifdef ALLCURVED
    allcurved = inmesh.allcurved;
    pcrdtocht = inmesh.pcrdtocht;
    pcrdtocht_nhist = inmesh.pcrdtocht_nhist;
    pcrdtocht1d = inmesh.pcrdtocht1d;
    pcrdtocht1d_nhist = inmesh.pcrdtocht1d_nhist;
    
    if (allcurved) {
        /* Allocate solution vector storage */
        crv.s.resize(maxpst,sm0,ND);
        crv.i.resize(maxpst,im0,ND);
        
        /* For ease of access have level 0 in time history reference ug */
        crvbd.resize(gbl->nhist+1);
        crvbd(0).s.reference(crv.s);
        crvbd(0).i.reference(crv.i);
        
        for(int i=1;i<gbl->nhist+1;++i) {
            crvbd(i).s.resize(maxpst,sm0,ND);
            crvbd(i).i.resize(maxpst,im0,ND);
        }
    }
#endif

	output_type = inmesh.output_type;

	/* ALLOCATE WORK ARRAYS */
	u.resize(NV);
	res.resize(NV);
	du.resize(NV,ND);
	uht.resize(NV);
	lf.resize(MAX(NV,ND));
	bdwk.resize(gbl->nhist+1,MAX(NV,ND));

	/* Allocate solution vector storage */
	ug.v.resize(maxpst,NV);
	ug.s.resize(maxpst,sm0,NV);
	ug.i.resize(maxpst,im0,NV);

	/* For ease of access have level 0 in time history reference ug */
	ugbd.resize(gbl->nhist+1);
	vrtxbd.resize(gbl->nhist+1);
	ugbd(0).v.reference(ug.v);
	ugbd(0).s.reference(ug.s);
	ugbd(0).i.reference(ug.i);
	vrtxbd(0).reference(pnts); 
	
	nebd = inmesh.nebd;
	nvbd = inmesh.nvbd;
	hp_ebdry.resize(nebd);
	hp_vbdry.resize(nvbd);
	for(int i=0;i<nebd;++i) 
		hp_ebdry(i) = inmesh.hp_ebdry(i)->create(*this,*ebdry(i));
	for(int i=0;i<nvbd;++i)
		hp_vbdry(i) = inmesh.hp_vbdry(i)->create(*this,*vbdry(i));
	helper = inmesh.helper->create(*this);


	switch (why) {
		case multigrid: {
			/* STUFF FOR MULTIGRID SOURCE TERMS */
#ifdef DIRK
            if (gbl->time_scheme > block_global::AM1) {
                *gbl->log << "Need to recompile to run BD schemes" << std::endl;
                sim::finalize(__LINE__,__FILE__,gbl->log);
            }
			ugbd(1).v.resize(maxpst,NV);
			vrtxbd(1).resize(maxpst); 
#else
			for (int i=1;i<gbl->nhist;++i) {
				ugbd(i).v.resize(maxpst,NV);
				vrtxbd(i).resize(maxpst);
			}
#endif
			vug_frst.resize(maxpst,NV);      
			dres.resize(1);
			dres(0).v.resize(maxpst,NV); 
			fadd = inmesh.fadd;

			/* GET MESH MOVEMENT FUNCTION */
			helper = inmesh.helper->create(*this); 

			/* UNSTEADY SOURCE TERMS */
			dugdt.resize(1);
			dxdt.resize(1);
			dugdt(0).resize(maxpst,NV,basis::tri(0)->gpx(),basis::tri(0)->gpn());
			dxdt(0).resize(maxpst,ND,basis::tri(0)->gpx(),basis::tri(0)->gpn());
            
#ifdef ALLCURVED
            pcrdtocht = crdtocht_standard;
            pcrdtocht_nhist = crdtocht_nhist_standard;
            pcrdtocht1d = crdtocht1d_standard;
            pcrdtocht1d_nhist = crdtocht1d_nhist_standard;
#endif
			break;
		}
		case adapt_storage: {
			for(int i=1;i<gbl->nadapt;++i) {
				ugbd(i).v.resize(maxpst,NV);
				ugbd(i).s.resize(maxpst,sm0,NV);
				ugbd(i).i.resize(maxpst,im0,NV);
				vrtxbd(i).resize(maxpst);
			}
			break;
		}
		default: {
			break;
		}
	}

	return;
}

tri_hp::~tri_hp() {
	
	delete helper;

	for(int i=0;i<nvbd;++i)
		delete hp_vbdry(i);

	for(int i=0;i<nebd;++i)
		delete hp_ebdry(i);
}

#ifndef petsc
FLT tri_hp::maxres() {
	int i,n;
	Array<FLT,1> mxr(NV);
	FLT flowerror;

	/* THIS ROUTINE WILL HAVE TO BE OVERWRITTEN TO GET CORRECT DIMENSIONAL NORM FOR EACH SYSTEM */
	if (mmovement == coupled_deformable) r_tri_mesh::maxres();

	for(n=0;n<NV;++n)
		mxr(n) = 0.0;

	for(i=0;i<npnt;++i) {
		for(n=0;n<NV;++n) {
			mxr(n) = MAX(mxr(n),fabs(hp_gbl->res.v(i,n)));
		}
	}

	for(n=0;n<NV;++n)
		*gbl->log << ' ' << mxr(n) << ' ';

	for(i=0;i<nebd;++i)
		hp_ebdry(i)->maxres();

	flowerror = 0.0;
	for(n=0;n<NV;++n)
		flowerror = MAX(flowerror,mxr(n));

	return(flowerror);

}
#else
FLT tri_hp::maxres() {
//	FLT petsc_norm;
//	VecNorm(petsc_f,NORM_2,&petsc_norm);
//	*gbl->log << ' ' << petsc_norm << ' ';
//	return(petsc_norm);
	
	*gbl->log << ' ' << max_residual << ' ';
	return(max_residual);
}
#endif
