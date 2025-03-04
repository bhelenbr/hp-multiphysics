/*
 *  allocate.cpp
 *  planar++
 *
 *  Created by helenbrk on Thu Oct 11 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */


#include "tet_hp.h"
#include <input_map.h>
#include <boundary.h>
#include "hp_boundary.h"

const int nmovetypes = 5;
const char movetypes[nmovetypes][80] = {"fixed","uncoupled_rigid","coupled_rigid","uncoupled_deformable","coupled_deformable"};

void tet_hp::init(input_map& inmap, std::shared_ptr<block_global> gin) {
	int i,ival,p;
	std::string keyword, line;
	std::istringstream data;
	std::string filename;
	
	gbl = gin;
	
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
	
	/* Initialize stuff for r_mesh */
//   if ((mmovement == coupled_deformable) || (mmovement == uncoupled_deformable)) temporary
//     r_tet_mesh::init(inmap,gin); 
//   else
		tet_mesh::init(inmap,gin);

	keyword = gbl->idprefix + "_nvariable";
	
	inmap.getwdefault(keyword,NV,1);
	
	keyword = gbl->idprefix + "_log2p";
	if (!inmap.get(keyword,log2p)) {
		inmap.getwdefault("log2p",log2p,0);
	}
	if (coarse_level) log2p = 0;
	

	/* Check that global basis has been allocated & allocate if necessary */
	if (basis::tet.extent(firstDim) < log2p +1) {
		basis::tet.resize(log2p+1);

		p = 1;
		for(i=0;i<log2p+1;++i) {
			basis::tet(i).initialize(p,p+1);
			p = p+1;
		}
	}
	
	p0 = basis::tet(log2p).p;
	em0 = basis::tet(log2p).em;
	fm0 = basis::tet(log2p).fm;
	im0 = basis::tet(log2p).im;
	log2pmax = log2p;
	
	TinyVector<std::string,3> output_purposes;
	TinyVector<int,3> defaults;
	output_purposes(0) = "display_type";
	defaults(0) = tet_hp::tecplot;
	output_purposes(1) = "restart_type";
	defaults(1) = tet_hp::binary;
	output_purposes(2) = "debug_type";
	defaults(2) = tet_hp::tecplot;
	for(int i=0;i<3;++i) {
		if (!inmap.get(gbl->idprefix + "_" + output_purposes(i),ival)) inmap.getwdefault(output_purposes(i),ival,defaults(i));
		output_type(i) = static_cast<filetype>(ival);
	}
	if (!inmap.get(gbl->idprefix + "_reload_type",ival)) inmap.getwdefault("reload_type",ival,static_cast<int>(tet_hp::binary));
	reload_type = static_cast<filetype>(ival); 
		
	/* Check that static work arrays are big enough */
	if (u.extent(firstDim) < NV) {
		u.resize(NV);
		u2d.resize(NV);
		u1d.resize(NV);
		res.resize(NV);
		res2d.resize(NV);
		res1d.resize(NV);
		du.resize(NV,ND);
		uht.resize(NV);
		lf.resize(MAX(NV,ND));
		//bdwk.resize(gbl->nhist+1,MAX(NV,ND));//temporary doesnt match .h
	}
	
		
	/* Allocate solution vector storage */
	ug.v.resize(maxvst,NV);
	ug.e.resize(maxvst,em0,NV);
	ug.f.resize(maxvst,fm0,NV);
	ug.i.resize(maxvst,im0,NV);
	
	/* For ease of access have level 0 in time history reference ug */
	ugbd.resize(gbl->nhist+1);
	vrtxbd.resize(gbl->nhist+1);
	ugbd(0).v.reference(ug.v);
	ugbd(0).e.reference(ug.e);
	ugbd(0).f.reference(ug.f);
	ugbd(0).i.reference(ug.i);
	vrtxbd(0).reference(pnts); 
	
	for(i=1;i<gbl->nhist+1;++i) {
		ugbd(i).v.resize(maxvst,NV);
		ugbd(i).e.resize(maxvst,em0,NV);
		ugbd(i).f.resize(maxvst,fm0,NV);
		ugbd(i).i.resize(maxvst,im0,NV);
		vrtxbd(i).resize(maxvst);
	} 
	
	/* GET INITIAL CONDITION FUNCTION */
	keyword = gbl->idprefix + "_ibc";
	std::string ibcname;
	if (!inmap.get(keyword,ibcname)) {
		keyword = "ibc";
		if (!inmap.get(keyword,ibcname)) {
			*gbl->log << "couldn't find ibc" << std::endl;
		}
	}
	hp_gbl->ibc = getnewibc(ibcname);
	hp_gbl->ibc->init(inmap,keyword);
	
	hp_fbdry.resize(nfbd);
	hp_ebdry.resize(nebd);
	hp_vbdry.resize(nvbd);
	for(i=0;i<nfbd;++i) {
		keyword =  fbdry(i)->idprefix + "_hp_type";
		std::string val;
		inmap.getwdefault(keyword,val,string("plain"));
		hp_fbdry(i) = getnewfaceobject(i,val);
	}
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
	for(i=0;i<nfbd;++i) hp_fbdry(i)->init(inmap);
	for(i=0;i<nebd;++i) hp_ebdry(i)->init(inmap);
	for(i=0;i<nvbd;++i) hp_vbdry(i)->init(inmap);
	setinfo();
	
	inmap.getwdefault("hp_fadd",fadd,1.0);

	/* GET HELPER FUNCTION */
	std::string helpername;
	if (!inmap.get(gbl->idprefix + "_helper",helpername))
		inmap.getwdefault("helper",helpername,std::string("plain"));
	helper = getnewhelper(helpername);
	helper->init(inmap,helpername);
	
	/* UNSTEADY SOURCE TERMS */
	dugdt.resize(log2p+1,maxvst,NV);
	dxdt.resize(log2p+1,maxvst,ND);
	/// COARSE MESH STOPS HERE */

	/* Multigrid Storage all except highest order (log2p+1)*/
	dres.resize(log2p);
	for(i=0;i<log2p;++i) {
		dres(i).v.resize(maxvst,NV);
		dres(i).e.resize(maxvst,basis::tet(i).em,NV);
		dres(i).f.resize(maxvst,basis::tet(i).fm,NV);
		dres(i).i.resize(maxvst,basis::tet(i).im,NV);
	}
		
	/* Allocate block stuff */
	hp_gbl->ug0.v.resize(maxvst,NV);
	hp_gbl->ug0.e.resize(maxvst,em0,NV);
	hp_gbl->ug0.f.resize(maxvst,fm0,NV);
	hp_gbl->ug0.i.resize(maxvst,im0,NV);
		
	hp_gbl->res.v.resize(maxvst,NV);
	hp_gbl->res.e.resize(maxvst,em0,NV);
	hp_gbl->res.f.resize(maxvst,fm0,NV);
	hp_gbl->res.i.resize(maxvst,im0,NV);
	
	hp_gbl->res_r.v.resize(maxvst,NV);
	hp_gbl->res_r.e.resize(maxvst,em0,NV);
	hp_gbl->res_r.f.resize(maxvst,fm0,NV);
	hp_gbl->res_r.i.resize(maxvst,im0,NV);
	hp_gbl->res_r.v = 0.0;
	hp_gbl->res_r.e = 0.0;
	hp_gbl->res_r.f = 0.0;
	hp_gbl->res_r.i = 0.0;

	hp_gbl->res0.v.resize(maxvst,NV);
	hp_gbl->res0.e.resize(maxvst,basis::tet(log2p).em,NV);
	hp_gbl->res0.f.resize(maxvst,basis::tet(log2p).fm,NV);
	hp_gbl->res0.i.resize(maxvst,basis::tet(log2p).im,NV); 

#ifdef JACOBI
	hp_gbl->jacob_diag.v.resize(maxvst,NV);
	hp_gbl->jacob_diag.e.resize(maxvst,em0,NV);
	hp_gbl->jacob_diag.f.resize(maxvst,fm0,NV);
	hp_gbl->jacob_diag.i.resize(maxvst,im0,NV);
#endif
	
	inmap.getwdefault("diagonal_preconditioner",hp_gbl->diagonal_preconditioner,true);
	if (hp_gbl->diagonal_preconditioner) {
		hp_gbl->vprcn.resize(maxvst,NV);
		hp_gbl->eprcn.resize(maxvst,NV);
		hp_gbl->fprcn.resize(maxvst,NV);
		hp_gbl->iprcn.resize(maxvst,NV);
	} else {
		hp_gbl->vprcn_ut.resize(maxvst,NV,NV);
		hp_gbl->eprcn_ut.resize(maxvst,NV,NV);
		hp_gbl->fprcn_ut.resize(maxvst,NV,NV);
		hp_gbl->iprcn_ut.resize(maxvst,NV,NV);
	}

	double CFLdflt[4] = {2.2, 1.5, 1.0, 0.5};
	if (!inmap.get(gbl->idprefix +"_cfl",hp_gbl->cfl.data(),log2pmax+1)) inmap.getwdefault("cfl",hp_gbl->cfl.data(),log2pmax+1,CFLdflt); 

	/***************************************************/
	/* RESTART SEQUENCE OR INITIAL CONDITION SEQUENCE */
	/**************************************************/
	int restartfile;
	if (inmap.get("restart",restartfile)) {
		std::ostringstream nstr;
		std::string fname;
		nstr << restartfile << std::flush;
		fname = "rstrt" +nstr.str() +"_" +gbl->idprefix;
		input(fname);
	} 
	else {
		for(i=0;i<nebd;++i)
			hp_ebdry(i)->curv_init();  /* FIXME WILL NEED TO CHANGE THIS TO "tobasis" */
		for(i=0;i<nfbd;++i)
			hp_fbdry(i)->curv_init();  /* FIXME WILL NEED TO CHANGE THIS TO "tobasis" */
		
		/* USE TOBASIS TO INITALIZE SOLUTION */
		tobasis(hp_gbl->ibc);
	}
	
	/* fix volume of tets after they are curved */
	if(basis::tet(log2p).p > 1) {
		for(i=0;i<ntet;++i) {
			tet(i).vol = volume_curved(i); /* which is actually 8*volume */
		}
	}
	
	if(basis::tet(log2p).p == 3){
		spkmass.resize(npnt);
		spklink.resize(npnt);
		spkpiv.resize(npnt);
		int maxnspk = 0;
		int nspk;
		for(int i=0; i < npnt; ++i){
			nspk=pnt(i).nspk;
			if (nspk > maxnspk)
				maxnspk = nspk;		
			spkmass(i).resize(nspk*(basis::tet(log2p).em-1),nspk*(basis::tet(log2p).em-1));
			spklink(i).resize(nspk);
			spkpiv(i).resize(nspk*(basis::tet(log2p).em-1));
		}
		*gbl->log << "max number of spokes "<<  maxnspk << endl;
		spkres.resize(maxnspk*(basis::tet(log2p).em-1));
		wkseg.resize(nseg,basis::tet(log2p).em-1,NV);
		//spoke();
	}
	

//#ifndef petsc
//	size_sparse_matrix = (npnt+nseg*em0+ntri*fm0+ntet*im0)*NV;
//
//	/* sparse matrix allocation */
//	ija.resize(MXTM*NV*ntet);//too much storage resize later
//	sa.resize(MXTM*NV*ntet);
//	number_sparse_elements = size_sparse_matrix;
//	sa = 0.0;
//	/* creates sparse matrix with zeros on diagonal */
//	for(int i = 0; i < number_sparse_elements+1; ++i)
//		ija(i) = size_sparse_matrix+1;
//#endif	
	
//#ifdef petsc
//	petsc_initialize();
//#else
//	initialize_sparse();
//#endif

	//test();
	
	/*********************************/
	/* ALLOCATE ADAPTATION STORAGE   */
	/*********************************/
	if (gbl->adapt_interval) {
		inmap.getwdefault("curvature_sensitivity",hp_gbl->curvature_sensitivity,20.0);
		hp_gbl->pstr = create();
		hp_gbl->pstr->init(*this,adapt_storage);

		/* LET EACH BOUNDARY CONDITION DIRECTLY FIND ITS ADAPTATION STORAGE */
		for(i=0;i<nebd;++i)
		hp_ebdry(i)->adapt_storage = hp_gbl->pstr->hp_ebdry(i);
	}

	return;
}


void tet_hp::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {
	std::string keyword;
	std::istringstream data;
	std::string filename;
	
	const tet_hp& inmesh = dynamic_cast<const tet_hp &>(in);
	gbl = inmesh.gbl;
	
	/* Initialize stuff for r_mesh */
//   mmovement = inmesh.mmovement;
//   if ((mmovement == coupled_deformable) || (mmovement == uncoupled_deformable)) temporary
//     r_tet_mesh::init(in,sizereduce1d); 
//   else
		tet_mesh::init(in,why,sizereduce1d);

	NV = inmesh.NV;

	if (why == multigrid) {
		log2p = 0; 
		p0 = 1;
		em0 = 0;
		fm0 = 0;
		im0 = 0;
		log2pmax = 0;
		coarse_flag = true;
		isfrst = true;
	}
	else {
		log2p = inmesh.log2p;
		p0 = inmesh.p0;
		em0 = inmesh.em0;
		fm0 = inmesh.fm0;
		im0 = inmesh.im0;
		log2pmax = inmesh.log2pmax;
		coarse_flag = false;
		isfrst = true;
	}

	
	output_type = inmesh.output_type;

	/* ALLOCATE WORK ARRAYS */
	u.resize(NV);
	u2d.resize(NV);
	u1d.resize(NV);
	res.resize(NV);
	res2d.resize(NV);
	res1d.resize(NV);
	du.resize(NV,ND);
	uht.resize(NV);
	lf.resize(MAX(NV,ND));
	// bdwk.resize(gbl->nhist+1,MAX(NV,ND)); TEMPORARY
	
	/* Allocate solution vector storage */
	ug.v.resize(maxvst,NV);
	ug.e.resize(maxvst,em0,NV);
	ug.f.resize(maxvst,fm0,NV);
	ug.i.resize(maxvst,im0,NV);
	
	/* For ease of access have level 0 in time history reference ug */
	ugbd.resize(gbl->nhist+1);
	vrtxbd.resize(gbl->nhist+1);
	ugbd(0).v.reference(ug.v);
	ugbd(0).e.reference(ug.e);
	ugbd(0).f.reference(ug.f);
	ugbd(0).i.reference(ug.i);
	vrtxbd(0).reference(pnts); 

	nfbd = inmesh.nfbd;
	nebd = inmesh.nebd;
	nvbd = inmesh.nvbd;
	hp_fbdry.resize(nfbd);
	hp_ebdry.resize(nebd);
	hp_vbdry.resize(nvbd);
	for(int i=0;i<nfbd;++i) hp_fbdry(i) = inmesh.hp_fbdry(i)->create(*this,*fbdry(i));
	for(int i=0;i<nebd;++i) hp_ebdry(i) = inmesh.hp_ebdry(i)->create(*this,*ebdry(i));
	for(int i=0;i<nvbd;++i) hp_vbdry(i) = inmesh.hp_vbdry(i)->create(*this,*vbdry(i));

	switch (why) {
		case multigrid: {
			/* STUFF FOR MULTIGRID SOURCE TERMS */
#ifdef DIRK
			ugbd(1).v.resize(maxvst,NV);
			vrtxbd(1).resize(maxvst); 
#else
			for (int i=1;i<gbl->nhist;++i) {
				ugbd(i).v.resize(maxvst,NV);
				vrtxbd(i).resize(maxvst);
			}
#endif
			vug_frst.resize(maxvst,NV);    
			dres.resize(1);
			dres(0).v.resize(maxvst,NV); 
			fadd = inmesh.fadd;
			
			/* GET MESH MOVEMENT FUNCTION */
			helper = inmesh.helper->create(*this); 
			
			/* UNSTEADY SOURCE TERMS */
			dugdt.resize(log2p+1,maxvst,NV);
			dxdt.resize(log2p+1,maxvst,ND);
			break;
		}
		case adapt_storage: {
			for(int i=1;i<gbl->nadapt;++i) {
				ugbd(i).v.resize(maxvst,NV);
				ugbd(i).e.resize(maxvst,em0,NV);
				ugbd(i).f.resize(maxvst,fm0,NV);
				ugbd(i).i.resize(maxvst,im0,NV);
				vrtxbd(i).resize(maxvst);
			}
			break;
		}
		default: {
			break;
		}
	}
	
	return;
}

tet_hp::~tet_hp() {
	
	for(int i=0;i<nvbd;++i)
		delete hp_vbdry(i);

	for(int i=0;i<nebd;++i)
		delete hp_ebdry(i);
		
	for(int i=0;i<nfbd;++i)
		delete hp_fbdry(i);
}

/* pnt.info marked with bdry number */
/* face boundary number in Tri.tet(1) & tet.tet(X)*/
/* seg.info marked with bdry number  */
/* tet.info marked with 0/-1 */

void tet_hp::setinfo() {
	int i;
	
	/* SET UP VRTX BC INFORMATION FOR OUTPUT */
	for(i=0;i<npnt;++i)
		pnt(i).info = -1;

	/* SET UP SIDE BC INFORMATION FOR CURVED SIDES OUTPUT */
	for(i=0;i<nseg;++i)
		seg(i).info = -1;
	
	for(i=0;i<ntri;++i)
		tri(i).info = -1;
		
	for(i=0;i<ntet;++i)
		tet(i).info = -1;

	for(i=0;i<nvbd;++i)
		pnt(vbdry(i)->pnt).info = i;
	
	if(log2p > 0){		
		for(i=0;i<nfbd;++i){
			if(hp_fbdry(i)->is_curved()){
				for(int j=0;j<fbdry(i)->ntri;++j){
					tri(fbdry(i)->tri(j).gindx).info = 0; /* tag segment as curved on face */
					tri(fbdry(i)->tri(j).gindx).bnum = i;
					tri(fbdry(i)->tri(j).gindx).lnum = j;
				}
				for(int j=0;j<fbdry(i)->nseg;++j){
					seg(fbdry(i)->seg(j).gindx).info = 0; /* tag segment as curved on face */
					seg(fbdry(i)->seg(j).gindx).bnum = i;
					seg(fbdry(i)->seg(j).gindx).lnum = j;
				}
			}
		}
		
		for(i=0;i<nebd;++i){
			if(hp_ebdry(i)->is_curved()){
				for(int j=0;j<ebdry(i)->nseg;++j){
					seg(ebdry(i)->seg(j).gindx).info = 1; /* overide face curvature tag */
					seg(ebdry(i)->seg(j).gindx).bnum = i;
					seg(ebdry(i)->seg(j).gindx).lnum = j;
				}
			}
		}
		
		/* tag any tet that has something curvy on it */
		for(i=0;i<ntet;++i){
			for(int j = 0; j < 4; ++j){
				int find = tet(i).tri(j);
				if(tri(find).info != -1){
					tet(i).info = 1;
					break;
				}
			}
			for(int j = 0; j < 6; ++j){
				int eind = tet(i).seg(j);
				if(seg(eind).info != -1){
					tet(i).info = 1;
					break;
				}
			}
		}
	}
	
	return;
}

FLT tet_hp::maxres() {
	int i,n;
	Array<FLT,1> mxr(NV);
	FLT flowerror;
	
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
		
	for(i=0;i<nfbd;++i)
		hp_fbdry(i)->maxres();

	flowerror = 0.0;
	for(n=0;n<NV;++n)
		flowerror = MAX(flowerror,mxr(n));

	return(flowerror);
	
}


/* calculate jacobian including curvature */
FLT tet_hp::jacobian(int tind) {
	
	TinyVector<int,4> v;
	int lgpx = basis::tet(log2p).gpx, lgpy = basis::tet(log2p).gpy, lgpz = basis::tet(log2p).gpz;
	int stridey = MXGP;
	int stridex = MXGP*MXGP; 
	FLT cjcb = 0.0;
	
	
	/* LOAD INDICES OF VERTEX POINTS */
	v = tet(tind).pnt;
	
	if (tet(tind).info < 0) {		
		for(int i=0;i<lgpx;++i) {
			for(int j=0;j<lgpy;++j) {
				for(int k=0;k<lgpz;++k) {
					for(int n=0;n<ND;++n) {
						dcrd(n)(0)(i)(j)(k) = 0.5*(pnts(tet(tind).pnt(3))(n) -pnts(tet(tind).pnt(2))(n));
						dcrd(n)(1)(i)(j)(k) = 0.5*(pnts(tet(tind).pnt(1))(n) -pnts(tet(tind).pnt(2))(n));
						dcrd(n)(2)(i)(j)(k) = 0.5*(pnts(tet(tind).pnt(0))(n) -pnts(tet(tind).pnt(2))(n));
					}
				}
			}
		}
	}
	else {
		crdtocht(tind);
		for(int n=0;n<ND;++n)
			basis::tet(log2p).proj_bdry(&cht(n)(0), &crd(n)(0)(0)(0), &dcrd(n)(0)(0)(0)(0), &dcrd(n)(1)(0)(0)(0),&dcrd(n)(2)(0)(0)(0),stridex,stridey);
	}
	
	for(int i=0;i<lgpx;++i) {
		for(int j=0;j<lgpy;++j) {
			for(int k=0;k<lgpz;++k) {
				cjcb += basis::tet(log2p).wtx(i)*basis::tet(log2p).wty(j)*basis::tet(log2p).wtz(k)*(dcrd(0)(0)(i)(j)(k)*(dcrd(1)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(1)(i)(j)(k))-dcrd(0)(1)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(0)(i)(j)(k))+dcrd(0)(2)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)-dcrd(1)(1)(i)(j)(k)*dcrd(2)(0)(i)(j)(k)));

			}
		}
	} 
		
	return(cjcb);
}

