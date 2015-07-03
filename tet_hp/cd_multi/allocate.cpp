/*
 *  allocatecd.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp_cd_multi.h"
#include <iostream> 

void tet_hp_cd_multi::init(input_map& inmap, void *gin) {

	gbl = static_cast<global *>(gin);
	tet_hp_cd::init(inmap,gin);
	
	
	marks.resize(maxvst);
	
	std::string gridname, filename;
	if (!inmap.get(gbl->idprefix + "_mesh",gridname)) {
		if (inmap.get("mesh",gridname)) {
			gridname = gridname +"_" +gbl->idprefix;
		}
		else {
			*gbl->log << "no mesh name" << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
	}

	size_t dotloc;
	dotloc = gridname.find_last_of('.');
	
	if (dotloc != string::npos) {
		/* Found and ending */
		filename = gridname.substr(0,dotloc) +".marks";
	}
	else {
		filename = gridname +".marks";
	}
	
	ifstream marks_file;
	marks_file.open(filename.c_str());
	if (marks_file.good()) {
		for (int i=0;i<ntet;++i) {
			marks_file.ignore(80,':');
			marks_file >> marks(i);
		}
	}
	
	std::ostringstream nstr;
	for(gbl->nmaterials=0, nstr.str(""), nstr << "Material" << gbl->nmaterials << "_conductivity" << std::flush; inmap.find(nstr.str()) != inmap.end(); ++gbl->nmaterials, nstr.str(""), nstr << "Material" << gbl->nmaterials << "_conductivity" << std::flush);
	
	gbl->kcond.resize(gbl->nmaterials);
	gbl->rhocv.resize(gbl->nmaterials);

	for(int n=0;n<gbl->nmaterials;++n) {
		nstr.str("");
		nstr << "Material" << n << "_conductivity";
		if (!inmap.get(nstr.str(),gbl->kcond(n))) {
			*gbl->log << "Couldn't load " << nstr.str() << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
		
		nstr.str("");
		nstr << "Material" << n << "_rho";
		FLT rho;
		if (!inmap.get(nstr.str(),rho)) {
			*gbl->log << "Couldn't load " << nstr.str() << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
		
		nstr.str("");
		nstr << "Material" << n << "_cv";
		FLT cv;
		if (!inmap.get(nstr.str(),cv)) {
			*gbl->log << "Couldn't load " << nstr.str() << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
		}
		
		gbl->rhocv(n) = rho*cv;
	}

	return;
}

void tet_hp_cd_multi::init(const multigrid_interface& in, init_purpose why, FLT sizereduce1d) {

	const tet_hp_cd_multi& inmesh = dynamic_cast<const tet_hp_cd_multi &>(in);
	gbl = inmesh.gbl;

	tet_hp_cd::init(in,why,sizereduce1d);
	marks.resize(maxvst);
	
	return;
}

/* A GENERIC CALCULATION OF SOURCES FOR AN AUTONOMOUS SYSTEM IN STANDARD FORM */
/* WILL NEED TO BE OVERRIDDEN FOR SPECIAL CASES */
void tet_hp_cd_multi::calculate_unsteady_sources() {
	int lgpx = basis::tet(log2p).gpx, lgpy = basis::tet(log2p).gpy, lgpz = basis::tet(log2p).gpz;
	int i,j,k,n,tind;
	int stridey = MXGP;
	int stridex = MXGP*MXGP;
	TinyVector<int,4> v;
	
	for (log2p=0;log2p<=log2pmax;++log2p) {
		for(tind=0;tind<ntet;++tind) {
			FLT rhocv = gbl->rhocv(marks(tind));
			v = tet(tind).pnt;
			
			if (tet(tind).info > -1) {
				crdtocht(tind,1);
				for(n=0;n<ND;++n)
					basis::tet(log2p).proj_bdry(&cht(n)(0), &crd(n)(0)(0)(0), &dcrd(n)(0)(0)(0)(0), &dcrd(n)(1)(0)(0)(0),&dcrd(n)(2)(0)(0)(0),stridex,stridey);
			}
			else {
				for(n=0;n<ND;++n)
					basis::tet(log2p).proj(vrtxbd(1)(v(0))(n),vrtxbd(1)(v(1))(n),vrtxbd(1)(v(2))(n),vrtxbd(1)(v(3))(n),&crd(n)(0)(0)(0),stridex,stridey);
				
				for(i=0;i<lgpx;++i) {
					for(j=0;j<lgpy;++j) {
						for(k=0;k<lgpz;++k) {
							for(n=0;n<ND;++n) {
								dcrd(n)(0)(i)(j)(k) = 0.5*(vrtxbd(1)(v(3))(n) -vrtxbd(1)(v(2))(n));
								dcrd(n)(1)(i)(j)(k) = 0.5*(vrtxbd(1)(v(1))(n) -vrtxbd(1)(v(2))(n));
								dcrd(n)(2)(i)(j)(k) = 0.5*(vrtxbd(1)(v(0))(n) -vrtxbd(1)(v(2))(n));
							}
						}
					}
				}
			}
			
			ugtouht(tind,1);
			basis::tet(log2p).proj(&uht(0)(0),&u(0)(0)(0)(0),stridex, stridey);
			
			for(i=0;i<lgpx;++i) {
				for(j=0;j<lgpy;++j) {
					for(k=0;k<lgpz;++k) {
						cjcb(i)(j)(k) = -gbl->bd(0)*(dcrd(0)(0)(i)(j)(k)*(dcrd(1)(1)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(1)(i)(j)(k))-dcrd(0)(1)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(2)(i)(j)(k)-dcrd(1)(2)(i)(j)(k)*dcrd(2)(0)(i)(j)(k))+dcrd(0)(2)(i)(j)(k)*(dcrd(1)(0)(i)(j)(k)*dcrd(2)(1)(i)(j)(k)-dcrd(1)(1)(i)(j)(k)*dcrd(2)(0)(i)(j)(k)));
						
						dugdt(log2p,tind,0)(i)(j)(k) = rhocv*u(0)(i)(j)(k)*cjcb(i)(j)(k);
						for(n=0;n<ND;++n)
							dxdt(log2p,tind,n)(i)(j)(k) = crd(n)(i)(j)(k);
					}
				}
			}
		}
		
	}
	log2p = log2pmax;
	
	return;
}


