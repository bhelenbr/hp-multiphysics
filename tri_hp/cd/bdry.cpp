/*
 *  cd_bdry.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 1/13/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

/*
 *  cd_bdry.h
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/7/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "bdry_cd.h"
#include "myblas.h"

//#define MPDEBUG

using namespace bdry_cd;

void generic::output(const std::string& filename, tri_hp::filetype typ,int tlvl) {
	int i,m,n,ind,sind,tind,seg;
	FLT visc[tri_mesh::ND];
	TinyVector<FLT,tri_mesh::ND> norm, mvel;
	FLT l,conv,diff,jcb;
	
	std::string fname;
	fname = filename +"_" +base.idprefix;
	
	switch(typ) {
		case(tri_hp::text): case(tri_hp::netcdf): case(tri_hp::binary): {
			hp_edge_bdry::output(filename,typ,tlvl);
			break;
		}
		case(tri_hp::tecplot): {
			if (!report_flag) return;
			hp_edge_bdry::output(filename,typ,tlvl);
			break;
			
			fname += ".dat";
			std::ofstream fout;
			fout.open(fname.c_str());
			
			fout << "VARIABLES=\"S\",\"X\",\"Y\",\"CFLUX\",\"DFLUX\"\nTITLE = " << base.idprefix << '\n'<< "ZONE\n";
			
			conv_total = 0.0;
			diff_total = 0.0;
			circumference = 0.0;
			ind = 0; 
			do { 
				sind = base.seg(ind);
				tind = x.seg(sind).tri(0);        
				
				for(seg=0;seg<3;++seg)
					if (x.tri(tind).seg(seg) == sind) break;
				assert(seg != 3);
				
				x.crdtocht(tind);
				for(m=basis::tri(x.log2p)->bm();m<basis::tri(x.log2p)->tm();++m)
					for(n=0;n<tri_mesh::ND;++n)
						x.cht(n,m) = 0.0;
				
				for(n=0;n<tri_mesh::ND;++n)
					basis::tri(x.log2p)->proj_side(seg,&x.cht(n,0), &x.crd(n)(0,0), &x.dcrd(n,0)(0,0), &x.dcrd(n,1)(0,0));
				
				x.ugtouht(tind);
				for(n=0;n<x.NV;++n)
					basis::tri(x.log2p)->proj_side(seg,&x.uht(n)(0),&x.u(n)(0,0),&x.du(n,0)(0,0),&x.du(n,1)(0,0));
				
				for (i=0;i<basis::tri(x.log2p)->gpx();++i) {
					jcb =  basis::tri(x.log2p)->wtx(i)*RAD(x.crd(0)(0,i))*sqrt(x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i) +x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i));
					circumference += jcb;
					
					x.cjcb(0,i) = x.hp_cd_gbl->kcond/(x.dcrd(0,0)(0,i)*x.dcrd(1,1)(0,i) -x.dcrd(1,0)(0,i)*x.dcrd(0,1)(0,i));
										
					/* DIFFUSIVE FLUXES ( FOR EXTRA VARIABLES) */
					visc[0] =  x.cjcb(0,i)*(x.dcrd(1,1)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,1)(0,i)*x.dcrd(0,0)(0,i));
					visc[1] = -x.cjcb(0,i)*(x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i));
					l = sqrt(x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i)  +x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i));
					norm(0) = x.dcrd(1,0)(0,i)/l;
					norm(1) = -x.dcrd(0,0)(0,i)/l; 
					for(n=0;n<tri_mesh::ND;++n) {
						mvel(n) = x.gbl->bd(0)*(x.crd(n)(0,i) -dxdt(x.log2p,ind)(n,i));
#ifdef MESH_REF_VEL
						mvel(n) += x.hp_gbl->mesh_ref_vel(n);
#endif
					}					
#ifdef CONST_A
					conv = x.u(0)(0,i)*((x.hp_cd_gbl->ax -mvel(0))*norm(0) +(x.hp_cd_gbl->ay -mvel(1))*norm(1));
#else
					TinyVector<FLT,tri_mesh::ND> pt;
					pt(0) = x.crd(0)(0,i);
					pt(1) = x.crd(1)(0,i);
					conv = x.u(0)(0,i)*((x.hp_cd_gbl->a->f(0,pt,x.gbl->time) -mvel(0))*norm(0) +(x.hp_cd_gbl->a->f(1,pt,x.gbl->time) -mvel(1))*norm(1));
#endif
					conv_total -= basis::tri(x.log2p)->wtx(i)*RAD(x.crd(0)(0,i))*conv*l;
					
					diff = -visc[0]*x.du(0,0)(0,i) -visc[1]*x.du(0,1)(0,i)/l;
					
					diff_total -= basis::tri(x.log2p)->wtx(i)*RAD(x.crd(0)(0,i))*diff*l;
					
					fout << circumference << ' ' << x.crd(0)(0,i) << ' ' << x.crd(1)(0,i) << ' ' << conv << ' ' << diff << std::endl;

				}	
			} while (++ind < base.nseg);
			streamsize oldprecision = (*x.gbl->log).precision(10);
			*x.gbl->log << base.idprefix << " circumference: " << circumference << std::endl;
			*x.gbl->log << base.idprefix << " total diffusive flux: " << diff_total << std::endl;
			*x.gbl->log << base.idprefix << " total convective flux: " << conv_total << std::endl;
			(*x.gbl->log).precision(oldprecision);
			
			fout.close();
			break;
		}
		default: 
			break;
	}
	
	return;
}



