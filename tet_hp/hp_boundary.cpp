/*
 *  hp_boundary.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/2/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp.h"
#include "hp_boundary.h"
#include <blitz/tinyvec-et.h>
#include <libbinio/binwrap.h>
#include <myblas.h>

//#define MPDEBUG

hp_vrtx_bdry* tet_hp::getnewvrtxobject(int bnum, input_map &bdrydata) {
	hp_vrtx_bdry *temp = new hp_vrtx_bdry(*this,*vbdry(bnum));   
	gbl->vbdry_gbls(bnum) = temp->create_global_structure();
	return(temp);
}

void hp_vrtx_bdry::setvalues(init_bdry_cndtn *ibc, Array<int,1> & dirichlets, int ndirichlets) {
	
	/* UPDATE BOUNDARY CONDITION VALUES */
	int v0 = base.pnt;
	for(int n=0;n<ndirichlets;++n)
		x.ug.v(v0,dirichlets(n)) = ibc->f(dirichlets(n),x.pnts(v0),x.gbl->time);		
		
	return;
}


hp_edge_bdry* tet_hp::getnewedgeobject(int bnum, input_map &bdrydata) {
	hp_edge_bdry *temp = new hp_edge_bdry(*this,*ebdry(bnum));
	gbl->ebdry_gbls(bnum) = temp->create_global_structure();
	return(temp);
}

void hp_edge_bdry::copy(const hp_edge_bdry &bin) {
	
	if (!curved || !x.em0) return;
	
	for(int i=0;i<x.gbl->nadapt; ++i)
		crvbd(i)(Range(0,base.nseg-1),Range::all()) = bin.crvbd(i)(Range(0,base.nseg-1),Range::all());
}

void hp_edge_bdry::init(input_map& inmap,void* gbl_in) {
	int i;
	std::string keyword;
	std::istringstream data;
	std::string filename;
		if (inmap.find(base.idprefix +"_ibc") != inmap.end()) {
		ibc = x.getnewibc(base.idprefix+"_ibc",inmap);
	}
	keyword = base.idprefix + "_curved";
	inmap.getwdefault(keyword,curved,false);

	keyword = base.idprefix + "_coupled";
	inmap.getwdefault(keyword,coupled,false);
	
	if (curved && !x.coarse_level) {
		crvbd.resize(x.gbl->nhist+1);
		crv.resize(base.maxseg,x.em0);
		for(i=1;i<x.gbl->nhist+1;++i)
			crvbd(i).resize(base.maxseg,x.em0);
		crvbd(0).reference(crv);
	}
	
	dxdt.resize(x.log2pmax+1,base.maxseg);
		
	base.resize_buffers(base.maxseg*(x.em0+2)*x.NV);

	return;
}

void hp_edge_bdry::output(std::ostream& fout, tet_hp::filetype typ,int tlvl) {
	int j,m,n;

	switch(typ) {
		case(tet_hp::text):
			fout << base.idprefix << " " << mytype << std::endl;
			if (curved) {
			fout << "p0: " << x.p0 << std::endl;

			for(j=0;j<base.nseg;++j) {
				for(m=0;m<x.em0;++m) {
					for(n=0;n<tet_mesh::ND;++n)
						fout << crvbd(tlvl)(j,m)(n) << ' ';
					fout << std::endl;
				}
			}
			}
			break;
			
		case(tet_hp::binary):
			if (curved) {
				binowstream bout(&fout);
				bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::BigEndian)),1);
				bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::FloatIEEE)),1);
				bout.writeInt(x.p0,sizeof(int));

			for(j=0;j<base.nseg;++j) {
				for(m=0;m<x.em0;++m) {
					for(n=0;n<tet_mesh::ND;++n)
						bout.writeFloat(crvbd(tlvl)(j,m)(n),binio::Double);
				}
			}
			}
			break;
			
		default:
			break;
	}
	return;
}

void hp_edge_bdry::input(ifstream& fin,tet_hp::filetype typ,int tlvl) {
	int j,m,n,pmin;
	std::string idin, mytypein;   
	switch(typ) {
		case(tet_hp::text):
			fin >> idin >> mytypein;
			if (curved) { 
			fin.ignore(80,':');
			fin >>  pmin;
			pmin = x.p0;
			for(j=0;j<base.nseg;++j) {
				for(m=0;m<pmin -1;++m) {
					for(n=0;n<tet_mesh::ND;++n)
						fin >> crvbd(tlvl)(j,m)(n);
				}
				for(m=pmin-1;m<x.em0;++m) {
					for(n=0;n<tet_mesh::ND;++n)
						crvbd(tlvl)(j,m)(n) = 0.0;
				}
			}
			}
			break;
		case(tet_hp::binary):
			if (curved) {
				biniwstream bin(&fin);
				
				/* HEADER INFORMATION */
				bin.setFlag(binio::BigEndian,bin.readInt(1));
				bin.setFlag(binio::FloatIEEE,bin.readInt(1));

				pmin = bin.readInt(sizeof(int));
			for(j=0;j<base.nseg;++j) {
				for(m=0;m<pmin-1;++m) {
					for(n=0;n<tet_mesh::ND;++n)
						crvbd(tlvl)(j,m)(n) = bin.readFloat(binio::Double);
				}
	            for(m=pmin-1;m<x.em0;++m) {
					for(n=0;n<tet_mesh::ND;++n)
						crvbd(tlvl)(j,m)(n) = 0.0;
				}				
			}
			}
			break;
			
		default:
			break;
	}
}

//void hp_edge_bdry::setvalues(init_bdry_cndtn *ibc, Array<int,1>& dirichlets, int ndirichlets) {
//    int j,k,m,n,v0,v1,sind,indx,info;
//	TinyVector<FLT,tri_mesh::ND> pt;
//    char uplo[] = "U";
//	
//    /* UPDATE BOUNDARY CONDITION VALUES */
//    for(j=0;j<base.nseg;++j) {
//        sind = base.seg(j);
//        v0 = x.seg(sind).pnt(0);
//        for(n=0;n<ndirichlets;++n)
//            x.ug.v(v0,dirichlets(n)) = ibc->f(dirichlets(n),x.pnts(v0),x.gbl->time);
//    }
//    v0 = x.seg(sind).pnt(1);
//    for(n=0;n<ndirichlets;++n)
//        x.ug.v(v0,dirichlets(n)) = ibc->f(dirichlets(n),x.pnts(v0),x.gbl->time);
//	
//    /*******************/    
//    /* SET SIDE VALUES */
//    /*******************/
//    for(j=0;j<base.nseg;++j) {
//        sind = base.seg(j);
//        v0 = x.seg(sind).pnt(0);
//        v1 = x.seg(sind).pnt(1);
//        
//        if (is_curved()) {
//            x.crdtocht1d(sind);
//            for(n=0;n<tri_mesh::ND;++n)
//                basis::tri(x.log2p).proj1d(&x.cht(n,0),&x.crd(n)(0,0),&x.dcrd(n,0)(0,0));
//        }
//        else {
//            for(n=0;n<tri_mesh::ND;++n) {
//                basis::tri(x.log2p).proj1d(x.pnts(v0)(n),x.pnts(v1)(n),&x.crd(n)(0,0));
//                
//                for(k=0;k<basis::tri(x.log2p).gpx;++k)
//                    x.dcrd(n,0)(0,k) = 0.5*(x.pnts(v1)(n)-x.pnts(v0)(n));
//            }
//        }
//		
//        if (basis::tri(x.log2p).sm) {
//            for(n=0;n<ndirichlets;++n)
//                basis::tri(x.log2p).proj1d(x.ug.v(v0,dirichlets(n)),x.ug.v(v1,dirichlets(n)),&x.res(dirichlets(n))(0,0));
//			
//            for(k=0;k<basis::tri(x.log2p).gpx; ++k) {
//                pt(0) = x.crd(0)(0,k);
//                pt(1) = x.crd(1)(0,k);
//                for(n=0;n<ndirichlets;++n)
//                    x.res(dirichlets(n))(0,k) -= ibc->f(dirichlets(n),pt,x.gbl->time);
//            }
//            for(n=0;n<ndirichlets;++n)
//                basis::tri(x.log2p).intgrt1d(&x.lf(dirichlets(n))(0),&x.res(dirichlets(n))(0,0));
//			
//            indx = sind*x.sm0;
//            for(n=0;n<ndirichlets;++n) {
//                PBTRS(uplo,basis::tri(x.log2p).sm,basis::tri(x.log2p).sbwth,1,&basis::tri(x.log2p).sdiag1d(0,0),basis::tri(x.log2p).sbwth+1,&x.lf(dirichlets(n))(2),basis::tri(x.log2p).sm,info);
//                for(m=0;m<basis::tri(x.log2p).sm;++m) 
//                    x.ug.s(sind,m,dirichlets(n)) = -x.lf(dirichlets(n))(2+m);
//            }
//        }
//    }    return;
//}

void hp_edge_bdry::curv_init(int tlvl) {
	int i,j,m,n,v0,v1,sind;
	TinyVector<FLT,tet_mesh::ND> pt;
	TinyMatrix<FLT,tet_mesh::ND,MXGP> crd;
	
	j = 0;
	do {
		sind = base.seg(j).gindx;
		v0 = x.seg(sind).pnt(0);
		base.mvpttobdry(j,-1.0, x.pnts(v0));
	} while (++j < base.nseg);
	v0 = x.seg(sind).pnt(1);
	base.mvpttobdry(base.nseg-1,1.0, x.pnts(v0));
	
	if (!curved || basis::tet(x.log2p).p == 1) return;

	/*****************************/
	/* SET UP HIGHER ORDER MODES */
	/*****************************/
	for(j=0;j<base.nseg;++j) {
		sind = base.seg(j).gindx;

		v0 = x.seg(sind).pnt(0);
		v1 = x.seg(sind).pnt(1);
		
		for(n=0;n<tet_mesh::ND;++n) 
			basis::tet(x.log2p).proj1d(x.pnts(v0)(n),x.pnts(v1)(n),&crd(n,0));
	
		for(i=0;i<basis::tet(x.log2p).gpx;++i) {
			pt(0) = crd(0,i);
			pt(1) = crd(1,i);
			pt(2) = crd(2,i);
			base.mvpttobdry(j,basis::tet(x.log2p).xp(i),pt);
			crd(0,i) -= pt(0);
			crd(1,i) -= pt(1);
			crd(2,i) -= pt(2);
		}
		
		for(n=0;n<tet_mesh::ND;++n) {
			basis::tet(x.log2p).intgrt1d(&x.cf(n,0),&crd(n,0));
		
			for(m=0;m<basis::tet(x.log2p).em;++m)
				crvbd(tlvl)(j,m)(n) = -x.cf(n,m+2)*basis::tet(x.log2p).diag1d(m);
		}
	}
	return;
}

void hp_edge_bdry::findandmovebdrypt(TinyVector<FLT,tet_mesh::ND>& xp,int &bel,FLT &psi) const {
	int sind,v0,v1,iter;
	FLT dx,dy,dz,ol,roundoff,dpsi;
	TinyVector<FLT,tet_mesh::ND> pt;
		
	base.findbdrypt(xp,bel,psi);
	if (!curved) {
		base.edge_bdry::mvpttobdry(bel,psi,xp);
		basis::tet(x.log2p).ptvalues1d(psi);
		return;
	}
	
	sind = base.seg(bel).gindx;
	v0 = x.seg(sind).pnt(0);
	v1 = x.seg(sind).pnt(1);
	dx = x.pnts(v1)(0) - x.pnts(v0)(0);
	dy = x.pnts(v1)(1) - x.pnts(v0)(1);
	dz = x.pnts(v1)(2) - x.pnts(v0)(2);
	ol = 2./(dx*dx +dy*dy +dz*dz);
	dx *= ol;
	dy *= ol;
	dz *= ol;
	
	/* FIND PSI SUCH THAT TANGENTIAL POSITION ALONG LINEAR SIDE STAYS THE SAME */
	/* THIS WAY, MULTIPLE CALLS WILL NOT GIVE DIFFERENT RESULTS */ 
	x.crdtocht1d(sind);
	
	iter = 0;
	roundoff = 10.0*EPSILON*(1.0 +(fabs(xp(0)*dx) +fabs(xp(1)*dy) +fabs(xp(2)*dz)));
	do {
		basis::tet(x.log2p).ptprobe1d(x.ND,pt.data(),psi,&x.cht(0)(0),MXTM);
		dpsi = (pt(0) -xp(0))*dx +(pt(1) -xp(1))*dy +(pt(2)-xp(2))*dz;
		psi -= dpsi;
		if (iter++ > 100) {
			*x.gbl->log << "#Warning: max iterations for curved side in bdry_locate type: " << base.idnum << " seg: " << bel << " sind: " << sind << " loc: " << xp << " dpsi: " << dpsi << std::endl;
			break;
		}  
	} while (fabs(dpsi) > roundoff);
	xp = pt;
}

void hp_edge_bdry::mvpttobdry(int bel,FLT psi,TinyVector<FLT,tet_mesh::ND> &xp) {

	/* SOLUTION IS BEING ADAPTED MUST GET INFO FROM ADAPT STORAGE */
	/* FIRST GET LINEAR APPROXIMATION TO LOCATION */
	cout << "awww does nothing edge bdry mvpttobdry " << endl;
	//base.edge_bdry::mvpttobdry(bel, psi, xp);
	//adapt_storage->findandmovebdrypt(xp,bel,psi);
	
	return;
}


void hp_edge_bdry::tadvance() {
	int stage = x.gbl->substep +x.gbl->esdirk;  
		
	if (x.p0 > 1 && curved) {
		if (stage) {
			/* BACK CALCULATE K TERM */
			for(int j=0;j<base.nseg;++j) {
			for(int m=0;m<basis::tet(x.log2p).em;++m) {
				for(int n=0;n<tet_mesh::ND;++n)
					crvbd(stage+1)(j,m)(n) = (crvbd(0)(j,m)(n)-crvbd(1)(j,m)(n))*x.gbl->adirk(stage-1,stage-1);
			}
			}
		}
		
		if (x.gbl->substep == 0) {
			/* STORE TILDE W */
			for(int j=0;j<base.nseg;++j) {
			for(int m=0;m<basis::tet(x.log2p).em;++m) {
				for(int n=0;n<tet_mesh::ND;++n)
					crvbd(1)(j,m)(n) = crv(j,m)(n);
			}
			}
		}
		
		/* UPDATE TILDE W */
		for (int s=0;s<stage;++s) {         
			for(int j=0;j<base.nseg;++j) {
			for(int m=0;m<basis::tet(x.log2p).em;++m) {
				for(int n=0;n<tet_mesh::ND;++n) {
					crvbd(1)(j,m)(n) += x.gbl->adirk(stage,s)*crvbd(s+2)(j,m)(n);
				}
			}
			}
		}

		/* EXTRAPOLATE GUESS? */
		if (stage && x.gbl->dti > 0.0) {
			FLT constant =  x.gbl->cdirk(x.gbl->substep);
			crvbd(0)(Range(0,base.nseg-1),Range(0,basis::tet(x.log2p).em-1)) += constant*crvbd(stage+1)(Range(0,base.nseg-1),Range(0,basis::tet(x.log2p).em-1));
		}
	}

	calculate_unsteady_sources();
		
	/* EXTRAPOLATE SOLUTION OR MOVE TO NEXT TIME POSITION */
	if (!coupled && curved) curv_init();
		
	return;
}

void hp_edge_bdry::calculate_unsteady_sources() {
	int i,j,n,sind;
	
	for(i=0;i<=x.log2pmax;++i) {
		for(j=0;j<base.nseg;++j) {
			sind = base.seg(j).gindx;
			x.crdtocht1d(sind,1);
			for(n=0;n<tet_mesh::ND;++n)
			basis::tet(i).proj1d(&x.cht(n)(0),&dxdt(i,j)(n,0));
		}
	}
	
	return;
}


void hp_edge_bdry::setvalues(init_bdry_cndtn *ibc, Array<int,1> & dirichlets, int ndirichlets) {
	int j,k,m,n,v0,v1,sind;
	TinyVector<FLT,tet_mesh::ND> pt;
	
	/* UPDATE BOUNDARY CONDITION VALUES */
	for(j=0;j<base.nseg;++j) {
		sind = base.seg(j).gindx;
		v0 = x.seg(sind).pnt(0);
		for(n=0;n<ndirichlets;++n)
			x.ug.v(v0,dirichlets(n)) = ibc->f(dirichlets(n),x.pnts(v0),x.gbl->time);		
	}
	v0 = x.seg(sind).pnt(1);
	for(n=0;n<ndirichlets;++n)
		x.ug.v(v0,dirichlets(n)) = ibc->f(dirichlets(n),x.pnts(v0),x.gbl->time);
	
	/*******************/    
	/* SET SIDE VALUES */
	/*******************/
	for(j=0;j<base.nseg;++j) {
		sind = base.seg(j).gindx;
		v0 = x.seg(sind).pnt(0);
		v1 = x.seg(sind).pnt(1);
		
		if (is_curved()) {
			x.crdtocht1d(sind);
			for(n=0;n<tet_mesh::ND;++n)
				basis::tet(x.log2p).proj1d(&x.cht(n)(0),&x.crd1d(n)(0));
		}
		else {
			for(n=0;n<tet_mesh::ND;++n) {
				basis::tet(x.log2p).proj1d(x.pnts(v0)(n),x.pnts(v1)(n),&x.crd1d(n)(0));
			}
		}
		if (basis::tet(x.log2p).em) {
			for(n=0;n<ndirichlets;++n)
				basis::tet(x.log2p).proj1d(x.ug.v(v0,dirichlets(n)),x.ug.v(v1,dirichlets(n)),&x.res1d(dirichlets(n))(0));
			
			for(k=0;k<basis::tet(x.log2p).gpx; ++k) {
				pt(0) = x.crd1d(0)(k);
				pt(1) = x.crd1d(1)(k);
				pt(2) = x.crd1d(2)(k);
				
				for(n=0;n<ndirichlets;++n)
					x.res1d(dirichlets(n))(k) -= ibc->f(dirichlets(n),pt,x.gbl->time);
			}
			for(n=0;n<ndirichlets;++n)
				basis::tet(x.log2p).intgrt1d(&x.lf(dirichlets(n))(0),&x.res1d(dirichlets(n))(0));
			
			for(n=0;n<ndirichlets;++n) 
				for(m=0;m<basis::tet(x.log2p).em;++m) 
					x.ug.e(sind,m,dirichlets(n)) = -x.lf(dirichlets(n))(2+m)*basis::tet(x.log2p).diag1d(m);
			
		}
	}
	
	return;
}


hp_face_bdry* tet_hp::getnewfaceobject(int bnum, input_map &bdrydata) {
	hp_face_bdry *temp = new hp_face_bdry(*this,*fbdry(bnum));
	gbl->fbdry_gbls(bnum) = temp->create_global_structure();
	return(temp);
}

void hp_face_bdry::copy(const hp_face_bdry &bin) {
	
	if (!curved || !x.em0) return;
	
	for(int i=0;i<x.gbl->nadapt; ++i) {
		ecrvbd(i)(Range(0,base.nseg-1),Range::all()) = bin.ecrvbd(i)(Range(0,base.nseg-1),Range::all());
		fcrvbd(i)(Range(0,base.nseg-1),Range::all()) = bin.fcrvbd(i)(Range(0,base.nseg-1),Range::all());
	}
}

void hp_face_bdry::init(input_map& inmap,void* gbl_in) {
	int i;
	std::string keyword;
	std::istringstream data;
	std::string filename;
	
	if (inmap.find(base.idprefix +"_ibc") != inmap.end()) {
		ibc = x.getnewibc(base.idprefix+"_ibc",inmap);
	}
	
	keyword = base.idprefix + "_curved";
	inmap.getwdefault(keyword,curved,false);

	keyword = base.idprefix + "_coupled";
	inmap.getwdefault(keyword,coupled,false);
	
	if (curved && !x.coarse_level) {
		ecrvbd.resize(x.gbl->nhist+1);
		ecrv.resize(base.maxpst,x.em0);
		for(i=1;i<x.gbl->nhist+1;++i)
			ecrvbd(i).resize(base.maxpst,x.em0);
		ecrvbd(0).reference(ecrv);
		
		fcrvbd.resize(x.gbl->nhist+1);
		fcrv.resize(base.maxpst,x.fm0);
		for(i=1;i<x.gbl->nhist+1;++i)
			fcrvbd(i).resize(base.maxpst,x.fm0);
		fcrvbd(0).reference(fcrv);
	}
	
	dxdt.resize(x.log2pmax+1,base.maxpst);
		
	base.resize_buffers(base.maxpst*(x.fm0+2)*x.NV);

	return;
}

void hp_face_bdry::output(std::ostream& fout, tet_hp::filetype typ,int tlvl) {
	int j,m,n;

	switch(typ) {
		case(tet_hp::text):
			fout << base.idprefix << " " << mytype << std::endl;
			if (curved) {
			fout << "p0: " << x.p0 << std::endl;

			for(j=0;j<base.nseg;++j) {
				for(m=0;m<x.em0;++m) {
					for(n=0;n<tet_mesh::ND;++n)
						fout << ecrvbd(tlvl)(j,m)(n) << ' ';
					fout << std::endl;
				}
			}

			for(j=0;j<base.ntri;++j) {
				for(m=0;m<x.fm0;++m) {
					for(n=0;n<tet_mesh::ND;++n)
						fout << fcrvbd(tlvl)(j,m)(n) << ' ';
					fout << std::endl;
				}
			}                       
			}
			break;
			
		case(tet_hp::binary):
			if (curved) {
				binowstream bout(&fout);
				bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::BigEndian)),1);
				bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::FloatIEEE)),1);
				bout.writeInt(x.p0,sizeof(int));

				for(j=0;j<base.nseg;++j) {
					for(m=0;m<x.em0;++m) {
						for(n=0;n<tet_mesh::ND;++n)
							bout.writeFloat(ecrvbd(tlvl)(j,m)(n),binio::Double);
					}
				}
				for(j=0;j<base.ntri;++j) {
					for(m=0;m<x.fm0;++m) {
						for(n=0;n<tet_mesh::ND;++n)
							bout.writeFloat(fcrvbd(tlvl)(j,m)(n),binio::Double);
					}
				}
			}
			break;
			
		default:
			break;
	}
	return;
}

void hp_face_bdry::input(ifstream& fin,tet_hp::filetype typ,int tlvl) {
	int j,m,n,pmin;
	std::string idin, mytypein;   
	switch(typ) {
		case(tet_hp::text):
			fin >> idin >> mytypein;
			if (curved) { 
			fin.ignore(80,':');
			fin >>  pmin;
			pmin = x.p0;
			for(j=0;j<base.nseg;++j) {
				for(m=0;m<pmin -1;++m) {
					for(n=0;n<tet_mesh::ND;++n)
						fin >> ecrvbd(tlvl)(j,m)(n);
				}
				for(m=pmin-1;m<x.em0;++m) {
					for(n=0;n<tet_mesh::ND;++n)
						ecrvbd(tlvl)(j,m)(n) = 0.0;
				}
			}

			// FIXME: DOESN'T WORK FOR CHANGING ORDERS
			for(j=0;j<base.ntri;++j) {
				for(m=0;m<x.fm0;++m) {
					for(n=0;n<tet_mesh::ND;++n)
						fin >> fcrvbd(tlvl)(j,m)(n);
				}
			}


			}
			break;
		case(tet_hp::binary):
			if (curved) {
				biniwstream bin(&fin);
				
				/* HEADER INFORMATION */
				bin.setFlag(binio::BigEndian,bin.readInt(1));
				bin.setFlag(binio::FloatIEEE,bin.readInt(1));

				pmin = bin.readInt(sizeof(int));
			for(j=0;j<base.nseg;++j) {
				for(m=0;m<pmin-1;++m) {
					for(n=0;n<tet_mesh::ND;++n)
						ecrvbd(tlvl)(j,m)(n) = bin.readFloat(binio::Double);
				}
	            for(m=pmin-1;m<x.em0;++m) {
					for(n=0;n<tet_mesh::ND;++n)
						ecrvbd(tlvl)(j,m)(n) = 0.0;
				}				
			}

			// FIXME: DOESN'T WORK FOR CHANGING ORDERS
			for(j=0;j<base.ntri;++j) {
				for(m=0;m<x.fm0;++m) {
					for(n=0;n<tet_mesh::ND;++n)
						fcrvbd(tlvl)(j,m)(n) = bin.readFloat(binio::Double);
				}
			}
			}
			break;
			
		default:
			break;
	}
}

void hp_face_bdry::setvalues(init_bdry_cndtn *ibc, Array<int,1> & dirichlets, int ndirichlets) {
	int i,j,k,m,n,v0,v1,v2,sind,find;
	TinyVector<FLT,tet_mesh::ND> pt;
	
	/* UPDATE BOUNDARY CONDITION VALUES */
	for(j=0;j<base.npnt;++j) {
		v0 = base.pnt(j).gindx;
		for(n=0;n<ndirichlets;++n)
			x.ug.v(v0,dirichlets(n)) = ibc->f(dirichlets(n),x.pnts(v0),x.gbl->time);		
	}
	
	/*******************/    
	/* SET SIDE VALUES */
	/*******************/
	for(j=0;j<base.nseg;++j) {
		sind = base.seg(j).gindx;
		v0 = x.seg(sind).pnt(0);
		v1 = x.seg(sind).pnt(1);
		
		if (is_curved()) {
			x.crdtocht1d(sind);
			for(n=0;n<tet_mesh::ND;++n)
				basis::tet(x.log2p).proj1d(&x.cht(n)(0),&x.crd1d(n)(0));
		}
		else {
			for(n=0;n<tet_mesh::ND;++n) {
				basis::tet(x.log2p).proj1d(x.pnts(v0)(n),x.pnts(v1)(n),&x.crd1d(n)(0));
			}
		}
		if (basis::tet(x.log2p).em) {
			for(n=0;n<ndirichlets;++n)
				basis::tet(x.log2p).proj1d(x.ug.v(v0,dirichlets(n)),x.ug.v(v1,dirichlets(n)),&x.res1d(dirichlets(n))(0));
			
			for(k=0;k<basis::tet(x.log2p).gpx; ++k) {
				pt(0) = x.crd1d(0)(k);
				pt(1) = x.crd1d(1)(k);
				pt(2) = x.crd1d(2)(k);
				
				for(n=0;n<ndirichlets;++n)
					x.res1d(dirichlets(n))(k) -= ibc->f(dirichlets(n),pt,x.gbl->time);
			}
			for(n=0;n<ndirichlets;++n)
				basis::tet(x.log2p).intgrt1d(&x.lf(dirichlets(n))(0),&x.res1d(dirichlets(n))(0));
			
			for(n=0;n<ndirichlets;++n) 
				for(m=0;m<basis::tet(x.log2p).em;++m) 
					x.ug.e(sind,m,dirichlets(n)) = -x.lf(dirichlets(n))(2+m)*basis::tet(x.log2p).diag1d(m);
			
		}
	}
	
	/*******************/    
	/* SET FACE VALUES */
	/*******************/
	for(j=0;j<base.ntri;++j) {
		find = base.tri(j).gindx;
		v0 = x.tri(find).pnt(0);
		v1 = x.tri(find).pnt(1);
		v2 = x.tri(find).pnt(2);
		
		if (is_curved()) {
			x.crdtocht2d(find);
			for(n=0;n<tet_mesh::ND;++n)
				basis::tet(x.log2p).proj2d_bdry(&x.cht(n)(0),&x.crd2d(n)(0)(0),MXGP);
		}
		else {
			for(n=0;n<tet_mesh::ND;++n) 
				basis::tet(x.log2p).proj2d(x.pnts(v0)(n),x.pnts(v1)(n),x.pnts(v2)(n),&x.crd2d(n)(0)(0),MXGP);                

				//basis::tet(x.log2p).proj2d(x.vrtxbd(0)(x.tri(find).pnt(0))(n),x.vrtxbd(0)(x.tri(find).pnt(1))(n),x.vrtxbd(0)(x.tri(find).pnt(2))(n),&x.crd2d(n)(0)(0),MXGP);

		}
		if (basis::tet(x.log2p).fm) {
//			for(n=0;n<ndirichlets;++n)
//				basis::tet(x.log2p).proj2d(x.ug.v(v0,dirichlets(n)),x.ug.v(v1,dirichlets(n)),x.ug.v(v2,dirichlets(n)),&x.res2d(dirichlets(n))(0)(0),MXGP);
//			cout << " there may be bug in face_bdry::setvalues for high order"<< endl;
			x.ugtouht2d_bdry(find);//temp may need to fix this for high order
			for(n=0;n<ndirichlets;++n)
				basis::tet(x.log2p).proj2d_bdry(&x.uht(dirichlets(n))(0),&x.res2d(dirichlets(n))(0)(0),MXGP);
			
			for(i=0;i<basis::tet(x.log2p).gpx; ++i) {
				for(k=0;k<basis::tet(x.log2p).gpy; ++k) {					
					pt(0) = x.crd2d(0)(i)(k);
					pt(1) = x.crd2d(1)(i)(k);
					pt(2) = x.crd2d(2)(i)(k);
					
					for(n=0;n<ndirichlets;++n)
						x.res2d(dirichlets(n))(i)(k) -= ibc->f(dirichlets(n),pt,x.gbl->time);
				}
			}
			for(n=0;n<ndirichlets;++n)
				basis::tet(x.log2p).intgrt2d(&x.lf(dirichlets(n))(0),&x.res2d(dirichlets(n))(0)(0),MXGP);
			
			for(n=0;n<ndirichlets;++n) {
				for(m=0;m<basis::tet(x.log2p).fm;++m) 
					x.ug.f(find,m,dirichlets(n)) = -x.lf(dirichlets(n))(3+3*basis::tet(x.log2p).em+m)*basis::tet(x.log2p).diag2d(m);
			}
		}
	}
	
	
	return;
}

void hp_face_bdry::curv_init(int tlvl) {
	int i,j,m,n,v0,v1,sind,info;
	TinyVector<FLT,tet_mesh::ND> pt;
	char uplo[] = "U";
	TinyMatrix<FLT,tet_mesh::ND,MXGP> crd;
	int stridey = MXGP;

	j = 0;
	do {
		sind = base.seg(j).gindx;
		v0 = x.seg(sind).pnt(0);
		base.mvpttobdry(j,-1.0, x.pnts(v0));
	} while (++j < base.nseg);
	v0 = x.seg(sind).pnt(1);
	base.mvpttobdry(base.nseg-1,1.0, x.pnts(v0));
	
	if (!curved || basis::tet(x.log2p).p == 1) return;
	
	
//	
//	/* SKIP END VERTICES */
//	for(j=1;j<base.npnt;++j) {
//		v0 = base.pnt(j).gindx;
//		base.mvpttobdry(base.pnt(j).tri,-1.0,-1.0,x.pnts(v0));  // FIXME WRONG
//	}
//
//	if (basis::tet(x.log2p).p == 1) return;
		
	/*****************************/
	/* SET UP HIGHER ORDER MODES */
	/*****************************/
	for(j=0;j<base.nseg;++j) {
		sind = base.seg(j).gindx;

		v0 = x.seg(sind).pnt(0);
		v1 = x.seg(sind).pnt(1);
		
		for(n=0;n<tet_mesh::ND;++n) 
			basis::tet(x.log2p).proj1d(x.pnts(v0)(n),x.pnts(v1)(n),&crd(n,0));
	
		for(i=0;i<basis::tet(x.log2p).gpx;++i) {
			pt(0) = crd(0,i);
			pt(1) = crd(1,i);
			pt(2) = crd(2,i);
			//base.mvpttobdry(base.seg(j).tri(0),-1.0,basis::tet(x.log2p).xp(i),pt); // FIXME
			base.mvpttobdry(base.seg(j).gindx,basis::tet(x.log2p).xp(i),pt); // FIXME
			crd(0,i) -= pt(0);
			crd(1,i) -= pt(1);
			crd(2,i) -= pt(2);
		}
		
		for(n=0;n<tet_mesh::ND;++n) {
			basis::tet(x.log2p).intgrt1d(&x.cf(n,0),&crd(n,0));
		
			for(m=0;m<basis::tet(x.log2p).em;++m)
			ecrvbd(tlvl)(j,m)(n) = -x.cf(n,m+2)*basis::tet(x.log2p).diag1d(m);
		}
	}
	
	if (basis::tet(x.log2p).fm <= 0) return;  
	
	for(j=0;j<base.ntri;++j) {
		int tind = base.tri(j).gindx;
		// x.crdtocht2d_bdry(tind,tlvl);  FIXME
		
		for(n=0;n<tet_mesh::ND;++n)
			basis::tet(x.log2p).proj2d_bdry(&x.cht(n)(0),&x.crd2d(n)(0)(0),stridey);
		
		for (i=0; i < basis::tet(x.log2p).gpx; ++i ) {
			for (j=0; j < basis::tet(x.log2p).gpy; ++j ) {
			pt(0) = x.crd2d(0)(i)(j);
			pt(1) = x.crd2d(1)(i)(j);
			pt(2) = x.crd2d(2)(i)(j);
			base.mvpttobdry(j,basis::tet(x.log2p).xp(i),basis::tet(x.log2p).yp(j),pt);
			x.crd2d(0)(i)(j) -= pt(0);
			x.crd2d(1)(i)(j) -= pt(1);
			x.crd2d(2)(i)(j) -= pt(2);
			}
		}
		
		for(n=0;n<tet_mesh::ND;++n) {
			basis::tet(x.log2p).intgrt2d(&x.lf(n)(0),&x.crd2d(n)(0)(0),stridey);
			for(i=0;i<basis::tet(x.log2p).fm;++i){
			fcrvbd(tlvl)(j,m)(n) = -x.lf(n)(3+3*basis::tet(x.log2p).em+i)*basis::tet(x.log2p).diag2d(i);
			}
		}
	}
	
	return;
}

void hp_face_bdry::tmatchsolution_snd(FLT *fdata, int bgnmode, int endmode, int modestride) {
	if (base.is_frst()) {
		base.tloadbuff(boundary::all,fdata,bgnmode*x.NV,(endmode+1)*x.NV-1,x.NV*modestride);
	}
	else {
		if (!base.in_group(boundary::all)) return;

		int count = 0;
		for (int j=0;j<base.ntri;++j) {
			int tind = base.tri(j).gindx;
			int offset = tind*modestride*x.NV;
			
			int msgn = 1;
			int indx = 0;
			for(int m = 1; m <= basis::tet(x.log2p).em-1; ++m) {
			for(int k = 1; k <= basis::tet(x.log2p).em-m; ++k) {
				if (indx >= bgnmode && indx <= endmode) {
					for(int n = 0; n < x.NV; ++n)
						base.fsndbuf(count++) = msgn*fdata[offset++];
				}
				++indx;
			}
			msgn *= -1;
			}
		}  
		
		base.sndsize() = count;
		base.sndtype() = boundary::flt_msg;
	}
	return;
}


void hp_face_bdry::tmatchsolution_rcv(FLT *fdata, int bgnmode, int endmode, int modestride) {
	if (base.is_frst()) {
		base.tfinalrcv(boundary::all,0,boundary::symmetric,boundary::average,fdata,bgnmode*x.NV,(endmode+1)*x.NV-1,x.NV*modestride);
	}
	else {
		bool reload = base.comm_finish(boundary::all,0,boundary::symmetric,boundary::average);
		if (!reload) return;
	
		int count = 0;
		for (int j=0;j<base.ntri;++j) {
			int tind = base.tri(j).gindx;
			int offset = tind*modestride*x.NV;
			
			int msgn = 1;
			int indx = 0;
			for(int m = 1; m <= basis::tet(x.log2p).em-1; ++m) {
			for(int k = 1; k <= basis::tet(x.log2p).em-m; ++k) {
				if (indx >= bgnmode && indx <= endmode) {
					for(int n = 0; n < x.NV; ++n)
						fdata[offset++] = msgn*base.fsndbuf(count++);
				}
				++indx;
			}
			msgn *= -1;
			}
		}    
	}
	return;
}


//
//void hp_face_bdry::findandmovebdrypt(TinyVector<FLT,tet_mesh::ND>& xp,int &bel,FLT &psi) const {
//   int sind,v0,v1,iter;
//   FLT dx,dy,dz,ol,roundoff,dpsi;
//   TinyVector<FLT,tet_mesh::ND> pt;
//      
//   base.findbdrypt(xp,bel,psi);
//   if (!curved) {
//      base.face_bdry::mvpttobdry(bel,psi,xp);
//      basis::tet(x.log2p).ptvalues1d(psi);
//      return;
//   }
//   
//   sind = base.seg(bel);
//   v0 = x.seg(sind).pnt(0);
//   v1 = x.seg(sind).pnt(1);
//   dx = x.pnts(v1)(0) - x.pnts(v0)(0);
//   dy = x.pnts(v1)(1) - x.pnts(v0)(1);
//   dz = x.pnts(v1)(2) - x.pnts(v0)(2);
//   ol = 2./(dx*dx +dy*dy +dz*dz);
//   dx *= ol;
//   dy *= ol;
//   dz *= ol;
//   
//   /* FIND PSI SUCH THAT TANGENTIAL POSITION ALONG LINEAR SIDE STAYS THE SAME */
//   /* THIS WAY, MULTIPLE CALLS WILL NOT GIVE DIFFERENT RESULTS */ 
//   x.crdtocht1d(sind);
//   
//   iter = 0;
//   roundoff = 10.0*EPSILON*(1.0 +(fabs(xp(0)*dx) +fabs(xp(1)*dy) +fabs(xp(2)*dz)));
//   do {
//      basis::tet(x.log2p).ptprobe1d(x.ND,pt.data(),psi,&x.cht(0)(0),MXTM);
//      dpsi = (pt(0) -xp(0))*dx +(pt(1) -xp(1))*dy +(pt(2)-xp(2))*dz;
//      psi -= dpsi;
//      if (iter++ > 100) {
//         *x.gbl->log << "#Warning: max iterations for curved side in bdry_locate type: " << base.idnum << " seg: " << bel << " sind: " << sind << " loc: " << xp << " dpsi: " << dpsi << std::endl;
//         break;
//      }  
//   } while (fabs(dpsi) > roundoff);
//   xp = pt;
//}

void hp_face_bdry::mvpttobdry(int bel,FLT psi,TinyVector<FLT,tet_mesh::ND> &xp) {

   /* SOLUTION IS BEING ADAPTED MUST GET INFO FROM ADAPT STORAGE */
   /* FIRST GET LINEAR APPROXIMATION TO LOCATION */
	cout << "awww does nothing face boundary mvpttobdry 1" << endl;
  // base.face_bdry::mvpttobdry(bel, psi, xp);
   //adapt_storage->findandmovebdrypt(xp,bel,psi);
   
   return;
}

void hp_face_bdry::mvpttobdry(int bel, FLT r, FLT s, TinyVector<FLT, tet_mesh::ND> &xp) {
	
	/* SOLUTION IS BEING ADAPTED MUST GET INFO FROM ADAPT STORAGE */
	/* FIRST GET LINEAR APPROXIMATION TO LOCATION */
	cout << "awww does nothing face boundary mvpttobdry 2" << endl;
	// base.face_bdry::mvpttobdry(bel, psi, xp);
	//adapt_storage->findandmovebdrypt(xp,bel,psi);
	
	return;
}

void hp_face_bdry::tadvance() {
	int stage = x.gbl->substep +x.gbl->esdirk;  
		
	if (x.p0 > 1 && curved) {
		if (stage) {
			/* BACK CALCULATE K TERM */
			for(int j=0;j<base.nseg;++j) {
			for(int m=0;m<basis::tet(x.log2p).em;++m) {
				for(int n=0;n<tet_mesh::ND;++n)
					ecrvbd(stage+1)(j,m)(n) = (ecrvbd(0)(j,m)(n)-ecrvbd(1)(j,m)(n))*x.gbl->adirk(stage-1,stage-1);
			}
			}
			
			/* BACK CALCULATE K TERM */
			for(int j=0;j<base.ntri;++j) {
			for(int m=0;m<basis::tet(x.log2p).fm;++m) {
				for(int n=0;n<tet_mesh::ND;++n)
					fcrvbd(stage+1)(j,m)(n) = (fcrvbd(0)(j,m)(n)-fcrvbd(1)(j,m)(n))*x.gbl->adirk(stage-1,stage-1);
			}
			}
		}
		
		if (x.gbl->substep == 0) {
			/* STORE TILDE W */
			for(int j=0;j<base.nseg;++j) {
			for(int m=0;m<basis::tet(x.log2p).em;++m) {
				for(int n=0;n<tet_mesh::ND;++n)
					ecrvbd(1)(j,m)(n) = ecrv(j,m)(n);
			}
			}
			
			for(int j=0;j<base.ntri;++j) {
			for(int m=0;m<basis::tet(x.log2p).fm;++m) {
				for(int n=0;n<tet_mesh::ND;++n)
					fcrvbd(1)(j,m)(n) = fcrv(j,m)(n);
			}
			}
		}
		
		/* UPDATE TILDE W */
		for (int s=0;s<stage;++s) {         
			for(int j=0;j<base.nseg;++j) {
			for(int m=0;m<basis::tet(x.log2p).em;++m) {
				for(int n=0;n<tet_mesh::ND;++n) {
					ecrvbd(1)(j,m)(n) += x.gbl->adirk(stage,s)*ecrvbd(s+2)(j,m)(n);
				}
			}
			}
			
			for(int j=0;j<base.ntri;++j) {
			for(int m=0;m<basis::tet(x.log2p).fm;++m) {
				for(int n=0;n<tet_mesh::ND;++n) {
					fcrvbd(1)(j,m)(n) += x.gbl->adirk(stage,s)*fcrvbd(s+2)(j,m)(n);
				}
			}
			}
		}

		/* EXTRAPOLATE GUESS? */
		if (stage && x.gbl->dti > 0.0) {
			FLT constant =  x.gbl->cdirk(x.gbl->substep);
			ecrvbd(0)(Range(0,base.nseg-1),Range(0,basis::tet(x.log2p).em-1)) += constant*ecrvbd(stage+1)(Range(0,base.nseg-1),Range(0,basis::tet(x.log2p).em-1));
			fcrvbd(0)(Range(0,base.ntri-1),Range(0,basis::tet(x.log2p).fm-1)) += constant*fcrvbd(stage+1)(Range(0,base.ntri-1),Range(0,basis::tet(x.log2p).fm-1));

		}
	}

	calculate_unsteady_sources();
		
	/* EXTRAPOLATE SOLUTION OR MOVE TO NEXT TIME POSITION */
	if (!coupled && curved) curv_init();
		
	return;
}

void hp_face_bdry::calculate_unsteady_sources() {
	int i,j,n;
	
	for(i=0;i<=x.log2pmax;++i) {
		for(j=0;j<base.ntri;++j) {
			int tind = base.tri(j).gindx;
			x.crdtocht2d(tind,1);
		
			for(n=0;n<tet_mesh::ND;++n)
				basis::tet(x.log2p).proj2d(&x.cht(n)(0),&dxdt(i,j)(n)(0)(0),MXGP);
		}
	}
	
	return;
}


void tet_hp::pc0load(int phase, FLT *pdata, int vrtstride) {
	int i;
			
	/* SEND COMMUNICATIONS TO ADJACENT MESHES */
	for(i=0;i<nfbd;++i) 
		hp_fbdry(i)->pmatchsolution_snd(phase,pdata,vrtstride);
	for(i=0;i<nebd;++i) 
		hp_ebdry(i)->pmatchsolution_snd(phase,pdata,vrtstride);
	for(i=0;i<nvbd;++i)
		hp_vbdry(i)->pmatchsolution_snd(phase,pdata,vrtstride);

	for(i=0;i<nfbd;++i)
		fbdry(i)->comm_prepare(boundary::all_phased,phase,boundary::symmetric);
	for(i=0;i<nebd;++i)
		ebdry(i)->comm_prepare(boundary::all_phased,phase,boundary::symmetric);
	for(i=0;i<nvbd;++i)
		vbdry(i)->comm_prepare(boundary::all_phased,phase,boundary::symmetric);
	
	return;
}

int tet_hp::pc0wait_rcv(int phase, FLT *pdata, int vrtstride) {
	int stop = 1;
	int i;

	for(i=0;i<nfbd;++i) {
		stop &= fbdry(i)->comm_wait(boundary::all_phased,phase,boundary::symmetric);
	}
	for(i=0;i<nebd;++i) {
		stop &= ebdry(i)->comm_wait(boundary::all_phased,phase,boundary::symmetric);
	}
	for(i=0;i<nvbd;++i) {
		stop &= vbdry(i)->comm_wait(boundary::all_phased,phase,boundary::symmetric);
	}
		
	for(i=0;i<nfbd;++i) 
		hp_fbdry(i)->pmatchsolution_rcv(phase,pdata,vrtstride);
	for(i=0;i<nebd;++i) 
		hp_ebdry(i)->pmatchsolution_rcv(phase,pdata,vrtstride);
	for(i=0;i<nvbd;++i)
		hp_vbdry(i)->pmatchsolution_rcv(phase,pdata,vrtstride);
		
	return(stop);
}

int tet_hp::pc0rcv(int phase, FLT *pdata, int vrtstride) {
	int stop = 1,i;

	for(i=0;i<nfbd;++i)
		stop &= fbdry(i)->comm_nowait(boundary::all_phased,phase,boundary::symmetric);
	for(i=0;i<nebd;++i)
		stop &= ebdry(i)->comm_nowait(boundary::all_phased,phase,boundary::symmetric);
	for(i=0;i<nvbd;++i)
		stop &= vbdry(i)->comm_nowait(boundary::all_phased,phase,boundary::symmetric);

	for(i=0;i<nfbd;++i) 
		hp_fbdry(i)->pmatchsolution_rcv(phase,pdata,vrtstride);      
	for(i=0;i<nebd;++i) 
		hp_ebdry(i)->pmatchsolution_rcv(phase,pdata,vrtstride);
	for(i=0;i<nvbd;++i)
		hp_vbdry(i)->pmatchsolution_rcv(phase,pdata,vrtstride);
		
	return(stop);
}

void tet_hp::sc0load(int phase, FLT *sdata, int bgn, int end, int stride) {
	int i;
	
	/* SEND COMMUNICATIONS TO ADJACENT MESHES */
	for(i=0;i<nfbd;++i) 
		hp_fbdry(i)->smatchsolution_snd(phase,sdata,bgn,end,stride);
	for(i=0;i<nebd;++i) 
		hp_ebdry(i)->smatchsolution_snd(phase,sdata,bgn,end,stride);

	for(i=0;i<nfbd;++i)
		fbdry(i)->comm_prepare(boundary::all_phased,phase,boundary::symmetric);
	for(i=0;i<nebd;++i)
		ebdry(i)->comm_prepare(boundary::all_phased,phase,boundary::symmetric);
	
	return;
}

int tet_hp::sc0wait_rcv(int phase, FLT *sdata, int bgn, int end, int stride) {
	int stop = 1;
	int i;

	for(i=0;i<nfbd;++i) {
		stop &= fbdry(i)->comm_wait(boundary::all_phased,phase,boundary::symmetric);
	}
	for(i=0;i<nebd;++i) {
		stop &= ebdry(i)->comm_wait(boundary::all_phased,phase,boundary::symmetric);
	}
		
	for(i=0;i<nfbd;++i) 
		hp_fbdry(i)->smatchsolution_rcv(phase,sdata,bgn,end,stride);
	for(i=0;i<nebd;++i) 
		hp_ebdry(i)->smatchsolution_rcv(phase,sdata,bgn,end,stride);
		
	return(stop);
}

int tet_hp::sc0rcv(int phase, FLT *sdata, int bgn, int end, int stride) {
	int stop = 1,i;

	for(i=0;i<nfbd;++i)
		stop &= fbdry(i)->comm_nowait(boundary::all_phased,phase,boundary::symmetric);
	for(i=0;i<nebd;++i)
		stop &= ebdry(i)->comm_nowait(boundary::all_phased,phase,boundary::symmetric);

	for(i=0;i<nfbd;++i) 
		hp_fbdry(i)->smatchsolution_rcv(phase,sdata,bgn,end,stride);      
	for(i=0;i<nebd;++i) 
		hp_ebdry(i)->smatchsolution_rcv(phase,sdata,bgn,end,stride);
		
	return(stop);
}   

void tet_hp::tc0load(FLT *tdata, int bgn, int end, int stride) {
	int i;
			
	/* SEND COMMUNICATIONS TO ADJACENT MESHES */
	for(i=0;i<nfbd;++i) 
		hp_fbdry(i)->tmatchsolution_snd(tdata,bgn,end,stride);

	for(i=0;i<nfbd;++i)
		fbdry(i)->comm_prepare(boundary::all,0,boundary::symmetric);
	
	return;
}

int tet_hp::tc0wait_rcv(FLT *tdata, int bgn, int end, int stride) {
	int stop = 1;
	int i;

	for(i=0;i<nfbd;++i) {
		stop &= fbdry(i)->comm_wait(boundary::all,0,boundary::symmetric);
	}
		
	for(i=0;i<nfbd;++i) 
		hp_fbdry(i)->tmatchsolution_rcv(tdata,bgn,end,stride);
		
	return(stop);
}

int tet_hp::tc0rcv(FLT *tdata, int bgn, int end, int stride) {
	int stop = 1,i;

	for(i=0;i<nfbd;++i)
		stop &= fbdry(i)->comm_nowait(boundary::all,0,boundary::symmetric);

	for(i=0;i<nfbd;++i) 
		hp_fbdry(i)->tmatchsolution_rcv(tdata,bgn,end,stride);      
		
	return(stop);
}  


void tet_hp::matchboundaries() {
	int i, m, n, msgn, bnum, count;
	int last_phase, mp_phase;
		
	/* Match boundary vertices */
	tet_mesh::matchboundaries();
			
	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		pc0load(mp_phase,ug.v.data());
		pmsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
		last_phase = true;
		last_phase &= pc0wait_rcv(mp_phase,ug.v.data());
	}
	
	if (!em0) return;
	
	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		sc0load(mp_phase,ug.e.data(),0,em0-1,ug.e.extent(secondDim));
		smsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
		last_phase = true;
		last_phase &= sc0wait_rcv(mp_phase,ug.e.data(),0,em0-1,ug.e.extent(secondDim));      
	}	

	/* Match curved sides */
	for(bnum=0;bnum<nebd;++bnum) {
		if (ebdry(bnum)->is_comm() && hp_ebdry(bnum)->is_curved()) {            
			count = 0;
			for(i=0;i<ebdry(bnum)->nseg;++i) {
			for(m=0;m<basis::tet(log2p).em;++m) {
				for(n=0;n<ND;++n)
					ebdry(bnum)->fsndbuf(count++) = hp_ebdry(bnum)->crde(i,m,n);
			}
			}
			ebdry(bnum)->sndsize() = count;
			ebdry(bnum)->sndtype() = boundary::flt_msg;
			ebdry(bnum)->comm_prepare(boundary::all,0,boundary::master_slave);
		}
	}
	
	for(bnum=0;bnum<nebd;++bnum) {
		if (ebdry(bnum)->is_comm() && hp_ebdry(bnum)->is_curved()) {            
			ebdry(bnum)->comm_exchange(boundary::all,0,boundary::master_slave);
		}
	}
	
	
	for(bnum=0;bnum<nebd;++bnum) {
		if (ebdry(bnum)->is_comm() && hp_ebdry(bnum)->is_curved()) {            
			ebdry(bnum)->comm_wait(boundary::all,0,boundary::master_slave);
			
			if (!ebdry(bnum)->is_frst()) {
			count = 0;
			for(i=ebdry(bnum)->nseg-1;i>=0;--i) {
				msgn = 1;
				for(m=0;m<basis::tet(log2p).em;++m) {
					for(n=0;n<ND;++n)
						hp_ebdry(bnum)->crde(i,m,n) = msgn*ebdry(bnum)->frcvbuf(0,count++);
					msgn *= -1;
				}
			}
			}
		}
	}
	
	if (!fm0) return;
	
	tc0load(ug.f.data(),0,fm0-1,ug.f.extent(secondDim));
	tmsgpass(boundary::all,0,boundary::symmetric);
	tc0wait_rcv(ug.f.data(),0,fm0-1,ug.f.extent(secondDim));      

	
	// need to match curved edges on face boundaries too fix me temp
	
	
	
	
//   /* Match curved sides */  // MAJOR FIXME HERE
//   for(bnum=0;bnum<nfbd;++bnum) {
//      if (fbdry(bnum)->is_comm() && hp_fbdry(bnum)->is_curved()) {            
//         count = 0;
//         for(i=0;i<fbdry(bnum)->nseg;++i) {
//            for(m=0;m<basis::tet(log2p).fm;++m) {
//               for(n=0;n<ND;++n)
//                  fbdry(bnum)->fsndbuf(count++) = hp_fbdry(bnum)->crds(i,m,n);
//            }
//         }
//         fbdry(bnum)->sndsize() = count;
//         fbdry(bnum)->sndtype() = boundary::flt_msg;
//         fbdry(bnum)->comm_prepare(boundary::all,0,boundary::master_slave);
//      }
//   }
//   
//   for(bnum=0;bnum<nfbd;++bnum) {
//      if (fbdry(bnum)->is_comm() && hp_fbdry(bnum)->is_curved()) {            
//         fbdry(bnum)->comm_exchange(boundary::all,0,boundary::master_slave);
//      }
//   }
//   
//   
//   for(bnum=0;bnum<nfbd;++bnum) {
//      if (fbdry(bnum)->is_comm() && hp_fbdry(bnum)->is_curved()) {            
//         fbdry(bnum)->comm_wait(boundary::all,0,boundary::master_slave);
//         
//         if (!fbdry(bnum)->is_frst()) {
//            count = 0;
//            for(i=fbdry(bnum)->nseg-1;i>=0;--i) {
//               msgn = 1;
//               for(m=0;m<basis::tet(log2p).fm;++m) {
//                  for(n=0;n<ND;++n)
//                     hp_fbdry(bnum)->crds(i,m,n) = msgn*fbdry(bnum)->frcvbuf(0,count++);
//                  msgn *= -1;
//               }
//            }
//         }
//      }
//   }      
	
	return;
}

