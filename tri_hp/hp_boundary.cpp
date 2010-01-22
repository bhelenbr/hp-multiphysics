/*
 *  hp_boundary.cpp
 *  spectral_hp
 *
 *  Created by Brian Helenbrook on 11/2/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp.h"
#include "hp_boundary.h"
#include <blitz/tinyvec-et.h>
#include <libbinio/binwrap.h>
#include <myblas.h>

//#define MPDEBUG

hp_vrtx_bdry* tri_hp::getnewvrtxobject(int bnum, input_map &bdrydata) {
	hp_vrtx_bdry *temp = new hp_vrtx_bdry(*this,*vbdry(bnum));  
	gbl->vbdry_gbls(bnum) = temp->create_global_structure();
	return(temp);
}

hp_edge_bdry* tri_hp::getnewsideobject(int bnum, input_map &bdrydata) {
	hp_edge_bdry *temp = new hp_edge_bdry(*this,*ebdry(bnum));
	gbl->ebdry_gbls(bnum) = temp->create_global_structure();
	return(temp);
}

void hp_edge_bdry::smatchsolution_snd(FLT *sdata, int bgn, int end, int stride) {
	/* CAN'T USE sfinalrcv BECAUSE OF CHANGING SIGNS */
	int j,k,n,count,offset,sind,sign;

	if (!base.is_comm()) return;

	int ebp1 = end-bgn+1;
	if (base.is_frst()) {
		count = 0;
		for(j=0;j<base.nseg;++j) {
			sind = base.seg(j);
			offset = (sind*stride +bgn)*x.NV;
			for(k=0;k<ebp1;++k) {
				for(n=0;n<x.NV;++n) {
					base.fsndbuf(count++) = sdata[offset++];
				}
			}
		}
	}
	else {
		int bgnsign = (bgn % 2 ? -1 : 1);
		count = 0;

		for(j=base.nseg-1;j>=0;--j) {
			sind = base.seg(j);
			offset = (sind*stride +bgn)*x.NV;
			sign = bgnsign;
			for (k=bgn;k<=end;++k) {
				for(n=0;n<x.NV;++n) {
					base.fsndbuf(count++) = sign*sdata[offset++];
				}
				sign *= -1;
			}
		}    
	}
	base.sndsize() = count;
	base.sndtype() = boundary::flt_msg;
}

void hp_edge_bdry::smatchsolution_rcv(FLT *sdata, int bgn, int end, int stride) {
	/* CAN'T USE sfinalrcv BECAUSE OF CHANGING SIGNS */
	int j,k,n,count,offset,sind,sign;

	bool reload = base.comm_finish(boundary::all,0,boundary::symmetric,boundary::average);
	if (!reload) return;

	int ebp1 = end -bgn +1;
	if (base.is_frst()) {
		count = 0;
		for(j=0;j<base.nseg;++j) {
			sind = base.seg(j);
			offset = (sind*stride +bgn)*x.NV;
			for(k=0;k<ebp1;++k) {
				for(n=0;n<x.NV;++n) {
					sdata[offset++] = base.fsndbuf(count++);
				}
			}
		}
	}
	else {
		int bgnsign = (bgn % 2 ? -1 : 1);
		int count = 0;

		for(j=base.nseg-1;j>=0;--j) {
			sind = base.seg(j);
			offset = (sind*stride +bgn)*x.NV;
			sign = bgnsign;
			for (k=bgn;k<=end;++k) {
				for(n=0;n<x.NV;++n) {
					sdata[offset++] = sign*base.fsndbuf(count++);
				}
				sign *= -1;
			}
		}    
	}
}


void hp_edge_bdry::copy(const hp_edge_bdry &bin) {

	if (!curved || !x.sm0) return;

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
		crv.resize(base.maxseg,x.sm0);
		for(i=1;i<x.gbl->nhist+1;++i)
			crvbd(i).resize(base.maxseg,x.sm0);
		crvbd(0).reference(crv);
	}

	dxdt.resize(x.log2pmax+1,base.maxseg);

	base.resize_buffers(base.maxseg*(x.sm0+2)*x.NV);

	return;
}

void hp_edge_bdry::output(std::ostream& fout, tri_hp::filetype typ,int tlvl) {
	int j,m,n;

	switch(typ) {
		case(tri_hp::text):
			fout << base.idprefix << " " << mytype << std::endl;
			if (curved) {
				fout << "p0: " << x.p0 << std::endl;

				for(j=0;j<base.nseg;++j) {
					for(m=0;m<x.sm0;++m) {
						for(n=0;n<tri_mesh::ND;++n)
							fout << crvbd(tlvl)(j,m)(n) << ' ';
						fout << std::endl;
					}
				}
			}
			break;

		case(tri_hp::binary):
			if (curved) {
				binowstream bout(&fout);
				bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::BigEndian)),1);
				bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::FloatIEEE)),1);
				bout.writeInt(x.p0,sizeof(int));

				for(j=0;j<base.nseg;++j) {
					for(m=0;m<x.sm0;++m) {
						for(n=0;n<tri_mesh::ND;++n)
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

void hp_edge_bdry::input(ifstream& fin,tri_hp::filetype typ,int tlvl) {
	int j,m,n,pmin;
	std::string idin, mytypein;    
	switch(typ) {
		case(tri_hp::text):
			fin >> idin >> mytypein;
			if (curved) { 
				fin.ignore(80,':');
				fin >>  pmin;
				pmin = x.p0;
				for(j=0;j<base.nseg;++j) {
					for(m=0;m<pmin -1;++m) {
						for(n=0;n<tri_mesh::ND;++n)
							fin >> crvbd(tlvl)(j,m)(n);
					}
					for(m=pmin-1;m<x.sm0;++m) {
						for(n=0;n<tri_mesh::ND;++n)
							crvbd(tlvl)(j,m)(n) = 0.0;
					}
				}
			}
			break;
		case(tri_hp::binary):
			if (curved) {
				biniwstream bin(&fin);

				/* HEADER INFORMATION */
				bin.setFlag(binio::BigEndian,bin.readInt(1));
				bin.setFlag(binio::FloatIEEE,bin.readInt(1));

				pmin = bin.readInt(sizeof(int));
				for(j=0;j<base.nseg;++j) {
					for(m=0;m<pmin-1;++m) {
						for(n=0;n<tri_mesh::ND;++n)
							crvbd(tlvl)(j,m)(n) = bin.readFloat(binio::Double);
					}
	                for(m=pmin-1;m<x.sm0;++m) {
						for(n=0;n<tri_mesh::ND;++n)
							crvbd(tlvl)(j,m)(n) = 0.0;
					}				
				}
			}
			break;

		default:
			break;
	}
}

void hp_edge_bdry::setvalues(init_bdry_cndtn *ibc, Array<int,1>& dirichlets, int ndirichlets) {
	int j,k,m,n,v0,v1,sind,indx,info;
	TinyVector<FLT,tri_mesh::ND> pt;
	char uplo[] = "U";

	/* UPDATE BOUNDARY CONDITION VALUES */
	j = 0;
	do {
		sind = base.seg(j);
		v0 = x.seg(sind).pnt(0);
		for(n=0;n<ndirichlets;++n)
			x.ug.v(v0,dirichlets(n)) = ibc->f(dirichlets(n),x.pnts(v0),x.gbl->time);
	} while(++j < base.nseg);
	v0 = x.seg(sind).pnt(1);
	for(n=0;n<ndirichlets;++n)
		x.ug.v(v0,dirichlets(n)) = ibc->f(dirichlets(n),x.pnts(v0),x.gbl->time);

	/*******************/    
	/* SET SIDE VALUES */
	/*******************/
	for(j=0;j<base.nseg;++j) {
		sind = base.seg(j);
		v0 = x.seg(sind).pnt(0);
		v1 = x.seg(sind).pnt(1);

		if (is_curved()) {
			x.crdtocht1d(sind);
			for(n=0;n<tri_mesh::ND;++n)
				basis::tri(x.log2p)->proj1d(&x.cht(n,0),&x.crd(n)(0,0),&x.dcrd(n,0)(0,0));
		}
		else {
			for(n=0;n<tri_mesh::ND;++n) {
				basis::tri(x.log2p)->proj1d(x.pnts(v0)(n),x.pnts(v1)(n),&x.crd(n)(0,0));

				for(k=0;k<basis::tri(x.log2p)->gpx();++k)
					x.dcrd(n,0)(0,k) = 0.5*(x.pnts(v1)(n)-x.pnts(v0)(n));
			}
		}

		if (basis::tri(x.log2p)->sm()) {
			for(n=0;n<ndirichlets;++n)
				basis::tri(x.log2p)->proj1d(x.ug.v(v0,dirichlets(n)),x.ug.v(v1,dirichlets(n)),&x.res(dirichlets(n))(0,0));

			for(k=0;k<basis::tri(x.log2p)->gpx(); ++k) {
				pt(0) = x.crd(0)(0,k);
				pt(1) = x.crd(1)(0,k);
				for(n=0;n<ndirichlets;++n)
					x.res(dirichlets(n))(0,k) -= ibc->f(dirichlets(n),pt,x.gbl->time);
			}
			for(n=0;n<ndirichlets;++n)
				basis::tri(x.log2p)->intgrt1d(&x.lf(dirichlets(n))(0),&x.res(dirichlets(n))(0,0));

			indx = sind*x.sm0;
			for(n=0;n<ndirichlets;++n) {
				PBTRS(uplo,basis::tri(x.log2p)->sm(),basis::tri(x.log2p)->sbwth(),1,(double *) &basis::tri(x.log2p)->sdiag1d(0,0),basis::tri(x.log2p)->sbwth()+1,&x.lf(dirichlets(n))(2),basis::tri(x.log2p)->sm(),info);
				for(m=0;m<basis::tri(x.log2p)->sm();++m) 
					x.ug.s(sind,m,dirichlets(n)) = -x.lf(dirichlets(n))(2+m);
			}
		}
	}
	return;
}


void hp_edge_bdry::curv_init(int tlvl) {
	int i,j,m,n,v0,v1,sind,info;
	TinyVector<FLT,2> pt;
	char uplo[] = "U";


	/* SKIP END VERTICES */
	for(j=1;j<base.nseg;++j) {
		sind = base.seg(j);
		v0 = x.seg(sind).pnt(0);
		base.mvpttobdry(j,-1.0, x.pnts(v0));
	}
//    v0 = x.seg(sind).pnt(1);
//    base.mvpttobdry(base.nseg-1,1.0, x.pnts(v0));

	if (!curved || basis::tri(x.log2p)->p() == 1) return;

	/*****************************/
	/* SET UP HIGHER ORDER MODES */
	/*****************************/
	for(j=0;j<base.nseg;++j) {
		sind = base.seg(j);

		v0 = x.seg(sind).pnt(0);
		v1 = x.seg(sind).pnt(1);

		for(n=0;n<tri_mesh::ND;++n) 
			basis::tri(x.log2p)->proj1d(x.pnts(v0)(n),x.pnts(v1)(n),&x.crd(n)(0,0));

		for(i=0;i<basis::tri(x.log2p)->gpx();++i) {
			pt(0) = x.crd(0)(0,i);
			pt(1) = x.crd(1)(0,i);
			base.mvpttobdry(j,basis::tri(x.log2p)->xp(i),pt);
			x.crd(0)(0,i) -= pt(0);
			x.crd(1)(0,i) -= pt(1);
		}

		for(n=0;n<tri_mesh::ND;++n) {
			basis::tri(x.log2p)->intgrt1d(&x.cf(n,0),&x.crd(n)(0,0));
			DPBTRS(uplo,basis::tri(x.log2p)->sm(),basis::tri(x.log2p)->sbwth(),1,(double *) &basis::tri(x.log2p)->sdiag1d(0,0),basis::tri(x.log2p)->sbwth()+1,&x.cf(n,2),basis::tri(x.log2p)->sm(),info);

			for(m=0;m<basis::tri(x.log2p)->sm();++m)
				crvbd(tlvl)(j,m)(n) = -x.cf(n,m+2);
		}

		/* TEST FOR A CIRCLE 
		basis::tri(x.log2p)->ptvalues1d(0.0);
		basis::tri(x.log2p)->ptprobe1d(tri_mesh::ND,&pt(0),&x.cht(0,0),MXTM);
		*gbl->log << pt << ' ' << pt(0)*pt(0) +pt(1)*pt(1) << std::endl;
		*/

	}
	return;
}

void hp_edge_bdry::findandmovebdrypt(TinyVector<FLT,2>& xp,int &bel,FLT &psi) const {
	int sind,v0,v1,iter;
	FLT dx,dy,ol,roundoff,dpsi;
	TinyVector<FLT,2> pt;

	base.findbdrypt(xp,bel,psi);
	if (!curved) {
		base.edge_bdry::mvpttobdry(bel,psi,xp);
		basis::tri(x.log2p)->ptvalues1d(psi);
		return;
	}

	sind = base.seg(bel);
	v0 = x.seg(sind).pnt(0);
	v1 = x.seg(sind).pnt(1);
	dx = x.pnts(v1)(0) - x.pnts(v0)(0);
	dy = x.pnts(v1)(1) - x.pnts(v0)(1);
	ol = 2./(dx*dx +dy*dy);
	dx *= ol;
	dy *= ol;

	/* FIND PSI SUCH THAT TANGENTIAL POSITION ALONG LINEAR SIDE STAYS THE SAME */
	/* THIS WAY, MULTIPLE CALLS WILL NOT GIVE DIFFERENT RESULTS */ 
	x.crdtocht1d(sind);

	iter = 0;
	roundoff = 10.0*EPSILON*(1.0 +(fabs(xp(0)*dx) +fabs(xp(1)*dy)));
	do {
		basis::tri(x.log2p)->ptprobe1d(x.ND,pt.data(),psi,&x.cht(0,0),MXTM);
		dpsi = (pt(0) -xp(0))*dx +(pt(1) -xp(1))*dy;
		psi -= dpsi;
		if (iter++ > 100) {
			*x.gbl->log << "#Warning: max iterations for curved side in bdry_locate type: " << base.idnum << " seg: " << bel << " sind: " << sind << " loc: " << xp << " dpsi: " << dpsi << std::endl;
			break;
		}  
	} while (fabs(dpsi) > roundoff);
	xp = pt;
}

void hp_edge_bdry::mvpttobdry(int bel,FLT psi,TinyVector<FLT,2> &xp) {

	/* SOLUTION IS BEING ADAPTED MUST GET INFO FROM ADAPT STORAGE */
	/* FIRST GET LINEAR APPROXIMATION TO LOCATION */
	base.edge_bdry::mvpttobdry(bel, psi, xp);
	adapt_storage->findandmovebdrypt(xp,bel,psi);

	return;
}

#ifdef DIRK
void hp_edge_bdry::tadvance() {
	int stage = x.gbl->substep +x.gbl->esdirk;  

	if (x.p0 > 1 && curved) {
		if (stage) {
			/* BACK CALCULATE K TERM */
			for(int j=0;j<base.nseg;++j) {
				for(int m=0;m<basis::tri(x.log2p)->sm();++m) {
					for(int n=0;n<tri_mesh::ND;++n)
						crvbd(stage+1)(j,m)(n) = (crvbd(0)(j,m)(n)-crvbd(1)(j,m)(n))*x.gbl->adirk(stage-1,stage-1);
				}
			}
		}

		if (x.gbl->substep == 0) {
			/* STORE TILDE W */
			for(int j=0;j<base.nseg;++j) {
				for(int m=0;m<basis::tri(x.log2p)->sm();++m) {
					for(int n=0;n<tri_mesh::ND;++n)
						crvbd(1)(j,m)(n) = crv(j,m)(n);
				}
			}
		}

		/* UPDATE TILDE W */
		for (int s=0;s<stage;++s) {            
			for(int j=0;j<base.nseg;++j) {
				for(int m=0;m<basis::tri(x.log2p)->sm();++m) {
					for(int n=0;n<tri_mesh::ND;++n) {
						crvbd(1)(j,m)(n) += x.gbl->adirk(stage,s)*crvbd(s+2)(j,m)(n);
					}
				}
			}
		}

		/* EXTRAPOLATE GUESS? */
		if (stage && x.gbl->dti > 0.0) {
			FLT constant =  x.gbl->cdirk(x.gbl->substep);
			crvbd(0)(Range(0,base.nseg-1),Range(0,basis::tri(x.log2p)->sm()-1)) += constant*crvbd(stage+1)(Range(0,base.nseg-1),Range(0,basis::tri(x.log2p)->sm()-1));
		}
	}

	calculate_unsteady_sources();

	/* EXTRAPOLATE SOLUTION OR MOVE TO NEXT TIME POSITION */
	if (!coupled) curv_init();

	return;
}

void hp_edge_bdry::calculate_unsteady_sources() {
	int i,j,n,sind;

	for(i=0;i<=x.log2pmax;++i) {
		for(j=0;j<base.nseg;++j) {
			sind = base.seg(j);
			x.crdtocht1d(sind,1);
			for(n=0;n<tri_mesh::ND;++n)
				basis::tri(i)->proj1d(&x.cht(n,0),&dxdt(i,j)(n,0));
		}
	}

	return;
}
#else
/* BACKWARDS DIFFERENCE STUFF */
void hp_edge_bdry::tadvance() {    
	if (x.p0 > 1 && curved) {
		for (int i=x.gbl->nhist-2;i>=0; --i) {            

			/* BACK CALCULATE K TERM */
			for(int j=0;j<base.nseg;++j) {
				for(int m=0;m<basis::tri(x.log2p)->sm();++m) {
					for(int n=0;n<tri_mesh::ND;++n)
						crvbd(i+1)(j,m)(n) = crvbd(i)(j,m)(n);
				}
			}
		}

		/* EXTRAPOLATE GUESS? */
		if (x.gbl->dti > 0.0 && x.gbl->tstep > 1) {
			crvbd(0)(Range(0,base.nseg-1),Range(0,basis::tri(x.log2p)->sm()-1)) += crvbd(1)(Range(0,base.nseg-1),Range(0,basis::tri(x.log2p)->sm()-1)) -crvbd(2)(Range(0,base.nseg-1),Range(0,basis::tri(x.log2p)->sm()-1));
		}
	}

	calculate_unsteady_sources();

	/* EXTRAPOLATE SOLUTION OR MOVE TO NEXT TIME POSITION */
	if (!coupled && curved) curv_init();

	return;
}

void hp_edge_bdry::calculate_unsteady_sources() {
	int j,n,sind;
	TinyMatrix<FLT,tri_mesh::ND,MXGP> crd;

	for (int level=1;level<min(x.gbl->nhist,x.gbl->tstep+1);++level) {
		for(int p=0;p<=x.log2pmax;++p) {
			for(j=0;j<base.nseg;++j) {
				dxdt(p,j) = 0.0;
				sind = base.seg(j);
				x.crdtocht1d(sind,level);
				for(n=0;n<tri_mesh::ND;++n) {
					basis::tri(p).proj1d(&x.cht(n,0),&crd(n,0));
					for (int i=0;i<basis::tri(p).gpx;++i)
						dxdt(p,j)(n,i) += x.gbl->bd(level)/x.gbl->bd(0)*crd(n,i);
				}
			}
		}
	}

	return;
}
#endif


void tri_hp::vc0load(int phase, FLT *pdata, int vrtstride) {
	int i;

	/* SEND COMMUNICATIONS TO ADJACENT MESHES */\
	for(i=0;i<nebd;++i) 
		hp_ebdry(i)->pmatchsolution_snd(phase,pdata,vrtstride);
	for(i=0;i<nvbd;++i)
		hp_vbdry(i)->pmatchsolution_snd(phase,pdata,vrtstride);

	for(i=0;i<nebd;++i)
		ebdry(i)->comm_prepare(boundary::all_phased,phase,boundary::symmetric);
	for(i=0;i<nvbd;++i)
		vbdry(i)->comm_prepare(boundary::all_phased,phase,boundary::symmetric);

	return;
}
int tri_hp::vc0wait_rcv(int phase, FLT *pdata, int vrtstride) {
	int stop = 1;
	int i;

	for(i=0;i<nebd;++i) {
		stop &= ebdry(i)->comm_wait(boundary::all_phased,phase,boundary::symmetric);
	}

	for(i=0;i<nvbd;++i) {
		stop &= vbdry(i)->comm_wait(boundary::all_phased,phase,boundary::symmetric);
	}


	for(i=0;i<nebd;++i) 
		hp_ebdry(i)->pmatchsolution_rcv(phase,pdata,vrtstride);
	for(i=0;i<nvbd;++i)
		hp_vbdry(i)->pmatchsolution_rcv(phase,pdata,vrtstride);

	return(stop);
}

int tri_hp::vc0rcv(int phase, FLT *pdata, int vrtstride) {
	int stop = 1,i;

	for(i=0;i<nebd;++i)
		stop &= ebdry(i)->comm_nowait(boundary::all_phased,phase,boundary::symmetric);
	for(i=0;i<nvbd;++i)
		stop &= vbdry(i)->comm_nowait(boundary::all_phased,phase,boundary::symmetric);

	for(i=0;i<nebd;++i) 
		hp_ebdry(i)->pmatchsolution_rcv(phase,pdata,vrtstride);
	for(i=0;i<nvbd;++i)
		hp_vbdry(i)->pmatchsolution_rcv(phase,pdata,vrtstride);

	return(stop);
}

void tri_hp::sc0load(FLT *sdata, int bgnmode, int endmode, int modestride) {
	int i;

	/* SEND COMMUNICATIONS TO ADJACENT MESHES */\
	for(i=0;i<nebd;++i) 
		hp_ebdry(i)->smatchsolution_snd(sdata,bgnmode,endmode,modestride);

	for(i=0;i<nebd;++i)
		ebdry(i)->comm_prepare(boundary::all,0,boundary::symmetric);

	return;
}

int tri_hp::sc0wait_rcv(FLT *sdata, int bgnmode, int endmode, int modestride) {
	int stop = 1;
	int i;

	for(i=0;i<nebd;++i) {
		stop &= ebdry(i)->comm_wait(boundary::all,0,boundary::symmetric);
	}

	for(i=0;i<nebd;++i) 
		hp_ebdry(i)->smatchsolution_rcv(sdata,bgnmode,endmode,modestride);

	return(stop);
}

int tri_hp::sc0rcv(FLT *sdata, int bgnmode, int endmode, int modestride) {
	int stop = 1,i;

	for(i=0;i<nebd;++i)
		stop &= ebdry(i)->comm_nowait(boundary::all,0,boundary::symmetric);

	for(i=0;i<nebd;++i) 
		hp_ebdry(i)->smatchsolution_rcv(sdata,bgnmode,endmode,modestride);

	return(stop);
}


void tri_hp::matchboundaries() {
	int i, m, n, msgn, bnum, count;
	int last_phase, mp_phase;

	/* Match boundary vertices */
	tri_mesh::matchboundaries();

	for(last_phase = false, mp_phase = 0; !last_phase; ++mp_phase) {
		vc0load(mp_phase,ug.v.data());
		pmsgpass(boundary::all_phased,mp_phase,boundary::symmetric);
		last_phase = true;
		last_phase &= vc0wait_rcv(mp_phase,ug.v.data());
	}

	if (!sm0) return;

	sc0load(ug.s.data(),0,sm0-1,ug.s.extent(secondDim));
	smsgpass(boundary::all,0,boundary::symmetric);
	sc0wait_rcv(ug.s.data(),0,sm0-1,ug.s.extent(secondDim));    

	/* Match curved sides */
	for(bnum=0;bnum<nebd;++bnum) {
		if (ebdry(bnum)->is_comm() && hp_ebdry(bnum)->is_curved()) {                
			count = 0;
			for(i=0;i<ebdry(bnum)->nseg;++i) {
				for(m=0;m<basis::tri(log2p)->sm();++m) {
					for(n=0;n<ND;++n)
						ebdry(bnum)->fsndbuf(count++) = hp_ebdry(bnum)->crds(i,m,n);
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
					for(m=0;m<basis::tri(log2p)->sm();++m) {
						for(n=0;n<ND;++n)
							hp_ebdry(bnum)->crds(i,m,n) = msgn*ebdry(bnum)->frcvbuf(0,count++);
						msgn *= -1;
					}
				}
			}
		}
	}

	return;
}

void hp_edge_bdry::findmax(FLT (*fxy)(TinyVector<FLT,2> &x)) {
	FLT ddpsi1, ddpsi2, psil, psir;
	TinyVector<FLT,2> xp, dx, maxloc, minloc;
	FLT max,min;
	int v0, sind;


	/* CALCULATE SLOPE AT ENDPOINT & TRANSMIT TO NEXT SURFACE */
	sind = base.seg(base.nseg-1);
	x.crdtocht1d(sind);
	basis::tri(x.log2p)->ptprobe1d(tri_mesh::ND,&xp(0),&dx(0),1.0,&x.cht(0,0),MXTM);
	ddpsi2 = (*fxy)(dx);
	if (base.vbdry(1) >= 0) {
		x.vbdry(base.vbdry(1))->vloadbuff(boundary::manifolds,&ddpsi2,0,1,1);
		x.vbdry(base.vbdry(1))->comm_prepare(boundary::manifolds,0,boundary::master_slave);
	}
	if (base.vbdry(0) >= 0) {
		x.vbdry(base.vbdry(0))->comm_prepare(boundary::manifolds,0,boundary::master_slave);
	}


	if (base.vbdry(1) >= 0) 
		x.vbdry(base.vbdry(1))->comm_exchange(boundary::manifolds,0,boundary::master_slave);
	if (base.vbdry(0) >= 0)
		x.vbdry(base.vbdry(0))->comm_exchange(boundary::manifolds,0,boundary::master_slave);

	if (base.vbdry(1) >= 0) 
		x.vbdry(base.vbdry(1))->comm_wait(boundary::manifolds,0,boundary::master_slave);
	if (base.vbdry(0) >= 0) {
		x.vbdry(base.vbdry(0))->comm_wait(boundary::manifolds,0,boundary::master_slave);
		if (x.vbdry(base.vbdry(0))->is_comm()) 
			ddpsi2 = x.vbdry(base.vbdry(0))->frcvbuf(0,0);
		else
			ddpsi2 = 0.0;
	}

	max = -1.0e99;
	min = 1.0e99;
	for(int indx=0;indx<base.nseg;++indx) {
		sind = base.seg(indx);
		x.crdtocht1d(sind);
		basis::tri(x.log2p)->ptprobe1d(tri_mesh::ND, &xp(0), &dx(0), -1.0, &x.cht(0,0), MXTM);
		ddpsi1 = (*fxy)(dx);
		if (ddpsi1 * ddpsi2 <= 0.0) {
			v0 = x.seg(base.seg(indx)).pnt(0);
			if ((*fxy)(x.pnts(v0)) > max) {
				maxloc[0] = x.pnts(v0)(0);
				maxloc[1] = x.pnts(v0)(1);
				max = (*fxy)(x.pnts(v0));
			}
			if ((*fxy)(x.pnts(v0)) < min) {
				minloc[0] = x.pnts(v0)(0);
				minloc[1] = x.pnts(v0)(1);
				min = (*fxy)(x.pnts(v0));
			}
			*x.gbl->log << "#LOCAL EXTREMA: " << x.pnts(v0)(0) << ' ' << x.pnts(v0)(1) << ' ' <<(*fxy)(x.pnts(v0)) << std::endl;
		}
		basis::tri(x.log2p)->ptprobe1d(tri_mesh::ND, &xp(0), &dx(0), 1.0, &x.cht(0,0), MXTM);
		ddpsi2 = (*fxy)(dx);
		if (ddpsi1 *ddpsi2 <= 0.0) {
			/* INTERIOR MAXIMUM */
			psil = -1.0;
			psir = 1.0;
			while (psir-psil > 1.0e-10) {
				basis::tri(x.log2p)->ptprobe1d(tri_mesh::ND, &xp(0), &dx(0), 0.5*(psil +psir), &x.cht(0,0), MXTM);
				if ((*fxy)(dx)*ddpsi1 < 0.0) 
					psir = 0.5*(psil+psir);
				else
					psil = 0.5*(psil+psir);
			}
			if ((*fxy)(xp) > max) {
				maxloc[0] = xp[0];
				maxloc[1] = xp[1];
				max = (*fxy)(xp);
			}
			if ((*fxy)(xp) < min) {
				minloc[0] = xp[0];
				minloc[1] = xp[1];
				min = (*fxy)(xp);
			}
			*x.gbl->log << "#LOCAL EXTREMA: " << xp[0] << ' ' << xp[1] << ' ' << (*fxy)(xp) << std::endl;
		}  
	}
	*x.gbl->log << "#MAX EXTREMA: " << maxloc[0] << ' ' << maxloc[1] << ' ' << max << std::endl;
	*x.gbl->log << "#MIN EXTREMA: " << minloc[0] << ' ' << minloc[1] << ' ' << min << std::endl;

	return;
}

 void hp_edge_bdry::findintercept(FLT (*fxy)(TinyVector<FLT,2> &x)) {
	FLT psil, psir;
	TinyVector<FLT,2> xp, dx;
	int v0, sind;
	FLT vl, vr;

	sind = base.seg(0);
	x.crdtocht1d(sind);
	v0 = x.seg(sind).pnt(0);
	vl = (*fxy)(x.pnts(v0));

	for(int indx=0;indx<base.nseg;++indx) {
		sind = base.seg(indx);
		x.crdtocht1d(sind);
		v0 = x.seg(sind).pnt(1);
		vr = (*fxy)(x.pnts(v0));

		if (vl*vr <= 0.0) {
			/* INTERIOR INTERCEPT */
			psil = -1.0;
			psir = 1.0;
			while (psir-psil > 1.0e-10) {
				basis::tri(x.log2p)->ptprobe1d(tri_mesh::ND,&xp(0),&dx(0),0.5*(psil+psir),&x.cht(0,0),MXTM);
				if ((*fxy)(xp)*vl < 0.0) 
					psir = 0.5*(psil+psir);
				else
					psil = 0.5*(psil+psir);
			}
			*x.gbl->log << "#INTERSECTION: " << xp[0] << ' ' << xp[1] << ' ' << (*fxy)(xp) << std::endl;
		}
		vl = vr; 
	}

	return;
}

void hp_edge_bdry::rsdl(int stage) {
	x.lf = 0.0;
	for(int j=0;j<base.nseg;++j) {
		int sind = base.seg(j);
		int v0 = x.seg(sind).pnt(0);
		int v1 = x.seg(sind).pnt(1);
		
		x.ugtouht1d(sind);
		element_rsdl(j,stage);
		
		for(int n=0;n<x.NV;++n)
			x.gbl->res.v(v0,n) += x.lf(n)(0);
		
		for(int n=0;n<x.NV;++n)
			x.gbl->res.v(v1,n) += x.lf(n)(1);
		
		for(int k=0;k<basis::tri(x.log2p)->sm();++k) {
			for(int n=0;n<x.NV;++n)
				x.gbl->res.s(sind,k,n) += x.lf(n)(k+2);
		}
	}
}


void hp_edge_bdry::element_jacobian(int indx, Array<FLT,2>& K) {
	int sm = basis::tri(x.log2p)->sm();	
	FLT dw = 1.0e-4;// fix me temp make a global value?
	Array<FLT,2> Rbar(x.NV,sm+2);
	Array<int,1> loc_to_glo(x.NV*(sm+2));

	/* Calculate and store initial residual */
	int sind = base.seg(indx);
	x.ugtouht1d(sind);
	
	x.lf = 0.0;
	element_rsdl(indx,0);

	for(int k=0;k<sm+2;++k) {
		for(int n=0;n<x.NV;++n) {
			Rbar(n,k) = x.lf(n)(k);
		}
	}
	
	if (x.mmovement != x.coupled_deformable) {
		/* Numerically create Jacobian */
		int kcol = 0;
		for(int mode = 0; mode < sm+2; ++mode){
			for(int var = 0; var < x.NV; ++var){
				x.uht(var)(mode) += dw;
				
				x.lf = 0.0;
				element_rsdl(indx,0);
				
				int krow = 0;
				for(int k=0;k<sm+2;++k)
					for(int n=0;n<x.NV;++n)
						K(krow++,kcol) = (x.lf(n)(k)-Rbar(n,k))/dw;
						
						++kcol;
				
				x.uht(var)(mode) -= dw;
			}
		}
	}
	else {
		/* Numerically create Jacobian */
		int kcol = 0;
		for(int mode = 0; mode < 2; ++mode) {
			for(int var = 0; var < x.NV; ++var){
				x.uht(var)(mode) += dw;
				
				x.lf = 0.0;
				element_rsdl(indx,0);
				
				int krow = 0;
				for(int k=0;k<sm+2;++k)
					for(int n=0;n<x.NV;++n)
						K(krow++,kcol) = (x.lf(n)(k)-Rbar(n,k))/dw;
						
						++kcol;
				
				x.uht(var)(mode) -= dw;
			}
			
			for(int var = 0; var < tri_mesh::ND; ++var){
				x.pnts(x.seg(sind).pnt(mode))(var) += dw;
				
				x.lf = 0.0;
				element_rsdl(indx,0);
				
				int krow = 0;
				for(int k=0;k<sm+2;++k)
					for(int n=0;n<x.NV;++n)
						K(krow++,kcol) = (x.lf(n)(k)-Rbar(n,k))/dw;
						
						++kcol;
				
				x.pnts(x.seg(sind).pnt(mode))(var) -= dw;
			}
		}
			
		for(int mode = 2; mode < sm+2; ++mode){
			for(int var = 0; var < x.NV; ++var){
				x.uht(var)(mode) += dw;
				
				x.lf = 0.0;
				element_rsdl(indx,0);
				
				int krow = 0;
				for(int k=0;k<sm+2;++k)
					for(int n=0;n<x.NV;++n)
						K(krow++,kcol) = (x.lf(n)(k)-Rbar(n,k))/dw;
						
						++kcol;
				
				x.uht(var)(mode) -= dw;
			}
		}
	}
}

#ifdef petsc
void hp_edge_bdry::petsc_jacobian() {
	int sm = basis::tri(x.log2p)->sm();	

	
	int vdofs;
	if (x.mmovement != x.coupled_deformable)
		vdofs = x.NV;
	else
		vdofs = x.NV+x.ND;
		
	/* If moving mesh, but not coupled, then vertices could slide along side */
	Array<FLT,2> K(x.NV*(sm+2),vdofs*2 +x.NV*sm);
	Array<int,1> rows(x.NV*(sm+2));
	Array<int,1> cols(vdofs*2 +x.NV*sm);
		
	for(int j=0;j<base.nseg;++j) {
		int sind = base.seg(j);
			
		int rind = 0;
		int cind = 0;
		for(int k=0;k<2;++k) {
			int gindx = vdofs*x.seg(sind).pnt(k);
			for(int n=0;n<x.NV;++n) {
				rows(rind++) = gindx;
				cols(cind++) = gindx++;
			}
			for(int n=x.NV;n<vdofs;++n) {
				cols(cind++) = gindx++;
			}
		}
		
		/* EDGE MODES */
		if (sm > 0) {
			int gbl_eind = x.npnt*vdofs + sind*sm*x.NV;
			for (int m = 0; m < sm; ++m) {
				for(int n = 0; n < x.NV; ++n) {
					rows(rind++) = gbl_eind;
					cols(cind++) = gbl_eind++;
				}
			}
		}
		
		element_jacobian(j,K);
		
		MatSetValues(x.petsc_J,x.NV*(sm+2),rows.data(),vdofs*2 +x.NV*sm,cols.data(),K.data(),ADD_VALUES);
	}
}

#include <r_tri_boundary.h>
void hp_edge_bdry::petsc_jacobian_dirichlet() {
	/* DO R_MESH B.C.'s */
	if (x.mmovement == x.coupled_deformable) {
		/* This is ugly, but I need to be able to control this for different boundary types */
		for(int i=0;i<x.nebd;++i) {
			if (x.hp_ebdry(i) == this) {
				x.r_sbdry(i)->jacobian_dirichlet();
				break;
			}
		}
	}
}
		
#endif
