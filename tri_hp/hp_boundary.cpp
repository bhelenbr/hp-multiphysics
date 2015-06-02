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
#include <libbinio/binwrap.h>
#include <myblas.h>

//#define MPDEBUG
//#define OLDWAY1
//#define OLDWAY2
//#define OLDWAY3

#ifdef DEBUG_JAC
#define SKIP_EXCHANGE
#endif

#ifdef OLDWAY1
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
#else
void hp_edge_bdry::smatchsolution_snd(FLT *sdata, int bgn, int end, int stride) {
    /* CAN'T USE sfinalrcv BECAUSE OF CHANGING SIGNS */
    int j,k,count,offset,sind,sign;
    
    if (!base.is_comm()) return;
    
    int ebp1 = end-bgn+1;
    if (base.is_frst()) {
        count = 0;
        for(j=0;j<base.nseg;++j) {
            sind = base.seg(j);
            offset = (sind*stride +bgn)*x.NV;
            for(k=0;k<ebp1;++k) {
                for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
                    base.fsndbuf(count++) = sdata[offset +*n];
                }
                offset+=x.NV;
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
                for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
                    base.fsndbuf(count++) = sign*sdata[offset +*n];
                }
                sign *= -1;
                offset+=x.NV;
            }
        }
    }
    base.sndsize() = count;
    base.sndtype() = boundary::flt_msg;
}

void hp_edge_bdry::smatchsolution_rcv(FLT *sdata, int bgn, int end, int stride) {
    /* CAN'T USE sfinalrcv BECAUSE OF CHANGING SIGNS */
    int j,k,count,offset,sind,sign;
    
    bool reload = base.comm_finish(boundary::all,0,boundary::symmetric,boundary::average);
    if (!reload) return;
    
    int ebp1 = end -bgn +1;
    if (base.is_frst()) {
        count = 0;
        for(j=0;j<base.nseg;++j) {
            sind = base.seg(j);
            offset = (sind*stride +bgn)*x.NV;
            for(k=0;k<ebp1;++k) {
                for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
                    sdata[offset+*n] = base.fsndbuf(count++);
                }
                offset+=x.NV;
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
                for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
                    sdata[offset+*n] = sign*base.fsndbuf(count++);
                }
                sign *= -1;
                offset+=x.NV;
            }
        }    
    }
}
#endif

void hp_edge_bdry::copy(const hp_edge_bdry &bin) {

	if (!curved || !x.sm0) return;

	for(int i=0;i<x.gbl->nadapt; ++i)
		crvbd(i)(Range(0,base.nseg-1),Range::all()) = bin.crvbd(i)(Range(0,base.nseg-1),Range::all());
}

void hp_vrtx_bdry::init(input_map& inmap,void* gbl_in) {
	std::string keyword,ibcname;
	
	/* FIND INITIAL CONDITION TYPE */
	keyword = base.idprefix + "_ibc";
	if (inmap.get(keyword,ibcname)) {
		ibc = x.getnewibc(ibcname);
		ibc->init(inmap,keyword);
	}
	
	keyword = base.idprefix + "_coupled";
	coupled = false;
	inmap.get(keyword,coupled);
	
	keyword = base.idprefix + "_frozen";
	frozen = false;
	inmap.get(keyword,frozen);
	
	keyword = base.idprefix +"_report";
	report_flag = false;
	inmap.get(keyword,report_flag);
	
	Array<int,1> atemp(x.NV);
	if (inmap.get(base.idprefix+"_hp_typelist", atemp.data(), x.NV)) {
		for (int n=0;n<x.NV;++n) {
			type[n] = static_cast<bctypes>(atemp(n));
			if (type[n] == essential) {
				essential_indices.push_back(n);
			}
		}
	}
	
	std::string val;
	if (!inmap.getline(base.idprefix +"_c0_indices",val)) {
		/* Default is that all variables are continuous across boundary */
		for(int n=0;n<x.NV;++n) {
			c0_indices.push_back(n);
		}
	}
	else {
		istringstream data;
		data.str(val);
		int ind;
		while (data >> ind)
			c0_indices.push_back(ind);
		data.clear();
	}
	c0_indices_xy = c0_indices;
	if (x.mmovement == x.coupled_deformable) {
		c0_indices_xy.push_back(x.NV);
		c0_indices_xy.push_back(x.NV+1);
	}
	
#ifdef petsc
	base.resize_buffers((x.NV+x.ND)*60*(3 +3*x.sm0+x.im0));  // Allows for 10 elements of jacobian entries to be sent 
#endif	
}

void hp_edge_bdry::init(input_map& inmap,void* gbl_in) {
	int i;
	std::string keyword,ibcname;
	std::istringstream data;
	std::string filename;
	
	/* FIND INITIAL CONDITION TYPE */
	keyword = base.idprefix + "_ibc";
	if (inmap.get(keyword,ibcname)) {
		ibc = x.getnewibc(ibcname);
		ibc->init(inmap,keyword);
	}

	if (base.is_comm() || inmap.find(base.idprefix+"_type") == inmap.end() || inmap[base.idprefix+"_type"] == "plain") {
		curved = false;
	}
	else {
		curved = true;
	}
	keyword = base.idprefix + "_curved";
	inmap.get(keyword,curved);
	
	keyword = base.idprefix + "_coupled";
	coupled = false;
	inmap.get(keyword,coupled);
	
	keyword = base.idprefix + "_frozen";
	frozen = false;
	inmap.get(keyword,frozen);
	
	keyword = base.idprefix +"_report";
	report_flag = false;
	inmap.get(keyword,report_flag);
	
	Array<int,1> atemp(x.NV);
	if (inmap.get(base.idprefix+"_hp_typelist", atemp.data(), x.NV)) {
		for (int n=0;n<x.NV;++n) {
			type[n] = static_cast<bctypes>(atemp(n));
			if (type[n] == essential) {
				essential_indices.push_back(n);
			}
		}
	}
	
	std::string val;
	if (!inmap.getline(base.idprefix +"_c0_indices",val)) {
		/* Default is that all variables are continuous across boundary */
		for(int n=0;n<x.NV;++n) {
			c0_indices.push_back(n);
		}
	}
	else {
		data.str(val);
		int ind;
		while (data >> ind)
			c0_indices.push_back(ind);
		data.clear();
	}
	c0_indices_xy = c0_indices;
	if (x.mmovement == x.coupled_deformable) {
		c0_indices_xy.push_back(x.NV);
		c0_indices_xy.push_back(x.NV+1);
	}

	if (curved && !x.coarse_level) {
		crvbd.resize(x.gbl->nhist+1);
		crv.resize(base.maxseg,x.sm0);
		for(i=1;i<x.gbl->nhist+1;++i)
			crvbd(i).resize(base.maxseg,x.sm0);
		crvbd(0).reference(crv);
	}
	dxdt.resize(x.log2pmax+1,base.maxseg);
	dxdt = 0.0;  // Do this to eliminate warning in valgrind when outputting before tadvance is called.

	if (inmap.find(base.idprefix+"_norm") == inmap.end())	{
		inmap[base.idprefix+"_norm"] = "0.0";
	}
	l2norm.init(inmap,base.idprefix+"_norm");

#ifndef petsc
	base.resize_buffers(base.maxseg*(x.sm0+2)*x.NV);
#else
	base.resize_buffers(base.maxseg*(x.sm0+2)*(x.NV+x.ND)*16*(3 +3*x.sm0+x.im0));  // Allows for 4 elements of jacobian entries to be sent 
#endif
	
	/* Set-up boundary flux functions */
	fluxes.resize(x.NV);
	Array<string,1> names(4);
	Array<int,1> dims(4);
	dims = x.ND;
	names(0) = "u";
	dims(0) = x.NV;
	names(1) = "x";
	names(2) = "xt";
	names(3) = "n";
	for(int n=0;n<x.NV;++n) {
		fluxes[n].set_arguments(4,dims,names);
	}
	
	input_map zeromap;
	zeromap["zero"] = "0.0";
	std::ostringstream nstr;
	for (int n=0;n<x.NV;++n) {
		if (type[n] == natural) {
			nstr.str("");
			nstr << base.idprefix << "_flux" << n << std::flush;
			if (inmap.find(nstr.str()) != inmap.end()) {
				fluxes[n].init(inmap,nstr.str());
			}
			else {
				fluxes[n].init(zeromap,"zero");
			}
		}
		else {
			fluxes[n].init(zeromap,"zero");
		}
	}
	
	return;
}

void hp_edge_bdry::output(const std::string& filename, tri_hp::filetype typ,int tlvl) {
	int j,m,n;
	
	switch(typ) {
		case(tri_hp::text): {
			std::string fname;
			fname = filename +"_" +x.gbl->idprefix +".txt";
			ofstream fout;
			fout.open(fname,std::ofstream::out | std::ofstream::app);
			
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
			fout.close();
			
			break;
		}

		case(tri_hp::binary): {
			if (curved) {
				std::string fname;
				fname = filename +"_" +x.gbl->idprefix +".bin";
				ofstream fout;
				fout.open(fname,std::ofstream::out | std::ofstream::app);
				
				binowstream bout(&fout);
				bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::BigEndian)),1);
				bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::FloatIEEE)),1);
				bout.writeInt(x.p0,sizeof(int));

				for(j=0;j<base.nseg;++j) {
					for(m=0;m<x.sm0;++m) {
						for(n=0;n<tri_mesh::ND;++n) {
							bout.writeFloat(crvbd(tlvl)(j,m)(n),binio::Double);
						}
					}
				}
				fout.close();
			}
			break;
		}
		
		case(tri_hp::tecplot): {
			if (!report_flag) break;
			
			std::string fname;
			fname = filename +"_" +base.idprefix +".dat";
			std::ofstream fout;
			fout.open(fname);
			
			fout << "VARIABLES=\"S\",\"X\",\"Y\",\"DXDT\",\"DYDT\",";
			for(n=0;n<x.NV;++n)
				fout << "\"V" << n << "\",\"DVTANG" << n << "\",\"DVNORM" << n << "\",";
			fout << "\nTITLE = " << base.idprefix << '\n'<< "ZONE\n";
			
			FLT circumference = 0.0;
			FLT l2error = 0.0;
			TinyVector<FLT,tri_mesh::ND> xpt,norm,mvel,visc;

			int ind = 0;
			do { 
				int sind = base.seg(ind);
				int tind = x.seg(sind).tri(0);        
				
				int seg;
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
				
				for (int i=0;i<basis::tri(x.log2p)->gpx();++i) {
					FLT arclength = sqrt(x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i) +x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i));
					FLT jcb =  basis::tri(x.log2p)->wtx(i)*RAD(x.crd(0)(0,i))*arclength;
					circumference += jcb;
	
					norm(0) = x.dcrd(1,0)(0,i);
					norm(1) = -x.dcrd(0,0)(0,i);                
					for(n=0;n<tri_mesh::ND;++n) {
						mvel(n) = x.gbl->bd(0)*(x.crd(n)(0,i) -dxdt(x.log2p,ind)(n,i));
#ifdef MESH_REF_VEL
						mvel(n) += x.gbl->mesh_ref_vel(n);
#endif
					}
					
					xpt(0) = x.crd(0)(0,i);
					xpt(1) = x.crd(1)(0,i);
					l2error += jcb*l2norm.Eval(xpt,x.gbl->time);
					
					jcb = 1./(arclength*(x.dcrd(0,0)(0,i)*x.dcrd(1,1)(0,i) -x.dcrd(1,0)(0,i)*x.dcrd(0,1)(0,i)));
									
					/* For normal derivative */
					visc[0] =  jcb*(x.dcrd(1,1)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,1)(0,i)*x.dcrd(0,0)(0,i));
					visc[1] = -jcb*(x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i));
					
					fout << circumference << ' ' << xpt(0) << ' ' << xpt(1) << ' ' << mvel(0) << ' ' << mvel(1) << ' ';
					
					for (n=0;n<x.NV;++n) {
						/* Output value, tangent, and normal derivatives */
						fout << x.u(n)(0,i) << ' ' << x.du(n,0)(0,i)/arclength 
							<< ' ' << -visc[0]*x.du(n,0)(0,i) -visc[1]*x.du(n,1)(0,i) << ' ';

					}
					fout << std::endl;
				}	
			} while (++ind < base.nseg);
			fout.close();
			
			streamsize oldprecision = (*x.gbl->log).precision(10);
			*x.gbl->log << base.idprefix << " circumference: " << circumference << std::endl;
			/* OUTPUT AUXILIARY FLUXES */
			*x.gbl->log << base.idprefix << " l2error: " << sqrt(l2error) << std::endl;
			(*x.gbl->log).precision(oldprecision);
			
			break;
	}

		default:
			break;
	}
	return;
}

void hp_edge_bdry::find_matching_boundary_name(input_map& inmap, std::string& matching_block, std::string& side_id) {

	std::ostringstream nstr;
	nstr << base.idnum << std::flush;
	side_id = "s" +nstr.str();
	
	if (!base.is_comm()) return;

	std::vector<std::string> blocks;
	if (!inmap.keys_with_ending(side_id+"_hp_type",blocks)) {
		*x.gbl->log << "Something screwy with surface blocks" << std::endl;
		sim::abort(__LINE__, __FILE__, x.gbl->log);
	}

	if (blocks.size() != 2) {
		*x.gbl->log << "Something screwy with surface blocks" << std::endl;
		sim::abort(__LINE__, __FILE__, x.gbl->log);
	}
	matching_block = blocks[0].substr(0,blocks[0].length()-9-side_id.length());
	if (matching_block == x.gbl->idprefix)
		matching_block = blocks[1].substr(0,blocks[1].length()-9-side_id.length());
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
						for(n=0;n<tri_mesh::ND;++n) {
							crvbd(tlvl)(j,m)(n) = bin.readFloat(binio::Double);
						}
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

void hp_edge_bdry::setvalues(init_bdry_cndtn *ibc, const std::vector<int>& indices) {
	int j,k,m,n,v0,v1,sind,indx,info;
	TinyVector<FLT,tri_mesh::ND> pt;
	char uplo[] = "U";

	/* UPDATE BOUNDARY CONDITION VALUES */
	j = 0;
	do {
		sind = base.seg(j);
		v0 = x.seg(sind).pnt(0);
		for(std::vector<int>::const_iterator n=indices.begin();n != indices.end();++n) 
			x.ug.v(v0,*n) = ibc->f(*n,x.pnts(v0),x.gbl->time);
	} while(++j < base.nseg);
	v0 = x.seg(sind).pnt(1);
	for(std::vector<int>::const_iterator n=indices.begin();n != indices.end();++n) 
		x.ug.v(v0,*n) = ibc->f(*n,x.pnts(v0),x.gbl->time);

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
			for(std::vector<int>::const_iterator n=indices.begin();n != indices.end();++n) 
				basis::tri(x.log2p)->proj1d(x.ug.v(v0,*n),x.ug.v(v1,*n),&x.res(*n)(0,0));

			for(k=0;k<basis::tri(x.log2p)->gpx(); ++k) {
				pt(0) = x.crd(0)(0,k);
				pt(1) = x.crd(1)(0,k);
				for(std::vector<int>::const_iterator n=indices.begin();n != indices.end();++n) 
					x.res(*n)(0,k) -= ibc->f(*n,pt,x.gbl->time);
			}
			for(std::vector<int>::const_iterator n=indices.begin();n != indices.end();++n) 
				basis::tri(x.log2p)->intgrt1d(&x.lf(*n)(0),&x.res(*n)(0,0));

			indx = sind*x.sm0;
			for(std::vector<int>::const_iterator n=indices.begin();n != indices.end();++n) {
				PBTRS(uplo,basis::tri(x.log2p)->sm(),basis::tri(x.log2p)->sbwth(),1,(double *) &basis::tri(x.log2p)->sdiag1d(0,0),basis::tri(x.log2p)->sbwth()+1,&x.lf(*n)(2),basis::tri(x.log2p)->sm(),info);
				for(m=0;m<basis::tri(x.log2p)->sm();++m) 
					x.ug.s(sind,m,*n) = -x.lf(*n)(2+m);
			}
		}
	}
	return;
}


void hp_edge_bdry::curv_init(int tlvl) {
	int i,j,m,n,v0,v1,sind,info;
	TinyVector<FLT,2> pt;
	char uplo[] = "U";

	j = 0;
	do {
		sind = base.seg(j);
		v0 = x.seg(sind).pnt(0);
		base.mvpttobdry(j,-1.0, x.pnts(v0));
	} while (++j < base.nseg);
	v0 = x.seg(sind).pnt(1);
	base.mvpttobdry(base.nseg-1,1.0, x.pnts(v0));

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
			FLT constant =  x.gbl->cdirk(x.gbl->substep)*x.gbl->extrapolate;
			crvbd(0)(Range(0,base.nseg-1),Range(0,basis::tri(x.log2p)->sm()-1)) += constant*crvbd(stage+1)(Range(0,base.nseg-1),Range(0,basis::tri(x.log2p)->sm()-1));
		}
	}

	calculate_unsteady_sources();

	/* EXTRAPOLATE SOLUTION OR MOVE TO NEXT TIME POSITION */
	if (!coupled) curv_init();
	
	if (!frozen) setvalues(ibc,essential_indices);

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
	
	setvalues(ibc,essential_indices);

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
	
	minloc = 0.0;
	maxloc = 0.0;


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

void hp_vrtx_bdry::rsdl(int stage) {
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;

	Array<FLT,1> lf(vdofs);
	element_rsdl(lf);

	for(int n=0;n<x.NV;++n)
		x.gbl->res.v(base.pnt,n) += lf(n);
	
}

void hp_edge_bdry::rsdl(int stage) {

	for(int j=0;j<base.nseg;++j) {
		int sind = base.seg(j);
		int v0 = x.seg(sind).pnt(0);
		int v1 = x.seg(sind).pnt(1);
		
		x.ugtouht1d(sind);
		element_rsdl(j,x.lf);
		
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

#ifdef petsc
#ifdef OLDWAY2
void hp_vrtx_bdry::non_sparse_snd(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi) {
	
	if (!base.is_comm()) return;

	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;

	/* Going to send all jacobian entries,  Diagonal entries for matching DOF's will be merged together not individual */
	/* Send number of non-zeros to matches */
	base.sndsize() = 0;
	base.sndtype() = boundary::int_msg;
	
	int pind = base.pnt*vdofs;
	for(int n=0;n<vdofs;++n) {
		base.isndbuf(base.sndsize()++) = nnzero(pind++);
	}
}


void hp_vrtx_bdry::non_sparse_rcv(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi) {
	
	if (!base.is_comm()) return;
	
	const int NV = x.NV;
	const int ND = tri_mesh::ND;
	
	int vdofs;
	if (x.mmovement != tri_hp::coupled_deformable) 
		vdofs = NV;
	else
		vdofs = ND+NV;
	
	/* Reload to avoid overlap with sides */
	int pind = base.pnt*vdofs;
	for(int n=0;n<vdofs;++n) {
		nnzero(pind) = base.isndbuf(n);
		nnzero_mpi(pind++) = 0;
	}
		
	for(int m=0;m<base.nmatches();++m) { 
		if (base.is_local(m)) {
			int count = 0;
			int pind = base.pnt*vdofs;
			for(int n=0;n<vdofs;++n) {
				nnzero(pind++) += base.ircvbuf(m,count++);
			}
		}
		else {
			int count = 0;
			int pind = base.pnt*vdofs;
			for(int n=0;n<vdofs;++n) {
				nnzero_mpi(pind++) += base.ircvbuf(m,count++);
			}
		}
	}
}

void hp_edge_bdry::non_sparse_snd(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi) {
	
	if (!base.is_comm()) return;
	
	const int sm=basis::tri(x.log2p)->sm();
	const int NV = x.NV;
	const int ND = tri_mesh::ND;
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;

	int begin_seg = x.npnt*vdofs;
			
	/* Going to send all jacobian entries,  Diagonal entries for matching DOF's will be merged together not individual */
	/* Send number of non-zeros to matches */
	base.sndsize() = 0;
	base.sndtype() = boundary::int_msg;
	for (int i=0;i<base.nseg;++i) {
		int sind = base.seg(i);
		int pind = x.seg(sind).pnt(0)*vdofs;
		for(int n=0;n<vdofs;++n)
			base.isndbuf(base.sndsize()++) = nnzero(pind++);
	}
	int sind = base.seg(base.nseg-1);
	int pind = x.seg(sind).pnt(1)*vdofs;
	for(int n=0;n<vdofs;++n)
		base.isndbuf(base.sndsize()++) = nnzero(pind++);
	
	/* Last thing to send is nnzero for edges */
	if (sm) {
		for (int i=0;i<base.nseg;++i) {
			int sind = base.seg(i);
			for (int m=0;m<sm;++m) {
				for(int n=0;n<NV;++n) {
					base.isndbuf(base.sndsize()++) = nnzero(begin_seg +sind*sm*NV +m*NV +n);
				}
			}
		}
	}
}

void hp_edge_bdry::non_sparse_rcv(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi) {
	
	if (!base.is_comm()) return;
	
	const int sm=basis::tri(x.log2p)->sm();
	const int NV = x.NV;
	const int ND = tri_mesh::ND;
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;

	int begin_seg = x.npnt*vdofs;
	
	for(int m=0;m<base.nmatches();++m) {
		if (base.is_local(m)) {
			int count = 0;
			for (int i=base.nseg-1;i>=0;--i) {
				int sind = base.seg(i);
				int pind = x.seg(sind).pnt(1)*vdofs;
				for(int n=0;n<vdofs;++n) {
					nnzero(pind++) += base.ircvbuf(m,count++);
				}
			}
			int sind = base.seg(0);
			int pind = x.seg(sind).pnt(0)*vdofs;
			for(int n=0;n<vdofs;++n) {
				nnzero(pind++) += base.ircvbuf(m,count++);
			}
			
			/* Now add to side degrees of freedom */
			if (sm) {
				for (int i=base.nseg-1;i>=0;--i) {
					int sind = base.seg(i);
					for (int mode=0;mode<sm;++mode) {
						for(int n=0;n<NV;++n) {
							nnzero(begin_seg +sind*sm*NV +mode*NV +n) += base.ircvbuf(m,count++);
						}
					}
				}
			}
		}
		else {
			int count = 0;
			for (int i=base.nseg-1;i>=0;--i) {
				int sind = base.seg(i);
				int pind = x.seg(sind).pnt(1)*vdofs;
				for(int n=0;n<vdofs;++n) {
					nnzero_mpi(pind++) = base.ircvbuf(m,count++);
				}
			}
			int sind = base.seg(0);
			int pind = x.seg(sind).pnt(0)*vdofs;
			for(int n=0;n<vdofs;++n) {
				nnzero_mpi(pind++) = base.ircvbuf(m,count++);
			}
			
			
			/* Now add to side degrees of freedom */
			if (sm) {
				for (int i=base.nseg-1;i>=0;--i) {
					int sind = base.seg(i);
					for (int mode=0;mode<sm;++mode) {
						for(int n=0;n<NV;++n) {
							nnzero_mpi(begin_seg +sind*sm*NV +mode*NV +n) += base.ircvbuf(m,count++);
						}
					}
				}
			}
		}
	}
}

#else
void hp_vrtx_bdry::non_sparse_snd(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi) {
    
    if (!base.is_comm()) return;
    
    const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
    
    /* Going to send all jacobian entries,  Diagonal entries for matching DOF's will be merged together not individual */
    /* Send number of non-zeros to matches */
    base.sndsize() = 0;
    base.sndtype() = boundary::int_msg;
    
    int pind = base.pnt*vdofs;
    for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
        base.isndbuf(base.sndsize()++) = nnzero(pind +*n);
    }
}


void hp_vrtx_bdry::non_sparse_rcv(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi) {
    
    if (!base.is_comm()) return;
    
    const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
    
    /* Reload to avoid overlap with sides */
    int pind = base.pnt*vdofs;
    int count = 0;
    for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
        nnzero(pind +*n) = base.isndbuf(count++);
        nnzero_mpi(pind +*n) = 0;
    }
    
    for(int m=0;m<base.nmatches();++m) {
        Array<int,1> target;
        if (base.is_local(m))  // only 1 matching boundary
            target.reference(nnzero);
        else
            target.reference(nnzero_mpi);
        
        int count = 0;
        int pind = base.pnt*vdofs;
        for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
            target(pind +*n) += base.ircvbuf(m,count++);
        }
    }
}


void hp_edge_bdry::non_sparse_snd(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi) {
	
	if (!base.is_comm()) return;
	
	const int sm=basis::tri(x.log2p)->sm();
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	const int begin_seg = x.npnt*vdofs;
	
	/* Going to send all jacobian entries,  Diagonal entries for matching DOF's will be merged together not individual */
	/* Send number of non-zeros to matches */
	base.sndsize() = 0;
	base.sndtype() = boundary::int_msg;
	for (int i=0;i<base.nseg;++i) {
		int sind = base.seg(i);
		int pind = x.seg(sind).pnt(0)*vdofs;
		for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
			base.isndbuf(base.sndsize()++) = nnzero(pind +*n);
		}
	}
	int sind = base.seg(base.nseg-1);
	int pind = x.seg(sind).pnt(1)*vdofs;
	for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
		base.isndbuf(base.sndsize()++) = nnzero(pind+*n);
	}
		
	/* Last thing to send is nnzero for edges */
	if (sm) {
		for (int i=0;i<base.nseg;++i) {
			int sind = base.seg(i);
			for (int m=0;m<sm;++m) {
				for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
					base.isndbuf(base.sndsize()++) = nnzero(begin_seg +sind*sm*x.NV +m*x.NV +*n);
				}
			}
		}
	}
	
	return;
}

void hp_edge_bdry::non_sparse_rcv(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi) {
	
	if (!base.is_comm()) return;
	
	const int sm=basis::tri(x.log2p)->sm();
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	const int begin_seg = x.npnt*vdofs;

	/* if local match then entries will go in to local jacobian */
	Array<int,1> target;
	if (base.is_local(0))  // only 1 matching boundary
		target.reference(nnzero);
	else
		target.reference(nnzero_mpi);
	
	int count = 0;
	for (int i=base.nseg-1;i>=0;--i) {
		int sind = base.seg(i);
		int pind = x.seg(sind).pnt(1)*vdofs;
		for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
			target(pind +*n) += base.ircvbuf(0,count++);
		}
	}
	int sind = base.seg(0);
	int pind = x.seg(sind).pnt(0)*vdofs;
	for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
		target(pind+*n) += base.ircvbuf(0,count++);
	}
	
	/* Now add to side degrees of freedom */
	if (sm) {
		for (int i=base.nseg-1;i>=0;--i) {
			int sind = base.seg(i);
			for (int mode=0;mode<sm;++mode) {
				for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
					target(begin_seg +sind*sm*x.NV +mode*x.NV +*n) += base.ircvbuf(0,count++);
				}
			}
		}
	}
}
#endif

#endif

void hp_vrtx_bdry::element_jacobian(Array<FLT,2>& K) {
	
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	Array<FLT,1> Rbar(vdofs),lf(vdofs);;
	
	/* Calculate and store initial residual */
	lf = 0.0;
	element_rsdl(lf);
	Rbar = lf;
	
	Array<FLT,1> dw(x.NV);
	dw = 0.0;
	for(int i=0;i<2;++i)
		for(int n=0;n<x.NV;++n)
			dw(n) = dw(n) + fabs(x.ug.v(base.pnt,n));
	
	dw = blitz::sum(dw)*eps_r;
	dw += eps_a;
	
	int kcol = 0;
	for(int var = 0; var < x.NV; ++var){
		x.ug.v(base.pnt,var) += dw(var);
		
		lf = 0.0;
		element_rsdl(lf);
		
		int krow = 0;
		for(int n=0;n<vdofs;++n)
			K(krow++,kcol) = (lf(n)-Rbar(n))/dw(var);
		
		++kcol;
		
		x.ug.v(base.pnt,var) -= dw(var);
	}

	int sind = x.ebdry(base.ebdry(1))->seg(0);
	FLT dx = eps_r*x.distance(x.seg(sind).pnt(0),x.seg(sind).pnt(1)) +eps_a;

	for(int var = 0; var < vdofs -x.NV; ++var){
		x.pnts(base.pnt)(var) += dx;
		
		lf = 0.0;
		element_rsdl(lf);
		
		int krow = 0;
		for(int n=0;n<vdofs;++n)
			K(krow++,kcol) = (lf(n)-Rbar(n))/dx;
		
		++kcol;
		
		x.pnts(base.pnt)(var) -= dx;
	}
}


void hp_edge_bdry::element_jacobian(int indx, Array<FLT,2>& K) {
	int sm = basis::tri(x.log2p)->sm();	
	Array<FLT,2> Rbar(x.NV,sm+2);
	Array<int,1> loc_to_glo(x.NV*(sm+2));
	
	/* Calculate and store initial residual */
	int sind = base.seg(indx);
	x.ugtouht1d(sind);
	
	x.lf = 0.0;
	element_rsdl(indx,x.lf);

	for(int k=0;k<sm+2;++k) {
		for(int n=0;n<x.NV;++n) {
			Rbar(n,k) = x.lf(n)(k);
		}
	}
	
	Array<FLT,1> dw(x.NV);
	dw = 0.0;
	for(int i=0;i<2;++i)
		for(int n=0;n<x.NV;++n)
			dw(n) = dw(n) + fabs(x.uht(n)(i));
	
	dw = blitz::sum(dw)*eps_r;
	dw += eps_a;
	
	
	if (x.mmovement != x.coupled_deformable) {
		/* Numerically create Jacobian */
		int kcol = 0;
		for(int mode = 0; mode < sm+2; ++mode){
			for(int var = 0; var < x.NV; ++var){
				x.uht(var)(mode) += dw(var);
				
				x.lf = 0.0;
				element_rsdl(indx,x.lf);
				
				int krow = 0;
				for(int k=0;k<sm+2;++k)
					for(int n=0;n<x.NV;++n)
						K(krow++,kcol) = (x.lf(n)(k)-Rbar(n,k))/dw(var);
						
						++kcol;
				
				x.uht(var)(mode) -= dw(var);
			}
		}
	}
	else {
		FLT dx = eps_r*x.distance(x.seg(sind).pnt(0),x.seg(sind).pnt(1)) +eps_a;

		/* Numerically create Jacobian */
		int kcol = 0;
		for(int mode = 0; mode < 2; ++mode) {
			for(int var = 0; var < x.NV; ++var){
				x.uht(var)(mode) += dw(var);
				
				x.lf = 0.0;
				element_rsdl(indx,x.lf);
				
				int krow = 0;
				for(int k=0;k<sm+2;++k)
					for(int n=0;n<x.NV;++n)
						K(krow++,kcol) = (x.lf(n)(k)-Rbar(n,k))/dw(var);
						
				++kcol;
				x.uht(var)(mode) -= dw(var);
			}
			
			for(int var = 0; var < tri_mesh::ND; ++var){
				x.pnts(x.seg(sind).pnt(mode))(var) += dx;
				
				x.lf = 0.0;
				element_rsdl(indx,x.lf);
				
				int krow = 0;
				for(int k=0;k<sm+2;++k)
					for(int n=0;n<x.NV;++n)
						K(krow++,kcol) = (x.lf(n)(k)-Rbar(n,k))/dx;
						
				++kcol;
				x.pnts(x.seg(sind).pnt(mode))(var) -= dx;
			}
		}
			
		for(int mode = 2; mode < sm+2; ++mode){
			for(int var = 0; var < x.NV; ++var){
				x.uht(var)(mode) += dw(var);
				
				x.lf = 0.0;
				element_rsdl(indx,x.lf);
				
				int krow = 0;
				for(int k=0;k<sm+2;++k)
					for(int n=0;n<x.NV;++n)
						K(krow++,kcol) = (x.lf(n)(k)-Rbar(n,k))/dw(var);
						
				++kcol;
				x.uht(var)(mode) -= dw(var);
			}
		}
	}
}

#ifdef petsc
int hp_edge_bdry::petsc_to_ug(PetscScalar *array) {
	int ind = 0;

	if (curved && coupled) {
		for(int j = 0; j < base.nseg;++j) {
			for(int m = 0; m < basis::tri(x.log2p)->sm(); ++m) {
				for(int n = 0; n < x.ND; ++n) { 
					crds(j,m,n) = array[ind++];
				}
			}
		}
	}
	
	return(ind);
}

void hp_edge_bdry::ug_to_petsc(int& ind) {
	if (curved && coupled) {
	
		for(int j = 0; j < base.nseg;++j) {
			for(int m = 0; m < basis::tri(x.log2p)->sm(); ++m) {
				for(int n = 0; n < x.ND; ++n) {
					VecSetValues(x.petsc_u,1,&ind,&crds(j,m,n),INSERT_VALUES);
					++ind;					
				}
			}
		}
	}
}

void hp_vrtx_bdry::petsc_jacobian() {
	
	/* Generic method for adding Jacobian terms */
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	Array<FLT,2> K(vdofs,vdofs);
	Array<int,1> rows(vdofs);
	Array<int,1> cols(vdofs);
	
	const int v0 = base.pnt;
	
	int rind = 0;
	int cind = 0;
	int gindx = vdofs*v0;
	for(int n=0;n<vdofs;++n) {
		rows(rind++) = gindx;
		cols(cind++) = gindx++;
	}

	element_jacobian(K);

#ifdef MY_SPARSE
	x.J.add_values(vdofs,rows,vdofs,cols,K);
#else
	MatSetValuesLocal(x.petsc_J,vdofs,rows.data(),vdofs,cols.data(),K.data(),ADD_VALUES);
#endif
}
	
void hp_edge_bdry::petsc_jacobian() {

	/* Generic method for adding Jacobian terms */
	const int sm = basis::tri(x.log2p)->sm();
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;

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
#ifdef MY_SPARSE
		x.J.add_values(x.NV*(sm+2),rows,vdofs*2 +x.NV*sm,cols,K);
#else
		MatSetValuesLocal(x.petsc_J,x.NV*(sm+2),rows.data(),vdofs*2 +x.NV*sm,cols.data(),K.data(),ADD_VALUES);
#endif
	}
}

#ifdef OLDWAY3
void hp_vrtx_bdry::petsc_matchjacobian_snd() {
	
	if (!base.is_comm())
		return;
	
	int vdofs;
	if (x.mmovement != x.coupled_deformable)
		vdofs = x.NV;
	else
		vdofs = x.NV+x.ND;
	
	int row;
	
#ifdef MY_SPARSE
	/* I am cheating here and sending floats and int's together */
	/* Send Jacobian entries */
	base.sndsize() = 0;
	base.sndtype() = boundary::flt_msg;
	/* First send number of entries for each vertex row */
	/* then append column numbers & values */
	row = base.pnt*vdofs; 
	/* attach diagonal column # to allow continuity enforcement */
	base.fsndbuf(base.sndsize()++) = row;
	base.fsndbuf(base.sndsize()++) = x.jacobian_start;
	for (int n=0;n<vdofs;++n) {
		base.fsndbuf(base.sndsize()++) = x.J._cpt(row+1) -x.J._cpt(row) +0.1;
#ifdef MPDEBUG
		*x.gbl->log << "sending " << x.J._cpt(row+1) -x.J._cpt(row) << " jacobian entries for vertex " << row/vdofs << " and variable " << n << std::endl;
#endif
		for (int col=x.J._cpt(row);col<x.J._cpt(row+1);++col) {
#ifdef MPDEBUG
			*x.gbl->log << x.J._col(col) << ' ';
#endif
			base.fsndbuf(base.sndsize()++) = x.J._col(col) +0.1;
			base.fsndbuf(base.sndsize()++) = x.J._val(col);
		}
#ifdef MPDEBUG
		*x.gbl->log << std::endl;
#endif
		++row;
	}
}

void hp_vrtx_bdry::petsc_matchjacobian_rcv(int phase)	{
	
	if (!base.is_comm() || base.matchphase(boundary::all_phased,0) != phase) return;
	
	sparse_row_major *pJ_mpi;

	int vdofs;
	if (x.mmovement != x.coupled_deformable)
		vdofs = x.NV;
	else
		vdofs = x.NV+x.ND;
	
	
	/* Reload values of Jacobian so don't have interference with sides */
	int count = 0;
	int row = base.pnt*vdofs;	

	assert(row == static_cast<int>(base.fsndbuf(count++)));
	assert(x.jacobian_start == static_cast<int>(base.fsndbuf(count++)));
	for(int n=0;n<vdofs;++n) {
		int ncol = static_cast<int>(base.fsndbuf(count++));
		for (int k = 0;k<ncol;++k) {
			int col = static_cast<int>(base.fsndbuf(count++));
			FLT val = base.fsndbuf(count++);
			/* This is a check that we aren't receiving unusued local degrees of freedom */
			if (col < INT_MAX-10 && col > -1) {
				x.J.set_values(row+n,col,val);
			}
		}
		x.J_mpi.multiply_row(row+n, 0.0);
	}
	
	
	for (int m=0;m<base.nmatches();++m) {
		count = 0;
		int row_mpi = static_cast<int>(base.frcvbuf(m,count++));
		int jstart_mpi = static_cast<int>(base.frcvbuf(m,count++));
		if (base.is_local(m)) {
			pJ_mpi = &x.J;
			jstart_mpi = 0;
		}
		else {
			pJ_mpi = &x.J_mpi;
		}
		row_mpi += jstart_mpi;
		
		for(int n=0;n<vdofs;++n) {
			int ncol = static_cast<int>(base.frcvbuf(m,count++));
#ifdef MPDEBUG
			*x.gbl->log << "receiving " << ncol << " jacobian entries for vertex " << row/vdofs << " and variable " << n << std::endl;
#endif
			for (int k = 0;k<ncol;++k) {
				int col = static_cast<int>(base.frcvbuf(m,count++));
				FLT val = base.frcvbuf(m,count++);
				
				if (col < INT_MAX-10 && col > -1) {
					col += jstart_mpi;
#ifdef MPDEBUG
					*x.gbl->log  << col << ' ';
#endif
					(*pJ_mpi).add_values(row+n,col,val);
				}
			}
#ifdef MPDEBUG
			*x.gbl->log << std::endl;
#endif
			/* Shift all entries for this vertex */
			for(int n_mpi=0;n_mpi<vdofs;++n_mpi) {
#ifndef SKIP_EXCHANGE
				FLT dval = (*pJ_mpi)(row+n,row_mpi+n_mpi);
				(*pJ_mpi)(row+n,row_mpi+n_mpi) = 0.0;				
				x.J(row+n,row+n_mpi) += dval;
#endif
			}
		}  
	}
	for(int n=0;n<vdofs;++n) {
		x.J.multiply_row(row+n,1.0/(base.nmatches()+1.0));
		x.J_mpi.multiply_row(row+n,1.0/(base.nmatches()+1.0));
	}
	
#else
	This Part not working
#endif
}

void hp_edge_bdry::petsc_matchjacobian_snd() {

	if (!base.is_comm())
		return;
		
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;

	int row,sind=-2;
	
#ifdef MY_SPARSE
	/* I am cheating here and sending floats and int's together */
	/* Send Jacobian entries */
	base.sndsize() = 0;
	base.sndtype() = boundary::flt_msg;
	base.fsndbuf(base.sndsize()++) = x.jacobian_start +0.1;
	
	/* First send number of entries for each vertex row */
	/* then append column numbers & values */
	for(int i=0;i<base.nseg;++i) {
		sind = base.seg(i);
		row = x.seg(sind).pnt(0)*vdofs; 
		/* attach diagonal column # to allow continuity enforcement */
		base.fsndbuf(base.sndsize()++) = row +0.1;;
		for (int n=0;n<vdofs;++n) {
			base.fsndbuf(base.sndsize()++) = x.J._cpt(row+1) -x.J._cpt(row) +0.1;
#ifdef MPDEBUG
			*x.gbl->log << "sending " << x.J._cpt(row+1) -x.J._cpt(row) << " jacobian entries for vertex " << row/vdofs << " and variable " << n << std::endl;
#endif
			for (int col=x.J._cpt(row);col<x.J._cpt(row+1);++col) {
#ifdef MPDEBUG
				*x.gbl->log << x.J._col(col) << ' ';
#endif
				base.fsndbuf(base.sndsize()++) = x.J._col(col) +0.1;
				base.fsndbuf(base.sndsize()++) = x.J._val(col);
			}
#ifdef MPDEBUG
			*x.gbl->log << std::endl;
#endif
			++row;
		}

		/* Send Side Information */
		row = x.npnt*vdofs +sind*x.NV*x.sm0;
		/* attach diagonal column # to allow continuity enforcement */
		base.fsndbuf(base.sndsize()++) = row +0.1;
		for(int mode=0;mode<x.sm0;++mode) {
			for (int n=0;n<x.NV;++n) {
				base.fsndbuf(base.sndsize()++) = x.J._cpt(row+1) -x.J._cpt(row) +0.1;
#ifdef MPDEBUG
				*x.gbl->log << "sending " << x.J._cpt(row+1) -x.J._cpt(row) << " jacobian entries for side " << sind << " and variable " << n << std::endl;
#endif
				for (int col=x.J._cpt(row);col<x.J._cpt(row+1);++col) {
#ifdef MPDEBUG
					*x.gbl->log << x.J._col(col) << ' ';
#endif
					base.fsndbuf(base.sndsize()++) = x.J._col(col) +0.1;
					base.fsndbuf(base.sndsize()++) = x.J._val(col);
				}

#ifdef MPDEBUG
				*x.gbl->log << std::endl;
#endif
				++row;
			}
		}
	}
	
	/* LAST POINT */
	row = x.seg(sind).pnt(1)*vdofs; 
	/* attach diagonal # to allow continuity enforcement */
	base.fsndbuf(base.sndsize()++) = row +0.1;
	for (int n=0;n<vdofs;++n) {
		base.fsndbuf(base.sndsize()++) = x.J._cpt(row+1) -x.J._cpt(row) +0.1;
#ifdef MPDEBUG
		*x.gbl->log << "sending " << x.J._cpt(row+1) -x.J._cpt(row) << " jacobian entries for vertex " << row/vdofs << " and variable " << n << std::endl;
#endif
		for (int col=x.J._cpt(row);col<x.J._cpt(row+1);++col) {
#ifdef MPDEBUG
			*x.gbl->log << x.J._col(col) << ' ';
#endif
			base.fsndbuf(base.sndsize()++) = x.J._col(col) +0.1;
			base.fsndbuf(base.sndsize()++) = x.J._val(col);
		}
#ifdef MPDEBUG
		*x.gbl->log << std::endl;
#endif
		++row;
	}
}
	
void hp_edge_bdry::petsc_matchjacobian_rcv(int phase)	{

	if (!base.is_comm() || base.matchphase(boundary::all_phased,0) != phase) return;
	
	int count = 0;
	int Jstart_mpi = static_cast<int>(base.frcvbuf(0, count++));
		
	sparse_row_major *pJ_mpi;
	if (base.is_local(0)) {
		pJ_mpi = &x.J;
		Jstart_mpi = 0;
	}
	else {
		pJ_mpi = &x.J_mpi;
	}
		
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;

	int row,sind;
	for (int i=base.nseg-1;i>=0;--i) {
		int sind = base.seg(i);
		int row = x.seg(sind).pnt(1)*vdofs; 
		int row_mpi = static_cast<int>(base.frcvbuf(0,count++)) +Jstart_mpi;
		for(int n=0;n<vdofs;++n) {
			int ncol = static_cast<int>(base.frcvbuf(0,count++));
#ifdef MPDEBUG
			*x.gbl->log << "receiving " << ncol << " jacobian entries for vertex " << row/vdofs << " and variable " << n << std::endl;
#endif
			for (int k = 0;k<ncol;++k) {
				int col = static_cast<int>(base.frcvbuf(0,count++));
				FLT val = base.frcvbuf(0,count++);
				if (col < INT_MAX-10 && col > -1) {
					 col += Jstart_mpi;
#ifdef MPDEBUG
					*x.gbl->log  << col << ' ';
#endif
					(*pJ_mpi).add_values(row+n,col,val);
				}
			}
#ifdef MPDEBUG
			*x.gbl->log << std::endl;
#endif
			/* Shift all entries for this vertex */
			for(int n_mpi=0;n_mpi<vdofs;++n_mpi) {
#ifdef MPDEBUG
				*x.gbl->log << "vertex swapping " << row+n << ',' << row_mpi+n_mpi << " for " << row+n << ',' << row+n_mpi << std::endl;
#endif
#ifndef SKIP_EXCHANGE
				FLT dval = (*pJ_mpi)(row+n,row_mpi+n_mpi);
				(*pJ_mpi)(row+n,row_mpi+n_mpi) = 0.0;				
				x.J(row+n,row+n_mpi) += dval;
#endif
			}
			x.J.multiply_row(row+n,0.5);
			x.J_mpi.multiply_row(row+n,0.5);
		}  
		
		/* Now receive side Jacobian information */
		row = x.npnt*vdofs +sind*x.NV*x.sm0;
		row_mpi = static_cast<int>(base.frcvbuf(0,count++)) +Jstart_mpi;
		int sgn = 1;
		int mcnt = 0;
		for(int mode=0;mode<x.sm0;++mode) {
			for(int n=0;n<x.NV;++n) {
				int ncol = static_cast<int>(base.frcvbuf(0,count++));
#ifdef MPDEBUG
				*x.gbl->log << "receiving " << ncol << " jacobian entries for side " << sind << " and variable " << n << std::endl;
#endif
				for (int k = 0;k<ncol;++k) {
					int col = static_cast<int>(base.frcvbuf(0,count++));
					FLT val = sgn*base.frcvbuf(0,count++);
					if (col < INT_MAX-10 && col > -1) {
						col += Jstart_mpi;
#ifdef MPDEBUG
						*x.gbl->log  << col << ' ';
#endif				
						(*pJ_mpi).add_values(row+mcnt,col,val);
					}
				}
#ifdef MPDEBUG
				*x.gbl->log << std::endl;
#endif
				/* Shift all modes in equation */
				int mcnt_mpi = 0;
				int sgn_mpi = 1;
				for(int mode_mpi=0;mode_mpi<x.sm0;++mode_mpi) {
					for(int n_mpi = 0;n_mpi<x.NV;++n_mpi) {
#ifdef MPDEBUG
						*x.gbl->log << "side swapping " << row+mcnt << ',' << row_mpi+mcnt_mpi << " for " << row+mcnt << ',' << row +mcnt_mpi << std::endl;
#endif
#ifndef SKIP_EXCHANGE
						FLT dval = (*pJ_mpi)(row+mcnt,row_mpi+mcnt_mpi);
						(*pJ_mpi)(row+mcnt,row_mpi+mcnt_mpi) = 0.0;				
						x.J(row+mcnt,row+mcnt_mpi) += sgn_mpi*dval;
#endif
						++mcnt_mpi;
					}
					sgn_mpi *= -1;
				}
				x.J.multiply_row(row+mcnt,0.5);
				x.J_mpi.multiply_row(row+mcnt,0.5);
				++mcnt;
			}
			sgn *= -1;
		}
	}
	sind = base.seg(0);
	row = x.seg(sind).pnt(0)*vdofs; 
	int row_mpi = static_cast<int>(base.frcvbuf(0,count++)) +Jstart_mpi;
	for(int n=0;n<vdofs;++n) {
		int ncol = static_cast<int>(base.frcvbuf(0,count++));
#ifdef MPDEBUG
		*x.gbl->log << "receiving " << ncol << " jacobian entries for vertex " << row/vdofs << " and variable " << n << std::endl;
#endif
		for (int k = 0;k<ncol;++k) {
			int col = static_cast<int>(base.frcvbuf(0,count++));
			FLT val = base.frcvbuf(0,count++);
			if (col < INT_MAX-10 && col > -1) {
				col += Jstart_mpi;
#ifdef MPDEBUG
				*x.gbl->log  << col << ' ';
#endif							
				(*pJ_mpi).add_values(row+n,col,val);
			}
		}
#ifdef MPDEBUG
		*x.gbl->log << std::endl;
#endif
		/* Shift all entries for this vertex */
		for(int n_mpi=0;n_mpi<vdofs;++n_mpi) {
#ifdef MPDEBUG
			*x.gbl->log << "vertex swapping " << row+n << ',' << row_mpi+n_mpi << " for " << row+n << ',' << row+n_mpi << std::endl;
#endif
#ifndef SKIP_EXCHANGE
			FLT dval = (*pJ_mpi)(row+n,row_mpi+n_mpi);
			(*pJ_mpi)(row+n,row_mpi+n_mpi) = 0.0;
			x.J(row+n,row+n_mpi) += dval;
#endif
		}
		x.J.multiply_row(row+n,0.5);
		x.J_mpi.multiply_row(row+n,0.5);
	} 
#else
    This Part not working
#endif
}
#else
void hp_vrtx_bdry::petsc_matchjacobian_snd() {
    
    if (!base.is_comm()) return;
    
    const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;

#ifdef MY_SPARSE
    /* I am cheating here and sending floats and int's together */
    /* Send Jacobian entries for continous variables  */
    base.sndsize() = 0;
    base.sndtype() = boundary::flt_msg;
    /* Send index of start of jacobian */
    base.fsndbuf(base.sndsize()++) = x.jacobian_start +0.1;
    /* Send index of start of surface unknowns (if they exist) */
    base.fsndbuf(base.sndsize()++) = jacobian_start+0.1;
    
    int rowbase = base.pnt*vdofs;
    /* Send continuous variables */
    for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
        int row = rowbase + *n;
        /* attach diagonal column # to allow continuity enforcement */
        base.fsndbuf(base.sndsize()++) = row +0.1;
        base.fsndbuf(base.sndsize()++) = x.J._cpt(row+1) -x.J._cpt(row) +0.1;
#ifdef MPDEBUG
        *x.gbl->log << "vertex sending " << x.J._cpt(row+1) -x.J._cpt(row) << " jacobian entries for vertex " << row/vdofs << " and variable " << *n << std::endl;
#endif
        for (int col=x.J._cpt(row);col<x.J._cpt(row+1);++col) {
#ifdef MPDEBUG
            *x.gbl->log << x.J._col(col) << ' ';
#endif
            base.fsndbuf(base.sndsize()++) = x.J._col(col) +0.1;
            base.fsndbuf(base.sndsize()++) = x.J._val(col);
        }
#ifdef MPDEBUG
        *x.gbl->log << std::endl;
#endif
    }
}

void hp_vrtx_bdry::petsc_matchjacobian_rcv(int phase)	{
	
	if (!base.is_comm() || base.matchphase(boundary::all_phased,0) != phase) return;
	
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	const int rowbase = base.pnt*vdofs;
	
	int count = 0;
	assert(x.jacobian_start == static_cast<int>(base.fsndbuf(count++)));
	assert(jacobian_start == static_cast<int>(base.fsndbuf(count++)));
	for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
		int row = static_cast<int>(base.fsndbuf(count++));
		int ncol = static_cast<int>(base.fsndbuf(count++));
		for (int k = 0;k<ncol;++k) {
			int col = static_cast<int>(base.fsndbuf(count++));
			FLT val = base.fsndbuf(count++);
			/* This is a check that we aren't receiving unusued local degrees of freedom */
			if (col < INT_MAX-10 && col > -1) {
				x.J.set_values(row,col,val);
			}
		}
		x.J_mpi.multiply_row(row, 0.0);
	}
	
	for (int m=0;m<base.nmatches();++m) {
		int count = 0;
		int Jstart_mpi = static_cast<int>(base.frcvbuf(m, count++)); // Start of jacobian on matching block
		int Jstart_mpi_vrtx_unknowns = static_cast<int>(base.frcvbuf(m, count++)); // Start of vertex unknowns on mathcing block (not used typically)
		
		sparse_row_major *pJ_mpi;
    if (base.is_local(m)) {
        pJ_mpi = &x.J;
        Jstart_mpi = 0;
    }
    else {
        pJ_mpi = &x.J_mpi;
    }
    
		vector<int> row_mpi_storage;
		for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
			int row = rowbase + *n;
			int row_mpi = static_cast<int>(base.frcvbuf(m,count++)) +Jstart_mpi;
			row_mpi_storage.push_back(row_mpi);
			int ncol = static_cast<int>(base.frcvbuf(m,count++));
#ifdef MPDEBUG
			*x.gbl->log << "vertex receiving " << ncol << " jacobian entries for vertex " << row/vdofs << " and variable " << *n << std::endl;
#endif
			for (int k = 0;k<ncol;++k) {
				int col = static_cast<int>(base.frcvbuf(m,count++));
				FLT val = base.frcvbuf(m,count++);
				if (col < INT_MAX-10 && col > -1) {
					col += Jstart_mpi;
#ifdef MPDEBUG
					*x.gbl->log  << col << ' ';
#endif
					(*pJ_mpi).add_values(row,col,val);
				}
			}
#ifdef MPDEBUG
			*x.gbl->log << std::endl;
#endif
		}
		
		/* Shift all diagonal block entries for this vertex */
		for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
			int row = rowbase + *n;
			std::vector<int>::iterator row_mpi=row_mpi_storage.begin();
			for(std::vector<int>::iterator n_mpi=c0_indices_xy.begin();n_mpi != c0_indices_xy.end();++n_mpi) {
#ifdef MPDEBUG
				*x.gbl->log << "vertex swapping " << row << ',' << *row_mpi << " for " << row << ',' << rowbase +*n_mpi << std::endl;
#endif
#ifndef SKIP_EXCHANGE
				FLT dval = (*pJ_mpi)(row,*row_mpi);
				(*pJ_mpi)(row,*row_mpi) = 0.0;
				x.J(row,rowbase+*n_mpi) += dval;
#endif
				++row_mpi;
			}
		}
		row_mpi_storage.clear();
	}
	
	for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
		x.J.multiply_row(rowbase+*n,1.0/(1.0+base.nmatches()));
		x.J_mpi.multiply_row(rowbase+*n,1.0/(1.0+base.nmatches()));
	}
	
#else
    This Part not working
#endif
}

void hp_edge_bdry::petsc_matchjacobian_snd() {
	
	if (!base.is_comm()) return;
	
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	
	/* Now do stuff for communication boundaries */
	int row,sind=-2;
	
	/* I am cheating here and sending floats and int's together */
#ifdef MY_SPARSE
	/* Send Jacobian entries for continous variables  */
	base.sndsize() = 0;
	base.sndtype() = boundary::flt_msg;
	/* Send index of start of jacobian */
	base.fsndbuf(base.sndsize()++) = x.jacobian_start +0.1;
	/* Send index of start of surface unknowns (if they exist) */
	base.fsndbuf(base.sndsize()++) = jacobian_start+0.1;
	
	for(int i=0;i<base.nseg;++i) {
		sind = base.seg(i);
		int rowbase = x.seg(sind).pnt(0)*vdofs;
		
		/* Send continuous variables */
		for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
			row = rowbase + *n;
			/* attach diagonal column # to allow continuity enforcement */
			base.fsndbuf(base.sndsize()++) = row +0.1;
			base.fsndbuf(base.sndsize()++) = x.J._cpt(row+1) -x.J._cpt(row) +0.1;
#ifdef MPDEBUG
			*x.gbl->log << "sending " << x.J._cpt(row+1) -x.J._cpt(row) << " jacobian entries for vertex " << row/vdofs << " and variable " << *n << std::endl;
#endif
			for (int col=x.J._cpt(row);col<x.J._cpt(row+1);++col) {
#ifdef MPDEBUG
				*x.gbl->log << x.J._col(col) << ' ';
#endif
				base.fsndbuf(base.sndsize()++) = x.J._col(col) +0.1;
				base.fsndbuf(base.sndsize()++) = x.J._val(col);
			}
#ifdef MPDEBUG
			*x.gbl->log << std::endl;
#endif
		}
		
		/* Send Side Information */
		rowbase = x.npnt*vdofs +sind*x.NV*x.sm0;
		for(int mode=0;mode<x.sm0;++mode) {
			for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
				row = rowbase + *n;
				/* attach diagonal column # to allow continuity enforcement */
				base.fsndbuf(base.sndsize()++) = row +0.1;
				base.fsndbuf(base.sndsize()++) = x.J._cpt(row+1) -x.J._cpt(row) +0.1;
#ifdef MPDEBUG
				*x.gbl->log << "sending " << x.J._cpt(row+1) -x.J._cpt(row) << " jacobian entries for side " << sind << " and variable " << *n << std::endl;
#endif
				for (int col=x.J._cpt(row);col<x.J._cpt(row+1);++col) {
#ifdef MPDEBUG
					*x.gbl->log << x.J._col(col) << ' ';
#endif
					base.fsndbuf(base.sndsize()++) = x.J._col(col) +0.1;
					base.fsndbuf(base.sndsize()++) = x.J._val(col);
				}
				
#ifdef MPDEBUG
				*x.gbl->log << std::endl;
#endif
			}
			rowbase += x.NV;
		}
	}
	
	/* LAST POINT */
	int rowbase = x.seg(sind).pnt(1)*vdofs;
	for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
		row = rowbase + *n;
		/* attach diagonal # to allow continuity enforcement */
		base.fsndbuf(base.sndsize()++) = row +0.1;
		base.fsndbuf(base.sndsize()++) = x.J._cpt(row+1) -x.J._cpt(row) +0.1;
#ifdef MPDEBUG
		*x.gbl->log << "sending " << x.J._cpt(row+1) -x.J._cpt(row) << " jacobian entries for vertex " << row/vdofs << " and variable " << *n << std::endl;
#endif
		for (int col=x.J._cpt(row);col<x.J._cpt(row+1);++col) {
#ifdef MPDEBUG
			*x.gbl->log << x.J._col(col) << ' ';
#endif
			base.fsndbuf(base.sndsize()++) = x.J._col(col) +0.1;
			base.fsndbuf(base.sndsize()++) = x.J._val(col);
		}
		
#ifdef MPDEBUG
		*x.gbl->log << std::endl;
#endif
	}
}


void hp_edge_bdry::petsc_matchjacobian_rcv(int phase) {
	
	if (!base.is_comm() || base.matchphase(boundary::all_phased,0) != phase) return;
	
	int count = 0;
	int Jstart_mpi = static_cast<int>(base.frcvbuf(0, count++)); // Start of jacobian on matching block
    count++; // Skip index of boundary unknowns on mathcing block.   Not used in this routine
	
	sparse_row_major *pJ_mpi;
	if (base.is_local(0)) {
		pJ_mpi = &x.J;
		Jstart_mpi = 0;
	}
	else {
		pJ_mpi = &x.J_mpi;
	}
	
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	
	/* Now do stuff for communication boundaries */
	int row;
	/* Now Receive Information */
	for (int i=base.nseg-1;i>=0;--i) {
		int sind = base.seg(i);
		int rowbase = x.seg(sind).pnt(1)*vdofs;
		vector<int> row_mpi_storage;
		for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
			row = rowbase + *n;
			int row_mpi = static_cast<int>(base.frcvbuf(0,count++)) +Jstart_mpi;
			row_mpi_storage.push_back(row_mpi);
			int ncol = static_cast<int>(base.frcvbuf(0,count++));
#ifdef MPDEBUG
			*x.gbl->log << "receiving " << ncol << " jacobian entries for vertex " << row/vdofs << " and variable " << *n << std::endl;
#endif
			for (int k = 0;k<ncol;++k) {
				int col = static_cast<int>(base.frcvbuf(0,count++));
				FLT val = base.frcvbuf(0,count++);
				if (col < INT_MAX-10 && col > -1) {
					col += Jstart_mpi;
#ifdef MPDEBUG
					*x.gbl->log  << col << ' ';
#endif
					(*pJ_mpi).add_values(row,col,val);
				}
			}
#ifdef MPDEBUG
			*x.gbl->log << std::endl;
#endif
		}
		
		/* Shift all diagonal block entries for this vertex */
		for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
			row = rowbase + *n;
			std::vector<int>::iterator row_mpi=row_mpi_storage.begin();
			for(std::vector<int>::iterator n_mpi=c0_indices_xy.begin();n_mpi != c0_indices_xy.end();++n_mpi) {
#ifdef MPDEBUG
				*x.gbl->log << "vertex swapping " << row << ',' << *row_mpi << " for " << row << ',' << rowbase +*n_mpi << std::endl;
#endif
#ifndef SKIP_EXCHANGE
				FLT dval = (*pJ_mpi)(row,*row_mpi);
				(*pJ_mpi)(row,*row_mpi) = 0.0;
				x.J(row,rowbase+*n_mpi) += dval;
#endif
				++row_mpi;
			}
			x.J.multiply_row(row,0.5);
			x.J_mpi.multiply_row(row,0.5);
		}
		row_mpi_storage.clear();
		
		/* Now receive side Jacobian information */
		rowbase = x.npnt*vdofs +sind*x.NV*x.sm0;
		int sgn = 1;
		for(int mode=0;mode<x.sm0;++mode) {
			for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
				row = rowbase +*n;
				int row_mpi = static_cast<int>(base.frcvbuf(0,count++)) +Jstart_mpi;
				row_mpi_storage.push_back(row_mpi);
				int ncol = static_cast<int>(base.frcvbuf(0,count++));
#ifdef MPDEBUG
				*x.gbl->log << "receiving " << ncol << " jacobian entries for side " << sind << " and variable " << *n << std::endl;
#endif
				for (int k = 0;k<ncol;++k) {
					int col = static_cast<int>(base.frcvbuf(0,count++));
					FLT val = sgn*base.frcvbuf(0,count++);
					if (col < INT_MAX-10 && col > -1) {
						col += Jstart_mpi;
#ifdef MPDEBUG
						*x.gbl->log  << col << ' ';
#endif
						(*pJ_mpi).add_values(row,col,val);
					}
				}
#ifdef MPDEBUG
				*x.gbl->log << std::endl;
#endif
			}
			rowbase += x.NV;
			sgn *= -1;
		}
		
		/* Shift diagonal block of modes */
		rowbase = x.npnt*vdofs +sind*x.NV*x.sm0;
		for(int mode=0;mode<x.sm0;++mode) {
			for(std::vector<int>::iterator n=c0_indices.begin();n != c0_indices.end();++n) {
				row = rowbase +mode*x.NV +*n;
				int sgn_mpi = 1;
				std::vector<int>::iterator row_mpi=row_mpi_storage.begin();
				for(int mode_mpi=0;mode_mpi<x.sm0;++mode_mpi) {
					for(std::vector<int>::iterator n_mpi=c0_indices.begin();n_mpi != c0_indices.end();++n_mpi) {
#ifdef MPDEBUG
						*x.gbl->log << "side swapping " << row << ',' << *row_mpi << " for " << row << ',' << rowbase +mode_mpi*x.NV +*n_mpi << std::endl;
#endif
#ifndef SKIP_EXCHANGE
						FLT dval = (*pJ_mpi)(row,*row_mpi);
						(*pJ_mpi)(row,*row_mpi) = 0.0;
						x.J(row,rowbase +mode_mpi*x.NV +*n_mpi) += sgn_mpi*dval;
#endif
						++row_mpi;
					}
					sgn_mpi *= -1;
				}
				x.J.multiply_row(row,0.5);
				x.J_mpi.multiply_row(row,0.5);
			}
		}
		row_mpi_storage.clear();
	}
	int sind = base.seg(0);
	int rowbase = x.seg(sind).pnt(0)*vdofs;
	vector<int> row_mpi_storage;
	for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
		row = rowbase + *n;
		int row_mpi = static_cast<int>(base.frcvbuf(0,count++)) +Jstart_mpi;
		row_mpi_storage.push_back(row_mpi);
		int ncol = static_cast<int>(base.frcvbuf(0,count++));
#ifdef MPDEBUG
		*x.gbl->log << "receiving " << ncol << " jacobian entries for vertex " << row/vdofs << " and variable " << *n << std::endl;
#endif
		for (int k = 0;k<ncol;++k) {
			int col = static_cast<int>(base.frcvbuf(0,count++));
			FLT val = base.frcvbuf(0,count++);
			if (col < INT_MAX-10 && col > -1) {
				col += Jstart_mpi;
#ifdef MPDEBUG
				*x.gbl->log  << col << ' ';
#endif
				(*pJ_mpi).add_values(row,col,val);
			}
		}
#ifdef MPDEBUG
		*x.gbl->log << std::endl;
#endif
	}
	
	/* Shift all diagonal block entries for this vertex */
	for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
		row = rowbase + *n;
		std::vector<int>::iterator row_mpi=row_mpi_storage.begin();
		for(std::vector<int>::iterator n_mpi=c0_indices_xy.begin();n_mpi != c0_indices_xy.end();++n_mpi) {
#ifdef MPDEBUG
			*x.gbl->log << "vertex swapping " << row << ',' << *row_mpi << " for " << row << ',' << rowbase +*n_mpi << std::endl;
#endif
#ifndef SKIP_EXCHANGE
			FLT dval = (*pJ_mpi)(row,*row_mpi);
			(*pJ_mpi)(row,*row_mpi) = 0.0;
			x.J(row,rowbase+*n_mpi) += dval;
#endif
			++row_mpi;
		}
		x.J.multiply_row(row,0.5);
		x.J_mpi.multiply_row(row,0.5);
	}
	row_mpi_storage.clear();
#endif
}
#endif
#endif



void hp_edge_bdry::element_rsdl(int eind, Array<TinyVector<FLT,MXTM>,1> lf) {
	int k,n,sind;
	TinyVector<FLT,2> pt,mvel,nrm;
	Array<FLT,1> u(x.NV),flx(x.NV);
	
	lf = 0.0;
	sind = base.seg(eind);
	
	/* Load coordinates */
	x.crdtocht1d(sind);
	
	/* Project coordinates to Gauss points & coordinate derivates */
	for(n=0;n<tri_mesh::ND;++n)
		basis::tri(x.log2p)->proj1d(&x.cht(n,0),&x.crd(n)(0,0),&x.dcrd(n,0)(0,0));
		
	/* Project solution to Gauss points */
	for(n=0;n<x.NV;++n)
		basis::tri(x.log2p)->proj1d(&x.uht(n)(0),&x.u(n)(0,0));
	
	/* Integrate flux across element boundary times basis functions */
	for(k=0;k<basis::tri(x.log2p)->gpx();++k) {
		/* Calculate boundary normal */
		nrm(0) = x.dcrd(1,0)(0,k);
		nrm(1) = -x.dcrd(0,0)(0,k);
		FLT jcb = sqrt(nrm(0)*nrm(0) +nrm(1)*nrm(1));
		nrm(0) /= jcb;
		nrm(1) /= jcb;

		/* Calculate the mesh velocity */
		for(n=0;n<tri_mesh::ND;++n) {
			pt(n) = x.crd(n)(0,k);
			mvel(n) = x.gbl->bd(0)*(x.crd(n)(0,k) -dxdt(x.log2p,eind)(n,k));
#ifdef MESH_REF_VEL
			mvel(n) += x.gbl->mesh_ref_vel(n);
#endif
			
		}
		
		for(n=0;n<x.NV;++n)
			u(n) = x.u(n)(0,k);
		
		/* Call flux function */
		flux(u,pt,mvel,nrm,jcb,flx);
		
		for(n=0;n<x.NV;++n)
			x.res(n)(0,k) = RAD(x.crd(0)(0,k))*flx(n)*jcb;
	}
	/* Integrate fluxes * basis */
	for(n=0;n<x.NV;++n)
		basis::tri(x.log2p)->intgrt1d(&lf(n)(0),&x.res(n)(0,0));
	
	return;
}

void hp_edge_bdry::flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx) {
	
	Array<FLT,1> axpt(tri_mesh::ND), amv(tri_mesh::ND), anorm(tri_mesh::ND);
	axpt(0) = xpt(0); axpt(1) = xpt(1);
	amv(0) = mv(0); amv(1) = mv(1);
	anorm(0)= norm(0); anorm(1) = norm(1);
	
	for (int n=0;n<x.NV;++n) {
		switch(type[n]) {
			case(essential): {
				flx(n) = 0.0;
				break;
			}
			case(natural): {
				flx(n) = fluxes[n].Eval(u,axpt,amv,anorm,x.gbl->time);
				break;
			}
		}
	}
	
	return;
}




void hp_edge_bdry::vdirichlet() {
	int sind,j,v0;
	
	j = 0;
	do {
		sind = base.seg(j);
		v0 = x.seg(sind).pnt(0);
		for(std::vector<int>::iterator n=essential_indices.begin();n != essential_indices.end();++n) 
			x.gbl->res.v(v0,*n)= 0.0;
	} while (++j < base.nseg);
	v0 = x.seg(sind).pnt(1);
	for(std::vector<int>::iterator n=essential_indices.begin();n != essential_indices.end();++n) 
		x.gbl->res.v(v0,*n) = 0.0;
}

void hp_edge_bdry::sdirichlet(int mode) {
	int sind;
	
	for(int j=0;j<base.nseg;++j) {
		sind = base.seg(j);
		for(std::vector<int>::iterator n=essential_indices.begin();n != essential_indices.end();++n) 
			x.gbl->res.s(sind,mode,*n) = 0.0;
	}
}

#ifdef petsc
void hp_vrtx_bdry::petsc_jacobian_dirichlet() {
	const int nessentials = essential_indices.size();
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	const int v0 = base.pnt;
	const int gind = v0*vdofs;
	Array<int,1> indices(nessentials);

	int counter = 0;
	for(std::vector<int>::const_iterator n=essential_indices.begin();n != essential_indices.end();++n) {
		indices(counter++)=gind +*n;
	}
	
#ifdef MY_SPARSE
	x.J.zero_rows(counter,indices);
	x.J_mpi.zero_rows(counter,indices);
	x.J.set_diag(counter,indices,1.0);
#else
	MatZeroRows(x.petsc_J,counter,indices.data(),1.0);
#endif
}

void hp_edge_bdry::petsc_jacobian_dirichlet() {	
	int sm=basis::tri(x.log2p)->sm();
	int nessentials = essential_indices.size();
	Array<int,1> indices((base.nseg+1)*nessentials +base.nseg*sm*nessentials);
	
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;

	
	/* only works if pressure is 4th variable */
	int gind,v0,sind;
	int counter = 0;
	
	int j = 0;
	do {
		sind = base.seg(j);
		v0 = x.seg(sind).pnt(0);
		gind = v0*vdofs;
		for(std::vector<int>::const_iterator n=essential_indices.begin();n != essential_indices.end();++n) {
			indices(counter++)=gind +*n;
		}
	} while (++j < base.nseg);
	v0 = x.seg(sind).pnt(1);
	gind = v0*vdofs;
	for(std::vector<int>::iterator n=essential_indices.begin();n != essential_indices.end();++n) {
		indices(counter++)=gind +*n;
	}
	
	for(int i=0;i<base.nseg;++i) {
		gind = x.npnt*vdofs+base.seg(i)*sm*x.NV;
		for(int m=0; m<sm; ++m) {
			for(std::vector<int>::iterator n=essential_indices.begin();n != essential_indices.end();++n) {
				indices(counter++)=gind+m*x.NV +*n;
			}
		}
	}	
	
#ifdef MY_SPARSE
	x.J.zero_rows(counter,indices);
	x.J_mpi.zero_rows(counter,indices);
	x.J.set_diag(counter,indices,1.0);
#else
	MatZeroRows(x.petsc_J,counter,indices.data(),1.0);
#endif
}
#endif


void symbolic_with_integration_by_parts::init(input_map& inmap,void* gbl_in) {
	std::string keyword;
	std::ostringstream nstr;
	input_map zeromap;
	zeromap["zero"] = "0.0";
	
	hp_edge_bdry::init(inmap,gbl_in);
	for (int n=0;n<x.NV;++n) {
		if (type[n] == natural) {
			nstr.str("");
			nstr << base.idprefix << "_dflux" << n << std::flush;
			if (inmap.find(nstr.str()) != inmap.end()) {
				derivative_fluxes[n].init(inmap,nstr.str());
			}
			else {
				derivative_fluxes[n].init(zeromap,"zero");
			}
		}
	}
	
	return;
}

void symbolic_with_integration_by_parts::element_rsdl(int indx, Array<TinyVector<FLT,MXTM>,1> lf) {
	int i,n,sind,v0,v1;
	TinyVector<FLT,tri_mesh::ND> norm, rp;
	Array<FLT,1> ubar(x.NV),flx(x.NV);
	FLT jcb;
	Array<TinyVector<FLT,MXGP>,1> u(x.NV);
	TinyMatrix<FLT,tri_mesh::ND,MXGP> crd, dcrd, mvel;
	Array<FLT,2> cflux(x.NV,MXGP),dflux(x.NV,MXGP);
	Array<FLT,1> axpt(tri_mesh::ND), amv(tri_mesh::ND), anorm(tri_mesh::ND),au(x.NV);
	
	lf = 0.0;
	
	sind = base.seg(indx);
	v0 = x.seg(sind).pnt(0);
	v1 = x.seg(sind).pnt(1);
	
	x.crdtocht1d(sind);
	for(n=0;n<tri_mesh::ND;++n)
		basis::tri(x.log2p)->proj1d(&x.cht(n,0),&crd(n,0),&dcrd(n,0));
	
	for(n=0;n<x.NV;++n)
		basis::tri(x.log2p)->proj1d(&x.uht(n)(0),&u(n)(0));    
	
	for(i=0;i<basis::tri(x.log2p)->gpx();++i) {
		norm(0) =  dcrd(1,i);
		norm(1) = -dcrd(0,i);
		jcb = sqrt(norm(0)*norm(0) +norm(1)*norm(1));
		
		/* RELATIVE VELOCITY STORED IN MVEL(N)*/
		for(n=0;n<tri_mesh::ND;++n) {
			mvel(n,i) = x.gbl->bd(0)*(crd(n,i) -dxdt(x.log2p,indx)(n,i));
#ifdef MESH_REF_VEL
			mvel(n,i) += x.gbl->mesh_ref_vel(n);
#endif
		}
		
		/* Evaluate Fluxes */
		axpt(0) = crd(0,i); axpt(1) = crd(1,i);
		amv(0) = mvel(0,i); amv(1) = mvel(1,i);
		anorm(0)= norm(0)/jcb; anorm(1) = norm(1)/jcb;
		for(int n=0;n<x.NV;++n)
			au(n) = u(n)(i);
		
		for(int n=0;n<x.NV;++n) {
			switch(type[n]) {
				case(essential): {
					cflux(n,i) = 0.0;
					dflux(n,i) = 0.0;
					break;
				}
				case(natural): {
					cflux(n,i) = RAD(crd(0,i))*fluxes[n].Eval(au,axpt,amv,anorm,x.gbl->time)*jcb;
					dflux(n,i) = RAD(crd(0,i))*derivative_fluxes[n].Eval(au,axpt,amv,anorm,x.gbl->time);
					break;
				}
			}
		}
	}
	
	for(int n=0;n<x.NV;++n) {
		basis::tri(x.log2p)->intgrt1d(&lf(n)(0),&cflux(n,0));
		basis::tri(x.log2p)->intgrtx1d(&lf(n)(0),&dflux(n,0));
	}
	
	return;
}
