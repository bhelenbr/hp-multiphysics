/*
 *  output.cpp
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 15 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp.h"
#include "hp_boundary.h"
#include <cstring>
#include <assert.h>
#include <stdlib.h>
#include <utilities.h>
#include <myblas.h>
#include <libbinio/binwrap.h>
#include <libbinio/binfile.h>
//#define DATATANK

#ifdef DATATANK
#include <DTSource.h>
#endif


void tri_hp::output(const std::string& fname, block::output_purpose why) {
	int i,j;
	std::string fnmapp, namewdot;
	std::ostringstream nstr;
	ofstream out;
	out.setf(std::ios::scientific, std::ios::floatfield);
	out.precision(4);

	switch(why) {
		case(block::display): {
			output(fname,output_type(0));
			helper->output();
			return;
		}
		case(block::restart): {
			namewdot = fname +"_d";
			for(i=0;i<gbl->nadapt;++i) {
				nstr.str("");
				nstr << i << std::flush;
				fnmapp = namewdot +nstr.str();
				output(fnmapp,output_type(1),i);
			}
			if (mmovement != fixed || gbl->adapt_interval) {
				namewdot = fname +"_v";
				if (output_type(1) == tri_hp::binary) {
					tri_mesh::output(fname,tri_mesh::binary);
					binofstream bout;
					for(i=1;i<gbl->nadapt;++i) {
						nstr.str("");
						nstr << i << std::flush;
						fnmapp = namewdot +nstr.str() +"_" +gbl->idprefix +".bin";
						bout.open(fnmapp.c_str());
						bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::BigEndian)),1);
						bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::FloatIEEE)),1);
						for (j=0;j<npnt;++j) { 
							bout.writeFloat(vrtxbd(i)(j)(0),binio::Double);
							bout.writeFloat(vrtxbd(i)(j)(1),binio::Double);
						}
						bout.close();
					}
				}
				else {
					tri_mesh::output(fname,tri_mesh::grid);
					tri_mesh::output(fname,tri_mesh::vlength);
					for(i=1;i<gbl->nadapt;++i) {
						nstr.str("");
						nstr << i << std::flush;
						fnmapp = namewdot +nstr.str() +"_" +gbl->idprefix +".txt";
						out.open(fnmapp.c_str());
						out << npnt << std::endl;
						for (j=0;j<npnt;++j) 
						out << j << ':' << vrtxbd(i)(j)(0) << ' ' << vrtxbd(i)(j)(1) << std::endl;
						out.close();
					}
				}
			}			
			return;
		}
		case(block::debug): {
			output(fname,output_type(2));
			helper->output();
		}
	}
	return;
}

 void tri_hp::output(const std::string& filename, filetype typ, int tlvl) {
	ofstream out;
	std::string fname, fnmapp;
	int i,j,k,m,n,v0,v1,sind,tind,indx,sgn;
	int ijind[MXTM][MXTM];

	out.setf(std::ios::scientific, std::ios::floatfield);
	out.precision(8);
	 
	fname = filename +"_" +gbl->idprefix;

	switch (typ) {
		case (text):
			fnmapp = fname +".txt";
			out.open(fnmapp.c_str());
			if (!out) {
				*gbl->log << "couldn't open text output file " << fnmapp;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}

			/* HEADER INFORMATION */
			out << "p0 = " << p0 << std::endl;
			out << "npnt = " << npnt << ", nseg = " << nseg << ", ntri = " << ntri << std::endl;
			out << "END OF HEADER" << std::endl;

			for(i=0;i<npnt;++i) {
				for(n=0;n<NV;++n)
					out << ugbd(tlvl).v(i,n) << '\t';
				out << std::endl;
			}

			for(i=0;i<nseg;++i) {
				for(m=0;m<sm0;++m) {
					for(n=0;n<NV;++n)
						out << ugbd(tlvl).s(i,m,n) << '\t';
					out << std::endl;
				}
			}

			for(i=0;i<ntri;++i) {
				for(m=0;m<im0;++m) {
					for(n=0;n<NV;++n)
						out << ugbd(tlvl).i(i,m,n) << '\t';
					out << std::endl;
				}
			}
			out.close();

			break;

		case (binary): {
			fnmapp = fname +".bin";
			out.open(fnmapp.c_str(),std::ios::binary);
			if (!out) {
				*gbl->log << "couldn't open text output file " << fnmapp;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			binowstream bout(&out);
			bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::BigEndian)),1);
			bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::FloatIEEE)),1);

			/* HEADER INFORMATION */
			bout.writeInt(p0,sizeof(int));
			bout.writeInt(npnt,sizeof(int));
			bout.writeInt(nseg,sizeof(int));
			bout.writeInt(ntri,sizeof(int));

			for(i=0;i<npnt;++i) {
				for(n=0;n<NV;++n)
					bout.writeFloat(ugbd(tlvl).v(i,n),binio::Double);
			}

			for(i=0;i<nseg;++i) {
				for(m=0;m<sm0;++m) {
					for(n=0;n<NV;++n)
						bout.writeFloat(ugbd(tlvl).s(i,m,n),binio::Double);
				}
			}

			for(i=0;i<ntri;++i) {
				for(m=0;m<im0;++m) {
					for(n=0;n<NV;++n)
						bout.writeFloat(ugbd(tlvl).i(i,m,n),binio::Double);
				}
			}
			out.close();

			break;
		}

		case (vtk):
			fnmapp = fname +".vtk";
			out.open(fnmapp.c_str());
			if (!out) {
				*gbl->log << "couldn't open vtk output file " << fnmapp;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			
			out << "# vtk DataFile Version 1.0" << endl;
			out << "add title here" << endl;
			out << "ASCII" << endl << endl;
			
			out << "DATASET UNSTRUCTURED_GRID" << endl;
			out << "POINTS " << npnt+basis::tri(log2p)->sm()*nseg+basis::tri(log2p)->im()*ntri << " float" << endl;
			
			for(i=0;i<npnt;++i) {
				for(n=0;n<ND;++n)
					out << vrtxbd(tlvl)(i)(n) << ' ';
				out << 0.0 << endl;
			}
			
			if (basis::tri(log2p)->p() > 1) {
				/* SIDE MODES */
				for(sind=0;sind<nseg;++sind) {
					if (seg(sind).info < 0) {
						v0 = seg(sind).pnt(0);
						v1 = seg(sind).pnt(1);
						for(n=0;n<ND;++n)
							basis::tri(log2p)->proj1d_leg(vrtxbd(tlvl)(v0)(n),vrtxbd(tlvl)(v1)(n),&crd(n)(0,0));
					}
					else {
						crdtocht1d(sind,tlvl);
						
						for(n=0;n<ND;++n)
							basis::tri(log2p)->proj1d_leg(&cht(n,0),&crd(n)(0,0));
					}
					
					for(i=1;i<basis::tri(log2p)->sm()+1;++i) {
						for(n=0;n<ND;++n)
							out << crd(n)(0,i) << ' ';                    
						out << 0.0 << endl;
					}
				}
				
				/* INTERIOR MODES */
				if (basis::tri(log2p)->p() > 2) {
					for(tind = 0; tind < ntri; ++tind) {
						if (tri(tind).info < 0) {
							for(n=0;n<ND;++n)
								basis::tri(log2p)->proj_leg(vrtxbd(tlvl)(tri(tind).pnt(0))(n),vrtxbd(tlvl)(tri(tind).pnt(1))(n),vrtxbd(tlvl)(tri(tind).pnt(2))(n),&crd(n)(0,0),MXGP);
						}
						else {
							crdtocht(tind,tlvl);
							for(n=0;n<ND;++n)
								basis::tri(log2p)->proj_bdry_leg(&cht(n,0),&crd(n)(0,0),MXGP);
						}
						
						for(i=1;i<basis::tri(log2p)->sm();++i) {
							for(j=1;j<basis::tri(log2p)->sm()-(i-1);++j) {
								for(n=0;n<ND;++n)
									out << crd(n)(i,j) << ' ';
								out << 0.0 << endl;

							}
						}
					}
				}
			}
			
			out << "CELLS " << ntri*(basis::tri(log2p)->sm()+1)*(basis::tri(log2p)->sm()+1) << ' ' << 4*ntri*(basis::tri(log2p)->sm()+1)*(basis::tri(log2p)->sm()+1) << endl;
			
			for(tind=0;tind<ntri;++tind) {
				
				/* VERTICES */
				ijind[0][basis::tri(log2p)->sm()+1] = tri(tind).pnt(0);
				ijind[0][0] = tri(tind).pnt(1);
				ijind[basis::tri(log2p)->sm()+1][0] = tri(tind).pnt(2);
				
				/* SIDES */
				indx = tri(tind).seg(0);
				sgn = tri(tind).sgn(0);
				if (sgn < 0) {
					for(i=0;i<basis::tri(log2p)->sm();++i)
						ijind[i+1][0] = npnt +(indx+1)*basis::tri(log2p)->sm() -(i+1);
				}
				else {
					for(i=0;i<basis::tri(log2p)->sm();++i)
						ijind[i+1][0] = npnt +indx*basis::tri(log2p)->sm() +i;
				}
				
				indx = tri(tind).seg(1);
				sgn = tri(tind).sgn(1);
				if (sgn > 0) {
					for(i=0;i<basis::tri(log2p)->sm();++i)
						ijind[basis::tri(log2p)->sm()-i][i+1] = npnt +indx*basis::tri(log2p)->sm() +i;
				}
				else {
					for(i=0;i<basis::tri(log2p)->sm();++i)
						ijind[basis::tri(log2p)->sm()-i][i+1] = npnt +(indx+1)*basis::tri(log2p)->sm() -(i+1);
				}
				
				indx = tri(tind).seg(2);
				sgn = tri(tind).sgn(2);
				if (sgn > 0) {
					for(i=0;i<basis::tri(log2p)->sm();++i)
						ijind[0][i+1] = npnt +(indx+1)*basis::tri(log2p)->sm() -(i+1);
				}
				else {
					for(i=0;i<basis::tri(log2p)->sm();++i)
						ijind[0][i+1] = npnt +indx*basis::tri(log2p)->sm() +i;
				}
				
				/* INTERIOR VERTICES */
				k = 0;
				for(i=1;i<basis::tri(log2p)->sm();++i) {
					for(j=1;j<basis::tri(log2p)->sm()-(i-1);++j) {
						ijind[i][j] = npnt +nseg*basis::tri(log2p)->sm() +tind*basis::tri(log2p)->im() +k;
						++k;
					}
				}
				
				/* OUTPUT CONNECTION LIST */        
				for(i=0;i<basis::tri(log2p)->sm()+1;++i) {
					for(j=0;j<basis::tri(log2p)->sm()-i;++j) {
						out << 3 << ' ' << ijind[i][j] << ' ' << ijind[i+1][j] << ' ' << ijind[i][j+1] << std::endl;
						out << 3 << ' ' << ijind[i+1][j] << ' ' << ijind[i+1][j+1] << ' ' << ijind[i][j+1] << std::endl;
					}
					out << 3 << ' ' << ijind[i][basis::tri(log2p)->sm()-i] << ' ' << ijind[i+1][basis::tri(log2p)->sm()-i] << ' ' << ijind[i][basis::tri(log2p)->sm()+1-i] << std::endl;
				}
			}
			
			
			out << "CELL_TYPES " << ntri*(basis::tri(log2p)->sm()+1)*(basis::tri(log2p)->sm()+1) << endl;
			
			for(i=0;i<ntri*(basis::tri(log2p)->sm()+1)*(basis::tri(log2p)->sm()+1);++i)
				out << 5 << endl;
			
			out << "POINT_DATA " << npnt+basis::tri(log2p)->sm()*nseg+basis::tri(log2p)->im()*ntri << endl;
			out << "SCALARS u float " << NV << endl;
			out << "LOOKUP_TABLE default" << endl;
			
			/* VERTEX MODES */
			for(i=0;i<npnt;++i) {
				for(n=0;n<NV;++n)
					out << ugbd(tlvl).v(i,n) << ' ';                    
				out << std::endl;
			}
			
			if (basis::tri(log2p)->p() > 1) {
				/* SIDE MODES */
				for(sind=0;sind<nseg;++sind) {
					ugtouht1d(sind,tlvl);
					for(n=0;n<NV;++n)
						basis::tri(log2p)->proj1d_leg(&uht(n)(0),&u(n)(0,0));
					
					for(i=1;i<basis::tri(log2p)->sm()+1;++i) {
						for(n=0;n<NV;++n)
							out << u(n)(0,i) << ' ';                    
						out << std::endl;
					}
				}
				
				/* INTERIOR MODES */
				if (basis::tri(log2p)->p() > 2) {
					for(tind = 0; tind < ntri; ++tind) {
						ugtouht(tind,tlvl);
						for(n=0;n<NV;++n)
							basis::tri(log2p)->proj_leg(&uht(n)(0),&u(n)(0,0),MXGP);
							
						for(i=1;i<basis::tri(log2p)->sm();++i) {
							for(j=1;j<basis::tri(log2p)->sm()-(i-1);++j) {
								for(n=0;n<NV;++n)
									out << u(n)(i,j) << ' ';                    
								out << std::endl;
							}
						}
					}
				}
			}
			
			
			out.close();
			break;
			
		case(tecplot):
			fnmapp = fname +".dat";
			out.open(fnmapp.c_str());
			if (!out) {
				*gbl->log << "couldn't open tecplot output file " << fnmapp;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}

			out << "ZONE F=FEPOINT, ET=TRIANGLE, N = " << npnt+basis::tri(log2p)->sm()*nseg+basis::tri(log2p)->im()*ntri << ", E = " << ntri*(basis::tri(log2p)->sm()+1)*(basis::tri(log2p)->sm()+1) << std::endl;

			/* VERTEX MODES */
			for(i=0;i<npnt;++i) {
				for(n=0;n<ND;++n)
					out << vrtxbd(tlvl)(i)(n) << ' ';
				for(n=0;n<NV;++n)
					out << ugbd(tlvl).v(i,n) << ' ';                    
				out << std::endl;
			}

			if (basis::tri(log2p)->p() > 1) {
				/* SIDE MODES */
				for(sind=0;sind<nseg;++sind) {
					if (seg(sind).info < 0) {
						v0 = seg(sind).pnt(0);
						v1 = seg(sind).pnt(1);
						for(n=0;n<ND;++n)
							basis::tri(log2p)->proj1d_leg(vrtxbd(tlvl)(v0)(n),vrtxbd(tlvl)(v1)(n),&crd(n)(0,0));
					}
					else {
						crdtocht1d(sind,tlvl);

						for(n=0;n<ND;++n)
							basis::tri(log2p)->proj1d_leg(&cht(n,0),&crd(n)(0,0));
					}
					ugtouht1d(sind,tlvl);
					for(n=0;n<NV;++n)
						basis::tri(log2p)->proj1d_leg(&uht(n)(0),&u(n)(0,0));

					for(i=1;i<basis::tri(log2p)->sm()+1;++i) {
						for(n=0;n<ND;++n)
							out << crd(n)(0,i) << ' ';
						for(n=0;n<NV;++n)
							out << u(n)(0,i) << ' ';                    
						out << std::endl;
					}
				}

				/* INTERIOR MODES */
				if (basis::tri(log2p)->p() > 2) {
					for(tind = 0; tind < ntri; ++tind) {
						ugtouht(tind,tlvl);
						for(n=0;n<NV;++n)
							basis::tri(log2p)->proj_leg(&uht(n)(0),&u(n)(0,0),MXGP);

						if (tri(tind).info < 0) {
							for(n=0;n<ND;++n)
								basis::tri(log2p)->proj_leg(vrtxbd(tlvl)(tri(tind).pnt(0))(n),vrtxbd(tlvl)(tri(tind).pnt(1))(n),vrtxbd(tlvl)(tri(tind).pnt(2))(n),&crd(n)(0,0),MXGP);
						}
						else {
							crdtocht(tind,tlvl);
							for(n=0;n<ND;++n)
								basis::tri(log2p)->proj_bdry_leg(&cht(n,0),&crd(n)(0,0),MXGP);
						}

						for(i=1;i<basis::tri(log2p)->sm();++i) {
							for(j=1;j<basis::tri(log2p)->sm()-(i-1);++j) {
								for(n=0;n<ND;++n)
									out << crd(n)(i,j) << ' ';
								for(n=0;n<NV;++n)
									out << u(n)(i,j) << ' ';                    
								out << std::endl;
							}
						}
					}
				}
			}

			/* OUTPUT CONNECTIVY INFO */
			out << std::endl << "#CONNECTION DATA#" << std::endl;

			for(tind=0;tind<ntri;++tind) {

				/* VERTICES */
				ijind[0][basis::tri(log2p)->sm()+1] = tri(tind).pnt(0);
				ijind[0][0] = tri(tind).pnt(1);
				ijind[basis::tri(log2p)->sm()+1][0] = tri(tind).pnt(2);

				/* SIDES */
				indx = tri(tind).seg(0);
				sgn = tri(tind).sgn(0);
				if (sgn < 0) {
					for(i=0;i<basis::tri(log2p)->sm();++i)
						ijind[i+1][0] = npnt +(indx+1)*basis::tri(log2p)->sm() -(i+1);
				}
				else {
					for(i=0;i<basis::tri(log2p)->sm();++i)
						ijind[i+1][0] = npnt +indx*basis::tri(log2p)->sm() +i;
				}

				indx = tri(tind).seg(1);
				sgn = tri(tind).sgn(1);
				if (sgn > 0) {
					for(i=0;i<basis::tri(log2p)->sm();++i)
						ijind[basis::tri(log2p)->sm()-i][i+1] = npnt +indx*basis::tri(log2p)->sm() +i;
				}
				else {
					for(i=0;i<basis::tri(log2p)->sm();++i)
						ijind[basis::tri(log2p)->sm()-i][i+1] = npnt +(indx+1)*basis::tri(log2p)->sm() -(i+1);
				}

				indx = tri(tind).seg(2);
				sgn = tri(tind).sgn(2);
				if (sgn > 0) {
					for(i=0;i<basis::tri(log2p)->sm();++i)
						ijind[0][i+1] = npnt +(indx+1)*basis::tri(log2p)->sm() -(i+1);
				}
				else {
					for(i=0;i<basis::tri(log2p)->sm();++i)
						ijind[0][i+1] = npnt +indx*basis::tri(log2p)->sm() +i;
				}

				/* INTERIOR VERTICES */
				k = 0;
				for(i=1;i<basis::tri(log2p)->sm();++i) {
					for(j=1;j<basis::tri(log2p)->sm()-(i-1);++j) {
						ijind[i][j] = npnt +nseg*basis::tri(log2p)->sm() +tind*basis::tri(log2p)->im() +k;
						++k;
					}
				}

				/* OUTPUT CONNECTION LIST */        
				for(i=0;i<basis::tri(log2p)->sm()+1;++i) {
					for(j=0;j<basis::tri(log2p)->sm()-i;++j) {
						out << ijind[i][j]+1 << ' ' << ijind[i+1][j]+1 << ' ' << ijind[i][j+1]+1 << std::endl;
						out << ijind[i+1][j]+1 << ' ' << ijind[i+1][j+1]+1 << ' ' << ijind[i][j+1]+1 << std::endl;
					}
					out << ijind[i][basis::tri(log2p)->sm()-i]+1 << ' ' << ijind[i+1][basis::tri(log2p)->sm()-i]+1 << ' ' << ijind[i][basis::tri(log2p)->sm()+1-i]+1 << std::endl;
				}
			}
			out.close();

			break; 

		case (datatank): {
#ifdef DATATANK
			fnmapp = fname +".dt";
			int dt_pnts = npnt+basis::tri(log2p)->sm()*nseg+basis::tri(log2p)->im()*ntri;
			int dt_tris = ntri*(basis::tri(log2p)->sm()+1)*(basis::tri(log2p)->sm()+1);

			DTMutableIntArray dt_tvrtx(3,dt_tris);
			DTMutableDoubleArray dt_vrtx(2,dt_pnts);
			DTMutableDoubleArray dt_vals(dt_pnts);


			/* VERTEX MODES */
			for(i=0;i<npnt;++i) {
				for(n=0;n<ND;++n)
					dt_vrtx(n,i) =  vrtxbd(tlvl)(i)(n);
				for(n=0;n<NV;++n)
					dt_vals(i) = ugbd(tlvl).v(i,n);                    
			}
			
			dt_pnts = npnt;
			if (basis::tri(log2p)->p() > 1) {
				/* SIDE MODES */
				for(sind=0;sind<nseg;++sind) {
					if (seg(sind).info < 0) {
						v0 = seg(sind).pnt(0);
						v1 = seg(sind).pnt(1);
						for(n=0;n<ND;++n)
							basis::tri(log2p)->proj1d_leg(vrtxbd(tlvl)(v0)(n),vrtxbd(tlvl)(v1)(n),&crd(n)(0,0));
					}
					else {
						crdtocht1d(sind,tlvl);
						
						for(n=0;n<ND;++n)
							basis::tri(log2p)->proj1d_leg(&cht(n,0),&crd(n)(0,0));
					}
					ugtouht1d(sind,tlvl);
					for(n=0;n<NV;++n)
						basis::tri(log2p)->proj1d_leg(&uht(n)(0),&u(n)(0,0));
					
					for(i=1;i<basis::tri(log2p)->sm()+1;++i) {
						for(n=0;n<ND;++n)
							dt_vrtx(n,dt_pnts) = crd(n)(0,i);
						for(n=0;n<NV;++n)
							dt_vals(dt_pnts) = u(n)(0,i);                    
						++dt_pnts;
					}
				}
				
				/* INTERIOR MODES */
				if (basis::tri(log2p)->p() > 2) {
					for(tind = 0; tind < ntri; ++tind) {
						ugtouht(tind,tlvl);
						for(n=0;n<NV;++n)
							basis::tri(log2p)->proj_leg(&uht(n)(0),&u(n)(0,0),MXGP);
						
						if (tri(tind).info < 0) {
							for(n=0;n<ND;++n)
								basis::tri(log2p)->proj_leg(vrtxbd(tlvl)(tri(tind).pnt(0))(n),vrtxbd(tlvl)(tri(tind).pnt(1))(n),vrtxbd(tlvl)(tri(tind).pnt(2))(n),&crd(n)(0,0),MXGP);
						}
						else {
							crdtocht(tind,tlvl);
							for(n=0;n<ND;++n)
								basis::tri(log2p)->proj_bdry_leg(&cht(n,0),&crd(n)(0,0),MXGP);
						}
						
						for(i=1;i<basis::tri(log2p)->sm();++i) {
							for(j=1;j<basis::tri(log2p)->sm()-(i-1);++j) {
								for(n=0;n<ND;++n)
									dt_vrtx(n,dt_pnts) = crd(n)(i,j);
								for(n=0;n<NV;++n)
									dt_vals(dt_pnts) = u(n)(i,j);                    
								++dt_pnts;
							}
						}
					}
				}
			}

			dt_tris = 0;
			for(tind=0;tind<ntri;++tind) {
				
				/* VERTICES */
				ijind[0][basis::tri(log2p)->sm()+1] = tri(tind).pnt(0);
				ijind[0][0] = tri(tind).pnt(1);
				ijind[basis::tri(log2p)->sm()+1][0] = tri(tind).pnt(2);
				
				/* SIDES */
				indx = tri(tind).seg(0);
				sgn = tri(tind).sgn(0);
				if (sgn < 0) {
					for(i=0;i<basis::tri(log2p)->sm();++i)
						ijind[i+1][0] = npnt +(indx+1)*basis::tri(log2p)->sm() -(i+1);
				}
				else {
					for(i=0;i<basis::tri(log2p)->sm();++i)
						ijind[i+1][0] = npnt +indx*basis::tri(log2p)->sm() +i;
				}
				
				indx = tri(tind).seg(1);
				sgn = tri(tind).sgn(1);
				if (sgn > 0) {
					for(i=0;i<basis::tri(log2p)->sm();++i)
						ijind[basis::tri(log2p)->sm()-i][i+1] = npnt +indx*basis::tri(log2p)->sm() +i;
				}
				else {
					for(i=0;i<basis::tri(log2p)->sm();++i)
						ijind[basis::tri(log2p)->sm()-i][i+1] = npnt +(indx+1)*basis::tri(log2p)->sm() -(i+1);
				}
				
				indx = tri(tind).seg(2);
				sgn = tri(tind).sgn(2);
				if (sgn > 0) {
					for(i=0;i<basis::tri(log2p)->sm();++i)
						ijind[0][i+1] = npnt +(indx+1)*basis::tri(log2p)->sm() -(i+1);
				}
				else {
					for(i=0;i<basis::tri(log2p)->sm();++i)
						ijind[0][i+1] = npnt +indx*basis::tri(log2p)->sm() +i;
				}
				
				/* INTERIOR VERTICES */
				k = 0;
				for(i=1;i<basis::tri(log2p)->sm();++i) {
					for(j=1;j<basis::tri(log2p)->sm()-(i-1);++j) {
						ijind[i][j] = npnt +nseg*basis::tri(log2p)->sm() +tind*basis::tri(log2p)->im() +k;
						++k;
					}
				}
				
				/* OUTPUT CONNECTION LIST */        
				for(i=0;i<basis::tri(log2p)->sm()+1;++i) {
					for(j=0;j<basis::tri(log2p)->sm()-i;++j) {
						dt_tvrtx(0,dt_tris) = ijind[i][j];
						dt_tvrtx(1,dt_tris) = ijind[i+1][j];
						dt_tvrtx(2,dt_tris) = ijind[i][j+1];
						++dt_tris;
						
						dt_tvrtx(0,dt_tris) = ijind[i+1][j];
						dt_tvrtx(1,dt_tris) = ijind[i+1][j+1];
						dt_tvrtx(2,dt_tris) = ijind[i][j+1];
						++dt_tris;
					}
					dt_tvrtx(0,dt_tris) = ijind[i][basis::tri(log2p)->sm()-i];
					dt_tvrtx(1,dt_tris) = ijind[i+1][basis::tri(log2p)->sm()-i];
					dt_tvrtx(2,dt_tris) = ijind[i][basis::tri(log2p)->sm()+1-i];
					++dt_tris;
				}
			}
			
			DTDataFile outputFile(fnmapp.c_str(),DTFile::NewReadWrite);
			outputFile.Save("Group","Seq_Var");
			outputFile.Save("mesh","SeqInfo_Var_1N");
			outputFile.Save("TriangularGrid2D","SeqInfo_Var_1T");
			
			outputFile.Save("var1","SeqInfo_Var_2N");
			outputFile.Save("TriangularMesh2D","SeqInfo_Var_2T");
			
			outputFile.Save("var2","SeqInfo_Var_3N");
			outputFile.Save("TriangularMesh2D","SeqInfo_Var_3T");
			
			outputFile.Save(3,"SeqInfo_Var_N");
			outputFile.Save("Group","SeqInfo_Var");

			// Structures for shared grids.
			DTTriangularGrid2D_SaveInfo Info_DTTriangularGrid2D;
			
			DTTriangularGrid2D dt_grid(dt_tvrtx,dt_vrtx);
			Write(outputFile,"Var_mesh",dt_grid,Info_DTTriangularGrid2D);
			DTTriangularMesh2D dt_mesh(dt_grid,dt_vals);
			Write(outputFile,"Var_var1",dt_mesh,Info_DTTriangularGrid2D);
			Write(outputFile,"Var_var2",dt_mesh,Info_DTTriangularGrid2D);
			Write(outputFile,"Var",DTDoubleArray()); // So that DataTank can see the variable.
#else
			*gbl->log << "Not supported on this platform" << std::endl;
#endif
			break;
		}

		case(adapt_diagnostic): {            
			fnmapp = fname+ "_length.dat";
			out.open(fnmapp.c_str());
			if (!out) {
				*gbl->log << "couldn't open tecplot output file " << fnmapp;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}

			out << "ZONE F=FEPOINT, ET=TRIANGLE, N = " << npnt << ", E = " << ntri << std::endl;
			for(i=0;i<npnt;++i) {
				for(n=0;n<ND;++n)
					out << pnts(i)(n) << ' ';
				out << lngth(i) << '\n';                    
			}

			/* OUTPUT CONNECTIVY INFO */
			out << "\n#CONNECTION DATA# \n";

			for(tind=0;tind<ntri;++tind)
				out << tri(tind).pnt(0)+1 << ' ' << tri(tind).pnt(1)+1 << ' ' << tri(tind).pnt(2)+1 << std::endl;

			out.close();
			
			fnmapp = fname +"_truncation.dat";
			out.open(fnmapp.c_str());
			if (!out) {
				*gbl->log << "couldn't open tecplot output file " << fnmapp << '\n';
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			
			/* OUTPUTS DISCONNECTED TRIS & ERROR FOR TRI */
			out << "ZONE F=FEPOINT, ET=TRIANGLE, N = " << 3*ntri << ", E = " << ntri << std::endl;
			for(tind=0;tind<ntri;++tind) {
				for(j=0;j<3;++j) {
					for(n=0;n<ND;++n)
						out << pnts(tri(tind).pnt(j))(n) << ' ';
					out << gbl->fltwk(tind) << '\n';
				}
			}

			/* OUTPUT CONNECTIVY INFO */
			out << "\n#CONNECTION DATA# \n";

			for(tind=0;tind<ntri;++tind)
				out << 3*tind+1 << ' ' << 3*tind+2 << ' ' << 3*tind+3 << std::endl;

			out.close();
			break;

		}

		default:
			*gbl->log << "can't output a tri_hp to that filetype" << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
			break;
	}
	 
	 /* BOUNDARY INFO */
	 for(i=0;i<nebd;++i) {
		 hp_ebdry(i)->output(filename,typ,tlvl);
	 }
	 
	 for(i=0;i<nvbd;++i) {
		 hp_vbdry(i)->output(filename,typ,tlvl);
	 }

	return;
}

void tri_hp::input(const std::string& filename) {
	int i,j;
	std::string fname,fnmapp;
	std::ostringstream nstr;
	ifstream fin;
	binifstream bin;
	
	fname = filename +"_" +gbl->idprefix;

	if (reload_type == tri_hp::binary) {
		fnmapp = fname +".bin";
		fin.open(fnmapp.c_str(),ios::in);
		if(fin.is_open()) {
			fin.close();
			input_map blank;
			tri_mesh::input(fnmapp,tri_mesh::binary,1,blank);
			setinfo();
			
			for(i=1;i<gbl->nadapt;++i) {
				nstr.str("");
				nstr << i << std::flush;
				fnmapp = filename +"_v" +nstr.str() +"_" +gbl->idprefix +".bin";
				bin.open(fnmapp.c_str());
				if (bin.error()) {
					*gbl->log << "#couldn't open input file " << fnmapp << std::endl;
					vrtxbd(i)(Range(0,npnt-1)) = pnts(Range(0,npnt-1));
				}
				else {
					bin.setFlag(binio::BigEndian,bin.readInt(1));
					bin.setFlag(binio::FloatIEEE,bin.readInt(1));
					for (j=0;j<npnt;++j) {
						vrtxbd(i)(j)(0) = bin.readFloat(binio::Double);
						vrtxbd(i)(j)(1) = bin.readFloat(binio::Double);
					}
				}
				bin.close();
			}
		}
		else {
			for(i=1;i<gbl->nadapt;++i)
				vrtxbd(i)(Range(0,npnt-1)) = pnts(Range(0,npnt-1));
		}

		for(i=0;i<gbl->nadapt;++i) {
			nstr.str("");
			nstr << i << std::flush;
			fnmapp = filename +"_d" +nstr.str();
			input(fnmapp,reload_type,i);
		}
	}
	else {
		fnmapp = fname +"_" +gbl->idprefix +".grd";
		fin.open(fnmapp.c_str(),ios::in);
		if(fin.is_open()) {
			fin.close();
			input_map blank;
			tri_mesh::input(fnmapp,tri_mesh::grid,1,blank);
			setinfo();
			for(i=1;i<gbl->nadapt;++i) {
				nstr.str("");
				nstr << i << std::flush;
				fnmapp = filename +"_v" +nstr.str() +"_" +gbl->idprefix +".txt";
				fin.open(fnmapp.c_str());
				if (!fin.is_open()) {
					*gbl->log << "couldn't open input file " << fnmapp << std::endl;
					sim::abort(__LINE__,__FILE__,gbl->log);
				}
				fin.ignore(80,'\n');  // SKIP NUMBER OF VERTICES
				for (j=0;j<npnt;++j) {
					fin.ignore(80,':');
					fin >> vrtxbd(i)(j)(0) >> vrtxbd(i)(j)(1);
				}
				fin.close();
			}
		}
		else {
			for(i=1;i<gbl->nadapt;++i)
				vrtxbd(i)(Range(0,npnt-1)) = pnts(Range(0,npnt-1));
		}


		for(i=0;i<gbl->nadapt;++i) {
			nstr.str("");
			nstr << i << std::flush;
			fnmapp = filename +"_d" +nstr.str();
			input(fnmapp,reload_type,i);
		}
	}

	return;

}

void tri_hp::input(const std::string& filename, filetype typ, int tlvl) {
	int i,j,k,m,n,pin,pmin,indx;
	int bnum;
	std::string fnapp;
	char buffer[80];
	char trans[] = "T";
	ifstream in;
	FLT fltskip;
	
	std::string fname;
	fname = filename +"_" +gbl->idprefix;

	switch(typ) {

		case (text):
			fnapp = fname +".txt";
			in.open(fnapp.c_str());
			if (!in) {
				*gbl->log << "couldn't open text input file " << fnapp << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}

			/* HEADER INFORMATION */
			/* INPUT # OF SIDE MODES (ONLY THING THAT CAN BE DIFFERENT) THEN SKIP THE REST */
			in.ignore(80,'=');
			in >> pin;
			pmin = MIN(p0,pin);

			do {
				in.ignore(80,'\n');
				buffer[0] = in.get();
			} while (buffer[0] != 'E');  // END OF HEADER
			in.ignore(80,'\n');

			for(i=0;i<npnt;++i) {
				for(n=0;n<NV;++n)
					in >> ugbd(tlvl).v(i,n);
			}

			for(i=0;i<nseg;++i) {
				for(m=0;m<(pmin-1);++m) {
					for(n=0;n<NV;++n)
						in >> ugbd(tlvl).s(i,m,n);
				}

				for(m=0;m<(pin-p0);++m) {
					for(n=0;n<NV;++n)
						in >> fltskip;
				}              
			}

			for(i=0;i<ntri;++i) {
				indx = 0;
				for(m=1;m<pmin-1;++m) {
					for(k=0;k<pmin-1-m;++k) {
						for(n=0;n<NV;++n) 
							in >> ugbd(tlvl).i(i,indx,n);
						++indx;
					}
					indx += p0 -pmin;

					for(k=0;k<pin-p0;++k) {
						for(n=0;n<NV;++n) 
							in >> fltskip;
					}
				}

				for(m=pmin-1;m<pin-1;++m) {
					for(k=0;k<pin-1-m;++k) {
						for(n=0;n<NV;++n) 
							in >> fltskip; 
					}
				}              
			}
			in.ignore(80,'\n');


			/* BOUNDARY INFO */
			for(i=0;i<nebd;++i)
				hp_ebdry(i)->input(in,typ,tlvl);

			for(i=0;i<nvbd;++i)
				hp_vbdry(i)->input(in,typ,tlvl);

			in.close();
			break;

		case (binary): {
			fnapp = fname +".bin";
			in.open(fnapp.c_str());
			if (!in) {
				*gbl->log << "couldn't open binary input file " << fnapp << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			biniwstream bin(&in);

			/* HEADER INFORMATION */
			bin.setFlag(binio::BigEndian,bin.readInt(1));
			bin.setFlag(binio::FloatIEEE,bin.readInt(1));

			/* INPUT # OF SIDE MODES (ONLY THING THAT CAN BE DIFFERENT) THEN SKIP THE REST */
			pin = bin.readInt(sizeof(int));
			pmin = MIN(p0,pin);

			if (bin.readInt(sizeof(int))  != npnt) {
				*gbl->log << "mismatched pnt counts?" << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			if (bin.readInt(sizeof(int))  != nseg) {
				*gbl->log << "mismatched seg counts" << std::endl;;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			if (bin.readInt(sizeof(int))  != ntri) {
				*gbl->log << "mismatched tri counts?" << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}

			for(i=0;i<npnt;++i) {
				for(n=0;n<NV;++n)
					ugbd(tlvl).v(i,n) = bin.readFloat(binio::Double);
			}

			for(i=0;i<nseg;++i) {
				for(m=0;m<(pmin-1);++m) {
					for(n=0;n<NV;++n)
						ugbd(tlvl).s(i,m,n) = bin.readFloat(binio::Double);
				}

				for(m=0;m<(pin-p0);++m) {
					for(n=0;n<NV;++n)
						bin.readFloat(binio::Double);
				}              
			}

			for(i=0;i<ntri;++i) {
				indx = 0;
				for(m=1;m<pmin-1;++m) {
					for(k=0;k<pmin-1-m;++k) {
						for(n=0;n<NV;++n) 
							ugbd(tlvl).i(i,indx,n) = bin.readFloat(binio::Double);
						++indx;
					}
					indx += p0 -pmin;

					for(k=0;k<pin-p0;++k) {
						for(n=0;n<NV;++n) 
							bin.readFloat(binio::Double);
					}
				}

				for(m=pmin-1;m<pin-1;++m) {
					for(k=0;k<pin-1-m;++k) {
						for(n=0;n<NV;++n) 
							bin.readFloat(binio::Double); 
					}
				}              
			}

			/* BOUNDARY INFO */
			for(i=0;i<nebd;++i)
				hp_ebdry(i)->input(in,typ,tlvl);

			for(i=0;i<nvbd;++i)
				hp_vbdry(i)->input(in,typ,tlvl);

			in.close();
			break;
		}

		case(tecplot):
			/* CAN ONLY DO THIS IF HAVE MESH FILE */
			fnapp = fname +".dat";
			in.open(fnapp.c_str());
			if (!in) {
				*gbl->log << "couldn't open tecplot input file " << fnapp << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			in.ignore(80,'\n');

			for(i=0;i<npnt;++i) {
				in >> fltskip >> fltskip;
				for(n=0;n<NV;++n)
					in >> ugbd(tlvl).v(i,n);
			}

			if (basis::tri(log2p)->sm() < 1) return;

			FLT temp,matrix[MXTM][MXTM];
			int sind,info,ipiv[2*MXTM];

			/* REVERSE OUTPUTING PROCESS */
			for(m=0;m<basis::tri(log2p)->sm();++m)
				for(n=0;n<basis::tri(log2p)->sm();++n)
					matrix[n][m] = basis::tri(log2p)->lgrnge1d(m+2,n+1);

			GETRF(basis::tri(log2p)->sm(),basis::tri(log2p)->sm(),matrix[0],MXTM,ipiv,info);
			if (info != 0) {
				*gbl->log << "DGETRF FAILED FOR INPUTING TECPLOT SIDES" << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}

			for(sind=0;sind<nseg;++sind) {
				if (seg(sind).info < 0) {
					ugtouht1d(sind,tlvl);
					for(n=0;n<NV;++n)
						basis::tri(log2p)->proj1d_leg(&uht(n)(0),&u(n)(0,0));

					for(m=0;m<basis::tri(log2p)->sm();++m) {
						in >> fltskip >> fltskip;
						for(n=0;n<NV;++n) {
							in >> temp;
							u(n)(0,m+1) -= temp;
						}
					}
					for(n=0;n<NV;++n) {
						GETRS(trans,basis::tri(log2p)->sm(),1,matrix[0],MXTM,ipiv,&u(n)(0,1),MXTM,info);
						for(m=0;m<basis::tri(log2p)->sm();++m)
							ugbd(tlvl).s(sind,m,n) = -u(n)(0,1+m);
					}
				}
				else {
					crdtocht1d(sind,tlvl);
					for(n=0;n<ND;++n)
						basis::tri(log2p)->proj1d_leg(&cht(n,0),&crd(n)(0,0));

					ugtouht1d(sind,tlvl);
					for(n=0;n<NV;++n)
						basis::tri(log2p)->proj1d_leg(&uht(n)(0),&u(n)(0,0));

					for(m=0;m<basis::tri(log2p)->sm();++m) {
						for(n=0;n<ND;++n) {
							in >> temp;
							crd(n)(0,m+1) -= temp;
						}

						for(n=0;n<NV;++n) {
							in >> temp;
							u(n)(0,m+1) -= temp;
						}
					}
					for(n=0;n<NV;++n) {
						GETRS(trans,basis::tri(log2p)->sm(),1,matrix[0],MXTM,ipiv,&u(n)(0,1),MXTM,info);
						for(m=0;m<basis::tri(log2p)->sm();++m)
							ugbd(tlvl).s(sind,m,n) = -u(n)(0,1+m);
					}

					bnum = getbdrynum(seg(sind).tri(1));
					indx = getbdryseg(seg(sind).tri(1));
					for(n=0;n<ND;++n) {
						GETRS(trans,basis::tri(log2p)->sm(),1,matrix[0],MXTM,ipiv,&crd(n)(0,1),MXTM,info);
						for(m=0;m<basis::tri(log2p)->sm();++m)
							hp_ebdry(bnum)->crdsbd(tlvl,indx,m,n) = -crd(n)(0,1+m);
					}
				}
			}

			if (basis::tri(log2p)->im() < 1) return;

			int tind;

			ugbd(tlvl).i(Range(0,ntri-1),Range::all(),Range::all()) = 0.0;

			/* REVERSE OUTPUTING PROCESS */
			for(m=0;m<basis::tri(log2p)->im();++m) {
				n = 0;
				for(i=1;i<basis::tri(log2p)->sm();++i) {
					for(j=1;j<basis::tri(log2p)->sm()-(i-1);++j) {
						matrix[n++][m] = basis::tri(log2p)->lgrnge(m+basis::tri(log2p)->bm(),i,j);
					}
				}
			}

			GETRF(basis::tri(log2p)->im(),basis::tri(log2p)->im(),matrix[0],MXTM,ipiv,info);
			if (info != 0) {
				*gbl->log << "DGETRF FAILED FOR INPUTING TECPLOT SIDES" << std::endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}

			for(tind=0;tind<ntri;++tind) {
				ugtouht(tind,tlvl);
				for(n=0;n<NV;++n)
					basis::tri(log2p)->proj_leg(&uht(n)(0),&u(n)(0,0),MXGP);

				m = 0;
				for(i=1;i<basis::tri(log2p)->sm();++i) {
					for(j=1;j<basis::tri(log2p)->sm()-(i-1);++j) {
						for(n=0;n<NV;++n)
							uht(n)(m) = u(n)(i,j);
						in >> fltskip >> fltskip;
						for(n=0;n<NV;++n) {
							in >> temp;    
							uht(n)(m) -= temp;
						}
						++m;
					}
				}
				for(n=0;n<NV;++n) {
					GETRS(trans,basis::tri(log2p)->im(),1,matrix[0],MXTM,ipiv,&uht(n)(0),MXTM,info);
					for(m=0;m<basis::tri(log2p)->im();++m)
						ugbd(tlvl).i(tind,m,n) = -uht(n)(m);
				}
			}                

			break;

		default:
			*gbl->log << "can't input a tri_hp from that filetype" << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
			break;
	}

	return;
}
