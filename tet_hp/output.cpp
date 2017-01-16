/*
 *  output.cpp
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 15 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp.h"
#include "hp_boundary.h"
#include <cstring>
#include <assert.h>
#include <stdlib.h>
#include <myblas.h>
#include <libbinio/binwrap.h>
#include <libbinio/binfile.h>

// #define DATATANK
#ifdef DATATANK
#include <DTSource.h>
#endif


void tet_hp::output(const std::string& fname, block::output_purpose why) {
	int i,j;
	std::string fnmapp, namewdot;
	std::ostringstream nstr;
	ofstream out;
	out.setf(std::ios::scientific, std::ios::floatfield);
	out.precision(10);
	
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
			if (output_type(1) == tet_hp::binary) {
				tet_mesh::output(fname,tet_mesh::binary);
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
						bout.writeFloat(vrtxbd(i)(j)(2),binio::Double);
					}
					bout.close();
				}
			}
			else {
				tet_mesh::output(fname,tet_mesh::grid);
				tet_mesh::output(fname,tet_mesh::vlength);
				for(i=1;i<gbl->nadapt;++i) {
					nstr.str("");
					nstr << i << std::flush;
					fnmapp = namewdot +nstr.str() +"_" +gbl->idprefix +".txt";
					out.open(fnmapp.c_str());
					out << npnt << std::endl;
					for (j=0;j<npnt;++j) 
						out << j << ':' << vrtxbd(i)(j)(0) << ' ' << vrtxbd(i)(j)(1) << ' ' << vrtxbd(i)(j)(2) << std::endl;
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

void tet_hp::output(const std::string& filename, filetype typ, int tlvl) {
	ofstream out;
	std::string fname, fnmapp;
	int i,j,k,m,n,v0,v1,eind,find,tind,indx,sgn,rot;
	int ijind[basis::tet(log2p).em+2][basis::tet(log2p).em+2][basis::tet(log2p).em+2];
	int stridey = MXGP;
	int stridex = MXGP*MXGP;
	TinyVector<double,3> pt;
	
	out.setf(std::ios::scientific, std::ios::floatfield);
	out.precision(10);
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
			out << "npnt = " << npnt << ", nseg = " << nseg << ", ntri = " << ntri << ", ntet = " << ntet << std::endl;
			out << "END OF HEADER" << std::endl;

			for(i=0;i<npnt;++i) {
				for(n=0;n<NV;++n)
					out << ugbd(tlvl).v(i,n) << '\t';
				out << std::endl;
			}
			
			for(i=0;i<nseg;++i) {
				for(m=0;m<em0;++m) {
					for(n=0;n<NV;++n)
						out << ugbd(tlvl).e(i,m,n) << '\t';
					out << std::endl;
				}
			}
			
			for(i=0;i<ntri;++i) {
				for(m=0;m<fm0;++m) {
					for(n=0;n<NV;++n)
						out << ugbd(tlvl).f(i,m,n) << '\t';
					out << std::endl;
				}
			}
			
			for(i=0;i<ntet;++i) {
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
				*gbl->log << "couldn't open binary output file " << fnmapp;
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
			bout.writeInt(ntet,sizeof(int));
			
			for(i=0;i<npnt;++i) {
				for(n=0;n<NV;++n)
					bout.writeFloat(ugbd(tlvl).v(i,n),binio::Double);
			}
			
			for(i=0;i<nseg;++i) {
				for(m=0;m<em0;++m) {
					for(n=0;n<NV;++n)
						bout.writeFloat(ugbd(tlvl).e(i,m,n),binio::Double);
				}
			}
			
			for(i=0;i<ntri;++i) {
				for(m=0;m<fm0;++m) {
					for(n=0;n<NV;++n)
						bout.writeFloat(ugbd(tlvl).f(i,m,n),binio::Double);
				}
			}
			
			for(i=0;i<ntet;++i) {
				for(m=0;m<im0;++m) {
					for(n=0;n<NV;++n)
						bout.writeFloat(ugbd(tlvl).i(i,m,n),binio::Double);
				}
			}
			
			out.close();
			break;
		}

		case (vtu): {
			fnmapp = fname +".vtu";
			out.open(fnmapp.c_str());
			if (!out) {
				*gbl->log<< "couldn't open vtu output file " << fnmapp;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			if (basis::tet(log2p).p > 2) {
				*gbl->log << "vtu output routine does not work for p > 2 " << endl;
				exit(2);
			}
			
			int numpnts,numtets;
			numpnts = npnt+basis::tet(log2p).em*nseg+basis::tet(log2p).fm*ntri+basis::tet(log2p).im*ntet;
			numtets = ntet*(basis::tet(log2p).em+1)*(basis::tet(log2p).em+1)*(basis::tet(log2p).em+1);

			out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
			out << "	<UnstructuredGrid>" << endl;
			out << "		<Piece NumberOfPoints=\"" << numpnts << "\" NumberOfCells=\"" << numtets << "\">" << endl;
			if(NV == 1) {
				out << "			<PointData Scalars=\"Temperature\">" << endl;
				out << "				<DataArray type=\"Float32\" Name=\"Temperature\" format=\"ascii\">" << endl;
				
			} 
			else {
				out << "			<PointData Scalars=\"Pressure\" Vectors=\"Velocity\">" << endl;
				out << "				<DataArray type=\"Float32\" Name=\"Pressure\" format=\"ascii\">" << endl;				
			}

			/* VERTEX MODES */
			for(i=0;i<npnt;++i) 
				out << ugbd(tlvl).v(i,0)<< std::endl;
			
			if (basis::tet(log2p).p > 1) {
				/* EDGE MODES */
				for(eind=0;eind<nseg;++eind) {
					ugtouht1d(eind,tlvl);
					basis::tet(log2p).proj1d_leg(&uht(0)(0),&u1d(0)(0));
					
					for(i=1;i<basis::tet(log2p).em+1;++i) 
						out << u1d(0)(i)<< std::endl;               
						         					
				}
			}
			
			out << "				</DataArray>" << endl;
			
			if (NV > 1) {			
				out << "				<DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
				
				/* VERTEX MODES */
				for(i=0;i<npnt;++i) {
					for(n=1;n<NV-1;++n) 
						out << ugbd(tlvl).v(i,n)<< ' ';
					out << std::endl;
				}
				
				if (basis::tet(log2p).p > 1) {
					/* EDGE MODES */
					for(eind=0;eind<nseg;++eind) {
						ugtouht1d(eind,tlvl);
						for(n=1;n<NV-1;++n)
							basis::tet(log2p).proj1d_leg(&uht(n)(0),&u1d(n)(0));
						
						for(i=1;i<basis::tet(log2p).em+1;++i) {
							for(n=1;n<NV-1;++n){
								out << u1d(n)(i)<< ' ';               
							}            
							out << std::endl; 
						}
					}
				}
				
				out << "				</DataArray>" << endl;
			}
			if (NV > 4) {
				out << "				<DataArray type=\"Float32\" Name=\"Temperature\" format=\"ascii\">" << endl;
				
				/* VERTEX MODES */
				for(i=0;i<npnt;++i) 
					out << ugbd(tlvl).v(i,NV-1)<< std::endl;
				
				if (basis::tet(log2p).p > 1) {
					/* EDGE MODES */
					for(eind=0;eind<nseg;++eind) {
						ugtouht1d(eind,tlvl);
						basis::tet(log2p).proj1d_leg(&uht(NV-1)(0),&u1d(NV-1)(0));
						
						for(i=1;i<basis::tet(log2p).em+1;++i) 
							out << u1d(NV-1)(i)<< std::endl;               
						
					}
				}
				
				out << "				</DataArray>" << endl;
			}
			
			out << "			</PointData>" << endl;
			
			out << "			<CellData>" << endl;
			out << "			</CellData>" << endl;
			
			out << "			<Points>" << endl;
			out << "				<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
			for(i=0;i<npnt;++i) {
				for(n=0;n<ND;++n) 
					out << vrtxbd(tlvl)(i)(n) << ' ';
				          
				out << std::endl;
			}
			
			if (basis::tet(log2p).p > 1) {
				/* EDGE MODES */
				for(eind=0;eind<nseg;++eind) {
					if (seg(eind).info < 0) {
						v0 = seg(eind).pnt(0);
						v1 = seg(eind).pnt(1);
						for(n=0;n<ND;++n)
							basis::tet(log2p).proj1d_leg(vrtxbd(tlvl)(v0)(n),vrtxbd(tlvl)(v1)(n),&crd1d(n)(0));
					}
					else {
						crdtocht1d(eind,tlvl);
						
						for(n=0;n<ND;++n)
							basis::tet(log2p).proj1d_leg(&cht(n)(0),&crd1d(n)(0));
					}
					
					for(i=1;i<basis::tet(log2p).em+1;++i) {
						for(n=0;n<ND;++n) 
							out << crd1d(n)(i) << ' ';
						            
						out << std::endl; 
					}
				}
			}
			
			out << "				</DataArray>" << endl;
			out << "			</Points>" << endl;
			
			out << "			<Cells>" << endl;
			out << "				<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << endl;
			/* OUTPUT CONNECTIVY INFO */
			for(tind=0;tind<ntet;++tind) {
				
				/* VERTICES */
				ijind[0][0][basis::tet(log2p).em+1] = tet(tind).pnt(0);
				ijind[0][basis::tet(log2p).em+1][0] = tet(tind).pnt(1);
				ijind[0][0][0] = tet(tind).pnt(2);
				ijind[basis::tet(log2p).em+1][0][0] = tet(tind).pnt(3);
				
				/* EDGES */
				indx = tet(tind).seg(0);
				sgn = tet(tind).sgn(0);
				if (sgn < 0) {  
					for(i=0;i<basis::tet(log2p).em;++i)
						ijind[i+1][0][0] = npnt +(indx+1)*basis::tet(log2p).em -(i+1);
				}
				else {
					for(i=0;i<basis::tet(log2p).em;++i)
						ijind[i+1][0][0] = npnt +indx*basis::tet(log2p).em +i;
				}
				
				indx = tet(tind).seg(1);
				sgn = tet(tind).sgn(1);
				if (sgn > 0) { 
					for(i=0;i<basis::tet(log2p).em;++i)
						ijind[basis::tet(log2p).em-i][i+1][0] = npnt +indx*basis::tet(log2p).em +i;
				}
				else {
					for(i=0;i<basis::tet(log2p).em;++i)
						ijind[basis::tet(log2p).em-i][i+1][0] = npnt +(indx+1)*basis::tet(log2p).em -(i+1);
				}
				
				indx = tet(tind).seg(2);
				sgn = tet(tind).sgn(2);
				if (sgn < 0) { 
					for(i=0;i<basis::tet(log2p).em;++i)
						ijind[0][i+1][0] = npnt +(indx+1)*basis::tet(log2p).em -(i+1);
				}
				else {
					for(i=0;i<basis::tet(log2p).em;++i)
						ijind[0][i+1][0] = npnt +indx*basis::tet(log2p).em +i;
				}
				
				indx = tet(tind).seg(3);
				sgn = tet(tind).sgn(3);
				if (sgn > 0) {  
					for(i=0;i<basis::tet(log2p).em;++i)
						ijind[0][i+1][basis::tet(log2p).em-i] = npnt +(indx+1)*basis::tet(log2p).em -(i+1);
				}
				else {
					for(i=0;i<basis::tet(log2p).em;++i)
						ijind[0][i+1][basis::tet(log2p).em-i] = npnt +indx*basis::tet(log2p).em +i;
				}
				
				indx = tet(tind).seg(4);
				sgn = tet(tind).sgn(4);
				if (sgn > 0) { 
					for(i=0;i<basis::tet(log2p).em;++i)
						ijind[0][0][i+1] = npnt +indx*basis::tet(log2p).em +i;
				}
				else {
					for(i=0;i<basis::tet(log2p).em;++i)
						ijind[0][0][i+1] = npnt +(indx+1)*basis::tet(log2p).em -(i+1);
				}
				
				indx = tet(tind).seg(5);
				sgn = tet(tind).sgn(5);
				if (sgn > 0) { 
					for(i=0;i<basis::tet(log2p).em;++i)
						ijind[i+1][0][basis::tet(log2p).em-i] = npnt +(indx+1)*basis::tet(log2p).em -(i+1);
				}
				else {
					for(i=0;i<basis::tet(log2p).em;++i)
						ijind[i+1][0][basis::tet(log2p).em-i] = npnt +indx*basis::tet(log2p).em +i;
				}
				
				/* FACES */
				k = 0;
				indx = tet(tind).tri(0);
				rot = -tet(tind).rot(0);
				if (rot == -1){
					for(i=1;i<basis::tet(log2p).em;++i) { 
						for(j=1;j<basis::tet(log2p).em-(i-1);++j) {
							ijind[basis::tet(log2p).em-i-j+1][j][0] = npnt +nseg*basis::tet(log2p).em +indx*basis::tet(log2p).fm +k;
							++k;
						}
					}
				}
				if (rot == 1) {
					for(i=1;i<basis::tet(log2p).em;++i) { 
						for(j=1;j<basis::tet(log2p).em-(i-1);++j) {
							ijind[i][j][0] = npnt +nseg*basis::tet(log2p).em +indx*basis::tet(log2p).fm +k;
							++k;
						}
					}
				}
				
				k = 0;
				indx = tet(tind).tri(1);
				rot = -tet(tind).rot(1);
				if(rot == 1){
					for(i=1;i<basis::tet(log2p).em;++i) {
						for(j=1;j<basis::tet(log2p).em-(i-1);++j) {
							ijind[basis::tet(log2p).em-i-j+1][0][j] = npnt +nseg*basis::tet(log2p).em +indx*basis::tet(log2p).fm +k;
							++k;
						}
					}
				}
				if (rot == -1) {
					for(i=1;i<basis::tet(log2p).em;++i) {
						for(j=1;j<basis::tet(log2p).em-(i-1);++j) {
							ijind[i][0][j] = npnt +nseg*basis::tet(log2p).em +indx*basis::tet(log2p).fm +k;
							++k;
						}
					}
				}
				
				k = 0;
				indx = tet(tind).tri(2);
				rot = -tet(tind).rot(2);
				if (rot == 1){
					for(i=1;i<basis::tet(log2p).em;++i) {
						for(j=1;j<basis::tet(log2p).em-(i-1);++j) {
							ijind[i][basis::tet(log2p).em-i-j+1][j] = npnt +nseg*basis::tet(log2p).em +indx*basis::tet(log2p).fm +k;
							++k;
						}
					}
				}
				if (rot == -1) {
					for(i=1;i<basis::tet(log2p).em;++i) {
						for(j=1;j<basis::tet(log2p).em-(i-1);++j) {
							ijind[basis::tet(log2p).em-i-j+1][i][j] = npnt +nseg*basis::tet(log2p).em +indx*basis::tet(log2p).fm +k;
							++k;
						}
					}			
				}
				
				k = 0;
				indx = tet(tind).tri(3);
				rot = -tet(tind).rot(3);
				if (rot == -1){
					for(i=1;i<basis::tet(log2p).em;++i) {
						for(j=1;j<basis::tet(log2p).em-(i-1);++j) {
							ijind[0][basis::tet(log2p).em-i-j+1][j] = npnt +nseg*basis::tet(log2p).em +indx*basis::tet(log2p).fm +k;
							++k;
						}
					}
				}
				if (rot == 1) {
					for(i=1;i<basis::tet(log2p).em;++i) {
						for(j=1;j<basis::tet(log2p).em-(i-1);++j) {
							ijind[0][i][j] = npnt +nseg*basis::tet(log2p).em +indx*basis::tet(log2p).fm +k;
							++k;
						}
					}
				}
				
				/* INTERIOR */
				m = 0;
				for(i=1;i<basis::tet(log2p).em;++i) {
					for(j=1;j<basis::tet(log2p).em-(i-1);++j) {
						for(k=1;k<basis::tet(log2p).em-i-j+1;++k){
							ijind[i][j][k] = npnt +nseg*basis::tet(log2p).em +ntri*basis::tet(log2p).fm +tind*basis::tet(log2p).im +m;
							++m;
						}
					}
				}
				
				
				/* OUTPUT CONNECTION LIST */     
				int em = basis::tet(log2p).em;
				for(i=0;i<em;++i) {
					for(j=0;j<em-i;++j) {
						for(k=0;k<em-i-j-1;++k) { 
							out << ijind[i][j][k] << ' ' << ijind[i][j][k+1] << ' ' << ijind[i+1][j][k] <<  ' ' << ijind[i][j+1][k] << std::endl;
							out << ijind[i][j][k+1] << ' ' << ijind[i+1][j][k+1] << ' ' << ijind[i+1][j][k] <<  ' ' << ijind[i][j+1][k] << std::endl;
							out << ijind[i][j][k+1] << ' ' << ijind[i+1][j][k+1] << ' ' << ijind[i][j+1][k] <<  ' ' << ijind[i][j+1][k+1] << std::endl;
							out << ijind[i+1][j][k+1] << ' ' << ijind[i][j+1][k] << ' ' << ijind[i+1][j+1][k] <<  ' ' << ijind[i+1][j][k] << std::endl;
							out << ijind[i+1][j][k+1] << ' ' << ijind[i][j+1][k] << ' ' << ijind[i][j+1][k+1] <<  ' ' << ijind[i+1][j+1][k] << std::endl;
							out << ijind[i][j+1][k+1] << ' ' << ijind[i+1][j][k+1] << ' ' << ijind[i+1][j+1][k] <<  ' ' << ijind[i+1][j+1][k+1] << std::endl;
							
						}			
						// five tet new way
						out << ijind[i][j][em-j-1-i] << ' ' << ijind[i][j][em-j-i] << ' ' << ijind[i+1][j][em-j-1-i] <<  ' ' << ijind[i][j+1][em-j-1-i] << std::endl;
						out << ijind[i][j][em-j-i] << ' ' << ijind[i+1][j][em-j-i] << ' ' << ijind[i+1][j][em-j-1-i] <<  ' ' << ijind[i][j+1][em-j-1-i] << std::endl;
						out << ijind[i][j][em-j-i] << ' ' << ijind[i+1][j][em-j-i] << ' ' << ijind[i][j+1][em-j-1-i] <<  ' ' << ijind[i][j+1][em-j-i] << std::endl;
						out << ijind[i+1][j][em-j-i] << ' ' << ijind[i][j+1][em-j-1-i] << ' ' << ijind[i+1][j+1][em-j-1-i] <<  ' ' << ijind[i+1][j][em-j-1-i] << std::endl;
						out << ijind[i+1][j][em-j-i] << ' ' << ijind[i][j+1][em-j-1-i] << ' ' << ijind[i][j+1][em-j-i] <<  ' ' << ijind[i+1][j+1][em-j-1-i] << std::endl;
						
						// single tet
						out << ijind[i][j][em-j-i] << ' ' << ijind[i+1][j][em-j-i] << ' ' << ijind[i][j+1][em-j-i] <<  ' ' << ijind[i][j][em-j+1-i] << std::endl;
						
					}
					//single tet
					out << ijind[i][em-i][0] << ' ' << ijind[i+1][em-i][0] << ' ' << ijind[i][em+1-i][0] <<  ' ' << ijind[i][em-i][1] << std::endl;
				}
				//single tet
				out << ijind[em][0][0] << ' ' << ijind[em+1][0][0] << ' ' << ijind[em][1][0] <<  ' ' << ijind[em][0][1] << std::endl;			
			}
			
			
			out << "				</DataArray>" << endl;
			out << "				<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << endl;
			out << "					";
			/* offsets */
			for(int i = 0; i < numtets; ++i)
				out << (i+1)*4 << ' ';
			out << endl;
			
			out << "				</DataArray>" << endl;
			out << "				<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
			out << "					";
			/* vtk file type 10 for tetrahedrals 24 for quadratic tet */
			for(int i = 0; i < numtets; ++i)
				out << 10 << ' '; 
			out << endl;
			
			out << "				</DataArray>" << endl;
			out << "			</Cells>" << endl;
			out << "		</Piece>" << endl;
			out << "	</UnstructuredGrid>" << endl;
			out << "</VTKFile>" << endl;
			
			out.close();
			if(gbl->idnum==0){
				std::ostringstream nstr;
				nstr.str("");
				int tstep = gbl->tstep;
				if (tstep == -1) tstep = 0; 
				
				nstr << tstep << std::flush;
				fnmapp = "data" +nstr.str()+".pvtu";				
				out.open(fnmapp.c_str());
				
				fnmapp = "data" +nstr.str();
				
				out << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << endl;
				out << "	<PUnstructuredGrid GhostLevel=\"0\">" << endl;
				if(NV == 1) {
					out << "		<PPointData Scalars=\"Temperature\">" << endl;
					out << "			<PDataArray type=\"Float32\" Name=\"Temperature\"/>" << endl;
				}
				else{
					out << "		<PPointData Scalars=\"Pressure\" Vectors=\"Velocity\">" << endl;
					out << "			<PDataArray type=\"Float32\" Name=\"Pressure\"/>" << endl;
				}
				if(NV > 1) out << "			<PDataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\"/>" << endl;
				if(NV > 4) out << "			<PDataArray type=\"Float32\" Name=\"Temperature\"/>" << endl;
				out << "		</PPointData>" << endl;
				out << "		<PCellData>" << endl;
				out << "		</PCellData>" << endl;
				out << "		<PPoints>" << endl;
				out << "			<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>" << endl;
				out << "		</PPoints>" << endl;
				
				for(int i=0;i<sim::blks.nblock; ++i)
					out << "		<Piece Source=\"" << fnmapp << "_b" << i << ".vtu\"/>" << endl;

				out << "	</PUnstructuredGrid>" << endl; 
				out << "</VTKFile>" << endl;
				
				out.close();
			}
			
			break;
		}
			
		case(tecplot): {
			fnmapp = fname +".dat";
			out.open(fnmapp.c_str());
			if (!out) {
				*gbl->log<< "couldn't open tecplot output file " << fnmapp;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			out << "ZONE F=FEPOINT, ET=TETRAHEDRON, N = " << npnt+basis::tet(log2p).em*nseg+basis::tet(log2p).fm*ntri+basis::tet(log2p).im*ntet << ", E = "  << ntet*(basis::tet(log2p).em+1)*(basis::tet(log2p).em+1)*(basis::tet(log2p).em+1) << std::endl;
			
			/* VERTEX MODES */
			for(i=0;i<npnt;++i) {
				for(n=0;n<ND;++n) {
					out << vrtxbd(tlvl)(i)(n) << ' ';
				}
				for(n=0;n<NV;++n) {
					out << ugbd(tlvl).v(i,n)<< ' ';
				}
				
				out << std::endl;
			}
			
			if (basis::tet(log2p).p > 1) {
				/* EDGE MODES */
				for(eind=0;eind<nseg;++eind) {
					if (seg(eind).info < 0) {
						v0 = seg(eind).pnt(0);
						v1 = seg(eind).pnt(1);
						for(n=0;n<ND;++n)
							basis::tet(log2p).proj1d_leg(vrtxbd(tlvl)(v0)(n),vrtxbd(tlvl)(v1)(n),&crd1d(n)(0));
					}
					else {
						crdtocht1d(eind,tlvl);
						for(n=0;n<ND;++n)
							basis::tet(log2p).proj1d_leg(&cht(n)(0),&crd1d(n)(0));
					}
					ugtouht1d(eind,tlvl);
					for(n=0;n<NV;++n)
						basis::tet(log2p).proj1d_leg(&uht(n)(0),&u1d(n)(0));
					
					for(i=1;i<basis::tet(log2p).em+1;++i) {
						for(n=0;n<ND;++n) {
							out << crd1d(n)(i) << ' ';
						}
						for(n=0;n<NV;++n){
							out << u1d(n)(i)<< ' ';
						}
						out << std::endl;
					}
				}
			}
			
			/* FACE MODES */
			if (basis::tet(log2p).p > 2) {
				for(find = 0; find < ntri; ++find) {
					ugtouht2d(find,tlvl);
					for(n=0;n<NV;++n)
						basis::tet(log2p).proj2d_leg(&uht(n)(0),&u2d(n)(0)(0),MXGP);
					
					if (tri(find).info < 0) {
						for(n=0;n<ND;++n)
							basis::tet(log2p).proj2d_leg(vrtxbd(tlvl)(tri(find).pnt(0))(n),vrtxbd(tlvl)(tri(find).pnt(1))(n),vrtxbd(tlvl)(tri(find).pnt(2))(n),&crd2d(n)(0)(0),MXGP);
					}
					else {
						crdtocht2d(find,tlvl);
						for(n=0;n<ND;++n)
							basis::tet(log2p).proj2d_bdry_leg(&cht(n)(0),&crd2d(n)(0)(0),MXGP);
					}
					
					for(i=1;i<basis::tet(log2p).em;++i) {
						for(j=1;j<basis::tet(log2p).em-(i-1);++j) {
							for(n=0;n<ND;++n) {
								out << crd2d(n)(i)(j) << ' ';
							}
							//							u2d(0)(i)(j)-=1.0;
							for(n=0;n<NV;++n) {
								out << u2d(n)(i)(j) << ' ';
							}
							out << std::endl;
						}
					}
				}
			}
			
			/* INTERIOR MODES */
			if (basis::tet(log2p).p > 3) {
				for(tind = 0; tind < ntet; ++tind) {
					ugtouht(tind,tlvl);
					for(n=0;n<NV;++n)
						basis::tet(log2p).proj_leg(&uht(n)(0),&u(n)(0)(0)(0),stridex,stridey);
					
					if (tri(find).info < 0) {
						for(n=0;n<ND;++n)
							basis::tet(log2p).proj_leg(vrtxbd(tlvl)(tet(tind).pnt(0))(n),vrtxbd(tlvl)(tet(tind).pnt(1))(n),vrtxbd(tlvl)(tet(tind).pnt(2))(n),vrtxbd(tlvl)(tet(tind).pnt(3))(n),&crd(n)(0)(0)(0),stridex,stridey);
					}
					else {
						crdtocht(tind,tlvl);
						for(n=0;n<ND;++n)
							basis::tet(log2p).proj_bdry_leg(&cht(n)(0),&crd(n)(0)(0)(0),stridex,stridey);
					}
					
					for(i=1;i<basis::tet(log2p).em;++i) {
						for(j=1;j<basis::tet(log2p).em+1-i;++j) {
							for(int k = 1; k < basis::tet(log2p).em+1-i-j; ++k) {
								for(n=0;n<ND;++n) {
									out << crd(n)(i)(j)(k) << ' ';
								}
								for(n=0;n<NV;++n) {
									out << u(n)(i)(j)(k)<< ' ';
								}
								out << std::endl;
							}
						}
					}
				}
			}
			
			/* OUTPUT CONNECTIVY INFO */
			for(tind=0;tind<ntet;++tind) {
				
				/* VERTICES */
				ijind[0][0][basis::tet(log2p).em+1] = tet(tind).pnt(0);
				ijind[0][basis::tet(log2p).em+1][0] = tet(tind).pnt(1);
				ijind[0][0][0] = tet(tind).pnt(2);
				ijind[basis::tet(log2p).em+1][0][0] = tet(tind).pnt(3);
				
				/* EDGES */
				indx = tet(tind).seg(0);
				sgn = tet(tind).sgn(0);
				if (sgn < 0) {
					for(i=0;i<basis::tet(log2p).em;++i)
						ijind[i+1][0][0] = npnt +(indx+1)*basis::tet(log2p).em -(i+1);
				}
				else {
					for(i=0;i<basis::tet(log2p).em;++i)
						ijind[i+1][0][0] = npnt +indx*basis::tet(log2p).em +i;
				}
				
				indx = tet(tind).seg(1);
				sgn = tet(tind).sgn(1);
				if (sgn > 0) {
					for(i=0;i<basis::tet(log2p).em;++i)
						ijind[basis::tet(log2p).em-i][i+1][0] = npnt +indx*basis::tet(log2p).em +i;
				}
				else {
					for(i=0;i<basis::tet(log2p).em;++i)
						ijind[basis::tet(log2p).em-i][i+1][0] = npnt +(indx+1)*basis::tet(log2p).em -(i+1);
				}
				
				indx = tet(tind).seg(2);
				sgn = tet(tind).sgn(2);
				if (sgn < 0) {
					for(i=0;i<basis::tet(log2p).em;++i)
						ijind[0][i+1][0] = npnt +(indx+1)*basis::tet(log2p).em -(i+1);
				}
				else {
					for(i=0;i<basis::tet(log2p).em;++i)
						ijind[0][i+1][0] = npnt +indx*basis::tet(log2p).em +i;
				}
				
				indx = tet(tind).seg(3);
				sgn = tet(tind).sgn(3);
				if (sgn > 0) {
					for(i=0;i<basis::tet(log2p).em;++i)
						ijind[0][i+1][basis::tet(log2p).em-i] = npnt +(indx+1)*basis::tet(log2p).em -(i+1);
				}
				else {
					for(i=0;i<basis::tet(log2p).em;++i)
						ijind[0][i+1][basis::tet(log2p).em-i] = npnt +indx*basis::tet(log2p).em +i;
				}
				
				indx = tet(tind).seg(4);
				sgn = tet(tind).sgn(4);
				if (sgn > 0) {
					for(i=0;i<basis::tet(log2p).em;++i)
						ijind[0][0][i+1] = npnt +indx*basis::tet(log2p).em +i;
				}
				else {
					for(i=0;i<basis::tet(log2p).em;++i)
						ijind[0][0][i+1] = npnt +(indx+1)*basis::tet(log2p).em -(i+1);
				}
				
				indx = tet(tind).seg(5);
				sgn = tet(tind).sgn(5);
				if (sgn > 0) {
					for(i=0;i<basis::tet(log2p).em;++i)
						ijind[i+1][0][basis::tet(log2p).em-i] = npnt +(indx+1)*basis::tet(log2p).em -(i+1);
				}
				else {
					for(i=0;i<basis::tet(log2p).em;++i)
						ijind[i+1][0][basis::tet(log2p).em-i] = npnt +indx*basis::tet(log2p).em +i;
				}
				
				/* FACES */
				k = 0;
				indx = tet(tind).tri(0);
				rot = -tet(tind).rot(0);
				if (rot == -1){
					for(i=1;i<basis::tet(log2p).em;++i) {
						for(j=1;j<basis::tet(log2p).em-(i-1);++j) {
							ijind[basis::tet(log2p).em-i-j+1][j][0] = npnt +nseg*basis::tet(log2p).em +indx*basis::tet(log2p).fm +k;
							++k;
						}
					}
				}
				if (rot == 1) {
					for(i=1;i<basis::tet(log2p).em;++i) {
						for(j=1;j<basis::tet(log2p).em-(i-1);++j) {
							ijind[i][j][0] = npnt +nseg*basis::tet(log2p).em +indx*basis::tet(log2p).fm +k;
							++k;
						}
					}
				}
				
				k = 0;
				indx = tet(tind).tri(1);
				rot = -tet(tind).rot(1);
				if(rot == 1){
					for(i=1;i<basis::tet(log2p).em;++i) {
						for(j=1;j<basis::tet(log2p).em-(i-1);++j) {
							ijind[basis::tet(log2p).em-i-j+1][0][j] = npnt +nseg*basis::tet(log2p).em +indx*basis::tet(log2p).fm +k;
							++k;
						}
					}
				}
				if (rot == -1) {
					for(i=1;i<basis::tet(log2p).em;++i) {
						for(j=1;j<basis::tet(log2p).em-(i-1);++j) {
							ijind[i][0][j] = npnt +nseg*basis::tet(log2p).em +indx*basis::tet(log2p).fm +k;
							++k;
						}
					}
				}
				
				k = 0;
				indx = tet(tind).tri(2);
				rot = -tet(tind).rot(2);
				if (rot == 1){
					for(i=1;i<basis::tet(log2p).em;++i) {
						for(j=1;j<basis::tet(log2p).em-(i-1);++j) {
							ijind[i][basis::tet(log2p).em-i-j+1][j] = npnt +nseg*basis::tet(log2p).em +indx*basis::tet(log2p).fm +k;
							++k;
						}
					}
				}
				if (rot == -1) {
					for(i=1;i<basis::tet(log2p).em;++i) {
						for(j=1;j<basis::tet(log2p).em-(i-1);++j) {
							ijind[basis::tet(log2p).em-i-j+1][i][j] = npnt +nseg*basis::tet(log2p).em +indx*basis::tet(log2p).fm +k;
							++k;
						}
					}
				}
				
				k = 0;
				indx = tet(tind).tri(3);
				rot = -tet(tind).rot(3);
				if (rot == -1){
					for(i=1;i<basis::tet(log2p).em;++i) {
						for(j=1;j<basis::tet(log2p).em-(i-1);++j) {
							ijind[0][basis::tet(log2p).em-i-j+1][j] = npnt +nseg*basis::tet(log2p).em +indx*basis::tet(log2p).fm +k;
							++k;
						}
					}
				}
				if (rot == 1) {
					for(i=1;i<basis::tet(log2p).em;++i) {
						for(j=1;j<basis::tet(log2p).em-(i-1);++j) {
							ijind[0][i][j] = npnt +nseg*basis::tet(log2p).em +indx*basis::tet(log2p).fm +k;
							++k;
						}
					}
				}
				
				/* INTERIOR */
				m = 0;
				for(i=1;i<basis::tet(log2p).em;++i) {
					for(j=1;j<basis::tet(log2p).em-(i-1);++j) {
						for(k=1;k<basis::tet(log2p).em-i-j+1;++k){
							ijind[i][j][k] = npnt +nseg*basis::tet(log2p).em +ntri*basis::tet(log2p).fm +tind*basis::tet(log2p).im +m;
							++m;
						}
					}
				}
				
				
				/* OUTPUT CONNECTION LIST */
				int em = basis::tet(log2p).em;
				for(i=0;i<em;++i) {
					for(j=0;j<em-i;++j) {
						for(k=0;k<em-i-j-1;++k) {
							out << ijind[i][j][k]+1 << ' ' << ijind[i][j][k+1]+1 << ' ' << ijind[i+1][j][k]+1 <<  ' ' << ijind[i][j+1][k]+1 << std::endl;
							out << ijind[i][j][k+1]+1 << ' ' << ijind[i+1][j][k+1]+1 << ' ' << ijind[i+1][j][k]+1 <<  ' ' << ijind[i][j+1][k]+1 << std::endl;
							out << ijind[i][j][k+1]+1 << ' ' << ijind[i+1][j][k+1]+1 << ' ' << ijind[i][j+1][k]+1 <<  ' ' << ijind[i][j+1][k+1]+1 << std::endl;
							out << ijind[i+1][j][k+1]+1 << ' ' << ijind[i][j+1][k]+1 << ' ' << ijind[i+1][j+1][k]+1 <<  ' ' << ijind[i+1][j][k]+1 << std::endl;
							out << ijind[i+1][j][k+1]+1 << ' ' << ijind[i][j+1][k]+1 << ' ' << ijind[i][j+1][k+1]+1 <<  ' ' << ijind[i+1][j+1][k]+1 << std::endl;
							out << ijind[i][j+1][k+1]+1 << ' ' << ijind[i+1][j][k+1]+1 << ' ' << ijind[i+1][j+1][k]+1 <<  ' ' << ijind[i+1][j+1][k+1]+1 << std::endl;
							
						}
						// five tet new way
						out << ijind[i][j][em-j-1-i]+1 << ' ' << ijind[i][j][em-j-i]+1 << ' ' << ijind[i+1][j][em-j-1-i]+1 <<  ' ' << ijind[i][j+1][em-j-1-i]+1 << std::endl;
						out << ijind[i][j][em-j-i]+1 << ' ' << ijind[i+1][j][em-j-i]+1 << ' ' << ijind[i+1][j][em-j-1-i]+1 <<  ' ' << ijind[i][j+1][em-j-1-i]+1 << std::endl;
						out << ijind[i][j][em-j-i]+1 << ' ' << ijind[i+1][j][em-j-i]+1 << ' ' << ijind[i][j+1][em-j-1-i]+1 <<  ' ' << ijind[i][j+1][em-j-i]+1 << std::endl;
						out << ijind[i+1][j][em-j-i]+1 << ' ' << ijind[i][j+1][em-j-1-i]+1 << ' ' << ijind[i+1][j+1][em-j-1-i]+1 <<  ' ' << ijind[i+1][j][em-j-1-i]+1 << std::endl;
						out << ijind[i+1][j][em-j-i]+1 << ' ' << ijind[i][j+1][em-j-1-i]+1 << ' ' << ijind[i][j+1][em-j-i]+1 <<  ' ' << ijind[i+1][j+1][em-j-1-i]+1 << std::endl;
						
						// single tet
						out << ijind[i][j][em-j-i]+1 << ' ' << ijind[i+1][j][em-j-i]+1 << ' ' << ijind[i][j+1][em-j-i]+1 <<  ' ' << ijind[i][j][em-j+1-i]+1 << std::endl;
						
					}
					//single tet
					out << ijind[i][em-i][0]+1 << ' ' << ijind[i+1][em-i][0]+1 << ' ' << ijind[i][em+1-i][0]+1 <<  ' ' << ijind[i][em-i][1]+1 << std::endl;
				}
				//single tet
				out << ijind[em][0][0]+1 << ' ' << ijind[em+1][0][0]+1 << ' ' << ijind[em][1][0]+1 <<  ' ' << ijind[em][0][1]+1 << std::endl;
			}
			out.close();
			break; 
		}

			
		default: {
			*gbl->log << "can't output a tet_hp to that filetype" << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
			break;
		}
	}
	
	/* BOUNDARY INFO */
	for(i=0;i<nfbd;++i) {
		hp_fbdry(i)->output(filename,typ,tlvl);
	}
	
	for(i=0;i<nebd;++i) {
		hp_ebdry(i)->output(filename,typ,tlvl);
	}
	
	for(i=0;i<nvbd;++i) {
		hp_vbdry(i)->output(filename,typ,tlvl);
	}
	
	return;
}

void tet_hp::input(const std::string& filename) {
	int i,j;
	std::string fname,fnmapp;
	std::ostringstream nstr;
	ifstream fin;
	binifstream bin;
	
	fname = filename +"_" +gbl->idprefix;
	
	if (reload_type == tet_hp::binary) {
		fnmapp = fname +".bin";
		fin.open(fnmapp.c_str(),ios::in);
		if(fin.is_open()) {
			fin.close();
			input_map blank;
			tet_mesh::input(fnmapp,tet_mesh::binary,1,blank);
			for(i=1;i<gbl->nadapt;++i) {
				nstr.str("");
				nstr << i << std::flush;
				fnmapp = filename +"_v" +nstr.str() +"_" +gbl->idprefix +".bin";
				bin.open(fnmapp.c_str());
				if (bin.error()) {
					*gbl->log << "couldn't open input file " << fnmapp << std::endl;
					vrtxbd(i)(Range(0,npnt-1)) = pnts(Range(0,npnt-1));
				}
				bin.setFlag(binio::BigEndian,bin.readInt(1));
				bin.setFlag(binio::FloatIEEE,bin.readInt(1));
				for (j=0;j<npnt;++j) {
					vrtxbd(i)(j)(0) = bin.readFloat(binio::Double);
					vrtxbd(i)(j)(1) = bin.readFloat(binio::Double);
					vrtxbd(i)(j)(2) = bin.readFloat(binio::Double);
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
			fnmapp = fname +".v";
			input_map blank;
			tet_mesh::input(fnmapp,tet_mesh::grid,1,blank);
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
					fin >> vrtxbd(i)(j)(0) >> vrtxbd(i)(j)(1) >> vrtxbd(i)(j)(2);
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
			fnmapp = filename +".d" +nstr.str();
			input(fnmapp,reload_type,i);
		}
	}

	return;

}

void tet_hp::input(const std::string& filename, filetype typ, int tlvl) {
	int i,k,m,n,pin,pmin,indx;
	std::string fnapp;
	char buffer[80];
	ifstream in;
	FLT fltskip;
	
	std::string fname;
	fname = filename +"_" +gbl->idprefix;

	switch(typ) {
		case (text): {
			fnapp = filename +".txt";
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
						in >> ugbd(tlvl).e(i,m,n);
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
							in >> ugbd(tlvl).f(i,indx,n);
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

			/* FIXME: NOT WORKING FOR CHANGES OF ORDER */
			for(i=0;i<ntet;++i) {
				for(m=0;m<im0;++i) {
					for(n=0;n<NV;++n) 
						in >> ugbd(tlvl).i(i,m,n);
				}
			}
			in.ignore(80,'\n');

			/* BOUNDARY INFO */
			for(i=0;i<nfbd;++i)
				hp_fbdry(i)->input(in,typ,tlvl);

			for(i=0;i<nebd;++i)
				hp_ebdry(i)->input(in,typ,tlvl);

			for(i=0;i<nvbd;++i)
				hp_vbdry(i)->input(in,typ,tlvl);

			in.close();
			break;
		}
		
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
				*gbl->log << "mismatched pnt counts?\n";
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			if (bin.readInt(sizeof(int))  != nseg) {
				*gbl->log << "mismatched seg counts?\n";
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			if (bin.readInt(sizeof(int))  != ntri) {
				*gbl->log << "mismatched tri counts?\n";
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			if (bin.readInt(sizeof(int))  != ntet) {
				*gbl->log << "mismatched tet counts?\n";
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			
			for(i=0;i<npnt;++i) {
				for(n=0;n<NV;++n)
					ugbd(tlvl).v(i,n) = bin.readFloat(binio::Double);
			}

			for(i=0;i<nseg;++i) {
				for(m=0;m<(pmin-1);++m) {
					for(n=0;n<NV;++n)
						ugbd(tlvl).e(i,m,n) = bin.readFloat(binio::Double);
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
							ugbd(tlvl).f(i,indx,n) = bin.readFloat(binio::Double);
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
			
			/* FIXME: NOT WORKING FOR CHANGES OF ORDER */
			for(i=0;i<ntet;++i) {
				for(m=0;m<im0;++i) {
					for(n=0;n<NV;++n) 
						ugbd(tlvl).i(i,m,n) = bin.readFloat(binio::Double);
				}
			}

			/* BOUNDARY INFO */
			for(i=0;i<nfbd;++i)
				hp_fbdry(i)->input(in,typ,tlvl);

			for(i=0;i<nebd;++i)
				hp_ebdry(i)->input(in,typ,tlvl);

			for(i=0;i<nvbd;++i)
				hp_vbdry(i)->input(in,typ,tlvl);

			in.close();
			break;
		}                    
		default:
			*gbl->log << "can't input a tet_hp from that filetype" << std::endl;
			sim::abort(__LINE__,__FILE__,gbl->log);
			break;
	}
	
	return;
}
