#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include "tet_mesh.h"
#ifdef USING_MADLIB
#include "MAdLibInterface.h"
#endif
// #define DATATANK

#ifdef DATATANK
#include <DTSource.h>
#endif

using namespace std;

void tet_mesh::output(const std::string &filename, tet_mesh::filetype filetype) const {
	std::string fnmapp;
	ofstream out;
	//bostream bout;

	
	out.setf(std::ios::scientific, std::ios::floatfield);
	out.precision(10);
			
	switch (filetype) {    
		case (easymesh): 
			fnmapp = filename +".p";
			out.open(fnmapp.c_str());
			if (!out) {
				*gbl->log << "couldn't open output file " << fnmapp << "for output" << endl;
				exit(1);
			}
			
			out << npnt << endl;
			for(int i=0;i<npnt;++i) {
				out << i << ": ";
				for(int n=0;n<ND;++n)
					out << pnts(i)(n) << ' ';
				out << pnt(i).info << endl;
			}
			out.close();

			/* SIDE FILE */
			fnmapp = filename +".s";
			out.open(fnmapp.c_str());
			if (!out) {
				*gbl->log << "couldn't open output file " << fnmapp << "for output" << endl;
				exit(1);
			}
			out << nseg << endl;
			for(int i=0;i<nseg;++i) {
				out << i << ": " << seg(i).pnt(0) << ' ' << seg(i).pnt(1) << ' ' << seg(i).info <<endl;
			}
			out.close();

			fnmapp = filename +".f";
			out.open(fnmapp.c_str());
			if (!out) {
				*gbl->log << "couldn't open output file " << fnmapp << "for output" << endl;
				exit(1);
			}
			out << ntri << endl;
			for(int i=0;i<ntri;++i) {
				out << i << ": " << tri(i).pnt(0) << ' ' << tri(i).pnt(1) << ' ' << tri(i).pnt(2);
				out << ' ' << tri(i).seg(0) << ' ' << tri(i).seg(1) << ' ' << tri(i).seg(2);
				out << ' ' << tri(i).sgn(0) << ' ' << tri(i).sgn(1) << ' ' << tri(i).sgn(2);
				out << ' ' << tri(i).tet(0) << ' ' << tri(i).tet(1) << ' ' << tri(i).info << endl;
			}
			out.close();

			fnmapp = filename +".t";
			out.open(fnmapp.c_str());
			if (!out) {
				*gbl->log << "couldn't open output file " << fnmapp << "for output" << endl;
				exit(1);
			}
			out << ntet << endl;
			for(int i=0;i<ntet;++i) {
				out << i << ": " << tet(i).pnt(0) << ' ' << tet(i).pnt(1) << ' ' << tet(i).pnt(2) << ' ' << tet(i).pnt(3);
				out << ' ' << tet(i).seg(0) << ' ' << tet(i).seg(1) << ' ' << tet(i).seg(2) << ' ' << tet(i).seg(3) << ' ' << tet(i).seg(4) << ' ' << tet(i).seg(5);
				out << ' ' << tet(i).sgn(0) << ' ' << tet(i).sgn(1) << ' ' << tet(i).sgn(2) << ' ' << tet(i).sgn(3) << ' ' << tet(i).sgn(4) << ' ' << tet(i).sgn(5);
				out << ' ' << tet(i).tri(0) << ' ' << tet(i).tri(1) << ' ' << tet(i).tri(2) << ' ' << tet(i).tri(3);
				out << ' ' << tet(i).rot(0) << ' ' << tet(i).rot(1) << ' ' << tet(i).rot(2) << ' ' << tet(i).rot(3);				
				out << ' ' << tet(i).tet(0) << ' ' << tet(i).tet(1) << ' ' << tet(i).tet(2) << ' ' << tet(i).tet(3) << ' ' << tet(i).info << endl;
			}
			out.close();			
			
			
			
			break;


		case (tecplot):
			fnmapp = filename +".dat";
			out.open(fnmapp.c_str());
			if (!out) {
				*gbl->log << "couldn't open output file " << fnmapp << "for output" << endl;
				exit(1);
			}
			out << "VARIABLES = \"X\", \"Y\", \"Z\" " << endl;
			out << "ZONE F=FEPOINT, ET=TETRAHEDRON, N = " << npnt << ", E = " << ntet << endl;
			
			for(int i=0;i<npnt;++i) {
				for(int n=0;n<ND;++n)
					out << pnts(i)(n) << ' ';
				out << endl;
			}

			out << "\n#CONNECTION DATA#\n";
			
			//tet_mesh::vertexball(50);
//            for(i=0;i<pnt(50).nnbor;++i){
//                tind = i2wk(i);
//                out << tet(tind).pnt(0)+1 << ' ' << tet(tind).pnt(1)+1 << ' ' << tet(tind).pnt(2)+1 << ' ' << tet(tind).pnt(3)+1 << endl;
//            }
			for(int i=0;i<ntet;++i)
				out << tet(i).pnt(0)+1 << ' ' << tet(i).pnt(1)+1 << ' ' << tet(i).pnt(2)+1 << ' ' << tet(i).pnt(3)+1 << endl;
				
			out.close();
			break;

		case(text):
			/* JUST OUTPUT VERTEX POSITIONS FOR DEFORMING MESH */
			fnmapp = filename +".txt";
			out.open(fnmapp.c_str());
			if (!out) {
				*gbl->log << "couldn't open output file" << fnmapp << "for output" << endl;
				exit(1);
			}
			out << npnt << endl;
			for(int i=0;i<npnt;++i) {
				out << i << ":";
				for(int n=0;n<ND;++n)
					out << pnts(i)(n) << ' '; 
				out << "\n";
			}
				
			out.close();
			break;
			
		case (grid):
		
			fnmapp = filename +".grd";
			out.open(fnmapp.c_str());

			out << "npnt:" << ' ' << npnt << "  nseg:" << ' ' << nseg << "  ntri:" << ' ' << ntri << "  ntet:" << ' ' << ntet << endl;
			
			// output vertices
			for(int i = 0; i < npnt; ++i)
				out << i << ": " << pnts(i)(0) << ' ' << pnts(i)(1) << ' ' << pnts(i)(2) << endl;
			
			// output seg connection data 
			for(int i = 0; i < nseg; ++i){
				out << i << ": ";
				for(int j = 0; j < 2; ++j)
					out << seg(i).pnt(j) << ' ' ;
				out << endl;
			}			
			
			// output face connection data
			for(int i = 0; i < ntri; ++i){
				out << i << ": ";
				for(int j = 0; j < 3; ++j)
					out << tri(i).pnt(j) << ' ' ;
				out << endl;
			}

			// output tet connection data
			for(int i = 0; i < ntet; ++i){
				out << i << ": ";
				for(int j = 0; j < 4; ++j)
					out << tet(i).pnt(j) << ' ' ;
				out << endl;
			}
				
						
			/* VERTEX BOUNDARY INFO HEADER */
			out << "nvbd: " << nvbd << endl;
			for(int i=0;i<nvbd;++i) {
				out << "idnum: " << vbdry(i)->idnum << endl;
				out << "point: " << vbdry(i)->pnt << endl;
			}
			
			/* SIDE BOUNDARY INFO HEADER */
			out << "nebd: " << nebd << endl;
			for(int i=0;i<nebd;++i) {
				out << "idnum: " << ebdry(i)->idnum << endl;
				out << "number: " << ebdry(i)->nseg << endl;
				for(int j=0;j<ebdry(i)->nseg;++j)
					out << j << ": " << ebdry(i)->seg(j).gindx << std::endl;
			}    
			
			/* FACE BOUNDARY INFO HEADER */
			out << "nfbd: " << nfbd << endl ;
			for(int i = 0; i < nfbd; ++i){
				out << "idnum: " << fbdry(i)->idnum << endl;
				out << "npnt: " << fbdry(i)->npnt << " nseg: " << fbdry(i)->nseg<< " ntri: " << fbdry(i)->ntri << endl;
				//out << fbdry(i)->type
				for(int j = 0; j < fbdry(i)->npnt; ++j){
					out << j << ": " << fbdry(i)->pnt(j).gindx << endl;
				}
				for(int j = 0; j < fbdry(i)->nseg; ++j){
					out << j << ": " << fbdry(i)->seg(j).gindx << endl;
				}
				for(int j = 0; j < fbdry(i)->ntri; ++j){
					out << j << ": " << fbdry(i)->tri(j).gindx << endl; 
				}
			}

			
			out.close();
			break;

		case (datatank): {
#ifdef DATATANK
			DTMutableIntArray dt_tvrtx(4,ntet);
			DTMutableDoubleArray dt_vrtx(3,npnt);

			for (int tind=0;tind<ntet;++tind)
				for (int j=0;j<4;++j)
					dt_tvrtx(j,tind) = tet(tind).pnt(j);

			for (int vind=0;vind<npnt;++vind)
				for (int j=0;j<3;++j)
					dt_vrtx(j,vind) = pnts(vind)(j);
			// find out what to put here        
			DTTriangularGrid2D dt_grid(dt_tvrtx,dt_vrtx);
			
			std::string outputFilename(filename +".dt");
			DTDataFile outputFile(outputFilename.c_str(),DTFile::NewReadWrite);
			// Output from computation
			Write(outputFile,"grid",dt_grid);
			outputFile.Save("TetrahedralGrid3D","Seq_grid");
 #else
			*gbl->log << "Not supported on this platform\n";
#endif
			break;
		}
			
			
		case (gmsh): {
#ifdef USING_MADLIB
			MAdLib_output(filename);
#else
			*gbl->log << "gmsh Not supported on this platform\n";
#endif
			break;
		}
			
		case(vlength): {
			fnmapp = filename +".lngth";
			out.open(fnmapp.c_str());
			if (!out) {
				*gbl->log << "couldn't open output file" << fnmapp << "for output" << endl;
				exit(1);
			}

			for(int i=0;i<npnt;++i)
				out << lngth(i) << endl;
				
			break;
		}


		default:
			*gbl->log << "That filetype is not supported yet" << endl;
			break;
	}
	
	return;
}


/* partition requires this routine */
void tet_mesh::setinfo() {
	
	/* SET UP VRTX BC INFORMATION FOR OUTPUT */
	for(int i=0;i<npnt;++i)
		pnt(i).info = -1;
	
	for(int i=0;i<nvbd;++i)
		pnt(vbdry(i)->pnt).info = vbdry(i)->idnum;

	/* SET UP SIDE BC INFORMATION FOR EASYMESH OUTPUT */
	for(int i=0;i<nseg;++i)
		seg(i).info = -1;
	
	for(int i=0;i<nebd;++i)
		for(int j=0;j<ebdry(i)->nseg;++j)
			seg(ebdry(i)->seg(j).gindx).info = ebdry(i)->idnum;

	for(int i=0;i<ntri;++i)
		tri(i).info = -1;
		
	for(int i=0;i<nfbd;++i)
		for(int j=0;j<fbdry(i)->ntri;++j)
			tri(fbdry(i)->tri(j).gindx).info = fbdry(i)->idnum;
	
	/* SET UP TRI INFO FOR EASYMESH OUTPUT */            
	for(int i=0;i<ntet;++i)
		tet(i).info = -1;
			
	return;
}
