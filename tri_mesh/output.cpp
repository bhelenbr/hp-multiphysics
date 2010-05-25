#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <libbinio/binfile.h>
#include "tri_mesh.h"

// #define DATATANK

#ifdef DATATANK
#include <DTSource.h>
#endif

using namespace std;

void tri_mesh::output(const std::string &filename, tri_mesh::filetype filetype) const {
	std::string fnmapp, grd_nm;
	int i,j,n,tind,count;
	ofstream out;
	binofstream bout;

	out.setf(std::ios::scientific, std::ios::floatfield);
	out.precision(10);
	
	grd_nm = filename;
	
	/* Override filetype based on ending? */
	size_t dotloc;
	dotloc = grd_nm.find_last_of('.');
	string ending;
	ending = grd_nm.substr(dotloc+1);
	if (ending == "grd") {
		filetype = grid;
		grd_nm = grd_nm.substr(0,dotloc);
	}
	else if (ending == "d") {
		filetype = boundary;
		grd_nm = grd_nm.substr(0,dotloc);
	}
	else if (ending == "bin") {
		filetype = binary;
		grd_nm = grd_nm.substr(0,dotloc);
	}
	else if (ending == "dat") {
		filetype = tecplot;
		grd_nm = grd_nm.substr(0,dotloc);
	}		
	

	switch (filetype) {

		case (easymesh):
			/* CREATE EASYMESH OUTPUT FILES */
			fnmapp = grd_nm +".n";
			out.open(fnmapp.c_str());
			if (!out) {
				*gbl->log << "couldn't open output file " << fnmapp << "for output" << endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}

			out << npnt << endl;
			for(i=0;i<npnt;++i) {
				out << i << ": ";
				for(n=0;n<ND;++n)
					out << pnts(i)(n) << ' ';
				out << pnt(i).info << endl;
			}
			out.close();

			/* SIDE FILE */
			fnmapp = grd_nm +".s";
			out.open(fnmapp.c_str());
			if (!out) {
				*gbl->log << "couldn't open output file " << fnmapp << "for output" << endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			out << nseg << endl;
			for(i=0;i<nseg;++i) {
				out << i << ": " << seg(i).pnt(0) << ' ' << seg(i).pnt(1) << ' ';
				out << seg(i).tri(0) << ' ' << seg(i).tri(1) << ' ' << seg(i).info <<endl;
			}
			out.close();

			fnmapp = grd_nm +".e";
			out.open(fnmapp.c_str());
			if (!out) {
				*gbl->log << "couldn't open output file " << fnmapp << "for output" << endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			out << ntri << endl;
			for(i=0;i<ntri;++i) {
				out << i << ": " << tri(i).pnt(0) << ' ' << tri(i).pnt(1) << ' ' << tri(i).pnt(2);
				out << ' ' << tri(i).tri(0) << ' ' << tri(i).tri(1) << ' ' << tri(i).tri(2);
				out << ' ' << tri(i).seg(0) << ' ' << tri(i).seg(1) << ' ' << tri(i).seg(2);
				out << " 0.0 0.0 " << tri(i).info << endl;
			}
			out.close();
			break;

		case (tecplot):
			fnmapp = grd_nm +".dat";
			out.open(fnmapp.c_str());
			if (!out) {
				*gbl->log << "couldn't open output file " << fnmapp << "for output" << endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}

			out << "ZONE F=FEPOINT, ET=TRIANGLE, N = " << npnt << ", E = " << ntri << endl;

			for(i=0;i<npnt;++i) {
				for(n=0;n<ND;++n)
					out << pnts(i)(n) << ' ';
				out << endl;
			}

			out << "\n#CONNECTION DATA#\n";

			for(i=0;i<ntri;++i)
				out << tri(i).pnt(0)+1 << ' ' << tri(i).pnt(1)+1 << ' ' << tri(i).pnt(2)+1 << endl;

			out.close();
			break;

		case (mavriplis):
			fnmapp = grd_nm +".GDS";
			out.open(fnmapp.c_str());
			if (!out) {
				*gbl->log << "couldn't open output file " << fnmapp << "for output" << endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}

			out << "YIN " << npnt << endl;

			for(i=0;i<npnt;++i) {
				for(n=0;n<ND;++n)
					out << "         " << pnts(i)(n);
				if (ND == 2) out << "         " << 0.0;
				out << endl;
			}

			out << "0\n";
			out << ntri << '\n';

			for(i=0;i<ntri;++i) {
				out << "1 Default" << endl;
				out << 3 << ' ' << tri(i).pnt(0)+1 << ' ' << tri(i).pnt(1)+1 << ' ' << tri(i).pnt(2)+1 << endl;
			}
				out.close();
			break;

		case (gambit):
			fnmapp = grd_nm +".FDNEUT";
			out.open(fnmapp.c_str());
			if (!out) {
				*gbl->log << "couldn't open output file " << fnmapp << "for output" << endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}

			count = ntri;
			for(i=0;i<nebd;++i)
				count += ebdry(i)->nseg;

			out << "** FIDAP NEUTRAL FILE\n";
			out << grd_nm << "\n";
			out << "VERSION     8.01\n";
			out << "29 Nov 1999     13:23:58\n";
			out << "    NO. OF NODES    NO. ELEMENTS NO. ELT GROUPS             NDFCD             NDFVL\n";

			out << setw(15) << npnt << setw(15) << count << setw(15) << nebd+1 << setw(15) << ND << setw(15) << ND << endl;
			out << "    STEADY/TRANS      TURB. FLAG FREE SURF FLAG     COMPR. FLAG    RESULTS ONLY\n";
			out << setw(12) << 0 << setw(12) << 0 << setw(12) << 0 << setw(12) << 0 << setw(12) << 0 << endl;
			out << "TEMPERATURE/SPECIES FLAGS\n";
			out << " 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n";
			out << "PRESSURE FLAGS - IDCTS, IPENY MPDF\n";
			out << "            1            1            0\n";
			out << "NODAL COORDINATES\n";

			for(i=0;i<npnt;++i) {
				out << setw(10) << i+1;
				for(n=0;n<ND;++n)
					out << setw(20) << pnts(i)(n);
				out << '\n';
			}

			out << "BOUNDARY CONDITIONS\n";
			out << "            0            0            0      0.0\n";
			out << "ELEMENT GROUPS\n";
			out << "GROUP:" << setw(9) << 1;
			out << " ELEMENTS:" << setw(10) << ntri;
			out << " NODES:" << setw(13) << 3;
			out << " GEOMETRY:" << setw(5) << 2;
			out << " TYPE:" << setw(4) << 2 << endl;
			out << "ENTITY NAME:    fluid\n";

			for(tind=0;tind<ntri;++tind) {
				out << setw(8) << tind+1;
				out << setw(8) << tri(tind).pnt(0)+1;
				out << setw(8) << tri(tind).pnt(1)+1;
				out << setw(8) << tri(tind).pnt(2)+1;
				out << endl;
			}

			count = ntri+1;

			for(i=0;i<nebd;++i) {
				out << "GROUP:" << setw(9) << i+2;
				out << " ELEMENTS:" << setw(10) << ebdry(i)->nseg;
				out << " NODES:" << setw(13) << 2;
				out << " GEOMETRY:" << setw(5) << 0;
				out << " ID:" << setw(4) << ebdry(i)->idnum << endl;

				out << "ENTITY NAME:    surface.1\n";
				for(j=0;j<ebdry(i)->nseg;++j) {
					out << setw(8) << count++;
					out << setw(8) << seg(ebdry(i)->seg(j)).pnt(0)+1;
					out << setw(8) << seg(ebdry(i)->seg(j)).pnt(1)+1 << endl;
				}
			}
			out.close();
			break;

		case(text):
			/* JUST OUTPUT VERTEX POSITIONS FOR DEFORMING MESH */
			fnmapp = grd_nm +".txt";
			out.open(fnmapp.c_str());
			if (!out) {
				*gbl->log << "couldn't open output file" << fnmapp << "for output" << endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			out << npnt << endl;
			for(i=0;i<npnt;++i) {
				out << i << ":";
				for(n=0;n<ND;++n)
					out << pnts(i)(n) << ' ';
				out << "\n";
			}

			out.close();
			break;

		case(grid):
			fnmapp = grd_nm +".grd";
			out.open(fnmapp.c_str());
			if (!out) {
				*gbl->log << "couldn't open output file" << fnmapp << "for output" << endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}

			/* HEADER LINES */
			out << "npnt: " << npnt;
			out << " nseg: " << nseg;
			out << " ntri: " << ntri << endl;

			/* POINT INFO */
			for(i=0;i<npnt;++i) {
				out << i << ": ";
				for(n=0;n<ND;++n)
					out << pnts(i)(n) << ' ';
				out << "\n";
			}

			/* SIDE INFO */
			for(i=0;i<nseg;++i)
				out << i << ": " << seg(i).pnt(0) << ' ' << seg(i).pnt(1) << endl;

			/* THEN TRI INFO */
			for(i=0;i<ntri;++i)
				out << i << ": " << tri(i).pnt(0) << ' ' << tri(i).pnt(1) << ' ' << tri(i).pnt(2) << endl;

			/* SIDE BOUNDARY INFO HEADER */
			out << "nebd: " << nebd << endl;
			for(i=0;i<nebd;++i) {
				out << "idnum: " << ebdry(i)->idnum << endl;
				out << "number: " << ebdry(i)->nseg << endl;
				for(int j=0;j<ebdry(i)->nseg;++j)
					out << j << ": " << ebdry(i)->seg(j) << std::endl;
			}

			/* VERTEX BOUNDARY INFO HEADER */
			out << "nvbd: " << nvbd << endl;
			for(i=0;i<nvbd;++i) {
				out << "idnum: " << vbdry(i)->idnum << endl;
				out << "point: " << vbdry(i)->pnt << endl;
			}

			out.close();

			break;

		case(binary):
			fnmapp = grd_nm +".bin";
			bout.open(fnmapp.c_str());
			if (bout.error()) {
				*gbl->log << "couldn't open output file" << fnmapp << "for output" << endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}

			/* HEADER LINES */
			bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::BigEndian)),1);
			bout.writeInt(static_cast<unsigned char>(bout.getFlag(binio::FloatIEEE)),1);
			bout.writeInt(npnt,sizeof(int));
			bout.writeInt(nseg,sizeof(int));
			bout.writeInt(ntri,sizeof(int));
			/* POINT INFO */
			for(i=0;i<npnt;++i) {
				for(n=0;n<ND;++n)
					bout.writeFloat(pnts(i)(n),binio::Double);
				bout.writeFloat(lngth(i),binio::Double);
			}

			/* SIDE INFO */
			for(i=0;i<nseg;++i) {
				bout.writeInt(seg(i).pnt(0),sizeof(int));
				bout.writeInt(seg(i).pnt(1),sizeof(int));
			}

			/* THEN TRI INFO */
			for(i=0;i<ntri;++i) {
				bout.writeInt(tri(i).pnt(0),sizeof(int));
				bout.writeInt(tri(i).pnt(1),sizeof(int));
				bout.writeInt(tri(i).pnt(2),sizeof(int));
			}

			/* SIDE BOUNDARY INFO HEADER */
			bout.writeInt(nebd,sizeof(int));
			for(i=0;i<nebd;++i) {
				bout.writeInt(ebdry(i)->idnum,sizeof(int));
				bout.writeInt(ebdry(i)->nseg,sizeof(int));
				for(int j=0;j<ebdry(i)->nseg;++j)
					bout.writeInt(ebdry(i)->seg(j),sizeof(int));
			}

			/* VERTEX BOUNDARY INFO HEADER */
			bout.writeInt(nvbd,sizeof(int));
			for(i=0;i<nvbd;++i) {
				bout.writeInt(vbdry(i)->idnum,sizeof(int));
				bout.writeInt(vbdry(i)->pnt,sizeof(int));
			}

			bout.close();

			break;

		case (datatank): {
#ifdef DATATANK
			DTMutableIntArray dt_tvrtx(3,ntri);
			DTMutableDoubleArray dt_vrtx(2,npnt);

			for (int tind=0;tind<ntri;++tind)
				for (int j=0;j<3;++j)
					dt_tvrtx(j,tind) = tri(tind).pnt(j);

			for (int pind=0;pind<npnt;++pind)
				for (int j=0;j<2;++j)
					dt_vrtx(j,pind) = pnts(pind)(j);

			DTTriangularGrid2D dt_grid(dt_tvrtx,dt_vrtx);

			std::string outputFilename(grd_nm +".dt");
			DTDataFile outputFile(outputFilename.c_str(),DTFile::NewReadWrite);
			// Output from computation
			Write(outputFile,"grid",dt_grid);
			outputFile.Save("TriangularGrid2D","Seq_grid");
#else
			*gbl->log << "Not supported on this platform\n";
#endif
			break;
		}

		case(vlength): {
			fnmapp = grd_nm +".lngth";
			out.open(fnmapp.c_str());
			if (!out) {
				*gbl->log << "couldn't open output file" << fnmapp << "for output" << endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}

			for(i=0;i<npnt;++i)
				out << lngth(i) << endl;

			break;
		}

		case(debug_adapt):
			/* CREATE EASYMESH OUTPUT FILES */
			fnmapp = grd_nm +".n";
			out.open(fnmapp.c_str());
			if (!out) {
				*gbl->log << "couldn't open output file " << fnmapp << "for output" << endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}

			out << npnt << endl;
			for(i=0;i<npnt;++i) {
				out << i << ": ";
				for(n=0;n<ND;++n)
					out << pnts(i)(n) << ' ';
				out << (tri(i).info&PDLTE != 0 ? -1 : 0) +(tri(i).info&PTOUC != 0 ? 1 : 0) << endl;
			}
			out.close();

			/* SIDE FILE */
			fnmapp = grd_nm +".s";
			out.open(fnmapp.c_str());
			if (!out) {
				*gbl->log << "couldn't open output file " << fnmapp << "for output" << endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			out << nseg << endl;
			for(i=0;i<nseg;++i) {
				out << i << ": " << seg(i).pnt(0) << ' ' << seg(i).pnt(1) << ' ';
				out << seg(i).tri(0) << ' ' << seg(i).tri(1) << ' ' << (tri(i).info&SDLTE ? -1 : (tri(i).info&STOUC ? 2 : 0)) << endl;
			}
			out.close();

			fnmapp = grd_nm +".e";
			out.open(fnmapp.c_str());
			if (!out) {
				*gbl->log << "couldn't open output file " << fnmapp << "for output" << endl;
				sim::abort(__LINE__,__FILE__,gbl->log);
			}
			out << ntri << endl;
			for(i=0;i<ntri;++i) {
				out << i << ": " << tri(i).pnt(0) << ' ' << tri(i).pnt(1) << ' ' << tri(i).pnt(2);
				out << ' ' << tri(i).tri(0) << ' ' << tri(i).tri(1) << ' ' << tri(i).tri(2);
				out << ' ' << tri(i).seg(0) << ' ' << tri(i).seg(1) << ' ' << tri(i).seg(2);
				out << " 0.0 0.0 " << (tri(i).info&TDLTE != 0 ? -1 : 0) +(tri(i).info&TTOUC != 0 ? 1 : 0) << endl;
			}
			out.close();

			break;
		case boundary: {
			fnmapp = grd_nm +"_bdry.inpt";
			out.open(fnmapp.c_str());
			for(i=0;i<nvbd;++i) vbdry(i)->output(out);
			for(i=0;i<nebd;++i) ebdry(i)->output(out);
			out.close();
			break;
		}

		default:
			*gbl->log << "That filetype is not supported yet" << endl;
			break;
	}

	return;
}

void tri_mesh::setinfo() {
	int i,j;

	/* SET UP POINT BC INFORMATION FOR OUTPUT */
	for(i=0;i<npnt;++i)
		pnt(i).info = 0;

	for(i=0;i<nvbd;++i)
		pnt(vbdry(i)->pnt).info = vbdry(i)->idnum;

	/* SET UP SIDE BC INFORMATION FOR EASYMESH OUTPUT */
	for(i=0;i<nseg;++i)
		seg(i).info = 0;

	for(i=0;i<nebd;++i)
		for(j=0;j<ebdry(i)->nseg;++j)
			seg(ebdry(i)->seg(j)).info = ebdry(i)->idnum;

	/* SET UP TRI INFO FOR EASYMESH OUTPUT */
	for(i=0;i<ntri;++i)
		tri(i).info = 0;

	return;
}
