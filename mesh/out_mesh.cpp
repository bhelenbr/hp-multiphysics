#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include "mesh.h"

using namespace std;

template int mesh<2>::out_mesh(FLT (*vin)[2], const char *filename, FTYPE filetype) const;
template int mesh<3>::out_mesh(FLT (*vin)[3], const char *filename, FTYPE filetype) const;

template<int ND> int mesh<ND>::out_mesh(FLT (*vin)[ND], const char *filename, FTYPE filetype) const {
   char fnmapp[100];
   int i,j,n,tind,count;
   ofstream out;
   
   out.setf(std::ios::scientific, std::ios::floatfield);
   out.precision(10);
         
   switch (filetype) {
   
      case (easymesh):
         /* CREATE EASYMESH OUTPUT FILES */
         strcpy(fnmapp,filename);
         strcat(fnmapp,".n");
         out.open(fnmapp);
         if (!out) {
            *log << "couldn't open output file " << fnmapp << "for output" << endl;
            exit(1);
         }
         
         out << nvrtx << endl;
         for(i=0;i<nvrtx;++i) {
            out << i << ": ";
            for(n=0;n<ND;++n)
               out << vin[i][n] << ' ';
            out << vinfo[i] << endl;
         }            
         out.close();

         /* SIDE FILE */      
         strcpy(fnmapp,filename);
         strcat(fnmapp,".s");
         out.open(fnmapp);
         if (!out) {
            *log << "couldn't open output file " << fnmapp << "for output" << endl;
            exit(1);
         }
         out << nside << endl;   
         for(i=0;i<nside;++i) {
            out << i << ": " << svrtx[i][0] << ' ' << svrtx[i][1] << ' ';
            out << stri[i][0] << ' ' << stri[i][1] << ' ' << sinfo[i] <<endl;
         }
         out.close();
   
         strcpy(fnmapp,filename);
         strcat(fnmapp,".e");
         out.open(fnmapp);
         if (!out) {
            *log << "couldn't open output file " << fnmapp << "for output" << endl;
            exit(1);
         }
         out << ntri << endl;
         for(i=0;i<ntri;++i) {
            out << i << ": " << tvrtx[i][0] << ' ' << tvrtx[i][1] << ' ' << tvrtx[i][2];
            out << ' ' << ttri[i][0] << ' ' << ttri[i][1] << ' ' << ttri[i][2];
            out << ' ' << tside[i].side[0] << ' ' << tside[i].side[1] << ' ' << tside[i].side[2];
            out << " 0.0 0.0 " << tinfo[i] << endl;
         }   
         out.close();
         break;
      
      case (tecplot):
         strcpy(fnmapp,filename);
         strcat(fnmapp,".dat");
         out.open(fnmapp);
         if (out == NULL ) {
            *log << "couldn't open output file " << fnmapp << "for output" << endl;
            exit(1);
         }

         out << "ZONE F=FEPOINT, ET=TRIANGLE, N = " << nvrtx << " E = " << ntri << endl;
         
         for(i=0;i<nvrtx;++i) {
            for(n=0;n<ND;++n)
               out << vin[i][n] << ' ';
            out << endl;
         }

         out << "\n#CONNECTION DATA#\n";
         
         for(i=0;i<ntri;++i)
            out << tvrtx[i][0]+1 << ' ' << tvrtx[i][1]+1 << ' ' << tvrtx[i][2]+1 << endl;
            
         out.close();
         break;

      case (gambit):
         strcpy(fnmapp,filename);
         strcat(fnmapp,".FDNEUT");
         out.open(fnmapp);
         if (!out) {
            *log << "couldn't open output file " << fnmapp << "for output" << endl;
            exit(1);
         }

         count = ntri;
         for(i=0;i<nsbd;++i)
            count += sbdry[i]->nsd();
   
         out << "** FIDAP NEUTRAL FILE\n";
         out << filename << "\n";
         out << "VERSION    8.01\n";
         out << "29 Nov 1999    13:23:58\n";
         out << "   NO. OF NODES   NO. ELEMENTS NO. ELT GROUPS          NDFCD          NDFVL\n";
         
         out << setw(15) << nvrtx << setw(15) << count << setw(15) << nsbd+1 << setw(15) << ND << setw(15) << ND << endl;
         out << "   STEADY/TRANS     TURB. FLAG FREE SURF FLAG    COMPR. FLAG   RESULTS ONLY\n";
         out << setw(12) << 0 << setw(12) << 0 << setw(12) << 0 << setw(12) << 0 << setw(12) << 0 << endl;
         out << "TEMPERATURE/SPECIES FLAGS\n";
         out << " 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0\n";
         out << "PRESSURE FLAGS - IDCTS, IPENY MPDF\n";
         out << "         1         1         0\n";
         out << "NODAL COORDINATES\n";
         
         for(i=0;i<nvrtx;++i) {
            out << setw(10) << i+1;
            for(n=0;n<ND;++n)
               out << setw(20) << vin[i][n];
            out << '\n';
         }
            
         out << "BOUNDARY CONDITIONS\n";
         out << "         0         0         0     0.0\n";
         out << "ELEMENT GROUPS\n";
         out << "GROUP:" << setw(9) << 1;
         out << " ELEMENTS:" << setw(10) << ntri;
         out << " NODES:" << setw(13) << 3;
         out << " GEOMETRY:" << setw(5) << 2;
         out << " TYPE:" << setw(4) << 2 << endl;
         out << "ENTITY NAME:   fluid\n";
         
         for(tind=0;tind<ntri;++tind) {
            out << setw(8) << tind+1;
            out << setw(8) << tvrtx[tind][0]+1;
            out << setw(8) << tvrtx[tind][1]+1;
            out << setw(8) << tvrtx[tind][2]+1;
            out << endl;
         }
        
         count = ntri+1;
        
         for(i=0;i<nsbd;++i) {
            out << "GROUP:" << setw(9) << i+2;
            out << " ELEMENTS:" << setw(10) << sbdry[i]->nsd();
            out << " NODES:" << setw(13) << 2;
            out << " GEOMETRY:" << setw(5) << 0;
            out << " TYPE:" << setw(4) << sbdry[i]->idnty() << endl;
            out << "ENTITY NAME:   surface.1\n";
            for(j=0;j<sbdry[i]->nsd();++j) {
               out << setw(8) << count++;
               out << setw(8) << svrtx[sbdry[i]->sd(j)][0]+1;
               out << setw(8) << svrtx[sbdry[i]->sd(j)][1]+1 << endl;
            }
         }
         out.close();
         break;

      case(text):
         /* JUST OUTPUT VERTEX POSITIONS FOR DEFORMING MESH */
         strcpy(fnmapp,filename);
         strcat(fnmapp,".txt");
         out.open(fnmapp);
         if (!out) {
            *log << "couldn't open output file" << fnmapp << "for output" << endl;
            exit(1);
         }
         out << nvrtx << endl;
         for(i=0;i<nvrtx;++i) {
            out << i << ":";
            for(n=0;n<ND;++n)
               out << vin[i][n] << ' '; 
            out << "\n";
         }
            
         out.close();
         break;
         
      case(grid):
         strcpy(fnmapp,filename);
         strcat(fnmapp,".grd");
         out.open(fnmapp);
         if (!out) {
            *log << "couldn't open output file" << fnmapp << "for output" << endl;
            exit(1);
         }
         
         /* HEADER LINES */
         out << "nvrtx: " << nvrtx;
         out << " nside: " << nside;
         out << " ntri: " << ntri << endl;

         /* VRTX INFO */                        
         for(i=0;i<nvrtx;++i) {
            for(n=0;n<ND;++n)
               out << vin[i][n] << ' ';
            out << "\n";
         }
                    
         /* SIDE INFO */
         for(i=0;i<nside;++i)
            out << svrtx[i][0] << ' ' << svrtx[i][1] << endl;

         /* THEN TRI INFO */
         for(i=0;i<ntri;++i)
            out << tvrtx[i][0] << ' ' << tvrtx[i][1] << ' ' << tvrtx[i][2] << endl;

         /* SIDE BOUNDARY INFO HEADER */
         out << "nsbd: " << nsbd << endl;
         for(i=0;i<nsbd;++i)
            sbdry[i]->output(out);

         /* VERTEX BOUNDARY INFO HEADER */
         out << "nvbd: " << nvbd << endl;
         for(i=0;i<nvbd;++i)
            vbdry[i]->output(out);
         
         out.close();
         
         break;

      default:
         *log << "That filetype is not supported yet" << endl;
         break;
   }

    return(1);
}

template void mesh<2>::setbcinfo();
template void mesh<3>::setbcinfo();

template<int ND> void mesh<ND>::setbcinfo() {
   int i,j;
   
   /* SET UP VRTX BC INFORMATION FOR OUTPUT */
   for(i=0;i<nvrtx;++i)
      vinfo[i] = 0;
   
   for(i=0;i<nvbd;++i)
      vinfo[vbdry[i]->v()] = vbdry[i]->idnty();

   /* SET UP SIDE BC INFORMATION FOR EASYMESH OUTPUT */
   for(i=0;i<nside;++i)
      sinfo[i] = 0;
   
   for(i=0;i<nsbd;++i)
      for(j=0;j<sbdry[i]->nsd();++j)
         sinfo[sbdry[i]->sd(j)] = sbdry[i]->idnty();

   /* SET UP TRI INFO FOR EASYMESH OUTPUT */         
   for(i=0;i<ntri;++i)
      tinfo[i] = 0;
         
   return;
}
