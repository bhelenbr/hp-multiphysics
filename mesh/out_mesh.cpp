#include <string.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include "mesh.h"
#include "boundary.h"

using namespace std;

int mesh::output(const char *filename, ftype::name filetype) const {
   char fnmapp[100];
   int i,j,n,tind,count;
   ofstream out;
   
   out.setf(std::ios::scientific, std::ios::floatfield);
   out.precision(10);
         
   switch (filetype) {
   
      case (ftype::easymesh):
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
               out << vrtx[i][n] << ' ';
            out << vd[i].info << endl;
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
            out << i << ": " << sd[i].vrtx[0] << ' ' << sd[i].vrtx[1] << ' ';
            out << sd[i].tri[0] << ' ' << sd[i].tri[1] << ' ' << sd[i].info <<endl;
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
            out << i << ": " << td[i].vrtx[0] << ' ' << td[i].vrtx[1] << ' ' << td[i].vrtx[2];
            out << ' ' << td[i].tri[0] << ' ' << td[i].tri[1] << ' ' << td[i].tri[2];
            out << ' ' << td[i].side[0] << ' ' << td[i].side[1] << ' ' << td[i].side[2];
            out << " 0.0 0.0 " << td[i].info << endl;
         }   
         out.close();
         break;
      
      case (ftype::tecplot):
         strcpy(fnmapp,filename);
         strcat(fnmapp,".dat");
         out.open(fnmapp);
         if (out == NULL ) {
            *log << "couldn't open output file " << fnmapp << "for output" << endl;
            exit(1);
         }

         out << "ZONE F=FEPOINT, ET=TRIANGLE, N = " << nvrtx << ", E = " << ntri << endl;
         
         for(i=0;i<nvrtx;++i) {
            for(n=0;n<ND;++n)
               out << vrtx[i][n] << ' ';
            out << endl;
         }

         out << "\n#CONNECTION DATA#\n";
         
         for(i=0;i<ntri;++i)
            out << td[i].vrtx[0]+1 << ' ' << td[i].vrtx[1]+1 << ' ' << td[i].vrtx[2]+1 << endl;
            
         out.close();
         break;

      case (ftype::mavriplis):
         strcpy(fnmapp,filename);
         strcat(fnmapp,".GDS");
         out.open(fnmapp);
         if (out == NULL ) {
            *log << "couldn't open output file " << fnmapp << "for output" << endl;
            exit(1);
         }

         out << "YIN " << nvrtx << endl;

         for(i=0;i<nvrtx;++i) {
            for(n=0;n<ND;++n)
               out << "       " << vrtx[i][n];
            if (ND == 2) out << "       " << 0.0;
            out << endl;
         }

         out << "0\n";
         out << ntri << '\n';

         for(i=0;i<ntri;++i) {
            out << "1 Default" << endl;
            out << 3 << ' ' << td[i].vrtx[0]+1 << ' ' << td[i].vrtx[1]+1 << ' ' << td[i].vrtx[2]+1 << endl;
         }
            out.close();
         break; 
         
      case (ftype::gambit):
         strcpy(fnmapp,filename);
         strcat(fnmapp,".FDNEUT");
         out.open(fnmapp);
         if (!out) {
            *log << "couldn't open output file " << fnmapp << "for output" << endl;
            exit(1);
         }

         count = ntri;
         for(i=0;i<nsbd;++i)
            count += sbdry[i]->nel;
   
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
               out << setw(20) << vrtx[i][n];
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
            out << setw(8) << td[tind].vrtx[0]+1;
            out << setw(8) << td[tind].vrtx[1]+1;
            out << setw(8) << td[tind].vrtx[2]+1;
            out << endl;
         }
        
         count = ntri+1;
        
         for(i=0;i<nsbd;++i) {
            out << "GROUP:" << setw(9) << i+2;
            out << " ELEMENTS:" << setw(10) << sbdry[i]->nel;
            out << " NODES:" << setw(13) << 2;
            out << " GEOMETRY:" << setw(5) << 0;
            out << " ID:" << setw(4) << sbdry[i]->idnum << endl;

            out << "ENTITY NAME:   surface.1\n";
            for(j=0;j<sbdry[i]->nel;++j) {
               out << setw(8) << count++;
               out << setw(8) << sd[sbdry[i]->el[j]].vrtx[0]+1;
               out << setw(8) << sd[sbdry[i]->el[j]].vrtx[1]+1 << endl;
            }
         }
         out.close();
         break;

      case(ftype::text):
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
               out << vrtx[i][n] << ' '; 
            out << "\n";
         }
            
         out.close();
         break;
         
      case(ftype::grid):
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
            out << i << ": ";
            for(n=0;n<ND;++n)
               out << vrtx[i][n] << ' ';
            out << "\n";
         }
                    
         /* SIDE INFO */
         for(i=0;i<nside;++i)
            out << i << ": " << sd[i].vrtx[0] << ' ' << sd[i].vrtx[1] << endl;

         /* THEN TRI INFO */
         for(i=0;i<ntri;++i)
            out << i << ": " << td[i].vrtx[0] << ' ' << td[i].vrtx[1] << ' ' << td[i].vrtx[2] << endl;

         /* SIDE BOUNDARY INFO HEADER */
         out << "nsbd: " << nsbd << endl;
         for(i=0;i<nsbd;++i) {
            out << "idnum: " << sbdry[i]->idnum << endl;
         	out << "number: " << sbdry[i]->nel << endl;
            for(int j=0;j<sbdry[i]->nel;++j)
               out << j << ": " << sbdry[i]->el[j] << std::endl;
         }

         /* VERTEX BOUNDARY INFO HEADER */
         out << "nvbd: " << nvbd << endl;
         for(i=0;i<nvbd;++i) {
            out << "idnum: " << vbdry[i]->idnum << endl;
            out << "point: " << vbdry[i]->v0 << endl;
         }
         
         out.close();
         
         break;

      default:
         *log << "That filetype is not supported yet" << endl;
         break;
   }
   
   return(1);
}

void mesh::bdry_output(const char *filename) const {
   char fnmapp[120];
   ofstream out;
   int i;
   
   strcpy(fnmapp,filename);
   strcat(fnmapp,"_bdry.inpt");
   out.open(fnmapp);
   for(i=0;i<nvbd;++i) vbdry[i]->output(out);
   for(i=0;i<nsbd;++i) sbdry[i]->output(out);
   out.close();
}

void mesh::setbcinfo() {
   int i,j;
   
   /* SET UP VRTX BC INFORMATION FOR OUTPUT */
   for(i=0;i<nvrtx;++i)
      vd[i].info = 0;
   
   for(i=0;i<nvbd;++i)
      vd[vbdry[i]->v0].info = vbdry[i]->idnum;

   /* SET UP SIDE BC INFORMATION FOR EASYMESH OUTPUT */
   for(i=0;i<nside;++i)
      sd[i].info = 0;
   
   for(i=0;i<nsbd;++i)
      for(j=0;j<sbdry[i]->nel;++j)
         sd[sbdry[i]->el[j]].info = sbdry[i]->idnum;

   /* SET UP TRI INFO FOR EASYMESH OUTPUT */         
   for(i=0;i<ntri;++i)
      td[i].info = 0;
         
   return;
}
