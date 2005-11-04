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

 void tri_hp::output(char *name, filetype typ, int tlvl) {
   ofstream out;
   char fnmapp[100];
   int i,j,k,m,n,v0,v1,sind,tind,indx,sgn;
   int ijind[MXTM][MXTM];
   
   out.setf(std::ios::scientific, std::ios::floatfield);
   out.precision(10);
         
   switch (typ) {
      case (text):
         strcpy(fnmapp,name);
         strcat(fnmapp,".txt");
         out.open(fnmapp);
         if (!out) {
            *sim::log << "couldn't open text output file " << fnmapp;
            exit(1);
         }
         
         /* HEADER INFORMATION */
         out << "p0 = " << p0 << std::endl;
         out << "nvrtx = " << nvrtx << ", nside = " << nside << ", ntri = " << ntri << std::endl;
         out << "END OF HEADER" << std::endl;

         for(i=0;i<nvrtx;++i) {
            for(n=0;n<NV;++n)
               out << ugbd(tlvl).v(i,n) << '\t';
            out << std::endl;
         }
         
         for(i=0;i<nside;++i) {
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

         /* BOUNDARY INFO */
         for(i=0;i<nsbd;++i) {
            hp_sbdry(i)->output(out,typ,tlvl);
         }
         
         out.close();
         break;
      
      case(tecplot):
         strcpy(fnmapp,name);
         strcat(fnmapp,".dat");
         out.open(fnmapp);
         if (!out) {
            *sim::log << "couldn't open tecplot output file " << fnmapp;
            exit(1);
         }
      
         out << "ZONE F=FEPOINT, ET=TRIANGLE, N = " << nvrtx+basis::tri(log2p).sm*nside+basis::tri(log2p).im*ntri << ", E = " << ntri*(basis::tri(log2p).sm+1)*(basis::tri(log2p).sm+1) << std::endl;

         /* VERTEX MODES */
         for(i=0;i<nvrtx;++i) {
            for(n=0;n<ND;++n)
               out << vrtxbd(tlvl)(i)(n) << ' ';
            for(n=0;n<NV;++n)
               out << ugbd(tlvl).v(i,n) << ' ';               
            out << std::endl;
         }
         
         if (basis::tri(log2p).p > 1) {
            /* SIDE MODES */
            for(sind=0;sind<nside;++sind) {
               if (sd(sind).info < 0) {
                  v0 = sd(sind).vrtx(0);
                  v1 = sd(sind).vrtx(1);
                  for(n=0;n<ND;++n)
                     basis::tri(log2p).proj1d_leg(vrtxbd(tlvl)(v0)(n),vrtxbd(tlvl)(v1)(n),&crd(n)(0,0));
               }
               else {
                  crdtocht1d(sind,tlvl);
   
                  for(n=0;n<ND;++n)
                  	basis::tri(log2p).proj1d_leg(&cht(n,0),&crd(n)(0,0));
               }
               ugtouht1d(sind,tlvl);
               for(n=0;n<NV;++n)
                  basis::tri(log2p).proj1d_leg(&uht(n)(0),&u(n)(0,0));

               for(i=1;i<basis::tri(log2p).sm+1;++i) {
                  for(n=0;n<ND;++n)
                     out << crd(n)(0,i) << ' ';
                  for(n=0;n<NV;++n)
                     out << u(n)(0,i) << ' ';               
                  out << std::endl;
               }
            }
   
            /* INTERIOR MODES */
            if (basis::tri(log2p).p > 2) {
               for(tind = 0; tind < ntri; ++tind) {
                  ugtouht(tind,tlvl);
                  for(n=0;n<NV;++n)
                     basis::tri(log2p).proj_leg(&uht(n)(0),&u(n)(0,0),MXGP);
                     
                  if (td(tind).info < 0) {
                     for(n=0;n<ND;++n)
                        basis::tri(log2p).proj_leg(vrtxbd(tlvl)(td(tind).vrtx(0))(n),vrtxbd(tlvl)(td(tind).vrtx(1))(n),vrtxbd(tlvl)(td(tind).vrtx(2))(n),&crd(n)(0,0),MXGP);
                  }
                  else {
                     crdtocht(tind,tlvl);
                     for(n=0;n<ND;++n)
                        basis::tri(log2p).proj_bdry_leg(&cht(n,0),&crd(n)(0,0),MXGP);
                  }
                  
                  for(i=1;i<basis::tri(log2p).sm;++i) {
                     for(j=1;j<basis::tri(log2p).sm-(i-1);++j) {
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
            ijind[0][basis::tri(log2p).sm+1] = td(tind).vrtx(0);
            ijind[0][0] = td(tind).vrtx(1);
            ijind[basis::tri(log2p).sm+1][0] = td(tind).vrtx(2);

            /* SIDES */
            indx = td(tind).side(0);
            sgn = td(tind).sign(0);
            if (sgn < 0) {
               for(i=0;i<basis::tri(log2p).sm;++i)
                  ijind[i+1][0] = nvrtx +(indx+1)*basis::tri(log2p).sm -(i+1);
            }
            else {
               for(i=0;i<basis::tri(log2p).sm;++i)
                  ijind[i+1][0] = nvrtx +indx*basis::tri(log2p).sm +i;
            }
            
            indx = td(tind).side(1);
            sgn = td(tind).sign(1);
            if (sgn > 0) {
               for(i=0;i<basis::tri(log2p).sm;++i)
                  ijind[basis::tri(log2p).sm-i][i+1] = nvrtx +indx*basis::tri(log2p).sm +i;
            }
            else {
               for(i=0;i<basis::tri(log2p).sm;++i)
                  ijind[basis::tri(log2p).sm-i][i+1] = nvrtx +(indx+1)*basis::tri(log2p).sm -(i+1);
            }

            indx = td(tind).side(2);
            sgn = td(tind).sign(2);
            if (sgn > 0) {
               for(i=0;i<basis::tri(log2p).sm;++i)
                  ijind[0][i+1] = nvrtx +(indx+1)*basis::tri(log2p).sm -(i+1);
            }
            else {
               for(i=0;i<basis::tri(log2p).sm;++i)
                  ijind[0][i+1] = nvrtx +indx*basis::tri(log2p).sm +i;
            }
   
            /* INTERIOR VERTICES */
            k = 0;
            for(i=1;i<basis::tri(log2p).sm;++i) {
               for(j=1;j<basis::tri(log2p).sm-(i-1);++j) {
                  ijind[i][j] = nvrtx +nside*basis::tri(log2p).sm +tind*basis::tri(log2p).im +k;
                  ++k;
               }
            }
   
            /* OUTPUT CONNECTION LIST */      
            for(i=0;i<basis::tri(log2p).sm+1;++i) {
               for(j=0;j<basis::tri(log2p).sm-i;++j) {
                  out << ijind[i][j]+1 << ' ' << ijind[i+1][j]+1 << ' ' << ijind[i][j+1]+1 << std::endl;
                  out << ijind[i+1][j]+1 << ' ' << ijind[i+1][j+1]+1 << ' ' << ijind[i][j+1]+1 << std::endl;
               }
               out << ijind[i][basis::tri(log2p).sm-i]+1 << ' ' << ijind[i+1][basis::tri(log2p).sm-i]+1 << ' ' << ijind[i][basis::tri(log2p).sm+1-i]+1 << std::endl;
            }
         }
         out.close();
         break; 
               
      case(adapt_diagnostic): {         
         strcpy(fnmapp,name);
         strcat(fnmapp,".dat");
         out.open(fnmapp);
         if (!out) {
            *sim::log << "couldn't open tecplot output file " << fnmapp;
            exit(1);
         }
      
         out << "ZONE F=FEPOINT, ET=TRIANGLE, N = " << nvrtx << ", E = " << ntri << std::endl;
         for(i=0;i<nvrtx;++i) {
            for(n=0;n<ND;++n)
               out << vrtx(i)(n) << ' ';
            out << vlngth(i) << ' ' << log10(fscr1(i)) << std::endl;               
         }

         /* OUTPUT CONNECTIVY INFO */
         out << std::endl << "#CONNECTION DATA#" << std::endl;
         
         for(tind=0;tind<ntri;++tind)
            out << td(tind).vrtx(0)+1 << ' ' << td(tind).vrtx(1)+1 << ' ' << td(tind).vrtx(2)+1 << std::endl;
         
         break;
      }
      
      default:
         *sim::log << "can't output a tri_hp to that filetype" << std::endl;
         exit(1);
         break;
     }
   
    return;
}

 void tri_hp::input(char *name, filetype typ, int tlvl) {
   int i,j,k,m,n,pin,pmin,indx;
   int bnum;
   char buffer[200];
   char trans[] = "T";
   ifstream in;
   FLT fltskip;

   switch(typ) {
   
      case (text):
         strcpy(buffer,name);
         strcat(buffer,".txt");
         in.open(buffer);
         if (!in ) {
            printf("couldn't open text input file %s\n",name);
            exit(1);
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

         for(i=0;i<nvrtx;++i) {
            for(n=0;n<NV;++n)
               in >> ugbd(tlvl).v(i,n);
         }
         
         for(i=0;i<nside;++i) {
            for(m=0;m<(pmin-1);++m) {
               for(n=0;n<NV;++n)
                  in >> ugbd(tlvl).s(i,m,n);
            }
            indx += p0 -pmin;

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
         

         /* BOUNDARY INFO */
         for(i=0;i<nsbd;++i)
            hp_sbdry(i)->input(in,typ,tlvl);

         in.close();
         break;

      case(tecplot):
         /* CAN ONLY DO THIS IF HAVE MESH FILE */
         strcpy(buffer,name);
         strcat(buffer,".dat");
         in.open(buffer);
         if (!in) {
            *sim::log << "couldn't open tecplot input file " << name << std::endl;
            exit(1);
         }
         in.ignore(80,'\n');

         for(i=0;i<nvrtx;++i) {
            in >> fltskip >> fltskip;
            for(n=0;n<NV;++n)
               in >> ugbd(tlvl).v(i,n);
         }
         
         if (basis::tri(log2p).sm < 1) return;
            
         FLT temp,matrix[MXTM][MXTM];
         int sind,info,ipiv[2*MXTM];
            
         /* REVERSE OUTPUTING PROCESS */
         for(m=0;m<basis::tri(log2p).sm;++m)
            for(n=0;n<basis::tri(log2p).sm;++n)
               matrix[n][m] = basis::tri(log2p).lgrnge1d(m+2,n+1);
            
         GETRF(basis::tri(log2p).sm,basis::tri(log2p).sm,matrix[0],MXTM,ipiv,info);
         if (info != 0) {
            printf("DGETRF FAILED FOR INPUTING TECPLOT SIDES\n");
            exit(1);
         }
            
         for(sind=0;sind<nside;++sind) {
            if (sd(sind).info < 0) {
               ugtouht1d(sind,tlvl);
               for(n=0;n<NV;++n)
                  basis::tri(log2p).proj1d_leg(&uht(n)(0),&u(n)(0,0));
               
               for(m=0;m<basis::tri(log2p).sm;++m) {
                  in >> fltskip >> fltskip;
                  for(n=0;n<NV;++n) {
                     in >> temp;
                     u(n)(0,m+1) -= temp;
                  }
               }
               for(n=0;n<NV;++n) {
                  GETRS(trans,basis::tri(log2p).sm,1,matrix[0],MXTM,ipiv,&u(n)(0,1),MXTM,info);
                  for(m=0;m<basis::tri(log2p).sm;++m)
                     ugbd(tlvl).s(sind,m,n) = -u(n)(0,1+m);
               }
            }
            else {
               crdtocht1d(sind,tlvl);
               for(n=0;n<ND;++n)
                  basis::tri(log2p).proj1d_leg(&cht(n,0),&crd(n)(0,0));
               
               ugtouht1d(sind,tlvl);
               for(n=0;n<NV;++n)
                  basis::tri(log2p).proj1d_leg(&uht(n)(0),&u(n)(0,0));
               
               for(m=0;m<basis::tri(log2p).sm;++m) {
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
                  GETRS(trans,basis::tri(log2p).sm,1,matrix[0],MXTM,ipiv,&u(n)(0,1),MXTM,info);
                  for(m=0;m<basis::tri(log2p).sm;++m)
                     ugbd(tlvl).s(sind,m,n) = -u(n)(0,1+m);
               }
               
               bnum = getbdrynum(sd(sind).tri(1));
               indx = getbdryel(sd(sind).tri(1));
               for(n=0;n<ND;++n) {
                  GETRS(trans,basis::tri(log2p).sm,1,matrix[0],MXTM,ipiv,&crd(n)(0,1),MXTM,info);
                  for(m=0;m<basis::tri(log2p).sm;++m)
                     hp_sbdry(bnum)->crdsbdin(indx,m,n,tlvl) = -crd(n)(0,1+m);
               }
            }
         }
            
         if (basis::tri(log2p).im < 1) return;
         
         int tind;
         
         ugbd(tlvl).i(Range(0,ntri),Range::all(),Range::all()) = 0.0;
         
         /* REVERSE OUTPUTING PROCESS */
         for(m=0;m<basis::tri(log2p).im;++m) {
            n = 0;
            for(i=1;i<basis::tri(log2p).sm;++i) {
               for(j=1;j<basis::tri(log2p).sm-(i-1);++j) {
                  matrix[n++][m] = basis::tri(log2p).lgrnge(m+basis::tri(log2p).bm,i,j);
               }
            }
         }
    
         GETRF(basis::tri(log2p).im,basis::tri(log2p).im,matrix[0],MXTM,ipiv,info);
         if (info != 0) {
            printf("DGETRF FAILED FOR INPUTING TECPLOT SIDES\n");
            exit(1);
         }

         for(tind=0;tind<ntri;++tind) {
            ugtouht(tind,tlvl);
            for(n=0;n<NV;++n)
               basis::tri(log2p).proj_leg(&uht(n)(0),&u(n)(0,0),MXGP);
            
            m = 0;
            for(i=1;i<basis::tri(log2p).sm;++i) {
               for(j=1;j<basis::tri(log2p).sm-(i-1);++j) {
                  for(n=0;n<ND;++n)
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
               GETRS(trans,basis::tri(log2p).im,1,matrix[0],MXTM,ipiv,&uht(n)(0),MXTM,info);
               for(m=0;m<basis::tri(log2p).im;++m)
                  ugbd(tlvl).i(tind,m,n) = -uht(n)(m);
            }
         }            
               
         break;
               
      default:
         *sim::log << "can't input a tri_hp from that filetype" << std::endl;
         exit(1);
         break;
   }
   
   return;
}