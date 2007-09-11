/*
 *  l2error.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on Tue Jun 11 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_hp.h"
#include <utilities.h>
#include <boundary.h>
#include "hp_boundary.h"

void tri_hp::l2error(init_bdry_cndtn *comparison) {
	int i,j,n,tind,loc[NV];
	FLT err,mxr[NV],l2r[NV];
    TinyVector<FLT,2> pt;
    
#ifdef TAYLOR
    extern FLT ppipi;
//    
//    ptprobe(0.5,0.5,l2r);
//    ppipi = l2r[2];
    
/* MATCH PRESSURE AT ONE POINT */
    ppipi = 0.0;
    ppipi = -comparsion->f(2,vrtx(0)(0),vrtx(0)(1))+ug.v(0,2);
#endif
    
	for(n=0;n<NV;++n) {
		mxr[n] = 0.0;
		l2r[n] = 0.0;
	}
	
	for(tind=0;tind<ntri;++tind) {
        
        if (td(tind).info > -1) {
            crdtocht(tind);
            for(n=0;n<ND;++n)
                basis::tri(log2p).proj_bdry(&cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), &dcrd(n,1)(0,0),MXGP);
        }
        else {
            for(n=0;n<ND;++n)
                basis::tri(log2p).proj(vrtx(td(tind).vrtx(0))(n),vrtx(td(tind).vrtx(1))(n),vrtx(td(tind).vrtx(2))(n),&crd(n)(0,0),MXGP);

            for(i=0;i<basis::tri(log2p).gpx;++i) {
                for(j=0;j<basis::tri(log2p).gpn;++j) {
                    for(n=0;n<ND;++n) {
                        dcrd(n,0)(i,j) = 0.5*(vrtx(td(tind).vrtx(1))(n) -vrtx(td(tind).vrtx(0))(n));
                        dcrd(n,1)(i,j) = 0.5*(vrtx(td(tind).vrtx(2))(n) -vrtx(td(tind).vrtx(0))(n));
                    }
                }
            }
        }

        ugtouht(tind);
		for(n=0;n<NV;++n)
			basis::tri(log2p).proj(&uht(n)(0),&u(n)(0,0),MXGP);
		
 		for (i=0;i<basis::tri(log2p).gpx;++i) {	
			for (j=0;j<basis::tri(log2p).gpn;++j) {
                cjcb(i,j) = (dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j));
                pt(0) = crd(0)(i,j);
                pt(1) = crd(1)(i,j);
                for(n=0;n<NV;++n) {
                    err =  fabs(u(n)(i,j)-comparison->f(n,pt));
                    if (err > mxr[n]) {
                        mxr[n] = err;
                        loc[n] = tind;
                    }
                    l2r[n] += err*err*basis::tri(log2p).wtx(i)*basis::tri(log2p).wtn(j)*cjcb(i,j);
                }
            }
        }	
	}
    
	for(n=0;n<NV;++n) {
		l2r[n] = sqrt(l2r[n]); 
		printf("#L_2: %.3e L_inf %.3e %4d ",l2r[n],mxr[n],loc[n]);
	}
	printf("\n");
		
	return;
}

/* CALCULATE AREA/CIRCUMFERENCE/YBAR */
void tri_hp::integrated_averages(Array<FLT,1> a) {
    int i,j,n,tind;

    /* a(0) = area */
    /* a(1) = xbar */
    /* a(2) = ybar */
    /* a(3-...) variable averages */
    a = 0.0;
    
    for(tind=0;tind<ntri;++tind) {
        if (td(tind).info > -1) {
            crdtocht(tind);
            for(n=0;n<mesh::ND;++n)
                basis::tri(log2p).proj_bdry(&cht(n,0), &crd(n)(0,0), &dcrd(n,0)(0,0), &dcrd(n,1)(0,0),MXGP);
        }
        else {
            for(n=0;n<mesh::ND;++n)
                basis::tri(log2p).proj(vrtx(td(tind).vrtx(0))(n),vrtx(td(tind).vrtx(1))(n),vrtx(td(tind).vrtx(2))(n),&crd(n)(0,0),MXGP);

            for(i=0;i<basis::tri(log2p).gpx;++i) {
                for(j=0;j<basis::tri(log2p).gpn;++j) {
                    for(n=0;n<mesh::ND;++n) {
                        dcrd(n,0)(i,j) = 0.5*(vrtx(td(tind).vrtx(2))(n) -vrtx(td(tind).vrtx(1))(n));
                        dcrd(n,1)(i,j) = 0.5*(vrtx(td(tind).vrtx(0))(n) -vrtx(td(tind).vrtx(1))(n));
                    }
                }
            }
        }
        
        ugtouht(tind);
        for(n=0;n<NV;++n)
            basis::tri(log2p).proj(&uht(n)(0),&u(n)(0,0),MXGP);

        for(i=0;i<basis::tri(log2p).gpx;++i) {
            for(j=0;j<basis::tri(log2p).gpn;++j) {
                cjcb(i,j) = dcrd(0,0)(i,j)*dcrd(1,1)(i,j) -dcrd(1,0)(i,j)*dcrd(0,1)(i,j);
                a(0) += RAD(crd(0)(i,j))*basis::tri(log2p).wtx(i)*basis::tri(log2p).wtn(j)*cjcb(i,j);
                a(1) += crd(0)(i,j)*RAD(crd(0)(i,j))*basis::tri(log2p).wtx(i)*basis::tri(log2p).wtn(j)*cjcb(i,j);
                a(2) += crd(1)(i,j)*RAD(crd(0)(i,j))*basis::tri(log2p).wtx(i)*basis::tri(log2p).wtn(j)*cjcb(i,j);
                for(n=0;n<NV;++n) 
                    a(3+n) += u(n)(i,j)*RAD(crd(0)(i,j))*basis::tri(log2p).wtx(i)*basis::tri(log2p).wtn(j)*cjcb(i,j);
            }
        }
    }
    a(Range(1,2+NV)) /= a(0);
    
    return;
}

/* UTILITY ROUTINE TO OUTPUT TRUNCATION ERROR ASSUMING STORED IN FSCR1 */
void tri_hp::output_error() {
   int i,n,tind;
   ostringstream fname;
   ofstream out;
    
   fname << idprefix << +"_truncation" << sim::time << ".dat";
   out.open(fname.str().c_str());
   if (!out) {
       *sim::log << "couldn't open tecplot output file " << fname.str();
       exit(1);
   }

   out << "ZONE F=FEPOINT, ET=TRIANGLE, N = " << nvrtx << ", E = " << ntri << std::endl;
   
   /* VERTEX MODES */
   for(i=0;i<nvrtx;++i) {
      for(n=0;n<ND;++n)
         out << vrtx(i)(n) << ' ';
      out << fscr1(i) << '\n';                    
   }
   
   

//   /* TO RENORMALIZE */
//   for(i=0;i<nvrtx;++i)
//      fscr1(i) = log10(fscr1(i)/(vd(i).nnbor*trncerr));

   /* OUTPUT CONNECTIVY INFO */
   out << std::endl << "#CONNECTION DATA#" << std::endl;

   for(tind=0;tind<ntri;++tind)
      out << td(tind).vrtx(0)+1 << ' ' << td(tind).vrtx(1)+1 << ' ' << td(tind).vrtx(2)+1 << '\n';
 
   out.close();
   
   return;
}
