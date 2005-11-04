#include "tri_hp.h"
#include "boundary.h"
#include "hp_boundary.h"
#include <assert.h>

 void tri_hp::ptprobe(TinyVector<FLT,2> xp, Array<FLT,1> uout, int tlvl) {
   FLT r,s;
   int tind;
   
   findinteriorpt(xp,tind,r,s);
   ugtouht(tind,tlvl);  
   basis::tri(log2p).ptprobe(NV,uout.data(),&uht(0)(0),MXTM);
}

 void tri_hp::ptprobe_bdry(int bnum, TinyVector<FLT,2> xp, Array<FLT,1> uout,int tlvl) {
   FLT psi;
   int sind;
   
   hp_sbdry(bnum)->findbdrypt(xp,sind,psi);
   ugtouht1d(sind,tlvl);  
   basis::tri(log2p).ptprobe1d(NV,uout.data(),&uht(0)(0),MXTM);
}

 void tri_hp::findandmvptincurved(TinyVector<FLT,2> xp, int &tind, FLT &r, FLT &s) {
   TinyVector<FLT,3> wgt;
   int v0;
   
   qtree.nearpt(xp.data(),v0);
   tind = findtri(xp,v0);
   getwgts(wgt);
   
   if (tind < 0) {
      *sim::log << "#Warning: couldn't find tri " << xp << " nearpt " << v0 << " neartri " << tind << std::endl;
      tind = abs(tind);
   }

   /* TRIANGLE COORDINATES */   
   s = wgt(0)*2 -1.0;
   r = wgt(2)*2 -1.0;
   
   if (td(tind).info < 0) {
      basis::tri(log2p).ptvalues(r,s);
      return;
   }
   
   *sim::log << "#In find and move " << xp << std::endl;

   /* MOVE POINT WITH SIDE CURVATURE */
   crdtocht(tind);
   basis::tri(log2p).ptprobe_bdry(ND,xp.data(),r,s,&cht(0,0),MXTM);

   return;
}


 void tri_hp::findinteriorpt(TinyVector<FLT,ND> xp, int &tind, FLT &r, FLT &s) {
   FLT dr,ds,dx,dy,x[ND],ddr[3],dds[3],det,roundoff,dxmax[ND];
   TinyVector<FLT,3> wgt;
   int n,iter,v0;

   qtree.nearpt(xp.data(),v0);
   tind = findtri(xp,v0);
   getwgts(wgt);

   if (tind < 0) {
      *sim::log << "#Warning: couldn't find boundary tri " << xp << " nearpt " << v0 << " neartri " << tind << std::endl;
      tind = abs(tind);
   }

   /* TRIANGLE COORDINATES */   
   s = wgt[0]*2 -1.0;
   r = wgt[2]*2 -1.0;
   
   if (td(tind).info < 0) {
      basis::tri(log2p).ptvalues(r,s);
      return;
   }

   /* DEAL WITH CURVED SIDES */
   crdtocht(tind);
   
   for(n=0;n<ND;++n)
      dxmax[n] = fabs(cht(n,0)-cht(n,1)) +fabs(cht(n,1)-cht(n,2));
   roundoff = 10.0*EPSILON*(1.0 +(fabs(xp(0))*dxmax[1] +fabs(xp(1))*dxmax[0])/area(tind));
   
   iter = 0;
   do {
      basis::tri(log2p).ptprobe_bdry(ND,x,ddr,dds,r,s,&cht(0,0),MXTM);
      det = 1.0/(fabs(ddr[0]*dds[1] - ddr[1]*dds[0]) +10.0*EPSILON);
      dx = xp(0)-x[0];
      dy = xp(1)-x[1];
      dr =  (dds[1]*dx -dds[0]*dy)*det;
      ds = -(ddr[1]*dx -ddr[0]*dy)*det;

      r += dr;
      s += ds;
      if (iter++ > 100) {
         *sim::log << "#Warning: max iterations for curved triangle " << tind << " loc: " << xp << std::endl;
         break;
      }
   } while (fabs(dr) +fabs(ds) > roundoff);

   return;
}