/*
 *  spectral_hp.h
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 01 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#include <r_mesh.h>
#include <float.h>
#include <hpbasis.h>
#include "defines.h"

/* SOLUTION VECTOR */
struct vsi {
   FLT (*v)[NV];
   FLT (*s)[NV];
   FLT (*i)[NV];
};

/* struct vsi {
   Array<TinyVector<FLT,NV>,1> v;
   Array<TinyVector<FLT,NV>,1> s;
   Array<TinyVector<FLT,NV>,1> i;
}; */

/* BOUNDARY INFORMATION */
struct bistruct {
   FLT flx[NV];
   FLT curv[ND];
};

class spectral_hp : public r_mesh  {
   protected:
      int size;
      
      hpbasis *b;
      int p0, sm0, im0;  // INITIALIZATION VALUES 
      
      /* SOLUTION INFORMATION */
      struct vsi ug;
      
      /* BOUNDARY INFORMATION */
      /* NOT ALL BOUNDARIES NEED THIS INFO, BUT EXTRA STORAGE IS MINISCULE */
   public:
      struct bistruct *binfo[MAXSB];
      
   protected:
      /* STATIC WORK ARRAYS */
      static FLT u[NV][MXGP][MXGP],du[NV][ND][MXGP][MXGP],res[NV][MXGP][MXGP];
      static FLT crd[ND][MXGP][MXGP], dcrd[ND][ND][MXGP][MXGP], cjcb[MXGP][MXGP];
      static FLT uht[NV][MXTM], lf[NV][MXTM];
      static FLT cht[ND][MXTM], cf[ND][MXTM];
      static int pmax;
      
      /* FUNCTIONS FOR MOVING GLOBAL TO LOCAL */
      void ugtouht(int tind);
      void ugtouht(int tind, struct vsi ug);
      void ugtouht_bdry(int tind);
      void ugtouht_bdry(int tind, struct vsi ug);
      void ugtouht1d(int sind);
      void ugtouht1d(int sind, struct vsi ug);
      void crdtocht(int tind);
      void crdtocht(int tind, FLT (*vrtx)[ND], struct bistruct **binfo);
      void crdtocht1d(int sind);
      void crdtocht1d(int sind, FLT (*vrtx)[ND], struct bistruct **binfo);

      /* THIS FUNCTION ADDS LF TO GLOBAL VECTORS */
      void lftog(int tind, struct vsi);

      /* SETUP V/S/T INFO */
      void setbcinfo();
                  
   public:
      spectral_hp() : r_mesh() , size(0) {}
      void copy(const spectral_hp& copy);
      void allocate(class hpbasis *bas);
      inline void loadbasis(class hpbasis *bas) { b = bas;}
      void tobasis(struct vsi g, FLT (*func)(int, FLT, FLT));
      void l2error(FLT (*func)(int, FLT, FLT));
      FLT findmax(int type, FLT (*fxy)(FLT x[ND]));
      FLT findmaxx(int type);
      FLT findmaxy(int type);
      inline void tobasis(FLT (*func)(int, FLT, FLT)) {tobasis(ug,func);}
      void curvinit(int MASK = ~0);
      void input(struct vsi g, FLT (*vin)[ND], struct bistruct **bin, char *name, FILETYPE type);
      inline void input(char *name, FILETYPE type) {input(ug,vrtx,binfo,name,type);}
      inline void input(struct vsi g, char *name, FILETYPE type) {input(g,vrtx,binfo,name,type);}
      void output(struct vsi g, FLT (*vin)[ND], struct bistruct **bin, char *name, FILETYPE typ = tecplot);
      inline void output(struct vsi g, char *name, FILETYPE type=tecplot) {output(g,vrtx,binfo,name,type);}
      inline void output(char *name, FILETYPE type=tecplot) {output(ug,vrtx,binfo,name,type);}
      void ptprobe(FLT xp, FLT yp, FLT u[NV]);
      void ptprobe1d(int typ, FLT xp, FLT yp, FLT uout[NV]);
      int findinteriorpt(FLT xp, FLT yp, FLT &r, FLT &s);
      int findandmvptincurved(FLT &xp, FLT &yp, FLT &r, FLT &s);
      int findbdrypt(int typ, FLT &x, FLT &y, FLT &psi);
};

