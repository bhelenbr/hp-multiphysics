/*
 *  spectral_hp.h
 *  planar++
 *
 *  Created by helenbrk on Mon Oct 01 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#include<mesh.h>
#include<float.h>
#include"hpbasis.h"

#define NV 3

/* SOLUTION VECTOR */
struct vsi {
   FLT (*v)[NV];
   FLT (*s)[NV];
   FLT (*i)[NV];
};

/* BOUNDARY INFORMATION */
struct bistruct {
   FLT flx[NV];
   FLT curv[ND];
};

class spectral_hp : public mesh  {
   protected:
      int size;
      
      hpbasis b;
      int p0, sm0, im0;  // INITIALIZATION VALUES 
      
/*		SOLUTION INFORMATION */
      struct vsi ug;
      
/*		BOUNDARY INFORMATION */
/*		NOT ALL BOUNDARIES NEED THIS INFO, BUT EXTRA STORAGE IS MINISCULE */
      struct bistruct *binfo[MAXSB];
      
/*		STATIC WORK ARRAYS */
      static FLT **u[NV],**du[NV][ND],**res[NV];
      static FLT **crd[ND], **dcrd[ND][ND], **cjcb;
      static FLT *uht[NV], *lf[NV], *lf1[NV];
      static int pmax;
      
/*		FUNCTIONS FOR MOVING GLOBAL TO LOCAL */
      void ugtouht(int tind);
      void ugtouht(int tind, struct vsi ug);
      void ugtouht_bdry(int tind);
      void ugtouht_bdry(int tind, struct vsi ug);
      void ugtouht1d(int sind);
      void ugtouht1d(int sind, struct vsi ug);
      void crdtouht(int tind);
      void crdtouht(int tind, FLT (*vrtx)[ND], struct bistruct **binfo);
      void crdtouht1d(int sind);
      void crdtouht1d(int sind, FLT (*vrtx)[ND], struct bistruct **binfo);

/*		THIS FUNCTION ADDS LF TO GLOBAL VECTORS */
      void lftog(int tind, struct vsi);

/*		SETUP V/S/T INFO */
      void setbcinfo();
                  
   public:
      spectral_hp() : mesh() , size(0) {}
      void copy(const spectral_hp& copy);
      void allocate(class hpbasis& bas);
      inline void loadbasis(class hpbasis& bas) { b = bas;}
      void tobasis(struct vsi g, FLT (*func)(int, FLT, FLT));
      inline void tobasis(FLT (*func)(int, FLT, FLT)) {tobasis(ug,func);}
      void curvinit();
      void input(struct vsi g, char *name, FILETYPE type);
      inline void input(char *name, FILETYPE type) {input(ug,name,type);}
      void output(struct vsi g, char *name, FILETYPE type);
      inline void output(char *name, FILETYPE type) {output(ug,name,type);}
      void ptprobe(FLT xp, FLT yp, FLT u[NV]);
      void ptprobe1d(int typ, FLT xp, FLT yp, FLT uout[NV]);
      int findinteriorpt(FLT xp, FLT yp, FLT &r, FLT &s);
      int findbdrypt(int typ, FLT &x, FLT &y, FLT &psi);
};