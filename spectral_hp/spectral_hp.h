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

class spectral_hp : public mesh  {
   protected:
      int size;
      
      hpbasis b;
      int p0, sm0, im0;  // INITIALIZATION VALUES 
      
/*		SOLUTION INFORMATION */
      FLT (*vug)[NV], (*sug)[NV], (*iug)[NV];
      
/*		BOUNDARY INFORMATION */
/*		NOT ALL BOUNDARIES NEED THIS INFO, BUT EXTRA STORAGE IS MINISCULE */
      struct bistruct {
         FLT flx[NV];
         FLT curv[ND];
         FLT sfct;
      };
      struct bistruct *binfo[MAXSB];
      
/*		STATIC WORK ARRAYS */
      static FLT **u[NV],**du[NV][ND],**res[NV];
      static FLT **crd[ND], **dcrd[ND][ND], **cjcb;
      static FLT *uht[NV], *lf[NV], *lf1[NV];
      static int pmax;
      
/*		FUNCTIONS FOR MOVING GLOBAL TO LOCAL */
      void ugtouht(int tind);
      void ugtouht_bdry(int tind);
      void ugtouht1d(int sind);
      void crdtouht1d(int sind);
      void crdtouht(int tind);
/*		THIS FUNCTION ADDS LF TO GLOBAL VECTORS */
      void lftog(int tind, FLT (*v)[NV], FLT (*s)[NV], FLT (*i)[NV]);

/*		SETUP V/S/T INFO */
      void setbcinfo();
                  
   public:
      spectral_hp() : mesh() , size(0) {}
      spectral_hp(const spectral_hp& copy) : mesh(copy), size(0) { *this = copy;}
      spectral_hp& operator=(const spectral_hp& copy);
      void allocate(class hpbasis& bas);
      inline void loadbasis(class hpbasis& bas) { b = bas;}
      void tobasis(FLT (*func)(int, FLT, FLT));
      void curvinit();
      void input(char *name, FILETYPE type);
      void output(char *name, FILETYPE type);
      void density(FLT terr, FLT min, FLT max);
      void adapt(class spectral_hp& bgn, FLT tolerance);
      int ptprobe(FLT xp, FLT yp, FLT u[NV]);
      int ptprobe1d(int typ, FLT xp, FLT yp, FLT uout[NV]);
      FLT bdry_locate(int type, FLT& x, FLT &y, int& sind);
};