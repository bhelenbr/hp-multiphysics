/*
 *  gtol.cpp
 *  planar++
 *
 *  Created by helenbrk on Sun Oct 14 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

#include "spectral_hp.h"

void spectral_hp::ugtouht(int tind) {
	static int i,k,m,n,indx,cnt;
	static int sign, msgn;

/*	THIS IS FOR FLOW VARIABLES ON ANY MESH */
/*	VERTICES */	
	for (i=0; i<3; ++i) {
		indx = tvrtx[tind][i];
		for(n=0; n<NV; ++n)
			uht[n][i] = ug.v[indx][n];
	}
	

/* SIDES */
/* SIDE ORDERING OF MESH IS DIFFERENT THAN BASIS BASIS 0,1,2 = MESH 2,0,1 */
   cnt = 3;
   indx = tside[tind].side[2]*sm0;
   sign = tside[tind].sign[2];
   msgn = 1;
   for (m = 0; m < b.sm; ++m) {
      for(n=0; n<NV; ++n)
         uht[n][cnt] = msgn*ug.s[indx +m][n];
      msgn *= sign;
      ++cnt;
   }
   
   indx = tside[tind].side[0]*sm0;
   sign = tside[tind].sign[0];
   msgn = 1;
   for (m = 0; m < b.sm; ++m) {
      for(n=0; n<NV; ++n)
         uht[n][cnt] = msgn*ug.s[indx +m][n];
      msgn *= sign;
      ++cnt;
   }
   
   indx = tside[tind].side[1]*sm0;
   sign = tside[tind].sign[1];
   msgn = 1;
   for (m = 0; m < b.sm; ++m) {
      for(n=0; n<NV; ++n)
         uht[n][cnt] = msgn*ug.s[indx +m][n];
      msgn *= sign;
      ++cnt;
   }

/*	INTERIORS */	
   if (b.im > 0) {	
      indx = tind*im0;
      for(m = 1; m < b.sm; ++m) {
         for(k = 0; k < b.sm-m; ++k) {
            for(n=0; n<NV; ++n)
               uht[n][cnt] = ug.i[indx][n];
            ++cnt; ++indx;
         }
         indx += sm0 -b.sm;
      }
   }
	
	return;
}

void spectral_hp::ugtouht(int tind, struct vsi ug) {
	static int i,k,m,n,indx,cnt;
	static int sign, msgn;

/*	THIS IS FOR FLOW VARIABLES ON ANY MESH */
/*	VERTICES */	
	for (i=0; i<3; ++i) {
		indx = tvrtx[tind][i];
		for(n=0; n<NV; ++n)
			uht[n][i] = ug.v[indx][n];
	}
	

/* SIDES */
/* SIDE ORDERING OF MESH IS DIFFERENT THAN BASIS BASIS 0,1,2 = MESH 2,0,1 */
   cnt = 3;
   indx = tside[tind].side[2]*sm0;
   sign = tside[tind].sign[2];
   msgn = 1;
   for (m = 0; m < b.sm; ++m) {
      for(n=0; n<NV; ++n)
         uht[n][cnt] = msgn*ug.s[indx +m][n];
      msgn *= sign;
      ++cnt;
   }
   
   indx = tside[tind].side[0]*sm0;
   sign = tside[tind].sign[0];
   msgn = 1;
   for (m = 0; m < b.sm; ++m) {
      for(n=0; n<NV; ++n)
         uht[n][cnt] = msgn*ug.s[indx +m][n];
      msgn *= sign;
      ++cnt;
   }
   
   indx = tside[tind].side[1]*sm0;
   sign = tside[tind].sign[1];
   msgn = 1;
   for (m = 0; m < b.sm; ++m) {
      for(n=0; n<NV; ++n)
         uht[n][cnt] = msgn*ug.s[indx +m][n];
      msgn *= sign;
      ++cnt;
   }

/*	INTERIORS */	
   if (b.im > 0) {	
      indx = tind*im0;
      for(m = 1; m < b.sm; ++m) {
         for(k = 0; k < b.sm-m; ++k) {
            for(n=0; n<NV; ++n)
               uht[n][cnt] = ug.i[indx][n];
            ++cnt; ++indx;
         }
         indx += sm0 -b.sm;
      }
   }
	
	return;
}

void spectral_hp::ugtouht_bdry(int tind) {
	static int i,m,n,indx,cnt;
	static int sign, msgn;
	
	for (i=0; i<3; ++i) {
		indx = tvrtx[tind][i];
		for(n=0; n<NV; ++n)
			uht[n][i] = ug.v[indx][n];
	}

/* SIDE ORDERING OF MESH IS DIFFERENT THAN BASIS BASIS 0,1,2 = MESH 2,0,1 */
   cnt = 3;
   indx = tside[tind].side[2]*sm0;
   sign = tside[tind].sign[2];
   msgn = 1;
   for (m = 0; m < b.sm; ++m) {
      for(n=0; n<NV; ++n)
         uht[n][cnt] = msgn*ug.s[indx +m][n];
      msgn *= sign;
      ++cnt;
   }
   
   indx = tside[tind].side[0]*sm0;
   sign = tside[tind].sign[0];
   msgn = 1;
   for (m = 0; m < b.sm; ++m) {
      for(n=0; n<NV; ++n)
         uht[n][cnt] = msgn*ug.s[indx +m][n];
      msgn *= sign;
      ++cnt;
   }
   
   indx = tside[tind].side[1]*sm0;
   sign = tside[tind].sign[1];
   msgn = 1;
   for (m = 0; m < b.sm; ++m) {
      for(n=0; n<NV; ++n)
         uht[n][cnt] = msgn*ug.s[indx +m][n];
      msgn *= sign;
      ++cnt;
   }
	
	return;
}

void spectral_hp::ugtouht1d(int sind) {
   static int m,n,indx,v0,v1;
   
   v0 = svrtx[sind][0];
   v1 = svrtx[sind][1];
   for(n=0;n<NV;++n) {
      uht[n][0] = ug.v[v0][n];
      uht[n][1] = ug.v[v1][n];
   }
 	
   indx = sind*sm0;
   for(m=0;m<b.sm;++m)
    for(n=0;n<NV;++n) 
      uht[n][m+2] = ug.s[indx+m][n];
         
   return;
}

void spectral_hp::crdtouht(int tind) {
   static int i,m,n,cnt,bnum,sind,indx;
   
/*	VERTICES */	
	for (i=0; i<3; ++i) {
		indx = tvrtx[tind][i];
		for(n=0; n<NV; ++n)
			uht[n][i] = vrtx[indx][n];
	}

	if (b.sm == 0) return;
   
/* SIDES */
   cnt = 3;
   for (i=0; i<3;++i) {	
      sind = tside[tind].side[(i+2)%3];
      if (sinfo[sind] < 0) {
         for(m=0;m<b.sm;++m) {
            for(n=0;n<NV;++n)
               uht[n][cnt] = 0.0;
            ++cnt;
         }
      }
      else {
         bnum = -stri[sind][1]/maxsbel -1;
         indx = (-stri[sind][1] -(bnum+1)*maxsbel)*sm0;
         for (m = 0; m < b.sm; ++m) {
            for(n=0; n<NV; ++n)
               uht[n][cnt] = binfo[bnum][indx+m].curv[n];
            ++cnt;
         }
      }
   }
   
   return;
}

void spectral_hp::crdtouht(int tind, FLT (*vrtx)[ND], struct bistruct **binfo) {
   static int i,m,n,cnt,bnum,sind,indx;
   
/*	VERTICES */	
	for (i=0; i<3; ++i) {
		indx = tvrtx[tind][i];
		for(n=0; n<NV; ++n)
			uht[n][i] = vrtx[indx][n];
	}

	if (b.sm == 0) return;
   
/* SIDES */
   cnt = 3;
   for (i=0; i<3;++i) {	
      sind = tside[tind].side[(i+2)%3];
      if (sinfo[sind] < 0) {
         for(m=0;m<b.sm;++m) {
            for(n=0;n<NV;++n)
               uht[n][cnt] = 0.0;
            ++cnt;
         }
      }
      else {
         bnum = -stri[sind][1]/maxsbel -1;
         indx = (-stri[sind][1] -(bnum+1)*maxsbel)*sm0;
         for (m = 0; m < b.sm; ++m) {
            for(n=0; n<NV; ++n)
               uht[n][cnt] = binfo[bnum][indx+m].curv[n];
            ++cnt;
         }
      }
   }
   
   return;
}

void spectral_hp::crdtouht1d(int sind) {
   static int m,n,bnum,indx,v0,v1;
   
   v0 = svrtx[sind][0];
   v1 = svrtx[sind][1];
   for(n=0;n<ND;++n) {
      uht[n][0] = vrtx[v0][n];
      uht[n][1] = vrtx[v1][n];
   }
   
   bnum = -stri[sind][1]/maxsbel -1;
   indx = (-stri[sind][1] -(bnum+1)*maxsbel)*sm0;
   
   for(m=0;m<b.sm;++m)
    for(n=0;n<ND;++n) 
         uht[n][m+2] = binfo[bnum][indx+m].curv[n];
         
   return;
}

void spectral_hp::crdtouht1d(int sind,FLT (*vrtx)[ND], struct bistruct **binfo) {
   static int m,n,bnum,indx,v0,v1;
   
   v0 = svrtx[sind][0];
   v1 = svrtx[sind][1];
   for(n=0;n<ND;++n) {
      uht[n][0] = vrtx[v0][n];
      uht[n][1] = vrtx[v1][n];
   }
   
   bnum = -stri[sind][1]/maxsbel -1;
   indx = (-stri[sind][1] -(bnum+1)*maxsbel)*sm0;
   
   for(m=0;m<b.sm;++m)
    for(n=0;n<ND;++n) 
         uht[n][m+2] = binfo[bnum][indx+m].curv[n];
         
   return;
}

void spectral_hp::lftog(int tind, struct vsi g) {
   static int m,n,indx,gindx,sgn,msgn;
   
/*	VERTEX MODES */
   for (m=0; m<3; ++m) {
      gindx = tvrtx[tind][m];
      for(n=0;n<NV;++n)
         g.v[gindx][n] += lf[n][m];
   }


   if (b.p > 1) {
/* 	SIDE MODES */
/* 	SIDE ORDERING OF MESH IS DIFFERENT THAN BASIS BASIS 0,1,2 = MESH 2,0,1 */
      indx = 3;
      gindx = tside[tind].side[2]*b.sm;
      sgn = tside[tind].sign[2];
      msgn = 1;
      for (m = 0; m < b.sm; ++m) {
         for(n=0;n<NV;++n)
            g.s[gindx][n] += msgn*lf[n][indx];
         msgn *= sgn;
         ++gindx;
         ++indx;
      }
      
      gindx = tside[tind].side[0]*b.sm;
      sgn = tside[tind].sign[0];
      msgn = 1;
      for (m = 0; m < b.sm; ++m) {
         for(n=0;n<NV;++n)
            g.s[gindx][n] += msgn*lf[n][indx];
         msgn *= sgn;
         ++gindx;
         ++indx;
      }
      
      gindx = tside[tind].side[1]*b.sm;
      sgn = tside[tind].sign[1];
      msgn = 1;
      for (m = 0; m < b.sm; ++m) {
         for(n=0;n<NV;++n)
            g.s[gindx][n] += msgn*lf[n][indx];
         msgn *= sgn;
         ++gindx;
         ++indx;
      }

      gindx = tind*b.im;
      indx = b.bm;
      for(m=0;m<b.im;++m) {
         for(n=0;n<NV;++n)
            g.i[gindx][n] += lf[n][indx];
         ++gindx;
         ++indx;
      }
   }
   
   return;
}
