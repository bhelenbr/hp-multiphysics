/*
 *  gtol.cpp
 *  planar++
 *
 *  Created by helenbrk on Sun Oct 14 2001.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */
 
#include "tet_hp.h"
#include "hp_boundary.h"
#include <assert.h>

//#define NODAL			

 /* Global to Local */
 void tet_hp::ugtouht(int tind) {
    int i,k,m,n,indx,eind,cnt,find;
    int sign, msgn;


   /* THIS IS FOR FLOW VARIABLES ON ANY MESH */
   /* VERTICES */   
   for (i = 0; i < 4; ++i) {
      indx = tet(tind).pnt(i);
      for(n = 0; n < NV; ++n)
         uht(n)(i) = ug.v(indx,n);
   }   

   /* EDGES */
   if (basis::tet(log2p).em > 0){
	   cnt = 4;
	   for(i = 0; i < 6; ++i) {
		  eind = tet(tind).seg(i);
		  sign = tet(tind).sgn(i);
		  msgn = 1;
		  for (m = 0; m < basis::tet(log2p).em; ++m) {
			 for(n = 0; n < NV; ++n)
				uht(n)(cnt) = msgn*ug.e(eind,m,n);      
			 msgn *= sign;
			 ++cnt;
		  }
	   }
   }
 
    /* FACES */   
	if (basis::tet(log2p).fm > 0) {
		for(int i = 0; i < 4; ++i) { 
			find = tet(tind).tri(i);
			sign = tet(tind).rot(i);
			msgn = 1;
			indx = 0;
			for(m = 1; m <= basis::tet(log2p).em-1; ++m) {
				for(k = 1; k <= basis::tet(log2p).em-m; ++k) {
					for(n = 0; n < NV; ++n)
						uht(n)(cnt) = msgn*ug.f(find,indx,n);
					++cnt; ++indx;
				}
				msgn *= sign;
				indx += em0 -basis::tet(log2p).em;// index shift for p-mg
			}
		}
	}
	
	/* INTERIOR */
	if (basis::tet(log2p).im > 0) {
		indx = 0;
		for(m = 1; m <= basis::tet(log2p).em-1; ++m) {
			for(k = 1; k <= basis::tet(log2p).em-m; ++k) {
				for(i = 1; i <= basis::tet(log2p).em-m-k; ++i){
					for(n = 0; n < NV; ++n){
						uht(n)(cnt) = ug.i(tind,indx,n);
					}
					++cnt; ++indx;
				}
				indx += em0 -basis::tet(log2p).em;				
			}

		}
	}
	
	
   
   return;
}

 void tet_hp::ugtouht(int tind, int tlvl) {
    int i,k,m,n,eind,indx,cnt,find;
    int sign, msgn;
    vefi &ug = ugbd(tlvl);
	
    
   /* THIS IS FOR FLOW VARIABLES ON ANY MESH */
   /* VERTICES */   
   for (i = 0; i < 4; ++i) {
      indx = tet(tind).pnt(i);
      for(n = 0; n < NV; ++n)
         uht(n)(i) = ug.v(indx,n);
   }   

   /* EDGES */
   if (basis::tet(log2p).em > 0){
	   cnt = 4;
	   for(i = 0; i < 6; ++i) {
		  eind = tet(tind).seg(i);
		  sign = tet(tind).sgn(i);
		  msgn = 1;
		  for (m = 0; m < basis::tet(log2p).em; ++m) {
			 for(n = 0; n < NV; ++n)
				uht(n)(cnt) = msgn*ug.e(eind,m,n);
			 msgn *= sign;
			 ++cnt;
		  }
	   }
   }
 
    /* FACES */   
	if (basis::tet(log2p).fm > 0) {
		for(int i = 0; i < 4; ++i) { 
			find = tet(tind).tri(i);
			sign = tet(tind).rot(i);
			msgn = 1;
			indx = 0;
			for(m = 1; m <= basis::tet(log2p).em-1; ++m) {
				for(k = 1; k <= basis::tet(log2p).em-m; ++k) {
					for(n = 0; n < NV; ++n)
						uht(n)(cnt) = msgn*ug.f(find,indx,n);
					++cnt; ++indx;
				}
				msgn *= sign;
				indx += em0 -basis::tet(log2p).em;				
			}
		}
	}
	
	/* INTERIOR */
	if (basis::tet(log2p).im > 0) {
		indx = 0;
		for(m = 1; m <= basis::tet(log2p).em-1; ++m) {
			for(k = 1; k <= basis::tet(log2p).em-m; ++k) {
				for(i = 1; i <= basis::tet(log2p).em-m-k; ++i){
					for(n = 0; n < NV; ++n)
						uht(n)(cnt) = ug.i(tind,indx,n);
					++cnt; ++indx;
				}
				indx += em0 -basis::tet(log2p).em;				
			}
		}
	}
   
   return;
}

 void tet_hp::ugtouht_bdry(int tind) {
    int i,m,n,k,indx,eind,find,cnt;
    int sign, msgn;
	
   
   /* VERTICES */   
   for (i = 0; i < 4; ++i) {
      indx = tet(tind).pnt(i);
      for(n = 0; n < NV; ++n)
         uht(n)(i) = ug.v(indx,n);
   }
   

   /* EDGES */
   if (basis::tet(log2p).em > 0){
	   cnt = 4;
	   for(i = 0; i < 6; ++i) {
		  eind = tet(tind).seg(i);
		  sign = tet(tind).sgn(i);
		  msgn = 1;
		  for (m = 0; m < basis::tet(log2p).em; ++m) {
			 for(n = 0; n < NV; ++n)
				uht(n)(cnt) = msgn*ug.e(eind,m,n);
			 msgn *= sign;
			 ++cnt;
		  }
	   }
   }
 
    /* FACES */   
	if (basis::tet(log2p).fm > 0) {
		for(int i = 0; i < 4; ++i) { 
			find = tet(tind).tri(i);
			sign = tet(tind).rot(i);
			msgn = 1;
			indx = 0;
			for(m = 1; m <= basis::tet(log2p).em-1; ++m) {
				for(k = 1; k <= basis::tet(log2p).em-m; ++k) {
					for(n = 0; n < NV; ++n)
						uht(n)(cnt) = msgn*ug.f(find,indx,n);
					++cnt; ++indx;
				}
				msgn *= sign;
				indx += em0 -basis::tet(log2p).em;				
			}
		}
	}
	
   return;
}

 void tet_hp::ugtouht_bdry(int tind, int tlvl) {
    int i,m,n,k,indx,eind,find,cnt;
    int sign, msgn;
    vefi &ug = ugbd(tlvl);
	
   
   /* VERTICES */   
   for (i = 0; i < 4; ++i) {
      indx = tet(tind).pnt(i);
      for(n = 0; n < NV; ++n)
         uht(n)(i) = ug.v(indx,n);
   }   

   /* EDGES */
   if (basis::tet(log2p).em > 0){
	   cnt = 4;
	   for(i = 0; i < 6; ++i) {
		  eind = tet(tind).seg(i);
		  sign = tet(tind).sgn(i);
		  msgn = 1;
		  for (m = 0; m < basis::tet(log2p).em; ++m) {
			 for(n = 0; n < NV; ++n)
				uht(n)(cnt) = msgn*ug.e(eind,m,n);
			 msgn *= sign;
			 ++cnt;
		  }
	   }
   }
 
    /* FACES */   
	if (basis::tet(log2p).fm > 0) {
		for(int i = 0; i < 4; ++i) { 
			find = tet(tind).tri(i);
			sign = tet(tind).rot(i);		
			msgn = 1;
			indx = 0;
			for(m = 1; m <= basis::tet(log2p).em-1; ++m) {
				for(k = 1; k <= basis::tet(log2p).em-m; ++k) {
					for(n = 0; n < NV; ++n)
						uht(n)(cnt) = msgn*ug.f(find,indx,n);
					++cnt; ++indx;
				}
				msgn *= sign;
				indx += em0 -basis::tet(log2p).em;				
			}
		}
	}
	
	return;
}


 void tet_hp::ugtouht2d(int find) {
   int i,m,n,k,v1,v2,v3,indx,cnt,sign,eind,msgn;
   
   v1 = tri(find).pnt(0);
   v2 = tri(find).pnt(1);
   v3 = tri(find).pnt(2);
   for(n=0;n<NV;++n) {
	  uht(n)(0) = ug.v(v1,n);
      uht(n)(1) = ug.v(v2,n);
      uht(n)(2) = ug.v(v3,n);
   }
    
   /* EDGES */
   cnt = 3;
   if (basis::tet(log2p).em > 0) {
	   for(i = 0; i < 3; ++i) {
		  eind = tri(find).seg(i);
		  sign = tri(find).sgn(i);
		  msgn = 1;
		  for (m = 0; m < basis::tet(log2p).em; ++m) {
			 for(n = 0; n < NV; ++n)
				uht(n)(cnt) = msgn*ug.e(eind,m,n);
			 msgn *= sign;
			 ++cnt;
		  }
	   }
   }
 
    /* FACE */   
	if (basis::tet(log2p).fm > 0) {
		indx = 0;
		for(m = 1; m <= basis::tet(log2p).em-1; ++m) {
			for(k = 1; k <= basis::tet(log2p).em-m; ++k) {
				for(n = 0; n < NV; ++n)
					uht(n)(cnt) = ug.f(find,indx,n);
				++cnt; ++indx;
			}
			indx += em0 -basis::tet(log2p).em;				
		}
	}
         
   return;
}

 void tet_hp::ugtouht2d(int find, int tlvl) {
   int i,k,m,n,v1,v2,v3,indx,cnt,eind,sign,msgn;
   vefi &ug = ugbd(tlvl);
   
   v1 = tri(find).pnt(0);
   v2 = tri(find).pnt(1);
   v3 = tri(find).pnt(2);
   for(n=0;n<NV;++n) {
	  uht(n)(0) = ug.v(v1,n);
      uht(n)(1) = ug.v(v2,n);
      uht(n)(2) = ug.v(v3,n);
   }
    
   /* EDGES */
   cnt = 3;
   if (basis::tet(log2p).em > 0) {
	   for(i = 0; i < 3; ++i) {
		  eind = tri(find).seg(i);
		  sign = tri(find).sgn(i);
		  msgn = 1;
		  for (m = 0; m < basis::tet(log2p).em; ++m) {
			 for(n = 0; n < NV; ++n)
				uht(n)(cnt) = msgn*ug.e(eind,m,n);
			 msgn *= sign;
			 ++cnt;
		  }
	   }
   }
 
    /* FACE */   
	if (basis::tet(log2p).fm > 0) {
		indx = 0;
		for(m = 1; m <= basis::tet(log2p).em-1; ++m) {
			for(k = 1; k <= basis::tet(log2p).em-m; ++k) {
				for(n = 0; n < NV; ++n)
					uht(n)(cnt) = ug.f(find,indx,n);
				++cnt; ++indx;
			}
			indx += em0 -basis::tet(log2p).em;				
		}
	}   return;
}


 void tet_hp::ugtouht2d_bdry(int find) {
   int i,m,n,v1,v2,v3,cnt,sign,eind,msgn;
   
   v1 = tri(find).pnt(0);
   v2 = tri(find).pnt(1);
   v3 = tri(find).pnt(2);
   for(n=0;n<NV;++n) {
	  uht(n)(0) = ug.v(v1,n);
      uht(n)(1) = ug.v(v2,n);
      uht(n)(2) = ug.v(v3,n);
   }
    
   /* EDGES */
   cnt = 3;
   if (basis::tet(log2p).em > 0) {
	   for(i = 0; i < 3; ++i) {
		  eind = tri(find).seg(i);
		  sign = tri(find).sgn(i);
		  msgn = 1;
		  for (m = 0; m < basis::tet(log2p).em; ++m) {
			 for(n = 0; n < NV; ++n)
				uht(n)(cnt) = msgn*ug.e(eind,m,n);
			 msgn *= sign;
			 ++cnt;
		  }
	   }
   }
          
   return;
}

 void tet_hp::ugtouht2d_bdry(int find, int tlvl) {
   int i,m,n,v1,v2,v3,eind,cnt,sign,msgn;
   vefi &ug = ugbd(tlvl);
   
   v1 = tri(find).pnt(0);
   v2 = tri(find).pnt(1);
   v3 = tri(find).pnt(2);
   for(n=0;n<NV;++n) {
	  uht(n)(0) = ug.v(v1,n);
      uht(n)(1) = ug.v(v2,n);
      uht(n)(2) = ug.v(v3,n);
   }
    
   /* EDGES */
   cnt = 3;
   if (basis::tet(log2p).em > 0) {
	   for(i = 0; i < 3; ++i) {
		  eind = tri(find).seg(i);
		  sign = tri(find).sgn(i);
		  msgn = 1;
		  for (m = 0; m < basis::tet(log2p).em; ++m) {
			 for(n = 0; n < NV; ++n)
				uht(n)(cnt) = msgn*ug.e(eind,m,n);
			 msgn *= sign;
			 ++cnt;
		  }
	   }
   }
  return;
}



 void tet_hp::ugtouht1d(int eind) {
   int m,n,v0,v1;
   
   v0 = seg(eind).pnt(0);
   v1 = seg(eind).pnt(1);
   for(n=0;n<NV;++n) {
      uht(n)(0) = ug.v(v0,n);
      uht(n)(1) = ug.v(v1,n);
   }
    
   for(m=0;m<basis::tet(log2p).em;++m)
    for(n=0;n<NV;++n) 
      uht(n)(m+2) = ug.e(eind,m,n);
         
   return;
}

 void tet_hp::ugtouht1d(int eind, int tlvl) {
   int m,n,v0,v1;
   vefi &ug = ugbd(tlvl);
   
   v0 = seg(eind).pnt(0);
   v1 = seg(eind).pnt(1);
   for(n=0;n<NV;++n) {
      uht(n)(0) = ug.v(v0,n);
      uht(n)(1) = ug.v(v1,n);
   }
    
   for(m=0;m<basis::tet(log2p).em;++m)
    for(n=0;n<NV;++n) 
      uht(n)(m+2) = ug.e(eind,m,n);
       
   return;
}

 void tet_hp::crdtocht(int tind) {
   int i,m,n,cnt,bnum,eind,find,indx;
   
   /* VERTICES */   
   for (i=0; i < 4; ++i) {
      indx = tet(tind).pnt(i);
      for(n=0; n<ND; ++n)
         cht(n)(i) = pnts(indx)(n);
	}
	
	if (basis::tet(log2p).em == 0) return;
	
	/* EDGES */
	cnt = 4;
	for (i = 0; i < 6; ++i){
		eind=tet(tind).seg(i);
		if(seg(eind).info < 0){
			for(m = 0; m < basis::tet(log2p).em; ++m) {
				for(n = 0; n < ND; ++n)
				   cht(n)(cnt) = 0.0;
				++cnt;
			 }
		}
		else {
			bnum = getbdrynum(seg(eind).info);
			indx = getbdryseg(seg(eind).info);
			for(m = 0; m < basis::tet(log2p).em; ++m) {
				for(n = 0; n < ND; ++n)
					cht(n)(cnt) = hp_ebdry(bnum)->crde(indx,m,n);/// temporary fix crde(indx,m,n)
				++cnt;
			 }		
		}
	}
	
	if (basis::tet(log2p).fm == 0) return;

	/* FACES */
	for (i = 0; i < 4; ++i){
		find=tet(tind).tri(i);
		if(tri(find).info < 0){
			for(m = 0; m < basis::tet(log2p).fm; ++m) {
				for(n = 0; n < ND; ++n)
				   cht(n)(cnt) = 0.0;
				++cnt;
			 }
		}
		else {
			bnum = getbdrynum(tri(find).tet(1));
			indx = getbdryseg(tri(find).tet(1));
			for(m = 0; m < basis::tet(log2p).fm; ++m) {
				for(n = 0; n < ND; ++n)
					cht(n)(cnt) = hp_fbdry(bnum)->crdf(indx,m,n);/// temporary fix crdf(indx,m,n)
				++cnt;
			 }		
		}
	}
   
   return;
}


 void tet_hp::crdtocht(int tind, int tlvl) {
   int i,m,n,cnt,bnum,eind,find,indx;
      
      /* VERTICES */   
   for (i=0; i < 4; ++i) {
      indx = tet(tind).pnt(i);
      for(n=0; n<ND; ++n)
         cht(n)(i) = vrtxbd(tlvl)(indx)(n);
	}
	
	if (basis::tet(log2p).em == 0) return;
	
	/* EDGES */
	cnt = 4;
	for (i = 0; i < 6; ++i){
		eind=tet(tind).seg(i);
		if(seg(eind).info < 0){
			for(m = 0; m < basis::tet(log2p).em; ++m) {
				for(n = 0; n < ND; ++n)
				   cht(n)(cnt) = 0.0;
				++cnt;
			 }
		}
		else {
			bnum = getbdrynum(seg(eind).info);
			indx = getbdryseg(seg(eind).info);
			for(m = 0; m < basis::tet(log2p).em; ++m) {
				for(n = 0; n < ND; ++n)
					cht(n)(cnt) = hp_ebdry(bnum)->crdebd(tlvl,indx,m,n);/// fix crdebd(tlvl,indx,m,n)
				++cnt;
			 }		
		}
	}
	
	if (basis::tet(log2p).fm == 0) return;

	/* FACES */
	for (i = 0; i < 4; ++i){
		find=tet(tind).tri(i);
		//cout << endl << tri(find).info << endl;
		if(tri(find).info < 0){
			for(m = 0; m < basis::tet(log2p).fm; ++m) {
				for(n = 0; n < ND; ++n)
				   cht(n)(cnt) = 0.0;
				++cnt;
			 }
		}
		else {
			bnum = getbdrynum(tri(find).tet(1));
			indx = getbdryseg(tri(find).tet(1));
			for(m = 0; m < basis::tet(log2p).fm; ++m) {
				for(n = 0; n < ND; ++n)
					cht(n)(cnt) = hp_fbdry(bnum)->crdfbd(tlvl,indx,m,n);/// fix crdf(indx,m,n)
				++cnt;
			 }		
		}
	}
	
	return;
}


 void tet_hp::crdtocht2d(int find){
   int i,m,n,cnt,bnum,indx,eind;
   
   /* VERTICES */   
   for (i=0; i < 3; ++i) {
      indx = tri(find).pnt(i);
      for(n=0; n<ND; ++n)
         cht(n)(i) = pnts(indx)(n);
   }
   
   	if (basis::tet(log2p).em == 0) return;
	
	/* EDGES */
	cnt = 3;
	for (i = 0; i < 3; ++i){
		eind=tri(find).seg(i);
		if(seg(eind).info < 0){
			for(m = 0; m < basis::tet(log2p).em; ++m) {
				for(n = 0; n < ND; ++n)
				   cht(n)(cnt) = 0.0;
				++cnt;
			 }
		}
		else {
			bnum = getbdrynum(seg(eind).info);
			indx = getbdryseg(seg(eind).info);
			for(m = 0; m < basis::tet(log2p).em; ++m) {
				for(n = 0; n < ND; ++n)
					cht(n)(cnt) = hp_ebdry(bnum)->crde(indx,m,n);/// fix crde(indx,m,n)
				++cnt;
			 }		
		}
	}
	
  
   return;
}


 void tet_hp::crdtocht2d(int find, int tlvl) {
   int i,m,n,cnt,bnum,eind,indx;
      
   /* VERTICES */   
   for (i=0; i < 3; ++i) {
		//cout << tri(find).pnt(i) << endl;
      indx = tri(find).pnt(i);
      for(n=0; n<ND; ++n)
         cht(n)(i) = vrtxbd(tlvl)(indx)(n);
   }
   
   	if (basis::tet(log2p).em == 0) return;
	
	/* EDGES */
	cnt = 3;
	for (i = 0; i < 3; ++i){
		eind=tri(find).seg(i);
		// if (eind != 0) cout << endl << eind << endl;
		if(seg(eind).info < 0){
			for(m = 0; m < basis::tet(log2p).em; ++m) {
				for(n = 0; n < ND; ++n)
				   cht(n)(cnt) = 0.0;
				++cnt;
			 }
		}
		else {
			bnum = getbdrynum(seg(eind).info);
			indx = getbdryseg(seg(eind).info);
			for(m = 0; m < basis::tet(log2p).em; ++m) {
				for(n = 0; n < ND; ++n)
					cht(n)(cnt) = hp_ebdry(bnum)->crdebd(tlvl,indx,m,n);/// fix crde(indx,m,n)
				++cnt;
			 }		
		}
	}	

	return;
}


 void tet_hp::crdtocht1d(int eind) {
   int m,n,bnum,indx,v0,v1;
   
   v0 = seg(eind).pnt(0);
   v1 = seg(eind).pnt(1);
   for(n=0;n<ND;++n) {
      cht(n)(0) = pnts(v0)(n);
      cht(n)(1) = pnts(v1)(n);
   }
   
   if (seg(eind).info < 0) {
      for(m=0;m<basis::tet(log2p).em;++m)
         for(n=0;n<ND;++n) 
            cht(n)(m+2) = 0.0;
   }
   else {
        bnum = getbdrynum(seg(eind).info);
        indx = getbdryseg(seg(eind).info);      
		for(m=0;m<basis::tet(log2p).em;++m)
			for(n=0;n<ND;++n) 
				cht(n)(m+2) = hp_ebdry(bnum)->crde(indx,m,n);     
 
   }
         
   return;
}

 void tet_hp::crdtocht1d(int eind,int tlvl) {
   int m,n,bnum,indx,v0,v1;
   
   v0 = seg(eind).pnt(0);
   v1 = seg(eind).pnt(1);
   for(n=0;n<ND;++n) {
      cht(n)(0) = vrtxbd(tlvl)(v0)(n);
      cht(n)(1) = vrtxbd(tlvl)(v1)(n);
   }
   //cout << seg(eind).info << endl;

   if (seg(eind).info < 0) {
      for(m=0;m<basis::tet(log2p).em;++m)
         for(n=0;n<ND;++n) 
            cht(n)(m+2) = 0.0;
   }
   else {
        bnum = getbdrynum(seg(eind).info);
        indx = getbdryseg(seg(eind).info);      
		for(m=0;m<basis::tet(log2p).em;++m)
			for(n=0;n<ND;++n) 
				cht(n)(m+2) = hp_ebdry(bnum)->crdebd(tlvl,indx,m,n);     
 
   }         
   return;
}

#ifndef NODAL 
 /* Local to Global */
 void tet_hp::lftog(int tind, struct vefi g) {
   int i,m,n,indx,gindx,eind,find,sgn,msgn;
   
   /* VERTEX MODES */
   for (m = 0; m < 4; ++m) {
      gindx = tet(tind).pnt(m);
      for (n = 0; n < NV; ++n)
         g.v(gindx,n) += lf(n)(m);
   }

	/* EDGE MODES */
	if (basis::tet(log2p).p > 1) {
		indx = 4;
		for(i = 0; i < 6; ++i) {
			eind = tet(tind).seg(i);
			sgn = tet(tind).sgn(i);
			msgn = 1;
			for (m = 0; m < basis::tet(log2p).em; ++m) {
				for(n = 0; n < NV; ++n)
					g.e(eind,m,n) += msgn*lf(n)(indx);
				msgn *= sgn;
				++indx;
			}
		}
	}
	
	/* FACE MODES */
	if (basis::tet(log2p).p > 2) {
		for(i = 0; i < 4; ++i){
			sgn = tet(tind).rot(i);
			find = tet(tind).tri(i);
			gindx = 0;
			msgn = 1;		
			for(m = 1; m <= basis::tet(log2p).em-1; ++m) {
				for(int j = 1; j <= basis::tet(log2p).em-m; ++j){
					for(n = 0; n < NV; ++n)
						g.f(find,gindx,n) += msgn*lf(n)(indx);
					++gindx;
					++indx;
				}
				msgn *= sgn;
			}
		}		
	}
	

	
	/* INTERIOR MODES */
	gindx = 0;
	for(m = 0; m < basis::tet(log2p).im; ++m) {
		for(n = 0; n < NV; ++n)
			g.i(tind,gindx,n) += lf(n)(indx);
		++gindx;
		++indx;
	}

   
   return;
}
#endif


#ifdef NODAL 
 void tet_hp::lftog(int tind, struct vefi g) {
   int i,m,n,indx,gindx,eind,find,sgn,msgn;
   
   /* VERTEX MODES */
   for (m = 0; m < 4; ++m) {
      gindx = tet(tind).pnt(m);
      for (n = 0; n < NV; ++n)
         g.v(gindx,n) += lf(n)(m);
   }

	/* EDGE MODES */
	if (basis::tet(log2p).p > 1) {
		indx = 4;
		for(i = 0; i < 6; ++i) {
			eind = tet(tind).seg(i);
			sgn = tet(tind).sgn(i);
			if(sgn > 0){
				for (m = 0; m < basis::tet(log2p).em; ++m) {
					for(n = 0; n < NV; ++n)
						g.e(eind,m,n) += lf(n)(indx);
					++indx;
				}
			}
			if(sgn < 0){
				for (m = 0; m < basis::tet(log2p).em; ++m) {
					for(n = 0; n < NV; ++n)
						g.e(eind,basis::tet(log2p).em-m-1,n) += lf(n)(indx);
					++indx;
				}
			}
		}
	}
	
	/* FACE MODES */
	if (basis::tet(log2p).p > 2) {
		for(i = 0; i < 4; ++i){
			sgn = tet(tind).rot(i);
			find = tet(tind).tri(i);
			gindx = 0;
			msgn = 1;		
			for(m = 1; m <= basis::tet(log2p).em-1; ++m) {
				for(int j = 1; j <= basis::tet(log2p).em-m; ++j){
					for(n = 0; n < NV; ++n)
						g.f(find,gindx,n) += msgn*lf(n)(indx);
					++gindx;
					++indx;
				}
				msgn *= sgn;
			}
		}		
	}
	

	
	/* INTERIOR MODES */
	gindx = 0;
	for(m = 0; m < basis::tet(log2p).im; ++m) {
		for(n = 0; n < NV; ++n)
			g.i(tind,gindx,n) += lf(n)(indx);
		++gindx;
		++indx;
	}

   
   return;
}

#endif

