#include "mesh.h"
#include <utilities.h>
#include<assert.h>

void mesh::copy(const mesh& tgt) {
   int i,n;
      
   if (!initialized) {
      allocate_duplicate(1.0,tgt);
   }
   else {
      /* CHECK IF BIG ENOUGH */
      if (tgt.nside > maxvst) {
         *sim::log << "mesh is too big to copy" << std::endl;
         exit(1);
      }
   }
   
   /* COPY VERTEX INFO OVER */
   nvrtx = tgt.nvrtx;
   for(i=0;i<nvrtx;++i)
      for(n=0;n<ND;++n)
         vrtx(i)(n) = tgt.vrtx(i)(n);

   for(i=0;i<nvrtx;++i) {
      vlngth(i) = tgt.vlngth(i);
      vd(i).info = tgt.vd(i).info;
      vd(i).nnbor = tgt.vd(i).nnbor;
      vd(i).tri = tgt.vd(i).tri;
   }
      
   /* COPY VERTEX BOUNDARY INFO */
   for(i=0;i<nvbd;++i)
      vbdry(i)->copy(*tgt.vbdry(i));
         
/* COPY SIDE INFORMATION */
   nside = tgt.nside;
   for(i=0;i<nside;++i) {
      for(n=0;n<2;++n)
         sd(i).vrtx(n) = tgt.sd(i).vrtx(n);
      for(n=0;n<2;++n)
         sd(i).tri(n) = tgt.sd(i).tri(n);
      sd(i).info = tgt.sd(i).info;
   }
      
   /* COPY SIDE BOUNDARY INFO */
   for(i=0;i<nsbd;++i)
      sbdry(i)->copy(*tgt.sbdry(i));
   
   /* COPY ELEMENT DATA */
   ntri = tgt.ntri;
   for(i=0;i<ntri;++i) {
      for(n=0;n<3;++n) {
         td(i).vrtx(n) = tgt.td(i).vrtx(n);
         td(i).tri(n) = tgt.td(i).tri(n);
         td(i).side(n) = tgt.td(i).side(n);
         td(i).sign(n) = tgt.td(i).sign(n);
      }
      td(i).info = tgt.td(i).info;
   }
      
   qtree.copy(tgt.qtree);
   qtree.change_vptr((FLT (*)[ND]) vrtx(0).data() );
   
   return;  
}

void mesh::append(const mesh &z) {
	int i,j,k,n,nel,vrt,sind,flip;
   int nvrtxold, nsideold,ntriold;
   int sind1,tind1,v1a,v1b;
   int sind2,tind2,v2a,v2b;
         
	for(i=0;i<z.nvrtx;++i) {
		for(n=0;n<ND;++n) 
			vrtx(i+nvrtx)(n) = z.vrtx(i)(n);
      qtree.addpt(i+nvrtx);
   }
      
   for(i=0;i<z.nvrtx;++i)
      vlngth(i+nvrtx) = z.vlngth(i);

   for(i=0;i<z.nvrtx;++i) {
      vd(i+nvrtx).nnbor = z.vd(i).nnbor;
      vd(i+nvrtx).tri = z.vd(i).tri +ntri;
   }
      
   /* MOVE BOUNDARY INFO */
   vbdry.resizeAndPreserve(nvbd+z.nvbd);
   for(i=0;i<z.nvbd;++i) {
      vbdry(nvbd) = z.vbdry(i)->create(*this);
      vbdry(nvbd)->alloc(4);
      vbdry(nvbd)->v0 = z.vbdry(i)->v0 +nvrtx;
      ++nvbd;
   }
	
	for(i=0;i<z.nside;++i) {
		for(n=0;n<2;++n) {
			sd(i+nside).vrtx(n) = z.sd(i).vrtx(n) +nvrtx;
         sd(i+nside).tri(n) = z.sd(i).tri(n) +ntri;
      }
   }
   
   /* MOVE BOUNDARY INFO */
   sbdry.resizeAndPreserve(nsbd+z.nsbd);
   for(i=0;i<z.nsbd;++i) {
		sbdry(nsbd++) = z.sbdry(i)->create(*this);
      sbdry(nsbd-1)->alloc(z.sbdry(i)->maxel);
      sbdry(nsbd-1)->nel = z.sbdry(i)->nel;
		for(j=0;j<z.sbdry(i)->nel;++j)
			sbdry(nsbd-1)->el(j) = z.sbdry(i)->el(j) +nside;
	}

   /* MOVE TRI INFO */
   for(i=0;i<z.ntri;++i) {
		for(n=0;n<3;++n) {
			td(i+ntri).vrtx(n) = z.td(i).vrtx(n) +nvrtx;
         td(i+ntri).side(n) = z.td(i).side(n)+nside;
         td(i+ntri).sign(n) = z.td(i).sign(n);
         td(i+ntri).tri(n) = z.td(i).tri(n)+ntri;
      }
   }
	
   nvrtxold = nvrtx;
   nsideold = nside;
   ntriold = ntri;
   
	nvrtx += z.nvrtx;
	nside += z.nside;
	ntri += z.ntri;		
   
   bdrylabel();  // CHANGES STRI / TTRI ON BOUNDARIES TO POINT TO GROUP/ELEMENT   
   
   /* KEEPS TRACK OF DELETED SIDES = -3, TOUCHED SIDES =-2, UNTOUCHED SIDES =-1 */
   for(i=nsideold;i<nside;++i)
      sd(i).info = -1;

   /* VINFO TO KEEP TRACK OF DELETED VERTICES (-3) */
   for(i=nvrtxold;i<nvrtx;++i)
      vd(i).info = -1;
   
   /* FIND MATCHING COMMUNICATION BOUNDARIES */
   for(i=0;i<nsbd -z.nsbd;++i) {
      if (!sbdry(i)->is_comm()) continue;
      
      for(j=0;j<z.nsbd;++j) {
         if (sbdry(i)->idnum == z.sbdry(j)->idnum) {
            nel = sbdry(i)->nel;
            
            /* FIRST POINT DONE REVERSE */
            sind1 = sbdry(i)->el(nel-1);
            tind1 = sd(sind1).tri(0);            
            v1b = sd(sind1).vrtx(1);
            sind2 = z.sbdry(j)->el(0);
            tind2 = z.sd(sind2).tri(0) +ntriold;
            v2a = z.sd(sind2).vrtx(0) +nvrtxold;
            vd(v1b).nnbor += z.vd(v2a-nvrtxold).nnbor-1;
            qtree.dltpt(v2a);
            vd(v2a).info = -3;
            do {
               for(vrt=0;vrt<3;++vrt) 
                  if (td(tind2).vrtx(vrt) == v2a) break; 
               assert(vrt != 3);
               
               td(tind2).vrtx(vrt) = v1b;
               sind = td(tind2).side((vrt +1)%3);
               flip = (1 +td(tind2).sign((vrt +1)%3))/2;
               assert(sd(sind).vrtx(flip) == v2a);
               sd(sind).vrtx(flip) = v1b;
               tind2 = td(tind2).tri((vrt +1)%3);
               
            } while(tind2 > 0); 
            
            for(k=0;k<nel;++k) {
               sind1 = sbdry(i)->el(nel-k-1);
               tind1 = sd(sind1).tri(0);
               v1a = sd(sind1).vrtx(0);
               sind2 = z.sbdry(j)->el(k);
               tind2 = z.sd(sind2).tri(0) +ntriold;
               v2b = z.sd(sind2).vrtx(1) +nvrtxold;
               sind2 += nsideold;               
               vd(v1a).nnbor += z.vd(v2b-nvrtxold).nnbor-2;

               for(vrt=0;vrt<3;++vrt) 
                  if (td(tind1).side(vrt) == sind1) break;
               assert(vrt < 3);
               td(tind1).tri(vrt) = tind2;
               sd(sind1).tri(1) = tind2;
               
               for(vrt=0;vrt<3;++vrt) 
                  if (td(tind2).vrtx(vrt) == v2b) break;
               assert(vrt < 3);
               td(tind2).side((vrt+1)%3) = sind1;
               td(tind2).sign((vrt+1)%3) = -1;
               td(tind2).tri((vrt+1)%3) = tind1;
               sd(sind2).info = -3;
               vd(v2b).info = -3;
               qtree.dltpt(v2b);
               
               for(;;) {
                  td(tind2).vrtx(vrt) = v1a;
                  sind = td(tind2).side((vrt +2)%3);
                  flip = (1 -td(tind2).sign((vrt +2)%3))/2;
                  assert(sd(sind).vrtx(flip) == v2b);
                  sd(sind).vrtx(flip) = v1a;
                  tind2 = td(tind2).tri((vrt +2)%3);
                  
                  if (tind2 < 0) break;
                  
                  for(vrt=0;vrt<3;++vrt) 
                     if (td(tind2).vrtx(vrt) == v2b) break; 
                  assert(vrt != 3);
               }
            }
            /* LAST VRTX ONLY SHARES ONE EDGE */
            vd(v1a).nnbor += 1;
            
            delete sbdry(i);
            for(k=i;k<nsbd-1;++k)
               sbdry(k) = sbdry(k+1);
            delete sbdry(nsbd-z.nsbd+j-1);
            for(k=nsbd-z.nsbd+j-1;k<nsbd-2;++k)
               sbdry(k) = sbdry(k+1);
            nsbd -= 2;
            break;
         }
      }
   }
   
      /* DELETE LEFTOVER VERTICES */
   /* VINFO > NVRTX STORES VRTX MOVEMENT HISTORY */
   for(i=0;i<nvrtx;++i) 
      if (vd(i).info == -3) 
         dltvrtx(i);
   
   /* FIX BOUNDARY CONDITION POINTERS */
   for(i=nvbd-z.nvbd;i<nvbd;++i)
      if (vbdry(i)->v0 >= nvrtx) 
         vbdry(i)->v0 = vd(vbdry(i)->v0).info;  
                        
   /* CLEAN UP SIDES */
   /* SINFO WILL END UP STORING -1 UNTOUCHED, -2 TOUCHED, or INITIAL INDEX OF UNTOUCHED SIDE */
   /* SINFO > NSIDE WILL STORE MOVEMENT HISTORY */
   for(i=nsideold;i<nside;++i) 
      if (sd(i).info == -3) dltsd(i);
   
   /* FIX BOUNDARY CONDITION POINTERS */
   for(i=nsbd-z.nsbd+1;i<nsbd;++i)
      for(j=0;j<sbdry(i)->nel;++j) 
         if (sbdry(i)->el(j) >= nside) 
            sbdry(i)->el(j) = sd(sbdry(i)->el(j)).info; 

   for (i=0;i<nsbd;++i) {
      sbdry(i)->reorder();
   }
   
  	bdrylabel();  // CHANGES STRI / TTRI ON BOUNDARIES TO POINT TO GROUP/ELEMENT               
      
   return;
}


