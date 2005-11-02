#include "mesh.h"
#include "boundary.h"
#include <utilities.h>
#include <assert.h>
#include <float.h>

void mesh::collapse(int sind, int delt) {
   int ntsrnd, nssrnd, nperim;
   int i,j,vn,vnear,prev,tind,tind1,sind1,stoptri,dir;
   int v0,v1,pt,sd1,sd2,sd3,t1,t2;

   /* i2wk_lst1 = triangles surrounding delete vertex */
   /* i2wk_lst2 = sides surrounding delete vertex */
   /* -1 indx lists how many */
   
   /* FIND TRIANGLES / SIDES SURROUNDING ENDPOINT */
   /* EXCLUDES TRIANGLES ADJACENT TO DELETED SIDE */
   /* EXCLUDES SIDES OF TRIANGLES ADJACENT TO DELETED SIDE */
   vnear = sd(sind).vrtx(delt);
   tind = vd(vnear).tri;
   if (tind != sd(sind).tri(0) && tind != sd(sind).tri(1))
      prev = 0;
   else 
      prev = 1;
   stoptri = tind;
   dir = 1;
   ntsrnd = 0;
   nssrnd = 0;
   nperim = 0;
   do {
      for(vn=0;vn<3;++vn) 
         if (td(tind).vrtx(vn) == vnear) break;
               
      tind1 = td(tind).tri((vn +dir)%3);
      if (tind1 < 0) {
         if (dir > 1) break;
         /* REVERSE DIRECTION AND GO BACK TO START */
         ++dir;
         tind1 = vd(vnear).tri;
         prev = 1;
         stoptri = -1;
      }
      
      if (tind1 != sd(sind).tri(0) && tind1 != sd(sind).tri(1)) {
         i2wk_lst1(ntsrnd++) = tind1;
         if (!prev) {
            i2wk_lst2(nssrnd++) = td(tind).side((vn +dir)%3);
         }
         prev = 0;
      }
      else {
         prev = 1;
      }

      tind = tind1;

   } while(tind != stoptri); 
   i2wk_lst1(-1) = ntsrnd;
   i2wk_lst2(-1) = nssrnd;
   
   /* UPDATE TVRTX & SVRTX */
   v0 = sd(sind).vrtx(1-delt);
   v1 = sd(sind).vrtx(delt);
   
   for(i=0;i<ntsrnd;++i) {
      tind = i2wk_lst1(i);
      td(tind).info |= TTOUC;
      for(j=0;j<3;++j) {
         if (td(tind).vrtx(j) == v1) {
            td(tind).vrtx(j) = v0;
            sd3 = td(tind).side((j+1)%3);
            td(sd3).info |= STOUC;
            pt = (1 +td(tind).sign((j+1)%3))/2;
            assert(sd(sd3).vrtx(pt) == v1 || sd(sd3).vrtx(pt) == v0);
            sd(sd3).vrtx(pt) = v0;
            sd3 = td(tind).side((j+2)%3);
            td(sd3).info |= STOUC;
            pt = (1 -td(tind).sign((j+2)%3))/2;
            assert(sd(sd3).vrtx(pt) == v1 || sd(sd3).vrtx(pt) == v0);
            sd(sd3).vrtx(pt) = v0;
            i2wk_lst3(nperim++) = td(tind).side(j);
            break;
         }
      }
   }

   /* MARK SIDE AS DELETED */ 
   td(sind).info |= SDLTE; 

   /* CLOSE THE GAP */
   for(i=0;i<2;++i) {
      tind = sd(sind).tri(i);
      if (tind < 0) continue;

      /* MARK TRI AS DELETED TO BE REMOVED LATER */
      td(tind).info |= TDLTE;  
      
      for(sd1=0;sd1<3;++sd1)
         if(td(tind).vrtx(sd1) == v1) break;
         
      assert(sd1 != 3);
               
      if (td(tind).side((sd1+1)%3) == sind)
         sd2 = (sd1+2)%3;
      else
         sd2 = (sd1+1)%3;
               
      sind1 = td(tind).side(sd2);
      if (sd(sind1).tri(1) < 0) {
         /* BOUNDARY SIDE SO DELETE OTHER ONE */
         /* ANOMALY OF TRIANGLE WITH 2 EDGES ON BOUNDARY */
         td(sind1).info |= STOUC;
         if (sd(sind1).vrtx(0) == v1) 
            sd(sind1).vrtx(0) = v0;
         if (sd(sind1).vrtx(1) == v1) 
            sd(sind1).vrtx(1) = v0;   
         sind1 = sd2;
         sd2 = sd1;
         sd1 = sind1;
         sind1 = td(tind).side(sd2);
      }
      td(sind1).info |= SDLTE;
      i2wk_lst2(nssrnd++) = sind1;
            
      t1 = td(tind).tri(sd1);
      t2 = td(tind).tri(sd2);

      /* UPDATE TTRI FOR T1 */  
      if (t1 > -1) {
         for(j=0;j<3;++j) {
            if(td(t1).tri(j) == tind) {
               td(t1).tri(j) = t2;
               break;
            }
         }
         for(j=0;j<3;++j)
            vd(td(tind).vrtx(j)).tri = t1;
      }
      /* UPDATE STRI FOR KEPT SIDE */
      sind1 = td(tind).side(sd1);
      pt = (1 -td(tind).sign(sd1))/2;
      sd(sind1).tri(pt) = t2;
      i2wk_lst3(nperim++) = sind1;

      /* UPDATE TTRI/TSIDE FOR T2 */
      if (t2 > -1) {
         for(j=0;j<3;++j) {
            if(td(t2).tri(j) == tind) {
               td(t2).tri(j) = t1;
               td(t2).side(j) = sind1;
               td(t2).sign(j) = td(tind).sign(sd1);
               break;
            }
         }
         for(j=0;j<3;++j)
            vd(td(tind).vrtx(j)).tri = t2;
      }
   }
   
   /* NEED TO REMOVE LEFTOVERS */
   qtree.dltpt(v1);
   td(v1).info |= VDLTE;
      
   /* SWAP AFFECTED SIDES */      
   swap(i2wk_lst2(-1),&i2wk_lst2(0));
   i2wk_lst2(-1) = nssrnd;
   i2wk_lst3(-1) = nperim;

   return;
}

/* DELETE UNREFERENCED TRIANGLE */
void mesh::dlttri(int tind) {
   int i,j,v0,t1,sind,flip;
   
   
   /* MOVE UP FROM BOTTOM UNTIL FIND UNDELETED TRI */
   while(td(ntri-1).info&TDLTE)
      --ntri;
   
   if (ntri <=  tind)  {
      ntri = tind;
      return;
   }
   
   --ntri;

   for (j=0;j<3;++j) {
      v0 = td(ntri).vrtx(j);
      td(tind).vrtx(j) = v0;
      vd(v0).tri = tind;
      td(tind).side(j) = td(ntri).side(j);
      td(tind).sign(j) = td(ntri).sign(j);
      sind = td(ntri).side(j);
      flip = (1 -td(ntri).sign(j))/2;
      sd(sind).tri(flip) = tind;
      
      t1 = td(ntri).tri(j);
      td(tind).tri(j) = t1;
      if (t1 > -1) {
         for(i=0;i<3;++i) {
            if(td(t1).tri(i) == ntri) {
               td(t1).tri(i) = tind;
               break;
            }
         }
      }
   }
   
   if (td(ntri).info&TTOUC) {
      td(tind).info = -2;
      updatetdata(tind);
   }
   else {
      td(tind).info = ntri;
      movetdata(ntri,tind);
   }
   /* THIS IS TO PREVENT ROLL BACK NEAR END */
   td(tind).info &= ~TDLTE;

   td(ntri).info = tind;
   
   return;
}

/* DELETE UNREFERENCED SIDE */
void mesh::dltsd(int sind) {
   int j,k,tind;

   /* MOVE UP FROM BOTTOM UNTIL FIND UNDELETED SIDE */
   while(td(nside-1).info&SDLTE)
      --nside;
   
   if (nside <= sind) {
      nside = sind;
      return;
   }
   
   /* DELETE SIDE */
   --nside;
   
   sd(sind).vrtx(0) = sd(nside).vrtx(0);
   sd(sind).vrtx(1) = sd(nside).vrtx(1);
   sd(sind).tri(0) = sd(nside).tri(0);
   sd(sind).tri(1) = sd(nside).tri(1);

   for(j=0;j<2;++j) {
      tind = sd(nside).tri(j);
      if (tind > -1) {
         for(k=0;k<3;++k)
            if (td(tind).side(k) == nside) break;
         td(tind).side(k) = sind;
      }
   }
   
   if (td(nside).info&STOUC) {
      sd(sind).info = -2;
      updatesdata(sind);
   }
   else {
      sd(sind).info = nside;
      movesdata(nside,sind);
   }
   sd(nside).info = sind;
   
   return;
}

/* DELETE UNREFERENCED VERTEX */
void mesh::dltvrtx(int v0) {
   int vn,stoptri,dir;
   int tind, sind, flip;
            
   /* MOVE UP FROM BOTTOM UNTIL FIND UNDELETED VERTEX */
   while(td(nvrtx-1).info&VDLTE) {
      --nvrtx;
   }
      
   if (nvrtx <= v0) {
      nvrtx = v0;
      return;
   }
   
   --nvrtx;
   
   vrtx(v0)(0) = vrtx(nvrtx)(0);
   vrtx(v0)(1) = vrtx(nvrtx)(1);
   vlngth(v0) = vlngth(nvrtx);
   vd(v0).nnbor = vd(nvrtx).nnbor;
   vd(v0).tri = vd(nvrtx).tri;

   qtree.movept(nvrtx,v0);
   
   tind = vd(nvrtx).tri;   
   stoptri = tind;
   dir = 1;
   do {
      for(vn=0;vn<3;++vn) 
         if (td(tind).vrtx(vn) == nvrtx) break; 
         
      assert(vn != 3);
      
      td(tind).vrtx(vn) = v0;
      sind = td(tind).side((vn +dir)%3);
      flip = (1 +(3-2*dir)*td(tind).sign((vn +dir)%3))/2;
      assert(sd(sind).vrtx(flip) == nvrtx);
      sd(sind).vrtx(flip) = v0;

      tind = td(tind).tri((vn +dir)%3);
      if (tind < 0) {
         if (dir > 1) break;
         /* REVERSE DIRECTION AND GO BACK TO START */
         ++dir;
         tind = vd(nvrtx).tri;
         for(vn=0;vn<3;++vn) {
            if (td(tind).vrtx(vn) == v0) {
               td(tind).vrtx(vn) = nvrtx;
               break;
            }
         }
         stoptri = -1;
      }
      
   } while(tind != stoptri); 
   
   if (td(nvrtx).info&VTOUC) {
      vd(v0).info = -2;
      updatevdata(v0);
   }
   else {
      vd(v0).info = nvrtx;
      movevdata(nvrtx,v0);
   }
   vd(nvrtx).info = v0;

   return;
}
