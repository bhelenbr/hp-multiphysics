/*
 *  insert.cpp
 *  mblock
 *
 *  Created by helenbrk on Fri Aug 31 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#include "mesh.h"
#include "boundary.h"
#include <math.h>
#include <float.h>   
#include <assert.h>
#include <stdlib.h>

int mesh::insert(const TinyVector<FLT,ND> &x) {
   int n,tind,vnear,err;
   
   for(n=0;n<ND;++n)
      vrtx(nvrtx)(n) = x(n);
   
   qtree.addpt(nvrtx);
   qtree.nearpt(nvrtx,vnear);

   /* FIND TRIANGLE CONTAINING POINT */      
   tind = findtri(x,vnear);
   if (tind < 0) {
      std::cerr << "couldn't find triangle for point: " << x(0) << ' ' << x(1) << " vnear: " << vnear << std::endl;
      std::cerr << "maxsrch: " << maxsrch << "vtri: " << vd(vnear).tri << std::endl;
      output("error",ftype::grid);
      exit(1);
   }    
   if (nvrtx >= maxvst) {
      *log << "need to use larger growth factor: too many vertices" << std::endl;
      exit(1);
   }
   err = insert(nvrtx,tind);
   nvrtx += 1 -err;
   
   return(err);
}

int mesh::insert(int vnum, int tnum) {
   int ntdel, nskeep, nsdel;
   int i,j,tind,tin,tnxt,v0,sstart,dir,rstrt;
   int sind,sind1,snum;
   
   /* SEARCH SURROUNDING FOR NONDELAUNEY TRIANGLES */
   /* SIDES ON BOUNDARY OF HOLE (SKEEP) */
   /* SIDES IN INTERIOR OF HOLE (SDEL) (STARTS AT HALFWAY THROUGH LIST) */
   ntdel = 0;
   nskeep = 0;
   nsdel = 0;

   i2wk_lst1(ntdel++) = tnum;
   td(tnum).info |= TSRCH;
   
   tind = tnum;
   v0 = td(tnum).vrtx(0);
   /* FIRST SIDE DONE SEPARATELY SO KNOW WHEN TO STOP */
   snum = 2;
   tnxt = td(tind).tri(snum);
   sind = td(tind).side(snum);
   dir = (1+td(tind).sign(snum))/2;
   sstart = sind;
   rstrt = 1;
   
   if (tnxt < 0) {
      i2wk_lst2(nskeep++) = sind;
      i2wk_lst2(nskeep++) = dir;
      v0 = sd(sind).vrtx(dir);
      snum = (snum+1)%3;
   }
   else if (incircle(tnxt,vrtx(vnum)) > 0.0) {
      if (dir > 0) i2wk_lst3(nsdel++) = sind;
      rstrt = 2;
      i2wk_lst1(ntdel++) = tnxt;
      td(tnxt).info |= TSRCH;
      tind = tnxt;
      for(snum=0;snum<3;++snum) {
         if (td(tind).vrtx(snum) == v0) break;
      }
      snum = (snum+2)%3;
   }
   else {
      i2wk_lst2(nskeep++) = sind;
      i2wk_lst2(nskeep++) = dir;
      v0 = sd(sind).vrtx(dir);
      snum = (snum+1)%3;
   }
   tnxt = td(tind).tri(snum);
   sind = td(tind).side(snum);
   dir = (1+td(tind).sign(snum))/2;   
   
   /* GO COUNTER-CLOCKWISE AROUND VERTICES OF HOLE */
   /* IF START SIDE IS IN THE INTERIOR MUST HIT IT TWICE */
   for(j=0;j<rstrt;++j) {
      do  {
         if (tnxt < 0) {
            i2wk_lst2(nskeep++) = sind;
            i2wk_lst2(nskeep++) = dir;
            v0 = sd(sind).vrtx(dir);
            snum = (snum+1)%3;
         }
         else if (td(tnxt).info&TSRCH) {
            if (dir > 0) i2wk_lst3(nsdel++) = sind;
            tind = tnxt;
            for(snum=0;snum<3;++snum) {
               if (td(tind).vrtx(snum) == v0) break;
            }
            snum = (snum+2)%3;
         }
         else if (incircle(tnxt,vrtx(vnum)) > 0.0) {
            if (dir > 0) i2wk_lst3(nsdel++) = sind;
            i2wk_lst1(ntdel++) = tnxt;
            td(tnxt).info |= TSRCH;
            tind = tnxt;
            for(snum=0;snum<3;++snum) {
               if (td(tind).vrtx(snum) == v0) break;
            }
            snum = (snum+2)%3;
         }
         else {
            i2wk_lst2(nskeep++) = sind;
            i2wk_lst2(nskeep++) = dir;
            v0 = sd(sind).vrtx(dir);
            snum = (snum+1)%3;
         }
         tnxt = td(tind).tri(snum);
         sind = td(tind).side(snum);
         dir = (1+td(tind).sign(snum))/2;   
      } while (sind != sstart);
   }
      
   /* RESET TSRCH FLAGS */
   for(i=0;i<ntdel;++i)
      td(i2wk_lst1(i)).info &= ~TSRCH;
      
   nskeep = nskeep >> 1;
         
   /*	CHECK THAT WE AREN'T INSERTING POINT VERY CLOSE TO BOUNDARY */
   for(i=0;i<nskeep;++i) {
      sind = i2wk_lst2(i);
      if(fabs(minangle(vnum, sd(sind).vrtx(0) , sd(sind).vrtx(1))) < 1.0e-3*M_PI/180.0) {
         *log << "#Warning: inserting too close to boundary" << std::endl;
         return(1);
      }
   }
   
   /* APPEND NEW INDICES TO LIST OF INDICES TO USE FOR NEW SIDES & TRIS*/
   for(i=nsdel;i<nskeep;++i)
      i2wk_lst3(i) = nside +(i-nsdel);
   
   for(i=ntdel;i<nskeep;++i) 
      i2wk_lst1(i) = ntri +(i-ntdel);
   
   /* PERIODIC TRIANGLE */
   i2wk_lst1(nskeep) = i2wk_lst1(0);

   ntri += 2;
   nside += 3;
   
   if (ntri > maxvst || nside > maxvst) {
      *log << "need to use bigger growth factor: too many sides/tris" << nside << ntri << std::endl;
      exit(1);
   }
   
   for(i=0;i<nskeep;++i) {
      tind = i2wk_lst1(i);
      td(tind).info &= TTOUC;

      tnxt = i2wk_lst1(i+1);
      sind = i2wk_lst2(i<<1);
      dir = i2wk_lst2((i<<1) +1);
         
      /* CREATE NEW INFO */
      v0 = sd(sind).vrtx(1-dir);
      td(tind).vrtx(0) = v0;
      vd(v0).tri = tind;
      v0 = sd(sind).vrtx(dir);
      td(tind).vrtx(1) = v0;
      vd(v0).tri = tind;
      td(tind).vrtx(2) = vnum;
      vd(vnum).tri = tind;

      /* SIDE 2 INFO */      
      td(tind).side(2) = sind;
      td(tind).sign(2) = -1 +2*dir;
      
      sd(sind).tri(1-dir) = tind;
      tin = sd(sind).tri(dir);
      td(tind).tri(2) = tin;
      if (tin > -1) {
         for(j=0;j<3;++j) {
            if (td(tin).side(j) == sind) {
               td(tin).tri(j) = tind;
               break;
            }
         }
      }
         
      /* CREATE SIDE 0 */
      v0 = sd(sind).vrtx(dir);
      sind1 = i2wk_lst3(i);
      sd(sind1).tri(0) = tind;
      sd(sind1).vrtx(0) = v0;
      sd(sind1).vrtx(1) = vnum;
      td(sind1).info &= STOUC;


      td(tind).side(0) = sind1;
      td(tind).sign(0) = 1;
      td(tind).tri(0) = tnxt;

      sd(sind1).tri(1) = tnxt;
      td(tnxt).side(1) = sind1;
      td(tnxt).sign(1) = -1;
      td(tnxt).tri(1) = tind;
   }
   
   i2wk_lst1(-1) = ntdel;
   i2wk_lst2(-1) = nskeep;
   i2wk_lst3(-1) = nsdel;
      
   return(0);
}

void mesh::bdry_insert(int vnum, int sind, int endpt) {
   int ntdel, nsdel, nskeep;
   int sstart;
   int i,j,tin,tind,tnxt,v0,dir;
   int snum, sind1;
   
   /* ADD POINT TO QUADTREE */
   qtree.addpt(nvrtx);
   td(vnum).info &= VTOUC;
   
   /* SEARCH SURROUNDING FOR NONDELAUNEY TRIANGLES */
   /* SIDES ON BOUNDARY OF HOLE (SKEEP) */
   /* SIDES IN INTERIOR OF HOLE (SDEL) (STARTS AT HALFWAY THROUGH LIST) */
   ntdel = 0;
   nskeep = 0;
   nsdel = 0;

   tind = sd(sind).tri(0);
   i2wk_lst1(ntdel++) = tind;
   td(tind).info |= TSRCH;
   for(snum=0;snum<3;++snum)
      if (td(tind).side(snum) == sind) break;

   sstart = sind;
   v0 = sd(sind).vrtx(1);
   snum = (snum+1)%3;
   tnxt = td(tind).tri(snum);
   sind = td(tind).side(snum);
   dir = (1+td(tind).sign(snum))/2;
   
   /* GO COUNTER-CLOCKWISE AROUND VERTICES OF HOLE */
   do  {
      if (tnxt < 0) {
         i2wk_lst2(nskeep++) = sind;
         i2wk_lst2(nskeep++) = dir;
         v0 = sd(sind).vrtx(dir);
         snum = (snum+1)%3;
      }
      else if (td(tnxt).info&TSRCH) {
         if (dir > 0) i2wk_lst3(nsdel++) = sind;
         tind = tnxt;
         for(snum=0;snum<3;++snum) {
            if (td(tind).vrtx(snum) == v0) break;
         }
         snum = (snum+2)%3;
      }
      else if (incircle(tnxt,vrtx(vnum)) > 0.0) {
         if (dir > 0) i2wk_lst3(nsdel++) = sind;
         i2wk_lst1(ntdel++) = tnxt;
         td(tnxt).info |= TSRCH;
         tind = tnxt;
         for(snum=0;snum<3;++snum) {
            if (td(tind).vrtx(snum) == v0) break;
         }
         snum = (snum+2)%3;
      }
      else {
         i2wk_lst2(nskeep++) = sind;
         i2wk_lst2(nskeep++) = dir;
         v0 = sd(sind).vrtx(dir);
         snum = (snum+1)%3;
      }
      tnxt = td(tind).tri(snum);
      sind = td(tind).side(snum);
      dir = (1+td(tind).sign(snum))/2;   
   } while (sind != sstart);

   /* RESET TSRCH FLAGS */
   for(i=0;i<ntdel;++i)
      td(i2wk_lst1(i)).info &= ~TSRCH;

   /* ALTER OLD BOUNDARY SIDE & CREATE NEW SIDE */
   sd(nside).vrtx(endpt) = vnum;
   sd(nside).vrtx(1-endpt) = sd(sind).vrtx(1-endpt);
   sd(sind).vrtx(1-endpt) = vnum;
   td(sind).info &= STOUC;
   td(nside).info &= STOUC;
   
   /* ADD NEW SIDE TO BOUNDARY GROUP */
   /* NEED TO REORDER WHEN FINISHED */
   i = getbdrynum(sd(sind).tri(1));
   sd(nside).tri(1) = trinumatbdry(i,sbdry(i)->nel);
   sbdry(i)->el(sbdry(i)->nel++) = nside;
   ++nside;
   
   nskeep = nskeep>>1;
      
   /* APPEND NEW INDICES TO LIST OF INDICES TO USE FOR NEW SIDES & TRIS*/
   for(i=nsdel;i<nskeep;++i)
      i2wk_lst3(i) = nside +(i-nsdel);
   
   for(i=ntdel;i<nskeep;++i) 
      i2wk_lst1(i) = ntri +(i-ntdel);
         
   ++ntri;
   ++nside;  
   
   if (ntri > maxvst || nside > maxvst) {
      *log << "need to use bigger growth factor: too many sides/tris:" << nside << ntri << std::endl;
      exit(1);
   }
   
   /* FIX FIRST AND LAST TRIANGLE BOUNDARY SIDE */
   if (endpt) {
      tind = i2wk_lst1(0);
      td(tind).side(1) = sind;
      td(tind).sign(1) = 1;
      td(tind).tri(1) = sd(sind).tri(1);
      sd(sind).tri(0) = tind;
      
      tind = i2wk_lst1(nskeep-1);
      td(tind).side(0) = nside-2;
      td(tind).sign(0) = 1;
      td(tind).tri(0) = sd(nside-2).tri(1);
      sd(nside-2).tri(0) = tind;
 
   }
   else {
      tind = i2wk_lst1(0);
      td(tind).side(1) = nside-2;
      td(tind).sign(1) = 1;
      td(tind).tri(1) = sd(nside-2).tri(1);  
      sd(nside-2).tri(0) = tind;
      
      tind = i2wk_lst1(nskeep-1);
      td(tind).side(0) = sind;
      td(tind).sign(0) = 1;
      td(tind).tri(0) = sd(sind).tri(1);
      sd(sind).tri(0) = tind;
   }
   
   
   for(i=0;i<nskeep-1;++i) {
      tind = i2wk_lst1(i);
      td(tind).info &= TTOUC;

      tnxt = i2wk_lst1(i+1);
      sind = i2wk_lst2(i<<1);
      dir = i2wk_lst2((i<<1) +1);
         
      /* CREATE NEW INFO */
      v0 = sd(sind).vrtx(1-dir);
      td(tind).vrtx(0) = v0;
      vd(v0).tri = tind;
      v0 = sd(sind).vrtx(dir);
      td(tind).vrtx(1) = v0;
      vd(v0).tri = tind;
      td(tind).vrtx(2) = vnum;
      vd(vnum).tri = tind;

      /* SIDE 2 INFO */      
      td(tind).side(2) = sind;
      td(tind).sign(2) = -1 +2*dir;
      
      sd(sind).tri(1-dir) = tind;
      tin = sd(sind).tri(dir);
      td(tind).tri(2) = tin;
      if (tin > -1) {
         for(j=0;j<3;++j) {
            if (td(tin).side(j) == sind) {
               td(tin).tri(j) = tind;
               break;
            }
         }
      }
         
      /* CREATE SIDE 0 */
      v0 = sd(sind).vrtx(dir);
      sind1 = i2wk_lst3(i);
      sd(sind1).tri(0) = tind;
      sd(sind1).vrtx(0) = v0;
      sd(sind1).vrtx(1) = vnum;
      td(sind1).info &= STOUC;


      td(tind).side(0) = sind1;
      td(tind).sign(0) = 1;
      td(tind).tri(0) = tnxt;

      sd(sind1).tri(1) = tnxt;
      td(tnxt).side(1) = sind1;
      td(tnxt).sign(1) = -1;
      td(tnxt).tri(1) = tind;
   }
   
   /* LAST TRIANGLE */
   i = nskeep-1;
   tind = i2wk_lst1(i);
   td(tind).info &= TTOUC;

   sind = i2wk_lst2(i<<1);
   dir = i2wk_lst2((i<<1) +1);
      
   /* CREATE NEW INFO */
   v0 = sd(sind).vrtx(1-dir);
   td(tind).vrtx(0) = v0;
   vd(v0).tri = tind;
   v0 = sd(sind).vrtx(dir);
   td(tind).vrtx(1) = v0;
   vd(v0).tri = tind;
   td(tind).vrtx(2) = vnum;
   vd(vnum).tri = tind;

   /* SIDE 2 INFO */      
   td(tind).side(2) = sind;
   td(tind).sign(2) = -1 +2*dir;
   
   sd(sind).tri(1-dir) = tind;
   tin = sd(sind).tri(dir);
   td(tind).tri(2) = tin;
   if (tin > -1) {
      for(j=0;j<3;++j) {
         if (td(tin).side(j) == sind) {
            td(tin).tri(j) = tind;
            break;
         }
      }
   }
   
   return;
}

int mesh::findtri(const TinyVector<FLT,ND> x, int vnear) const {
   int i,j,vn,dir,stoptri,tin,tind;
   int ntdel;
   int tclose,nsurround;
   FLT minclosest,closest;
   
#if (((-1)^TSRCH)&TSRCH)
#define NOT 
#else
#define NOT !
#endif
   
   /* HERE WE USE i1wk & i2wk THIS MUST BE -1 BEFORE USING */
   tind = vd(vnear).tri;
   stoptri = tind;
   dir = 1;
   ntdel = 0;
   do {
      if (intri(tind,x) < area(tind)*10.*EPSILON) goto FOUND;
      i1wk(tind) ^= TSRCH;
      i2wk(ntdel++) = tind;
   
      for(vn=0;vn<3;++vn) 
         if (td(tind).vrtx(vn) == vnear) break;
      
      tind = td(tind).tri((vn +dir)%3);
      if (tind < 0) {
         if (dir > 1) break;
         /* REVERSE DIRECTION AND GO BACK TO START */
         ++dir;
         tind = vd(vnear).tri;
         for(vn=0;vn<3;++vn) 
            if (td(tind).vrtx(vn) == vnear) break;

         tind = td(tind).tri((vn +dir)%3);
         if (tind < 0) break;
         stoptri = -1;
      }
   } while(tind != stoptri); 
   
   nsurround = ntdel;
      
   /* DIDN'T FIND TRIANGLE */
   /* NEED TO SEARCH SURROUNDING TRIANGLES */
   for(i=0;i<ntdel;++i) {
      tin = i2wk(i);
      for(j=0;j<3;++j) {
         tind = td(tin).tri(j);
         if (tind < 0) continue;
         if (NOT(i1wk(tind)&TSRCH)) continue;
         i1wk(tind) ^= TSRCH;
         i2wk(ntdel++) = tind;         
         if (intri(tind,x) < area(tind)*10.*EPSILON) goto FOUND;
      }
      if (ntdel >= maxsrch-4) break;
   }
//   std::cerr << "couldn't find tri for point " << x[0] << ' ' << x[1] << ' ' << vnear << std::endl;
   tind = i2wk(0);
   minclosest = intri(tind,x)/area(tind);
   tclose = tind;
   for (i=1;i<nsurround;++i) {
      tind = i2wk(i);
      if ((closest = intri(tind,x)/area(tind)) < minclosest) {
         minclosest = closest;
         tclose = tind;
      }
   }
   intri(tclose,x);
   tind = -tclose;
      
FOUND:
   /* RESET INTWKW1 TO -1 */
   for(i=0;i<ntdel;++i) {
      i1wk(i2wk(i)) = -1;
   }
 
   return(tind);
}

