/*
 *  yabers.cpp
 *  mblock
 *
 *  Created by helenbrk on Fri Sep 14 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */

#include "mesh.h"
#include "boundary.h"
#include <utilities.h>
#include <assert.h>

/* THIS IS SUPPOSED TO DO THE REVERSE OF THE REBAY ROUTINE I HOPE */
/* THUS THE NAME YABER -> REBAY */

#define NO_DEBUG_ADAPT

#ifdef DEBUG_ADAPT
int adapt_count;
static std::string adapt_file;
#endif

extern int nlst; 

void mesh::yaber(FLT tolsize) {
   int i,j,tind,sind,sind1,v0,cnt,endpt,sum;
   FLT x,y,a,asum,dx,dy,l0,l1;
   int ntsrnd, nssrnd, nperim;
   int vn,vnear,prev,tind1,stoptri,dir;
   int v1;
   int snum;
   FLT sratio;
   TinyVector<int,3> badside;

   FLT minvl;

   /* TO ADJUST FOR INSCRIBED RADIUS */
   tolsize *= 6./sqrt(3.);
   
   /* SET UP FLTWK */
   nlst = 0;
   for(i=0;i<ntri;++i) {
      if (td(i).info&TDLTE) continue;
      minvl = vlngth(td(i).vrtx(0));
      minvl = MIN(minvl,vlngth(td(i).vrtx(1)));
      minvl = MIN(minvl,vlngth(td(i).vrtx(2)));
      fscr1(i) = minvl/inscribedradius(i);
      if (fscr1(i) > tolsize) putinlst(i);
   }
   
   /* MARK BOUNDARY VERTEX POINTS */
   /* THESE SHOULD NOT BE DELETED */
   for(i=0;i<nsbd;++i) {
      for(j=0;j<sbdry(i)->nel;++j) {
         sind = sbdry(i)->el(j);
         v0 = sd(sind).vrtx(0);
         td(v0).info |= VSPEC;
      }      
   }
   
   cnt = 0;
   /* BEGIN COARSENING ALGORITHM */
   while (nlst > 0) {
      for(i=nlst-1;i>=0;--i) {  // START WITH LARGEST TGT TO ACTUAL RATIO
         for(j=0;j<3;++j) {
            tind = td(sd(i).info).tri(j);
            if (tind < 0 || vd(tind).info == -1)  goto TFOUND;
         }
      }
      
      TFOUND:
      tind = sd(i).info;
      tkoutlst(tind);
      
      /* FIND SIDE ON TRIANGLE WITH LARGEST TARGET TO LENGTH RATIO */
      minvl = 0.0;
      for(j=0;j<3;++j) {
         sind = td(tind).side(j);
         if (sd(sind).tri(1) < 0) {
            badside(j) = false;
            continue;
         }
         sratio = MIN(vlngth(sd(sind).vrtx(0)),vlngth(sd(sind).vrtx(1)))/distance(sd(sind).vrtx(0),sd(sind).vrtx(1));
         if (sratio > tolsize*sqrt(3.)/6.) {
            badside(j) = true;
         }
         else {
            badside(j) = false;
         }
         if (sratio > minvl) {
            sind1 = sind;
            snum = j;
            minvl =sratio;
         }
      }

      sind = sind1;

      /* REMOVE TRIANGLES THAT WILL BE DELETED */
      tind = sd(sind).tri(0);
      if (vd(tind).info > -1) tkoutlst(tind);
      tind = sd(sind).tri(1);
      if (vd(tind).info > -1) tkoutlst(tind);

      /* DON'T DELETE BOUNDARY POINT */
      sum = (td(sd(sind).vrtx(0)).info&VSPEC) +(td(sd(sind).vrtx(1)).info&VSPEC);
      if (sum > 0) {
         if (sum > VSPEC) {
            std::cout << "Big oops near point" << vrtx(sd(sind).vrtx(0)) << std::endl;
            continue;
         }
         if (td(sd(sind).vrtx(0)).info&VSPEC) endpt = 1;
         else endpt = 0;
      }
      else {
         /* TRY TO FIND DIRECTION OF ACCEPTED TRIS AND DELETE AWAY FROM THAT DIRECTION */
//         if (vd(td(tind).tri((snum+1)%3)).info == -1  && vd(td(tind).tri((snum+2)%3)).info > -1) {
//            endpt = (1-td(tind).sign(snum))/2;
//         }
//         else if (vd(td(tind).tri((snum+1)%3)).info > -1  && vd(td(tind).tri((snum+2)%3)).info == -1) {
//            endpt = (1+td(tind).sign(snum))/2;
//         }
//         else {
         {
            /* KEEP POINT WHICH IS CLOSEST TO CENTER OF AREA */
            /* THIS WAY WORKS BEST BUT COSTS MORE */
            /* TEMPORARY:: NEED TO ELIMINATE DUPLICATION OF WORK BETWEEN THIS & COLLAPSE */
            x = 0.0;
            y = 0.0;
            asum = 0.0;
            for(endpt=0;endpt<2;++endpt) {
               vnear = sd(sind).vrtx(endpt);
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
               
               tind = sd(sind).tri(endpt);
               if (tind > -1) {
                  a = area(tind);
                  asum += a;
                  for(vn=0;vn<3;++vn) {
                     x += a*vrtx(td(tind).vrtx(vn))(0);
                     y += a*vrtx(td(tind).vrtx(vn))(1);
                  }
               }            
               for(j=0;j<ntsrnd;++j) {
                  tind = i2wk_lst1(j);
                  a = area(tind);
                  asum += a;
                  for(vn=0;vn<3;++vn) {
                     x += a*vrtx(td(tind).vrtx(vn))(0);
                     y += a*vrtx(td(tind).vrtx(vn))(1);
                  }
               }
            }

            asum = 1./(3.*asum);
            x = x*asum;
            y = y*asum;
            v0 = sd(sind).vrtx(0);
            v1 = sd(sind).vrtx(1);
            dx = vrtx(v0)(0) -x;
            dy = vrtx(v0)(1) -y;
            l0 = dx*dx +dy*dy;   
            dx = vrtx(v1)(0) -x;
            dy = vrtx(v1)(1) -y;
            l1 = dx*dx +dy*dy;
            
            endpt = (l0 > l1 ? 0 : 1);
            
         }
      }

#ifdef DEBUG_ADAPT
      std::cout << "collapsing interior" << cnt << ' ' << adapt_count << ' ' << sind << " endpt " << endpt << std::endl;
#endif
      
      /* COLLAPSE EDGE */
      collapse(sind,endpt);
      ++cnt;
      
      /* RECLASSIFY AFFECTED TRIANGLES */
      for(i=0;i<i2wk_lst1(-1);++i) {
         tind = i2wk_lst1(i);
         if (vd(tind).info > -1) tkoutlst(tind);
         minvl = vlngth(td(tind).vrtx(0));
         minvl = MIN(minvl,vlngth(td(tind).vrtx(1)));
         minvl = MIN(minvl,vlngth(td(tind).vrtx(2)));
         fscr1(tind) = minvl/inscribedradius(tind);
         if (fscr1(tind) > tolsize) putinlst(tind);
      }
      
#ifdef DEBUG_ADAPT
      std::ostringstream nstr;
      nstr << adapt_count++ << std::flush;
      adapt_file = "adapt" +nstr.str();
      nstr.str("");
      output(adapt_file,grid);
#endif
   }

   *sim::log << "#Yaber finished: " << cnt << " sides coarsened" << std::endl;

   return;
}

void mesh::checkintegrity() {
   int i,j,sind,dir;
   
   for(i=0;i<ntri;++i) {
      if (td(i).info < 0) continue;
      
      if (area(i) < 0.0) *sim::log << "negative area" << i << std::endl;
      
      for(j=0;j<3;++j) {
         sind = td(i).side(j);
         dir = -(td(i).sign(j) -1)/2;
         
         if (sd(sind).info == -3) {
            *sim::log << "references deleted side" <<  i << sind << std::endl;
            for(i=0;i<nside;++i)
               sd(i).info += 2;
            output("error");
            exit(1);
         }

         if (sd(sind).vrtx(dir) != td(i).vrtx((j+1)%3) && sd(sind).vrtx(1-dir) != td(i).vrtx((j+2)%3)) {
            *sim::log << "failed vrtx check tind" << i << "sind" << sind << std::endl;
            for(i=0;i<nside;++i)
               sd(i).info += 2;
            output("error"); 
            exit(1);
         }    
         
         if (sd(sind).tri(dir) != i) {
            *sim::log << "failed side check tind" << i << "sind" << sind << std::endl;
            for(i=0;i<nside;++i)
               sd(i).info += 2;
            output("error"); 
            exit(1);
         }
         
         if (td(i).tri(j) != sd(sind).tri(1-dir)) {
            *sim::log << "failed ttri check tind" << i << "sind" << sind << std::endl;
            for(i=0;i<nside;++i)
               sd(i).info += 2;
            output("error"); 
            exit(1);
         }
         
         if (td(i).tri(j) > 0) {
            if(td(td(i).tri(j)).info < 0) {
               *sim::log << "references deleted tri" << std::endl;
               for(i=0;i<nside;++i)
                  sd(i).info += 2;
               output("error"); 
               exit(1);
            }
         }
      }
   }
   
   return;
}

void mesh::bdry_yaber(FLT tolsize) {
   int sind,endpt,v0,v1,count;
   int el, nel, pel, next, sindprev, sindnext, saffect;

   /* COARSEN FIRST BOUNDARIES */
   count = 0;
   for(int bnum=0;bnum<nsbd;++bnum) {
      if (!sbdry(bnum)->is_frst()) {
         sbdry(bnum)->master_slave_prepare();
         continue;
      }
      
      nlst = 0;
      for(int indx=0;indx<sbdry(bnum)->nel;++indx) {
         sind = sbdry(bnum)->el(indx);
         if (td(sind).info&SDLTE) continue;
         fscr1(sind) = MIN(vlngth(sd(sind).vrtx(0)),vlngth(sd(sind).vrtx(1)))/distance(sd(sind).vrtx(0),sd(sind).vrtx(1));
         if (fscr1(sind) > tolsize) {
            putinlst(sind);
         }
      }
      
      /* SKIP FIRST SPOT SO CAN SEND LENGTH FIRST */
      sbdry(bnum)->sndsize() = 1;
      while (nlst > 0) {
         // START WITH LARGEST SIDE LENGTH RATIO
         sind = sd(nlst-1).info;
         el = getbdryel(sd(sind).tri(1));
         
         /* FIND ADJACENT NON-DELETED SIDES */
         nel = -1;
         for(next = el+1;next<sbdry(bnum)->nel;++next) {
            if(!(td(sbdry(bnum)->el(next)).info&SDLTE)) {
               nel = next;
               sindnext = sbdry(bnum)->el(nel);
               break;
            }
         }
         pel = -1;
         for(next = el-1;next>=0;--next) {
            if(!(td(sbdry(bnum)->el(next)).info&SDLTE)) {
               pel = next;
               sindprev = sbdry(bnum)->el(pel);
               break;
            }
         }
         
         /* PICK ENDPT */
         v0 = sd(sind).vrtx(0);
         v1 = sd(sind).vrtx(1);

         if (td(v0).info&VSPEC && td(v1).info&VSPEC) {
            tkoutlst(sind);
            continue;
         }
         else if (td(v0).info&VSPEC) {
            endpt = 1;
            saffect = sindnext;
         }
         else if (td(v1).info&VSPEC) {
            endpt = 0;
            saffect = sindprev;
         }
         else {
            /* PICK MORE ARPROPRIATE VERTEX TO DELETE */
            if (fscr1(sindprev) > fscr1(sindnext)) {
               endpt = 0;
               saffect = sindprev;
            }
            else {
               endpt = 1;
               saffect = sindnext;
            }
         }

#ifdef DEBUG_ADAPT
      std::cout << "collapsing boundary"  << adapt_count << ' ' << sind << " endpt " << endpt << std::endl;
#endif
         collapse(sind,endpt);
         tkoutlst(sind);
         
         /* STORE ADAPTATION INFO FOR COMMUNICATION BOUNDARIES */
         sbdry(bnum)->isndbuf(sbdry(bnum)->sndsize()++) = el;
         sbdry(bnum)->isndbuf(sbdry(bnum)->sndsize()++) = endpt;
         
         /* UPDATE AFFECTED SIDE */
         if (vd(saffect).info > -1) tkoutlst(saffect);
         fscr1(saffect) = MIN(vlngth(sd(saffect).vrtx(0)),vlngth(sd(saffect).vrtx(1)))/distance(sd(saffect).vrtx(0),sd(saffect).vrtx(1));
         if (fscr1(saffect) > tolsize) putinlst(saffect);
         ++count;
#ifdef DEBUG_ADAPT
         number_str(adapt_file,"adapt",adapt_count++,5);
         output(adapt_file,grid);
#endif            
      }
      sbdry(bnum)->isndbuf(0) = sbdry(bnum)->sndsize();
   }
   *sim::log << "#Boundary coarsening finished, " << count << " sides coarsened" << std::endl;
   return;
}

void mesh::bdry_yaber1() {
   int i,el,endpt,sind,sndsize;
   
   for(int bnum=0;bnum<nsbd;++bnum) {
      
      if (sbdry(bnum)->is_frst() || !sbdry(bnum)->is_comm()) continue;
      
      sbdry(bnum)->master_slave_wait();
      
      sndsize = sbdry(bnum)->ircvbuf(0,0);
      
      for(i=1;i<sndsize;i+=2) {
         el = sbdry(bnum)->nel -1 -sbdry(bnum)->ircvbuf(0,i);
         endpt = 1 -sbdry(bnum)->ircvbuf(0,i+1);
         sind = sbdry(bnum)->el(el);
#ifdef DEBUG_ADAPT
         std::cout << "collapsing boundary" << adapt_count << ' ' << sind << " endpt " << endpt << std::endl;
#endif
         collapse(sind,endpt);
#ifdef DEBUG_ADAPT
         number_str(adapt_file,"adapt",adapt_count++,5);
         output(adapt_file,grid);
#endif
      }
   }
   return;
}


void mesh::checkintwk() const {
   int i;
   
   for(i=0;i<maxvst;++i)
      if (i1wk(i) != -1) *sim::log << "failed intwk1 check" << i << i1wk(i) << std::endl;
   
   return;
}
