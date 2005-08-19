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

#ifdef DEBUG_ADAPT
extern int adapt_count;
extern char adapt_file[100];
#endif

extern int nlst; 

void mesh::yaber(FLT tolsize) {
   int i,j,tind,sind,v0,cnt,endpt;
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
      if (td(td(tind).vrtx(j)).info&VSPEC) {
         tkoutlst(tind);
         continue;
      }

      sind = td(tind).side((j+1)%3);
      endpt = (1+td(tind).sign((j+1)%3))/2;
      
      /* REMOVE TRIANGLES THAT WILL BE DELETED */
      tind = sd(sind).tri(0);
      if (vd(tind).info > -1) tkoutlst(tind);
      tind = sd(sind).tri(1);
      if (vd(tind).info > -1) tkoutlst(tind);
   
      /* COLLAPSE EDGE */
      collapse(sind,endpt);
      ++cnt;
      
      /* RECLASSIFY AFFECTED TRIANGLES */
      for(i=0;i<i2wk_lst1(-1);++i) {
         tind = i2wk_lst1(i);
         if (vd(tind).info > -1) {
            tkoutlst(tind);
            minvl = vlngth(td(tind).vrtx(0));
            minvl = MIN(minvl,vlngth(td(tind).vrtx(1)));
            minvl = MIN(minvl,vlngth(td(tind).vrtx(2)));
            fscr1(tind) = minvl/inscribedradius(tind);
            if (fscr1(tind) > tolsize) putinlst(tind);
         }
      }
      
#ifdef DEBUG_ADAPT
      std::cout << "collapsed " << sind << " endpt " << endpt << std::endl;
      number_str(adapt_file,"adapt",adapt_count++,5);
      output(adapt_file,ftype::grid);
#endif
   }

   *log << "#Yaber finished: " << cnt << " sides coarsened" << std::endl;

   return;
}

void mesh::checkintegrity() {
   int i,j,sind,dir;
   
   for(i=0;i<ntri;++i) {
      if (td(i).info < 0) continue;
      
      if (area(i) < 0.0) *log << "negative area" << i << std::endl;
      
      for(j=0;j<3;++j) {
         sind = td(i).side(j);
         dir = -(td(i).sign(j) -1)/2;
         
         if (sd(sind).info == -3) {
            *log << "references deleted side" <<  i << sind << std::endl;
            for(i=0;i<nside;++i)
               sd(i).info += 2;
            output("error");
            exit(1);
         }

         if (sd(sind).vrtx(dir) != td(i).vrtx((j+1)%3) && sd(sind).vrtx(1-dir) != td(i).vrtx((j+2)%3)) {
            *log << "failed vrtx check tind" << i << "sind" << sind << std::endl;
            for(i=0;i<nside;++i)
               sd(i).info += 2;
            output("error"); 
            exit(1);
         }    
         
         if (sd(sind).tri(dir) != i) {
            *log << "failed side check tind" << i << "sind" << sind << std::endl;
            for(i=0;i<nside;++i)
               sd(i).info += 2;
            output("error"); 
            exit(1);
         }
         
         if (td(i).tri(j) != sd(sind).tri(1-dir)) {
            *log << "failed ttri check tind" << i << "sind" << sind << std::endl;
            for(i=0;i<nside;++i)
               sd(i).info += 2;
            output("error"); 
            exit(1);
         }
         
         if (td(i).tri(j) > 0) {
            if(td(td(i).tri(j)).info < 0) {
               *log << "references deleted tri" << std::endl;
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
         fscr1(sind) = stgt_to_actual(sind);
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
         
         collapse(sind,endpt);
         tkoutlst(sind);
         
         /* STORE ADAPTATION INFO FOR COMMUNICATION BOUNDARIES */
         sbdry(bnum)->isndbuf(sbdry(bnum)->sndsize()++) = el;
         sbdry(bnum)->isndbuf(sbdry(bnum)->sndsize()++) = endpt;
         
         /* UPDATE AFFECTED SIDE */
         if (vd(saffect).info > -1) tkoutlst(saffect);
         fscr1(saffect) = stgt_to_actual(saffect);
         if (fscr1(saffect) > tolsize) putinlst(saffect);
         ++count;
   #ifdef DEBUG_ADAPT
         std::cout << "collapsed boundary side " << sind << " endpt " << endpt << std::endl; 
         number_str(adapt_file,"adapt",adapt_count++,5);
         output(adapt_file,ftype::grid);
   #endif            
      }
      sbdry(bnum)->isndbuf(0) = sbdry(bnum)->sndsize();
   }
   *log << "#Boundary coarsening finished, " << count << " sides coarsened" << std::endl;
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
         collapse(sind,endpt);
#ifdef DEBUG_ADAPT
         number_str(adapt_file,"adapt",adapt_count++,5);
         output(adapt_file,ftype::grid);
#endif
      }
   }
   return;
}


void mesh::checkintwk() const {
   int i;
   
   for(i=0;i<maxvst;++i)
      if (i1wk(i) != -1) *log << "failed intwk1 check" << i << i1wk(i) << std::endl;
   
   return;
}
