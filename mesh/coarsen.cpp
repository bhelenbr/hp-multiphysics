#include "mesh.h"
#include "boundary.h"
#include <utilities.h>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <assert.h>
#include <typeinfo>


int mesh::coarsen(FLT factor, const class mesh& inmesh) {
   int i,j,n,sind;
   int v0, v1, odd;
   FLT mindist;
   
   if (!initialized) {
      /* VERTEX STORAGE ALLOCATION */
      allocate_duplicate(1.9,inmesh);
   }
   
   log = inmesh.log;
   
   /* PREPARE FOR COARSENING */
   for(i=0;i<inmesh.nvrtx;++i)
      fscr1(i) = 1.0e8;
   
   for(i=0;i<inmesh.nside;++i) {
      v0 = inmesh.sd(i).vrtx(0);
      v1 = inmesh.sd(i).vrtx(1);
      fscr1(v0) = MIN(inmesh.distance(v0,v1),fscr1(v0));
      fscr1(v1) = MIN(inmesh.distance(v0,v1),fscr1(v1));
   }

   for(i=0;i<inmesh.nvrtx;++i)
      fscr1(i) *= factor; /* SHOULD BE BETWEEN 1.5 and 2.0 */
      
   nvrtx = 0;
   nside = 0;
   ntri  = 0;

   /* USE I1WK TO KEEP TRACK OF INDICES */
   
   /* COARSEN SIDES   */
   for(i=0;i<nsbd;++i) {
      sbdry(i)->nel = 0;
      if (typeid(sbdry(i)) != typeid(inmesh.sbdry(i))) {
         *log << "can't coarsen into object with different boundaries" << std::endl;
         exit(1);
      }

      /* CHECK IF FIRST POINT INSERTED*/
      v0 = inmesh.sd(inmesh.sbdry(i)->el(0)).vrtx(0);
      if (i1wk(v0) < 0) {
         for(n=0;n<ND;++n)
            vrtx(nvrtx)(n) = inmesh.vrtx(v0)(n);
         i1wk(v0) = nvrtx;
         sd(nside).vrtx(0) = nvrtx;
         ++nvrtx;
      }
      else { 
         sd(nside).vrtx(0) = i1wk(v0);
      }

      odd = inmesh.sbdry(i)->nel%2;
      if (odd) {
         for(j=2;j<inmesh.sbdry(i)->nel/2;j+=2) {
            v0 = inmesh.sd(inmesh.sbdry(i)->el(j)).vrtx(0);
            for(n=0;n<ND;++n)
               vrtx(nvrtx)(n) = inmesh.vrtx(v0)(n);
            i1wk(v0) = nvrtx;
            sd(nside).vrtx(1) = nvrtx;
            sbdry(i)->el(sbdry(i)->nel) = nside;
            sd(nside).tri(1) = -1;
            ++nside; 
            ++sbdry(i)->nel;
            sd(nside).vrtx(0) = nvrtx;
            ++nvrtx;
         } 

         /* MIDDLE POINT OF ODD NUMBERED SIDE */
         if (inmesh.sbdry(i)->nel > 1) {
            j = inmesh.sbdry(i)->nel/2;
            v0 = inmesh.sd(inmesh.sbdry(i)->el(j)).vrtx(0);
            v1 = inmesh.sd(inmesh.sbdry(i)->el(j)).vrtx(1);
            for(n=0;n<ND;++n)
               vrtx(nvrtx)(n) = 0.5*(inmesh.vrtx(v0)(n) +inmesh.vrtx(v1)(n));
            i1wk(v0) = nvrtx;
            i1wk(v1)= nvrtx;
            sd(nside).vrtx(1) = nvrtx;
            sbdry(i)->el(sbdry(i)->nel) = nside;
            sd(nside).tri(1) = -1;
            ++nside; 
            ++sbdry(i)->nel;
            sd(nside).vrtx(0) = nvrtx;
            ++nvrtx;
         }
         
         for(j = inmesh.sbdry(i)->nel -((inmesh.sbdry(i)->nel-2)/4)*2;j<inmesh.sbdry(i)->nel;j+=2) {
            v0 = inmesh.sd(inmesh.sbdry(i)->el(j)).vrtx(0);
            for(n=0;n<ND;++n)
               vrtx(nvrtx)(n) = inmesh.vrtx(v0)(n);
            i1wk(v0) = nvrtx;
            sd(nside).vrtx(1) = nvrtx;
            sbdry(i)->el(sbdry(i)->nel) = nside;
            sd(nside).tri(1) = -1;
            ++nside; 
            ++sbdry(i)->nel;
            sd(nside).vrtx(0) = nvrtx;
            ++nvrtx;
         }
      }
      else {
         for(j=2;j<inmesh.sbdry(i)->nel;j+=2) {
            v0 = inmesh.sd(inmesh.sbdry(i)->el(j)).vrtx(0);
            for(n=0;n<ND;++n)
               vrtx(nvrtx)(n) = inmesh.vrtx(v0)(n);
            i1wk(v0) = nvrtx;
            sd(nside).vrtx(1) = nvrtx;
            sbdry(i)->el(sbdry(i)->nel) = nside;
            sd(nside).tri(1) = -1;
            ++nside; 
            ++sbdry(i)->nel;
            sd(nside).vrtx(0) = nvrtx;
            ++nvrtx;
         }
      }
      
      /* INSERT LAST POINT */
      v0 = inmesh.sd(inmesh.sbdry(i)->el(inmesh.sbdry(i)->nel-1)).vrtx(1);
      if (i1wk(v0) < 0) {
         for(n=0;n<ND;++n)
            vrtx(nvrtx)(n) = inmesh.vrtx(v0)(n);
         i1wk(v0) = nvrtx;
         sd(nside).vrtx(1) = nvrtx;
         sbdry(i)->el(sbdry(i)->nel) = nside;
         sd(nside).tri(1) = -1;
         ++nside;
         ++sbdry(i)->nel;
         ++nvrtx;
      }
      else {
         sd(nside).vrtx(1) = i1wk(v0);
         sbdry(i)->el(sbdry(i)->nel) = nside;
         sd(nside).tri(1) = -1;
         ++nside;
         ++sbdry(i)->nel;
      }
   }
   
   /* MOVE VERTEX BDRY INFORMATION */
   for(i=0;i<inmesh.nvbd;++i) {
      vbdry(i)->copy(*inmesh.vbdry(i));
      vbdry(i)->v0 = i1wk(inmesh.vbdry(i)->v0);
   }
   
   /* RESET i1wk */
   for(i=0;i<inmesh.nsbd;++i) {
      for(j=0;j<inmesh.sbdry(i)->nel;++j) {
         sind = inmesh.sbdry(i)->el(j);
         i1wk(inmesh.sd(sind).vrtx(0)) = -1;
         i1wk(inmesh.sd(sind).vrtx(1)) = -1;
      }
   }
   
   treeinit();

   for(i=0;i<nside;++i)
      i2wk_lst1(i) = i+1;
      
   triangulate(nside);

   /****************************************************/         
   /* Boyer-Watson Algorithm to insert interior points */
   /****************************************************/
   for(i=0;i<inmesh.nvrtx;++i) {
      if (i2wk(i) == 0) continue;
      
      mindist = qtree.nearpt(inmesh.vrtx(i).data(),j);
      if (sqrt(mindist) < fscr1(i)) continue;
            
      insert(inmesh.vrtx(i));
   }
   cnt_nbor();
   
   /* RESET i2wk */
   for(i=0;i<inmesh.nsbd;++i) {
      for(j=0;j<inmesh.sbdry(i)->nel;++j) {
         sind = inmesh.sbdry(i)->el(j);
         i2wk(inmesh.sd(sind).vrtx(0)) = -1;
         i2wk(inmesh.sd(sind).vrtx(1)) = -1;
      }
   }
   
   bdrylabel();
   initvlngth();
   
   /* PRINT SOME GENERAL DEBUGGING INFO */
   *log << "#" << std::endl << "#COARSE MESH " << std::endl;
   *log << "#MAXVST:" << maxvst << " VERTICES:" << nvrtx << " SIDES:" << nside << " ELEMENTS:" << ntri << std::endl;   
   /* PRINT BOUNDARY INFO */
   for(i=0;i<nsbd;++i)
      *log << "#" << sbdry(i)->idprefix << " " << sbdry(i)->mytype << " " << sbdry(i)->nel << std::endl;

   return(1);
}


int mesh::smooth_cofa(int niter) {
   int iter,sind,i,j,n,v0,v1;
      
   for(i=0;i<nvrtx;++i)
      vd(i).info = 0;
      
   for(i=0;i<nsbd;++i) {
      for(j=0;j<sbdry(i)->nel;++j) {
         sind = sbdry(i)->el(j);
         vd(sd(sind).vrtx(0)).info = -1;
         vd(sd(sind).vrtx(1)).info = -1;
      }
   }
   
   for(iter=0; iter< niter; ++iter) {
      /* SMOOTH POINT DISTRIBUTION X*/
      for(n=0;n<ND;++n) {
         for(i=0;i<nvrtx;++i)
            fscr1(i) = 0.0;
   
         for(i=0;i<nside;++i) {
            v0 = sd(i).vrtx(0);
            v1 = sd(i).vrtx(1);
            fscr1(v0) += vrtx(v1)(n);
            fscr1(v1) += vrtx(v0)(n);
         }
   
         for(i=0;i<nvrtx;++i) {
            if (vd(i).info == 0) {
               vrtx(i)(n) = fscr1(i)/vd(i).nnbor;
            }
         }
      }
   }
   
   return(1);
}

void mesh::coarsen2(FLT factor, const class mesh &inmesh, FLT size_reduce) {
   int i;
   
   if (!initialized) {
      allocate_duplicate(size_reduce,inmesh);
   }      
   copy(inmesh);
   initvlngth();
   for(i=0;i<nvrtx;++i)
      vlngth(i) = factor*vlngth(i);
      
   yaber(1.414);
   qtree.reinit();  // REMOVES UNUSED QUADS
   
   return;
}

void mesh::coarsen3() {
   int i,j,v0,tind,sind,node,cnt=0;
   
   /* SET-UP ADAPTION TRACKING STUFF */
   /* For vertices 0 = untouched, 1 = touched, 2 = deleted, 3 = special */
   /* For sides 0 = untouced, 1 = touched, 2 = deleted */
   /* For triangles 0 = untouched, 1 = touched, 2 = deleted, 3 = searched */
   for(i=0;i<maxvst;++i)
      td(i).info = 0;
   
   /* COARSEN SIDE EDGES FIRST */
   for(i=0;i<nsbd;++i) {
      *log << "coarsening boundary " << i << ": type " << sbdry(i)->mytype << " sides " << sbdry(i)->nel << std::endl;
      for(j=0;j<sbdry(i)->nel;++j) {
         sind = sbdry(i)->el(j);
         v0 = sd(sind).vrtx(0);
         if (vd(v0).info > 0) {
            /* COARSEN SIDE */
            collapse(sind,0);
            vd(v0).info = -1;
            ++cnt;
         }
      }
   }
   
   
   for(i=0;i<nvrtx;++i) {
      if (vd(i).info > 0) {
         tind = vd(i).tri;
         for(j=0;j<3;++j) {
            if (td(tind).vrtx(j) == i) {
               break;
            }
         }
         j = (j+1)%3;
         sind = td(tind).side(j);
         node = (1+td(tind).sign(j))/2;
         collapse(sind,node);
         vd(v0).info = -1;
         ++cnt;
      }
   }
   
   cleanup_after_adapt();
   
   *log << "#Coarsen finished: " << cnt << " sides coarsened" << std::endl;
   for(i=0;i<nsbd;++i)
      *log << "boundary " << i << ": type " << sbdry(i)->mytype << " sides " << sbdry(i)->nel << std::endl;
   
   return;
}
   
   

