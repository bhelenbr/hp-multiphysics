#include "mesh.h"
#include "boundary.h"
#include "utilities.h"
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
   /* USE GLOBAL i2wk TO KEEP TRACK OF WHETHER VERTICES ARE IN */   
   for(i=0;i<inmesh.nvrtx;++i)
      fwk[i] = 1.0e8;
   
   for(i=0;i<inmesh.nside;++i) {
      v0 = inmesh.sd[i].vrtx[0];
      v1 = inmesh.sd[i].vrtx[1];
      fwk[v0] = MIN(inmesh.distance(v0,v1),fwk[v0]);
      fwk[v1] = MIN(inmesh.distance(v0,v1),fwk[v1]);
   }

   for(i=0;i<inmesh.nvrtx;++i)
      fwk[i] *= factor; /* SHOULD BE BETWEEN 1.5 and 2.0 */

   nvrtx = 0;
   nside = 0;
   ntri  = 0;

   /* COARSEN SIDES   */
   for(i=0;i<nsbd;++i) {
      sbdry[i]->nel = 0;
      if (typeid(sbdry[i]) != typeid(inmesh.sbdry[i])) {
         *log << "can't coarsen into object with different boundaries" << std::endl;
         exit(1);
      }

      /* CHECK IF FIRST POINT INSERTED*/
      v0 = inmesh.sd[inmesh.sbdry[i]->el[0]].vrtx[0];
      if (i2wk[v0] < 0) {
         for(n=0;n<ND;++n)
            vrtx[nvrtx][n] = inmesh.vrtx[v0][n];
         i2wk[v0] = nvrtx;
         sd[nside].vrtx[0] = nvrtx;
         ++nvrtx;
         assert(nvrtx < maxvst-1);
      }
      else { 
         sd[nside].vrtx[0] = i2wk[v0];
      }

      odd = inmesh.sbdry[i]->nel%2;
      if (odd) {
         for(j=2;j<inmesh.sbdry[i]->nel/2;j+=2) {
            v0 = inmesh.sd[inmesh.sbdry[i]->el[j]].vrtx[0];
            for(n=0;n<ND;++n)
               vrtx[nvrtx][n] = inmesh.vrtx[v0][n];
            i2wk[v0] = nvrtx;
            sd[nside].vrtx[1] = nvrtx;
            sbdry[i]->el[sbdry[i]->nel] = nside;
            sd[nside].tri[1] = -1;
            ++nside; 
            ++sbdry[i]->nel;
            sd[nside].vrtx[0] = nvrtx;
            ++nvrtx;
            assert(sbdry[i]->nel < sbdry[i]->maxel -1);
            assert(nside < maxvst -1);
            assert(nvrtx < maxvst -1);
         } 

         /* MIDDLE POINT OF ODD NUMBERED SIDE */
         if (inmesh.sbdry[i]->nel > 1) {
            j = inmesh.sbdry[i]->nel/2;
            v0 = inmesh.sd[inmesh.sbdry[i]->el[j]].vrtx[0];
            v1 = inmesh.sd[inmesh.sbdry[i]->el[j]].vrtx[1];
            for(n=0;n<ND;++n)
               vrtx[nvrtx][n] = 0.5*(inmesh.vrtx[v0][n] +inmesh.vrtx[v1][n]);
            i2wk[v0] = nvrtx;
            i2wk[v1]= nvrtx;
            sd[nside].vrtx[1] = nvrtx;
            sbdry[i]->el[sbdry[i]->nel] = nside;
            sd[nside].tri[1] = -1;
            ++nside; 
            ++sbdry[i]->nel;
            sd[nside].vrtx[0] = nvrtx;
            ++nvrtx;
            assert(sbdry[i]->nel < sbdry[i]->maxel -1);
            assert(nside < maxvst -1);
            assert(nvrtx < maxvst -1); 
         }
         
         for(j = inmesh.sbdry[i]->nel -((inmesh.sbdry[i]->nel-2)/4)*2;j<inmesh.sbdry[i]->nel;j+=2) {
            v0 = inmesh.sd[inmesh.sbdry[i]->el[j]].vrtx[0];
            for(n=0;n<ND;++n)
               vrtx[nvrtx][n] = inmesh.vrtx[v0][n];
            i2wk[v0] = nvrtx;
            sd[nside].vrtx[1] = nvrtx;
            sbdry[i]->el[sbdry[i]->nel] = nside;
            sd[nside].tri[1] = -1;
            ++nside; 
            ++sbdry[i]->nel;
            sd[nside].vrtx[0] = nvrtx;
            ++nvrtx;
            assert(sbdry[i]->nel < sbdry[i]->maxel -1);
            assert(nside < maxvst -1);
            assert(nvrtx < maxvst -1);
         }
      }
      else {
         for(j=2;j<inmesh.sbdry[i]->nel;j+=2) {
            v0 = inmesh.sd[inmesh.sbdry[i]->el[j]].vrtx[0];
            for(n=0;n<ND;++n)
               vrtx[nvrtx][n] = inmesh.vrtx[v0][n];
            i2wk[v0] = nvrtx;
            sd[nside].vrtx[1] = nvrtx;
            sbdry[i]->el[sbdry[i]->nel] = nside;
            sd[nside].tri[1] = -1;
            ++nside; 
            ++sbdry[i]->nel;
            sd[nside].vrtx[0] = nvrtx;
            ++nvrtx;
            assert(sbdry[i]->nel < sbdry[i]->maxel -1);
            assert(nside < maxvst -1);
            assert(nvrtx < maxvst -1);
         }
      }
      
      /* INSERT LAST POINT */
      v0 = inmesh.sd[inmesh.sbdry[i]->el[inmesh.sbdry[i]->nel-1]].vrtx[1];
      if (i2wk[v0] < 0) {
         for(n=0;n<ND;++n)
            vrtx[nvrtx][n] = inmesh.vrtx[v0][n];
         i2wk[v0] = nvrtx;
         sd[nside].vrtx[1] = nvrtx;
         sbdry[i]->el[sbdry[i]->nel] = nside;
         sd[nside].tri[1] = -1;
         ++nside;
         ++sbdry[i]->nel;
         ++nvrtx;
      }
      else {
         sd[nside].vrtx[1] = i2wk[v0];
         sbdry[i]->el[sbdry[i]->nel] = nside;
         sd[nside].tri[1] = -1;
         ++nside;
         ++sbdry[i]->nel;
      }
      assert(sbdry[i]->nel < sbdry[i]->maxel -1);
      assert(nside < maxvst -1);
      assert(nvrtx < maxvst -1);      
   }
   
   /* MOVE VERTEX BDRY INFORMATION */
   for(i=0;i<inmesh.nvbd;++i) {
      vbdry[i]->copy(*inmesh.vbdry[i]);
      vbdry[i]->v0 = i2wk[inmesh.vbdry[i]->v0];
   }
   
   treeinit();

   for(i=0;i<nside;++i)
      i1wk[i] = i+1;
      
   triangulate(nside);

   /****************************************************/         
   /* Boyer-Watson Algorithm to insert interior points */
   /****************************************************/
   /* MARK ALL FINE MESH BOUNDARY NODES SO WE KNOW NOT TO INSERT */
   for(i=0;i<inmesh.nsbd;++i) {
      for(j=0;j<inmesh.sbdry[i]->nel;++j) {
         sind = inmesh.sbdry[i]->el[j];
         i2wk[inmesh.sd[sind].vrtx[0]] = 0;
         i2wk[inmesh.sd[sind].vrtx[1]] = 0;
      }
   }

   for(i=0;i<inmesh.nvrtx;++i) {
      if (i2wk[i] == 0) continue;
      
      mindist = qtree.nearpt(inmesh.vrtx[i],j);
      if (sqrt(mindist) < fwk[i]) continue;
            
      insert(inmesh.vrtx[i]);
   }
   cnt_nbor();
   
   /* RESET i2wk */
   for(i=0;i<inmesh.nsbd;++i) {
      for(j=0;j<inmesh.sbdry[i]->nel;++j) {
         sind = inmesh.sbdry[i]->el[j];
         i2wk[inmesh.sd[sind].vrtx[0]] = -1;
         i2wk[inmesh.sd[sind].vrtx[1]] = -1;
      }
   }
   
   bdrylabel();
   initvlngth();
   
   /* PRINT SOME GENERAL DEBUGGING INFO */
   *log << "#" << std::endl << "#COARSE MESH " << std::endl;
   *log << "#MAXVST:" << maxvst << " VERTICES:" << nvrtx << " SIDES:" << nside << " ELEMENTS:" << ntri << std::endl;   
   /* PRINT BOUNDARY INFO */
   for(i=0;i<nsbd;++i)
      *log << "#" << sbdry[i]->idprefix << " " << typeid(*sbdry[i]).name() << " " << sbdry[i]->nel << std::endl;

   return(1);
}


int mesh::smooth_cofa(int niter) {
   int iter,sind,i,j,n,v0,v1;
   
   for(i=0;i<nvrtx;++i)
      vd[i].info = 0;
      
   for(i=0;i<nsbd;++i) {
      for(j=0;j<sbdry[i]->nel;++j) {
         sind = sbdry[i]->el[j];
         vd[sd[sind].vrtx[0]].info = -1;
         vd[sd[sind].vrtx[1]].info = -1;
      }
   }
   
   for(iter=0; iter< niter; ++iter) {
      /* SMOOTH POINT DISTRIBUTION X*/
      for(n=0;n<ND;++n) {
         for(i=0;i<nvrtx;++i)
            fwk[i] = 0.0;
   
         for(i=0;i<nside;++i) {
            v0 = sd[i].vrtx[0];
            v1 = sd[i].vrtx[1];
            fwk[v0] += vrtx[v1][n];
            fwk[v1] += vrtx[v0][n];
         }
   
         for(i=0;i<nvrtx;++i) {
            if (vd[i].info == 0) {
               vrtx[i][n] = fwk[i]/vd[i].nnbor;
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
      vlngth[i] = factor*vlngth[i];
      
   yaber(1.414, 0);
   qtree.reinit();  // REMOVES UNUSED QUADS
   
   return;
}
   
   

