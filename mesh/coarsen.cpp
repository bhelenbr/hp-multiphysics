#include"mesh.h"
#include"utilities.h"
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<cfloat>
#include<assert.h>

int mesh::coarsen(const class mesh& inmesh) {
	int i,j,k,sind;
	int v0, v1, odd;
   int sideord[MAXSB], *sidelst[MAXSB], nsdloop[MAXSB];
   int nloop;
   FLT mindist;
   
   if (!initialized) {
/*		VERTEX STORAGE ALLOCATION */
      maxvst =  MAX((int) (inmesh.maxvst/3.5),10);
      maxsbel = MAX(inmesh.maxsbel/2,10);
      allocate(maxvst);
      bdryalloc(maxsbel);
      nsbd = inmesh.nsbd;
      nvbd = inmesh.nvbd;
      qtree.allocate(vrtx,maxvst);
      initialized = 1;
   }
   
/* PREPARE FOR COARSENING */
/*	USE GLOBAL INTWK2 TO KEEP TRACK OF WHETHER VERTICES ARE IN */   
	for(i=0;i<inmesh.nvrtx;++i)
		fltwk[i] = 1.0e8;
	
	for(i=0;i<inmesh.nside;++i) {
		v0 = inmesh.svrtx[i][0];
		v1 = inmesh.svrtx[i][1];
		fltwk[v0] = MIN(inmesh.distance(v0,v1),fltwk[v0]);
		fltwk[v1] = MIN(inmesh.distance(v0,v1),fltwk[v1]);
	}

	for(i=0;i<inmesh.nvrtx;++i)
		fltwk[i] *= 1.6; /* SHOULD BE BETWEEN 1.5 and 2.0 */

	nvrtx = 0;
	nside = 0;
	ntri  = 0;

/*	COARSEN SIDES	*/
	for(i=0;i<nsbd;++i) {
		sbdry[i].num = 0;
		sbdry[i].type = inmesh.sbdry[i].type;

/*		CHECK IF FIRST POINT INSERTED*/
		v0 = inmesh.svrtx[inmesh.sbdry[i].el[0]][0];
		if (intwk2[v0] < 0) {
			vrtx[nvrtx][0] = inmesh.vrtx[v0][0];
			vrtx[nvrtx][1] = inmesh.vrtx[v0][1];
			intwk2[v0] = nvrtx;
			svrtx[nside][0] = nvrtx;
			++nvrtx;
         assert(nvrtx < maxvst-1);
		}
		else { 
			svrtx[nside][0] = intwk2[v0];
		}

      odd = inmesh.sbdry[i].num%2;
      if (odd) {
         for(j=2;j<inmesh.sbdry[i].num/2;j+=2) {
            v0 = inmesh.svrtx[inmesh.sbdry[i].el[j]][0];
            vrtx[nvrtx][0] = inmesh.vrtx[v0][0];
            vrtx[nvrtx][1] = inmesh.vrtx[v0][1];
            intwk2[v0] = nvrtx;
            svrtx[nside][1] = nvrtx;
            sbdry[i].el[sbdry[i].num] = nside;
            stri[nside][1] = -1;
            ++nside; 
            ++sbdry[i].num;
            svrtx[nside][0] = nvrtx;
            ++nvrtx;
            assert(sbdry[i].num < maxsbel -1);
            assert(nside < maxvst -1);
            assert(nvrtx < maxvst -1);
         } 

/*			MIDDLE POINT OF ODD NUMBERED SIDE */
         if (inmesh.sbdry[i].num > 1) {
            j = inmesh.sbdry[i].num/2;
            v0 = inmesh.svrtx[inmesh.sbdry[i].el[j]][0];
            v1 = inmesh.svrtx[inmesh.sbdry[i].el[j]][1];
            vrtx[nvrtx][0] = 0.5*(inmesh.vrtx[v0][0] +inmesh.vrtx[v1][0]);
            vrtx[nvrtx][1] = 0.5*(inmesh.vrtx[v0][1] +inmesh.vrtx[v1][1]);
            intwk2[v0] = nvrtx;
            intwk2[v1]= nvrtx;
            svrtx[nside][1] = nvrtx;
            sbdry[i].el[sbdry[i].num] = nside;
            stri[nside][1] = -1;
            ++nside; 
            ++sbdry[i].num;
            svrtx[nside][0] = nvrtx;
            ++nvrtx;
            assert(sbdry[i].num < maxsbel -1);
            assert(nside < maxvst -1);
            assert(nvrtx < maxvst -1); 
         }
         
         for(j = inmesh.sbdry[i].num -((inmesh.sbdry[i].num-2)/4)*2;j<inmesh.sbdry[i].num;j+=2) {
            v0 = inmesh.svrtx[inmesh.sbdry[i].el[j]][0];
            vrtx[nvrtx][0] = inmesh.vrtx[v0][0];
            vrtx[nvrtx][1] = inmesh.vrtx[v0][1];
            intwk2[v0] = nvrtx;
            svrtx[nside][1] = nvrtx;
            sbdry[i].el[sbdry[i].num] = nside;
            stri[nside][1] = -1;
            ++nside; 
            ++sbdry[i].num;
            svrtx[nside][0] = nvrtx;
            ++nvrtx;
            assert(sbdry[i].num < maxsbel -1);
            assert(nside < maxvst -1);
            assert(nvrtx < maxvst -1);
         }
      }
      else {
         for(j=2;j<inmesh.sbdry[i].num;j+=2) {
            v0 = inmesh.svrtx[inmesh.sbdry[i].el[j]][0];
            vrtx[nvrtx][0] = inmesh.vrtx[v0][0];
            vrtx[nvrtx][1] = inmesh.vrtx[v0][1];
            intwk2[v0] = nvrtx;
            svrtx[nside][1] = nvrtx;
            sbdry[i].el[sbdry[i].num] = nside;
            stri[nside][1] = -1;
            ++nside; 
            ++sbdry[i].num;
            svrtx[nside][0] = nvrtx;
            ++nvrtx;
            assert(sbdry[i].num < maxsbel -1);
            assert(nside < maxvst -1);
            assert(nvrtx < maxvst -1);
         }
      }
      
/*		INSERT LAST POINT */
		v0 = inmesh.svrtx[inmesh.sbdry[i].el[inmesh.sbdry[i].num-1]][1];
		if (intwk2[v0] < 0) {
			vrtx[nvrtx][0] = inmesh.vrtx[v0][0];
			vrtx[nvrtx][1] = inmesh.vrtx[v0][1];
			intwk2[v0] = nvrtx;
			svrtx[nside][1] = nvrtx;
			sbdry[i].el[sbdry[i].num] = nside;
         stri[nside][1] = -1;
			++nside;
			++sbdry[i].num;
			++nvrtx;
		}
		else {
			svrtx[nside][1] = intwk2[v0];
			sbdry[i].el[sbdry[i].num] = nside;
         stri[nside][1] = -1;
			++nside;
			++sbdry[i].num;
		}
      assert(sbdry[i].num < maxsbel -1);
      assert(nside < maxvst -1);
      assert(nvrtx < maxvst -1);      
	}
   
/*	MOVE VERTEX BDRY INFORMATION */
   for(i=0;i<inmesh.nvbd;++i) {
      vbdry[i].type = inmesh.vbdry[i].type;
      vbdry[i].num = inmesh.vbdry[i].num;
      for(j=0;j<vbdry[i].num;++j)
         vbdry[i].el[j] = intwk2[inmesh.vbdry[i].el[j]];
   }
   
   treeinit();

/*	CREATE ORDERED LIST OF SIDES */
/*	STORAGE FOR SIDELST */
   for(i=0;i<nsbd;++i)
      sidelst[i] = new int[maxsbel*nsbd];

/*	TEMPORARY LIST OF BOUNDARIES */
   for(i=0;i<nsbd;++i) 
      sideord[i] = i;

   nloop = 0;
   nsdloop[0] = 0;
	for(i=0; i<nsbd; ++i) {
      for(j=0;j<sbdry[sideord[i]].num;++j)
         sidelst[nloop][nsdloop[nloop]++] = sbdry[sideord[i]].el[j] +1;

      v0 = svrtx[sbdry[sideord[i]].el[sbdry[sideord[i]].num-1]][1];
		for(j=i+1; j<nsbd;++j) {
			if (v0 ==  svrtx[sbdry[sideord[j]].el[0]][0]) {
            k = sideord[i+1];
				sideord[i+1] = sideord[j];
            sideord[j] = k;
				goto FINDNEXT;
			}
		}
/*		NEW LOOP */
      nsdloop[++nloop] = 0;       
FINDNEXT:
      continue;
	}

   if (nloop == 0) {
      printf("Problem with boundaries nloop: %d\n",nloop);
      exit(1);
   }
   
/*	CREATE INITIAL TRIANGULATION */            
   triangulate(sidelst,nsdloop,nloop);
         
   for(i=0;i<10;++i) 
      delete []sidelst[i];

/****************************************************/			
/*	Boyer-Watson Algorithm to insert interior points */
/****************************************************/
/*	MARK ALL FINE MESH BOUNDARY NODES SO WE KNOW NOT TO INSERT */
   for(i=0;i<inmesh.nsbd;++i) {
      for(j=0;j<inmesh.sbdry[i].num;++j) {
         sind = inmesh.sbdry[i].el[j];
         intwk2[inmesh.svrtx[sind][0]] = 0;
         intwk2[inmesh.svrtx[sind][1]] = 0;
      }
   }

   for(i=0;i<inmesh.nvrtx;++i) {
		if (intwk2[i] == 0) continue;
      
      mindist = qtree.nearpt(inmesh.vrtx[i][0],inmesh.vrtx[i][1],j);
      if (sqrt(mindist) < fltwk[i]) continue;
      
      insert(inmesh.vrtx[i][0],inmesh.vrtx[i][1]);
   }
   cnt_nbor();
   
/*	RESET INTWK2 */
   for(i=0;i<inmesh.nsbd;++i) {
      for(j=0;j<inmesh.sbdry[i].num;++j) {
         sind = inmesh.sbdry[i].el[j];
         intwk2[inmesh.svrtx[sind][0]] = -1;
         intwk2[inmesh.svrtx[sind][1]] = -1;
      }
   }
   
   bdrylabel();
   initvlngth();
   
/*	PRINT SOME GENERAL DEBUGGING INFO */
   printf("#\n#\n#COARSE MESH\n");
   printf("#MAXVST %d VERTICES %d SIDES %d ELEMENTS %d\n",maxvst,nvrtx,nside,ntri);
   
/* PRINT BOUNDARY INFO */
   for(i=0;i<nsbd;++i)
      printf("#MAX %d BDRY %d TYPE %d SIDES %d\n",maxsbel,i,sbdry[i].type,sbdry[i].num);

	return(1);
}


int mesh::smooth_cofa(int niter) {
	int iter,sind,i,j,n,v0,v1;
   
   for(i=0;i<nvrtx;++i)
      vinfo[i] = 0;
      
   for(i=0;i<nsbd;++i) {
      for(j=0;j<sbdry[i].num;++j) {
         sind = sbdry[i].el[j];
         vinfo[svrtx[sind][0]] = -1;
         vinfo[svrtx[sind][1]] = -1;
      }
   }
	
	for(iter=0; iter< niter; ++iter) {
/*		SMOOTH POINT DISTRIBUTION X*/
      for(n=0;n<ND;++n) {
         for(i=0;i<nvrtx;++i)
            fltwk[i] = 0.0;
   
         for(i=0;i<nside;++i) {
            v0 = svrtx[i][0];
            v1 = svrtx[i][1];
            fltwk[v0] += vrtx[v1][n];
            fltwk[v1] += vrtx[v0][n];
         }
   
         for(i=0;i<nvrtx;++i) {
            if (vinfo[i] == 0) {
               vrtx[i][n] = fltwk[i]/nnbor[i];
            }
         }
      }
	}

	return(1);
}