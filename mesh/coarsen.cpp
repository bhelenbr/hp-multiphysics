#include"mesh.h"
#include"utilities.h"
#include<cstdio>
#include<cstdlib>
#include<cstring>
#include<cmath>
#include<cfloat>
#include<assert.h>

int mesh::coarsen(const class mesh& fmesh) {
	int i,j,sind;
	int v0, v1;
   int nloop;
   FLT mindist;
   
   if (!initialized) {
/*		VERTEX STORAGE ALLOCATION */
      maxvst = (int) (fmesh.maxvst/3.5);
      vrtx = (FLT (*)[ND]) xmalloc(maxvst*ND*sizeof(FLT));
      vlngth = new FLT[maxvst];
      vinfo = new int[maxvst+1];
      ++vinfo;  //  ALLOWS US TO ACCES VINFO[-1]
      nnbor = new int[maxvst];
      vtri = new int[maxvst];
	
/*		VERTEX BOUNDARY STORAGE INFORMATION */		
      nvbd = fmesh.nvbd;
      maxvbel = fmesh.maxvbel;
      for(i=0;i<nvbd;++i)
         vbdry[i].el = new int[maxvbel];

/*		SIDE STORAGE ALLOCATION */	
      svrtx  = (int (*)[2]) xmalloc(maxvst*2*sizeof(int));
      stri   = (int (*)[2]) xmalloc(maxvst*2*sizeof(int));
      sinfo = new int[maxvst+1];
      ++sinfo; // ALLOWS US TO ACCESS SINFO[-1]

/*		SIDE BOUNDARY STORAGE ALLOCATION */
      nsbd = fmesh.nsbd;
      maxsbel = MAX(fmesh.maxsbel/2,2);
      for(i=0;i<nsbd;++i)
         sbdry[i].el = new int[maxsbel];      

/*		TRIANGLE ALLOCATION */			
      tvrtx = (int (*)[3]) xmalloc(maxvst*3*sizeof(int));
      ttri = (int (*)[3]) xmalloc(maxvst*3*sizeof(int));
      tside = new struct tsidedata[maxvst];
      tinfo = new int[maxvst+1];
      ++tinfo; // ALLOWS US TO ACCESS TINFO(-1)
   }
   
   initialized = 1;


/* PREPARE FOR COARSENING */
/*	USE GLOBAL INTWK2 TO KEEP TRACK OF WHETHER VERTICES ARE IN */   
	for(i=0;i<fmesh.nvrtx;++i)
		fltwk[i] = 1.0e8;
	
	for(i=0;i<fmesh.nside;++i) {
		v0 = fmesh.svrtx[i][0];
		v1 = fmesh.svrtx[i][1];
		fltwk[v0] = MIN(fmesh.distance(v0,v1),fltwk[v0]);
		fltwk[v1] = MIN(fmesh.distance(v0,v1),fltwk[v1]);
	}

	for(i=0;i<fmesh.nvrtx;++i)
		fltwk[i] *= 1.6; /* SHOULD BE BETWEEN 1.5 and 2.0 */

	nvrtx = 0;
	nside = 0;
	ntri  = 0;

/*	COARSEN SIDES	*/
	for(i=0;i<nsbd;++i) {
		sbdry[i].num = 0;
		sbdry[i].type = fmesh.sbdry[i].type;
      if (sbdry[i].type & ALLD_MP)
         recvbuf[i] = new FLT[8*maxsbel];

/*		CHECK IF FIRST POINT INSERTED*/
		v0 = fmesh.svrtx[fmesh.sbdry[i].el[0]][0];
		if (intwk2[v0] < 0) {
			vrtx[nvrtx][0] = fmesh.vrtx[v0][0];
			vrtx[nvrtx][1] = fmesh.vrtx[v0][1];
			intwk2[v0] = nvrtx;
			svrtx[nside][0] = nvrtx;
			++nvrtx;
         assert(nvrtx < maxvst-1);
		}
		else { 
			svrtx[nside][0] = intwk2[v0];
		}

		for(j=2;j<fmesh.sbdry[i].num;j+=2) {
			v0 = fmesh.svrtx[fmesh.sbdry[i].el[j]][0];
			vrtx[nvrtx][0] = fmesh.vrtx[v0][0];
			vrtx[nvrtx][1] = fmesh.vrtx[v0][1];
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
		
/*		INSERT LAST POINT */
		v0 = fmesh.svrtx[fmesh.sbdry[i].el[fmesh.sbdry[i].num-1]][1];
		if (intwk2[v0] < 0) {
			vrtx[nvrtx][0] = fmesh.vrtx[v0][0];
			vrtx[nvrtx][1] = fmesh.vrtx[v0][1];
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

/*	SMOOTH BOUNDARY POINT DISTRIBUTION FOR SIDES WITH UNEVEN NUMBER */
	int s0;
	int s1;
	for(i=0;i<nsbd;++i) {
		if (fmesh.sbdry[i].num%2 != 0 && fmesh.sbdry[i].num > 1) {
			s0 = sbdry[i].el[sbdry[i].num-2];
			s1 = sbdry[i].el[sbdry[i].num-1];
/*			CHECK COLINEARITY */			
			if (abs(area(s0,svrtx[s1][1])) < 10.*FLT_EPSILON) {
				vrtx[svrtx[s0][1]][0] = 
					0.5*(vrtx[svrtx[s0][0]][0] +vrtx[svrtx[s1][1]][0]);
				vrtx[svrtx[s0][1]][1] = 
					0.5*(vrtx[svrtx[s0][0]][1] +vrtx[svrtx[s1][1]][1]);					
			}
		}
	}
   
   treeinit();

/*	ALLOCATE TEMPORARY STORAGE FOR ORDERED BOUNDARY LOOPS */
   int **sidelst;
   int *nsdloop;
   
   sidelst = new int *[10];
   nsdloop = new int[10];
   for(i=0;i<10;++i) {
      sidelst[i] = new int[maxsbel*nsbd];
      nsdloop[i] = 0;
   }

/*	MAKE ORDERED LIST OF SIDES */
/*	EACH LOOP ASSIGNED TO DIFFERENT ARRAY */
/*	ASSUMES CCW ORDERING OF GROUPS/ELEMENTS */   
   nloop = 0;
   v0 = svrtx[sbdry[0].el[0]][0];
   for(i=0;i<nsbd;++i) {
      for(j=0;j<sbdry[i].num;++j)
         sidelst[nloop][nsdloop[nloop]++] = sbdry[i].el[j] +1;
         
      if (v0 == svrtx[sbdry[i].el[sbdry[i].num-1]][1]) {
         ++nloop;
         if (nloop > 10) {
            printf("Too many internal loops %d\n",nloop);
            exit(1);
         }
         if (i < nsbd -1) v0 = svrtx[sbdry[i+1].el[0]][0];
      }
   }

   if (nloop == 0) {
      printf("Problem with boundaries nloop: %d\n",nloop);
      exit(1);
   }
   
/*	CREATE INITIAL TRIANGULATION */            
   triangulate(sidelst,nsdloop,nloop);

/****************************************************/			
/*	Boyer-Watson Algorithm to insert interior points */
/****************************************************/
/*	MARK ALL FINE MESH BOUNDARY NODES SO WE KNOW NOT TO INSERT */
   for(i=0;i<fmesh.nsbd;++i) {
      for(j=0;j<fmesh.sbdry[i].num;++j) {
         sind = fmesh.sbdry[i].el[j];
         intwk2[fmesh.svrtx[sind][0]] = 0;
         intwk2[fmesh.svrtx[sind][1]] = 0;
      }
   }
   
   for(i=0;i<fmesh.nvrtx;++i) {
		if (intwk2[i] == 0) continue;
      
      mindist = qtree.nearpt(fmesh.vrtx[i][0],fmesh.vrtx[i][1],j);
      if (sqrt(mindist) < fltwk[i]) continue;
      
      insert(fmesh.vrtx[i][0],fmesh.vrtx[i][1]);
   }
   cnt_nbor();
   
/*	RESET INTWK2 */
   for(i=0;i<fmesh.nsbd;++i) {
      for(j=0;j<fmesh.sbdry[i].num;++j) {
         sind = fmesh.sbdry[i].el[j];
         intwk2[fmesh.svrtx[sind][0]] = -1;
         intwk2[fmesh.svrtx[sind][1]] = -1;
      }
   }
   
   bdrylabel();
   initvlngth();

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