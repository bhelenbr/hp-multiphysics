#include"mesh.h"
#include"utilities.h"
#include<cstdio>
#include<float.h>

/*	CREATE SIDELIST FROM TRIANGLE VERTEX LIST */
/*	USES VINFO TO STORE FIRST SIND FROM VERTEX */
/* USES SINFO TO STORE NEXT SIND FROM SIND */
/* TVRTX MUST BE COUNTERCLOCKWISE ORDERED */
void mesh::createsideinfo(void) {
   int i,j,tind,v1,v2,vout,temp,minv,maxv,order,sind,sindprev;
   
   for(i=0;i<nvrtx;++i)
      vinfo[i] = -1;
      
   nside = 0;
   for(tind=0;tind<ntri;++tind) {
      vout = tvrtx[tind][0];
      v1 = tvrtx[tind][1];
      v2 = tvrtx[tind][2];
      for(j=0;j<3;++j) {
/*			CHECK IF SIDE HAS BEEN CREATED ALREADY */
         if (v2 > v1) {
            minv = v1;
            maxv = v2;
            order = 0;
         }
         else {
            minv = v2;
            maxv = v1;
            order = 1;
         }
         sind = vinfo[minv];
         while (sind >= 0) {
            if (maxv == svrtx[sind][order]) {
               if (stri[sind][1] >= 0) {
                  printf("Error: side %d has been matched with Triangle %d 3 times\n",sind,tind);
                  exit(1);
               }
               else {
                  stri[sind][1] = tind;
                  tside[tind].side[j] = sind;
                  tside[tind].sign[j] = -1;
                  goto NEXTTRISIDE;
               }
            }
            sindprev = sind;
            sind = sinfo[sind];
         }
/*			NEW SIDE */
         svrtx[nside][0] = v1;
         svrtx[nside][1] = v2;
         stri[nside][0] = tind;
         stri[nside][1] = -1;
         tside[tind].side[j] = nside;
         tside[tind].sign[j] = 1;
         sinfo[nside] = -1;
         if (vinfo[minv] < 0)
            vinfo[minv] = nside;
         else 
            sinfo[sindprev] = nside;
         ++nside;
NEXTTRISIDE:
         temp = vout;
         vout = v1;
         v1 = v2;
         v2 = temp;
      }
   }

   return;
}

/*	CALCULATE NUMBER OF NEIGHBORS TO EACH CELL */
void mesh::createvtri(void) {
	int i,tind;
	
/*	THIS ALLOWS US TO GET TO LOCAL HIGHER ENTITIES FROM VERTEX NUMBER */
   for (tind=0;tind<ntri;++tind)
      for(i=0;i<3;++i)
         vtri[tvrtx[tind][i]] = tind;
	
	return;
}

/*	CALCULATE NUMBER OF NEIGHBORS TO EACH CELL */
void mesh::cnt_nbor(void) {
	int i;
	
	for (i=0;i<nvrtx;++i)
		nnbor[i] = 0;
	
	for(i=0;i<nside;++i) {
		++nnbor[svrtx[i][0]];
		++nnbor[svrtx[i][1]];
	}

	return;
}


/*	CREATES TRIANGLE TO TRIANGLE POINTER */
void mesh::createttri(void) {
   int tind,sind,j,flip;
   
   for(tind=0;tind<ntri;++tind) {
      for(j=0;j<3;++j) {
         sind = tside[tind].side[j];
         flip = (1 +tside[tind].sign[j])/2;
         ttri[tind][j] = stri[sind][flip];
      }
   }
   
   return;
}

void mesh::treeinit() {
   int i,j,sind,v0;
   FLT x1,y1,x2,y2;
   
   x1 = vrtx[0][0];
   y1 = vrtx[0][1];
   x2 = x1;
   y2 = y1;   
   
   for (i=0;i<nsbd;++i) {
      for(j=0;j<sbdry[i].num;++j) {
         sind = sbdry[i].el[j];
         v0 = svrtx[sind][0];
         x1 = MIN(x1,vrtx[v0][0]);
         y1 = MIN(y1,vrtx[v0][1]);
         x2 = MAX(x2,vrtx[v0][0]);
         y2 = MAX(y2,vrtx[v0][1]);
      }
   }

   qtree.init(vrtx,maxvst,x1-EPSILON,y1-EPSILON,x2+EPSILON,y2+EPSILON);
   
   for(i=0;i<nvrtx;++i) 
      qtree.addpt(i);

   return;
}

/*	REORDERS BOUNDARIES TO BE SEQUENTIAL */
/*	USES INTWK1 & INTWK2 AS WORK ARRAYS */
void mesh::bdrysidereorder(int bdrynum) {
	int i,j,count,total,sind,nsegs,first[MAXSB];
   
   total = sbdry[bdrynum].num;
   
/*	STORE SIDE INDICES BY VERTEX NUMBER */
	for(i=0; i < sbdry[bdrynum].num; ++i) {
      sind = sbdry[bdrynum].el[i];
      intwk1[svrtx[sind][1]] = sind;
      intwk2[svrtx[sind][0]] = sind;
   }

/*	FIND FIRST SIDE(S) */   
   nsegs = 0;
   first[0] = sbdry[bdrynum].el[0];
   for(i=0;i<sbdry[bdrynum].num;++i) {
      sind = sbdry[bdrynum].el[i];
      if (intwk1[svrtx[sind][0]] == -1) 
         first[nsegs++] = sind;
   }
   
/*	MAKE FIRST LIST */
   count = 0;
   sind = first[0];
   do {
      sbdry[bdrynum].el[count++] = sind;
      sind = intwk2[svrtx[sind][1]];
   } while (sind >= 0 && sind != first[0]);
   sbdry[bdrynum].num = count;
   
/*	RESET INTWK FOR FIRST SEG TO -1 */
	for(i=0; i < sbdry[bdrynum].num; ++i) {
      sind = sbdry[bdrynum].el[i];
      intwk1[svrtx[sind][1]] = -1;
      intwk2[svrtx[sind][0]] = -1;
   }
   
   
/*	MAKE NEW BOUNDARIES FOR DISCONNECTED SEGMENTS */
   for(i=1;i<nsegs;++i) {
      sbdry[nsbd].el = new int[maxsbel];
      sbdry[nsbd].num = 0;
      sbdry[nsbd].type = sbdry[bdrynum].type;
      sind = first[i];
      do {
         sbdry[nsbd].el[sbdry[nsbd].num++] = sind;
         sind = intwk2[svrtx[sind][1]];
      } while (sind >= 0 && sind != first[i]);
      count += sbdry[nsbd].num;
      
/*		RESET INTWK FOR SEGMENT */
      for(j=0; j < sbdry[nsbd].num; ++j) {
         sind = sbdry[nsbd].el[j];
         intwk1[svrtx[sind][1]] = -1;
         intwk2[svrtx[sind][0]] = -1;
      }

      ++nsbd;
   }
   
   if (count != total) {
      printf("Error: Multiple loops in one boundary group\n");
      exit(1);
   }

/* UPDATE STRI/TTRI POINTERS ON BOUNDARY */   
   bdrylabel();
      
	return;
}

/*	REARRANGE BOUNDARY GROUPS TO BE IN CCW ORDER */
/*	THIS IS A BIT BRUTE FORCE */

void mesh::bdrygroupreorder(void) {
   int i,j,k,sideord[MAXSB],nloop,strtloop[MAXSB],v0;
   int bnum,bgn,end,count;
   struct boundary tmpsb;
   FLT sum;
   
/*	FIND COUNTERCLOCKWISE ORDERING OF BOUNDARY */
   for(i=0;i<nsbd;++i) 
      sideord[i] = i;

   nloop = 0;
   strtloop[nloop++] = 0;
	for(i=1; i<nsbd; ++i) {
      v0 = svrtx[sbdry[sideord[i-1]].el[sbdry[sideord[i-1]].num-1]][1];
		for(j=i; j<nsbd;++j) {
			if (v0 ==  svrtx[sbdry[sideord[j]].el[0]][0]) {
            k = sideord[i];
				sideord[i] = sideord[j];
            sideord[j] = k;
				goto FINDNEXT;
			}
		}
      strtloop[nloop++] = i;       
FINDNEXT:
      continue;
	}
   strtloop[nloop] = nsbd;
         
/*	DETERMINE WHICH LOOP IS THE EXTERIOR LOOP */
   if (nloop > 1) {
      for(i=0;i<nloop;++i) {
         sum = 0.0;
         for(j=strtloop[i];j<strtloop[i+1];++j) {
            bnum = sideord[j];
            for(k=0;k<sbdry[bnum].num;++k) 
               sum += area(sbdry[bnum].el[k],0);
         }
         if (sum > 0.0) break;
      }
      bgn = strtloop[i];
      end = strtloop[i+1];
   }
   else {
      bgn = 0;
      end = nsbd;
   }
      
/*	REARRANGE GROUPS EXTERIOR GOES FIRST */
/* JUST CHANGES sbel POINTER DOESN'T MOVE WHOLE LIST */
   count = 0;
   for(i=bgn;i<end;++i) {
      bnum = sideord[i];
      
      tmpsb = sbdry[count];
      sbdry[count] = sbdry[bnum];
      sbdry[bnum] = tmpsb;
   
/*		KEEP SIDEORD UP TO DATE */
      for(j=0;j<nsbd;++j) {
         if (sideord[j] == count) {
            sideord[j] = bnum;
            break;
         }
      }
      sideord[i] = count;
      ++count;
   }
   
   for(i=end;i<nsbd;++i)
      sideord[bgn +i-end] = sideord[i];

/*	SHIFT INTERIOR BOUNDARIES */ 
   for(i=0;i<nsbd - (end-bgn);++i) {
      bnum = sideord[i];
      
      tmpsb = sbdry[count];
      sbdry[count] = sbdry[bnum];
      sbdry[bnum] = tmpsb;
      
/*		KEEP SIDEORD UP TO DATE */
      for(j=0;j<nsbd;++j) {
         if (sideord[j] == count) {
            sideord[j] = bnum;
            break;
         }
      }
      sideord[i] = count;
      ++count;
   }

/* UPDATE STRI/TTRI POINTERS ON BOUNDARY */   
   bdrylabel();
   
   return;
}

/* FIX STRI TTRI TO POINT TO GROUP/SIDE ON BOUNDARY */
void mesh::bdrylabel() {
   static int i,j,k,sind,tind;
   
   for(i=0;i<nsbd;++i) {
      for(j=0;j<sbdry[i].num;++j) {
         sind = sbdry[i].el[j];
         stri[sind][1] = -((i+1)*maxsbel +j);
         tind = stri[sind][0];
         for(k=0;k<3;++k)
            if (tside[tind].side[k] == sind) break;
            
         ttri[tind][k] = stri[sind][1];
      }
   }
   
   return;
}


void mesh::initvlngth() {
   int i,v0,v1;
   FLT l;
   
   for(i=0;i<nvrtx;++i)
      vlngth[i] = 0.0;
      
   for(i=0;i<nside;++i) {
      v0 = svrtx[i][0];
      v1 = svrtx[i][1];
      l = distance(svrtx[i][0],svrtx[i][1]);
      vlngth[v0] += l;
      vlngth[v1] += l;
   }
   
   for(i=0;i<nvrtx;++i)
      vlngth[i] /= nnbor[i];
   
   return;
}
      
   
   