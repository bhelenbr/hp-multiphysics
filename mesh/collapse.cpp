#include "mesh.h"
#include "utilities.h"
#include<stdio.h>
#include<assert.h>
#include<float.h>

/*	THESE TELL WHICH TRIS WERE AFFECTED & WHICH SIDES WERE DELETED */
extern int ntdel, tdel[MAXLST +1];
extern int nsdel, sdel[MAXLST +1];

int mesh::collapse(int sind) {
   int i,j,vn,vnear,prev,tind,tind1,stoptri,dir[2];
   int ntsrnd[2],tsrnd[2][MAXLST],nssrnd[2],ssrnd[2][MAXLST];
   int delt,v0,v1,sd,pt,sd1,sd2,t1,t2;
   FLT x,y,a,asum,dx,dy,l0,l1;
   
/*	FIND TRIANGLES / SIDES SURROUNDING BOTH ENDPOINTS */
   for(i=0;i<2;++i) {
      vnear = svrtx[sind][i];
      tind = vtri[vnear];
      if (tind != stri[sind][0] && tind != stri[sind][1])
         prev = 0;
      else 
         prev = 1;
      stoptri = tind;
      dir[i] = 1;
      ntsrnd[i] = 0;
      nssrnd[i] = 0;
      do {
         for(vn=0;vn<3;++vn) 
            if (tvrtx[tind][vn] == vnear) break;
         assert(vn != 3);
         
         tind1 = ttri[tind][(vn +dir[i])%3];
         if (tind1 < 0) {
            if (dir[i] > 1) break;
/*				REVERSE DIRECTION AND GO BACK TO START */
            ++dir[i];
            tind1 = vtri[vnear];
            prev = 1;
            stoptri = -1;
         }
         
         
         if (tind1 != stri[sind][0] && tind1 != stri[sind][1]) {
            tsrnd[i][ntsrnd[i]++] = tind1;
            if (!prev) {
               ssrnd[i][nssrnd[i]++] = tside[tind].side[(vn +dir[i])%3];
            }
            assert(ntsrnd[i] < MAXLST);
            assert(nssrnd[i] < MAXLST);
            prev = 0;
         }
         else {
            prev = 1;
         }

         tind = tind1;

      } while(tind != stoptri); 
   }

/*	DECIDE WHICH POINT TO COLLAPSE TOWARDS */
/*	REGARDLESS OF OUTCOME THIS SIDE SHOULD BE REMOVED FROM LIST */
   nsdel = 0;
   sdel[nsdel++] = sind;

   if (dir[0] + dir[1] > 2) {
/*		ONE OR BOTH POINTS ON BOUNDARY  */
      if (dir[0] + dir[1] == 3) {
/* 		ONE POINT ON BOUNDARY */
         if (vinfo[svrtx[sind][0]] == 1)
            delt = 1;
         else if (vinfo[svrtx[sind][1]] == 1) 
            delt = 0;
         else if (dir[0] == 1)
            delt = 0;
         else if (dir[1] == 1)
            delt = 1;
      }
      else {
/* 		BOTH ON BOUNDARY */
/*			IF NOT BOUNDARY EDGE OR TWO ENDPOINTS RETURN */
         if (stri[sind][1] > -1 || vinfo[svrtx[sind][0]] +vinfo[svrtx[sind][1]] == 2) return(1);
         if (vinfo[svrtx[sind][0]] == 1)
            delt = 1;
         else if (vinfo[svrtx[sind][1]] == 1) 
            delt = 0;
         else
            goto SEARCH;
      }
      goto DELETE;
   }

SEARCH:

/*	KEEP POINT WHICH IS CLOSEST TO CENTER OF AREA */
   x = 0.0;
   y = 0.0;
   asum = 0.0;
   for(i=0;i<2;++i) {
      tind = stri[sind][i];
      if (tind > -1) {
         a = area(tind);
         asum += a;
         for(vn=0;vn<3;++vn) {
            x += a*vrtx[tvrtx[tind][vn]][0];
            y += a*vrtx[tvrtx[tind][vn]][1];
         }
      }            
      for(j=0;j<ntsrnd[i];++j) {
         tind = tsrnd[i][j];
         a = area(tind);
         asum += a;
         for(vn=0;vn<3;++vn) {
            x += a*vrtx[tvrtx[tind][vn]][0];
            y += a*vrtx[tvrtx[tind][vn]][1];
         }
      }
   }
   asum = 1./(3.*asum);
   x = x*asum;
   y = y*asum;
   v0 = svrtx[sind][0];
   v1 = svrtx[sind][1];
   dx = vrtx[v0][0] -x;
   dy = vrtx[v0][1] -y;
   l0 = dx*dx +dy*dy;   
   dx = vrtx[v1][0] -x;
   dy = vrtx[v1][1] -y;
   l1 = dx*dx +dy*dy;
   
   delt = (l0 > l1 ? 0 : 1);

DELETE:
   sinfo[sind] = -3;
   
/*	UPDATE TVRTX & SVRTX */
   v0 = svrtx[sind][1-delt];
   v1 = svrtx[sind][delt];
   
   for(i=0;i<ntsrnd[delt];++i) {
      tind = tsrnd[delt][i];
      for(j=0;j<3;++j) {
         if (tvrtx[tind][j] == v1) {
            tvrtx[tind][j] = v0;
            sd = tside[tind].side[(j+1)%3];
            pt = (1 +tside[tind].sign[(j+1)%3])/2;
            assert(svrtx[sd][pt] == v1 || svrtx[sd][pt] == v0);
            svrtx[sd][pt] = v0;
            sd = tside[tind].side[(j+2)%3];
            pt = (1 -tside[tind].sign[(j+2)%3])/2;
            assert(svrtx[sd][pt] == v1 || svrtx[sd][pt] == v0);
            svrtx[sd][pt] = v0;           
            break;
         }
      }
      assert(j != 3);
   }

/*	CLOSE THE GAP */
   for(i=0;i<2;++i) {
      tind = stri[sind][i];
      if (tind < 0) continue;

/*		MARK TRI AS DELETED TO BE REMOVED LATER */
      tinfo[tind] = -1;  
      
      for(sd1=0;sd1<3;++sd1)
         if(tvrtx[tind][sd1] == v1) break;
         
      assert(sd1 != 3);
               
      if (tside[tind].side[(sd1+1)%3] == sind)
         sd2 = (sd1+2)%3;
      else
         sd2 = (sd1+1)%3;
         
/*		MARK SIDE AS DELETED TO BE REMOVED LATER */
      sd = tside[tind].side[sd2];
      sinfo[sd] = -3;
      sdel[nsdel++] = sd;
            
      t1 = ttri[tind][sd1];
      t2 = ttri[tind][sd2];

/*		UPDATE TTRI FOR T1 */  
      if (t1 > -1) {
         for(j=0;j<3;++j) {
            if(ttri[t1][j] == tind) {
               ttri[t1][j] = t2;
               break;
            }
         }
         assert(j != 3);
         vtri[tvrtx[tind][(sd1+1)%3]] = t1;
         vtri[tvrtx[tind][(sd1+2)%3]] = t1;
      }
/*		UPDATE STRI FOR KEPT SIDE */
      sd = tside[tind].side[sd1];
      pt = (1 -tside[tind].sign[sd1])/2;
      stri[sd][pt] = t2;

/*		UPDATE TTRI/TSIDE FOR T2 */
      if (t2 > -1) {
         for(j=0;j<3;++j) {
            if(ttri[t2][j] == tind) {
               ttri[t2][j] = t1;
               tside[t2].side[j] = sd;
               tside[t2].sign[j] = tside[tind].sign[sd1];
               break;
            }
         }
         assert(j != 3);
         vtri[tvrtx[tind][(sd1+1)%3]] = t2;
         vtri[tvrtx[tind][(sd1+2)%3]] = t2;
      }
   }
   
/* NEED TO REMOVE LEFTOVERS */
   qtree.dltpt(v1);
   vinfo[v1] = -1;
   
/*	STORE LIST OF AFFECTED TRIANGLES */
   ntdel = ntsrnd[delt];
   for(i=0;i<ntsrnd[delt];++i)
      tdel[i] = tsrnd[delt][i];

/* MARK AFFECTED SIDES AS TOUCHED */
   for(i=0;i<nssrnd[delt];++i)
      sinfo[ssrnd[delt][i]] = -2;
      
/*	SWAP AFFECTED SIDES */      
   swap(nssrnd[delt],ssrnd[delt]);

   return(0);
}

/* DELETE UNREFERENCED TRIANGLE */
void mesh::dlttri(int tind) {
   int i,j,v0,t1,sind,flip;
   
   --ntri;
   if (tind == ntri) return;
   
   for (j=0;j<3;++j) {
      v0 = tvrtx[ntri][j];
      tvrtx[tind][j] = v0;
      vtri[v0] = tind;
      tside[tind].side[j] = tside[ntri].side[j];
      tside[tind].sign[j] = tside[ntri].sign[j];
      sind = tside[ntri].side[j];
      flip = (1 -tside[ntri].sign[j])/2;
      stri[sind][flip] = tind;
      
      t1 = ttri[ntri][j];
      ttri[tind][j] = t1;
      if (t1 > -1) {
         for(i=0;i<3;++i) {
            if(ttri[t1][i] == ntri) {
               ttri[t1][i] = tind;
               break;
            }
         }
      }
   }
   
   return;
}

/* DELETE UNREFERENCED SIDE */
void mesh::dltside(int sind) {
   int j,k,tind;
   
/*	DELETE SIDE */
   --nside;
   
   if (sind == nside) return;
   
   svrtx[sind][0] = svrtx[nside][0];
   svrtx[sind][1] = svrtx[nside][1];
   stri[sind][0] = stri[nside][0];
   stri[sind][1] = stri[nside][1];

   for(j=0;j<2;++j) {
      tind = stri[nside][j];
      if (tind > -1) {
         for(k=0;k<3;++k)
            if (tside[tind].side[k] == nside) break;
         assert(k != 3);
         tside[tind].side[k] = sind;
      }
   }
   
   return;
}

/* DELETE UNREFERENCED VERTEX */
void mesh::dltvrtx(int v0) {
   int vn,stoptri,dir;
   int tind, sind, flip;
      
   --nvrtx;
   
   if (v0 == nvrtx) return;
   
   qtree.movept(nvrtx,v0);
   
   vrtx[v0][0] = vrtx[nvrtx][0];
   vrtx[v0][1] = vrtx[nvrtx][1];
   vtri[v0] = vtri[nvrtx];
   vlngth[v0] = vlngth[nvrtx];
   
   tind = vtri[nvrtx];   
   stoptri = tind;
   dir = 1;
   do {
      for(vn=0;vn<3;++vn) 
         if (tvrtx[tind][vn] == nvrtx) break; 
         
      assert(vn != 3);
      
      tvrtx[tind][vn] = v0;
      sind = tside[tind].side[(vn +dir)%3];
      flip = (1 +(3-2*dir)*tside[tind].sign[(vn +dir)%3])/2;
      assert(svrtx[sind][flip] == nvrtx);
      svrtx[sind][flip] = v0;

      tind = ttri[tind][(vn +dir)%3];
      if (tind < 0) {
         if (dir > 1) break;
/*			REVERSE DIRECTION AND GO BACK TO START */
         ++dir;
         tind = vtri[nvrtx];
         for(vn=0;vn<3;++vn) {
            if (tvrtx[tind][vn] == v0) {
               tvrtx[tind][vn] = nvrtx;
               break;
            }
         }
         stoptri = -1;
      }
      
   } while(tind != stoptri); 

   return;
}


