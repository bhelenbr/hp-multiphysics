#include "mesh.h"
#include "boundary.h"
#include <utilities.h>
#include <assert.h>
#include <float.h>

/* THESE TELL WHICH TRIS WERE AFFECTED & WHICH SIDES WERE DELETED */
/* ntdel, tdel[maxlst +1]; */
/* nsdel, sdel[maxlst +1]; */

int mesh::collapse(int sind,int& ntdel,int *tdel,int& nsdel,int *sdel) {
   int i,j,vn,vnear,prev,tind,tind1,sind1,stoptri,dir[2],bnum;
   int ntsrnd[2],tsrnd[2][maxlst],nssrnd[2],ssrnd[2][maxlst],bside[2][2];
   int delt,v0,v1,pt,sd1,sd2,sd3,t1,t2;
   FLT x,y,a,asum,dx,dy,l0,l1;

   
   /* FIND TRIANGLES / SIDES SURROUNDING BOTH ENDPOINTS */
   /* EXCLUDES TRIANGLES ADJACENT TO DELETED SIDE */
   for(i=0;i<2;++i) {
      vnear = sd(sind).vrtx(i);
      tind = vd(vnear).tri;
      if (tind != sd(sind).tri(0) && tind != sd(sind).tri(1))
         prev = 0;
      else 
         prev = 1;
      stoptri = tind;
      dir[i] = 1;
      ntsrnd[i] = 0;
      nssrnd[i] = 0;
      do {
         for(vn=0;vn<3;++vn) 
            if (td(tind).vrtx(vn) == vnear) break;
         assert(vn != 3);
         
         tind1 = td(tind).tri((vn +dir[i])%3);
         if (tind1 < 0) {
            bside[i][dir[i]-1] = td(tind).side((vn +dir[i])%3);
            if (dir[i] > 1) break;
            /* REVERSE DIRECTION AND GO BACK TO START */
            ++dir[i];
            tind1 = vd(vnear).tri;
            prev = 1;
            stoptri = -1;
         }
         
         
         if (tind1 != sd(sind).tri(0) && tind1 != sd(sind).tri(1)) {
            tsrnd[i][ntsrnd[i]++] = tind1;
            if (!prev) {
               ssrnd[i][nssrnd[i]++] = td(tind).side((vn +dir[i])%3);
            }
            assert(ntsrnd[i] < maxlst);
            assert(nssrnd[i] < maxlst);
            prev = 0;
         }
         else {
            prev = 1;
         }

         tind = tind1;

      } while(tind != stoptri); 
   }

   /* DECIDE WHICH POINT TO COLLAPSE TOWARDS */
   /* REGARDLESS OF OUTCOME THIS SIDE SHOULD BE REMOVED FROM LIST */
   nsdel = 0;
   sdel[nsdel++] = sind;

   if (dir[0] + dir[1] > 2) {
      /* ONE OR BOTH POINTS ON BOUNDARY  */
      if (dir[0] + dir[1] == 3) {
         /* ONE POINT ON BOUNDARY */
         if (dir[0] == 1)
            delt = 0;
         else if (dir[1] == 1)
            delt = 1;
      }
      else {         
         /* BOTH ON BOUNDARY */
         /* IF NOT BOUNDARY EDGE OR TWO ENDPOINTS RETURN */
         if (sd(sind).tri(1) > -1 || vd(sd(sind).vrtx(0)).info +vd(sd(sind).vrtx(1)).info == -4) return(-1);

         /* CHECK IF CORNER POINT */
         if (vd(sd(sind).vrtx(0)).info == -2) {
            delt = 1;
            goto DELETE;
         }
         if (vd(sd(sind).vrtx(1)).info == -2) {
            delt = 0;
            goto DELETE;
         }

         /* CHECK IF TWO EDGES OF TRIANGLE ARE ON BOUNDARY */         
         tind = sd(sind).tri(0);
         for (i=0;i<3;++i)
            if(td(tind).side(i) == sind) break;
         assert(i != 3);
         if (td(tind).tri((i+1)%3) < 0) {
            delt = 0;
            goto DELETE;
         }
         if (td(tind).tri((i+2)%3) < 0) {
            delt = 1;
            goto DELETE;
         }
         
         /* CAN PICK EITHER POINT KEEP ONE CLOSEST TO CENTER OF PREVIOUS/NEXT VERTICES ON BOUNDARY */
         /* KEEP COMMUNICATION BOUNDARIES COHERENT */
         sd1 = bside[0][0];
         sd2 = bside[1][1];
         
         v0 = sd(sd1).vrtx(0);
         v1 = sd(sd2).vrtx(1);
         x = 0.5*(vrtx(v0)(0) +vrtx(v1)(0));
         y = 0.5*(vrtx(v0)(1) +vrtx(v1)(1));
         
         v0 = sd(sind).vrtx(0);
         dx = vrtx(v0)(0)-x;
         dy = vrtx(v0)(1)-y;
         l0 = dx*dx +dy*dy;
         
         v0 = sd(sind).vrtx(1);
         dx = vrtx(v0)(0)-x;
         dy = vrtx(v0)(1)-y;
         l1 = dx*dx +dy*dy;
         delt = (l0 < l1 ? 1 : 0);

         if (fabs(l0 - l1)/(l0 +l1) < 10.*EPSILON) {
            /* CONSISTENT WAY TO PICK IN DEGENERATE CASE? */
            bnum = getbdrynum(sd(sind).tri(1));
            delt = sbdry(bnum)->is_frst();
               *log << "#Warning: degenerate case in edge collapse for bdry" << bnum << sbdry(bnum)->idnum << "picking "
                << delt << "(" << vrtx(sd(sind).vrtx(delt))(0) << "," << vrtx(sd(sind).vrtx(delt))(1) << ")" << std::endl;
         }
      }
   }
   else {
      /* THIS IS AN INTERIOR EDGE WITH NO CONNECTION TO BOUNDARY */
      /* KEEP POINT WHICH IS CLOSEST TO CENTER OF AREA */
      x = 0.0;
      y = 0.0;
      asum = 0.0;
      for(i=0;i<2;++i) {
         tind = sd(sind).tri(i);
         if (tind > -1) {
            a = area(tind);
            asum += a;
            for(vn=0;vn<3;++vn) {
               x += a*vrtx(td(tind).vrtx(vn))(0);
               y += a*vrtx(td(tind).vrtx(vn))(1);
            }
         }            
         for(j=0;j<ntsrnd[i];++j) {
            tind = tsrnd[i][j];
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
      
      delt = (l0 > l1 ? 0 : 1);
   }
   
DELETE:
   
   /* UPDATE TVRTX & SVRTX */
   v0 = sd(sind).vrtx(1-delt);
   v1 = sd(sind).vrtx(delt);
   
   for(i=0;i<ntsrnd[delt];++i) {
      tind = tsrnd[delt][i];
      for(j=0;j<3;++j) {
         if (td(tind).vrtx(j) == v1) {
            td(tind).vrtx(j) = v0;
            sd3 = td(tind).side((j+1)%3);
            sd(sd3).info = -2;
            pt = (1 +td(tind).sign((j+1)%3))/2;
            assert(sd(sd3).vrtx(pt) == v1 || sd(sd3).vrtx(pt) == v0);
            sd(sd3).vrtx(pt) = v0;
            sd3 = td(tind).side((j+2)%3);
            sd(sd3).info = -2;
            pt = (1 -td(tind).sign((j+2)%3))/2;
            assert(sd(sd3).vrtx(pt) == v1 || sd(sd3).vrtx(pt) == v0);
            sd(sd3).vrtx(pt) = v0;
            break;
         }
      }
      assert(j != 3);
   }

   /* MARK SIDE AS DELETED */      
   sd(sind).info = -3;

   /* CLOSE THE GAP */
   for(i=0;i<2;++i) {
      tind = sd(sind).tri(i);
      if (tind < 0) continue;

      /* MARK TRI AS DELETED TO BE REMOVED LATER */
      td(tind).info = -3;  
      
      for(sd1=0;sd1<3;++sd1)
         if(td(tind).vrtx(sd1) == v1) break;
         
      assert(sd1 != 3);
               
      if (td(tind).side((sd1+1)%3) == sind)
         sd2 = (sd1+2)%3;
      else
         sd2 = (sd1+1)%3;
         
      /* MARK SIDE AS DELETED TO BE REMOVED LATER */
      
      sind1 = td(tind).side(sd2);
      if (sd(sind1).tri(1) < 0) {
         /* BOUNDARY SIDE SO DELETE OTHER ONE */
         /* ANOMALY OF TRIANGLE WITH 2 EDGES ON BOUNDARY */
         sd(sind1).info = -2;
         if (sd(sind1).vrtx(0) == v1) 
            sd(sind1).vrtx(0) = v0;
         if (sd(sind1).vrtx(1) == v1) 
            sd(sind1).vrtx(1) = v0;   
         sind1 = sd2;
         sd2 = sd1;
         sd1 = sind1;
         sind1 = td(tind).side(sd2);
      }
      sd(sind1).info = -3;
      sdel[nsdel++] = sind1;
            
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
         assert(j != 3);
         for(j=0;j<3;++j)
            vd(td(tind).vrtx(j)).tri = t1;
      }
      /* UPDATE STRI FOR KEPT SIDE */
      sind1 = td(tind).side(sd1);
      pt = (1 -td(tind).sign(sd1))/2;
      sd(sind1).tri(pt) = t2;

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
         assert(j != 3);
         for(j=0;j<3;++j)
            vd(td(tind).vrtx(j)).tri = t2;
      }
   }
   
   /* NEED TO REMOVE LEFTOVERS */
   qtree.dltpt(v1);
   vd(v1).info = -3;
   
   /* STORE LIST OF AFFECTED TRIANGLES */
   ntdel = ntsrnd[delt];
   for(i=0;i<ntsrnd[delt];++i)
      tdel[i] = tsrnd[delt][i];
      
   /* SWAP AFFECTED SIDES */      
   swap(nssrnd[delt],ssrnd[delt]);

   return(delt);
}

int mesh::collapse1(int sind,int delt,int& ntdel,int *tdel,int& nsdel,int *sdel) {
   int i,j,vn,vnear,prev,tind,tind1,sind1,stoptri,dir[2];
   int ntsrnd[2],tsrnd[2][maxlst],nssrnd[2],ssrnd[2][maxlst],bside[2][2];
   int v0,v1,pt,sd1,sd2,sd3,t1,t2;

   
   /* FIND TRIANGLES / SIDES SURROUNDING BOTH ENDPOINTS */
   /* EXCLUDES TRIANGLES ADJACENT TO DELETED SIDE */
   for(i=0;i<2;++i) {
      vnear = sd(sind).vrtx(i);
      tind = vd(vnear).tri;
      if (tind != sd(sind).tri(0) && tind != sd(sind).tri(1))
         prev = 0;
      else 
         prev = 1;
      stoptri = tind;
      dir[i] = 1;
      ntsrnd[i] = 0;
      nssrnd[i] = 0;
      do {
         for(vn=0;vn<3;++vn) 
            if (td(tind).vrtx(vn) == vnear) break;
         assert(vn != 3);
         
         tind1 = td(tind).tri((vn +dir[i])%3);
         if (tind1 < 0) {
            bside[i][dir[i]-1] = td(tind).side((vn +dir[i])%3);
            if (dir[i] > 1) break;
            /* REVERSE DIRECTION AND GO BACK TO START */
            ++dir[i];
            tind1 = vd(vnear).tri;
            prev = 1;
            stoptri = -1;
         }
         
         
         if (tind1 != sd(sind).tri(0) && tind1 != sd(sind).tri(1)) {
            tsrnd[i][ntsrnd[i]++] = tind1;
            if (!prev) {
               ssrnd[i][nssrnd[i]++] = td(tind).side((vn +dir[i])%3);
            }
            assert(ntsrnd[i] < maxlst);
            assert(nssrnd[i] < maxlst);
            prev = 0;
         }
         else {
            prev = 1;
         }

         tind = tind1;

      } while(tind != stoptri); 
   }

   /* DECIDE WHICH POINT TO COLLAPSE TOWARDS */
   /* REGARDLESS OF OUTCOME THIS SIDE SHOULD BE REMOVED FROM LIST */
   nsdel = 0;
   sdel[nsdel++] = sind;
   
   if (dir[0] + dir[1] == 4) {
      /* BOTH POINTS OF EDGE ARE ON BOUNDARY */
      if (sd(sind).tri(1) > -1 || vd(sd(sind).vrtx(0)).info +vd(sd(sind).vrtx(1)).info == -4) return(-1);
   }
   
   /* UPDATE TVRTX & SVRTX */
   v0 = sd(sind).vrtx(1-delt);
   v1 = sd(sind).vrtx(delt);
   
   for(i=0;i<ntsrnd[delt];++i) {
      tind = tsrnd[delt][i];
      for(j=0;j<3;++j) {
         if (td(tind).vrtx(j) == v1) {
            td(tind).vrtx(j) = v0;
            sd3 = td(tind).side((j+1)%3);
            sd(sd3).info = -2;
            pt = (1 +td(tind).sign((j+1)%3))/2;
            assert(sd(sd3).vrtx(pt) == v1 || sd(sd3).vrtx(pt) == v0);
            sd(sd3).vrtx(pt) = v0;
            sd3 = td(tind).side((j+2)%3);
            sd(sd3).info = -2;
            pt = (1 -td(tind).sign((j+2)%3))/2;
            assert(sd(sd3).vrtx(pt) == v1 || sd(sd3).vrtx(pt) == v0);
            sd(sd3).vrtx(pt) = v0;
            break;
         }
      }
      assert(j != 3);
   }

   /* MARK SIDE AS DELETED */      
   sd(sind).info = -3;

   /* CLOSE THE GAP */
   for(i=0;i<2;++i) {
      tind = sd(sind).tri(i);
      if (tind < 0) continue;

      /* MARK TRI AS DELETED TO BE REMOVED LATER */
      td(tind).info = -3;  
      
      for(sd1=0;sd1<3;++sd1)
         if(td(tind).vrtx(sd1) == v1) break;
         
      assert(sd1 != 3);
               
      if (td(tind).side((sd1+1)%3) == sind)
         sd2 = (sd1+2)%3;
      else
         sd2 = (sd1+1)%3;
         
      /* MARK SIDE AS DELETED TO BE REMOVED LATER */
      sind1 = td(tind).side(sd2);
      if (sd(sind1).tri(1) < 0) {
         /* BOUNDARY SIDE SO DELETE OTHER ONE */
         /* ANOMALY OF TRIANGLE WITH 2 EDGES ON BOUNDARY */
         sd(sind1).info = -2;
         if (sd(sind1).vrtx(0) == v1) 
            sd(sind1).vrtx(0) = v0;
         if (sd(sind1).vrtx(1) == v1) 
            sd(sind1).vrtx(1) = v0;   
         sind1 = sd2;
         sd2 = sd1;
         sd1 = sind1;
         sind1 = td(tind).side(sd2);
      }
      sd(sind1).info = -3;
      sdel[nsdel++] = sind1;
            
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
         assert(j != 3);
         for(j=0;j<3;++j)
            vd(td(tind).vrtx(j)).tri = t1;
      }
      /* UPDATE STRI FOR KEPT SIDE */
      sind1 = td(tind).side(sd1);
      pt = (1 -td(tind).sign(sd1))/2;
      sd(sind1).tri(pt) = t2;

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
         assert(j != 3);
         for(j=0;j<3;++j)
            vd(td(tind).vrtx(j)).tri = t2;
      }
   }

   
   /* NEED TO REMOVE LEFTOVERS */
   qtree.dltpt(v1);
   vd(v1).info = -3;
   
   /* STORE LIST OF AFFECTED TRIANGLES */
   ntdel = ntsrnd[delt];
   for(i=0;i<ntsrnd[delt];++i)
      tdel[i] = tsrnd[delt][i];
      
   /* SWAP AFFECTED SIDES */      
   swap(nssrnd[delt],ssrnd[delt]);

   return(delt);
}









/* DELETE UNREFERENCED TRIANGLE */
void mesh::dlttri(int tind) {
   int i,j,v0,t1,sind,flip;
   
   
   /* MOVE UP FROM BOTTOM UNTIL FIND UNDELETED TRI */
   while(td(ntri-1).info == -3)
      --ntri;
   
   if (tind >= ntri) return;
   
   --ntri;

   td(tind).info = ntri;
   td(ntri).info = tind;
   
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
   
   return;
}

/* DELETE UNREFERENCED SIDE */
void mesh::dltsd(int sind) {
   int j,k,tind;

   /* MOVE UP FROM BOTTOM UNTIL FIND UNDELETED SIDE */
   while(sd(nside-1).info == -3)
      --nside;
   
   if (sind >= nside) return;
   
   /* DELETE SIDE */
   --nside;
   
   sd(sind).vrtx(0) = sd(nside).vrtx(0);
   sd(sind).vrtx(1) = sd(nside).vrtx(1);
   sd(sind).tri(0) = sd(nside).tri(0);
   sd(sind).tri(1) = sd(nside).tri(1);
   if (sd(nside).info == -1) {
      sd(sind).info = nside;
   }
   else {
      sd(sind).info = -2;
   }
   sd(nside).info = sind;

   for(j=0;j<2;++j) {
      tind = sd(nside).tri(j);
      if (tind > -1) {
         for(k=0;k<3;++k)
            if (td(tind).side(k) == nside) break;
         td(tind).side(k) = sind;
      }
   }
   
   return;
}

/* DELETE UNREFERENCED VERTEX */
void mesh::dltvrtx(int v0) {
   int vn,stoptri,dir;
   int tind, sind, flip;
         
   /* MOVE UP FROM BOTTOM UNTIL FIND UNDELETED VERTEX */
   while(vd(nvrtx-1).info == -3)
      --nvrtx;
      
   if (v0 >= nvrtx) return;   
   
   --nvrtx;

   vd(v0).info = nvrtx;
   vd(nvrtx).info = v0;
   
   qtree.movept(nvrtx,v0);
   
   vrtx(v0)(0) = vrtx(nvrtx)(0);
   vrtx(v0)(1) = vrtx(nvrtx)(1);
   vd(v0).nnbor = vd(nvrtx).nnbor;
   vd(v0).tri = vd(nvrtx).tri;
   vlngth(v0) = vlngth(nvrtx);
   
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

   return;
}


