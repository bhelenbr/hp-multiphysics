/*
 *  r_boundary.h
 *  mesh
 *
 *  Created by Brian Helenbrook on Mon Jun 10 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */
 #include "boundary.h"

class r_side_bdry {
   protected:
      std::string mytype;
      r_mesh &x;
      side_bdry &base;

   public:
      r_side_bdry(r_mesh &xin, side_bdry &bin) : x(xin), base(bin) {mytype="plain";}
      /* VIRTUAL FUNCTIONS FOR BOUNDARY DEFORMATION */
      virtual void output(std::ostream& fout) {
         fout << base.idprefix << ".r_type: " << mytype << std::endl;         
      }
      virtual void input(std::map<std::string,std::string>& bdrydata) {}
      virtual void tadvance() {}
      virtual void dirichlet() {}
      virtual void fixdx2() {}
      virtual ~r_side_bdry() {}
};

class r_fixed : public r_side_bdry {
   public:
      int dstart, dstop;
      
      r_fixed(r_mesh &xin, side_bdry &bin) : r_side_bdry(xin,bin), dstart(0), dstop(1) {mytype = "fixed";} 
         
      void input(std::map <std::string,std::string>& inmap) {
         r_side_bdry::input(inmap);
         
         std::istringstream data(inmap[base.idprefix+".r_dirfixed"]);
         if (!(data >> dstart >> dstop)) {dstart = 0; dstop = 1;}
         data.clear();
      }  
      
      void output(std::ostream& fout) {
         r_side_bdry::output(fout);
         fout << base.idprefix << ".r_dirfixed: " << dstart << '\t' << dstop << std::endl;
      }
           
      void dirichlet() {         
         for(int j=0;j<base.nel;++j) {
            int sind = base.el[j];
               for(int n=dstart;n<=dstop;++n) {
                  x.res[x.sd[sind].vrtx[0]][n] = 0.0;
                  x.res[x.sd[sind].vrtx[1]][n] = 0.0;
               }
         }
         return;
      }
};

class r_fixed4 : public r_fixed {
   public:
      int d2start, d2stop;

      r_fixed4(r_mesh &xin, side_bdry &bin, int dstr, int dstp, int d2str, int d2stp)
         : r_fixed(xin,bin), d2start(d2str), d2stop(d2stp) {mytype="fixed4";} 
      
      void fixdx2() {
         for(int j=0;j<base.nel;++j) {
            int sind = base.el[j];
               for(int n=d2start;n<=d2stop;++n) {
                  x.work[x.sd[sind].vrtx[0]][n] = 0.0;
                  x.work[x.sd[sind].vrtx[1]][n] = 0.0;
               }
         }
         return;
      }
};

class r_translating : public r_fixed {
   private:
      FLT dx[2];
   public:
      r_translating(r_mesh &xin, side_bdry &bin): r_fixed(xin,bin) { 
         for(int n=0;n<mesh::DIM;++n)
            dx[n] = 0.0;
         mytype = "translating";
      }
      void input(std::map <std::string,std::string>& inmap) {
         r_fixed::input(inmap);
         
         std::istringstream data(inmap[base.idprefix+".r_translate"]);
         if (!(data >> dx[0] >> dx[1])) {dx[0] = 0; dx[1] = 1;}
         data.clear();
      }  
      
      void output(std::ostream& fout) {
         r_fixed::output(fout);
         fout << base.idprefix << ".r_translate: " << dx[0] << '\t' << dx[1] << std::endl;
      }
      void tadvance() {
         int n,v0;
         for(int j=0;j<base.nel;++j) {
            v0 = x.sd[base.el[j]].vrtx[0];
            for(n=0;n<2;++n)
               x.vrtx[v0][n] += dx[n];
         }
      }
};

class r_oscillating : public r_fixed {
   public:
      FLT v0, amp, omega;
      FLT time;
      
      r_oscillating(r_mesh &xin, side_bdry &bin) : 
         r_fixed(xin,bin), v0(0.0), amp(0.0), omega(0.0), time(0.0) {mytype="oscillating";}
         
      void input(std::map <std::string,std::string>& inmap) {
         r_fixed::input(inmap);
         
         std::istringstream data(inmap[base.idprefix+".r_oscillate"]);
         if (!(data >> v0 >> amp >> omega)) {v0 = 0.0; amp = 0.0; omega = 0.0;}
         data.clear();
      }  
      
      void output(std::ostream& fout) {
         r_fixed::output(fout);
         fout << base.idprefix << ".r_oscillate: " << v0 << '\t' << amp << '\t' << omega << std::endl;
      }
         
      void tadvance() {
         int n,vrt;
         FLT center[2], center1[2], dx[2],xp[2];
         FLT theta, theta1,dtheta;
         
         center[0] = v0*time;
         center[1] = amp*(1-cos(omega*time));
         theta = atan(omega*amp*sin(omega*time)/v0);         

         time = time +1.0;
         
         center1[0] = v0*time;
         center1[1] = amp*(1-cos(omega*time));
         theta1 = atan(omega*amp*sin(omega*time)/v0);
         dtheta = theta1-theta;
                  
         for(int j=0;j<base.nel;++j) {
            vrt = x.sd[base.el[j]].vrtx[0];
            xp[0] = x.vrtx[vrt][0]-center[0];
            xp[1] = x.vrtx[vrt][1]-center[1];         
            dx[0] = center1[0]-center[0] +xp[0]*cos(dtheta)-xp[1]*sin(dtheta) -xp[0];
            dx[1] = center1[1]-center[1] +xp[0]*sin(dtheta)+xp[1]*cos(dtheta) -xp[1];
            for(n=0;n<2;++n)
               x.vrtx[vrt][n] += dx[n];
         }
      }
};

class r_stype {
   public:
      static const int ntypes = 4;
      enum ids {plain=1, fixed, translating, oscillating};
      const static char names[ntypes][40];
      static int getid(const char *nin) {
         for(int i=0;i<ntypes;++i) 
            if (!strcmp(nin,names[i])) return(i+1);
         return(-1);
      }
};



