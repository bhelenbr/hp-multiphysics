/*
 *  r_boundary.h
 *  mesh
 *
 *  Created by Brian Helenbrook on Mon Jun 10 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */
 
#include "mesh.h"

class r_side_bdry {
    protected:
        std::string mytype;
        r_mesh &x;
        side_bdry &base;

    public:
        r_side_bdry(r_mesh &xin, side_bdry &bin) : x(xin), base(bin) {mytype="plain";}
        /* VIRTUAL FUNCTIONS FOR BOUNDARY DEFORMATION */
        virtual void output(std::ostream& fout) {
            fout << base.idprefix << "_r_type: " << mytype << std::endl;            
        }
        virtual void input(input_map& bdrydata) {}
        virtual void tadvance() {}
        virtual void dirichlet() {}
        virtual void fixdx2() {}
        virtual ~r_side_bdry() {}
};

class r_fixed : public r_side_bdry {
    public:
        int dstart, dstop;
        
        r_fixed(r_mesh &xin, side_bdry &bin) : r_side_bdry(xin,bin), dstart(0), dstop(1) {mytype = "fixed";} 
            
        void input(input_map& inmap) {
            r_side_bdry::input(inmap);
            
            std::string line;
            inmap.getlinewdefault(base.idprefix+"_r_dir",line,"0 1");
            std::istringstream data(line);
            data >> dstart >> dstop;
            data.clear();
        }  
        
        void output(std::ostream& fout) {
            r_side_bdry::output(fout);
            fout << base.idprefix << "_r_dir: " << dstart << '\t' << dstop << std::endl;
        }
              
        void dirichlet() {            
            for(int j=0;j<base.nel;++j) {
                int sind = base.el(j);
                    for(int n=dstart;n<=dstop;++n) {
                        x.gbl_ptr->res(x.sd(sind).vrtx(0))(n) = 0.0;
                        x.gbl_ptr->res(x.sd(sind).vrtx(1))(n) = 0.0;
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
                int sind = base.el(j);
                    for(int n=d2start;n<=d2stop;++n) {
                        x.gbl_ptr->res1(x.sd(sind).vrtx(0))(n) = 0.0;
                        x.gbl_ptr->res1(x.sd(sind).vrtx(1))(n) = 0.0;
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
            for(int n=0;n<mesh::ND;++n)
                dx[n] = 0.0;
            mytype = "translating";
        }
        void input(input_map& inmap) {
            r_fixed::input(inmap);
            
            std::string line;
            if(inmap.getline(base.idprefix+"_r_translate",line)) {
                std::istringstream data(line);
                data >> dx[0] >> dx[1];
                data.clear();
            }
            else {
                dx[0] = 0;
                dx[1] = 1;
            }
        }  
        
        void output(std::ostream& fout) {
            r_fixed::output(fout);
            fout << base.idprefix << "_r_translate: " << dx[0] << '\t' << dx[1] << std::endl;
        }
        void tadvance() {
            int n,v0;
            for(int j=1;j<base.nel;++j) {
                v0 = x.sd(base.el(j)).vrtx(0);
                for(n=0;n<2;++n)
                    x.vrtx(v0)(n) += dx[n];
            }
            
            
            /* TEMPORARY */
                v0 = x.sd(base.el(0)).vrtx(0);
                for(n=0;n<2;++n)
                    x.vrtx(v0)(n) += dx[n]/2;
                    
                v0 = x.sd(base.el(base.nel-1)).vrtx(1);
                for(n=0;n<2;++n)
                    x.vrtx(v0)(n) += dx[n]/2;
                
        }
};

class r_oscillating : public r_fixed {
    public:
        FLT v0, amp, omega;
        
        r_oscillating(r_mesh &xin, side_bdry &bin) : 
            r_fixed(xin,bin), v0(0.0), amp(0.0), omega(0.0) {mytype="oscillating";}
            
        void input(input_map& inmap) {
            r_fixed::input(inmap);
            
            std::string line;
            if(inmap.getline(base.idprefix+"_r_oscillate",line)) {
                std::istringstream data(line);
                data >> v0 >> amp >> omega;
                data.clear();
            }
            else {
                v0 = 0.0;
                amp = 0.0;
                omega = 0.0;
            }
        }  
        
        void output(std::ostream& fout) {
            r_fixed::output(fout);
            fout << base.idprefix << "_r_oscillate: " << v0 << '\t' << amp << '\t' << omega << std::endl;
        }
            
        void tadvance() {
            int n,vrt;
            FLT center[2], center1[2], dx[2],xp[2];
            FLT theta, theta1,dtheta;
            
            center[0] = v0*sim::time;
            center[1] = amp*(1-cos(omega*sim::time));
            theta = atan(omega*amp*sin(omega*sim::time)/v0);            
            
            center1[0] = v0*sim::time;
            center1[1] = amp*(1-cos(omega*sim::time));
            theta1 = atan(omega*amp*sin(omega*sim::time)/v0);
            dtheta = theta1-theta;
                        
            for(int j=0;j<base.nel;++j) {
                vrt = x.sd(base.el(j)).vrtx(0);
                xp[0] = x.vrtx(vrt)(0)-center[0];
                xp[1] = x.vrtx(vrt)(1)-center[1];            
                dx[0] = center1[0]-center[0] +xp[0]*cos(dtheta)-xp[1]*sin(dtheta) -xp[0];
                dx[1] = center1[1]-center[1] +xp[0]*sin(dtheta)+xp[1]*cos(dtheta) -xp[1];
                for(n=0;n<2;++n)
                    x.vrtx(vrt)(n) += dx[n];
            }
        }
};

