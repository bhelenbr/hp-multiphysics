/*
 *  r_tri_boundary.h
 *  mesh
 *
 *  Created by Brian Helenbrook on Mon Jun 10 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */
 
#include "tri_mesh.h"

class r_side_bdry {
    protected:
        std::string mytype;
        r_tri_mesh &x;
        edge_bdry &base;

    public:
        r_side_bdry(r_tri_mesh &xin, edge_bdry &bin) : x(xin), base(bin) {mytype="plain";}
        r_side_bdry(const r_side_bdry &inbdry, r_tri_mesh &xin, edge_bdry &bin) : x(xin), base(bin) {mytype="plain";}
        virtual r_side_bdry* create(r_tri_mesh &xin, edge_bdry &bin) const {return new r_side_bdry(*this,xin,bin);}
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
        
        r_fixed(r_tri_mesh &xin, edge_bdry &bin) : r_side_bdry(xin,bin), dstart(0), dstop(1) {mytype = "fixed";} 
        r_fixed(const r_fixed &inbdry, r_tri_mesh &xin, edge_bdry &bin) : r_side_bdry(inbdry,xin,bin), dstart(inbdry.dstart), dstop(inbdry.dstop) {mytype = "fixed";}
        r_fixed* create(r_tri_mesh &xin, edge_bdry &bin) const {return new r_fixed(*this,xin,bin);}
                   
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
            for(int j=0;j<base.nseg;++j) {
                int sind = base.seg(j);
                    for(int n=dstart;n<=dstop;++n) {
                        x.gbl->res(x.seg(sind).pnt(0))(n) = 0.0;
                        x.gbl->res(x.seg(sind).pnt(1))(n) = 0.0;
                    }
            }
            return;
        }
};

class r_fixed4 : public r_fixed {
    public:
        int d2start, d2stop;

        r_fixed4(r_tri_mesh &xin, edge_bdry &bin, int dstr, int dstp, int d2str, int d2stp)
            : r_fixed(xin,bin), d2start(d2str), d2stop(d2stp) {mytype="fixed4";} 
        r_fixed4(const r_fixed4 &inbdry, r_tri_mesh &xin, edge_bdry &bin) : r_fixed(inbdry,xin,bin), d2start(inbdry.d2start), d2stop(inbdry.d2stop) {mytype="fixed4";}
        r_fixed4* create(r_tri_mesh &xin, edge_bdry &bin) const {return new r_fixed4(*this,xin,bin);}
        
        
        void fixdx2() {
            for(int j=0;j<base.nseg;++j) {
                int sind = base.seg(j);
                    for(int n=d2start;n<=d2stop;++n) {
                        x.gbl->res1(x.seg(sind).pnt(0))(n) = 0.0;
                        x.gbl->res1(x.seg(sind).pnt(1))(n) = 0.0;
                    }
            }
            return;
        }
};

class r_translating : public r_fixed {
    private:
        FLT dx[2];
    public:
        r_translating(r_tri_mesh &xin, edge_bdry &bin): r_fixed(xin,bin) { 
            for(int n=0;n<tri_mesh::ND;++n)
                dx[n] = 0.0;
            mytype = "translating";
        }
        r_translating(const r_translating &inbdry, r_tri_mesh &xin, edge_bdry &bin) : r_fixed(inbdry,xin,bin){ 
            for(int n=0;n<tri_mesh::ND;++n)
                dx[n] = inbdry.dx[n];
            mytype = "translating";
        }
        r_translating* create(r_tri_mesh &xin, edge_bdry &bin) const {return new r_translating(*this,xin,bin);}

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
            int n,p0;
            for(int j=1;j<base.nseg;++j) {
                p0 = x.seg(base.seg(j)).pnt(0);
                for(n=0;n<2;++n)
                    x.pnts(p0)(n) += dx[n];
            }
            
            
            /* TEMPORARY */
                p0 = x.seg(base.seg(0)).pnt(0);
                for(n=0;n<2;++n)
                    x.pnts(p0)(n) += dx[n]/2;
                    
                p0 = x.seg(base.seg(base.nseg-1)).pnt(1);
                for(n=0;n<2;++n)
                    x.pnts(p0)(n) += dx[n]/2;
                
        }
};

class r_oscillating : public r_fixed {
    public:
        FLT p0, amp, omega;
        
        r_oscillating(r_tri_mesh &xin, edge_bdry &bin) : r_fixed(xin,bin), p0(0.0), amp(0.0), omega(0.0) {mytype="oscillating";}
        r_oscillating(const r_oscillating &inbdry, r_tri_mesh &xin, edge_bdry &bin) : r_fixed(inbdry,xin,bin), p0(inbdry.p0), amp(inbdry.amp), omega(inbdry.omega) {mytype="oscillating";}
        r_oscillating* create(r_tri_mesh &xin, edge_bdry &bin) const {return new r_oscillating(*this,xin,bin);}            
        
            
        void input(input_map& inmap) {
            r_fixed::input(inmap);
            
            std::string line;
            if(inmap.getline(base.idprefix+"_r_oscillate",line)) {
                std::istringstream data(line);
                data >> p0 >> amp >> omega;
                data.clear();
            }
            else {
                p0 = 0.0;
                amp = 0.0;
                omega = 0.0;
            }
        }  
        
        void output(std::ostream& fout) {
            r_fixed::output(fout);
            fout << base.idprefix << "_r_oscillate: " << p0 << '\t' << amp << '\t' << omega << std::endl;
        }
            
        void tadvance() {
            int n,vrt;
            FLT center[2], center1[2], dx[2],xp[2];
            FLT theta, theta1,dtheta;
            
            center[0] = p0*x.gbl->time;
            center[1] = amp*(1-cos(omega*x.gbl->time));
            theta = atan(omega*amp*sin(omega*x.gbl->time)/p0);            
            
            center1[0] = p0*x.gbl->time;
            center1[1] = amp*(1-cos(omega*x.gbl->time));
            theta1 = atan(omega*amp*sin(omega*x.gbl->time)/p0);
            dtheta = theta1-theta;
                        
            for(int j=0;j<base.nseg;++j) {
                vrt = x.seg(base.seg(j)).pnt(0);
                xp[0] = x.pnts(vrt)(0)-center[0];
                xp[1] = x.pnts(vrt)(1)-center[1];            
                dx[0] = center1[0]-center[0] +xp[0]*cos(dtheta)-xp[1]*sin(dtheta) -xp[0];
                dx[1] = center1[1]-center[1] +xp[0]*sin(dtheta)+xp[1]*cos(dtheta) -xp[1];
                for(n=0;n<2;++n)
                    x.pnts(vrt)(n) += dx[n];
            }
        }
};
