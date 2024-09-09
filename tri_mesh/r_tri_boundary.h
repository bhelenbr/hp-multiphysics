/*
 *  r_tri_boundary.h
 *  mesh
 *
 *  Created by Brian Helenbrook on Mon Jun 10 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_mesh.h"

#ifdef MY_SPARSE
#include <myblas.h>
#else
#include <petscksp.h>
#endif

class r_side_bdry {
	protected:
		std::string mytype;
		r_tri_mesh &x;
		edge_bdry &base;

	public:
		r_side_bdry(r_tri_mesh &xin, edge_bdry &bin) : x(xin), base(bin) {mytype="plain";}
		r_side_bdry(const r_side_bdry &inbdry, r_tri_mesh &xin, edge_bdry &bin) : x(xin), base(bin) {mytype="plain";}
		virtual r_side_bdry* create(r_tri_mesh &xin, edge_bdry &bin) const {return new r_side_bdry(*this,xin,bin);}
		virtual void init(input_map& bdrydata) {}
		virtual void tadvance() {}
		virtual void dirichlet() {}
		virtual void fixdx2() {}
		virtual ~r_side_bdry() {}
#ifdef MY_SPARSE
		virtual void jacobian(sparse_row_major &J, sparse_row_major &J_mpi, int stride) {}
		virtual void jacobian_dirichlet(sparse_row_major &J, sparse_row_major &J_mpi, int stride) {}
#else
		virtual void jacobian(Mat petsc_J, int stride) {}
		virtual void jacobian_dirichlet(Mat petsc_J, int stride) {}
#endif
};

class r_fixed : public r_side_bdry {
	public:
		int dstart, dstop;

		r_fixed(r_tri_mesh &xin, edge_bdry &bin) : r_side_bdry(xin,bin), dstart(0), dstop(1) {mytype = "fixed";}
		r_fixed(const r_fixed &inbdry, r_tri_mesh &xin, edge_bdry &bin) : r_side_bdry(inbdry,xin,bin), dstart(inbdry.dstart), dstop(inbdry.dstop) {mytype = "fixed";}
		r_fixed* create(r_tri_mesh &xin, edge_bdry &bin) const {return new r_fixed(*this,xin,bin);}

		void init(input_map& inmap) {
			r_side_bdry::init(inmap);

			std::string line;
			inmap.getlinewdefault(base.idprefix+"_r_dir",line,"0 1");
			std::istringstream data(line);
			data >> dstart >> dstop;
			data.clear();
		}

		void dirichlet() {
			for(int j=0;j<base.nseg;++j) {
				int sind = base.seg(j);
				for(int n=dstart;n<=dstop;++n) {
					x.r_gbl->res(x.seg(sind).pnt(0))(n) = 0.0;
					x.r_gbl->res(x.seg(sind).pnt(1))(n) = 0.0;
				}
			}
			return;
		}
		
#ifdef MY_SPARSE
		void jacobian_dirichlet(sparse_row_major &J, sparse_row_major &J_mpi, int stride) {
//		void jacobian_dirichlet(Mat petsc_J, int stride) {
			int np = (base.nseg+1)*(dstop -dstart +1);
			Array<int,1> points(np);
			
			int j = 0;
			int cnt = 0;
			int sind;
			do {
				sind = base.seg(j);
				int gindx = x.seg(sind).pnt(0)*stride +stride -tri_mesh::ND +dstart;
				for (int n=dstart;n<=dstop;++n)
					points(cnt++) = gindx++;
			} while(++j < base.nseg);
			int gindx = x.seg(sind).pnt(1)*stride +stride -tri_mesh::ND +dstart;
			for (int n=dstart;n<=dstop;++n)
				points(cnt++) = gindx++;
			
			J.zero_rows(cnt,points);
			J_mpi.zero_rows(cnt,points);
			J.set_diag(cnt,points,1.0);
//			MatZeroRows(petsc_J,cnt,points.data(),1.0);
		}
#endif
};

class r_fixed_angled : public r_side_bdry {
	public:
		FLT theta;

		r_fixed_angled(r_tri_mesh &xin, edge_bdry &bin) : r_side_bdry(xin,bin), theta(0.0) {mytype = "fixed_angled";}
		r_fixed_angled(const r_fixed_angled &inbdry, r_tri_mesh &xin, edge_bdry &bin) : r_side_bdry(inbdry,xin,bin), theta(inbdry.theta) {mytype = "fixed_angled";}
		r_fixed_angled* create(r_tri_mesh &xin, edge_bdry &bin) const {return new r_fixed_angled(*this,xin,bin);}

		void init(input_map& inmap) {
			r_side_bdry::init(inmap);

			inmap.get(base.idprefix+"_r_theta",theta);
			theta *= M_PI/180.0;
		}

		void dirichlet() {
			/* SKIP ENDPOINTS */
			for(int j=1;j<base.nseg;++j) {
				int sind = base.seg(j);
				/* Tangent Residual */
				FLT res = -x.r_gbl->res(x.seg(sind).pnt(0))(0)*sin(theta) +x.r_gbl->res(x.seg(sind).pnt(0))(1)*cos(theta);
				x.r_gbl->res(x.seg(sind).pnt(0))(0) = -res*sin(theta);
				x.r_gbl->res(x.seg(sind).pnt(0))(1) = res*cos(theta);
			}
			
			return;
		}
		
#ifdef MY_SPARSE 
		void jacobian_dirichlet(sparse_row_major &J, sparse_row_major &J_mpi, int stride) {
			/* SKIP ENDPOINTS */
			for(int j=0;j<base.nseg-1;++j) {
				int sind = base.seg(j);
				int row = x.seg(sind).pnt(1)*stride +stride-x.ND;

				int nnz1 = J._cpt(row+1) -J._cpt(row);
				int nnz2 = J._cpt(row+2) -J._cpt(row+1);
				Array<int,1> cols(nnz1);
				Array<FLT,2> vals(2,nnz1);
				
				/* SOME ERROR CHECKING TO MAKE SURE ROW SPARSENESS PATTERN IS THE SAME */
				if (nnz1 != nnz2) {
					*x.gbl->log << "zeros problem in deforming mesh on angled boundary\n";
					sim::abort(__LINE__,__FILE__,x.gbl->log);
				}
				int row1 = J._cpt(row);
				int row2 = J._cpt(row+1);
				for(int col=0;col<nnz1;++col) {
					if (J._col(row1++) != J._col(row2++)) {
						*x.gbl->log << "zeros indexing problem in deforming mesh on angled boundary\n";
						sim::abort(__LINE__,__FILE__,x.gbl->log);
					}	
				}

				vals(0,Range(0,nnz1-1)) = -J._val(Range(J._cpt(row),J._cpt(row+1)-1))*sin(theta);
				vals(0,Range(0,nnz1-1)) += J._val(Range(J._cpt(row+1),J._cpt(row+2)-1))*cos(theta);
								
				/* Replace x equation with tangential position equation */
				/* Replacy y equation with normal displacement equation */
				/* Normal Equation */
				vals(1,Range::all()) = 0.0;
				for(int col=0;col<nnz1;++col) {
					if (J._col(row1+col) == row) {
						vals(1,col) = cos(theta);
						break;
					}
				}
				for(int col=0;col<nnz1;++col) {
					if (J._col(row1+col) == row+1) {
						vals(1,col) = sin(theta);
						break;
					}
				}
			
				/* tangent = -sin(theta) i +cos(theta) j */
				/* normal = cos(theta) i + sin(theta) j */
				/* Rotate equations for diagonal dominance to match what is done to residual */
				Array<FLT,2> temp(2,nnz1);
				temp(0,Range::all()) = -vals(0,Range::all())*sin(theta) +vals(1,Range::all())*cos(theta);
				temp(1,Range::all()) =  vals(0,Range::all())*cos(theta) +vals(1,Range::all())*sin(theta);
			
				J._val(Range(J._cpt(row),J._cpt(row+1)-1)) = temp(0,Range::all());
				J._val(Range(J._cpt(row+1),J._cpt(row+2)-1)) = temp(1,Range::all());
			}
		}
#endif
};

class r_curved : public r_side_bdry {
	public:
		r_curved(r_tri_mesh &xin, edge_bdry &bin) : r_side_bdry(xin,bin) {mytype = "curved";}
		r_curved(const r_curved &inbdry, r_tri_mesh &xin, edge_bdry &bin) : r_side_bdry(inbdry,xin,bin) {mytype = "curved";}
		r_curved* create(r_tri_mesh &xin, edge_bdry &bin) const {return new r_curved(*this,xin,bin);}
		
		void dirichlet() {
			TinyVector<FLT,tri_mesh::ND> pt, norm;
			/* SKIP ENDPOINTS */
			for(int j=0;j<base.nseg-1;++j) {
				int sind = base.seg(j);
				int pnt = x.seg(sind).pnt(1);
				base.bdry_normal(j,1.0,norm);
				pt = x.pnts(pnt);
				base.mvpttobdry(j,1.0,pt);
				
				/* Tangent Residual */
				FLT rest = -x.r_gbl->res(pnt)(0)*norm(1) +x.r_gbl->res(pnt)(1)*norm(0);
				/* Normal Residual */
				FLT resn = ((x.pnts(pnt)(0)-pt(0))*norm(0) +(x.pnts(pnt)(1)-pt(1))*norm(1))/x.r_gbl->diag(pnt);

				x.r_gbl->res(pnt)(0) = -rest*norm(1) +resn*norm(0);
				x.r_gbl->res(pnt)(1) = rest*norm(0) +resn*norm(1);
			}
			
			return;
		}
	
#ifdef MY_SPARSE 
		void jacobian_dirichlet(sparse_row_major &J, sparse_row_major &J_mpi, int stride) {
			TinyVector<FLT,tri_mesh::ND> norm;			
			
			/* SKIP ENDPOINTS? */
			for(int j=0;j<base.nseg-1;++j) {
				base.bdry_normal(j,1.0,norm);
				int sind = base.seg(j);
				int pnt = x.seg(sind).pnt(1);
				int row = x.seg(sind).pnt(1)*stride +stride-x.ND;
				
				int nnz1 = J._cpt(row+1) -J._cpt(row);
				int nnz2 = J._cpt(row+2) -J._cpt(row+1);
				Array<int,1> cols(nnz1);
				Array<FLT,2> vals(2,nnz1);
				
				/* SOME ERROR CHECKING TO MAKE SURE ROW SPARSENESS PATTERN IS THE SAME */
				if (nnz1 != nnz2) {
					*x.gbl->log << "zeros problem in deforming mesh on angled boundary\n";
					sim::abort(__LINE__,__FILE__,x.gbl->log);
				}
				int row1 = J._cpt(row);
				int row2 = J._cpt(row+1);
				for(int col=0;col<nnz1;++col) {
					if (J._col(row1++) != J._col(row2++)) {
						*x.gbl->log << "zeros indexing problem in deforming mesh on angled boundary\n";
						sim::abort(__LINE__,__FILE__,x.gbl->log);
					}	
				}
				
				vals(0,Range(0,nnz1-1)) = -J._val(Range(J._cpt(row),J._cpt(row+1)-1))*norm(1);
				vals(0,Range(0,nnz1-1)) += J._val(Range(J._cpt(row+1),J._cpt(row+2)-1))*norm(0);
				
				/* Replace x equation with tangential position equation */
				/* Replacy y equation with normal displacement equation */
				/* Normal Equation */
				vals(1,Range::all()) = 0.0;
				for(int col=0;col<nnz1;++col) {
					if (J._col(row1+col) == row) {
						vals(1,col) = norm(0)/x.r_gbl->diag(pnt);
						break;
					}
				}
				for(int col=0;col<nnz1;++col) {
					if (J._col(row1+col) == row+1) {
						vals(1,col) = norm(1)/x.r_gbl->diag(pnt);
						break;
					}
				}
				
				/* tangent = -norm(1) i +norm(0) j */
				/* normal = norm(0) i + norm(1) j */
				/* Rotate equations for diagonal dominance to match what is done to residual */
				Array<FLT,2> temp(2,nnz1);
				temp(0,Range::all()) = -vals(0,Range::all())*norm(1) +vals(1,Range::all())*norm(0);
				temp(1,Range::all()) =  vals(0,Range::all())*norm(0) +vals(1,Range::all())*norm(1);
				
				J._val(Range(J._cpt(row),J._cpt(row+1)-1)) = temp(0,Range::all());
				J._val(Range(J._cpt(row+1),J._cpt(row+2)-1)) = temp(1,Range::all());
			}
		}
#endif
	
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
						x.r_gbl->res1(x.seg(sind).pnt(0))(n) = 0.0;
						x.r_gbl->res1(x.seg(sind).pnt(1))(n) = 0.0;
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

		void init(input_map& inmap) {
			r_fixed::init(inmap);

			FLT dx_dflt[2] = {0.0, 1.0};
			inmap.getwdefault(base.idprefix+"_r_translate",dx,2,dx_dflt);
		}

		void tadvance() {
			int n,p0;
			for(int j=1;j<base.nseg;++j) {
				p0 = x.seg(base.seg(j)).pnt(0);
				for(n=0;n<2;++n)
					x.pnts(p0)(n) += dx[n];
			}

			/* TEMPORARY I think this was to get my mpi test right */
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


		void init(input_map& inmap) {
			r_fixed::init(inmap);

			std::string line;
			inmap.getlinewdefault(base.idprefix+"_r_oscillate",line,std::string("0.0 0.0 0.0"));
			std::istringstream data(line);
			data >> p0 >> amp >> omega;
			data.clear();
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


class r_deforming : public r_fixed {
	bool do_left, do_right;
	public:
		r_deforming(r_tri_mesh &xin, edge_bdry &bin) : r_fixed(xin,bin), do_left(false), do_right(false) {mytype="deforming";}
		r_deforming(const r_deforming &inbdry, r_tri_mesh &xin, edge_bdry &bin) : r_fixed(inbdry,xin,bin), do_left(inbdry.do_left), do_right(inbdry.do_right) {mytype="deforming";}
		r_deforming* create(r_tri_mesh &xin, edge_bdry &bin) const {return new r_deforming(*this,xin,bin);}

		void init(input_map& inmap) {
			r_fixed::init(inmap);
			inmap.getwdefault(base.idprefix+"_r_do_left",do_left,false);
			inmap.getwdefault(base.idprefix+"_r_do_right",do_right,false);
		}
		
		void tadvance() {
			int vrt;
			for(int j=1-do_left;j<base.nseg;++j) {
				vrt = x.seg(base.seg(j)).pnt(0);
				base.mvpttobdry(j,-1.0,x.pnts(vrt));
			}
			if (do_right) {
				vrt = x.seg(base.seg(base.nseg-1)).pnt(1);
				base.mvpttobdry(base.nseg-1,1.0,x.pnts(vrt));
			}
		}
};

class r_vrtx_bdry {
	protected:
		std::string mytype;
		r_tri_mesh &x;
		vrtx_bdry &base;

	public:
		r_vrtx_bdry(r_tri_mesh &xin, vrtx_bdry &bin) : x(xin), base(bin) {mytype="plain";}
		r_vrtx_bdry(const r_vrtx_bdry &inbdry, r_tri_mesh &xin, vrtx_bdry &bin) : x(xin), base(bin) {mytype="plain";}
		virtual r_vrtx_bdry* create(r_tri_mesh &xin, vrtx_bdry &bin) const {return new r_vrtx_bdry(*this,xin,bin);}
		virtual void init(input_map& bdrydata) {}
		virtual void tadvance() {}
		virtual void dirichlet() {}
		virtual void fixdx2() {}
		virtual ~r_vrtx_bdry() {}
#ifdef MY_SPARSE
		virtual void jacobian(sparse_row_major &J, sparse_row_major &J_mpi, int stride) {}
		virtual void jacobian_dirichlet(sparse_row_major &J, sparse_row_major &J_mpi, int stride) {}
#else
		virtual void jacobian(Mat petsc_J, int stride) {}
		virtual void jacobian_dirichlet(Mat petsc_J, int stride) {}
#endif
	
};

class r_vfixed : public r_vrtx_bdry {
	public:
		int dstart, dstop;

		r_vfixed(r_tri_mesh &xin, vrtx_bdry &bin) : r_vrtx_bdry(xin,bin), dstart(0), dstop(1) {mytype = "vfixed";}
		r_vfixed(const r_vfixed &inbdry, r_tri_mesh &xin, vrtx_bdry &bin) : r_vrtx_bdry(inbdry,xin,bin), dstart(inbdry.dstart), dstop(inbdry.dstop) {mytype = "vfixed";}
		r_vfixed* create(r_tri_mesh &xin, vrtx_bdry &bin) const {return new r_vfixed(*this,xin,bin);}

		void init(input_map& inmap) {
			r_vrtx_bdry::init(inmap);

			std::string line;
			inmap.getlinewdefault(base.idprefix+"_r_dir",line,"0 1");
			std::istringstream data(line);
			data >> dstart >> dstop;
			data.clear();
		}

		void dirichlet() {
			int pnt = base.pnt;
			for(int n=dstart;n<=dstop;++n) {
				x.r_gbl->res(pnt)(n) = 0.0;
			}
			return;
		}
		
#ifdef MY_SPARSE
		void jacobian_dirichlet(sparse_row_major &J, sparse_row_major &J_mpi, int stride) {
			//		void jacobian_dirichlet(Mat petsc_J, int stride) {
			int np = (dstop -dstart +1);
			Array<int,1> points(np);
			
			int cnt = 0;
			int gindx = base.pnt*stride +stride -tri_mesh::ND +dstart;
			for (int n=dstart;n<=dstop;++n)
				points(cnt++) = gindx++;
			J.zero_rows(cnt,points);
			J_mpi.zero_rows(cnt,points);
			J.set_diag(cnt,points,1.0);
			//			MatZeroRows(petsc_J,cnt,points.data(),1.0);
		}
#endif
	

};

class r_vmoving : public r_vfixed {
	public:
		r_vmoving(r_tri_mesh &xin, vrtx_bdry &bin) : r_vfixed(xin,bin) {mytype="moving";}
		r_vmoving(const r_vmoving &inbdry, r_tri_mesh &xin, vrtx_bdry &bin) : r_vfixed(inbdry,xin,bin) {mytype="moving";}
		r_vmoving* create(r_tri_mesh &xin, vrtx_bdry &bin) const {return new r_vmoving(*this,xin,bin);}

		void tadvance() {
			int pnt = base.pnt;
			base.mvpttobdry(x.pnts(pnt));
		}
};



