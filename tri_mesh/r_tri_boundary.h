/*
 *  r_tri_boundary.h
 *  mesh
 *
 *  Created by Brian Helenbrook on Mon Jun 10 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */

#include "tri_mesh.h"

#ifndef petsc
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
		virtual void jacobian() {}
		virtual void jacobian_dirichlet() {}
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

class r_fixed_angled : public r_side_bdry {
	public:
		FLT theta;

		r_fixed_angled(r_tri_mesh &xin, edge_bdry &bin) : r_side_bdry(xin,bin), theta(0.0) {mytype = "fixed_angled";}
		r_fixed_angled(const r_fixed_angled &inbdry, r_tri_mesh &xin, edge_bdry &bin) : r_side_bdry(inbdry,xin,bin), theta(inbdry.theta) {mytype = "fixed_angled";}
		r_fixed_angled* create(r_tri_mesh &xin, edge_bdry &bin) const {return new r_fixed_angled(*this,xin,bin);}

		void input(input_map& inmap) {
			r_side_bdry::input(inmap);

			inmap.get(base.idprefix+"_r_theta",theta);
			theta *= M_PI/180.0;
		}

		void output(std::ostream& fout) {
			r_side_bdry::output(fout);
			fout << base.idprefix << "_r_theta: " << theta*180.0/M_PI << std::endl;
		}

		void dirichlet() {
			int sind;
			int j = 0;
			do {
				sind = base.seg(j);
				/* Tangent Residual */
				FLT res = -x.gbl->res(x.seg(sind).pnt(0))(0)*sin(theta) +x.gbl->res(x.seg(sind).pnt(0))(1)*cos(theta);
				x.gbl->res(x.seg(sind).pnt(0))(0) = -res*sin(theta);
				x.gbl->res(x.seg(sind).pnt(0))(1) = res*cos(theta);
			} while (++j < base.nseg);
			FLT res = -x.gbl->res(x.seg(sind).pnt(1))(0)*sin(theta) +x.gbl->res(x.seg(sind).pnt(1))(1)*cos(theta);
			x.gbl->res(x.seg(sind).pnt(1))(0) = -res*sin(theta);
			x.gbl->res(x.seg(sind).pnt(1))(1) = res*cos(theta);
			return;
		}
};
#else

#include <../tri_hp/tri_hp.h>

class r_side_bdry {
protected:
	std::string mytype;
	tri_hp& x;
	edge_bdry &base;
	
public:
	r_side_bdry(r_tri_mesh &xin, edge_bdry &bin) : x(dynamic_cast<tri_hp &>(xin)), base(bin) {mytype="plain";}
	r_side_bdry(const r_side_bdry &inbdry, r_tri_mesh &xin, edge_bdry &bin) : x(dynamic_cast<tri_hp &>(xin)), base(bin) {mytype="plain";}
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
	virtual void jacobian() {}
	virtual void jacobian_dirichlet() {}
};

class r_fixed : public r_side_bdry {
	protected:
		int dstart, dstop;
	
	public:
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
					x.r_tri_mesh::gbl->res(x.seg(sind).pnt(0))(n) = 0.0;
					x.r_tri_mesh::gbl->res(x.seg(sind).pnt(1))(n) = 0.0;
				}
			}
			return;
		}
		
		
		void jacobian_dirichlet() {
			int stride = x.NV +tri_mesh::ND;
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

			points += x.jacobian_start;
#ifdef MY_SPARSE
			x.J.zero_rows(cnt,points);
			x.J_mpi.zero_rows(cnt,points);
			x.J.set_diag(cnt,points,1.0);
#else
			MatZeroRows(x.petsc_J,cnt,points.data(),1.0);
#endif
		}
};

class r_fixed_angled : public r_fixed {
	public:
		FLT theta;

		r_fixed_angled(r_tri_mesh &xin, edge_bdry &bin) : r_fixed(xin,bin), theta(0.0) {mytype = "fixed_angled";}
		r_fixed_angled(const r_fixed_angled &inbdry, r_tri_mesh &xin, edge_bdry &bin) : r_fixed(inbdry,xin,bin), theta(inbdry.theta) {mytype = "fixed_angled";}
		r_fixed_angled* create(r_tri_mesh &xin, edge_bdry &bin) const {return new r_fixed_angled(*this,xin,bin);}

		void input(input_map& inmap) {
			r_fixed::input(inmap);

			inmap.get(base.idprefix+"_r_theta",theta);
			theta *= M_PI/180.0;
		}

		void output(std::ostream& fout) {
			r_fixed::output(fout);
			fout << base.idprefix << "_r_theta: " << theta*180.0/M_PI << std::endl;
		}

		void dirichlet() {
			int sind;
			
			/* SKIP ENDPOINTS */
			for(int j=0;j<base.nseg-1;++j) {
				sind = base.seg(j);
				/* Tangent Residual */
				/* tangent = -sin(theta) i +cos(theta) j */
				/* normal = cos(theta) i + sin(theta) j */
				FLT res = -x.r_tri_mesh::gbl->res(x.seg(sind).pnt(1))(0)*sin(theta) +x.r_tri_mesh::gbl->res(x.seg(sind).pnt(1))(1)*cos(theta);
				x.r_tri_mesh::gbl->res(x.seg(sind).pnt(1))(0) = -res*sin(theta);
				x.r_tri_mesh::gbl->res(x.seg(sind).pnt(1))(1) = res*cos(theta);
			} 
			
			return;
		}
		
		void jacobian() {}
		
		void jacobian_dirichlet() {
			int stride = x.NV +tri_mesh::ND;
			
#ifdef MY_SPARSE
			/* SKIP ENDPOINTS? */
			for(int j=0;j<base.nseg-1;++j) {
				int sind = base.seg(j);
				int row = x.seg(sind).pnt(1)*stride +x.NV+x.jacobian_start;

				int nnz1 = x.J._cpt(row+1) -x.J._cpt(row);
				int nnz2 = x.J._cpt(row+2) -x.J._cpt(row+1);
				Array<int,1> cols(nnz1);
				Array<FLT,2> vals(2,nnz1);
				
				/* SOME ERROR CHECKING TO MAKE SURE ROW SPARSENESS PATTERN IS THE SAME */
				if (nnz1 != nnz2) {
					*x.gbl->log << "zeros problem in deforming mesh on angled boundary\n";
					sim::abort(__LINE__,__FILE__,x.gbl->log);
				}
				int row1 = x.J._cpt(row);
				int row2 = x.J._cpt(row+1);
				for(int col=0;col<nnz1;++col) {
					if (x.J._col(row1++) != x.J._col(row2++)) {
						*x.gbl->log << "zeros indexing problem in deforming mesh on angled boundary\n";
						sim::abort(__LINE__,__FILE__,x.gbl->log);
					}	
				}

				vals(0,Range(0,nnz1-1)) = -x.J._val(Range(x.J._cpt(row),x.J._cpt(row+1)-1))*sin(theta);
				vals(0,Range(0,nnz1-1)) += x.J._val(Range(x.J._cpt(row+1),x.J._cpt(row+2)-1))*cos(theta);
								
				/* Replace x equation with tangential position equation */
				/* Replacy y equation with normal displacement equation */
				/* Normal Equation */
				vals(1,Range::all()) = 0.0;
				for(int col=0;col<nnz1;++col) {
					if (x.J._col(row1+col) == row) {
						vals(1,col) = cos(theta);
						break;
					}
				}
				for(int col=0;col<nnz1;++col) {
					if (x.J._col(row1+col) == row+1) {
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
			
				x.J._val(Range(x.J._cpt(row),x.J._cpt(row+1)-1)) = temp(0,Range::all());
				x.J._val(Range(x.J._cpt(row+1),x.J._cpt(row+2)-1)) = temp(1,Range::all());
			}
	#else
			const PetscInt *cols1, *cols2;
			const PetscScalar *vals1, *vals2;
			
			/* SKIP ENDPOINTS? */
			for(int j=0;j<base.nseg-1;++j) {
				int sind = base.seg(j);
				int row = x.seg(sind).pnt(1)*stride +x.NV +x.jacobian_start;
				int nnz1, nnz2;
				MatGetRow(x.petsc_J,row,&nnz1,&cols1,&vals1);
				MatGetRow(x.petsc_J,row+1,&nnz2,&cols2,&vals2);
				if (nnz1 != nnz2) {
					*x.gbl->log << "zeros problem in deforming mesh on angled boundary\n";
				}
				
				Array<int,1> cols(nnz1);
				Array<FLT,2> vals(2,nnz1);
				for (int col=0;col<nnz1;++col) {
					cols(col) = cols1[col];
					vals(0,col) = -vals1[col]*sin(theta);
				}
	
				for (int col=0;col<nnz2;++col) {
					if (cols(col) != cols2[col]) {
						*x.gbl->log << "zeros problem in deforming mesh on angled boundary\n";
						*x.gbl->log << cols << ' ' << col << ' ' << cols2[col] << std::endl;
					}
					// cols(col+nnz1) = cols2[col];
					vals(0,col) += vals2[col]*cos(theta);
				}
				MatRestoreRow(x.petsc_J,row,&nnz1,&cols1,&vals1);
				MatRestoreRow(x.petsc_J,row+1,&nnz2,&cols2,&vals2);
				
				/* Replace x equation with tangential position equation */
				/* Replacy y equation with normal displacement equation */
				/* Normal Equation */
				vals(1,Range::all()) = 0.0;
				for(int col=0;col<nnz1;++col) {
					if (cols(col) == row) {
						vals(1,col) = cos(theta);
						break;
					}
				}
				for(int col=0;col<nnz1;++col) {
					if (cols(col) == row+1) {
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
				
				TinyVector<int,2> rows(row,row+1);
				MatSetValues(x.petsc_J,2,rows.data(),nnz1,cols.data(),temp.data(),INSERT_VALUES);
				MatAssemblyBegin(x.petsc_J,MAT_FINAL_ASSEMBLY);
				MatAssemblyEnd(x.petsc_J,MAT_FINAL_ASSEMBLY);	
			}
	#endif
		}
	
};
#endif


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

			FLT dx_dflt[2] = {0.0, 1.0};
			inmap.getwdefault(base.idprefix+"_r_translate",dx,2,dx_dflt);
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


		void input(input_map& inmap) {
			r_fixed::input(inmap);

			std::string line;
			inmap.getlinewdefault(base.idprefix+"_r_oscillate",line,std::string("0.0 0.0 0.0"));
			std::istringstream data(line);
			data >> p0 >> amp >> omega;
			data.clear();
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


class r_deforming : public r_fixed {
	public:
		r_deforming(r_tri_mesh &xin, edge_bdry &bin) : r_fixed(xin,bin) {mytype="deforming";}
		r_deforming(const r_deforming &inbdry, r_tri_mesh &xin, edge_bdry &bin) : r_fixed(inbdry,xin,bin) {mytype="deforming";}
		r_deforming* create(r_tri_mesh &xin, edge_bdry &bin) const {return new r_deforming(*this,xin,bin);}

		void tadvance() {
			int vrt;
			for(int j=1;j<base.nseg;++j) {
				vrt = x.seg(base.seg(j)).pnt(0);
				base.mvpttobdry(j,-1.0,x.pnts(vrt));
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
		/* VIRTUAL FUNCTIONS FOR BOUNDARY DEFORMATION */
		virtual void output(std::ostream& fout) {
			fout << base.idprefix << "_r_vtype: " << mytype << std::endl;
		}
		virtual void input(input_map& bdrydata) {}
		virtual void tadvance() {}
		virtual void dirichlet() {}
		virtual void fixdx2() {}
		virtual ~r_vrtx_bdry() {}
};

class r_vfixed : public r_vrtx_bdry {
	public:
		int dstart, dstop;

		r_vfixed(r_tri_mesh &xin, vrtx_bdry &bin) : r_vrtx_bdry(xin,bin), dstart(0), dstop(1) {mytype = "vfixed";}
		r_vfixed(const r_vfixed &inbdry, r_tri_mesh &xin, vrtx_bdry &bin) : r_vrtx_bdry(inbdry,xin,bin), dstart(inbdry.dstart), dstop(inbdry.dstop) {mytype = "vfixed";}
		r_vfixed* create(r_tri_mesh &xin, vrtx_bdry &bin) const {return new r_vfixed(*this,xin,bin);}

		void input(input_map& inmap) {
			r_vrtx_bdry::input(inmap);

			std::string line;
			inmap.getlinewdefault(base.idprefix+"_r_dir",line,"0 1");
			std::istringstream data(line);
			data >> dstart >> dstop;
			data.clear();
		}

		void output(std::ostream& fout) {
			r_vrtx_bdry::output(fout);
			fout << base.idprefix << "_r_dir: " << dstart << '\t' << dstop << std::endl;
		}

		void dirichlet() {
			int pnt = base.pnt;
			for(int n=dstart;n<=dstop;++n) {
				x.gbl->res(pnt)(n) = 0.0;
			}
			return;
		}

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



