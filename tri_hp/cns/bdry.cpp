#include "bdry_cns.h"
#include <myblas.h>

/*************************************************/
/* SET DIRICHLET BOUNDARY VALUES & FLUXES ********/
/* (THINGS THAT ARE INDEPENDENT OF THE SOLUTION) */
/* BUT NOT INDEPENDENT OF TIME/MESH POSITION     */
/*************************************************/
using namespace bdry_cns;

void generic::output(const std::string& filename, tri_hp::filetype typ,int tlvl) {
	int i,m,n,ind,sind,tind,seg;
	const int NV = 4;
	TinyMatrix<TinyMatrix<FLT,tri_mesh::ND,tri_mesh::ND>,NV-1,NV-1> visc;
	TinyMatrix<FLT,tri_mesh::ND,tri_mesh::ND> kcond;
	TinyVector<FLT,tri_mesh::ND> norm, mvel;
	FLT lkcond = x.gbl->kcond;
	FLT lmu = x.gbl->mu;
	FLT convect,jcb;
	FLT gam = x.gbl->gamma;
	FLT gm1 = gam-1.0;
	FLT ogm1 = 1.0/gm1;
	FLT gogm1 = gam*ogm1;
    
    FLT Pout, uout, RTout;
	
	std::string fname;
	fname = filename +"_" +base.idprefix;

	switch(typ) {
		case(tri_hp::tecplot): {
			if (!report_flag) return;
			
			Array<FLT,1> lconv_flux(x.NV),ldiff_flux(x.NV);
			fname += ".dat";
			
			std::ofstream file_out;
			file_out.open(fname);
			
			file_out << "VARIABLES=\"S\",\"X\",\"Y\",\"CFLUX0\",\"DLUX0\",\"CFLUX1\",\"DLUX1\",\"CFLUX2\",\"DLUX2\",\"DLUX3\",\"CFLUX3\",\"DLUX3\"\nTITLE = " << base.idprefix << '\n'<< "ZONE\n";

			conv_flux = 0.0;
			diff_flux = 0.0;
			moment = 0.0;
			circumference = 0.0;
            
            Pout = 0.0;
            uout = 0.0;
            RTout = 0.0;
#ifdef L2_ERROR
			FLT l2error = 0.0;
			TinyVector<FLT,2> xpt;
#endif

			for(ind=0; ind < base.nseg; ++ind) {
				sind = base.seg(ind);
				tind = x.seg(sind).tri(0);        

				for(seg=0;seg<3;++seg)
					if (x.tri(tind).seg(seg) == sind) break;
				assert(seg != 3);

				x.crdtocht(tind);
				for(m=basis::tri(x.log2p)->bm();m<basis::tri(x.log2p)->tm();++m)
					for(n=0;n<tri_mesh::ND;++n)
						x.cht(n,m) = 0.0;

				for(n=0;n<tri_mesh::ND;++n)
					basis::tri(x.log2p)->proj_side(seg,&x.cht(n,0), &x.crd(n)(0,0), &x.dcrd(n,0)(0,0), &x.dcrd(n,1)(0,0));

				x.ugtouht(tind);
				for(n=0;n<x.NV;++n)
					basis::tri(x.log2p)->proj_side(seg,&x.uht(n)(0),&x.u(n)(0,0),&x.du(n,0)(0,0),&x.du(n,1)(0,0));

				for (i=0;i<basis::tri(x.log2p)->gpx();++i) {
					jcb =  basis::tri(x.log2p)->wtx(i)*RAD(x.crd(0)(0,i))*sqrt(x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i) +x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i));
					circumference += jcb;

					double mujcbi = lmu*RAD(x.crd(0)(0,i))/(x.dcrd(0,0)(0,i)*x.dcrd(1,1)(0,i) -x.dcrd(1,0)(0,i)*x.dcrd(0,1)(0,i));
					double kcjcbi = lkcond*RAD(x.crd(0)(0,i))/x.gbl->R/(x.dcrd(0,0)(0,i)*x.dcrd(1,1)(0,i) -x.dcrd(1,0)(0,i)*x.dcrd(0,1)(0,i));

					/* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
					/* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
					visc(0,0)(0,0) = -mujcbi*(4./3.*x.dcrd(1,1)(0,i)*x.dcrd(1,1)(0,i) +x.dcrd(0,1)(0,i)*x.dcrd(0,1)(0,i));
					visc(0,0)(1,1) = -mujcbi*(4./3.*x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i));
					visc(0,0)(0,1) =  mujcbi*(4./3.*x.dcrd(1,1)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,1)(0,i)*x.dcrd(0,0)(0,i));
#define                 viscI0II0II1II0I visc(0,0)(0,1)

					visc(1,1)(0,0) = -mujcbi*(x.dcrd(1,1)(0,i)*x.dcrd(1,1)(0,i) +4./3.*x.dcrd(0,1)(0,i)*x.dcrd(0,1)(0,i));
					visc(1,1)(1,1) = -mujcbi*(x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i) +4./3.*x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i));
					visc(1,1)(0,1) =  mujcbi*(x.dcrd(1,1)(0,i)*x.dcrd(1,0)(0,i) +4./3.*x.dcrd(0,1)(0,i)*x.dcrd(0,0)(0,i));
#define                 viscI1II1II1II0I visc(1,1)(0,1)

					visc(0,1)(0,0) =  mujcbi*1./3.*x.dcrd(0,1)(0,i)*x.dcrd(1,1)(0,i);
					visc(0,1)(1,1) =  mujcbi*1./3.*x.dcrd(0,0)(0,i)*x.dcrd(1,0)(0,i);
					visc(0,1)(0,1) = -mujcbi*(x.dcrd(0,1)(0,i)*x.dcrd(1,0)(0,i)-2./3.*x.dcrd(0,0)(0,i)*x.dcrd(1,1)(0,i));
					visc(0,1)(1,0) = -mujcbi*(x.dcrd(0,0)(0,i)*x.dcrd(1,1)(0,i)-2./3.*x.dcrd(0,1)(0,i)*x.dcrd(1,0)(0,i));

					/* OTHER SYMMETRIES     */
#define				viscI1II0II0II0I visc(0,1)(0,0)
#define				viscI1II0II1II1I visc(0,1)(1,1)
#define				viscI1II0II0II1I visc(0,1)(1,0)
#define				viscI1II0II1II0I visc(0,1)(0,1)

					/* HEAT DIFFUSION TENSOR */
					/* DIFFUSIVE FLUXES ( FOR EXTRA VARIABLES) */
					kcond(1,1) = -kcjcbi*(x.dcrd(1,0)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,0)(0,i)*x.dcrd(0,0)(0,i));
					kcond(0,1) =  kcjcbi*(x.dcrd(1,1)(0,i)*x.dcrd(1,0)(0,i) +x.dcrd(0,1)(0,i)*x.dcrd(0,0)(0,i));
#define				kcondI1II0I kcond(0,1)

					ldiff_flux(0) = 0.0;
					ldiff_flux(1) =    basis::tri(x.log2p)->wtx(i)*(-x.u(0)(0,i)*RAD(x.crd(0)(0,i))*x.dcrd(1,0)(0,i)
									-viscI0II0II1II0I*x.du(1,0)(0,i) -visc(0,1)(1,0)*x.du(2,0)(0,i)
									-visc(0,0)(1,1)*x.du(1,1)(0,i) -visc(0,1)(1,1)*x.du(2,1)(0,i));
					ldiff_flux(2) =    basis::tri(x.log2p)->wtx(i)*( x.u(0)(0,i)*RAD(x.crd(0)(0,i))*x.dcrd(0,0)(0,i)
									-viscI1II0II1II0I*x.du(1,0)(0,i) -viscI1II1II1II0I*x.du(2,0)(0,i)
									-viscI1II0II1II1I*x.du(1,1)(0,i) -visc(1,1)(1,1)*x.du(2,1)(0,i));
					ldiff_flux(3) = basis::tri(x.log2p)->wtx(i)*(-kcond(0,1)*x.du(3,0)(0,i) -kcond(1,1)*x.du(3,1)(0,i));
					diff_flux -= ldiff_flux;
					ldiff_flux /= jcb;

					norm(0) = x.dcrd(1,0)(0,i);
					norm(1) = -x.dcrd(0,0)(0,i);
					for(n=0;n<tri_mesh::ND;++n) {
						mvel(n) = x.gbl->bd(0)*(x.crd(n)(0,i) -dxdt(x.log2p,ind)(n,i));
#ifdef MESH_REF_VEL
						mvel(n) += x.gbl->mesh_ref_vel(n);
#endif
					}
					double rho = x.u(0)(0,i)/x.u(NV-1)(0,i);
					double h = gogm1*x.u(NV-1)(0,i) +0.5*(x.u(1)(0,i)*x.u(1)(0,i)+x.u(2)(0,i)*x.u(2)(0,i));
					convect = basis::tri(x.log2p)->wtx(i)*RAD(x.crd(0)(0,i))*rho*((x.u(1)(0,i)-mvel(0))*norm(0) +(x.u(2)(0,i)-mvel(1))*norm(1));
					lconv_flux(0) = convect;
					lconv_flux(1) = x.u(1)(0,i)*convect;
					lconv_flux(2) = x.u(2)(0,i)*convect;
					lconv_flux(3) = h*convect;
					conv_flux -= lconv_flux;
					lconv_flux /= jcb;

#ifdef L2_ERROR
					xpt(0) = x.crd(0)(0,i);
					xpt(1) = x.crd(1)(0,i);
					l2error += jcb*l2norm->Eval(xpt,x.gbl->time);
#endif

					file_out << circumference << ' ' << x.crd(0)(0,i) << ' ' << x.crd(1)(0,i) << ' ';
					for(int n=0;n<x.NV;++n)
						file_out << x.u(n)(0,i) << ' ' << lconv_flux(n) << ' ' << ldiff_flux(n) << ' ';
					file_out << std::endl;
                    
                    
//                    Pout += basis::tri(x.log2p)->wtx(i)*x.u(0)(0,i)*sqrt(norm(0)*norm(0)+norm(1)*norm(1));
//                    uout += basis::tri(x.log2p)->wtx(i)*x.u(1)(0,i)*sqrt(norm(0)*norm(0)+norm(1)*norm(1));
//                    RTout += basis::tri(x.log2p)->wtx(i)*x.u(3)(0,i)*sqrt(norm(0)*norm(0)+norm(1)*norm(1));
                    
                    Pout += basis::tri(x.log2p)->wtx(i)*(x.u(0)(0,i)*norm(0)+x.u(0)(0,i)*norm(1));
                    uout += basis::tri(x.log2p)->wtx(i)*(x.u(1)(0,i)*norm(0)+x.u(2)(0,i)*norm(1));
                    RTout += basis::tri(x.log2p)->wtx(i)*(x.u(3)(0,i)*norm(0)+x.u(3)(0,i)*norm(1));
                    
				}

			}
			file_out.close();


			*x.gbl->log << setprecision(14) << base.idprefix << " circumference: " << circumference << std::endl;
			*x.gbl->log << base.idprefix << " viscous/pressure flux: " << diff_flux << std::endl;
			*x.gbl->log << base.idprefix << " convective flux: " << conv_flux << std::endl;
			*x.gbl->log << base.idprefix << " circulation: " << circulation << std::endl;

			/* OUTPUT AUXILIARY FLUXES */
			*x.gbl->log << base.idprefix << "total fluxes: " << total_flux << std::endl;
#ifdef L2_ERROR
			*x.gbl->log << base.idprefix << "l2error: " << sqrt(l2error) << std::endl;
#endif
            
            *x.gbl->log << base.idprefix << " Integrated Pressure: " << Pout << std::endl;
            *x.gbl->log << base.idprefix << " Integrated vel: " << uout << std::endl;
            *x.gbl->log << base.idprefix << " Integrated RT: " << RTout << std::endl;
            
			break;
		}
			
		default: {
			hp_edge_bdry::output(filename,typ,tlvl);
			break;
		}
	}

	return;
}

void generic::flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx) {
	
	//			TinyVector<FLT,4> ub,fluxtemp;
	//			FLT mag = sqrt(norm(0)*norm(0) + norm(1)*norm(1));
	//			norm /= mag;
	//			
	//			/* rotate coordinate system */
	//			FLT ul =  (u(1)-mv(0))*norm(0) +(u(2)-mv(1))*norm(1);
	//			FLT vl = -(u(1)-mv(0))*norm(1) +(u(2)-mv(1))*norm(0);
	//						
	//			/* far field */
	//			FLT pr = ibc->f(0,xpt,x.gbl->time);
	//			FLT RT = ibc->f(x.NV-1,xpt,x.gbl->time);
	//			//RT = u(x.NV-1);
	//
	//			FLT rho = pr*RT;
	//			
	//			/* CONTINUITY */
	//			fluxtemp(0) = rho*ul;
	//			fluxtemp(1) = fluxtemp(0)*ul+pr;
	//			fluxtemp(2) = fluxtemp(0)*vl;
	//			
	//			FLT h = x.gbl->gamma/(x.gbl->gamma-1.0)*RT +0.5*(ul*ul+vl*vl);
	//			fluxtemp(3) = fluxtemp(0)*h;
	//			
	//			/* CHANGE BACK TO X,Y COORDINATES */
	//			flx(0) = fluxtemp(0);
	//			flx(1) = fluxtemp(1)*norm(0) - fluxtemp(2)*norm(1);
	//			flx(2) = fluxtemp(1)*norm(1) + fluxtemp(2)*norm(0);
	//			flx(3) = fluxtemp(3);
	//			
	//			flx *= mag;
	
	flx(0) = ibc->f(0, xpt, x.gbl->time)/u(x.NV-1)*((u(1) -mv(0))*norm(0) +(u(2) -mv(1))*norm(1));
	//		flx(0) = ibc->f(0, xpt, x.gbl->time)/ibc->f(x.NV-1, xpt, x.gbl->time)*((u(1) -mv(0))*norm(0) +(u(2) -mv(1))*norm(1));
	
	/* X&Y MOMENTUM */
#ifdef INERTIALESS
	for (int n=1;n<tri_mesh::ND+1;++n)
		flx(n) = ibc->f(0, xpt, x.gbl->time)*norm(n-1);
#else
	for (int n=1;n<tri_mesh::ND+1;++n)
		flx(n) = flx(0)*u(n) +ibc->f(0, xpt, x.gbl->time)*norm(n-1);
#endif
	
	/* ENERGY EQUATION */
    //double E = u(x.NV-1)/(x.gbl->gamma-1.0)+0.5*(u(1)*u(1)+u(2)*u(2));
    double rho = u(0)/u(x.NV-1);
    double E = ibc->f(0, xpt, x.gbl->time)/(rho*(x.gbl->gamma-1.0))+0.5*(u(1)*u(1)+u(2)*u(2));
    flx(x.NV-1) = rho*E*((u(1)-mv(0))*norm(0)+(u(2)-mv(1))*norm(1))+ibc->f(0, xpt, x.gbl->time)*(u(1)*norm(0)+u(2)*norm(1));
	
	return;
}

void inflow::modify_boundary_residual() {
	int j,k,m,n,v0,v1,sind,info;
	TinyVector<FLT,tri_mesh::ND> pt;
	TinyVector<double,MXGP> res1d;
	TinyVector<double,MXTM> rescoef;
	char uplo[] = "U";
	
	FLT ogm1 = 1.0/(x.gbl->gamma-1.0);
	
	j = 0;
	do {
		sind = base.seg(j);	
		v0 = x.seg(sind).pnt(0);
		
		FLT KE = 0.5*(ibc->f(1, x.pnts(v0), x.gbl->time)*ibc->f(1, x.pnts(v0), x.gbl->time)+ibc->f(2, x.pnts(v0), x.gbl->time)*ibc->f(2, x.pnts(v0), x.gbl->time));
		
		x.gbl->res.v(v0,1) = x.gbl->res.v(v0,0)*ibc->f(1, x.pnts(v0), x.gbl->time);
		x.gbl->res.v(v0,2) = x.gbl->res.v(v0,0)*ibc->f(2, x.pnts(v0), x.gbl->time);
		x.gbl->res.v(v0,3) = x.gbl->res.v(v0,0)*(ibc->f(3, x.pnts(v0), x.gbl->time)*ogm1+KE);
		
	} while (++j < base.nseg);
	
	v0 = x.seg(sind).pnt(1);
	
	FLT KE = 0.5*(ibc->f(1, x.pnts(v0), x.gbl->time)*ibc->f(1, x.pnts(v0), x.gbl->time)+ibc->f(2, x.pnts(v0), x.gbl->time)*ibc->f(2, x.pnts(v0), x.gbl->time));
	
	x.gbl->res.v(v0,1) = x.gbl->res.v(v0,0)*ibc->f(1, x.pnts(v0), x.gbl->time);
	x.gbl->res.v(v0,2) = x.gbl->res.v(v0,0)*ibc->f(2, x.pnts(v0), x.gbl->time);
	x.gbl->res.v(v0,3) = x.gbl->res.v(v0,0)*(ibc->f(3, x.pnts(v0), x.gbl->time)*ogm1+KE);
	
	if(basis::tri(x.log2p)->sm()) {
		for(j=0;j<base.nseg;++j) {
			sind = base.seg(j);
			v0 = x.seg(sind).pnt(0);
			v1 = x.seg(sind).pnt(1);
			
			if (is_curved()) {
				x.crdtocht1d(sind);
				for(n=0;n<tri_mesh::ND;++n)
					basis::tri(x.log2p)->proj1d(&x.cht(n,0),&x.crd(n)(0,0),&x.dcrd(n,0)(0,0));
			}
			else {
				for(n=0;n<tri_mesh::ND;++n) {
					basis::tri(x.log2p)->proj1d(x.pnts(v0)(n),x.pnts(v1)(n),&x.crd(n)(0,0));
					
					for(k=0;k<basis::tri(x.log2p)->gpx();++k)
						x.dcrd(n,0)(0,k) = 0.5*(x.pnts(v1)(n)-x.pnts(v0)(n));
				}
			}
			
			/* take global coefficients and put into local vector */
			/*  only need res_rho */
			for (m=0; m<2; ++m) 
				rescoef(m) = x.gbl->res.v(x.seg(sind).pnt(m),0);					
			
			for (m=0;m<basis::tri(x.log2p)->sm();++m) 
				rescoef(m+2) = x.gbl->res.s(sind,m,0);					
			
			basis::tri(x.log2p)->proj1d(&rescoef(0),&res1d(0));
			
			for(n=1;n<x.NV;++n)
				basis::tri(x.log2p)->proj1d(x.gbl->res.v(v0,n),x.gbl->res.v(v1,n),&x.res(n)(0,0));
			
			for(k=0;k<basis::tri(x.log2p)->gpx(); ++k) {
				pt(0) = x.crd(0)(0,k);
				pt(1) = x.crd(1)(0,k);
				
				FLT KE = 0.5*(ibc->f(1, pt, x.gbl->time)*ibc->f(1, pt, x.gbl->time)+ibc->f(2, pt, x.gbl->time)*ibc->f(2, pt, x.gbl->time));
				x.res(1)(0,k) -= res1d(k)*ibc->f(1, pt, x.gbl->time);
				x.res(2)(0,k) -= res1d(k)*ibc->f(2, pt, x.gbl->time);
				x.res(3)(0,k) -= res1d(k)*(ibc->f(3, pt, x.gbl->time)*ogm1+KE);
				
			}
			
			for(n=1;n<x.NV;++n){
				basis::tri(x.log2p)->intgrt1d(&x.lf(n)(0),&x.res(n)(0,0));
#ifdef F2CFortran
				PBTRS(uplo,basis::tri(x.log2p)->sm(),basis::tri(x.log2p)->sbwth(),1,(double *) &basis::tri(x.log2p)->sdiag1d(0,0),basis::tri(x.log2p)->sbwth()+1,&x.lf(n)(2),basis::tri(x.log2p)->sm(),info);
#else
                const int sbwth = basis::tri(x.log2p)->sbwth(), one = 1, sm = basis::tri(x.log2p)->sm();
                const int sbp1 = sbwth +1;
                dpbtrs_(uplo,&sm,&sbwth,&one,(double *) &basis::tri(x.log2p)->sdiag1d(0,0),&sbp1,&x.lf(n)(2),&sm,&info);
#endif
                for(m=0;m<basis::tri(x.log2p)->sm();++m)
					x.gbl->res.s(sind,m,n) = -x.lf(n)(2+m);						
				
			}								
		}
	}
	
	
	
	return;
}

void adiabatic::flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx) {
	
	/* CONTINUITY */
	flx(0) = ibc->f(0, xpt, x.gbl->time)/ibc->f(x.NV-1, xpt, x.gbl->time)*((u(1) -mv(0))*norm(0) +(u(2) -mv(1))*norm(1));
	
	/* MOMENTUM */
	for (int n=1;n<x.NV-1;++n)
		flx(n) = 0.0;
	
	/* ENERGY EQUATION */
	double h = x.gbl->gamma/(x.gbl->gamma-1.0)*ibc->f(x.NV-1, xpt, x.gbl->time) +0.5*(u(1)*u(1)+u(2)*u(2));
	flx(x.NV-1) = h*flx(0);
	
	return;
}


void adiabatic::modify_boundary_residual() {
	int j,k,m,n,v0,v1,sind,info;
	TinyVector<FLT,tri_mesh::ND> pt;
	TinyVector<double,MXGP> res1d;
	TinyVector<double,MXTM> rescoef;
	char uplo[] = "U";
	
//	FLT ogm1 = 1.0/(x.gbl->gamma-1.0);
	
	j = 0;
	do {
		sind = base.seg(j);	
		v0 = x.seg(sind).pnt(0);		

		x.gbl->res.v(v0,1) = x.gbl->res.v(v0,0)*ibc->f(1, x.pnts(v0), x.gbl->time);
		x.gbl->res.v(v0,2) = x.gbl->res.v(v0,0)*ibc->f(2, x.pnts(v0), x.gbl->time);
		
	} while (++j < base.nseg);
	
	v0 = x.seg(sind).pnt(1);
	
	x.gbl->res.v(v0,1) = x.gbl->res.v(v0,0)*ibc->f(1, x.pnts(v0), x.gbl->time);
	x.gbl->res.v(v0,2) = x.gbl->res.v(v0,0)*ibc->f(2, x.pnts(v0), x.gbl->time);
	
	if(basis::tri(x.log2p)->sm()) {
		for(j=0;j<base.nseg;++j) {
			sind = base.seg(j);
			v0 = x.seg(sind).pnt(0);
			v1 = x.seg(sind).pnt(1);
			
			if (is_curved()) {
				x.crdtocht1d(sind);
				for(n=0;n<tri_mesh::ND;++n)
					basis::tri(x.log2p)->proj1d(&x.cht(n,0),&x.crd(n)(0,0),&x.dcrd(n,0)(0,0));
			}
			else {
				for(n=0;n<tri_mesh::ND;++n) {
					basis::tri(x.log2p)->proj1d(x.pnts(v0)(n),x.pnts(v1)(n),&x.crd(n)(0,0));
					
					for(k=0;k<basis::tri(x.log2p)->gpx();++k)
						x.dcrd(n,0)(0,k) = 0.5*(x.pnts(v1)(n)-x.pnts(v0)(n));
				}
			}
			
			/* take global coefficients and put into local vector */
			/*  only need res_rho */
			for (m=0; m<2; ++m) 
				rescoef(m) = x.gbl->res.v(x.seg(sind).pnt(m),0);					
			
			for (m=0;m<basis::tri(x.log2p)->sm();++m) 
				rescoef(m+2) = x.gbl->res.s(sind,m,0);					
			
			basis::tri(x.log2p)->proj1d(&rescoef(0),&res1d(0));
			
			for(n=1;n<x.NV;++n)
				basis::tri(x.log2p)->proj1d(x.gbl->res.v(v0,n),x.gbl->res.v(v1,n),&x.res(n)(0,0));
			
			for(k=0;k<basis::tri(x.log2p)->gpx(); ++k) {
				pt(0) = x.crd(0)(0,k);
				pt(1) = x.crd(1)(0,k);
				
				x.res(1)(0,k) -= res1d(k)*ibc->f(1, pt, x.gbl->time);
				x.res(2)(0,k) -= res1d(k)*ibc->f(2, pt, x.gbl->time);
				
			}
			
			for(n=1;n<x.NV-1;++n){
				basis::tri(x.log2p)->intgrt1d(&x.lf(n)(0),&x.res(n)(0,0));
#ifdef F2CFortran
				PBTRS(uplo,basis::tri(x.log2p)->sm(),basis::tri(x.log2p)->sbwth(),1,(double *) &basis::tri(x.log2p)->sdiag1d(0,0),basis::tri(x.log2p)->sbwth()+1,&x.lf(n)(2),basis::tri(x.log2p)->sm(),info);
#else
                const int sbwth = basis::tri(x.log2p)->sbwth(), one = 1, sm = basis::tri(x.log2p)->sm();
                const int sbp1 = sbwth +1;
                dpbtrs_(uplo,&sm,&sbwth,&one,(double *) &basis::tri(x.log2p)->sdiag1d(0,0),&sbp1,&x.lf(n)(2),&sm,&info);
#endif
				for(m=0;m<basis::tri(x.log2p)->sm();++m) 
					x.gbl->res.s(sind,m,n) = -x.lf(n)(2+m);						
				
			}								
		}
	}
	
	
	
	return;
}

void characteristic::flux(Array<FLT,1>& pvu, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx) {
	
	TinyVector<FLT,4> lambda,Rl,Rr,ub,Roe,fluxtemp,fluxleft, fluxright;
	Array<FLT,2> A(x.NV,x.NV),V(x.NV,x.NV),VINV(x.NV,x.NV),temp(x.NV,x.NV),P(x.NV,x.NV),Pinv(x.NV,x.NV),dpdc(x.NV,x.NV), dcdp(x.NV,x.NV);
	FLT gam = x.gbl->gamma;
	FLT gm1 = gam-1.0;
	FLT gogm1 = gam/gm1;
    
    Array<FLT,2> Aeigs(x.NV,x.NV);
//    Array<FLT,1> Aeigs(x.NV);
	
	/* Left */
	/* Rotate Coordinate System */
	FLT ul =  pvu(1)*norm(0) +pvu(2)*norm(1);
	FLT vl = -pvu(1)*norm(1) +pvu(2)*norm(0);
    pvu(1) = ul; pvu(2) = vl;
    
    FLT mv_n = mv(0)*norm(0)+mv(1)*norm(1);
	
	/* Roe Variables */
	Rl(0) = sqrt(pvu(0)/pvu(x.NV-1)); // sqrt(rho)
	Rl(1) = Rl(0)*ul; // sqrt(rho)*u
	Rl(2) = Rl(0)*vl; // sqrt(rho)*v
	Rl(3) = Rl(0)*(pvu(x.NV-1)/gm1+0.5*(ul*ul+vl*vl)); // sqrt(rho)*E
	
	/* Right */
	for(int n=0;n<x.NV;++n)
		ub(n) = ibc->f(n,xpt,x.gbl->time);
	
	/* Rotate Coordinate System */
	FLT ur =  ub(1)*norm(0) +ub(2)*norm(1);
	FLT vr = -ub(1)*norm(1) +ub(2)*norm(0);
    ub(1) = ur; ub(2) = vr;
	
	/* Roe Variables */
	Rr(0) = sqrt(ub(0)/ub(x.NV-1));
	Rr(1) = Rr(0)*ur;
	Rr(2) = Rr(0)*vr;
	Rr(3) = Rr(0)*(ub(x.NV-1)/gm1+0.5*(ur*ur+vr*vr));	
	
	/* Average Roe Variables */
	Roe = 0.5*(Rl+Rr);
	
	/* Calculate u,v,c Variables */
	FLT rho = Roe(0)*Roe(0);
	FLT u = Roe(1)/Roe(0);
	FLT v = Roe(2)/Roe(0);
	FLT ke = 0.5*(u*u+v*v);
	FLT E = Roe(3)/Roe(0);
	FLT rt = gm1*(E-ke);
	FLT pr = rho*rt;
	FLT c2 = gam*rt;
	FLT c = sqrt(c2);

	fluxleft(0) = (pvu(0)/pvu(x.NV-1))*(pvu(1)-mv_n);
	fluxleft(1) = fluxleft(0)*pvu(1)+pvu(0);
	fluxleft(2) = fluxleft(0)*pvu(2);
//	fluxleft(3) = fluxleft(0)*(gogm1*pvu(x.NV-1)+0.5*(pvu(1)*pvu(1)+pvu(2)*pvu(2)));
    FLT E_pvu = pvu(x.NV-1)/gm1+0.5*(pvu(1)*pvu(1)+pvu(2)*pvu(2));
    fluxleft(3) = fluxleft(0)*E_pvu+pvu(0)*pvu(1);

	fluxright(0) = (ub(0)/ub(x.NV-1))*(ub(1)-mv_n);
	fluxright(1) = fluxright(0)*ub(1)+ub(0);
	fluxright(2) = fluxright(0)*ub(2);
//	fluxright(3) = fluxright(0)*(gogm1*ub(x.NV-1)+0.5*(ub(1)*ub(1)+ub(2)*ub(2)));
    FLT E_ub = ub(x.NV-1)/gm1+0.5*(ub(1)*ub(1)+ub(2)*ub(2));
    fluxright(3) = fluxright(0)*E_ub+ub(0)*ub(1);
	
	fluxtemp = 0.5*(fluxleft+fluxright);
	
	FLT nu = x.gbl->mu/rho;
//	FLT cp = gogm1*x.gbl->R;
//	FLT alpha = x.gbl->kcond/(rho*cp);
	FLT h = side_length; 
	
//	FLT hdt = 0.5*h*x.gbl->bd(0)/c;
    //To steaady state
    FLT hdt = 0.0;
    
	FLT umag = sqrt(u*u+v*v);
	FLT M = MAX(umag/c,1.0e-5);
	FLT nuh = 4.0*nu/(h*c);
//    FLT nuh = 3.0*MAX(4.0*nu/(3.0*h*c),alpha/(h*c));
//	FLT alh = 2.0*alpha/(h*c);//maybe it should be smaller?
//
//	FLT b2 = MIN(M*M/(1.0-M*M) + hdt*hdt + nuh*nuh + alh*alh,1.0);
    
    FLT b2;
    if (M<0.8){
        b2 = MIN(M*M/(1.0-M*M) + hdt*hdt + nuh*nuh,1.0);
//        b2 = MIN(3.0*M*M + (hdt+nuh)*(hdt+nuh), 1.0);
    }
    else{
        b2 = 1.0;
    }
    
	FLT alph = 0.0;
//    FLT alph = 1.0+b2;
    
    //Turn preconditioner off
    b2 = 1.0;
	
	/* Inverse of Preconditioner */
	Pinv = 1.0/b2,					 0.0, 0.0, 0.0,
		   alph*u/(pr*gam*b2),		 1.0, 0.0, 0.0,
		   alph*v/(pr*gam*b2),		 0.0, 1.0, 0.0,
	       -(b2-1.0)/(gogm1*rho*b2), 0.0, 0.0, 1.0;
	
	/* jacobian of conservative wrt primitive */
	dcdp = 1.0/rt,               0.0,   0.0,   -rho/rt,
		   u/rt,                 rho,   0.0,   -rho*u/rt,
		   v/rt,                 0.0,   rho,   -rho*v/rt,
		   (rt+gm1*ke)/(gm1*rt), rho*u, rho*v, -rho*ke/rt;
    
    temp = 0.0;
    for(int i=0; i<x.NV; ++i)
        for(int j=0; j<x.NV; ++j)
            for(int k=0; k<x.NV; ++k)
                temp(i,j)+=dcdp(i,k)*Pinv(k,j);
    
    Pinv = temp;
	
	FLT temp1 = sqrt(u*u*(1.0-2.0*b2+b2*b2)+4.0*b2*c2);
    
    V = 0.5*(u*(b2-1.0)+temp1)*rho,         0.5*(u*(b2-1.0)-temp1)*rho,         0.0, 0.0,
        1.0,                                1.0,                                0.0, 0.0,
        0.0,                                0.0,                                1.0, 0.0,
        0.5*(u*(b2-1.0)*gm1+gm1*temp1)/gam, 0.5*(u*(b2-1.0)*gm1-gm1*temp1)/gam, 0.0, 1.0;
    
    VINV =  1.0/(temp1*rho), -0.5*(u*(b2-1.0)-temp1)/temp1, 0.0, 0.0,
            -1.0/(temp1*rho), 0.5*(u*(b2-1.0)+temp1)/temp1, 0.0, 0.0,
            0.0,              0.0,                            1.0, 0.0,
            -gm1/(gam*rho),   0.0,                            0.0, 1.0;
    
    Aeigs = 0.5*(u+u*b2+temp1)-mv_n, 0.0, 0.0, 0.0,
            0.0, 0.5*(u+u*b2-temp1)-mv_n, 0.0, 0.0,
            0.0, 0.0, u-mv_n, 0.0,
            0.0, 0.0, 0.0, u-mv_n;

    for(int i=0; i<x.NV; ++i)
        Aeigs(i,i) = fabs(Aeigs(i,i));

    temp = 0.0;
    for(int i=0; i<x.NV; ++i)
        for(int j=0; j<x.NV; ++j)
            for(int k=0; k<x.NV; ++k)
                temp(i,j) += Aeigs(i,k)*VINV(k,j);

    A = 0.0;
    for(int i=0; i<x.NV; ++i)
        for(int j=0; j<x.NV; ++j)
            for(int k=0; k<x.NV; ++k)
                A(i,j)+=V(i,k)*temp(k,j);


    temp = 0.0;
    for(int i=0; i<x.NV; ++i)
        for(int j=0; j<x.NV; ++j)
            for(int k=0; k<x.NV; ++k)
                temp(i,j)+=Pinv(i,k)*A(k,j);

    A = temp;

//    temp = 0.0;
//    for(int i=0; i<x.NV; ++i)
//        for(int j=0; j<x.NV; ++j)
//            for(int k=0; k<x.NV; ++k)
//                temp(i,j)+=dcdp(i,k)*A(k,j);
//    A = temp;
    
    
	for(int i = 0; i < x.NV; ++i)
		for(int j = 0; j < x.NV; ++j)
			fluxtemp(i) -= 0.5*A(i,j)*(ub(j)-pvu(j));
	
	/* CHANGE BACK TO X,Y COORDINATES */
	flx(0) = fluxtemp(0);
	flx(1) = fluxtemp(1)*norm(0) - fluxtemp(2)*norm(1);
	flx(2) = fluxtemp(1)*norm(1) + fluxtemp(2)*norm(0);
	flx(3) = fluxtemp(3);
	
	return;
}

void applied_stress::init(input_map& inmap,void* gbl_in) {
	std::string keyword;
	std::ostringstream nstr;
	
	generic::init(inmap,gbl_in);
	
	stress.resize(x.NV-1);
	
	for(int n=0;n<x.NV-1;++n) {
		nstr.str("");
		nstr << base.idprefix << "_stress" << n << std::flush;
		if (inmap.find(nstr.str()) != inmap.end()) {
			stress(n).init(inmap,nstr.str());
		}
		else {
			*x.gbl->log << "couldn't find stress function " << nstr.str() << std::endl;
			sim::abort(__LINE__,__FILE__,x.gbl->log);
		}
	}
	
	return;
}

void applied_stress::flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx) {
	
	/* CONTINUITY */
	flx(0) = ibc->f(0, xpt, x.gbl->time)/u(x.NV-1)*((u(1) -mv(0))*norm(0) +(u(2) -mv(1))*norm(1));
	
	/* X&Y MOMENTUM */
#ifdef INERTIALESS
	for (int n=0;n<tri_mesh::ND;++n)
		flx(n+1) = -stress(n).Eval(xpt,x.gbl->time) +ibc->f(0, xpt, x.gbl->time)*norm(n);
#else
	for (int n=0;n<tri_mesh::ND;++n)
		flx(n+1) = flx(0)*u(n+1) -stress(n).Eval(xpt,x.gbl->time) +ibc->f(0, xpt, x.gbl->time)*norm(n);
#endif
	
	/* ENERGY EQUATION */
	double h = x.gbl->gamma/(x.gbl->gamma-1.0)*u(x.NV-1) +0.5*(u(1)*u(1)+u(2)*u(2));
	flx(x.NV-1) = h*flx(0)-stress(2).Eval(xpt,x.gbl->time);
	
	return;
}

void symmetry::tadvance() {
	int j,m,v0,sind;
	TinyVector<FLT,tri_mesh::ND> pt;
	
	hp_edge_bdry::tadvance();
	
	/* UPDATE BOUNDARY CONDITION VALUES */
	j = 0;
	do {
		sind = base.seg(j);
		v0 = x.seg(sind).pnt(0);
		x.ug.v(v0,dir+1) = 0.0;
	}	while (++j < base.nseg);
	v0 = x.seg(sind).pnt(1);
	x.ug.v(v0,dir+1) = 0.0;
	
	/*******************/
	/* SET SIDE VALUES */
	/*******************/
	for(j=0;j<base.nseg;++j) {
		sind = base.seg(j);
		for(m=0;m<basis::tri(x.log2p)->sm();++m) {
			x.ug.s(sind,m,dir+1) = 0.0;
		}
	}
	
	return;
}


void outflow_supersonic::flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx) {
    
    flx(0) = u(0)/u(x.NV-1)*((u(1) -mv(0))*norm(0) +(u(2) -mv(1))*norm(1));
    
    /* X&Y MOMENTUM */
#ifdef INERTIALESS
    for (int n=1;n<tri_mesh::ND+1;++n)
        flx(n) = ibc->f(0, xpt, x.gbl->time)*norm(n-1);
#else
    for (int n=1;n<tri_mesh::ND+1;++n)
        flx(n) = flx(0)*u(n) +u(0)*norm(n-1);
#endif
    
    /* ENERGY EQUATION */
    double rho = u(0)/u(x.NV-1);
    double E = u(0)/(rho*(x.gbl->gamma-1.0))+0.5*(u(1)*u(1)+u(2)*u(2));
    flx(x.NV-1) = rho*E*((u(1)-mv(0))*norm(0)+(u(2)-mv(1))*norm(1))+u(0)*(u(1)*norm(0)+u(2)*norm(1));
    
    return;
}



void euler::flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx) {
    
    Array<FLT,1> ub(x.NV);
    for(int n=0;n<x.NV;++n)
        ub(n) = ibc->f(n,xpt,x.gbl->time);
    
    flx(0) = u(0)/u(x.NV-1)*((ub(1) -mv(0))*norm(0) +(ub(2) -mv(1))*norm(1));
    
    /* X&Y MOMENTUM */
#ifdef INERTIALESS
    for (int n=1;n<tri_mesh::ND+1;++n)
        flx(n) = ibc->f(0, xpt, x.gbl->time)*norm(n-1);
#else
    for (int n=1;n<tri_mesh::ND+1;++n)
        flx(n) = flx(0)*ub(n) +u(0)*norm(n-1);
#endif
    
    /* ENERGY EQUATION */
    double rho = u(0)/u(x.NV-1);
    double E = u(0)/(rho*(x.gbl->gamma-1.0))+0.5*(ub(1)*ub(1)+ub(2)*ub(2));
    flx(x.NV-1) = rho*E*((ub(1)-mv(0))*norm(0)+(ub(2)-mv(1))*norm(1))+u(0)*(ub(1)*norm(0)+ub(2)*norm(1));
    
    return;
}


void wall::flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, FLT side_length, Array<FLT,1>& flx) {
    
    //Zero nomral velocity
//    FLT u_x, u_y, u_tan, norm_x, norm_y, tan_x, tan_y;
//    norm_x = norm(0);
//    norm_y = norm(1);
//
//    //Tangent vector points clockwise
//    tan_x = -norm(1);
//    tan_y = norm(0);
//    u_tan = u(1)*tan_x+u(2)*tan_y;
//
//    //u_norm is zero
//    if(tan_x != 0){
//        u_y = -u_tan*norm_x/(norm_y*tan_x-tan_y*norm_x);
//        u_x = u_tan/tan_x-u_y*(tan_y/tan_x);
//    }
//    else{
//        u_y = u_tan;
//        u_x = 0.0;
//    }
//
//
//    flx(0) = u(0)/u(x.NV-1)*((u_x -mv(0))*norm(0) +(u_y -mv(1))*norm(1));
//
//    /* X&Y MOMENTUM */
//    flx(1) = flx(0)*u_x +u(0)*norm(0);
//    flx(2) = flx(0)*u_y +u(0)*norm(1);
//
//    /* ENERGY EQUATION */
//    double rho = u(0)/u(x.NV-1);
//    double E = u(x.NV-1)/(x.gbl->gamma-1.0)+0.5*(u_x*u_x+u_y*u_y);
//    flx(x.NV-1) = rho*E*((u_x-mv(0))*norm(0)+(u_y-mv(1))*norm(1))+u(0)*(u_x*norm(0)+u_y*norm(1));
    
    
    flx(0) = 0.0;

    /* X&Y MOMENTUM */
    flx(1) = u(0)*norm(0);
    flx(2) = u(0)*norm(1);

    /* ENERGY EQUATION */
    flx(x.NV-1) = 0.0;
    
    
    return;
}


