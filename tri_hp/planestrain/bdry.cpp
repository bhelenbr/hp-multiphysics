#include "bdry_ps.h"
#include <myblas.h>

/*************************************************/
/* SET DIRICHLET BOUNDARY VALUES & FLUXES ********/
/* (THINGS THAT ARE INDEPENDENT OF THE SOLUTION) */
/* BUT NOT INDEPENDENT OF TIME/MESH POSITION     */
/*************************************************/
//void chrctr(FLT rho, FLT gam, double wl[NV], double wr[NV], double norm[ND], double mv[ND]);


using namespace bdry_ps;

void friction_wall::rsdl(int stage) {
	int j,k,m,n,seg,v0,v1,sind,tind;
	TinyVector<FLT,2> pt,nrm;
	TinyVector<FLT,3> u,flx;
	TinyVector<FLT,tri_mesh::ND> stress;
	FLT visc[x.ND][x.ND][x.ND][x.ND];

	for(j=0;j<base.nseg;++j) {
		sind = base.seg(j);
		v0 = x.seg(sind).pnt(0);
		v1 = x.seg(sind).pnt(1);
		tind = x.seg(sind).tri(0);        

		for(seg=0;seg<3;++seg)
			if (x.tri(tind).seg(seg) == sind) break;

		x.crdtocht(tind);
		for(m=basis::tri(x.log2p)->bm();m<basis::tri(x.log2p)->tm();++m)
			for(n=0;n<x.ND;++n)
				x.cht(n,m) = 0.0;

		for(n=0;n<x.ND;++n)
			basis::tri(x.log2p)->proj_side(seg,&x.cht(n,0), &x.crd(n)(0,0), &x.dcrd(n,0)(0,0), &x.dcrd(n,1)(0,0));

		x.ugtouht(tind);

		for(n=0;n<x.NV;++n)
			basis::tri(x.log2p)->proj_side(seg,&x.uht(n)(0),&x.u(n)(0,0),&x.du(n,0)(0,0),&x.du(n,1)(0,0));

		for (k=0;k<basis::tri(x.log2p)->gpx();++k) {
			x.cjcb(0,k) = x.hp_ps_gbl->mu/(x.dcrd(0,0)(0,k)*x.dcrd(1,1)(0,k) -x.dcrd(1,0)(0,k)*x.dcrd(0,1)(0,k));

			/* BIG FAT UGLY VISCOUS TENSOR (LOTS OF SYMMETRY THOUGH)*/
			/* INDICES ARE 1: EQUATION U OR V, 2: VARIABLE (U OR V), 3: EQ. DERIVATIVE (R OR S) 4: VAR DERIVATIVE (R OR S)*/
			visc[0][0][0][0] =  x.cjcb(0,k)*(2.*x.dcrd(1,1)(0,k)*x.dcrd(1,1)(0,k) +x.dcrd(0,1)(0,k)*x.dcrd(0,1)(0,k));
			visc[0][0][1][1] =  x.cjcb(0,k)*(2.*x.dcrd(1,0)(0,k)*x.dcrd(1,0)(0,k) +x.dcrd(0,0)(0,k)*x.dcrd(0,0)(0,k));
			visc[0][0][0][1] = -x.cjcb(0,k)*(2.*x.dcrd(1,1)(0,k)*x.dcrd(1,0)(0,k) +x.dcrd(0,1)(0,k)*x.dcrd(0,0)(0,k));
#define      viscI0II0II1II0I visc[0][0][0][1]

			visc[1][1][0][0] =  x.cjcb(0,k)*(x.dcrd(1,1)(0,k)*x.dcrd(1,1)(0,k) +2.*x.dcrd(0,1)(0,k)*x.dcrd(0,1)(0,k));
			visc[1][1][1][1] =  x.cjcb(0,k)*(x.dcrd(1,0)(0,k)*x.dcrd(1,0)(0,k) +2.*x.dcrd(0,0)(0,k)*x.dcrd(0,0)(0,k));
			visc[1][1][0][1] = -x.cjcb(0,k)*(x.dcrd(1,1)(0,k)*x.dcrd(1,0)(0,k) +2.*x.dcrd(0,1)(0,k)*x.dcrd(0,0)(0,k));
#define      viscI1II1II1II0I visc[1][1][0][1]

			visc[0][1][0][0] = -x.cjcb(0,k)*x.dcrd(0,1)(0,k)*x.dcrd(1,1)(0,k);
			visc[0][1][1][1] = -x.cjcb(0,k)*x.dcrd(0,0)(0,k)*x.dcrd(1,0)(0,k);
			visc[0][1][0][1] =  x.cjcb(0,k)*x.dcrd(0,1)(0,k)*x.dcrd(1,0)(0,k);
			visc[0][1][1][0] =  x.cjcb(0,k)*x.dcrd(0,0)(0,k)*x.dcrd(1,1)(0,k);

			/* OTHER SYMMETRIES     */                
#define      viscI1II0II0II0I visc[0][1][0][0]
#define      viscI1II0II1II1I visc[0][1][1][1]
#define      viscI1II0II0II1I visc[0][1][1][0]
#define      viscI1II0II1II0I visc[0][1][0][1]

			stress(0) =    (-x.u(2)(0,k)*x.dcrd(1,0)(0,k) 
							-viscI0II0II1II0I*x.du(0,0)(0,k) -visc[0][1][1][0]*x.du(1,0)(0,k)
							-visc[0][0][1][1]*x.du(0,1)(0,k) -visc[0][1][1][1]*x.du(1,1)(0,k));															
			stress(1) =    ( x.u(2)(0,k)*x.dcrd(0,0)(0,k)
							-viscI1II0II1II0I*x.du(0,0)(0,k) -viscI1II1II1II0I*x.du(1,0)(0,k)
							-viscI1II0II1II1I*x.du(0,1)(0,k) -visc[1][1][1][1]*x.du(1,1)(0,k));

			nrm(0) = x.dcrd(1,0)(0,k);
			nrm(1) = -x.dcrd(0,0)(0,k);                
			for(n=0;n<tri_mesh::ND;++n) {
				pt(n) = x.crd(n)(0,k);
			}

			for(n=0;n<x.NV;++n)
				u(n) = x.u(n)(0,k);

			flux(u,stress,pt,nrm,flx);

			for(n=0;n<x.NV;++n)
				x.res(n)(0,k) = flx(n);

		}

		for(n=0;n<x.NV;++n)
			basis::tri(x.log2p)->intgrt1d(&x.lf(n)(0),&x.res(n)(0,0));

		for(n=0;n<x.NV;++n)
			x.hp_gbl->res.v(v0,n) += x.lf(n)(0);

		for(n=0;n<x.NV;++n)
			x.hp_gbl->res.v(v1,n) += x.lf(n)(1);

		for(k=0;k<basis::tri(x.log2p)->sm();++k) {
			for(n=0;n<x.NV;++n)
				x.hp_gbl->res.s(sind,k,n) += x.lf(n)(k+2);
		}
	}

	return;
}

void curve_edges::setvalues(init_bdry_cndtn *ibc, const std::vector<int>& indices) {
    int j,v0,v1,sind,info;
    TinyVector<FLT,tri_mesh::ND> pt;
    char uplo[] = "U";
    
    /* SET DISPLACEMENT OF VERTICES TO ZERO */
    j = 0;
    do {
        sind = base.seg(j);
        v0 = x.seg(sind).pnt(0);
        for(int n=0;n<tri_mesh::ND;++n)
            x.ug.v(v0,n) = 0.0;
    } while(++j < base.nseg);
    v0 = x.seg(sind).pnt(1);
    for(int n=0;n<tri_mesh::ND;++n)
        x.ug.v(v0,n) = 0.0;
    
    if (basis::tri(x.log2p)->p() == 1) return;

    /*******************/
    /* SET SIDE DISPLACEMENTS */
    /*******************/
    for(j=0;j<base.nseg;++j) {
        sind = base.seg(j);
        v0 = x.seg(sind).pnt(0);
        v1 = x.seg(sind).pnt(1);
        
        base.mvpttobdry(j,-1.0,pt);
        for (int n=0;n<tri_mesh::ND;++n)
            x.ug.v(v0,n) = pt(n) -x.pnts(v0)(n);
            
        base.mvpttobdry(j,1.0,pt);
        for (int n=0;n<tri_mesh::ND;++n)
            x.ug.v(v1,n) = pt(n) -x.pnts(v1)(n);
            
        for(int n=0;n<tri_mesh::ND;++n) {
            basis::tri(x.log2p)->proj1d(x.ug.v(v0,n),x.ug.v(v1,n),&x.res(n)(0,0));
            basis::tri(x.log2p)->proj1d(x.pnts(v0)(n),x.pnts(v1)(n),&x.crd(n)(0,0));
            
            for(int k=0;k<basis::tri(x.log2p)->gpx();++k)
                x.dcrd(n,0)(0,k) = 0.5*(x.pnts(v1)(n)-x.pnts(v0)(n));
        }
        
        for(int k=0;k<basis::tri(x.log2p)->gpx(); ++k) {
            pt(0) = x.crd(0)(0,k)+x.res(0)(0,k);
            pt(1) = x.crd(1)(0,k)+x.res(1)(0,k);
            base.mvpttobdry(j,basis::tri(x.log2p)->xp(k),pt);
            x.crd(0)(0,k) -= pt(0);
            x.crd(1)(0,k) -= pt(1);
        }
        
        for(int n=0;n<tri_mesh::ND;++n) {
            basis::tri(x.log2p)->intgrt1d(&x.cf(n,0),&x.crd(n)(0,0));
#ifdef F2CFortran
            DPBTRS(uplo,basis::tri(x.log2p)->sm(),basis::tri(x.log2p)->sbwth(),1,(double *) &basis::tri(x.log2p)->sdiag1d(0,0),basis::tri(x.log2p)->sbwth()+1,&x.cf(n,2),basis::tri(x.log2p)->sm(),info);
#else
            const int sbwth = basis::tri(x.log2p)->sbwth(), one = 1, sm = basis::tri(x.log2p)->sm();
            const int sbp1 = sbwth +1;
            dpbtrs_(uplo,&sm,&sbwth,&one,(double *) &basis::tri(x.log2p)->sdiag1d(0,0),&sbp1,&x.cf(n,2),&sm,&info);
#endif
            for(int m=0;m<basis::tri(x.log2p)->sm();++m)
                x.ug.s(sind,m,n) = -x.cf(n,m+2);
        }
    }
    
    return;
}
