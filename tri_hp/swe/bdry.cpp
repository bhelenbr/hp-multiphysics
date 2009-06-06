#include "bdry_swe.h"

/*************************************************/
/* SET DIRICHLET BOUNDARY VALUES & FLUXES ********/
/* (THINGS THAT ARE INDEPENDENT OF THE SOLUTION) */
/* BUT NOT INDEPENDENT OF TIME/MESH POSITION     */
/*************************************************/

using namespace bdry_swe;

void characteristic::flux(Array<FLT,1>& u, TinyVector<FLT,tri_mesh::ND> xpt, TinyVector<FLT,tri_mesh::ND> mv, TinyVector<FLT,tri_mesh::ND> norm, Array<FLT,1>& flx) {
    FLT ul,vl,hl,hr,hul,hur,hvl,hvr,hu,hv,h,w1,w2,w3;
    FLT um,vm,c,c2,lam0,lam1,lam2,mag,dxmax;
    FLT pre,qmax,fmax,alpha,alpha2,sigma;
    Array<FLT,1> ub(x.NV), uvp(x.NV);

    /* CHARACTERISTIC FAR-FIELD B.C. */ 
    mag = sqrt(norm(0)*norm(0) + norm(1)*norm(1));
    norm(0) /= mag;
    norm(1) /= mag;
    um =  mv(0)*norm(0) +mv(1)*norm(1);
    vm = -mv(0)*norm(1) +mv(1)*norm(0);

    hl = u(x.NV-1);    
    hul =  u(0)*norm(0) +u(1)*norm(1);
    hvl = -u(0)*norm(1) +u(1)*norm(0);

    /* FREESTREAM CONDITIONS */
    for(int n=0;n<x.NV;++n)
		ub(n) = ibc->f(n,xpt,x.gbl->time);

    hr = ub(x.NV-1);    
    hur =  ub(0)*norm(0) +ub(1)*norm(1);
    hvr = -ub(0)*norm(1) +ub(1)*norm(0);


    /* CALCULATE EIGENVALUES */
    ul = hul/hl -um;
    vl = hvl/hl -vm;
    qmax = ul*ul +vl*vl;
    fmax = fabs(x.gbl->f0 +x.gbl->cbeta*xpt(1));
    c2 = x.gbl->g*hl;
    dxmax = 2.*mag/(0.25*(basis::tri(x.log2p).p +1)*(basis::tri(x.log2p).p+1));
    alpha = x.gbl->cd*dxmax/(2*hl);
    alpha2 = alpha*alpha;
    sigma = MAX((qmax -c2)/qmax,0);
//    pre = x.gbl->ptest*(pow(dxmax*(gbl->bd(0)+fmax),2.0) +(3.+alpha2)*qmax)/(dxmax*dxmax*(gbl->bd(0)*gbl->bd(0) +sigma*fmax*fmax) +c2 +(3.+sigma*alpha2)*qmax);
    pre = x.gbl->ptest*(pow(dxmax*(x.gbl->bd(0)*.5 +fmax),2.0) +(1.+alpha2)*qmax)/(dxmax*dxmax*(x.gbl->bd(0)*x.gbl->bd(0)*.25 +sigma*fmax*fmax) +c2 +(1.+sigma*alpha2)*qmax);

    lam0 = ul;
    c = sqrt(ul*ul*(1-pre) +pre*c2);
    lam1 = ul+c;
    lam2 = ul-c;

    if (lam0 < 0.0) {        
		if (lam1 < 0.0) {
			/* SUPERCRITICAL CASE */
			hu = hur;
			hv = hvr;
			h = hr;
		}
		else {
			/* SUBCRITICAL */
			w1 = hvr;
			w2 = 1./(2.*c)*(pre*hul -lam2*hl);
			w3 = 1./(2.*c)*(-pre*hur +lam1*hr);

			hv = w1;
			hu = 1/pre*(lam1*w2 +lam2*w3);
			h = w2 +w3;
		}
    }
    else {
		for(int n=tri_mesh::ND;n<x.NV-1;++n)
			ub(n) = u(n);

		if (lam2 > 0.0) {
			/* SUPERCRITICAL CASE */
			hu = hul;
			hv = hvl;
			h = hl;
		}
		else {
			/* SUBCRITICAL */
			w1 = hvl;
			w2 = 1./(2.*c)*(pre*hul -lam2*hl);
			w3 = 1./(2.*c)*(-pre*hur +lam1*hr);

			hv = w1;
			hu = 1/pre*(lam1*w2 +lam2*w3);
			h = w2 +w3;


		}
    }

    /* CHANGE BACK TO X,Y COORDINATES */
    ub(0) =  hu*norm(0) -hv*norm(1);
    ub(1) =  hu*norm(1) +hv*norm(0);
    ub(x.NV-1) =  h;

    norm *= mag;

    flx(x.NV-1) = ((ub(0) -ub(x.NV-1)*mv(0))*norm(0) +(ub(1) -ub(x.NV-1)*mv(1))*norm(1));

    for(int n=0;n<tri_mesh::ND;++n)
		flx(n) = flx(x.NV-1)*ub(n)/ub(x.NV-1) +ub(x.NV-1)*ub(x.NV-1)*x.gbl->g*norm(n)/2.;

    for(int n=tri_mesh::ND;n<x.NV-1;++n)
		flx(n) = flx(x.NV-1)*ub(n);

    return;
}
