#include "bdry_swirl.h"

/*************************************************/
/* SET DIRICHLET BOUNDARY VALUES & FLUXES ********/
/* (THINGS THAT ARE INDEPENDENT OF THE SOLUTION) */
/* BUT NOT INDEPENDENT OF TIME/MESH POSITION     */
/*************************************************/

using namespace bdry_swirl;

void symmetry::tadvance() {
	int j,m,v0,sind,indx;
	TinyVector<FLT,tri_mesh::ND> pt;

	hp_edge_bdry::tadvance();

	/* UPDATE BOUNDARY CONDITION VALUES */
	for(j=0;j<base.nseg;++j) {
		sind = base.seg(j);
		v0 = x.seg(sind).pnt(0);
		x.ug.v(v0,0) = 0.0;
		x.ug.v(v0,2) = 0.0;
	}
	v0 = x.seg(sind).pnt(1);
	x.ug.v(v0,0) = 0.0;
	x.ug.v(v0,2) = 0.0;

	/*******************/    
	/* SET SIDE VALUES */
	/*******************/
	for(j=0;j<base.nseg;++j) {
		sind = base.seg(j);
		indx = sind*x.sm0;
		for(m=0;m<basis::tri(x.log2p)->sm();++m) {
			x.ug.s(sind,m,0) = 0.0;
			x.ug.s(sind,m,2) = 0.0;
		}
	}

	return;
}

void actuator_disc::output(const std::string& filename, tri_hp::filetype typ,int tlvl) {
    int n,ind,sind;
    TinyVector<FLT,tri_mesh::ND> nrm, mvel, pt;
    FLT powerP, powerQ;
    
    switch(typ) {
        case(tri_hp::tecplot): {
            if (!report_flag) return;
            
            powerP = 0.0;
            powerQ = 0.0;
            
            for (ind=0;ind<base.nseg;++ind) {
                sind = base.seg(ind);
                
                x.crdtocht1d(sind);
                for(n=0;n<tri_mesh::ND;++n)
                    basis::tri(x.log2p)->proj1d(&x.cht(n,0),&x.crd(n)(0,0),&x.dcrd(n,0)(0,0));
                
                x.ugtouht1d(sind);
                for(n=0;n<x.NV;++n)
                    basis::tri(x.log2p)->proj1d(&x.uht(n)(0),&x.u(n)(0,0));
                
                for(int k=0;k<basis::tri(x.log2p)->gpx();++k) {
                    nrm(0) = x.dcrd(1,0)(0,k);
                    nrm(1) = -x.dcrd(0,0)(0,k);
                    for(n=0;n<tri_mesh::ND;++n) {
                        pt(n) = x.crd(n)(0,k);
                        mvel(n) = x.gbl->bd(0)*(x.crd(n)(0,k) -dxdt(x.log2p,ind)(n,k));
#ifdef MESH_REF_VEL
                        mvel(n) += x.hp_gbl->mesh_ref_vel(n);
#endif
                    }
                    
                    FLT length = sqrt(nrm(0)*nrm(0) +nrm(1)*nrm(1));
                    FLT norm_vel = ((x.u(0)(0,k) -mvel(0))*nrm(0) +(x.u(1)(0,k) -mvel(1))*nrm(1))/length;
                    TinyVector<FLT,3> inpt(pt(0),pt(1),norm_vel);
                    FLT delta_p = dp.Eval(inpt,x.gbl->time);
                    FLT delta_q = dq.Eval(inpt,x.gbl->time);
                    
                    powerP +=  basis::tri(x.log2p)->wtx(k)*RAD(x.crd(0)(0,k))*length*delta_p*norm_vel;
                    powerQ +=  basis::tri(x.log2p)->wtx(k)*RAD(x.crd(0)(0,k))*length*delta_q*x.u(2)(0,k);

                }
            }
            streamsize oldprecision = (*x.gbl->log).precision(10);
            *x.gbl->log << base.idprefix << " power: " << powerP +powerQ << ' ' << powerP << ' ' << powerQ << std::endl;
            (*x.gbl->log).precision(oldprecision);
            break;
        }
        default: {
            break;
        }
    }
    
    generic::output(filename,typ,tlvl);
    
    return;
}
