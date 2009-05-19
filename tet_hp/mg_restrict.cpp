#include "tet_hp.h"
#include "hp_boundary.h"

void tet_hp::mg_restrict() {
    int i,j,k,m,n,tind,v0,indx,indx1;
    
    isfrst = true;
    coarse_flag = true;
    
    for(i=0;i<nfbd;++i)
        hp_fbdry(i)->mg_restrict();
		
	helper->mg_restrict();
        
    if(!coarse_level) {
        --log2p;
        
        /* TRANSFER IS ON FINE MESH */
        gbl->res0.v(Range(0,npnt-1),Range::all()) = gbl->res.v(Range(0,npnt-1),Range::all());

        if (basis::tet(log2p).p > 1) {
            gbl->res0.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all()) = gbl->res.e(Range(0,nseg-1),Range(0,basis::tet(log2p).em-1),Range::all());
            
            if (basis::tet(log2p).p > 2) {
          
                for(tind=0;tind<ntri;++tind) {
                    indx = 0;
                    indx1 = 0;
                    for(m=1;m<basis::tet(log2p).em;++m) {
                        for(k=0;k<basis::tet(log2p).em-m;++k) {
                            gbl->res0.f(tind,indx,Range::all()) = gbl->res.f(tind,indx1,Range::all());
                            ++indx;
                            ++indx1;
                        }
                        indx1 += basis::tet(log2p).p;
                    }
                }
				if(basis::tet(log2p).p > 3){
					for(tind=0;tind<ntet;++tind) {
						indx = 0;
						indx1 = 0;
						for(m=1;m<=basis::tet(log2p).em;++m) {
							for(k=1;k<=basis::tet(log2p).em-m;++k) {
								for(j=1;j<=basis::tet(log2p).em-m-k;++j){
									gbl->res0.i(tind,indx,Range::all()) = gbl->res.i(tind,indx1,Range::all());
									++indx;
									++indx1;
								}
								indx1 += basis::tet(log2p).p;
							}
						}
					}
				}
            }
        }
    }
    else {
//        if (mmovement == coupled_deformable) { 
//            r_tri_mesh::mg_restrict();
//        }

        tet_hp *fmesh = dynamic_cast<tet_hp *>(fine);
        
        gbl->res0.v(Range(0,npnt-1),Range::all()) = 0.0;
                
        /* LOOP THROUGH FINE VERTICES TO CALCULATE RESIDUAL  */
        for(i=0;i<fmesh->npnt;++i) {
            tind = fmesh->ccnnct(i).tet;
            for(j=0;j<4;++j) {
                v0 = tet(tind).pnt(j);
                for(n=0;n<NV;++n)
                    gbl->res0.v(v0,n) += fmesh->ccnnct(i).wt(j)*gbl->res.v(i,n);
            }
        }
        
        /* LOOP THROUGH COARSE VERTICES    */
        /* TO CALCULATE VUG ON COARSE MESH */
        for(i=0;i<npnt;++i) {
            tind = fcnnct(i).tet;

            ug.v(i,Range::all()) = 0.0;

            for(j=0;j<4;++j) {
                ug.v(i,Range::all()) += fcnnct(i).wt(j)*fmesh->ug.v(fmesh->tet(tind).pnt(j),Range::all());
            }

            vug_frst(i,Range::all()) = ug.v(i,Range::all());
        }
    }


    return;
}

