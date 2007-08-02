#include "mesh.h"
#include <utilities.h>

void mesh::coarsen_substructured(const class mesh &zx,int p) {
    int i,sind,sign,p2;

    copy(zx);
    
    p2 = p*p;
    ntri = zx.ntri/p2;
    nside = (ntri*3 +zx.sbdry(0)->nel/p)/2;
    nvrtx = zx.nvrtx -nside*(p-1) -ntri*((p-1)*(p-2))/2;
    
    for(i=0;i<zx.ntri;i+=p2) {
        td(i/p2).vrtx(0) = zx.td(i+2*p-2).vrtx(2);
        td(i/p2).vrtx(1) = zx.td(i).vrtx(0);
        td(i/p2).vrtx(2) = zx.td(i+p2-1).vrtx(1);
        sind = (zx.td(i).vrtx(1) -nvrtx)/(p-1);
        sign = ((zx.td(i).vrtx(1) -nvrtx)%(p-1) == 0 ? 1 : -1);
        td(i/p2).side(0) = sind; 
        td(i/p2).sign(0) = sign;
        sd(sind).vrtx((1-sign)/2) = td(i/p2).vrtx(1);
        sd(sind).vrtx((1+sign)/2) = td(i/p2).vrtx(2);
        
        sind = (zx.td(i+p2-1).vrtx(2) -nvrtx)/(p-1);
        sign = ((zx.td(i+p2-1).vrtx(2) -nvrtx)%(p-1) == 0 ? 1 : -1);
        td(i/p2).side(1) = sind;
        td(i/p2).sign(1) = sign;
        sd(sind).vrtx((1-sign)/2) = td(i/p2).vrtx(2);
        sd(sind).vrtx((1+sign)/2) = td(i/p2).vrtx(0);
        
        sind = (zx.td(i+2*p-2).vrtx(0) -nvrtx)/(p-1);
        sign = ((zx.td(i+2*p-2).vrtx(0) -nvrtx)%(p-1) == 0 ? 1 : -1);
        td(i/p2).side(2) = sind;
        td(i/p2).sign(2) = sign;
        sd(sind).vrtx((1-sign)/2) = td(i/p2).vrtx(0);
        sd(sind).vrtx((1+sign)/2) = td(i/p2).vrtx(1);
    }
    sbdry(0)->nel = zx.sbdry(0)->nel/p;
    i = 0;
    if (zx.sd(zx.sbdry(0)->el(0)).vrtx(0) < nvrtx) ++i;
    for(;i<zx.sbdry(0)->nel;i+=p) {
        sind = (zx.sd(zx.sbdry(0)->el(i)).vrtx(0) -nvrtx)/(p-1);
        sbdry(0)->el(i/p) = sind;
    }
    
    return;
}

void mesh::symmetrize() {
	int i,j,sind,vct,bnum;
	double vmin, vmax;
	
	for(i=0;i<ntri;++i) {
		vmin = 0.0;
		vmax = 0.0;
		for(vct=0;vct<3;++vct) {
			vmin = MIN(vmin,vrtx(td(i).vrtx(vct))(1));
			vmax = MAX(vmax,vrtx(td(i).vrtx(vct))(1));
		}
		if (vmax > fabs(vmin)) 
			td(i).info = 0;
		else
			td(i).info = 1;
	}
    
    output("testing",easymesh);
    	
	mesh zpart[2];
    zpart[0].allocate(maxvst*4);
	zpart[0].partition(*this,1);
    
    /* FIND NEW BOUNDARY */
    bnum = -1;
    for(i=0;i<zpart[0].nsbd;++i) {
        if (zpart[0].sbdry(i)->is_comm()) {
            for(j=0;j<nsbd;++j)
                if (sbdry(j)->idnum == zpart[0].sbdry(i)->idnum) goto next;
            bnum = i;
            break;
            next: continue;
        }
    }
    
    if (bnum > -1) {
        for(j=0;j<zpart[0].sbdry(bnum)->nel;++j) {
            sind = zpart[0].sbdry(bnum)->el(j);
            for(vct=0;vct<2;++vct) 
                zpart[0].vrtx(zpart[0].sd(sind).vrtx(vct))(1) = 0.0;
        }
        zpart[0].smooth_cofa(2);
    }
    zpart[0].output("cut",grid);

	zpart[1].copy(zpart[0]);
	for(i=0;i<zpart[1].nvrtx;++i)
		zpart[1].vrtx(i)(1) *= -1.0;
    
	for(i=0;i<zpart[1].nside;++i) {
		vct = zpart[1].sd(i).vrtx(0);
		zpart[1].sd(i).vrtx(0) = zpart[1].sd(i).vrtx(1);
		zpart[1].sd(i).vrtx(1) = vct;
	}

	for(i=0;i<zpart[1].ntri;++i) {
		vct = zpart[1].td(i).vrtx(0);
		zpart[1].td(i).vrtx(0) = zpart[1].td(i).vrtx(1);
		zpart[1].td(i).vrtx(1) = vct;
        vct = zpart[1].td(i).side(0);
        zpart[1].td(i).side(0) = zpart[1].td(i).side(1);
		zpart[1].td(i).side(1) = vct;
        vct = zpart[1].td(i).sign(0);
        zpart[1].td(i).sign(0) = zpart[1].td(i).sign(1);
		zpart[1].td(i).sign(1) = vct;
        vct = zpart[1].td(i).tri(0);
        zpart[1].td(i).tri(0) = zpart[1].td(i).tri(1);
		zpart[1].td(i).tri(1) = vct;
	}	
	
	for(i=0;i<zpart[1].nsbd;++i)
		zpart[1].sbdry(i)->reorder();
        
    zpart[1].output("reflected");

    if (bnum > -1) {
        zpart[0].append(zpart[1]);
        zpart[0].output("appended");
    }

    return;
}

void mesh::cut(Array<double,1> indicator) {
	int i,j,sind,vct,bnum;
	double vmin, vmax;
	
	for(i=0;i<ntri;++i) {
		vmin = 0.0;
		vmax = 0.0;
		for(vct=0;vct<3;++vct) {
			vmin = MIN(vmin,indicator(td(i).vrtx(vct)));
			vmax = MAX(vmax,indicator(td(i).vrtx(vct)));
		}
		if (vmax > fabs(vmin)) 
			td(i).info = 0;
		else
			td(i).info = 1;
	}
    
    output("testing",easymesh);
    	
    for (int m = 0; m < 2; ++m) {
        mesh zpart[2];
        zpart[m].allocate(maxvst*4);
        zpart[m].partition(*this,m);

        /* FIND NEW BOUNDARY */
        bnum = -1;
        for(i=0;i<zpart[m].nsbd;++i) {
            if (zpart[m].sbdry(i)->is_comm()) {
                for(j=0;j<nsbd;++j)
                    if (sbdry(j)->idnum == zpart[m].sbdry(i)->idnum) goto next;
                bnum = i;
                break;
                next: continue;
            }
        }
        
        if (bnum > -1) {
            TinyVector<FLT,ND> grad,dphi,dx;
            TinyMatrix<FLT,ND,ND> ldcrd;
            FLT jcbi,dnorm,mag;
            int n,tind,vind,tindold, vindold;
            TinyVector<int,3> v;
            
            for(j=0;j<zpart[m].sbdry(bnum)->nel;++j) {
                sind = zpart[m].sbdry(bnum)->el(j);
                tind = zpart[m].sd(sind).tri(0);
                tindold = zpart[m].td(tind).info;
                
                /* CALCULATE GRADIENT */
                v = td(tindold).vrtx;
                for(n=0;n<ND;++n) {
                    ldcrd(n,0) = 0.5*(vrtx(v(2))(n) -vrtx(v(1))(n));
                    ldcrd(n,1) = 0.5*(vrtx(v(0))(n) -vrtx(v(1))(n));
                }
                jcbi = 1./(ldcrd(0,0)*ldcrd(1,1) -ldcrd(0,1)*ldcrd(1,0));
                dphi(0) = 0.5*(indicator(v(2)) -indicator(v(1)));
                dphi(1) = 0.5*(indicator(v(0)) -indicator(v(1)));

                grad(0) =  dphi(0)*ldcrd(1,1) -dphi(1)*ldcrd(1,0);
                grad(1) = -dphi(0)*ldcrd(0,1) +dphi(1)*ldcrd(0,0);
                grad *= jcbi;
                mag = sqrt(grad(0)*grad(0) +grad(1)*grad(1));
                
                /* MOVE FIRST VERTEX OF EACH EDGE */
                /* THE OLD MESH STORES  */
                /* vd(vind).info = new vrtx index or -1 */
                /* THE NEW MESH STORES */
                /* td(tind).info = old tri index */
                vind = zpart[m].sd(sind).vrtx(0);
                for (vindold=0;vindold<3;++vindold) 
                    if (vd(td(tindold).vrtx(vindold)).info == vind) break;
                vindold = td(tindold).vrtx(vindold);
                
                dnorm = -indicator(vindold)/mag;
                dx(0) = dnorm*grad(0)/mag;
                dx(1) = dnorm*grad(1)/mag;
                
                zpart[m].vrtx(vind) += dx;
            }
        }
        char buff[100];
        sprintf(buff,"cut%d",m);
        zpart[m].output(buff,grid);
    }


    return;
}	
			
	

void mesh::shift(TinyVector<FLT,ND>& s) {
    int n;
    
    for(int i=0;i<nvrtx;++i)
        for(n=0;n<ND;++n)
            vrtx(i)(1) += s(n);

    return;
}

void mesh::scale(TinyVector<FLT,ND>& s) {
    int n;
    for(int i=0;i<nvrtx;++i)
        for(n=0;n<ND;++n)
            vrtx(i)(1) *= s(n);

    return;
}

int mesh::smooth_cofa(int niter) {
    int iter,sind,i,j,n,v0,v1;
        
    for(i=0;i<nvrtx;++i)
        vd(i).info = 0;
        
    for(i=0;i<nsbd;++i) {
        for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j);
            vd(sd(sind).vrtx(0)).info = -1;
            vd(sd(sind).vrtx(1)).info = -1;
        }
    }
    
    for(iter=0; iter< niter; ++iter) {
        /* SMOOTH POINT DISTRIBUTION X*/
        for(n=0;n<ND;++n) {
            for(i=0;i<nvrtx;++i)
                gbl_ptr->fltwk(i) = 0.0;
    
            for(i=0;i<nside;++i) {
                v0 = sd(i).vrtx(0);
                v1 = sd(i).vrtx(1);
                gbl_ptr->fltwk(v0) += vrtx(v1)(n);
                gbl_ptr->fltwk(v1) += vrtx(v0)(n);
            }
    
            for(i=0;i<nvrtx;++i) {
                if (vd(i).info == 0) {
                    vrtx(i)(n) = gbl_ptr->fltwk(i)/vd(i).nnbor;
                }
            }
        }
    }
    
    return(1);
}
