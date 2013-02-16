#include "tri_mesh.h"
#include <utilities.h>

void tri_mesh::coarsen_substructured(const class tri_mesh &zx,int p) {
	int i,sind,sgn,p2;

	copy(zx);

	p2 = p*p;
	ntri = zx.ntri/p2;
	nseg = (ntri*3 +zx.ebdry(0)->nseg/p)/2;
	npnt = zx.npnt -nseg*(p-1) -ntri*((p-1)*(p-2))/2;

	for(i=0;i<zx.ntri;i+=p2) {
		tri(i/p2).pnt(0) = zx.tri(i+2*p-2).pnt(2);
		tri(i/p2).pnt(1) = zx.tri(i).pnt(0);
		tri(i/p2).pnt(2) = zx.tri(i+p2-1).pnt(1);
		sind = (zx.tri(i).pnt(1) -npnt)/(p-1);
		sgn = ((zx.tri(i).pnt(1) -npnt)%(p-1) == 0 ? 1 : -1);
		tri(i/p2).seg(0) = sind;
		tri(i/p2).sgn(0) = sgn;
		seg(sind).pnt((1-sgn)/2) = tri(i/p2).pnt(1);
		seg(sind).pnt((1+sgn)/2) = tri(i/p2).pnt(2);

		sind = (zx.tri(i+p2-1).pnt(2) -npnt)/(p-1);
		sgn = ((zx.tri(i+p2-1).pnt(2) -npnt)%(p-1) == 0 ? 1 : -1);
		tri(i/p2).seg(1) = sind;
		tri(i/p2).sgn(1) = sgn;
		seg(sind).pnt((1-sgn)/2) = tri(i/p2).pnt(2);
		seg(sind).pnt((1+sgn)/2) = tri(i/p2).pnt(0);

		sind = (zx.tri(i+2*p-2).pnt(0) -npnt)/(p-1);
		sgn = ((zx.tri(i+2*p-2).pnt(0) -npnt)%(p-1) == 0 ? 1 : -1);
		tri(i/p2).seg(2) = sind;
		tri(i/p2).sgn(2) = sgn;
		seg(sind).pnt((1-sgn)/2) = tri(i/p2).pnt(0);
		seg(sind).pnt((1+sgn)/2) = tri(i/p2).pnt(1);
	}
	ebdry(0)->nseg = zx.ebdry(0)->nseg/p;
	i = 0;
	if (zx.seg(zx.ebdry(0)->seg(0)).pnt(0) < npnt) ++i;
	for(;i<zx.ebdry(0)->nseg;i+=p) {
		sind = (zx.seg(zx.ebdry(0)->seg(i)).pnt(0) -npnt)/(p-1);
		ebdry(0)->seg(i/p) = sind;
	}

	return;
}

void tri_mesh::symmetrize() {
	int i,j,sind,vct,bnum;
	double vmin, vmax;

	for(i=0;i<ntri;++i) {
		vmin = 0.0;
		vmax = 0.0;
		for(vct=0;vct<3;++vct) {
			vmin = MIN(vmin,pnts(tri(i).pnt(vct))(1));
			vmax = MAX(vmax,pnts(tri(i).pnt(vct))(1));
		}
		if (vmax > fabs(vmin))
			tri(i).info = 0;
		else
			tri(i).info = 1;
	}

	output("testing",easymesh);

	tri_mesh zpart[2];
	zpart[0].allocate(maxpst*4);
	zpart[0].partition(*this,1);

	/* FIND NEW BOUNDARY */
	bnum = -1;
	for(i=0;i<zpart[0].nebd;++i) {
		if (zpart[0].ebdry(i)->is_comm()) {
			for(j=0;j<nebd;++j)
				if (ebdry(j)->idnum == zpart[0].ebdry(i)->idnum) goto next;
			bnum = i;
			break;
			next: continue;
		}
	}

	if (bnum > -1) {
		for(j=0;j<zpart[0].ebdry(bnum)->nseg;++j) {
			sind = zpart[0].ebdry(bnum)->seg(j);
			for(vct=0;vct<2;++vct)
				zpart[0].pnts(zpart[0].seg(sind).pnt(vct))(1) = 0.0;
		}
		zpart[0].smooth_cofa(2);
	}
	zpart[0].output("cut",grid);

	zpart[1].copy(zpart[0]);
	for(i=0;i<zpart[1].npnt;++i)
		zpart[1].pnts(i)(1) *= -1.0;

	for(i=0;i<zpart[1].nseg;++i) {
		vct = zpart[1].seg(i).pnt(0);
		zpart[1].seg(i).pnt(0) = zpart[1].seg(i).pnt(1);
		zpart[1].seg(i).pnt(1) = vct;
	}

	for(i=0;i<zpart[1].ntri;++i) {
		vct = zpart[1].tri(i).pnt(0);
		zpart[1].tri(i).pnt(0) = zpart[1].tri(i).pnt(1);
		zpart[1].tri(i).pnt(1) = vct;
		vct = zpart[1].tri(i).seg(0);
		zpart[1].tri(i).seg(0) = zpart[1].tri(i).seg(1);
		zpart[1].tri(i).seg(1) = vct;
		vct = zpart[1].tri(i).sgn(0);
		zpart[1].tri(i).sgn(0) = zpart[1].tri(i).sgn(1);
		zpart[1].tri(i).sgn(1) = vct;
		vct = zpart[1].tri(i).tri(0);
		zpart[1].tri(i).tri(0) = zpart[1].tri(i).tri(1);
		zpart[1].tri(i).tri(1) = vct;
	}

	for(i=0;i<zpart[1].nebd;++i)
		zpart[1].ebdry(i)->reorder();

	zpart[1].output("reflected");

	if (bnum > -1) {
		zpart[0].append(zpart[1]);
		zpart[0].output("appended");
	}

	return;
}

void tri_mesh::cut() {
	int i,j,sind,vct,bnum;
	double vmin, vmax;

	for(i=0;i<ntri;++i) {
		vmin = 0.0;
		vmax = 0.0;
		for(vct=0;vct<3;++vct) {
			vmin = MIN(vmin,gbl->fltwk(tri(i).pnt(vct)));
			vmax = MAX(vmax,gbl->fltwk(tri(i).pnt(vct)));
		}
		if (vmax > fabs(vmin))
			tri(i).info = 0;
		else
			tri(i).info = 1;
	}

	output("testing",easymesh);

	for (int m = 0; m < 2; ++m) {
		tri_mesh zpart[2];
		zpart[m].allocate(maxpst*4);
		zpart[m].partition(*this,m);

		char buff[100];
		sprintf(buff,"begin%d",m);
		zpart[m].output(buff,tecplot);

		/* FIND NEW BOUNDARY */
		bnum = -1;
		for(i=0;i<zpart[m].nebd;++i) {
			if (zpart[m].ebdry(i)->is_comm()) {
				for(j=0;j<nebd;++j)
					if (ebdry(j)->idnum == zpart[m].ebdry(i)->idnum) goto next;
				bnum = i;
				break;
				next: continue;
			}
		}

		if (bnum > -1) {
			TinyVector<FLT,ND> grad,dphi,dx;
			TinyMatrix<FLT,ND,ND> ldcrd;
			FLT jcbi,dnorm,mag;
			int n,tind,pind,tindold, vindold;
			TinyVector<int,3> v;

			for(j=0;j<zpart[m].ebdry(bnum)->nseg;++j) {
				sind = zpart[m].ebdry(bnum)->seg(j);
				tind = zpart[m].seg(sind).tri(0);
				tindold = zpart[m].tri(tind).info;

				/* CALCULATE GRADIENT */
				v = tri(tindold).pnt;
				for(n=0;n<ND;++n) {
					ldcrd(n,0) = 0.5*(pnts(v(2))(n) -pnts(v(1))(n));
					ldcrd(n,1) = 0.5*(pnts(v(0))(n) -pnts(v(1))(n));
				}
				jcbi = 1./(ldcrd(0,0)*ldcrd(1,1) -ldcrd(0,1)*ldcrd(1,0));
				dphi(0) = 0.5*(gbl->fltwk(v(2)) -gbl->fltwk(v(1)));
				dphi(1) = 0.5*(gbl->fltwk(v(0)) -gbl->fltwk(v(1)));

				grad(0) =  dphi(0)*ldcrd(1,1) -dphi(1)*ldcrd(1,0);
				grad(1) = -dphi(0)*ldcrd(0,1) +dphi(1)*ldcrd(0,0);
				grad *= jcbi;
				mag = sqrt(grad(0)*grad(0) +grad(1)*grad(1));

				/* MOVE FIRST VERTEX OF EACH EDGE */
				/* THE OLD MESH STORES  */
				/* pnt(pind).info = new pnt index or -1 */
				/* THE NEW MESH STORES */
				/* tri(tind).info = old tri index */
				pind = zpart[m].seg(sind).pnt(0);
				for (vindold=0;vindold<3;++vindold)
					if (pnt(tri(tindold).pnt(vindold)).info == pind) break;
				vindold = tri(tindold).pnt(vindold);

				dnorm = -gbl->fltwk(vindold)/mag;
				dx(0) = dnorm*grad(0)/mag;
				dx(1) = dnorm*grad(1)/mag;

				zpart[m].pnts(pind) += dx;
			}
		}
		sprintf(buff,"cut%d",m);
		zpart[m].output(buff,grid);
	}


	return;
}

void tri_mesh::trim() {
	int i,j,n,bsd,tin,tind,nsrch,ntdel;

	/* ASSUMES gbl->fltwk HAS BEEN SET WITH VALUES TO DETERMINE HOW MUCH TO TRIM OFF OF BOUNDARIES */

	for(i=0;i<ntri;++i)
		tri(i).info = 0;

	ntdel = 0;

	for (bsd=0;bsd<ebdry(0)->nseg;++bsd) {
		tind = seg(ebdry(0)->seg(bsd)).tri(0);
		if (tri(tind).info > 0) continue;

		gbl->intwk(0) = tind;
		tri(tind).info = 1;
		nsrch = ntdel+1;

		/* NEED TO SEARCH SURROUNDING TRIANGLES */
		for(i=ntdel;i<nsrch;++i) {
			tin = gbl->intwk(i);
			for (n=0;n<3;++n)
				if (gbl->fltwk(tri(tin).pnt(n)) < 0.0) goto NEXT;

			gbl->intwk(ntdel++) = tin;

			for(j=0;j<3;++j) {
				tind = tri(tin).tri(j);
				if (tind < 0) continue;
				if (tri(tind).info > 0) continue;
				tri(tind).info = 1;
				gbl->intwk(nsrch++) = tind;
			}
			NEXT: continue;
		}
	}

	for(i=0;i<ntri;++i)
		tri(i).info = 0;

	for(i=0;i<ntdel;++i)
		tri(gbl->intwk(i)).info = 1;

	for(i=0;i<maxpst;++i)
		gbl->intwk(i) = -1;

	tri_mesh ztrim;
	ztrim.init(*this);
	ztrim.partition(*this,0);
	ztrim.output("trim",grid);
	ztrim.output("trim",tri_mesh::boundary);

	return;
}





void tri_mesh::shift(TinyVector<FLT,ND>& s) {
	int n;

	for(int i=0;i<npnt;++i)
		for(n=0;n<ND;++n)
			pnts(i)(n) += s(n);

	return;
}

void tri_mesh::scale(TinyVector<FLT,ND>& s) {
	int n;
	for(int i=0;i<npnt;++i)
		for(n=0;n<ND;++n)
			pnts(i)(n) *= s(n);

	return;
}

int tri_mesh::smooth_cofa(int niter) {
	int iter,sind,i,j,n,p0,p1;

	for(i=0;i<npnt;++i)
		pnt(i).info = 0;

	for(i=0;i<nebd;++i) {
		for(j=0;j<ebdry(i)->nseg;++j) {
			sind = ebdry(i)->seg(j);
			pnt(seg(sind).pnt(0)).info = -1;
			pnt(seg(sind).pnt(1)).info = -1;
		}
	}

	for(iter=0; iter< niter; ++iter) {
		/* SMOOTH POINT DISTRIBUTION X*/
		for(n=0;n<ND;++n) {
			for(i=0;i<npnt;++i)
				gbl->fltwk(i) = 0.0;

			for(i=0;i<nseg;++i) {
				p0 = seg(i).pnt(0);
				p1 = seg(i).pnt(1);
				gbl->fltwk(p0) += pnts(p1)(n);
				gbl->fltwk(p1) += pnts(p0)(n);
			}

			for(i=0;i<npnt;++i) {
				if (pnt(i).info == 0) {
					pnts(i)(n) = gbl->fltwk(i)/pnt(i).nnbor;
				}
			}
		}
	}

	return(1);
}

int tri_mesh::smooth_lngth(int niter) {
	int iter,i,p0,p1;
	int sind, j;

	/* Smooth distribution along each edge fixing endpoint */
//	for(iter=0; iter< niter; ++iter) {
//		for(i=0;i<nebd;++i) {
//			gbl->fltwk(0) = 0.0;
//			for(j=0;j<ebdry(i)->nseg;++j) {
//				gbl->fltwk(j) += 1.0/lngth(seg(ebdry(i)->seg(j)).pnt(1));
//				gbl->fltwk(j+1) = 1.0/lngth(seg(ebdry(i)->seg(j)).pnt(0));
//			}
//			
//			for(j=1;j<ebdry(i)->nseg;++j) {
//				lngth(seg(ebdry(i)->seg(j)).pnt(0)) = 2.0/gbl->fltwk(j);
//			}
//		}
//	}
//	
//	for(i=0;i<npnt;++i)
//		pnt(i).info = 0;
//
//	for(i=0;i<nebd;++i) {
//		for(j=0;j<ebdry(i)->nseg;++j) {
//			sind = ebdry(i)->seg(j);
//			pnt(seg(sind).pnt(0)).info = -1;
//			pnt(seg(sind).pnt(1)).info = -1;
//		}
//	}

	for(iter=0; iter< niter; ++iter) {
		/* SMOOTH POINT DISTRIBUTION X*/
			for(i=0;i<npnt;++i)
				gbl->fltwk(i) = 0.0;

			for(i=0;i<nseg;++i) {
				p0 = seg(i).pnt(0);
				p1 = seg(i).pnt(1);
				gbl->fltwk(p0) += 1./lngth(p1);
				gbl->fltwk(p1) += 1./lngth(p0);
			}

			for(i=0;i<npnt;++i) {
//				if (pnt(i).info == 0) {
					lngth(i) = 1./(gbl->fltwk(i)/pnt(i).nnbor);
//				}
			}
	}

	return(1);
}

