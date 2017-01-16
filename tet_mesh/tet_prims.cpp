#include "tet_mesh.h"
//#include <blitz/tinyvec-et.h>
#include <cmath>

bool tet_mesh::findtet(const TinyVector<FLT,3> xp, int seedvrtx, int &tind) {
	tind = pnt(seedvrtx).tet;
	int nbor = pnt(seedvrtx).nnbor;
	int ind = 0;
	int ind2 = nbor;
	bool found = true;
	
	/* ADD TO LIST & MARK ADDED */
	gbl->i2wk(ind++) = tind;
	gbl->i1wk(tind) = 0;
	
	/* SEARCH FOR POINT IN VERTEX BALL FIRST */
	for(int i = 0; i < nbor; ++i) {	
		tind = gbl->i2wk(i);
		if (intet(tind, xp) < 10.*EPSILON*tet(tind).vol)
			goto FOUNDTET;
				
		/* FIND NEXT TETS */
		for(int j = 0; j < 4; ++j) {			
			tind = tet(gbl->i2wk(i)).tet(j);
			if (tind > -1 && gbl->i1wk(tind) == -1) {			
				if ((tet(tind).pnt(0)-seedvrtx)*(tet(tind).pnt(1)-seedvrtx)*(tet(tind).pnt(2)-seedvrtx)*(tet(tind).pnt(3)-seedvrtx)  == 0) {
					gbl->i2wk(ind++) = tind;
				}
				else {
					/* ADD TO BE SEARCHED LATER */
					gbl->i2wk(ind2++) = tind;
				}
				gbl->i1wk(tind) = 0;  // MARK ADDED
			}
		}
	}
	if (ind != nbor) {
		*gbl->log << "problem in findtet " << seedvrtx << ' ' << nbor << '\n';
		*gbl->log << "tets with common vertex" << std::endl;
		for (int i=0;i<ind;++i) {
			*gbl->log << gbl->i2wk(i);
			tind = gbl->i2wk(i);
			for (int j=0;j<4;++j)
				*gbl->log << '\t' << tet(tind).tet(j) << ',' << gbl->i1wk(tet(tind).tet(j));
			*gbl->log << std::endl;
		}
		*gbl->log << std::endl;
		*gbl->log << "other tets" << std::endl;

		for (int i=nbor;i<ind2;++i) {
			*gbl->log << gbl->i2wk(i);
			tind = gbl->i2wk(i);
			for (int j=0;j<4;++j)
				*gbl->log << '\t' << tet(tind).tet(j) << ',' << gbl->i1wk(tet(tind).tet(j));
			*gbl->log << std::endl;
		}
			
		output("error"+gbl->idprefix,easymesh);
		output("error"+gbl->idprefix);
		
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
			
	/* SEARCH SURROUNDING TETS */
	for(int i = nbor; i < ntet; ++i) {
		tind = gbl->i2wk(i);
		if (intet(tind, xp) < 10.*EPSILON*tet(tind).vol)
			goto FOUNDTET;
		
		/* ADD NEIGHBORS */
		for(int j = 0; j < 4; ++j) {
			tind = tet(gbl->i2wk(i)).tet(j);
			if (tind > -1 && gbl->i1wk(tind) == -1) {
				gbl->i2wk(ind2++) = tind;
				gbl->i1wk(tind) = 0;
			}
		}		
	}
	
	/* SEARCHED EVERY TET AND FAILED: RETURN 0 */
//	*gbl->log << "Error: could not find tet that contains point at " << xp << " with seed vertex " << seedvrtx << std::endl;     
//	for (int i=0;i<ind2;++i) {
//		tind = gbl->i2wk(i);
//		*gbl->log << i << ' ' << gbl->i2wk(i) << ' ' << intet(tind, xp) << std::endl;
//	}
//	output("error",easymesh);
//	output("error");
//	sim::abort(__LINE__,__FILE__,gbl->log);

	found = false;
	
	FOUNDTET:;	
	for(int i = 0; i < ind; ++i)
		gbl->i1wk(gbl->i2wk(i)) = -1; /* reset i1wk to -1 */
	for(int i = nbor; i < ind2; ++i)
		gbl->i1wk(gbl->i2wk(i)) = -1; /* reset i1wk to -1 */

	return(found);
}

bool tet_mesh::findtet(const TinyVector<FLT,3> xp, int &tind) {
	bool found = 1;
	int ind = 0;
	
	gbl->i2wk(ind++) = tind;
	gbl->i1wk(tind) = 0;
	
	/* SEARCH TET */
	for(int i = 0; i < ntet; ++i) {
		tind = gbl->i2wk(i);
		if (intet(tind, xp) < 10.*EPSILON*tet(tind).vol)
			goto FOUNDTET;
		
		/* ADD NEIGHBORS */
		for(int j = 0; j < 4; ++j) {
			tind = tet(gbl->i2wk(i)).tet(j);
			if (tind > -1 && gbl->i1wk(tind) == -1) {
				gbl->i2wk(ind++) = tind;
				gbl->i1wk(tind) = 0;
			}
		}		
	}
	
	/* SEARCHED EVERY TET AND FAILED: RETURN -1 */
	*gbl->log << "Error: could not find tet that contains point at " << xp << std::endl;     
	found = 0;
	
	FOUNDTET:;
	for(int i = 0; i < ind; ++i)
		gbl->i1wk(gbl->i2wk(i)) = -1; /* reset i1wk to -1 */
	
	/* RETURN FOUND TET */
	return(found);
}

FLT tet_mesh::intet(int tind, const TinyVector<FLT,3> &xp) {
	
	/* FIND VOLUME OF TET CREATED FROM A FACE OF ORIGINAL TET AND A POINT */
	for(int i = 0; i < 4; ++i) {
		/* CALCULATE VOLUME */
		tet_wgt(i) = -tet(tind).rot(i)*volume(xp,tet(tind).tri(i)); 
	}
	
	return(fabs(tet_wgt(0)) +fabs(tet_wgt(1)) +fabs(tet_wgt(2))+fabs(tet_wgt(3)) - (tet_wgt(0) +tet_wgt(1) +tet_wgt(2)+tet_wgt(3)));
}

int tet_mesh::facematch(TinyVector<int,3> v1, TinyVector<int,3> v2) {
	/* Proof by Michael Felland Clarkson University 8/1/07
	 Face matching algorithm: The only way to match is if both sets are identical. 
	 Let a be the common value--call the sets {a,b,c} and {a,d,e}.  You want a*b*c=a*d*e 
	 and a+b+c=a+d+e; a is not 0.  This means b*c=d*e and b+c=d+e.  Since the sums are the 
	 same, one of the two pairs of values must lie between the other--say b<d<e<c. If m
	 is the common mean value, m=(b+c)/2=(d+e)/2, then there are x and y so that
	 d=m-x,e=m+x,b=m-y,c=m+y with x<y.  d*e = m^2-x^2 > m^2-y^2 = b*c (this is if
	 all are positive or all negative--cases where some are negative can be also handled). */
	
	int lcl0;
	long lcl1,lcl2;
	
	lcl0=v1(0)+v1(1)+v1(2);
	lcl1=(v1(0)+1)*(v1(1)+1)*(v1(2)+1);
	lcl2=abs(lcl0-v2(0)-v2(1)-v2(2));
	lcl2+=abs(lcl1-(v2(0)+1)*(v2(1)+1)*(v2(2)+1));	
	
	return(lcl2);
}

/* RETURNS WEIGHTS FROM INTRI FUNCTION */
void tet_mesh::getwgts(TinyVector<FLT,4> &wt) const {
	int i;
	FLT sum;
	
	sum = tet_wgt(0) +tet_wgt(1) +tet_wgt(2)+tet_wgt(3);
	for(i=0;i<4;++i) 
		wt(i) = tet_wgt(i)/sum;
	
	return;
}

FLT tet_mesh::volume(int ttind) {
	TinyVector<int,4> v = tet(ttind).pnt;
	TinyVector<FLT,ND> dx1, dx2, dx3, nrm;
	FLT vol;
	
	dx1 = pnts(v(3))-pnts(v(2));
	dx2 = pnts(v(1))-pnts(v(2));
	dx3 = pnts(v(0))-pnts(v(2));
	nrm = cross(dx1,dx2);
	vol = dot(nrm,dx3);
	return(vol);
}

FLT tet_mesh::volume(TinyVector<FLT,ND> xp, int tind) {
	TinyVector<int,3> v = tri(tind).pnt;
	TinyVector<FLT,ND> dx1, dx2, dx3, nrm;
	FLT vol;
	
	dx1 = pnts(v(2))-pnts(v(1));
	dx2 = pnts(v(0))-pnts(v(1));
	dx3 = xp-pnts(v(1));
	nrm = cross(dx1,dx2);
	vol = dot(nrm,dx3);
	return(vol);
}


void tet_mesh::vertexball(int vind) {
	int i,j,k,ind,tind,tind2;
	int nbor = pnt(vind).nnbor;        

	ind = 0;
	// known tet connected to vertex
	gbl->i2wk(ind) = pnt(vind).tet;
	gbl->i1wk(gbl->i2wk(ind)) = 0;
	
	for(i = 0; i < nbor; ++i) {    
		tind = gbl->i2wk(i);  
		for(j = 0; j < 4; ++j) { 
			tind2 = tet(tind).tet(j);
			if (tind2 < 0)
				goto NEXTFACE;
			if (gbl->i1wk(tind2) < 0) { 
				for(k = 0; k < 4; ++k) {
					if (tet(tind2).pnt(k) == vind) {
						gbl->i2wk(++ind) = tind2; // connected tet found
						gbl->i1wk(tind2) = 0;
						goto NEXTFACE;                            
					}
				}
			}
			NEXTFACE:;        
		}
	}

	for(i = 0; i < nbor; ++i) {    
		gbl->i1wk(gbl->i2wk(i))=-1; // reset i1wk to -1
	}

	return;
}

//void tet_mesh::spokes(int vind) {
//    int tind,sind;
//    int ind = 0;
//    
//    for(int i=0; i < pnt(vind).nnbor; ++i) {
//        tind = gbl->i2wk(i);
//        for(int j=0; j < 6; ++j) {
//            sind=tet(tind).seg(j);
//            if (gbl->i1wk(sind) < 0) {
//                for(int k=0; k < 2; ++k) {
//                    if (vind == seg(sind).pnt(k)) {
//                        gbl->i3wk(ind++)=sind;
//                        gbl->i1wk(sind) = 1;
//                    }
//                }
//            }
//            else {
//                ++gbl->i1wk(sind);
//            }
//        }
//    }
//    
//    nspk=ind;
//    for(int i = 0; i < nspk; ++i) {    
//        gbl->i2wk(i)=gbl->i1wk(gbl->i3wk(i));//number of times edge is hit
//        gbl->i1wk(gbl->i3wk(i))=-1; // reset i1wk to -1
//    }    
//
//    return;
//}

void tet_mesh::ring(int eind) {
	int i,j,k,ind,tind,tind2;
	int nbor = seg(eind).nnbor;        

	ind = 0;
	// known tet connected to edge
	gbl->i2wk(ind) = seg(eind).tet;
	gbl->i1wk(gbl->i2wk(ind)) = 0;
	
	for(i = 0; i < nbor; ++i) {    
		tind = gbl->i2wk(i);  
		for(j = 0; j < 4; ++j) {            
			tind2 = tet(tind).tet(j);
			if (tind2 < 0)
				goto NEXTFACE;
			if (gbl->i1wk(tind2) < 0) {            
				for(k = 0; k < 6; ++k) {
					if (tet(tind2).seg(k) == eind) {
						gbl->i2wk(++ind) = tind2; // connected tet found
						gbl->i1wk(tind2) = 0;
						goto NEXTFACE;                            
					}
				}
			}
			NEXTFACE:;        
		}
	}
	
	for(i = 0; i < nbor; ++i) {    
		gbl->i1wk(gbl->i2wk(i))=-1; // reset i1wk to -1
	}
	return;
}


void tet_mesh::switch_edge_sign(int eind) {
	int i,j,k,tind,find,sind;
	int nbor = seg(eind).nnbor; 
	
	ring(eind);
	
	for(i = 0; i < nbor; ++i) {
		tind = gbl->i2wk(i);
		for(j=0;j<6;++j) {
			sind=tet(tind).seg(j);
			if (eind == sind) {
				tet(tind).sgn(j) *= -1;
				//break;
			}
		}
		for(j=0;j<4;++j) {
			find = tet(tind).tri(j);
			for(k=0;k<3;++k) {
				sind = tri(find).seg(k);
				if (eind == sind) {
					tri(find).sgn(k) *= -1;
					//break;
				}	
			}
		}
	}
	return;
}

FLT tet_mesh::area(int find) {
	int N1 = tri(find).pnt(0);
	int N2 = tri(find).pnt(1);
	int N3 = tri(find).pnt(2);
	double ax, ay, az;
	
	/* Calculate normal to face */
	ax = (pnts(N2)(1)-pnts(N1)(1))*(pnts(N3)(2)-pnts(N1)(2))-(pnts(N2)(2)-pnts(N1)(2))*(pnts(N3)(1)-pnts(N1)(1));
	ay = (pnts(N2)(2)-pnts(N1)(2))*(pnts(N3)(0)-pnts(N1)(0))-(pnts(N2)(0)-pnts(N1)(0))*(pnts(N3)(2)-pnts(N1)(2));
	az = (pnts(N2)(0)-pnts(N1)(0))*(pnts(N3)(1)-pnts(N1)(1))-(pnts(N2)(1)-pnts(N1)(1))*(pnts(N3)(0)-pnts(N1)(0));
	
	return(sqrt(ax*ax +ay*ay +az*az));
}

void tet_mesh::face_normal(int find,TinyVector<FLT,ND> normal) {
	int N1 = tri(find).pnt(0);
	int N2 = tri(find).pnt(1);
	int N3 = tri(find).pnt(2);
	
	/* Calculate normal to face */
	normal(0) = (pnts(N2)(1)-pnts(N1)(1))*(pnts(N3)(2)-pnts(N1)(2))-(pnts(N2)(2)-pnts(N1)(2))*(pnts(N3)(1)-pnts(N1)(1));
	normal(1) = (pnts(N2)(2)-pnts(N1)(2))*(pnts(N3)(0)-pnts(N1)(0))-(pnts(N2)(0)-pnts(N1)(0))*(pnts(N3)(2)-pnts(N1)(2));
	normal(2) = (pnts(N2)(0)-pnts(N1)(0))*(pnts(N3)(1)-pnts(N1)(1))-(pnts(N2)(1)-pnts(N1)(1))*(pnts(N3)(0)-pnts(N1)(0));
	
	return;
}

 //       
//FLT tri_mesh::minangle(int p0, int p1, int p2) const {
//    int i, i1, i2;
//    TinyVector<FLT,3> l,dx,dy;
//    FLT crossprod;
//    
//    dx(0) = pnts(p2)(0) -pnts(p1)(0);
//    dy(0) = pnts(p2)(1) -pnts(p1)(1);
//    l(0) = dx(0)*dx(0) +dy(0)*dy(0);    
//
//    dx(1) = pnts(p0)(0) -pnts(p2)(0);
//    dy(1) = pnts(p0)(1) -pnts(p2)(1);
//    l(1) = dx(1)*dx(1) +dy(1)*dy(1);
//    
//    dx(2) = pnts(p1)(0) -pnts(p0)(0);
//    dy(2) = pnts(p1)(1) -pnts(p0)(1);
//    l(2) = dx(2)*dx(2) +dy(2)*dy(2);
//        
//    i = (l(0) < l(1) ? 0 : 1);
//    i = (l(i) < l(2) ? i : 2);
//    i1 = (i+1)%3;
//    i2 = (i+2)%3;
//    
//    crossprod = -dx(i2)*dy(i1) +dy(i2)*dx(i1);
//
//    return(crossprod/sqrt(l(i1)*l(i2)));
//    
//}
//    
//FLT tri_mesh::angle(int p0, int p1, int p2) const {
//    TinyVector<FLT,3> l;
//    FLT dx, dy;
//    
//    dx = pnts(p1)(0) -pnts(p0)(0);
//    dy = pnts(p1)(1) -pnts(p0)(1);
//    l(0) = dx*dx +dy*dy;
//
//    dx = pnts(p2)(0) -pnts(p1)(0);
//    dy = pnts(p2)(1) -pnts(p1)(1);
//    l(1) = dx*dx +dy*dy;    
//
//    dx = pnts(p0)(0) -pnts(p2)(0);
//    dy = pnts(p0)(1) -pnts(p2)(1);
//    l(2) = dx*dx +dy*dy;  
//            
//    return((l(0) +l(1) -l(2))/(2.*sqrt(l(0)*l(1))));
//    
//}
//
//FLT tri_mesh::circumradius(int tind) const {
//    FLT alpha,beta;
//    FLT xmid1,ymid1,xmid2,ymid2,xcen,ycen;
//    FLT dx1,dy1,dx2,dy2,area;
//    int p0, p1, p2;
//    
//    p0 = tri(tind).pnt(0);
//    p1 = tri(tind).pnt(1);
//    p2 = tri(tind).pnt(2);
//    
//    dx1 =  (pnts(p0)(0)-pnts(p2)(0));
//    dy1 =  (pnts(p0)(1)-pnts(p2)(1));
//    dx2 =  (pnts(p1)(0)-pnts(p0)(0));
//    dy2 =  (pnts(p1)(1)-pnts(p0)(1));
//
//    /* RELATIVE TO POINT 0 TO AVOID ROUNDOFF */
//    xmid1 = 0.5*(pnts(p2)(0) -pnts(p0)(0));
//    ymid1 = 0.5*(pnts(p2)(1) -pnts(p0)(1));    
//    xmid2 = 0.5*(pnts(p1)(0) -pnts(p0)(0));
//    ymid2 = 0.5*(pnts(p1)(1) -pnts(p0)(1));
//        
//    area          = 1.0/(dx1*dy2 -dy1*dx2);
//    alpha         = dx2*xmid2 +dy2*ymid2;
//    beta          = dx1*xmid1 +dy1*ymid1;
//    xcen = area*(beta*dy2 -alpha*dy1);
//    ycen = area*(alpha*dx1 -beta*dx2);
//
//    return(sqrt(xcen*xcen +ycen*ycen));
//}
//
//void tri_mesh::circumcenter(int tind, TinyVector<FLT,2> &x) const {
//    FLT alpha,beta;
//    FLT xmid1,ymid1,xmid2,ymid2;
//    FLT dx1,dy1,dx2,dy2,area;
//    int p0, p1, p2;
//    
//    p0 = tri(tind).pnt(0);
//    p1 = tri(tind).pnt(1);
//    p2 = tri(tind).pnt(2);
//    
//    dx1 =  (pnts(p0)(0)-pnts(p2)(0));
//    dy1 =  (pnts(p0)(1)-pnts(p2)(1));
//    dx2 =  (pnts(p1)(0)-pnts(p0)(0));
//    dy2 =  (pnts(p1)(1)-pnts(p0)(1));
//
//    /* RELATIVE TO POINT V0 */
//    xmid1 = 0.5*(pnts(p2)(0) -pnts(p0)(0));
//    ymid1 = 0.5*(pnts(p2)(1) -pnts(p0)(1));    
//    xmid2 = 0.5*(pnts(p1)(0) -pnts(p0)(0));
//    ymid2 = 0.5*(pnts(p1)(1) -pnts(p0)(1));
//        
//    area          = 1.0/(dx1*dy2 -dy1*dx2);
//    alpha         = dx2*xmid2 +dy2*ymid2;
//    beta          = dx1*xmid1 +dy1*ymid1;
//    x(0) = area*(beta*dy2 -alpha*dy1) +pnts(p0)(0);
//    x(1) = area*(alpha*dx1 -beta*dx2) +pnts(p0)(1);
//    
//    return;
//}
//
//FLT tri_mesh::inscribedradius(int tind) const {
//    int p0,p1,p2;
//    FLT dx1,dy1,dx2,dy2,area,perim;
//    
//    p0 = tri(tind).pnt(0);
//    p1 = tri(tind).pnt(1);
//    p2 = tri(tind).pnt(2);
//    
//    dx1 =  (pnts(p0)(0)-pnts(p2)(0));
//    dy1 =  (pnts(p0)(1)-pnts(p2)(1));
//    dx2 =  (pnts(p1)(0)-pnts(p0)(0));
//    dy2 =  (pnts(p1)(1)-pnts(p0)(1));
//    
//    area = (dx1*dy2 -dy1*dx2);
//    perim = sqrt(dx1*dx1 +dy1*dy1) +sqrt(dx2*dx2 +dy2*dy2)
//        +sqrt((dx1+dx2)*(dx1+dx2) +(dy1+dy2)*(dy1+dy2));
//        
//    return(area/perim);
//}
//
//
//FLT tri_mesh::aspect(int tind) const {
//    int p0,p1,p2;
//    FLT dx1,dy1,dx2,dy2,area,perim;
//    
//    p0 = tri(tind).pnt(0);
//    p1 = tri(tind).pnt(1);
//    p2 = tri(tind).pnt(2);
//    
//    dx1 =  (pnts(p0)(0)-pnts(p2)(0));
//    dy1 =  (pnts(p0)(1)-pnts(p2)(1));
//    dx2 =  (pnts(p1)(0)-pnts(p0)(0));
//    dy2 =  (pnts(p1)(1)-pnts(p0)(1));
//    
//    area = (dx1*dy2 -dy1*dx2)*9/sqrt(3.)*4;
//    perim = sqrt(dx1*dx1 +dy1*dy1) +sqrt(dx2*dx2 +dy2*dy2)
//        +sqrt((dx1+dx2)*(dx1+dx2) +(dy1+dy2)*(dy1+dy2));
//        
//    return(area/(perim*perim));
//}

