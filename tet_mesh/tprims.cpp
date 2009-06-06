#include"tet_mesh.h"
#include<cmath>

bool tet_mesh::findtet(const TinyVector<FLT,3> xp, int seedvrtx, int &tind) {
	tind = pnt(seedvrtx).tet;
	int nbor = pnt(seedvrtx).nnbor;
	int ind = 0;
	int ind2 = nbor;
	bool found = true;

	/* KNOWN TET CONNECTED TO VERTEX */
	if(intet(tind, xp) < tet(tind).vol*10.*EPSILON)
		goto FOUNDTET;
	
	gbl->i2wk(ind++) = tind;
	gbl->i1wk(tind) = 0;

	/* SEARCH FOR POINT IN VERTEX BALL FIRST */
	for(int i = 0; i < nbor; ++i) {	
		for(int j = 0; j < 4; ++j) {			
			tind = tet(gbl->i2wk(i)).tet(j);
			if (tind == -1)
				goto NEXTFACE;
			if (gbl->i1wk(tind) == -1) {			
				for(int k = 0; k < 4; ++k) {
					if(tet(tind).pnt(k) == seedvrtx) {
						if(intet(tind, xp) < tet(tind).vol*10.*EPSILON)
							goto FOUNDTET;
						gbl->i2wk(ind++) = tind; 
						gbl->i1wk(tind) = 0;
						goto NEXTFACE;							
					}
				}
				/* STORE SOME SURROUNDING TETS TO SEARCH LATER */
				gbl->i2wk(ind2++) = tind; 
				gbl->i1wk(tind) = 0;				
			}
			NEXTFACE:;		
		}
	}

	/* SEARCH SURROUNDING TETS */
	for(int i = nbor; i < ntet; ++i) {
		if(intet(tind, xp) < tet(tind).vol*10.*EPSILON)
			goto FOUNDTET;
		for(int j = 0; j < 4; ++j) {
			tind = tet(gbl->i2wk(i)).tet(j);
			if(tind > -1 && gbl->i1wk(tind) == -1) {
				if(intet(tind, xp) < tet(tind).vol*10.*EPSILON)
					goto FOUNDTET;
				gbl->i2wk(ind2++) = tind;
				gbl->i1wk(tind) = 0;
			}
		}		
	}
	
	/* SEARCHED EVERY TET AND FAILED: RETURN 0 */
	*gbl->log << "Error: could not find tet that contains point at " << xp << " with seed vertex " << seedvrtx << std::endl;     
	found = false;
	
	FOUNDTET:;
	for(int i = 0; i < ind; ++i)
		gbl->i1wk(gbl->i2wk(i)) = -1; /* reset i1wk to -1 */
	for(int i = nbor; i < ind2; ++i)
		gbl->i1wk(gbl->i2wk(i)) = -1; /* reset i1wk to -1 */	

	return(found);
}

bool tet_mesh::findtet(const TinyVector<FLT,3> xp, int &tind){
	bool found = 1;
	int ind = 0;

	/* KNOWN TET CONNECTED TO VERTEX */
	if(intet(tind, xp) < tet(tind).vol*10.*EPSILON)
		goto FOUNDTET;
	
	gbl->i2wk(ind++) = tind;
	gbl->i1wk(tind) = 0;

	/* SEARCH SURROUNDING TETS */
	for(int i = 0; i < ntet; ++i) {
		for(int j = 0; j < 4; ++j) {
			tind = tet(gbl->i2wk(i)).tet(j);
			if(tind > -1 && gbl->i1wk(tind) == -1) {
				if(intet(tind, xp) < tet(tind).vol*10.*EPSILON)
					goto FOUNDTET;
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
	
	double rnx,rny,rnz;
	TinyVector<int,4> v;	
	TinyMatrix<FLT,5,3> pnt;
	
	for(int i=0;i < 3; ++i){
		pnt(4,i)=xp(i);
	}
	
	/* LOAD IN VERTICES */
	v=tet(tind).pnt;
	for(int i=0;i < 4; ++i){
		for(int j=0;j < 3; ++j){
			pnt(i,j)=pnts(v(i))(j); 
		}
	}
	
	/* FIND VOLUME OF TET CREATED FROM A FACE ORIGINAL TET AND A POINT */
	for(int i = 0; i < 4; ++i){
		v(0)=0;
		v(1)=1;
		v(2)=2;
		v(3)=3;
		v(i)=4;
		rnx = (pnt(v(1),1)-pnt(v(0),1))*(pnt(v(2),2)-pnt(v(0),2))-(pnt(v(2),1)-pnt(v(0),1))*(pnt(v(1),2)-pnt(v(0),2));
		rny = (pnt(v(1),2)-pnt(v(0),2))*(pnt(v(2),0)-pnt(v(0),0))-(pnt(v(2),2)-pnt(v(0),2))*(pnt(v(1),0)-pnt(v(0),0));
		rnz = (pnt(v(1),0)-pnt(v(0),0))*(pnt(v(2),1)-pnt(v(0),1))-(pnt(v(2),0)-pnt(v(0),0))*(pnt(v(1),1)-pnt(v(0),1));
		/* CALCULATE VOLUME */
		tet_wgt(i)=-rnx*(pnt(v(3),0)-pnt(v(0),0))-rny*(pnt(v(3),1)-pnt(v(0),1))-rnz*(pnt(v(3),2)-pnt(v(0),2)); 
	}
	
	return(fabs(tet_wgt(0)) +fabs(tet_wgt(1)) +fabs(tet_wgt(2))+fabs(tet_wgt(3)) - (tet_wgt(0) +tet_wgt(1) +tet_wgt(2)+tet_wgt(3)));
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

