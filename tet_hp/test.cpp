/*
 *  test.cpp
 *  tet_hp
 *
 *  Created by michael brazell on 2/19/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_hp.h"
#include "hp_boundary.h"
#include <cstring>
#include <assert.h>
#include <stdlib.h>
#include <utilities.h>
#include <myblas.h>


void tet_hp::test() {

	
	
	create_jacobian();
	
/* insert sparse matrix test */
//	TinyVector<FLT,5> b,x;
//	x(0)=1,x(1)=4,x(2)=9,x(3)=2,x(4)=67;
	
//	insert_sparse(0,0,.1);
//	insert_sparse(1,1,.2);
//	insert_sparse(2,2,.3);
//	insert_sparse(3,3,.4);
//	insert_sparse(4,4,.5);
//
//	insert_sparse(0,2,200);
//	insert_sparse(0,4,400);
//	insert_sparse(0,3,300);
//	insert_sparse(0,1,100);
//
//
	cout << ija << sa << endl;
	
	//test matrix multiply
//	for(int i = 0; i < size_sparse_matrix; ++i){
//		b(i) = sa(i)*x(i);
//		for(int j = ija(i); j < ija(i+1); ++j)
//			b(i)+=sa(j)*x(ija(j));
//	}
//		
//	cout << b << endl;
	
	exit(4);
			
	
	
	tet_mesh::checkintegrity();
	output("b0",tecplot);
	particle();
	exit(5);

	TinyMatrix<FLT,3,11> pt;
	TinyVector<TinyMatrix<FLT,MXGP,MXGP>,4> f,g;
	TinyVector<FLT,3> pnt;
	float r,s,t;
	int ind=0;
	int find;

	
	cout.precision(10);
		//	cout << npnt << ' ' << nseg << ' ' << ntri << ' ' << ntet << endl;
			/* pt probe test */
				Array<double,1> uout(1);

//			// edge 1
//			pt(0,0) = -.345676;
//			pt(1,0) = -1;
//			pt(2,0) = -1;
//
//			// edge 2
//			pt(0,1) = -.2346;
//			pt(1,1) = -pt(0,1);
//			pt(2,1) = -1;
//			
//			// edge 3
//			pt(0,2) = -1;
//			pt(1,2) = -.4537;
//			pt(2,2) = -1;
//
//			// edge 4
//			pt(0,3) = -1;
//			pt(1,3) = -.457;
//			pt(2,3) = -pt(1,3);
//
//			// edge 5
//			pt(0,4) = -1;
//			pt(1,4) = -1;
//			pt(2,4) = -.3457;
//
//			// edge 6
//			pt(0,5) = -.2346;
//			pt(1,5) = -1;
//			pt(2,5) = -pt(0,5);
//			
//			// face 0
//			pt(0,6) = -.8457;
//			pt(1,6) = -.5678;
//			pt(2,6) = -1;
//
//			// face 1
//			pt(0,7) = -.7;
//			pt(1,7) = -1;
//			pt(2,7) = -.645687;
//
//			// face 2
//			pt(0,8) = -.7;
//			pt(1,8) = -.6;
//			pt(2,8) = -1-pt(0,8)-pt(1,8);
//
//			// face 3
//			pt(0,9) = -1;
//			pt(1,9) = -.6;
//			pt(2,9) = -.46578;
//
//			// interior
//			pt(0,10) = -.34567;
//			pt(1,10) = -.235;
//			pt(2,10) = -.46578;
//
//			TinyVector<FLT,MXTM> lin;
//			lin=0.0;
//			
//			for(int td = 0; td < ntet; ++td){
//				for(int j = 0; j < 11; ++j){
//					FLT x,y,z;
//					uht(0) = 0;
//					basis::tet(log2p).ptvalues_rst(pt(0,j),pt(1,j),pt(2,j));
//					ugtouht(td,0);
//
//					basis::tet(log2p).ptprobe(1,uout.data(),&uht(0)(0),MXTM);
//					
//					for (int i = 0; i < 4; ++i)
//						lin(i) = pnts(tet(td).pnt(i))(0);			
//					basis::tet(log2p).ptprobe(1, &x, pt(0,j),pt(1,j),pt(2,j), lin.data(), MXTM);
//					for (int i = 0; i < 4; ++i)
//						lin(i) = pnts(tet(td).pnt(i))(1);
//					basis::tet(log2p).ptprobe(1, &y, pt(0,j),pt(1,j),pt(2,j), lin.data(), MXTM);
//					for (int i = 0; i < 4; ++i)
//						lin(i) = pnts(tet(td).pnt(i))(2);
//					basis::tet(log2p).ptprobe(1, &z, pt(0,j),pt(1,j),pt(2,j), lin.data(), MXTM);
//					
//					pnt(0) = x;
//					pnt(1) = y;
//					pnt(2) = z;
//					//cout << pnt << endl;
//					lcl = gbl->ibc->f(0,pnt,gbl->time) -uout(0);
//					if(fabs(lcl) > 10e-14){
//						++ind;
//						cout << " approximated - actual = " << lcl << " mode " << j << " tet " << td <<  endl;
//					}
//				}
//			}
//				
//			cout << "total bad modes: " << ind << endl;
			
				
			//tet_mesh::checkintegrity();
					//cout << ug.e(Range(0,nseg-1),Range::all(),0) << endl;




//			const std::string filename="grdtest";
//			tet_mesh::filetype in = static_cast<tet_mesh::filetype>(grid);
//			tet_mesh::output(filename,in);
//			output("testtobasis",block::display);
//			l2error(gbl->ibc);
			
			
			
			
//		cout << ug.v(Range(0,npnt-1),Range::all())<< endl;
//		cout << ug.e(Range(0,nseg-1),Range::all(),0) << endl;
//		cout << ug.f(Range(0,ntri-1),Range::all(),0) << endl;
//		cout << ug.i(Range(0,ntet-1),Range::all(),0) << endl;

//			std::ostringstream filename;
//			filename.str("");
//			filename << "refineby2" << std::flush;
//			cout << "running refineby2" << endl;
//
//			tet_mesh::output(filename.str(),tet_mesh::refineby2);
//
//			exit(4);
		
		
		
		
		
			cout << "running minvrt_test" << endl;

			tet_hp::minvrt_test();
			std::ostringstream filename;
			filename.str("");
			filename << "minvrt_test"  << std::flush;
			output(filename.str(),block::display);
			exit(5);
		
			

			
//			int td = 0;
//			
//			cout << tri(td).pnt <<' ' << tet(tri(td).tet(0)).pnt << ' ' << tet(tri(td).tet(1)).pnt << endl;
		
//		for(int i = 0; i < 4; ++i){
//			find=tet(td).tri(i);
//			ugtouht2d(find,0);
//			basis::tet(log2p).proj2d(&uht(0)(0),&f(i)(0,0), MXGP);
//			//cout << "face: "<< i << " f = " << f(i) << endl;
//		}
//		
//		for(int m = 0; m < 4; ++m){
//			find=tet(td).tri(m);
//			ugtouht2d(find,0);	
//			for(int i = 0; i < MXGP; ++i){
//				for(int j = 0; j < MXGP; ++j){
//					t = -1.0;
//					s = basis::tet(log2p).yp(j);
//					r = 0.5*(basis::tet(log2p).xp(i)+1.0)*(-s-t)-1.0;
//					basis::tet(log2p).ptprobe2d(1,&f(m)(i,j),r,s,&uht(0)(0),3+3*basis::tet(log2p).em+basis::tet(log2p).fm);
//				}
//			}
//		}
//
//		ugtouht(td,0);		  
//		cout << "face 0" << endl;	
//		for(int i = 0; i < MXGP; ++i){
//			for(int j = 0; j < MXGP; ++j){
//				t = -1.0;
//				s = 0.5*(basis::tet(log2p).yp(j)+1.0)*(1.0-t)-1.0;
//				r = 0.5*(basis::tet(log2p).xp(i)+1.0)*(-s-t)-1.0;
//				basis::tet(log2p).ptprobe(1,&g(0)(i,j),r,s,t,&uht(0)(0),basis::tet(log2p).tm);
//				cout << f(0)(i,j)-g(0)(i,j) << endl;
//			}
//			//cout << basis::tet(log2p).xp(i) <<  basis::tet(log2p).yp(i) << basis::tet(log2p).zp(i) << endl;
//		}
//		
//		cout << "face 1" << endl;	
//		for(int i = 0; i < MXGP; ++i){
//			for(int j = 0; j < MXGP; ++j){
//				t = -1.0;
//				s = 0.5*(basis::tet(log2p).yp(j)+1.0)*(1.0-t)-1.0;
//				r = 0.5*(basis::tet(log2p).xp(i)+1.0)*(-s-t)-1.0;
//				basis::tet(log2p).ptprobe(1,&g(1)(i,j),r,t,s,&uht(0)(0),basis::tet(log2p).tm);
//				cout << f(1)(i,j)-g(1)(i,j) << endl;
//			}
//		}
//		
//		cout << "face 3" << endl;	
//		for(int i = 0; i < MXGP; ++i){
//			for(int j = 0; j < MXGP; ++j){
//				t = -1.0;
//				s = 0.5*(basis::tet(log2p).yp(j)+1.0)*(1.0-t)-1.0;
//				r = 0.5*(basis::tet(log2p).xp(i)+1.0)*(-s-t)-1.0;
//				basis::tet(log2p).ptprobe(1,&g(3)(i,j),t,r,s,&uht(0)(0),basis::tet(log2p).tm);
//				cout << f(3)(i,j)-g(3)(i,j) << endl;
//			}
//		}
		
		//cout << tet(0).sgn << endl<< tri(0).sgn << endl<< tri(1).sgn<< endl<< tri(2).sgn<< endl<<tri(3).sgn<< endl;
//		for(int i = 0; i < ntet; ++i)
//			cout << tet(i).sgn << endl;
//		for(int i = 0; i < ntri; ++i)
//			cout << tri(i).sgn << tri(i).seg << endl;

	return;
	
}