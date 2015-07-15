/*
 *  test.cpp
 *  mesh3d
 *
 *  Created by Michael Brazell on 8/2/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "tet_mesh.h"
#include <utilities.h>
#include <cstdlib>
#include <cstring>
#include <input_map.h>
#include <iostream>

void tet_mesh::test() {
//    int n;

//    const std::string filename="box2";
//    tet_mesh::filetype in = static_cast<tet_mesh::filetype>(baker);
//    input_map input;
//    
//    input["growth factor"] = "1.0";
//    input["filetype"] = "1";
//    input["mesh"] = "box2";
//    
//    tet_mesh::init(input,0);
		
	//in = static_cast<tet_mesh::filetype>(tecplot);
	// output tecplot
	//tet_mesh::output(filename, in);
	
	
//    // output all fields to screen    
//    n = 4;
//    
//    cout << endl << "vertex = " << npnt << endl << endl;
//
//    for(int i = 0; i < n; ++i)
//        cout << pnts(i)(0) << ' ' << pnts(i)(1) << ' ' << pnts(i)(2) << "  tet = " << pnt(i).tet << "  nnbor = " << pnt(i).nnbor <<  endl;
//        
//    cout << endl<< endl << "seg = "<< nseg << endl << endl;    
//    for(int i = 0; i < n; ++i) {
//        for(int j = 0; j < 2; ++j)
//            cout << seg(i).pnt(j) << ' ' ;
//        cout << "  tet = " << seg(i).tet << "  nnbor = " << seg(i).nnbor<< endl;
//    }
//    
//    cout << endl <<endl<< "tri = "<< ntri << endl<< endl;
//    
//    cout << "vertex data" << endl;
//    for(int i = 0; i < n; ++i) {
//        for(int j = 0; j < 3; ++j)
//            cout << tri(i).pnt(j) << ' ' ;
//        cout << endl;
//    }
//    cout << "seg data " << endl;
//    for(int i = 0; i < n; ++i) {
//        for(int j = 0; j < 3; ++j)
//            cout << tri(i).seg(j) << ' ' ;
//        cout << endl;
//    }
//    cout << "tet data" << endl;
//    for(int i = 0; i < n; ++i) {
//        for(int j = 0; j < 2; ++j)
//            cout << tri(i).tet(j) << ' ' ;
//        cout << endl;
//    }
//    cout << "sgn data" << endl;
//    for(int i = 0; i < n; ++i) {
//        for(int j = 0; j < 3; ++j)
//            cout << tri(i).sgn(j) << ' ' ;
//        cout << endl;
//    }
//
//    cout << endl <<endl<< "tet's = "<< ntet << endl<< endl;
//    
//    cout << "vertex data" << endl;
//    for(int i = 0; i < n; ++i) {
//        for(int j = 0; j < 4; ++j)
//            cout << tet(i).pnt(j) << ' ' ;
//        cout << endl;
//    }
//    
//
//    
//    cout << "seg data" << endl;
//    for(int i = 0; i < n; ++i) {
//        for(int j = 0; j < 6; ++j)
//            cout << tet(i).seg(j) << ' ' ;
//        cout << endl;
//    }
//    cout << "face data" << endl;
//    for(int i = 0; i < n; ++i) {
//        for(int j = 0; j < 4; ++j)
//            cout << tet(i).tri(j) << ' ' ;
//        cout << endl;
//    }
//
//    cout << "tet data" << endl;
//    for(int i = 0; i < n; ++i) {
//        for(int j = 0; j < 4; ++j)
//            cout << tet(i).tet(j) << ' ' ;
//        cout << endl;
//    }
//    cout << "sgn data" << endl;
//    for(int i = 0; i < n; ++i) {
//        for(int j = 0; j < 6; ++j)
//            cout << tet(i).sgn(j) << ' ' ;
//        cout << endl;
//    }    
//    
//    cout << endl << "face bdry data" << endl;
//    for(int i = 0; i < n; ++i) {
//        for(int j = 0; j < 3; ++j)
//            cout << fbdry(0)->tri(i).pnt(j) << ' ' ;
//        cout << endl;
//    }
//    
//    cout << endl << "max nnbor pnt = " ;
//    int maxnnbor = 0;
//    int vfind;
//    for( int i = 0 ; i < npnt; ++i) {
//        if (pnt(i).nnbor > maxnnbor) {
//            maxnnbor = pnt(i).nnbor;
//             vfind = i;
//        }
//    }
//        
//    cout << maxnnbor << ' ' << vfind;
//    
//    cout << endl << "face bdry data re-indexed" << endl;
//    for(int i = 0; i < n; ++i) {
//        for(int j = 0; j < 3; ++j)
//            cout << fbdry(0)->tri(i).pnt(j) << ' ' ;
//        cout << endl;
//    }
//    
//    cout << endl << "global index data" << endl;
//    for(int i = 0; i < n; ++i)        
//        cout << fbdry(0)->pnt(i).gindx << endl;
//    
//    
//    cout << endl << "bdry side data" << endl;
//    for(int i = 0; i < n; ++i) {
//        for(int j = 0; j < 2; ++j)
//            cout << fbdry(0)->seg(i).pnt(j) << ' ' ;
//        cout << endl;
//    }
//    
//    cout << endl << "bdry side tri data" << endl;
//    for(int i = 0; i < n; ++i) {
//        for(int j = 0; j < 2; ++j)
//            cout << fbdry(0)->seg(i).tri(j) << ' ' ;
//        cout << endl;
//    }
//    
//    cout << endl << "bdry tri side data" << endl;
//    for(int i = 0; i < n; ++i) {
//        for(int j = 0; j < 3; ++j)
//            cout << fbdry(0)->tri(i).seg(j) << ' ' ;
//        cout << endl;
//    }
//    
//    cout << endl << "bdry tri tri data" << endl;
//    for(int i = 0; i < n; ++i) {
//        for(int j = 0; j < 3; ++j)
//            cout << fbdry(0)->tri(i).tri(j) << ' ' ;
//        cout << endl;
//    }
//    
//    int vind = 0;
//    
//    cout << endl << "tet's surrounding gbl vertex: " << vind << endl;
//    tet_mesh::vertexball(vind);
//    for(int i = 0; i < pnt(vind).nnbor; ++i)
//        cout << gbl->i2wk(i) << endl;
//        
//    int eind = 0;
//        
//    cout << endl << "tet's surrounding gbl edge: " << eind << " total tet's: " << seg(eind).nnbor << endl;
//    tet_mesh::ring(eind);
//    for(int i = 0; i < seg(eind).nnbor; ++i)
//        cout << gbl->i2wk(i) << endl;
		
	
//    cout << endl << "tri global index data" << endl;
//    for(int i = 0; i < n; ++i)        
//        cout << fbdry(0)->tet(i).gindx << endl;
//        
//        
//    cout << endl << "side global index data" << endl;
//    for(int i = 0; i < n; ++i)        
//        cout << fbdry(0)->seg(i).gindx << endl;
//    
//    cout << endl << "tri's surrounding lcl vertex: " << vind << endl;
//    fbdry(0)->vertexcircle(vind);
//    for(int i = 0; i < fbdry(0)->pnt(vind).nnbor; ++i)
//        cout << i2wk(i) << endl;
		
		
	////output mesh
//    in = static_cast<tet_mesh::filetype>(tecplot);
//    tet_mesh::output(filename, in);


//for each face:
//The segment listed has the same vertices as the vertices of that edge of the face
//The sign is compatible with the definition of the segment
//
//for each tet:
//Definition of seg's are compatible (by checking vertex numbers)
//Definition of signs are compatible ((by checking orientations)
//Definition of rot is combatible and that the vertex definitions in the face match those in the tet.
//Neighboring tet info in face is compatible with tet
//The neigboring tet info in tet is compatible with neighboring tet info in face

//    TinyVector<int,3> v3;
//    TinyVector<int,2> v2;
//    TinyVector<int,4> v4;
//    
//    TinyMatrix<int,3,2> vf;
//    vf(0,0) = 1, vf(0,1) = 2;
//    vf(1,0) = 2, vf(1,1) = 0;
//    vf(2,0) = 1, vf(2,1) = 0;
//    
//    TinyMatrix<int,6,2> vt;
//    vt(0,0) = 1, vt(0,1) = 2;
//    vt(1,0) = 2, vt(1,1) = 0;
//    vt(2,0) = 1, vt(2,1) = 0;
//    vt(3,0) = 1, vt(3,1) = 0;
//
//    int sign;
//    
//    for(int i = 0; i < ntri; ++i) {
//        v3=tri(i).pnt;
//        for(int j = 0; j < 3; ++j) {
//            sign = tri(i).sgn(j);
//            v2 = seg(tri(i).seg(j)).pnt;
//            if (sign == 1) {
//                if (v3(vf(j,0)) != v2(0) || v3(vf(j,1)) != v2(1))
//                    cout << "bad segment on tri: " << i << endl;
//            }
//            else{
//                if (v3(vf(j,0)) != v2(1) || v3(vf(j,1)) != v2(0))
//                    cout << "bad segment on tri: " << i << endl;
//            }        
//        }
//    } 

//    int tet1,tet2,fac1,fac2;
//    
//    for(find = 0; find < ntri; ++find) {
//        tet1=tri(find).tet(0);
//        tet2=tri(find).tet(1);
//        if (tet2 > -1) {
//            for(int i = 0; i < 4; ++i) {
//                if (tet(tet1).tri(i)=find)
//                    fac1=i;
//                if (tet(tet2).tri(i)=find)
//                    fac2=i;                                    
//            }
//            
//        }
//        
//    }
	checkintegrity();
return;
}

void tet_mesh::checkintegrity() {
	int tind,i,j,sind,find;
	// vertices of each segment on tet
	TinyMatrix<int,6,2> vs;    
	vs(0,0)=2, vs(0,1)=3; 
	vs(1,0)=3, vs(1,1)=1; 
	vs(2,0)=2, vs(2,1)=1;
	vs(3,0)=1, vs(3,1)=0;
	vs(4,0)=2, vs(4,1)=0; 
	vs(5,0)=3, vs(5,1)=0;
	
	// vertices of each face on a tet
	TinyMatrix<int,4,3> vf1;
	vf1(0,0)=1, vf1(0,2)=2, vf1(0,1)=3;
	vf1(1,0)=0, vf1(1,2)=3, vf1(1,1)=2;
	vf1(2,0)=0, vf1(2,2)=1, vf1(2,1)=3;
	vf1(3,0)=0, vf1(3,2)=2, vf1(3,1)=1;
	TinyMatrix<int,4,3> vf2;
	vf2(0,0)=1, vf2(0,1)=2, vf2(0,2)=3;
	vf2(1,0)=0, vf2(1,1)=3, vf2(1,2)=2;
	vf2(2,0)=0, vf2(2,1)=1, vf2(2,2)=3;
	vf2(3,0)=0, vf2(3,1)=2, vf2(3,2)=1; 
	
	// segments of each face on a tet
	TinyMatrix<int,4,3> sf1;
	sf1(0,0)=0, sf1(0,2)=1, sf1(0,1)=2;
	sf1(1,0)=0, sf1(1,2)=4, sf1(1,1)=5;
	sf1(2,0)=1, sf1(2,2)=5, sf1(2,1)=3;
	sf1(3,0)=2, sf1(3,2)=3, sf1(3,1)=4;
	TinyMatrix<int,4,3> sf2;
	sf2(0,0)=0, sf2(0,1)=1, sf2(0,2)=2;
	sf2(1,0)=0, sf2(1,1)=4, sf2(1,2)=5;
	sf2(2,0)=1, sf2(2,1)=5, sf2(2,2)=3;
	sf2(3,0)=2, sf2(3,1)=3, sf2(3,2)=4;
	
	// vertices of each segment on a tri
	TinyMatrix<int,3,2> vst;
	vst(0,0)=1, vst(0,1)=2;
	vst(1,0)=2, vst(1,1)=0;
	vst(2,0)=0, vst(2,1)=1;
	
	
	for(i=0;i<maxvst;++i) {
		if (gbl->i1wk(i) != -1) {
			*gbl->log << "gbl->i1wk check failed" << std::endl;
			exit(1);
		}
	}
		
	for(tind=0;tind<ntet;++tind) {
		if (tet(tind).vol < 0) {
			*gbl->log << "negative tet volume on tet: " << tind << endl;
			exit(1);
		}
		
		for(i = 0; i < 6; ++i) {
			sind = tet(tind).seg(i);
			if (tet(tind).sgn(i) == 1) {
				if (seg(sind).pnt(0) != tet(tind).pnt(vs(i,0)) && seg(sind).pnt(1) != tet(tind).pnt(vs(i,1))) {
					*gbl->log << "tet seg vertex error type 1"<< endl;
					exit(1);
				}
			}
			else{
				if (seg(sind).pnt(1) != tet(tind).pnt(vs(i,0)) && seg(sind).pnt(0) != tet(tind).pnt(vs(i,1))) {
					*gbl->log << "tet seg vertex error type 2 " << sind << ' ' << tind << ' ' << i << endl;
					exit(1);
				}
			}
		}
		
		/* Test Tri Data */
		for(i = 0; i < 4; ++i) {
			find = tet(tind).tri(i);
			if (tri(find).tet(1) > -1) {
				if (tet(tind).tet(i) < 0) {
					*gbl->log << "tet neighbor data error type 1"<< endl;
					exit(1);
				}
				if (tri(find).tet(1) != tet(tind).tet(i) && tri(find).tet(0) != tet(tind).tet(i)) {
					*gbl->log << "tet neighbor data error type 2"<< endl;
					exit(1);
				}
			}

			if (tet(tind).rot(i) == 1) {
				for(j = 0; j < 3; ++j) {
					if (tet(tind).pnt(vf1(i,j)) != tri(find).pnt(j)) {
						*gbl->log << "face vertex error type 1"<< endl;
						exit(1);
					}
					if (tet(tind).seg(sf1(i,j)) != tri(find).seg(j)) {
						*gbl->log << "face seg error type 1"<< endl;
						exit(1);
					}
				}
				if (tind != tri(find).tet(0)) {
					*gbl->log << "tet rot error type 1"<< endl;
					exit(1);
				}
			}
			else{
				for(j = 0; j < 3; ++j) {
					if (tet(tind).pnt(vf2(i,j)) != tri(find).pnt(j)) {
						*gbl->log << "face vertex error type 2"<< endl;
						exit(1);
					}
					if (tet(tind).seg(sf2(i,j)) != tri(find).seg(j)) {
						*gbl->log << "face seg error type 2"<< endl;
						exit(1);
					}
				}
//                if (tri(find).tet(1) != -1) {
//                    if (tet(tind).tet(i) != tri(find).tet(1)) {
//                        *gbl->log << "tet rot error type 2 "<< tri(find).tet(1)<<' '<<tet(tind).tet(i)<< ' '<< tri(find).tet(1)<< ' ' << tet(tind).rot(i) << endl;                        
//                        exit(1);
//                    }
//                }
			}
							
			for(j = 0; j < 3; ++j) {
				sind = tri(find).seg(j);
				if (tri(find).sgn(j) == 1) {
					if (tri(find).pnt(vst(j,0)) != seg(sind).pnt(0) && tri(find).pnt(vst(j,1)) != seg(sind).pnt(1)) {
						*gbl->log << "seg vertex error type 1 "<< find << ' ' << j << ' ' << sind << endl;
						output("error",easymesh);
						output("error",grid);
						exit(1);
					}
				}
				else{
					if (tri(find).pnt(vst(j,1)) != seg(sind).pnt(0) && tri(find).pnt(vst(j,0)) != seg(sind).pnt(1)) {
						*gbl->log << "seg vertex error type 2 "<< find << ' ' << j << ' ' << sind << endl;
						output("error",easymesh);
						output("error",grid);
						exit(1);
					}
				}
			}
		}   
	}
	
	for(find=0;find<ntri;++find) {
		for(j = 0; j < 3; ++j) {
			sind = tri(find).seg(j);
			if (tri(find).sgn(j) == 1) {
				if (tri(find).pnt(vst(j,0)) != seg(sind).pnt(0) && tri(find).pnt(vst(j,1)) != seg(sind).pnt(1)) {
					*gbl->log << "seg vertex error type 3 "<< find << ' ' << j << ' ' << sind << endl;
					output("error",grid);
					output("error",easymesh);
					exit(1);
				}
			}
			else{
				if (tri(find).pnt(vst(j,1)) != seg(sind).pnt(0) && tri(find).pnt(vst(j,0)) != seg(sind).pnt(1)) {
					*gbl->log << "seg vertex error type 4 "<< find << ' ' << j << ' ' << sind << endl;
					output("error",grid);
					output("error",easymesh);
					exit(1);
				}
			}
		}
	}
	
	
	for(int i=0;i<nfbd;++i)
		fbdry(i)->checkintegrity();
	
	for(int i=0;i<nebd;++i)
		ebdry(i)->checkintegrity();
	
	for(int i=0;i<nvbd;++i)
		vbdry(i)->checkintegrity();
	
	return;
}

