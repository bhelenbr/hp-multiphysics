//
//  hp_getnewbdry.cpp
//  tet_hp
//
//  Created by Brian Helenbrook on 7/2/15.
//
//

//
//  hp_getnewbdry.cpp
//  tri_hp
//
//  Created by Brian Helenbrook on 3/8/15.
//
//

#include <stdio.h>
#include "hp_boundary.h"
#include <tet_boundary.h>

class tet_hp_vtype {
public:
	static const int ntypes = 1;
	enum ids {unknown=-1,plain};
	const static char names[ntypes][40];
	static int getid(const char *nin) {
		for(int i=0;i<ntypes;++i)
			if (!strcmp(nin,names[i])) return(i);
		return(unknown);
	}
};

const char tet_hp_vtype::names[ntypes][40] = {"plain"};

hp_vrtx_bdry* tet_hp::getnewvrtxobject(int bnum, std::string name) {
	std::string keyword,val;
	std::istringstream data;
	int type;
	hp_vrtx_bdry *temp;
	
	type = tet_hp_vtype::getid(name.c_str());
	if (type == tet_hp_vtype::unknown) {
		*gbl->log << "unknown vrtx type:" << keyword << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
	switch(type) {
		default: {
			temp = new hp_vrtx_bdry(*this,*vbdry(bnum));
		}
	}
	return(temp);
}

class tet_hp_etype {
public:
	static const int ntypes = 4;
	enum ids {unknown=-1,plain};
	static const char names[ntypes][40];
	static int getid(const char *nin) {
		for(int i=0;i<ntypes;++i)
			if (!strcmp(nin,names[i])) return(i);
		return(-1);
	}
};

const char tet_hp_etype::names[ntypes][40] = {"plain"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_edge_bdry* tet_hp::getnewedgeobject(int bnum, std::string name) {
	std::string keyword,val;
	std::istringstream data;
	int type;
	hp_edge_bdry *temp;
	
	
	type = tet_hp_etype::getid(name.c_str());
	if (type == tet_hp_etype::unknown) {
		*gbl->log << "unknown side type:" << keyword << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
	
	switch(type) {
		default: {
			temp = new hp_edge_bdry(*this,*ebdry(bnum));
		}
	}
	
	return(temp);
}


class tet_hp_ftype {
public:
	static const int ntypes = 4;
	enum ids {unknown=-1,plain};
	static const char names[ntypes][40];
	static int getid(const char *nin) {
		for(int i=0;i<ntypes;++i)
			if (!strcmp(nin,names[i])) return(i);
		return(-1);
	}
};

const char tet_hp_ftype::names[ntypes][40] = {"plain"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_face_bdry* tet_hp::getnewfaceobject(int bnum, std::string name) {
	std::string keyword,val;
	std::istringstream data;
	int type;
	hp_face_bdry *temp;
	
	
	type = tet_hp_ftype::getid(name.c_str());
	if (type == tet_hp_ftype::unknown) {
		*gbl->log << "unknown side type:" << keyword << std::endl;
		sim::abort(__LINE__,__FILE__,gbl->log);
	}
	
	switch(type) {
		default: {
			temp = new hp_face_bdry(*this,*fbdry(bnum));
		}
	}
	
	return(temp);
}

