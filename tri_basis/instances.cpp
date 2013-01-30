/*
 *  instances.cpp
 *  tri_basis
 *
 *  Created by Brian Helenbrook on 7/15/09.
 *  Copyright 2009 Clarkson University. All rights reserved.
 *
 */
 
#include "tri_basis.h"
#include "intgrt.cpp"
#include "initialize.cpp"
#include "probe.cpp"
#include "intgrt1d.cpp"
#include "proj.cpp"
#include "proj1d.cpp"
#include "ptvalues.cpp"
 
template class tri_basis<1,0>;
template class tri_basis<2,0>;
template class tri_basis<4,0>;
//template class tri_basis<8,0>;

template class tri_basis<1,1>;
template class tri_basis<2,1>;
template class tri_basis<4,1>;
//template class tri_basis<8,1>;

