/*
 *  box.cpp
 *  mblock
 *
 *  Created by helenbrk on Wed Aug 08 2001.
 *  Copyright (c) 2001 __CompanyName__. All rights reserved.
 *
 */
#include"quad_impl.h"

/* STATIC VARIABLES FOR SEARCHING */
class box<1> **quadtree<1>::srchlst = 0;
int quadtree<1>::maxsrch = 0;

class box<2> **quadtree<2>::srchlst = 0;
int quadtree<2>::maxsrch = 0;
 
class box<3> **quadtree<3>::srchlst = 0;
int quadtree<3>::maxsrch = 0;

template class box<1>;
template class box<2>;
template class box<3>;

template class quadtree<1>;
template class quadtree<2>;
template class quadtree<3>;

