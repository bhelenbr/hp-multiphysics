// -*- C++ -*-
// -------------------------------------------------------------------
// MAdLib - Copyright (C) 2008-2009 Universite catholique de Louvain
//
// See the Copyright.txt and License.txt files for license information. 
// You should have received a copy of these files along with MAdLib. 
// If not, see <http://www.madlib.be/license/>
//
// Please report all bugs and problems to <contrib@madlib.be>
// -------------------------------------------------------------------
// Author: Gaetan Compere
//
// This file provides an example of an interface to MAdLib as it 
// could be implemented in a physical solver requiring mesh adaptivity.
// -------------------------------------------------------------------

#ifndef __MADLIBINTERFACE_H
#define __MADLIBINTERFACE_H

#include <MAdLib/ModelInterface.h>
#include <MAdLib/MeshDataBaseInterface.h>
#include <MAdLib/MeshDataBase.h>
#include <MAdLib/MAdModel.h>
#include <MAdLib/ModelInterface.h>
#include <MAdLib/AdaptInterface.h>
#include <MAdLib/PWLinearSField.h>
#include <input_map.h>
#include <utility>
#include <set>

#ifdef SINGLE
#define FLT float
#define EPSILON FLT_EPSILON
#else
#ifndef FLT
#define FLT double
#define EPSILON DBL_EPSILON
#endif
#endif

class tet_mesh;
//-----------------------------------------------------------------------------
// Class interfacing MAdLib with 'Solver'
//-----------------------------------------------------------------------------

namespace MAdLibInterface {

	// Input gmsh format using MAdLib
	// void input(tet_mesh *, const std::string filename, FLT grwfac, input_map& input);

	// Mesh to mesh conversion
	void importFromMAdMesh(const MAd::pMesh, tet_mesh *);
	void exportToMAdMesh(const tet_mesh *, MAd::pMesh);
	void exportToMAdModel(const tet_mesh *, MAd::pGModel);

	// Size field construction
	void buildSizeField(MAd::PWLSField *);
	
	// To coarsen a tet_mesh
	void coarsenMesh(FLT factor, tet_mesh *);

};

//-----------------------------------------------------------------------------

#endif

