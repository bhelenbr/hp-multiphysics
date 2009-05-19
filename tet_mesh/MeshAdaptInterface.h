//-----------------------------------------------------------------------------
// -*- C++ -*-

#ifndef __MESHADAPTINTERFACE_H
#define __MESHADAPTINTERFACE_H


#include "MeshDataBaseInterface.h"
#include "AdaptInterface.h"
#include "PWLinearSField.h"
#include "IsoMeshSize.h"
#include "dolfin.h"

namespace dolfin
{

 class Mesh;


 //-----------------------------------------------------------------------------
 class MeshAdaptInterface {

 public:

   MeshAdaptInterface(Mesh *);
   ~MeshAdaptInterface();

 protected:

   void adaptMesh();

 private:

   void deleteIdMappings();
   void getBoundaryVertices();

 protected: // should become protected

   virtual void updateSizeField() = 0;

   // allocate deallocate solver datas
   virtual void deallocateData() = 0;
   virtual void allocateAndComputeData() = 0;

   // build/destroy meshes
   void importMDBMesh();
   void exportToMDBMesh();
   void deleteDMesh();
   void deleteMDBMesh();

   // move nodes
   void prescribeNodesMotion() {};
   void updateMDBCoordinates() {};

   // constrain entities not to be adapted
   void constrainExternalBoundaries();
   void constrainInternalBoundaries();

   // transfer solution
   void attachSolutionToMDBMesh();
   void getSolutionFromMDBMesh();
   void addFunction(string name, Function** f);
   void clearFunctions();
   virtual void deallocateSolution() const = 0;
   virtual void allocateSolution() const = 0;
   virtual void projectSolution() const = 0;

 public:

   Mesh * dMesh;
   Mesh * dMeshOld;

   std::map<string,Function**> functions;

 protected:
   pGModel model;
   pMesh MDBMesh;

   MeshAdapter * adapter;
   PWLSField * sizeField;

   int * locToMDBIds;
   std::map<int,int> MDBToLocIds;

   // cell id's correspondance between MDB mesh (id=order in iterator) and Dolfin mesh
   int * MDBToLocCellIds;

   // Dolfin order
   std::vector<bool> boundaryVertices;
   std::vector<bool> constrainedVertices;
 };
 //-----------------------------------------------------------------------------

}


#endif