#include "MeshAdaptInterface.h"
#include "NullSField.h"
#include "ModelInterface.h"
//#include "modeler.h"
#include "NullModel.h"
#include "dolfin/DiscreteFunction.h"

#include <map>

#include <boost/timer.hpp>

using namespace dolfin;

//-----------------------------------------------------------------------------
void DolfinCBFunction (pPList before, pPList after, void *data ,
                      operationType type , pEntity ppp) {

 MeshAdaptInterface * ai = static_cast<MeshAdaptInterface *>(data);

 std::map<string,Function**>::const_iterator fIter = ai->functions.begin();
 std::map<string,Function**>::const_iterator fLast = ai->functions.end();
 for (; fIter != fLast; fIter++) {

   string fName = (*fIter).first;
   pMeshDataId dataId = MD_lookupMeshDataId(fName.c_str());

   Function * function = *((*fIter).second);

//     Mesh& dMesh = function->mesh();

//     int N = dMesh.numVertices();
//     int d = dMesh.topology().dim();

//     std::cout << "CB " << std::endl;
//     std::cout << "N: " << N << std::endl;
//     std::cout << "d: " << d << std::endl;


       switch (type) {
       case LMM_ESPLIT: // assumes the edge is split in the middle
         {
	   if ( fName[0] == 'c' )
           {
               // find the old edge
               void *tmp=0;
               pEntity pE = (pEntity) PList_next(before,&tmp);

               // get coordinates and datas at old nodes
               double data0 = 0.;
               pVertex pV0 = E_vertex((pEdge)pE, 0);
               int gotit0 = EN_getDataDbl((pEntity)pV0, dataId,  &data0);

               double data1 = 0.;
               pVertex pV1 = E_vertex((pEdge)pE, 1);
               int gotit1 = EN_getDataDbl((pEntity)pV1, dataId,  &data1);

               if (!gotit0 || !gotit1) {
                 printf("Error: one of the nodes has no data attached to with name %s",
                        fName.c_str());
                 throw;
               }

               // Interpolate the data at the new node. 
               // Warning: it assumes the point is located at the middle of the edge
               double newData = 0.5 * ( data0 + data1 );

               // attach it
               EN_attachDataDbl(ppp, dataId, newData);
             }
	   else if ( fName[0] == 'd' )
	   {
	     //std::cout << "refining d marker" << std::endl;
	   }
         }
         break;
       case LMM_ECOLLAPSE:
         {
           //EN_deleteData(ppp, dataId);
         }
         break;
       case LMM_FSWAP:
         {
         }
         break;
       case LMM_ESWAP:
         {
         }
         break;
       default:
         printf("Error in callback function: this type of operation is not handled: %d",type);
         throw;
       }

 }

};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
MeshAdaptInterface::MeshAdaptInterface(Mesh * _dMesh):
 locToMDBIds(NULL), MDBToLocCellIds(NULL), dMeshOld(NULL),
 model(0), MDBMesh(0), sizeField(0), adapter(0)
{
 dMesh = _dMesh;

 //pGModel model = 0;

//   MDBMesh = MS_newMesh(model);
//   exportToMDBMesh();

//   sizeField = new PWLSField(MDBMesh);
//   sizeField->setCurrentSize();

//   adapter = new MeshAdapter(MDBMesh,sizeField);

//   adapter->addCallback(DolfinCBFunction,(void*)this);

 ////M_writeSMS(MDBMesh,"init.msh",2,NULL);
 ////MDBMesh->shrink();
}

//-----------------------------------------------------------------------------
MeshAdaptInterface::~MeshAdaptInterface()
{
 if (locToMDBIds)        delete [] locToMDBIds;
 if (MDBToLocCellIds) delete [] MDBToLocCellIds;

 if (adapter)   delete adapter;
 if (sizeField) delete sizeField;
 if (MDBMesh)   delete MDBMesh;
 if (model)     delete model;

 if (dMeshOld)  { delete dMeshOld; dMeshOld=NULL; }
}

//-----------------------------------------------------------------------------
void MeshAdaptInterface::getBoundaryVertices()
{
 boundaryVertices.resize(dMesh->numVertices());
 for(int i = 0; i < dMesh->numVertices(); i++)
 {
   boundaryVertices[i] = false;
 }

 MeshFunction<unsigned int> bnd_vertex_map;
 MeshFunction<unsigned int> bnd_cell_map;

 BoundaryMesh boundary(*dMesh,bnd_vertex_map,bnd_cell_map);

 for (VertexIterator v(boundary); !v.end(); ++v)
 {
   int bnd_vi = v->index();
   int vol_vi = bnd_vertex_map(*v);

   boundaryVertices[vol_vi] = true;

   Vertex vv(*dMesh, vol_vi);
   //cout << vv.point() << endl;
 }

}
//-----------------------------------------------------------------------------
void MeshAdaptInterface::importMDBMesh()
{
//   dMeshOld = dMesh;

 //dMesh = new Mesh;

 deleteIdMappings();

 int dim = getDim(MDBMesh);
 int numVertices = M_numVertices(MDBMesh);
 int numCells = 0;

 CellType * cell_type;
 switch (dim) {
 case 3: 
   cell_type = CellType::create(CellType::tetrahedron);
   numCells = M_numRegions(MDBMesh);
   break;
 case 2:
   cell_type = CellType::create(CellType::triangle);
   numCells = M_numFaces(MDBMesh);
   break;
 case 1:
   cell_type = CellType::create(CellType::interval);
   numCells = M_numEdges(MDBMesh);
   break;
 case 0:
   cell_type = CellType::create(CellType::point);
   numCells = M_numVertices(MDBMesh);
   break;
 }

 // --- create the mesh ---
 MeshEditor editor;
 editor.open(*dMesh, cell_type->cellType(),
	     dim, dim);

 editor.initVertices(numVertices);
 editor.initCells(numCells);


 // --- add vertices ---
 int current_vertex = 0;
 if ( locToMDBIds )
 {
   delete [] locToMDBIds;
   locToMDBIds = 0;
 }
 locToMDBIds = new int[numVertices];
 VIter vit = M_vertexIter(MDBMesh);
 while (pVertex pv = VIter_next(vit))
   {
     double xyz[3];
     V_coord(pv,xyz);
     Point p(xyz[0],xyz[1],xyz[2]);
     editor.addVertex(current_vertex, p);
     pPoint pp = V_point(pv);
     locToMDBIds[current_vertex] = P_id(pp);
     current_vertex++;
   }
 VIter_delete(vit);

 // --- create reverse vertex id mapping ---
 for (int i=0; i<numVertices; i++) {
   MDBToLocIds[locToMDBIds[i]] = i;
 }

 // --- add cells ---
 int current_cell = 0;
 if ( MDBToLocCellIds )
 {
   delete [] MDBToLocCellIds;
   MDBToLocCellIds = 0;
 }
 MDBToLocCellIds = new int[numCells];
 if (dim==3) {
   RIter rit = M_regionIter(MDBMesh);
   while (pRegion pr = RIter_next(rit))
     {
       Array<uint> cell_vertices(cell_type->numEntities(0));

       for (int iV=0; iV<R_numVertices(pr); iV++) {
         pVertex pv = R_vertex(pr, iV);
         pPoint pp = V_point(pv);
         cell_vertices[iV] = MDBToLocIds[P_id(pp)];
       }

       editor.addCell(current_cell, cell_vertices);
       MDBToLocCellIds[current_cell] = current_cell; // if dolfin numbers cells in order of creation -> ?
       current_cell++;
     }
   RIter_delete(rit);
 }
 else  if (dim==2) {
   FIter fit = M_faceIter(MDBMesh);
   while (pFace pf = FIter_next(fit))
     {
       Array<uint> cell_vertices(cell_type->numEntities(0));

       for (int iV=0; iV<F_numVertices(pf); iV++) {
         pVertex pv = F_vertex(pf, iV);
         pPoint pp = V_point(pv);
         cell_vertices[iV] = MDBToLocIds[P_id(pp)];
       }

       editor.addCell(current_cell, cell_vertices);
       MDBToLocCellIds[current_cell] = current_cell; // if dolfin numbers cells in order of creation -> ?
       current_cell++;
     }
   FIter_delete(fit);
 }
 else  if (dim==1) {
   EIter eit = M_edgeIter(MDBMesh);
   while (pEdge pe = EIter_next(eit))
     {
       Array<uint> cell_vertices(cell_type->numEntities(0));

       for (int iV=0; iV<2; iV++) {
         pVertex pv = E_vertex(pe, iV);
         pPoint pp = V_point(pv);
         cell_vertices[iV] = MDBToLocIds[P_id(pp)];
       }

       editor.addCell(current_cell, cell_vertices);
       MDBToLocCellIds[current_cell] = current_cell; // if dolfin numbers cells in order of creation -> ?
       current_cell++;
     }
   EIter_delete(eit);
 }
 else  {
   VIter vit = M_vertexIter(MDBMesh);
   while (pVertex pv = VIter_next(vit))
     {
       Array<uint> cell_vertices(cell_type->numEntities(0));

       pPoint pp = V_point(pv);
       cell_vertices[0] = MDBToLocIds[P_id(pp)];

       editor.addCell(current_cell, cell_vertices);
       MDBToLocCellIds[current_cell] = current_cell; // if dolfin numbers cells in order of creation -> ?
       current_cell++;
     }
   VIter_delete(vit);
 }

 editor.close();

}

//-----------------------------------------------------------------------------
void MeshAdaptInterface::exportToMDBMesh()
{
 if (!MDBMesh->model)  MDBMesh->model = new NullModel; 

 std::cout << "Num pts: "<<dMesh->numVertices()<<", num cells: "<<dMesh->numCells()<<std::endl;

 deleteIdMappings();

 uint* cells = dMesh->cells();
 int numCellV = dMesh->type().numEntities(0);
 GEntity *geom = 0;

 if ( locToMDBIds )
 {
   delete [] locToMDBIds;
   locToMDBIds = 0;
 }
 locToMDBIds = new int[dMesh->numVertices()];

 if ( MDBToLocCellIds )
 {
   delete [] MDBToLocCellIds;
   MDBToLocCellIds = 0;
 }
 MDBToLocCellIds = new int[dMesh->numCells()];

 switch(dMesh->type().cellType()) {
 case CellType::tetrahedron:
   {
     real* xyz = dMesh->coordinates();
     for (int iV=0; iV < dMesh->numVertices(); iV++) {
       MDBMesh->add_point(iV+1,xyz[3*iV],xyz[3*iV+1],xyz[3*iV+2]);
       locToMDBIds[iV] = iV+1;
     }

     geom = MDBMesh->model->regionByTag(1);
     for (int iC=0; iC < dMesh->numCells(); iC++) {
       MDBMesh->add_tet(cells[iC*numCellV]  +1,
                        cells[iC*numCellV+1]+1,
                        cells[iC*numCellV+2]+1,
                        cells[iC*numCellV+3]+1,
                        geom); 
       MDBToLocCellIds[iC] = iC; // ok if the id in dolfin is iC !
     }
   }
   break;
 case CellType::triangle:
   {
     real* xyz = dMesh->coordinates();
     for (int iV=0; iV < dMesh->numVertices(); iV++) {
       MDBMesh->add_point(iV+1,xyz[2*iV],xyz[2*iV+1],0.);
       locToMDBIds[iV] = iV+1;
     }

     geom = MDBMesh->model->faceByTag(1);
     for (int iC=0; iC < dMesh->numCells(); iC++) {
       MDBMesh->add_triangle(cells[iC*numCellV]  +1,
                             cells[iC*numCellV+1]+1,
                             cells[iC*numCellV+2]+1,
                             geom); 
       MDBToLocCellIds[iC] = iC; // ok if the id in dolfin is iC !
     }
   }
   break;
 case CellType::interval:
   {
     real* xyz = dMesh->coordinates();
     for (int iV=0; iV < dMesh->numVertices(); iV++) {
       MDBMesh->add_point(iV+1,xyz[iV],0.,0.);
       locToMDBIds[iV] = iV+1;
     }

     geom = MDBMesh->model->edgeByTag(1);
     for (int iC=0; iC < dMesh->numCells(); iC++) {
       MDBMesh->add_edge(cells[iC*numCellV]  +1,
                         cells[iC*numCellV+1]+1,
                         geom); 
       MDBToLocCellIds[iC] = iC; // ok if the id in dolfin is iC !
     }
   }
   break;
 case CellType::point:
   {
     real* xyz = dMesh->coordinates();
     for (int iV=0; iV < dMesh->numVertices(); iV++) {
       MDBMesh->add_point(iV+1,xyz[iV],0.,0.);
       locToMDBIds[iV] = iV+1;
     }

     geom = MDBMesh->model->vertexByTag(1);
     for (int iC=0; iC < dMesh->numCells(); iC++) {
       pPoint pp = MDBMesh->find_point(cells[iC*numCellV]+1);
       pp->g = geom;
       MDBToLocCellIds[iC] = iC; // ok if the id in dolfin is iC !
     }
   }
   break;
 }

 MDBMesh->classify_unclassified_entities();

 MDBMesh->destroyStandAloneEntities();

 // --- create reverse vertex id mapping ---
 for (int i=0; i<dMesh->numVertices(); i++) {
   MDBToLocIds[locToMDBIds[i]] = i;
 }

 getBoundaryVertices();
}
//-----------------------------------------------------------------------------
void MeshAdaptInterface::deleteIdMappings()
{
 if (locToMDBIds) delete [] locToMDBIds;
 locToMDBIds = 0;
 MDBToLocIds.clear();
}

//-----------------------------------------------------------------------------
void MeshAdaptInterface::deleteMDBMesh()
{
 if (MDBMesh) delete MDBMesh;
 MDBMesh = NULL;
}

//-----------------------------------------------------------------------------
void MeshAdaptInterface::addFunction(string name, Function** fct)
{
//   DiscreteFunction* df = static_cast<DiscreteFunction*>(fct->f);
//   string signature = df->finite_element->signature();

//   functions[name] = std::make_pair(signature,fct);

 functions[name] = fct;
}

//-----------------------------------------------------------------------------
void MeshAdaptInterface::clearFunctions()
{
 functions.clear();
}

//-----------------------------------------------------------------------------
// TODO: extend to solutions on cells
void MeshAdaptInterface::attachSolutionToMDBMesh()
{
 std::map<string,Function**>::const_iterator fIter = functions.begin();
 std::map<string,Function**>::const_iterator fLast = functions.end();

 int N = dMesh->numVertices();
 int M = dMesh->numCells();
 int d = dMesh->topology().dim();

 cout << "Attach " << endl;
 cout << "N: " << N << endl;
 cout << "d: " << d << endl;


 for (; fIter != fLast; fIter++) {

   string fName = (*fIter).first;
   pMeshDataId dataId = MD_lookupMeshDataId(fName.c_str());

   Function * function = *((*fIter).second);
   assert(function);
   real* dataArray = function->vector().vec().array();

   if ( fName[0] == 'c' )
   {
     VIter vit = M_vertexIter(MDBMesh);
     while (pVertex pv = VIter_next(vit))
     {
	// get dolfin point id
	pPoint pp = V_point(pv);
	int mdbId = P_id(pp);
	int dId = MDBToLocIds[mdbId];
	
	real data = 0.0;
	if(fName == "cUx")
	 data = dataArray[0 * N + dId];
	else if(fName == "cUy")
	 data = dataArray[1 * N + dId];
	else if(fName == "cUz")
	 data = dataArray[2 * N + dId];
	else if(fName == "ch0f")
	 data = dataArray[dId];
	else if(fName == "cphi")
	 data = dataArray[dId];
	else
	 dolfin::error("Unknown function");
	
	EN_attachDataDbl((pEntity)pv,dataId,(double)data);
     }
     VIter_delete(vit);
   }
   else if ( fName[0] == 'd' )
   {
     // phi
     // Assume dim == 2
     int mdbId = 0;
     RIter rit = M_regionIter(MDBMesh);
     while (pRegion pr = RIter_next(rit))
     {
	// get dolfin element id
	int dId = MDBToLocCellIds[mdbId];

	real data = 0.0;

	if(fName == "dphi")
	{
	 data = dataArray[dId];
	}
	else if(fName == "dS00")
	{
	 data = dataArray[0 * M + dId];
	}
	else if(fName == "dS01")
	{
	 data = dataArray[1 * M + dId];
	}
	else if(fName == "dS02")
	{
	 data = dataArray[2 * M + dId];
	}
	else if(fName == "dS10")
	{
	 data = dataArray[3 * M + dId];
	}
	else if(fName == "dS11")
	{
	 data = dataArray[4 * M + dId];
	}
	else if(fName == "dS12")
	{
	 data = dataArray[5 * M + dId];
	}
	else if(fName == "dS20")
	{
	 data = dataArray[6 * M + dId];
	}
	else if(fName == "dS21")
	{
	 data = dataArray[7 * M + dId];
	}
	else if(fName == "dS22")
	{
	 data = dataArray[8 * M + dId];
	}
	else
	 dolfin::error("Unknown function");

	EN_attachDataDbl((pEntity)pr,dataId,(double)data);
	
	mdbId++;
     }
     RIter_delete(rit);
   }
   else
   {
     error("unknown function");
   }
 }
}

//-----------------------------------------------------------------------------
// TODO: extend to solutions on cells
void MeshAdaptInterface::getSolutionFromMDBMesh()
{
 std::map<string,Function**>::const_iterator fIter = functions.begin();
 std::map<string,Function**>::const_iterator fLast = functions.end();

 int N = dMesh->numVertices();
 int M = dMesh->numCells();
 int d = dMesh->topology().dim();

 cout << "N: " << N << endl;
 cout << "M: " << M << endl;
 cout << "d: " << d << endl;

 for (; fIter != fLast; fIter++) {

   string fName = (*fIter).first;
   pMeshDataId dataId = MD_lookupMeshDataId(fName.c_str());

   Function * function = *((*fIter).second);

   real* dataArray = function->vector().vec().array();

   if ( fName[0] == 'c' )
   {
     //cout << "fc: " << function->name() << endl;

     VIter vit = M_vertexIter(MDBMesh);
     while (pVertex pv = VIter_next(vit))
     {
	// get dolfin point id
	pPoint pp = V_point(pv);
	int mdbId = P_id(pp);
	int dId = MDBToLocIds[mdbId];
	
	double data;
	EN_getDataDbl((pEntity)pv,dataId,&data);
	EN_deleteData((pEntity)pv,dataId);
	
	if(fName == "cUx")
	 data = dataArray[0 * N + dId] = (real)data;
	else if(fName == "cUy")
	 data = dataArray[1 * N + dId] = (real)data;
	else if(fName == "cUz")
	 data = dataArray[2 * N + dId] = (real)data;
	else if(fName == "ch0f")
	 data = dataArray[dId] = (real)data;
	else if(fName == "cphi")
	 data = dataArray[dId] = (real)data;
	else
	 dolfin::error("Unknown function");
     }
     VIter_delete(vit);
   }
   else if(fName[0] = 'd')
   {
     //cout << "fd: " << function->name() << endl;


     int dim = getDim(MDBMesh);

     // phi
     if(dim==2)
     {
     }
     else if (dim==3)
     {
	RIter rit = M_regionIter(MDBMesh);
	int counter = 0;
       int mdbId = 0;
	while (pRegion pr = RIter_next(rit))
	{
	 // get dolfin point id
	 int dId = MDBToLocCellIds[mdbId];

	 //cout << "dId: " << dId << endl;

	 double data;
	 EN_getDataDbl((pEntity)pr,dataId,&data);
	 EN_deleteData((pEntity)pr,dataId);
	 
	 if(fName == "dphi")
	 {
	   dataArray[dId] = (real)data;
	 }
	 else if(fName == "dS00")
	 {
	   dataArray[0 * M + dId] = (real)data;
	 }
	 else if(fName == "dS01")
	 {
	   dataArray[1 * M + dId] = (real)data;
	 }
	 else if(fName == "dS02")
	 {
	   dataArray[2 * M + dId] = (real)data;
	 }
	 else if(fName == "dS10")
	 {
	   dataArray[3 * M + dId] = (real)data;
	 }
	 else if(fName == "dS11")
	 {
	   dataArray[4 * M + dId] = (real)data;
	 }
	 else if(fName == "dS12")
	 {
	   dataArray[5 * M + dId] = (real)data;
	 }
	 else if(fName == "dS20")
	 {
	   dataArray[6 * M + dId] = (real)data;
	 }
	 else if(fName == "dS21")
	 {
	   dataArray[7 * M + dId] = (real)data;
	 }
	 else if(fName == "dS22")
	 {
	   dataArray[8 * M + dId] = (real)data;
	 }
	 else
	   dolfin::error("Unknown function");
	 
	 mdbId++;
	}
	//cout << "counter: " << mdbId << endl;
     }
     else
     {
	error("unknown function");
     }
   }

   MD_deleteMeshDataId(dataId);
 }
}

//-----------------------------------------------------------------------------
void MeshAdaptInterface::constrainExternalBoundaries()
{
 VIter vit = M_vertexIter(MDBMesh);
 while (pVertex pv = VIter_next(vit)) {

   // get dolfin point id
   pPoint pp = V_point(pv);
   int mdbId = P_id(pp);
   int dId = MDBToLocIds[mdbId];

   // if on boundary, constrain the point
   if ( boundaryVertices[dId] || constrainedVertices[dId])
   {
     adapter->setConstraint(pv);
     EN_constrain(pv);
   }
 }
 VIter_delete(vit);

 int dim = getDim(MDBMesh);
 switch (dim) {
 case 2: 
   {
     EIter eit = M_edgeIter(MDBMesh);
     while (pEdge pe = EIter_next(eit)) {
	
	pVertex pV0 = E_vertex((pEdge)pe, 0);
	pVertex pV1 = E_vertex((pEdge)pe, 1);

	// get dolfin point id
	pPoint pp0 = V_point(pV0);
	pPoint pp1 = V_point(pV1);
	int mdbId0 = P_id(pp0);
	int mdbId1 = P_id(pp1);
	int dId0 = MDBToLocIds[mdbId0];
	int dId1 = MDBToLocIds[mdbId1];

	if ( (boundaryVertices[dId0] || constrainedVertices[dId0]) &&
	    (boundaryVertices[dId1] || constrainedVertices[dId1]))
	{
	 adapter->setConstraint(pe);
	 EN_constrain(pe);
	}
     }
     EIter_delete(eit);
   }
   break;
 case 3:
   {
     EIter eit = M_edgeIter(MDBMesh);
     while (pEdge pe = EIter_next(eit)) {
	
	pVertex pV0 = E_vertex((pEdge)pe, 0);
	pVertex pV1 = E_vertex((pEdge)pe, 1);

	// get dolfin point id
	pPoint pp0 = V_point(pV0);
	pPoint pp1 = V_point(pV1);
	int mdbId0 = P_id(pp0);
	int mdbId1 = P_id(pp1);
	int dId0 = MDBToLocIds[mdbId0];
	int dId1 = MDBToLocIds[mdbId1];

	if ( (boundaryVertices[dId0] || constrainedVertices[dId0]) &&
	    (boundaryVertices[dId1] || constrainedVertices[dId1]))
	{
	 adapter->setConstraint(pe);
	 EN_constrain(pe);
	}
     }
     EIter_delete(eit);


     FIter fit = M_faceIter(MDBMesh);
     while (pFace pf = FIter_next(fit)) {
	
	pVertex pV0 = F_vertex((pFace)pf, 0);
	pVertex pV1 = F_vertex((pFace)pf, 1);
	pVertex pV2 = F_vertex((pFace)pf, 2);

	// get dolfin point id
	pPoint pp0 = V_point(pV0);
	pPoint pp1 = V_point(pV1);
	pPoint pp2 = V_point(pV2);
	int mdbId0 = P_id(pp0);
	int mdbId1 = P_id(pp1);
	int mdbId2 = P_id(pp2);
	int dId0 = MDBToLocIds[mdbId0];
	int dId1 = MDBToLocIds[mdbId1];
	int dId2 = MDBToLocIds[mdbId2];

	if ( (boundaryVertices[dId0] || constrainedVertices[dId0]) &&
	    (boundaryVertices[dId1] || constrainedVertices[dId1]) &&
	    (boundaryVertices[dId2] || constrainedVertices[dId2]))
	{
	 adapter->setConstraint(pf);
	 EN_constrain(pf);
	}
     }
     FIter_delete(fit);

   }
   break;
 }
//   int dim = getDim(MDBMesh);
//   switch (dim) {
//   case 2: 
//     {
//       EIter eit = M_edgeIter(MDBMesh);
//       while (pEdge pe = EIter_next(eit)) {

//         // if on boundary, constrain the edge
//         if ( E_numFaces(pe) == 1 ) {
//           adapter->setConstraint(pe);
//         }
//       }
//       EIter_delete(eit);
//     }
//     break;
//   case 3: 
//     {
//       FIter fit = M_faceIter(MDBMesh);
//       while (pFace pf = FIter_next(fit)) {

//         // if on boundary, constrain the edge
//         if ( F_numRegions(pf) == 1 ) {
//           adapter->setConstraint(pf);
//         }
//       }
//       FIter_delete(fit);
//     }
//     break;
//   }
}
//-----------------------------------------------------------------------------
void MeshAdaptInterface::constrainInternalBoundaries()
{
}
//-----------------------------------------------------------------------------
void MeshAdaptInterface::adaptMesh()
{
 // Warning: memory peak!
 // Could be avoided by the deletion of the mesh-dependent data here

 // build the MDB mesh
 //MDBMesh->expand();

 boost::timer local_timer;

 if(sizeField)
 {
   delete sizeField;
   sizeField = 0;
 }

 if(adapter)
 {
   delete adapter;
   adapter = 0;
 }

 if (MDBMesh)
 {
   delete MDBMesh;
   MDBMesh = 0;
 }

 if (model)
 {
   delete model;
   model = 0;
 }

 MDBMesh = MS_newMesh(model);
 exportToMDBMesh();

 sizeField = new PWLSField(MDBMesh);
 sizeField->setCurrentSize();

 adapter = new MeshAdapter(MDBMesh,sizeField);

 //adapter->addCallback(DolfinCBFunction,(void*)this);


 local_timer.restart();

 updateSizeField();

 cout << "MA update h timer: " << local_timer.elapsed() << endl;

 local_timer.restart();

 // constrain entities not to be adapted
 constrainExternalBoundaries();
 constrainInternalBoundaries();

 cout << "MA constrain timer: " << local_timer.elapsed() << endl;

 local_timer.restart();

 // move nodes
 prescribeNodesMotion(); // here all nodes, but could be interface nodes
 updateMDBCoordinates(); // or moveAndReposition in LMM

 cout << "MA move nodes timer: " << local_timer.elapsed() << endl;

 local_timer.restart();

 // if callback function is used:
 attachSolutionToMDBMesh();
 // deallocate solver datas
 deallocateData();

 cout << "MA attach/deallocate timer: " << local_timer.elapsed() << endl;

 //deleteDMesh();

 local_timer.restart();
 // run the adaptation
 adapter->run();

 //adapter->optimiseElementShape();
//   adapter->optimiseElementShape();

 adapter->printStatistics(std::cout);

 cout << "MA adapt timer: " << local_timer.elapsed() << endl;

 //M_writeSMS(MDBMesh,"test.msh",2,NULL);

 local_timer.restart();

 importMDBMesh();

 cout << "MA import timer: " << local_timer.elapsed() << endl;

 local_timer.restart();

 // allocate and compute solver datas
 allocateAndComputeData();

 cout << "MA allocate timer: " << local_timer.elapsed() << endl;

 local_timer.restart();

 // get the solution from the MDB mesh
 getSolutionFromMDBMesh();

 cout << "MA get solution timer: " << local_timer.elapsed() << endl;

 //  MDBMesh->shrink();
}