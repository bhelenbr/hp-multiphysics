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

#include "MAdLibInterface.h"
#include "tet_mesh.h"

using namespace MAd;

void tet_mesh::MAdLib_input(const std::string filename, FLT grwfac, input_map& input) {
	pGModel MAdModel = NULL;
	GM_create(&MAdModel,"theModel");
	GM_readFromMSH(MAdModel, filename.c_str());
	pMesh MAdMesh = M_new(MAdModel);
	M_load(MAdMesh,filename.c_str());
	
	if (!initialized) {
		maxvst = M_numTets(MAdMesh);
		maxvst = (maxvst > M_numEdges(MAdMesh) ? maxvst : M_numEdges(MAdMesh));
		maxvst = (maxvst > M_numFaces(MAdMesh) ? maxvst : M_numFaces(MAdMesh));
		maxvst = (maxvst > M_numVertices(MAdMesh) ? maxvst : M_numVertices(MAdMesh));
		allocate(grwfac*maxvst);	

		/* Now Generate Boundary Structures */
		int Ntotal = M_numGeomFeatures(MAdMesh);
		nvbd = 0;
		nebd = 0;
		nfbd = 0;
		for (int i=0;i<Ntotal;) {
			int thisID = M_geomFeatureId(MAdMesh,i);
			pPGList list = M_geomFeature(MAdMesh, thisID);
			void * tmp = NULL;
			pGEntity pge;
			while ( ( pge = PGList_next(list,&tmp) ) ) {
// 			std::vector<pGEntity>::iterator it = list->entities.begin();
// 			for(std::vector<pGEntity>::iterator it = list->entities.begin(); it != list->entities.end(); ++it) {
				int gDim = GEN_type(pge);
// 				int gDim = GEN_type(*it);						
				switch (gDim) {
					case(0):
						++nvbd;
						break;
					case(1):
						++nebd;
						break;
					case(2):
						++nfbd;
						break;
				}
				++i;
			}
			PGList_delete(list);
		} 
		
		vbdry.resize(nvbd);
		ebdry.resize(nebd);
		fbdry.resize(nfbd);

		nvbd = 0;
		nebd = 0;
		nfbd = 0;
		for (int i=0;i<Ntotal;) {
			int thisID = M_geomFeatureId(MAdMesh,i);
			pPGList list = M_geomFeature(MAdMesh, thisID);
			void * tmp = NULL;
			pGEntity pge;
			while ( ( pge = PGList_next(list,&tmp) ) ) {
// 			std::vector<pGEntity>::iterator it = list->entities.begin();
// 			for(std::vector<pGEntity>::iterator it = list->entities.begin(); it != list->entities.end(); ++it) {
				int gDim = GEN_type(pge);
				int gId = GEN_tag(pge);			
				switch (gDim) {
					case(0):
						vbdry(nvbd) = getnewvrtxobject(gId,input);
						vbdry(nvbd)->alloc(4);
						++nvbd;
						break;
					case(1):
						ebdry(nebd) = getnewedgeobject(gId,input);
						ebdry(nebd)->alloc(grwfac*M_numClassifiedEdges(MAdMesh,pge));
						++nebd;
						break;
					case(2):
						fbdry(nfbd) = getnewfaceobject(gId,input);
						fbdry(nfbd)->alloc(grwfac*2*M_numClassifiedEdges(MAdMesh,pge));
						++nfbd;
						break;
				}
				++i;
			}
			PGList_delete(list);
		}
	}
	
	MAdLibInterface::importFromMAdMesh(MAdMesh,this);
}
	


//-----------------------------------------------------------------------------
// Main routine for adaptation
void MAdLibInterface::coarsenMesh(FLT factor, tet_mesh* mesh)
{
	//-----------------------------------------------------
	// Step 1: Prepare for adaptation
	//-----------------------------------------------------

	// 1.A. (optional): delete mesh/solution dependent data in 
	//									the solver to free memory
	// solver->deleteData();

	// 1.B. Build the MAdLib geometrical model.
	//			(see note in MAdLibInterface.h about models)
	pGModel MAdModel = NULL;
	GM_create(&MAdModel,"theModel");
	exportToMAdModel(mesh, MAdModel);

	// 1.C. Build the MadLib mesh and delete the solver mesh.
	pMesh MAdMesh = M_new(MAdModel);
	exportToMAdMesh(mesh, MAdMesh);
	M_writeMsh(MAdMesh,"starting_mesh.msh",2);

	// solver->deleteMesh();

	// 1.D. Build the size field used in adaptation
	PWLSField * sizeField = new PWLSField(MAdMesh);
	sizeField->setCurrentSize();
	sizeField->scale(factor);
	
	// 1.E. Build the adaptation tool
	MeshAdapter * adapter = new MeshAdapter(MAdMesh,sizeField);
	adapter->setMaxIterationsNumber(10);
	adapter->setCollapseOnBoundary(true,1.0);
	
	// 1.F. Register the callback function(s) of the solver
	// adapter->addCallback(Solver_CBFunction,(void*)this);
	
	// 1.G. Transfer solution to the MAdLib mesh as an attached data
	// attachSolutionToMesh(MAdMesh);
	// solver->deallocateSolution();

	//-----------------------------------------------------
	// Step 2: Run the adaptation
	//-----------------------------------------------------

	adapter->run();

	// optional output
	adapter->printStatistics(std::cout);
	M_writeMsh(MAdMesh,"adapted_mesh.msh",2);

	delete adapter;
	delete sizeField;

	//-----------------------------------------------------
	// Step 3: Rebuild solver data and mesh
	//-----------------------------------------------------

	// 3.A. Rebuild the solver mesh
	importFromMAdMesh(MAdMesh, mesh);
	
	// 3.B. get the solution from the MDB mesh
//	solver->allocateSolution();
//	getSolutionFromMesh(MAdMesh);

	// 3.C. Clean up MAdLib stuff
	delete MAdMesh;
	delete MAdModel;

	// 3.D. (optional with 1.A.) build mesh/solution 
	//													dependent data in the solver
//	solver->allocateAndComputeData();
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Converts a MAdLib mesh into a 'Solver_mesh'
void MAdLibInterface::importFromMAdMesh(const MAd::pMesh MAdMesh, tet_mesh* mesh) {
	std::map<int,int> MAdToSolverIds;
	std::map<int,int> SolverToMAdIds;

	// --- Build vertices ---
	mesh->npnt = 0;
	VIter vit = M_vertexIter(MAdMesh);
	while (pVertex pv = VIter_next(vit)) {
		// get MAdLib mesh id
		int MAd_Id = EN_id((pEntity)pv);

		pGEntity pge = V_whatIn(pv);
		int gDim = GEN_type(pge);
		int gId = GEN_tag(pge);

		if (gDim == 0) {
			// Vertex Boundary */
			for (int i=0;i<mesh->nvbd;++i) {
				if (gId == mesh->vbdry(i)->idnum) {
					mesh->vbdry(i)->pnt = mesh->npnt;
					break;
				}
			}
		}

		// get coordinates
		TinyVector<FLT,3> xyz;
		V_coord(pv,xyz.data());

		// fill in id's tables
		MAdToSolverIds[MAd_Id] = mesh->npnt;
		SolverToMAdIds[mesh->npnt] = MAd_Id;

		// add node in solver mesh
		mesh->pnts(mesh->npnt++) = xyz;
	}
	VIter_delete(vit);
	
	
	// --- Build edges ---
	mesh->nseg = 0;
	for (int i=0;i<mesh->nebd;++i) {
		mesh->ebdry(i)->nseg = 0;
	}
	EIter eit = M_edgeIter(MAdMesh);
	while (pEdge pe = EIter_next(eit)) {
		
		pGEntity pge = E_whatIn(pe);
		int gDim = GEN_type(pge);
		int gId = GEN_tag(pge);

		if (gDim == 1) {
			// Edge Boundary */
			for (int i=0;i<mesh->nebd;++i) {
				if (gId == mesh->ebdry(i)->idnum) {
					mesh->ebdry(i)->seg(mesh->ebdry(i)->nseg++).gindx = mesh->nseg;
					break;
				}
			}
		}		

		
		// get list of node id's in the solver mesh
		TinyVector<int,2> nodes;
		for (int i=0;i<2;++i) {
				pVertex pv = E_vertex(pe,i);
				int MAd_Id = EN_id((pEntity)pv);
				nodes(i) = MAdToSolverIds[MAd_Id];
		}
		
		// add node in solver mesh
		mesh->seg(mesh->nseg++).pnt = nodes;
	}
	EIter_delete(eit);
	
	// --- Build faces ---
	mesh->ntri = 0;
	for (int i=0;i<mesh->nfbd;++i) {
		mesh->fbdry(i)->ntri = 0;
	}
	FIter fit = M_faceIter(MAdMesh);
	while (pFace pf = FIter_next(fit)) {
		pGEntity pge = F_whatIn(pf);
		int gDim = GEN_type(pge);
		int gId = GEN_tag(pge);
		
		if (gDim == 2) {
			// Face Boundary */
			for (int i=0;i<mesh->nfbd;++i) {
				if (gId == mesh->fbdry(i)->idnum) {
					mesh->fbdry(i)->tri(mesh->fbdry(i)->ntri++).gindx = mesh->ntri;
					break;
				}
			}
		}
	
		// get list of node id's in the solver mesh
		TinyVector<int,3> nodes;
		pPList fVerts = F_vertices(pf,0);
		void * temp = NULL;
		int iN = 0;
		while ( pVertex pv = (pVertex)PList_next(fVerts,&temp) ) {
				int MAd_Id = EN_id((pEntity)pv);
				nodes(iN++) = MAdToSolverIds[MAd_Id];
		}
		PList_delete(fVerts);


		// add the element to the solver mesh
		mesh->tri(mesh->ntri++).pnt = nodes;
	}
	
	// --- Build elements ---
	mesh->ntet = 0;
	RIter rit = M_regionIter(MAdMesh);
	while (pRegion pr = RIter_next(rit)) {
		pGEntity pge = R_whatIn(pr);
		int gDim = GEN_type(pge);
		int gId = GEN_tag(pge);
		
		if (gDim == 2) {
			// Edge Boundary */
			for (int i=0;i<mesh->nfbd;++i) {
				if (gId == mesh->fbdry(i)->idnum) {
					mesh->fbdry(i)->tri(mesh->fbdry(i)->ntri++).gindx = mesh->ntet;
					break;
				}
			}
		}
	
		// get list of node id's in the solver mesh
		TinyVector<int,4> nodes;
		pPList rVerts = R_vertices(pr);
		void * temp = NULL;
		int iN = 0;
		while ( pVertex pv = (pVertex)PList_next(rVerts,&temp) ) {
				int MAd_Id = EN_id((pEntity)pv);
				nodes(iN++) = MAdToSolverIds[MAd_Id];
		}
		PList_delete(rVerts);


		// add the element to the solver mesh
		mesh->tet(mesh->ntet++).pnt = nodes;
	}
	RIter_delete(rit);
	
	for (int i=0;i<mesh->nebd;++i) {
		mesh->ebdry(i)->setup_next_prev();
		mesh->ebdry(i)->reorder();
	}
	
	mesh->reorient_tets();
	mesh->match_all();
	
	for(int i = 0; i < mesh->nfbd; ++i) {
		/* Load global vertex info */
		for(int j=0;j<mesh->fbdry(i)->ntri;++j) {
			int gindx =mesh->fbdry(i)->tri(j).gindx;
			mesh->fbdry(i)->tri(j).pnt = mesh->tri(gindx).pnt;
		}
		mesh->fbdry(i)->create_from_tri();
	}
	
	mesh->checkintegrity();
}

//-----------------------------------------------------------------------------
// Converts a 'tet_mesh' into a MAdLib mesh
void MAdLibInterface::exportToMAdMesh(const tet_mesh* mesh, MAd::pMesh MAdMesh) {
	
	// --- Build the vertices ---
	int nVerts = mesh->npnt;
	for (int i=0; i < nVerts; ++i) {
		MAdMesh->add_point(i,mesh->pnts(i)(0),mesh->pnts(i)(1),mesh->pnts(i)(2));
	}
	
	for (int i=0;i<mesh->nvbd;++i) {
		int idnum = mesh->vbdry(i)->idnum;
		pGVertex geom = GM_vertexByTag(MAdMesh->model,idnum);
		pVertex pv = MAdMesh->find_point(mesh->vbdry(i)->pnt);
		V_setWhatIn(pv, geom);
	}
	
	// --- Build the edges
	for (int i=0;i<mesh->nseg;++i)
		MAdMesh->add_edge(mesh->seg(i).pnt(0), mesh->seg(i).pnt(1));
		
	for (int i=0;i<mesh->nebd;++i) {
		int idnum = mesh->ebdry(i)->idnum;
		pGEdge geom = GM_edgeByTag(MAdMesh->model,idnum);
		for (int j=0;j<mesh->ebdry(i)->nseg;++j) {
			int gindx = mesh->ebdry(i)->seg(j).gindx;
			MAdMesh->add_edge(mesh->seg(gindx).pnt(0), mesh->seg(gindx).pnt(1),geom);
		}
	}	
	
	// --- Build the faces ---
	for (int i=0;i<mesh->ntri;++i)
		MAdMesh->add_triangle(mesh->tri(i).pnt(0), mesh->tri(i).pnt(1),mesh->tri(i).pnt(2));
		
	for (int i=0;i<mesh->nfbd;++i) {
		int idnum = mesh->fbdry(i)->idnum;
		pGFace geom = GM_faceByTag(MAdMesh->model,idnum);
		for (int j=0;j<mesh->fbdry(i)->ntri;++j) {
			int gindx = mesh->fbdry(i)->tri(j).gindx;
			MAdMesh->add_triangle(mesh->tri(gindx).pnt(0), mesh->tri(gindx).pnt(1),mesh->tri(gindx).pnt(2),geom);
		}
	}

	// --- Build the elements ---
	pGRegion geom = GM_regionByTag(MAdMesh->model,mesh->gbl->idnum);
	for (int i=0; i < mesh->ntet; ++i) {
		MAdMesh->add_tet(mesh->tet(i).pnt(0), mesh->tet(i).pnt(1),
										mesh->tet(i).pnt(2), mesh->tet(i).pnt(3),
										(pGEntity)geom); 
	}

	MAdMesh->classify_unclassified_entities();
	MAdMesh->destroyStandAloneEntities();
}

//-----------------------------------------------------------------------------
// Create in MAdModel all geometrical entities listed in solverModel.
void MAdLibInterface::exportToMAdModel(const tet_mesh *mesh, MAd::pGModel MAdModel) {
	for (int i=0;i<mesh->nvbd;++i)
		GM_entityByTag(MAdModel,0,mesh->vbdry(i)->idnum);

	for (int i=0;i<mesh->nebd;++i)
		GM_entityByTag(MAdModel,1,mesh->ebdry(i)->idnum);

	for (int i=0;i<mesh->nfbd;++i)
		GM_entityByTag(MAdModel,2,mesh->fbdry(i)->idnum);
		
	GM_entityByTag(MAdModel,3,mesh->gbl->idnum);
}

//-----------------------------------------------------------------------------
// Build a field of prescribed edges lengths on the domain.
void MAdLibInterface::buildSizeField(MAd::PWLSField * sizeField)
{
	// First option: keep actual edges lengths
	sizeField->setCurrentSize();
	sizeField->scale(0.5);
	
//	// Second option: compute it from solver functions
//	VIter vit = M_vertexIter(sizeField->getMesh());
//	while (pVertex pv = VIter_next(vit))
//		{
//			// get solver point id
//			int MAd_Id = EN_id((pEntity)pv);
//			void setSize(pEntity, double);
//			sizeField->setSize((pEntity)pv, length);
//
//			
//			
//			int solver_Id = MAdToSolverIds[MAd_Id];
//			
//			// get the edge length prescribed by the solver
//			double length = solver->wishEdgeLength(solver_Id);
//
//			// fill in the size field
//			sizeField->setSize((pEntity)pv, length);
//		}
//	VIter_delete(vit);
}

//-----------------------------------------------------------------------------
//void MAdLibInterface::attachSolutionToMesh(MAd::pMesh MAdMesh)
//{
//	// Get the solution database. Here we assume that it is a nodal solution.
//	const Solver_solution * solution = solver->getSolution();
//
//	// The data id used to identify the data attached to mesh entities
//	pMeshDataId dataId = MD_lookupMeshDataId("SolutionTag");
//
//	VIter vit = M_vertexIter(MAdMesh);
//	while (pVertex pv = VIter_next(vit))
//		{
//			// get solver point id
//			int MAd_Id = EN_id((pEntity)pv);
//			int solver_Id = MAdToSolverIds[MAd_Id];
//			
//			double data = (*solution)[solver_Id];
//			
//			// attach data to the mesh vertex
//			EN_attachDataDbl((pEntity)pv,dataId,data);
//		}
//	VIter_delete(vit);
//}
//
////-----------------------------------------------------------------------------
//void MAdLibInterface::getSolutionFromMesh(MAd::pMesh MAdMesh)
//{
//	// Get the solution database. Here we assume that it is a nodal solution.
//	Solver_solution * solution = solver->getSolution();
//
//	// The data id used to identify the data attached to mesh entities
//	pMeshDataId dataId = MD_lookupMeshDataId("SolutionTag");
//
//	VIter vit = M_vertexIter(MAdMesh);
//	while (pVertex pv = VIter_next(vit))
//		{
//			// get solver point id
//			pPoint pp = V_point(pv);
//			int MAdId = P_id(pp);
//			int solver_Id = MAdToSolverIds[MAdId];
//			
//			// get attached data and delete it
//			double data;
//			EN_getDataDbl((pEntity)pv,dataId,&data);
//			EN_deleteData((pEntity)pv,dataId);
//			
//			*(*solution)[solver_Id] = data;
//		}
//	VIter_delete(vit);
//}

//-----------------------------------------------------------------------------