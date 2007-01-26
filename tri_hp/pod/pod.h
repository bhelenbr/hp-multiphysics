/*
 *  pod.h
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 7/26/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

template<class BASE> class pod : public BASE {
   public:
      typedef typename BASE::gbl gbl;

   protected:
      int nsnapshots;
      int nmodes;
      Array<FLT,1> scaling;
      Array<FLT,1> coeffs;
      typedef typename BASE::vsi vsi;
      Array<vsi,1> modes;
      Array<FLT,1> psimatrix,psimatrix_recv;
      bool modes_set;
      
   private:
      int excpt;
   
   public:
      void init(input_map& input, gbl *gin); 
      pod<BASE>* create() { return new pod<BASE>();}
      block::ctrl tadvance(bool coarse,block::ctrl ctrl_message,Array<mesh::transfer,1> &fv_to_ct,Array<mesh::transfer,1> &cv_to_ft, tri_hp *fmesh);

//      block::ctrl setup_preconditioner(block::ctrl ctrl_message);
//      block::ctrl rsdl(block::ctrl ctrl_message, int stage=sim::NSTAGE);
//      block::ctrl update(block::ctrl ctrl_message);
//      FLT maxres();      
//      
//      /* MGRID TRANSFER */
//      inline void setlog2p(int value) { BASE::log2p = value; } /* To switch polynomial degree */
//      block::ctrl mg_getfres(block::ctrl ctrl_message, Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, pod<BASE> *fmesh);
//      block::ctrl mg_getcchng(block::ctrl ctrl_message,Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, pod<BASE> *cmesh);
//
//      /* ADVANCE TIME SOLUTION */
//      block::ctrl tadvance(bool coarse,block::ctrl ctrl_message,Array<mesh::transfer,1> &fv_to_ct,Array<mesh::transfer,1> &cv_to_ft, pod<BASE> *fmesh);
//      virtual void calculate_unsteady_sources(bool coarse);
};

#include "pod.cpp"