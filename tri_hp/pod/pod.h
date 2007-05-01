/*
 *  pod.h
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 7/26/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

template<class BASE> class pod_generate : public BASE {
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
      
   private:
      int excpt;
   
   public:
      void init(input_map& input, gbl *gin); 
      pod_generate<BASE>* create() { return new pod_generate<BASE>();}
      block::ctrl tadvance(bool coarse,block::ctrl ctrl_message,Array<mesh::transfer,1> &fv_to_ct,Array<mesh::transfer,1> &cv_to_ft, tri_hp *fmesh);
};

template<class BASE> class pod_simulate : public BASE {
   public:
      typedef typename BASE::gbl gbl;

   protected:
      int nmodes;
      Array<FLT,1> coeffs, rsdls, rsdls_recv;
      typedef typename BASE::vsi vsi;
      Array<vsi,1> modes;
      Array<FLT,2> jacobian;
      Array<int,1> ipiv;
      
   private:
      int excpt,excpt1,modeloop;
   
   public:
      void init(input_map& input, gbl *gin); 
      pod_simulate<BASE>* create() { return new pod_simulate<BASE>();}
      block::ctrl rsdl(block::ctrl ctrl_message, int stage=sim::NSTAGE);
      block::ctrl setup_preconditioner(block::ctrl ctrl_message);
      block::ctrl update(block::ctrl ctrl_message);
      FLT maxres();      

};




#include "pod.cpp"