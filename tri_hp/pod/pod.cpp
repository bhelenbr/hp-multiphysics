/*
 *  pod.cpp
 *  tri_hp
 *
 *  Created by Brian Helenbrook on 7/26/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include <myblas.h>

template<class BASE> void pod<BASE>::init(input_map& input, gbl *gin) {
   bool coarse;
   std::string filename,keyword,linebuff;
   std::ostringstream nstr;
   std::istringstream instr;
   int i;

   /* Initialize base class */
   BASE::init(input,gin);
   
   keyword = BASE::idprefix + ".coarse";
   input.getwdefault(keyword,coarse,false);

   if (coarse) return;
   
   keyword = BASE::idprefix + ".snapshots";
   if (!input.get(keyword,nsnapshots)) {
      input.getwdefault("snapshots",nsnapshots,10);
   }
   
   if (!input.get(BASE::idprefix + ".podmodes",nsnapshots)) input.getwdefault("podmodes",nmodes,nsnapshots);   
   
   /* THIS IS TO CHANGE THE WAY SNAPSHOT MATRIX ENTRIES ARE FORMED */
   scaling.resize(BASE::NV);
   scaling = 1;
   if (input.getline(BASE::idprefix + ".scale_vector",linebuff) || input.getline("scale_vector",linebuff)) {
      instr.str(linebuff);
      for(i=0;i<BASE::NV;++i)
         instr >> scaling(i);
   }
   
   nmodes = MAX(nmodes,2);
   modes.resize(nmodes);
   for(i=0;i<nmodes;++i) {
      modes(i).v.resize(BASE::maxvst,BASE::NV);
      modes(i).s.resize(BASE::maxvst,BASE::sm0,BASE::NV);
      modes(i).i.resize(BASE::maxvst,BASE::im0,BASE::NV);
   }
   coeffs.resize(nmodes);
   modes_set = false;
   
   return;
}
   
template<class BASE> block::ctrl pod<BASE>::tadvance(bool coarse,block::ctrl ctrl_message,Array<mesh::transfer,1> &fv_to_ct,Array<mesh::transfer,1> &cv_to_ft, tri_hp *fmesh) {
   int i,j,k,l,n,tind,info;
   int lgpx = basis::tri(BASE::log2p).gpx, lgpn = basis::tri(BASE::log2p).gpn;
   std::string filename,keyword,linebuff;
   std::ostringstream nstr;
   std::istringstream instr;
   char jobz[2] = "V", uplo[2] = "L";
   FLT cjcb;
   block::ctrl state;
   
   
   if (ctrl_message == block::begin) excpt = 0;
#ifdef CTRL_DEBUG
         *sim::log << "in pod::tadvance: ctrl_message: " << ctrl_message << " excpt: " << excpt << std::endl;
#endif
   
   psimatrix.resize(nsnapshots*(nsnapshots+1)/2);
   psimatrix_recv.resize(nsnapshots*(nsnapshots+1)/2);

   switch (excpt) {
      case(0): {
         /* ADVANCE PHYSICAL PROBLEM */
         if (ctrl_message != block::advance1) {
            state = BASE::tadvance(coarse,ctrl_message,fv_to_ct,cv_to_ft,fmesh);
            if (state != block::stop) return(state);
            return(block::advance1);
         }
         excpt += 1 +modes_set*2;
      }

      case(1): {
         int psi1dcounter = 0;
         vsi ugstore;
         ugstore.v.reference(BASE::ugbd(0).v);
         ugstore.s.reference(BASE::ugbd(0).s);
         ugstore.i.reference(BASE::ugbd(0).i);
         
         /* GENERATE POD MODES SNAPSHOT COEFFICIENTS */
         psimatrix = 0.0;
         for (k=0;k<nsnapshots;++k) {
            nstr.str("");
            nstr << k+1 << std::flush;
            filename = BASE::idprefix +"_rstrt" +nstr.str() +".d0";
            BASE::ugbd(0).v.reference(modes(0).v);
            BASE::ugbd(0).s.reference(modes(0).s);
            BASE::ugbd(0).i.reference(modes(0).i);
            BASE::input(filename, BASE::text);
            
            for(l=k;l<nsnapshots;++l) {
               nstr.str("");
               nstr << l+1 << std::flush;
               filename = BASE::idprefix +"_rstrt" +nstr.str() +".d0";
               BASE::ugbd(0).v.reference(modes(1).v);
               BASE::ugbd(0).s.reference(modes(1).s);
               BASE::ugbd(0).i.reference(modes(1).i);
               BASE::input(filename, BASE::text);
               
               for(tind=0;tind<BASE::ntri;++tind) {        
                  /* LOAD ISOPARAMETRIC MAPPING COEFFICIENTS */
                  BASE::crdtocht(tind);

                  /* PROJECT COORDINATES AND COORDINATE DERIVATIVES TO GAUSS POINTS */
                  for(n=0;n<BASE::ND;++n)
                     basis::tri(BASE::log2p).proj_bdry(&BASE::cht(n,0), &BASE::crd(n)(0,0), &BASE::dcrd(n,0)(0,0), &BASE::dcrd(n,1)(0,0),MXGP);
                               
                  BASE::ugbd(0).v.reference(modes(0).v);
                  BASE::ugbd(0).s.reference(modes(0).s);
                  BASE::ugbd(0).i.reference(modes(0).i);
                  BASE::ugtouht(tind,0);
                  for(n=0;n<BASE::NV;++n)
                     basis::tri(BASE::log2p).proj(&BASE::uht(n)(0),&BASE::u(n)(0,0),MXGP);
                     
                  BASE::ugbd(0).v.reference(modes(1).v);
                  BASE::ugbd(0).s.reference(modes(1).s);
                  BASE::ugbd(0).i.reference(modes(1).i);
                  BASE::ugtouht(tind,0);
                  for(n=0;n<BASE::NV;++n)
                     basis::tri(BASE::log2p).proj(&BASE::uht(n)(0),&BASE::res(n)(0,0),MXGP);
                  
                  for(i=0;i<lgpx;++i) {
                     for(j=0;j<lgpn;++j) {
                        cjcb = RAD(BASE::crd(0)(i,j))*basis::tri(BASE::log2p).wtx(i)*basis::tri(BASE::log2p).wtn(j)*(BASE::dcrd(0,0)(i,j)*BASE::dcrd(1,1)(i,j) -BASE::dcrd(1,0)(i,j)*BASE::dcrd(0,1)(i,j));
                        for(n=0;n<BASE::NV;++n) {
                           psimatrix(psi1dcounter) += BASE::u(n)(i,j)*BASE::res(n)(i,j)*scaling(n)*cjcb;
                        }
                     }
                  }
               }
               ++psi1dcounter;
            }	
            BASE::ugbd(0).v.reference(ugstore.v);
            BASE::ugbd(0).s.reference(ugstore.s);
            BASE::ugbd(0).i.reference(ugstore.i);
         }
         sim::blks.allreduce1(psimatrix.data(),psimatrix_recv.data());
         
         return(block::advance);
      }
      case(2): {
         sim::blks.allreduce2(psimatrix.extent(firstDim),blocks::flt_msg,blocks::sum);
         
         vsi ugstore;
         ugstore.v.reference(BASE::ugbd(0).v);
         ugstore.s.reference(BASE::ugbd(0).s);
         ugstore.i.reference(BASE::ugbd(0).i);
         Array<FLT,1> eigenvalues(nsnapshots);
         Array<FLT,2> eigenvectors(nsnapshots,nsnapshots);
         Array<FLT,1> work(3*nsnapshots);
         
         DSPEV(jobz,uplo,nsnapshots,psimatrix_recv.data(),eigenvalues.data(),eigenvectors.data(),nsnapshots,work.data(),info);
         if (info != 0) {
            *sim::log << "Failed to find eigenmodes " << info << std::endl;
            exit(1);
         }
         
         *sim::log << "eigenvalues "<<  eigenvalues << std::endl;
         
         eigenvectors.transposeSelf(secondDim,firstDim);  // FORTRAN ROUTINE RETURNS TRANSPOSE

         //reconstruct POD MODES
         for(k=0;k<nmodes;++k)	{
            
            BASE::ugbd(0).v.reference(ugstore.v);
            BASE::ugbd(0).s.reference(ugstore.s);
            BASE::ugbd(0).i.reference(ugstore.i);
         
            modes(k).v = 0.0;
            modes(k).s = 0.0;
            modes(k).i = 0.0;
            
            for(l=0;l<nsnapshots;++l)	{
               nstr.str("");
               nstr << l +1 << std::flush;
               filename = BASE::idprefix +"_rstrt" +nstr.str() +".d0";
               BASE::input(filename, BASE::text);

               modes(k).v += eigenvectors(l,k)*BASE::ug.v;
               modes(k).s += eigenvectors(l,k)*BASE::ug.s;
               modes(k).i += eigenvectors(l,k)*BASE::ug.i;
            }
            
            nstr.str("");
            nstr << k << std::flush;
            filename = BASE::idprefix +"_mode" +nstr.str();
            BASE::ugbd(0).v.reference(modes(k).v);
            BASE::ugbd(0).s.reference(modes(k).s);
            BASE::ugbd(0).i.reference(modes(k).i);
            BASE::output(filename, BASE::text);
            BASE::output(filename, BASE::tecplot);
         }
         BASE::ugbd(0).v.reference(ugstore.v);
         BASE::ugbd(0).s.reference(ugstore.s);
         BASE::ugbd(0).i.reference(ugstore.i);
         modes_set = true;
         
         return(block::stop);
      }
   }
         
   return(block::stop);
}
   
   
//block::ctrl tri_hp_ins::rsdl(block::ctrl ctrl_message, int stage) {
//   
//
//
//   return(block::stop);
//}

   
   
   
   
   
   
   
   
   
   
   