#include "r_mesh.h"
#include "r_boundary.h"
#include "boundary.h"
#include "block.h"
#include <utilities.h>
#include <iostream>
#include <cmath>
#include <input_map.h>
#include <fstream>


sharedmem* r_mesh::init(bool coarse, std::map <std::string,std::string>& input, std::string prefix, r_mesh::gbl *rgin, sharedmem *wkin) {
   std::string keyword;
   std::istringstream data;
   std::string filename;
   
   keyword = prefix + ".fadd";
   data.str(input[keyword]);
   if (!(data >> fadd)) {
      data.clear();
      keyword = "fadd";
      data.str(input[keyword]);
      if (!(data >> fadd)) fadd = 1.0;
   }
   *log << "#fadd: " << fadd << std::endl;
   data.clear();

   keyword = prefix + ".vnn";
   data.str(input[keyword]);
   if (!(data >> vnn)) vnn = 0.5; 
   *log << "#vnn: " << vnn << std::endl;
   data.clear();
   
   /* local storage */   
   ksprg = new FLT[maxvst];
   kvol = new FLT[maxvst];
   src = (FLT (*)[ND]) xmalloc(maxvst*ND*sizeof(FLT));
   if (coarse) vrtx_frst = (FLT (*)[ND]) xmalloc(maxvst*ND*sizeof(FLT));
   isfrst = false;
   
   /* SHARED INFORMATION */
   rg = rgin;
   
   if (!coarse) {
      scratch.resize((2*ND+1)*maxvst*sizeof(FLT));
      scratch.checkout();
   }
   get_scratch_pointers();
      
   std::string bdryfile;
   std::map<std::string,std::string> bdrymap;
   keyword = prefix + ".bdryfile";
   data.str(input[keyword]);
   if (data >> bdryfile) {
      if (!(strncmp("${HOME}",bdryfile.c_str(), 7))) {
         filename = getenv("HOME");
         filename = filename + (bdryfile.c_str()+7);
      }
      else 
         filename = bdryfile;
   }
   else {
      keyword = prefix + ".mesh";
      data.str(input[keyword]);
      filename = filename +"_bdry.inpt";
   }
   input_map(bdrymap,filename.c_str());
   data.clear();

   r_sbdry.resize(nsbd);
   for(int i=0;i<nsbd;++i)
      r_sbdry(i) = getnewsideobject(i,&bdrymap);
   
   return(&scratch);
}

void r_mesh::get_scratch_pointers() {
   mesh::get_scratch_pointers();
   fscr2 = static_cast<FLT (*)[ND]>(scratch.data()) +maxvst;
   fscr3 = fscr2+maxvst;
}

void r_mesh::bdry_output(const char *filename) const {
   char fnmapp[120];
   std::ofstream bout;
   int i;
   
   strcpy(fnmapp,filename);
   strcat(fnmapp,"_bdry.inpt");
   bout.open(fnmapp);
   for(i=0;i<nvbd;++i) 
      vbdry(i)->output(bout);
      
   for(i=0;i<nsbd;++i) {
      sbdry(i)->output(bout);
      r_sbdry(i)->output(bout);
   }
   bout.close();
}


void r_mesh::rklaplace() {
   int sind,tind,v0,v1,k;
   FLT dx,dy,l;
   
   for(sind=0;sind<nside;++sind)
      ksprg[sind] = 0.0;

   /* COEFFICIENTS FOR LAPLACE EQUATION */
   /* THIS REQUIRES 2 EVALUATIONS OF SIDE LENGTH FOR EACH SIDE */
   /* BUT IS LOGISTICALLY SIMPLE        */         
   for(tind=0;tind<ntri;++tind) {
      for(k=0;k<3;++k) {
         sind = td(tind).side(k);
         v0 = sd(sind).vrtx(0);
         v1 = sd(sind).vrtx(1);
         dx = vrtx(v1)(0) -vrtx(v0)(0);
         dy = vrtx(v1)(1) -vrtx(v0)(1);      
         l  = (dx*dx +dy*dy)/area(tind);

         ksprg[sind] -= l;
         sind = td(tind).side((k+1)%3);
         ksprg[sind] += l;
         sind = td(tind).side((k+2)%3);
         ksprg[sind] += l;
      }
   }      

   calc_kvol();
   
   return;
}


void r_mesh::rksprg() {
   int sind,v0,v1;
   double dx,dy;

   /* 2D SPRING CONSTANTS FINE MESH*/
   for(sind=0;sind<nside;++sind) {
      v0 = sd(sind).vrtx(0);
      v1 = sd(sind).vrtx(1);
      dx = vrtx(v1)(0) -vrtx(v0)(0);
      dy = vrtx(v1)(1) -vrtx(v0)(1);
      ksprg[sind] = 1.0/(dx*dx +dy*dy);
   }

   calc_kvol();

   return;
}

void r_mesh::calc_kvol() {
   int i,tind;
   
   for(i=0;i<nvrtx;++i) 
      kvol[i] = 0.0;

   for(tind=0;tind<ntri;++tind) 
      for(i=0;i<3;++i) 
         kvol[td(tind).vrtx(i)] += area(tind);
           
   return;
}

void r_mesh::kvoli() {
   int i;

   for(i=0;i<nvrtx;++i)
      kvol[i] = 1./kvol[i];
}

void r_mesh::rkmgrid(Array<mesh::transfer,1> &fv_to_ct, r_mesh *fmesh) {
   int i,j,sind,tind,tind0,tind1,v0,v1;   

   /* TEMPORARILY USE DIAG TO STORE DIAGONAL SUM */
   for(i=0;i<fmesh->nvrtx;++i)
      fscr1(i) = 0.0;

   /* FORM KIJ SUM AT VERTICES */
   for(sind=0;sind<fmesh->nside;++sind) {
      v0 = fmesh->sd(sind).vrtx(0);
      v1 = fmesh->sd(sind).vrtx(1);
      fscr1(v0) += fmesh->ksprg[sind];
      fscr1(v1) += fmesh->ksprg[sind];
   }

   for(i=0;i<nside;++i)
      ksprg[i] = 0.0;

   /* LOOP THROUGH FINE VERTICES   */
   /* TO CALCULATE KSPRG ON COARSE MESH */   
   for(i=0;i<fmesh->nvrtx;++i) {
      tind = fv_to_ct(i).tri;
      for(j=0;j<3;++j) {
         sind = td(tind).side(j);
         ksprg[sind] -= fv_to_ct(i).wt(j)*fv_to_ct(i).wt((j+1)%3)*fscr1(i);
      }
   }

   /* LOOP THROUGH FINE SIDES */
   for(i=0;i<fmesh->nside;++i) {
      v0 = fmesh->sd(i).vrtx(0);
      v1 = fmesh->sd(i).vrtx(1);
      tind0 = fv_to_ct(v0).tri;
      tind1 = fv_to_ct(v1).tri;
                  
      /* TEMPORARILY STORE WEIGHTS FOR FINE POINTS (0,1) */
      /* FOR EACH COARSE VERTEX */
      for(j=0;j<3;++j)  {
         fscr2[td(tind1).vrtx(j)][0] = 0.0;
         fscr2[td(tind0).vrtx(j)][1] = 0.0;
      }

      for(j=0;j<3;++j)  {
         fscr2[td(tind0).vrtx(j)][0] = fv_to_ct(v0).wt(j);
         fscr2[td(tind1).vrtx(j)][1] = fv_to_ct(v1).wt(j);
      }
      
      /* LOOP THROUGH COARSE TRIANGLE 0 SIDES */
      for(j=0;j<3;++j) {
         sind = td(tind0).side(j);
         ksprg[sind] += fmesh->ksprg[i]*
            (fscr2[sd(sind).vrtx(0)][0]*fscr2[sd(sind).vrtx(1)][1]
            +fscr2[sd(sind).vrtx(1)][0]*fscr2[sd(sind).vrtx(0)][1]);
      }

      if (tind0 != tind1) {
         for(j=0;j<3;++j) {
            sind = td(tind1).side(j);
            if (sd(sind).tri(0) +sd(sind).tri(1) != tind0 +tind1) {
               ksprg[sind] += fmesh->ksprg[i]*
                  (fscr2[sd(sind).vrtx(0)][0]*fscr2[sd(sind).vrtx(1)][1]
                  +fscr2[sd(sind).vrtx(1)][0]*fscr2[sd(sind).vrtx(0)][1]);

            }
         }
      }
   }
   
   /* CALCULATE KVOL USING MULTIGRID I VOL I^T*/
   /* LOOP THROUGH FINE VERTICES TO CALCULATE COARSE VOLUME  */
   for(i=0;i<nvrtx;++i)
      kvol[i] = 0.0;

   for(i=0;i<fmesh->nvrtx;++i) {
      tind = fv_to_ct(i).tri;
      for(j=0;j<3;++j) {
         v0 = td(tind).vrtx(j);
         kvol[v0] += fv_to_ct(i).wt(j)/fmesh->kvol[i];
      }
   }

   return;
}

/*************************************/ 
/* SECOND ORDER MESH MOVEMENT SCHEME */
/*************************************/
#ifndef FOURTH
void r_mesh::rsdl() {
   int i,n,sind,v0,v1;
   FLT dx,dy;
   
   for(i=0;i<nvrtx;++i)
      for(n=0;n<ND;++n)
         fscr3[i][n] = 0.0;

   for(sind=0;sind<nside;++sind) {
      v0 = sd(sind).vrtx(0);
      v1 = sd(sind).vrtx(1);

      dx = ksprg[sind]*(vrtx(v1)(0)-vrtx(v0)(0));
      dy = ksprg[sind]*(vrtx(v1)(1)-vrtx(v0)(1));

      fscr3[v0][0] -= dx;
      fscr3[v0][1] -= dy;

      fscr3[v1][0] += dx;
      fscr3[v1][1] += dy;
   }
   
   /* CALCULATE DRIVING TERM ON FIRST ENTRY TO COARSE MESH */
   if (isfrst) {
      for(i=0;i<nvrtx;++i) 
         for(n=0;n<ND;++n)
            src[i][n] -= fscr3[i][n];

      isfrst = false;
   }

   /* ADD IN MULTIGRID SOURCE OR FMESH SOURCE */
   for(i=0;i<nvrtx;++i) 
      for(n=0;n<ND;++n)
         fscr3[i][n] += src[i][n];

   /* APPLY DIRICHLET BOUNDARY CONDITIONS */
   for(i=0;i<nsbd;++i)
      r_sbdry(i)->dirichlet();
      
         
   /* TESTING 
   for(i=0;i<nvrtx;++i) {
      fscr3[i][0] = 0.525*vrtx(i)(0)+0.7*i;
      fscr3[i][1] = 0.525*vrtx(i)(1)+0.7*i;
   }
   
   std::cout << "74 " << fscr3[nvrtx-17][0] << " " <<  fscr3[nvrtx-17][1] << std::endl;
   std::cout << "90 " << fscr3[nvrtx-1][0] << " " <<  fscr3[nvrtx-1][1] << std::endl;
   */
   
   return;
}
#endif

#ifndef FOURTH
void r_mesh::vddt(void)
{
   int i,v0,v1,sind;

   /**************************************************/
   /* DETERMINE MESH MOVEMENT TIME STEP           */
   /**************************************************/
   for(i=0;i<nvrtx;++i)
      fscr1(i) = 0.0;

   /* FORM TIME STEP FOR MV_UPDATE */
   for(sind=0;sind<nside;++sind) {
      v0 = sd(sind).vrtx(0);
      v1 = sd(sind).vrtx(1);
      fscr1(v0) += fabs(ksprg[sind]);
      fscr1(v1) += fabs(ksprg[sind]);
   }
   
   return;
}
#endif

void r_mesh::vddti() {
   int i;
      
   for(i=0;i<nvrtx;++i)
      fscr1(i) = vnn/fscr1(i);
}

block::ctrl r_mesh::update(int excpt) {
   int i,n;

   for(i=0;i<nvrtx;++i)
      for(n=0;n<ND;++n)
         vrtx(i)(n) -= fscr1(i)*fscr3[i][n];
      
   return(block::stop);
}

void r_mesh::zero_source() {
   int i,n;
   
   for(i=0;i<nvrtx;++i) 
      for(n=0;n<ND;++n) 
         src[i][n] = 0.0;
            
   return;
}

void r_mesh::sumsrc() {
   int i,n;

   for(i=0;i<nvrtx;++i) 
      for(n=0;n<ND;++n)
         src[i][n] = -1.0*fscr3[i][n];
   
   return;
}

block::ctrl r_mesh::mg_getfres(int excpt, Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, r_mesh *fmesh) {
   int i,j,n,tind,v0;
   
   FLT (*fscr3)[ND] = fmesh->fscr3;
      
   for(i=0;i<nvrtx;++i)
      for(n=0;n<ND;++n)
         src[i][n] = 0.0;
         
   /* LOOP THROUGH FINE VERTICES TO CALCULATE RESIDUAL  */
   for(i=0;i<fmesh->nvrtx;++i) {
      tind = fv_to_ct(i).tri;
      for(j=0;j<3;++j) {
         v0 = td(tind).vrtx(j);
         for(n=0;n<ND;++n)
            src[v0][n] += fadd*fv_to_ct(i).wt(j)*fscr3[i][n];
      }
   }
   
   /* LOOP THROUGH fv_to_ct VERTICES   */
   /* TO CALCULATE VRTX ON fv_to_ct MESH */
   for(i=0;i<nvrtx;++i) {
      tind = cv_to_ft(i).tri;

      for(n=0;n<ND;++n)
         vrtx(i)(n) = 0.0;
         
      for(j=0;j<3;++j) {
         for(n=0;n<ND;++n)
            vrtx(i)(n) += cv_to_ft(i).wt(j)*fmesh->vrtx(fmesh->td(tind).vrtx(j))(n);
      }
      
      for(n=0;n<ND;++n)
         vrtx_frst[i][n] = vrtx(i)(n);
   }

   isfrst = true;
   
   return(block::stop);
}

block::ctrl r_mesh::mg_getcchng(int excpt,Array<mesh::transfer,1> &fv_to_ct, Array<mesh::transfer,1> &cv_to_ft, r_mesh *cmesh) {
   int i,j,n,ind,tind;  
   int stop=1;
   
   
   switch (excpt) {
      case(0):
         mp_phase = -1;
         /* DETERMINE CORRECTIONS ON COARSE MESH   */   
         for(i=0;i<cmesh->nvrtx;++i)
            for(n=0;n<ND;++n) 
               cmesh->vrtx_frst[i][n] -= cmesh->vrtx(i)(n);

         /* LOOP THROUGH FINE VERTICES   */
         /* TO DETERMINE CHANGE IN SOLUTION */   
         for(i=0;i<nvrtx;++i) {
            
            for(n=0;n<ND;++n)
               fscr3[i][n] = 0.0;
            
            tind = fv_to_ct(i).tri;
            
            for(j=0;j<3;++j) {
               ind = cmesh->td(tind).vrtx(j);
               for(n=0;n<ND;++n) 
                  fscr3[i][n] -= fv_to_ct(i).wt(j)*cmesh->vrtx_frst[ind][n];
            }
         }
         return(block::advance);
            
      case(1):
         /* SEND COMMUNICATION PACKETS  */ 
         ++mp_phase;

         switch(mp_phase%3) {
            case(0):
               for(i=0;i<nsbd;++i)
                  if (sbdry(i)->group()&0x1) sbdry(i)->loadbuff((FLT *) fscr3,0,1,2);
                  
               for(i=0;i<nsbd;++i) 
                  if (sbdry(i)->group()&0x1) sbdry(i)->comm_prepare(mp_phase/3);
               
               return(block::stay);
            case(1):
               for(i=0;i<nsbd;++i) 
                  if (sbdry(i)->group()&0x1) sbdry(i)->comm_transmit(mp_phase/3);
               return(block::stay);
            case(2):
               stop = 1;
               for(i=0;i<nsbd;++i) {
                  if (sbdry(i)->group()&0x1) {
                     stop &= sbdry(i)->comm_wait(mp_phase/3);
                     sbdry(i)->finalrcv(mp_phase/3,(FLT *) fscr3,0,1,2);
                  }
               }
               return(static_cast<block::ctrl>(stop));
         }

      case(2):
         for(i=0;i<nvrtx;++i)
            for(n=0;n<ND;++n) 
               vrtx(i)(n) += fscr3[i][n];
               
         return(block::stop);
   }
   *log << "Flow control error\n" << std::endl;
   return(block::stop);
}

void r_mesh::maxres() {
   int i,n;
   FLT mxr[ND];

   for(n=0;n<ND;++n)
      mxr[n] = 0.0;

   for(i=0;i<nvrtx;++i)
      for(n=0;n<ND;++n)
         mxr[n] = MAX(mxr[n],fabs(fscr3[i][n]));
         
   for(n=0;n<ND;++n)
      *log << ' ' << mxr[n] << ' ';
         
   return;
}


#ifdef FOURTH
/*************************************/ 
/* FOURTH ORDER MESH MOVEMENT SCHEME */
/* THIS HAS TO BE DONE IN TWO PARTS  */
/* TO CACLUATE LAPLACIAN              */
/*************************************/
void r_mesh::rsdl1() {
   int i,n,v0,v1,sind;
   FLT dx,dy;
   
   for(i=0;i<nvrtx;++i)
      for(n=0;n<ND;++n)
         fscr2[i][n] = 0.0;

   for (sind = 0; sind < nside; ++sind) {
      v0 = sd(sind).vrtx(0);
      v1 = sd(sind).vrtx(1);

      dx = ksprg[sind]*(vrtx(v1)(0)-vrtx(v0)(0));
      dy = ksprg[sind]*(vrtx(v1)(1)-vrtx(v0)(1));

      fscr2[v0][0] -= dx;
      fscr2[v0][1] -= dy;

      fscr2[v1][0] += dx;
      fscr2[v1][1] += dy;      
   }
   
   /* APPLY DIRICHLET BOUNDARY CONDITIONS */
   for(i=0;i<nsbd;++i)
      r_sbdry(i)->fixdx2();

   return;
}

void r_mesh::rsdl() {
   int i,n,v0,v1,sind;
   FLT dx,dy;
  
   /* DIVIDE BY VOLUME FOR AN APPROXIMATION TO D^2/DX^2 */
   for(i=0;i<nvrtx;++i)
      for(n=0;n<ND;++n) 
         fscr2[i][n] *= kvol[i];
   
   for(i=0;i<nvrtx;++i)
      for(n=0;n<ND;++n)
         fscr3[i][n] = 0.0;

   for (sind = 0; sind < nside; ++sind) {
      v0 = sd(sind).vrtx(0);
      v1 = sd(sind).vrtx(1);

      dx = ksprg[sind]*(fscr2[v1][0]-fscr2[v0][0]);
      dy = ksprg[sind]*(fscr2[v1][1]-fscr2[v0][1]);

      fscr3[v0][0] -= dx;
      fscr3[v0][1] -= dy;

      fscr3[v1][0] += dx;
      fscr3[v1][1] += dy;      
   }

   /* CALCULATE DRIVING TERM ON FIRST ENTRY TO COARSE MESH */
   if (isfrst) {
      for(i=0;i<nvrtx;++i) 
         for(n=0;n<ND;++n)
            src[i][n] -= fscr3[i][n];

      isfrst = false;
   }
   
   /* ADD IN MULTIGRID SOURCE OR FMESH SOURCE */
   for(i=0;i<nvrtx;++i) 
      for(n=0;n<ND;++n)
         fscr3[i][n] += src[i][n];  

   /* APPLY DIRICHLET BOUNDARY CONDITIONS */
   for(i=0;i<nsbd;++i)
      r_sbdry(i)->dirichlet();
   
   return;
}

void r_mesh::vddt1(void)
{
   int i,v0,v1,sind;

   /**************************************************/
   /* DETERMINE MESH MOVEMENT TIME STEP           */
   /**************************************************/
   for(i=0;i<nvrtx;++i)
      fscr1(i) = 0.0;
      
   for(sind=0;sind<nside;++sind) {
      v0 = sd(sind).vrtx(0);
      v1 = sd(sind).vrtx(1);
      fscr1(v0) += fabs(ksprg[sind]);
      fscr1(v1) += fabs(ksprg[sind]);
   }

   return;
}

void r_mesh::vddt() {
   int i,v0,v1,sind;

   for(i=0;i<nvrtx;++i)
      fscr2[i][0] = fscr1(i)*kvol[i];
      
   for(i=0;i<nvrtx;++i)
      fscr1(i) = 0.0;
   
   for(sind=0;sind<nside;++sind) {
      v0 = sd(sind).vrtx(0);
      v1 = sd(sind).vrtx(1);
      fscr1(v0) += fabs(ksprg[sind])*(fscr2[v0][0] +fabs(ksprg[sind])*kvol[v1]);
      fscr1(v1) += fabs(ksprg[sind])*(fscr2[v1][0] +fabs(ksprg[sind])*kvol[v0]);
   }

   return;
}
#endif


block::ctrl r_mesh::tadvance(bool coarse,int execpoint,Array<mesh::transfer,1> &fv_to_ct,Array<mesh::transfer,1> &cv_to_ft, r_mesh *fmesh) {
   if (!coarse) {
      switch (execpoint) {
         case (0):
            /* SETUP SPRING CONSTANTS  */
            mp_phase = -1;
            rklaplace();
            return(block::advance);
            
         case (1):
           ++mp_phase;
            /* MESSAGE PASSING SEQUENCE */
            switch(mp_phase%3) {
               case(0):
                  msgload(mp_phase/3,kvol,0,0,1);
                  return(block::stay);
               case(1):
                  msgpass(mp_phase/3);
                  return(block::stay);
               case(2):
                  return(static_cast<block::ctrl>(msgwait_rcv(mp_phase/3,kvol,0,0,1)));
            }
            
         case (2):
            kvoli();
            return(block::advance);

#ifdef FOURTH
#define P2 2
         case (3):
            rsdl1();
            mp_phase = -1;
            return(block::advance);
            
         case (4):
            ++mp_phase;
            /* MESSAGE PASSING SEQUENCE */
            switch(mp_phase%3) {
               case(0):
                  msgload(mp_phase/3,(FLT *) fscr2,0,1,2);
                  return(block::stay);
               case(1):
                  msgpass(mp_phase/3);
                  return(block::stay);
               case(2):
                  return(static_cast<block::ctrl>(msgwait_rcv(mp_phase/3,(FLT *) fscr2,0,1,2));
            }
#else
#define P2 0
#endif            
         case (3+P2):
            mp_phase = -1;
            zero_source();
            rsdl();
            return(block::advance);

         case (4+P2):
            ++mp_phase;
            /* MESSAGE PASSING SEQUENCE */
            switch(mp_phase%3) {
               case(0):
                  msgload(mp_phase/3,(FLT *) fscr3,0,1,2);
                  return(block::stay);
               case(1):
                  msgpass(mp_phase/3);
                  return(block::stay);
               case(2):
                  return(static_cast<block::ctrl>(msgwait_rcv(mp_phase/3,(FLT *) fscr3,0,1,2)));
            }
         
         case (5+P2):
            sumsrc();
            moveboundaries();
            return(block::stop);

         default:
            *log << "error in control flow tadvance 1" << std::endl;
            exit(1);
      }
   }
   else {
#ifdef GEOMETRIC
      switch (execpoint) {
         case (0):
            /* SETUP SPRING CONSTANTS  */
            mp_phase = -1;
            rklaplace();
            return(block::advance);
            
         case (1):
            ++mp_phase;
            switch(mp_phase%3) {
               case(0):
                  msgload(mp_phase/3,kvol,0,0,1);
                  return(block::stay);
               case(1):
                  msgpass(mp_phase/3);
                  return(block::stay);
               case(2):
                  return(static_cast<block::ctrl>(msgwait_rcv(mp_phase/3,kvol,0,0,1)));
            }
         case (2):
            kvoli();
            return(block::stop);
         
         default:
            *log << "error in control flow tadvance 2" << std::endl;
            exit(1);
      }
   }
#else
   /* USE MULTIGRID INTERPOLATION (ALGEBRAIC) */
   /* MUST BE DONE THIS WAY FOR SPRING METHOD */
   /* SETUP FIRST MESH */
      switch (phase) {
         case(0):
            mp_phase = -1;
            rkmgrid(fv_to_ct,fmesh);
            return(block::advance);
         case(1):
            ++mp_phase
            switch(mp_phase%3) {
               case(0):
                  msgload(mp_phase/3,kvol,0,0,1);
                  return(block::stay);
               case(1):
                  msgpass(mp_phase/3);
                  return(block::stay);
               case(2):
                  return(static_cast<block::ctrl>(msgwait_rcv(mp_phase/3,kvol,0,0,1)));
            }
         case(2):
            grd[lvl].kvoli();
            return(block::stop);
      }
   }
#endif  
}    

void r_mesh::moveboundaries() {
   
   /* MOVE BOUNDARY POSITIONS */
   for(int i=0;i<nsbd;++i)
      r_sbdry(i)->tadvance();
   
   return;
}

block::ctrl r_mesh::rsdl(int excpt) {

   switch (excpt) {
#ifdef FOURTH
      case(0):
         mp_phase = -1;
         rsdl1();
         return(block::advance);
      case(1):
         ++mp_phase;
         /* MESSAGE PASSING SEQUENCE */
         switch(mp_phase%3) {
            case(0):
               msgload(mp_phase/3,(FLT *) fscr2,0,1,2);
               return(block::stay);
            case(1):
               msgpass(mp_phase/3);
               return(block::stay);
            case(2):
               return(static_cast<block::ctrl>(msgwait_rcv(mp_phase/3,(FLT *) fscr2,0,1,2));
         }
#endif
      case(0+P2):
         mp_phase = -1;
         rsdl();
         return(block::advance);

         
      case(1+P2):
         ++mp_phase;
         /* MESSAGE PASSING SEQUENCE */
         switch(mp_phase%3) {
            case(0):
               msgload(mp_phase/3,(FLT *) fscr3,0,1,2);
               return(block::stay);
            case(1):
               msgpass(mp_phase/3);
               return(block::stay);
            case(2):
               return(static_cast<block::ctrl>(msgwait_rcv(mp_phase/3,(FLT *) fscr3,0,1,2)));
         }
         
      case(2+P2):
         return(block::stop);
         
      default:
         *log << "flow control error, rsdl" << std::endl;
         exit(1);
   }
   
   return(block::stop);
}


block::ctrl r_mesh::setup_preconditioner(int excpt) {
   switch (excpt) {
#ifdef FOURTH
      case(0):
         mp_phase = -1;
         vddt1();
         return(block::advance);
      case(1):
         ++mp_phase;
         /* MESSAGE PASSING SEQUENCE */
         switch(mp_phase%3) {
            case(0):
               msgload(mp_phase/3,fscr1,0,0,1);
               return(block::stay);
            case(1):
               msgpass(mp_phase/3);
               return(block::stay);
            case(2):
               return(static_cast<block::ctrl>(msgwait_rcv(mp_phase/3,fscr1,0,0,1));
         }
#endif
      case(0+P2):
         mp_phase = -1;
         vddt();
         return(block::advance);
         
      case(1+P2):
         ++mp_phase;
         /* MESSAGE PASSING SEQUENCE */
         switch(mp_phase%3) {
            case(0):
               msgload(mp_phase/3,fscr1.data(),0,0,1);
               return(block::stay);
            case(1):
               msgpass(mp_phase/3);
               return(block::stay);
            case(2):
               return(static_cast<block::ctrl>(msgwait_rcv(mp_phase/3,fscr1.data(),0,0,1)));
         }
         
      case(2+P2):
         vddti();
         return(block::stop);
   }
   
   *log << "flow control error: vddt" << std::endl;
   exit(1);
   
   return(block::stop);
}
