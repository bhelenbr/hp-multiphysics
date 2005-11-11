#include "tri_hp_cd.h"
#include "cd_bdry.h"
#include "hp_boundary.h"
#include <myblas.h>
#include <input_map.h>

/** \brief Helper object for vrtx_bdry 
 *
 * \ingroup boundary
 * Contains list of all vrtx_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tri_hp_cd_vtype {
   public:
      static const int ntypes = 1;
      enum ids {plain=1};
      const static char names[ntypes][40];
      static int getid(const char *nin) {
         for(int i=0;i<ntypes;++i) 
            if (!strcmp(nin,names[i])) return(i+1);
         return(-1);
      }
};

const char tri_hp_cd_vtype::names[ntypes][40] = {"plain"};

hp_vrtx_bdry* tri_hp_cd::getnewvrtxobject(int bnum, std::map<std::string,std::string> *bdrydata) {
   std::string keyword;
   std::istringstream data;
   std::map<std::string,std::string>::const_iterator mi;
   char idntystring[10];
   int type;        
   hp_vrtx_bdry *temp;  
   int idnum = vbdry(bnum)->idnum;
   
   type = idnum&0xffff;

   if (bdrydata) {
      sprintf(idntystring,"v%d",idnum);
      keyword = std::string(idntystring) + ".type";
      mi = (*bdrydata).find(keyword);
      if (mi != (*bdrydata).end()) {
         type = tri_hp_cd_vtype::getid((*mi).second.c_str());
         if (type < 0)  {
            *sim::log << "unknown vertex type:" << (*mi).second << std::endl;
            exit(1);
         }
      }
   }
   
   switch(type) {
      case tri_hp_cd_vtype::plain: {
         temp = new hp_vrtx_bdry(*this,*vbdry(bnum));
         break;
      }
      default: {
         std::cout << "unrecognizable vrtx type: " <<  type << " idnum: " << idnum << std::endl;
         temp = new hp_vrtx_bdry(*this,*vbdry(bnum));
         break;
      }
   } 
   
   if (bdrydata) temp->input(*bdrydata);
   
   return(temp);
}


/** \brief Helper object for side_bdry 
 *
 * \ingroup boundary
 * Contains list of all side_bdys's by name 
 * and has routine to return integer so can
 * allocate by name rather than by number
 */
class tri_hp_cd_stype {
   public:
      static const int ntypes = 2;
      enum ids {dirichlet=1,curved_dirichlet};
      static const char names[ntypes][40];
      static int getid(const char *nin) {
         for(int i=0;i<ntypes;++i)
            if (!strcmp(nin,names[i])) return(i+1);
         return(-1);
      }
};

const char tri_hp_cd_stype::names[ntypes][40] = {"dirichlet","curved_dirichlet"};

/* FUNCTION TO CREATE BOUNDARY OBJECTS */
hp_side_bdry* tri_hp_cd::getnewsideobject(int bnum, std::map<std::string,std::string> *bdrydata) {
   std::string keyword;
   std::istringstream data;
   std::map<std::string,std::string>::const_iterator mi;
   char idntystring[10];
   int type;        
   hp_side_bdry *temp;  
   int idnum = sbdry(bnum)->idnum;

   type = idnum&0xff;

   if (bdrydata) {
      sprintf(idntystring,"s%d",idnum);
      keyword = std::string(idntystring) + ".cd_type";
      mi = (*bdrydata).find(keyword);
      if (mi != (*bdrydata).end()) {
         type = tri_hp_cd_stype::getid((*mi).second.c_str());
         if (type < 0)  {
            *sim::log << "unknown side type:" << (*mi).second << std::endl;
            exit(1);
         }
      }
      else {
         *sim::log << "#couldn't find " << keyword << std::endl;
         type = tri_hp_cd_stype::dirichlet;
      }
   }

   switch(type) {
      case tri_hp_cd_stype::dirichlet: {
         temp = new dirichlet_cd<hp_side_bdry>(*this,*sbdry(bnum));
         break;
      }
      case tri_hp_cd_stype::curved_dirichlet: {
         temp = new dirichlet_cd<hp_curved>(*this,*sbdry(bnum));
         break;
      }
   }
   
   if (bdrydata) temp->input(*bdrydata);
   
   return(temp);
}


#ifdef SKIP


extern FLT df1d(int,FLT,FLT);
void chrctr(FLT ax, FLT ay, double wl[NV], double wr[NV], double norm[ND], double mv[ND]);

void hp_mgrid::setinflow() {
    int i,j,k,m,n,indx,v0,v1,info;
    FLT x,y;
    int sind;
   char uplo[] = "U";

   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type&INFL_MASK) {
         /* INFLOW BOUNDARIES */
         /* SET VERTEX VALUES OF U,V */   
       }
      
      if (sbdry[i].type&OUTF_MASK) {
         for(n=0;n<NV;++n)
            binfo[i][0].flx[n] = 0.0;
            
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(0);
            v1 = sd(sind).vrtx(1);
            
           
            for(n=0;n<ND;++n) {
               basis::tri(log2p).proj1d(vrtx(v0)(n),vrtx(v1)(n),&crd(n)(0,0));
               
               for(k=0;k<basis::tri(log2p).gpx;++k)
                  dcrd(n,0)(0,k) = 0.5*(vrtx(v1)(n)-vrtx(v0)(n));
            }
            
            /* NOW SET FLUXES */
            for(k=0;k<basis::tri(log2p).gpx;++k) 
               res(0)(0,k) = -hp_gbl->mu*RAD1D(k)*(df1d(0,crd(0)(0,k),crd(1)(0,k))*dcrd(1,0)(0,k));
            
            basis::tri(log2p).intgrt1d(&lf(0)(0),&res(0)(0,0));
            
            indx = j*(basis::tri(log2p).sm +1);
            binfo[i][indx++].flx[0] += lf(0)(0);
            for(m=0;m<basis::tri(log2p).sm;++m)
               binfo[i][indx++].flx[0] = lf(0)(m+2);
            binfo[i][indx].flx[0] = lf(0)(1);
         }
      }
   }
   
   return;
}

void hp_mgrid::addbflux(int mgrid) {
    int i,j,k,n,indx,indx1;
    int sind,v0,v1;
    FLT nrm[ND], wl[NV], wr[NV];
   FLT mvel[ND] = {0.0, 0.0};
   
   /***********************************/
   /* ADD SOURCE TERMS ON FINEST MESH */
   /***********************************/
   if(!mgrid) {
//      setinflow();  //TEMPORARY

      for(i=0;i<nsbd;++i) {         
         if (sbdry[i].type&OUTF_MASK) {
            /* ALLOWS FOR APPLIED STRESS ON BOUNDARY */
            indx = 0;
            for(j=0;j<sbdry(i)->nel;++j) {
               sind=sbdry(i)->el(j);
               v0 = sd(sind).vrtx(0);
               indx1 = sind*basis::tri(log2p).sm;
               for(n=0;n<NV;++n)
                  hp_gbl->res.v(v0,n) += binfo[i][indx].flx[n];
               ++indx;
               for(k=0;k<basis::tri(log2p).sm;++k) {
                  for(n=0;n<NV;++n)
                     hp_gbl->res.s(indx1)(n) += binfo[i][indx].flx[n];
                  ++indx;
                  ++indx1;
               }
            }
            v0 = sd(sind).vrtx(1);
            for(n=0;n<NV;++n)
               hp_gbl->res.v(v0,n) += binfo[i][indx].flx[n];
         }
      }
   }
   
   
   /* THESE ARE SOURCE TERMS WHICH CHANGE WITH THE SOLUTION */
   /* MUST BE UPDATED DURING MGRID FOR GOOD CONVERGENCE */
   for(i=0;i<nsbd;++i) {

      /* OUTFLOW BOUNDARY CONDITION    */
      if (sbdry[i].type&OUTF_MASK) {
         indx = 0;
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(0);
            v1 = sd(sind).vrtx(1);
            
            if (sbdry[i].type&CURV_MASK) {
               crdtocht1d(sind);
               for(n=0;n<ND;++n)
                  basis::tri(log2p).proj1d(&cht(n,0),&crd(n)(0,0),&dcrd(n,0)(0,0));
               
               crdtocht1d(sind,dvrtdt,hp_gbl->dbinfodt);
               for(n=0;n<ND;++n)
                  basis::tri(log2p).proj1d(&cht(n,0),&crd(n)(1,0));
            }
            else {
               for(n=0;n<ND;++n) {
                  basis::tri(log2p).proj1d(vrtx(v0)(n),vrtx(v1)(n),&crd(n)(0,0));
                  
                  for(k=0;k<basis::tri(log2p).gpx;++k)
                     dcrd(n,0)(0,k) = 0.5*(vrtx(v1)(n)-vrtx(v0)(n));
               
                  basis::tri(log2p).proj1d(dvrtdt[v0][n],dvrtdt[v1][n],&crd(n)(1,0));
               }
            }
            
            ugtouht1d(sind);
            for(n=0;n<NV;++n)
               basis::tri(log2p).proj1d(&uht(n)(0),&u(n)(0,0));
            
            for(k=0;k<basis::tri(log2p).gpx;++k) {
               for(n=0;n<NV;++n) {
                  wl[n] = u(n)(0,k);
                  wr[n] = (hp_gbl->initfunc)(n,crd(0)(0,k),crd(1)(0,k));
               }
               nrm[0] = dcrd(1,0)(0,k);
               nrm[1] = -dcrd(0,0)(0,k);

               for(n=0;n<ND;++n)
                  mvel[n] = sim::bd[0]*crd(n)(0,k) +crd(n)(1,k);
                  
               chrctr(axext,ayext,wl,wr,nrm,mvel);
                                 
               res(0)(0,k) = wl[0]*RAD1D(k)*((axext -mvel[0])*nrm[0] +(ayext -mvel[1])*nrm[1]);
            }
            
            for(n=0;n<NV;++n)
               basis::tri(log2p).intgrt1d(&lf(n)(0),&res(n)(0,0));
            
            for(n=0;n<NV;++n)
               hp_gbl->res.v(v0,n) += lf(n)(0);

            for(n=0;n<NV;++n)
               hp_gbl->res.v(v1,n) += lf(n)(1);
            
            indx1 = sind*basis::tri(log2p).sm;
            indx = 2;
            for(k=0;k<basis::tri(log2p).sm;++k) {
               for(n=0;n<NV;++n)
                  hp_gbl->res.s(indx1)(n) += lf(n)(indx);
               ++indx1;
               ++indx;
            }
         }
      }
   }
   
   return;
}

void hp_mgrid::bdry_vsnd() {
   int i,j,n,sind,count,v0,bnum;
   class mesh *tgt;
   
   /* SEND VERTEX INFO FOR Y_DIR*/
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & COMY_MASK) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         /* SEND VERTEX INFO */
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(0);
            for (n=0;n<NV;++n) 
               tgt->sbuff[bnum][count++] = hp_gbl->res.v(v0,n);
         }
         v0 = sd(sind).vrtx(1);
         for (n=0;n<NV;++n) 
            tgt->sbuff[bnum][count++] = hp_gbl->res.v(v0,n);
      }
   }
   
   return;
}

void hp_mgrid::bdry_mp() {
   int i,j,n,sind,count,v0,bnum;
   class mesh *tgt;
   
   /* THIS PART IS TO RECEIVE AND ZERO FOR VERTICES */
   /* RECEIVE VRTX MESSAGES */
   /* CALCULATE AVERAGE RESIDUAL */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & COMY_MASK) {
         count = 0;
         /* RECV VERTEX INFO */
         for(j=sbdry(i)->nel-1;j>=0;--j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(1);
            for(n=0;n<NV;++n)
               hp_gbl->res.v(v0,n) = 0.5*(hp_gbl->res.v(v0,n) +sbuff[i][count++]);
         }
         v0 = sd(sind).vrtx(0);
         for(n=0;n<NV;++n)
            hp_gbl->res.v(v0,n) = 0.5*(hp_gbl->res.v(v0,n) +sbuff[i][count++]);
      } 
   }
   
   /* SEND VERTEX INFO FOR X_DIR*/
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & COMX_MASK) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         /* SEND VERTEX INFO */
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(0);
            for (n=0;n<NV;++n) 
               tgt->sbuff[bnum][count++] = hp_gbl->res.v(v0,n);
         }
         v0 = sd(sind).vrtx(1);
         for (n=0;n<NV;++n) 
            tgt->sbuff[bnum][count++] = hp_gbl->res.v(v0,n);
      }
   }
   
   return;
}


void hp_mgrid::bdry_vrcvandzero() {
    int i,j,n;
    int sind,v0,count;
   
   /* THIS PART IS TO RECEIVE AND ZERO FOR VERTICES */
   /* RECEIVE VRTX MESSAGES */
   /* CALCULATE AVERAGE RESIDUAL */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & COMX_MASK) {
         count = 0;
         /* RECV VERTEX INFO */
         for(j=sbdry(i)->nel-1;j>=0;--j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(1);
            for(n=0;n<NV;++n)
               hp_gbl->res.v(v0,n) = 0.5*(hp_gbl->res.v(v0,n) +sbuff[i][count++]);
         }
         v0 = sd(sind).vrtx(0);
         for(n=0;n<NV;++n)
            hp_gbl->res.v(v0,n) = 0.5*(hp_gbl->res.v(v0,n) +sbuff[i][count++]);
      }         
   }

   /* APPLY VRTX DIRICHLET CONDITIONS TO RES */
   for(i=0;i<nvbd;++i) {
      if (vbdry[i].type&INFL_MASK) {
         for(j=0;j<vbdry[i].num;++j) {
            v0 = vbdry[i].el[j];
            for(n=0;n<NV;++n)
               hp_gbl->res.v(v0,n) = 0.0;
         }
      }
   }

   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type&INFL_MASK) {
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(0);
            for(n=0;n<NV;++n)
               hp_gbl->res.v(v0,n) = 0.0;
         }
         v0 = sd(sind).vrtx(1);
         for(n=0;n<NV;++n)
            hp_gbl->res.v(v0,n) = 0.0;
      }
   }
   
   return;
}

void hp_mgrid::bdry_ssnd(int mode) {
   int i,j,n,count,indx,bnum;
   class mesh *tgt;
   
   /* SEND SIDE INFO */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & (COMX_MASK +COMY_MASK)) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         for(j=0;j<sbdry(i)->nel;++j) {
            indx = sbdry(i)->el(j)*basis::tri(log2p).sm +mode;
            for (n=0;n<NV;++n) 
               tgt->sbuff[bnum][count++] = hp_gbl->res.s(indx)(n);
         }
      } 
   }
   
   return;
}

   
void hp_mgrid::bdry_srcvandzero(int mode) {
    int i,j,n;
    int sind,count,indx,sign;
   
   sign = (mode % 2 ? -1 : 1);
   
   /* THIS PART TO RECIEVE AND ZERO FOR SIDES */
   /* RECEIVE P'TH SIDE MODE MESSAGES */
   /* CALCULATE AVERAGE RESIDUAL */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & (COMX_MASK +COMY_MASK)) {
         count = 0;
         for(j=sbdry(i)->nel-1;j>=0;--j) {
            indx = sbdry(i)->el(j)*basis::tri(log2p).sm +mode;
            for(n=0;n<NV;++n)
               hp_gbl->res.s(indx)(n) = 0.5*(hp_gbl->res.s(indx)(n) +sign*sbuff[i][count++]);
         }
      }
   }

   /* APPLY SIDE DIRICHLET CONDITIONS TO MODE */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type&INFL_MASK) {
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j)*basis::tri(log2p).sm +mode;
            for(n=0;n<NV;++n)
               hp_gbl->res.s(sind)(n) = 0.0;
         }
      }
   }
   
   return;
}


void chrctr(FLT ax, FLT ay,double wl[NV], double wr[NV], double norm[ND], double mv[ND]) {
   FLT ul;
   FLT um,lam0,mag;
   
   /* CHARACTERISTIC FAR-FIELD B.C. */   
   mag = sqrt(norm[0]*norm[0] + norm[1]*norm[1]);
   
   norm[0] /= mag;
   norm[1] /= mag;
   
   ul =  ax*norm[0] +ay*norm[1];      
   um = mv[0]*norm[0] +mv[1]*norm[1];

   lam0 = ul-um;

   if (lam0 > 0.0)
      wl[0] = wl[0];
   else
      wl[0] = wr[0];

   /* SHOULDN'T CHANGE NORM */   
   norm[0] *= mag;
   norm[1] *= mag;
   
   return;
 
}

void hp_mgrid::bdrycheck1() {
   int i,j,n,v0,count,bnum,sind;
   class mesh *tgt;
   
   /* SEND SIDE INFO */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & (COMX_MASK +COMY_MASK +IFCE_MASK)) {
         bnum = sbdry[i].adjbnum;
         tgt = sbdry[i].adjmesh;
         count = 0;
         for(j=0;j<sbdry(i)->nel;++j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(0);
            for (n=0;n<ND;++n) 
               tgt->sbuff[bnum][count++] = vrtx(v0)(n);
         }
         v0 = sd(sind).vrtx(1);
         for (n=0;n<ND;++n) 
            tgt->sbuff[bnum][count++] = vrtx(v0)(n); 
      }    
   }
   
   return;
}

void hp_mgrid::bdrycheck2() {
    int i,j,n,v0;
    int sind,count;
      
   /* THIS PART TO RECIEVE AND ZERO FOR SIDES */
   /* RECEIVE P'TH SIDE MODE MESSAGES */
   /* CALCULATE AVERAGE RESIDUAL */
   for(i=0;i<nsbd;++i) {
      if (sbdry[i].type & PRDX_MASK) {
         count = 1;
         for(j=sbdry(i)->nel-1;j>=0;--j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(1);
            vrtx(v0)(1) = 0.5*(vrtx(v0)(1) +sbuff[i][count++]);
            ++count;
         }
         v0 = sd(sind).vrtx(0);
         vrtx(v0)(1) = 0.5*(vrtx(v0)(1) +sbuff[i][count++]); 
         continue;   
      }
      
      if (sbdry[i].type & PRDY_MASK) {
         count = 0;
         for(j=sbdry(i)->nel-1;j>=0;--j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(1);
            vrtx(v0)(0) = 0.5*(vrtx(v0)(0) +sbuff[i][count++]);
            ++count;
         }
         v0 = sd(sind).vrtx(0);
         vrtx(v0)(0) = 0.5*(vrtx(v0)(0) +sbuff[i][count++]); 
         continue;   
      }
      
      if (sbdry[i].type & (COMX_MASK +COMY_MASK +IFCE_MASK)) {
         count = 0;
         for(j=sbdry(i)->nel-1;j>=0;--j) {
            sind = sbdry(i)->el(j);
            v0 = sd(sind).vrtx(1);
            for(n=0;n<ND;++n)
               vrtx(v0)(n) = 0.5*(vrtx(v0)(n) +sbuff[i][count++]); 
         }
      }
   }
   
   return;
}

#endif


