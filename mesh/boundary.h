/*
 *  boundary.h
 *  mesh
 *
 *  Created by Brian Helenbrook on Fri Jun 07 2002.
 *  Copyright (c) 2002 __MyCompanyName__. All rights reserved.
 *
 */
#include <utilities.h>
#include <stdio.h>
#include <input_map.h>

#ifdef MPISRC
#include <mpi.h>
#endif

#ifndef _boundary_h_
#define _boundary_h_

#ifdef SINGLE
#define FLT float
#define EPSILON FLT_EPSILON
#else
#ifndef FLT
#define FLT double
#define EPSILON DBL_EPSILON
#endif
#endif

/** \defgroup boundary Mesh Boundary Objects 
 *  This group contains object for dealing with mesh boundaries
 */
 
/** \brief Generic interface for a b.c. 
 *
 *  \ingroup boundary
 *  Contains all the basic functions for parallel communications
 */
class boundary {
    public:
        int idnum;
        std::string idprefix; 
        std::string mytype;

        boundary(int idin) : idnum(idin) {
            char buffer[100];
            std::string keyname;
            sprintf(buffer,"%d",idnum);
            idprefix = std::string(buffer);
            mytype = "boundary";
        }
        virtual void alloc(int n) {}
        virtual void output(std::ostream& fout) {
            fout << idprefix << "_type: " << mytype << std::endl;            
        }
        virtual void input(input_map& bdrydata) {}
        
        /* VIRTUAL FUNCTIONS FOR COMMUNICATION BOUNDARIES */
        enum msg_type {flt_msg, int_msg};
        enum groups {all,all_phased,partitions,manifolds};
        enum comm_type {symmetric,master_slave,slave_master};
        enum operation {average,sum,maximum,replace};
        union  {
            bool bdum;
            int idum;
            FLT fdum;
            msg_type mdum;
        } dummy;
        virtual bool is_comm() {return(false);}
        virtual bool& is_frst() {return(dummy.bdum=true);}
        virtual int& group() {return(dummy.idum=1);}
        virtual int matches() {return(0);}
        virtual int local_cnnct(boundary *bin, int msg_tag) {return 1;}
#ifdef MPISRC
        virtual int mpi_cnnct(int proc_tgt, int msg_tag) {return 1;}
#endif
        virtual int& matchphase(boundary::groups group, int matchnum) {return(dummy.idum=0);}
        virtual void resize_buffers(int size) {}
        virtual void *psndbuf() {return(&dummy);}
        virtual int& isndbuf(int indx) {return(dummy.idum);}
        virtual FLT& fsndbuf(int indx) {return(dummy.fdum);}
        virtual int& ircvbuf(int m,int indx) {return(dummy.idum);}
        virtual FLT& frcvbuf(int m,int indx) {return(dummy.fdum);}
        virtual int& sndsize() {return(dummy.idum=0);}
        virtual boundary::msg_type& sndtype() {return(dummy.mdum);}
        virtual void comm_prepare(boundary::groups group, int phase, comm_type type) {}
        virtual void comm_exchange(boundary::groups group, int phase, comm_type type) {}
        virtual int comm_wait(boundary::groups group, int phase, comm_type type) {return 1;}
        virtual int comm_nowait(boundary::groups group, int phase, comm_type type) {return 1;}
        virtual ~boundary() {}
    protected:
        int excpt;
};
#endif


