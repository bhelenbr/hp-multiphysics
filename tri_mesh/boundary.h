#ifndef _boundary_h_
#define _boundary_h_

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
#include <symbolic_function.h>
#include <float.h>
#include "blocks.h"

using namespace blitz;

#ifdef MPISRC
#include <mpi.h>
#endif

#ifdef SINGLE
#define FLT float
#define EPSILON FLT_EPSILON
#else
#ifndef FLT
#define FLT double
#define EPSILON DBL_EPSILON
#endif
#endif

//#define MPDEBUG

/** \brief Template class to make a communciation boundary
 *
 * \ingroup boundary
 * Contains variables and routines for communciation boundaries.
 * Can be applied to vrtx_bdry or edge_bdry to make specific type
 */
template<class BASE,class MESH> class comm_bdry : public BASE {
	protected:
		static const int maxmatch = 8;
		bool first; //!< For master-slave communication. Only one boundary of matching boundaries is first
		int maxgroup;
		int groupmask;  //!< To make groups that only communicate in restricted situations group 0 all, group 1 all phased, group 2 partitions, group 3 manifolds
		Array<int,1> maxphase; //!<  For phased symmetric message passing for each group
		Array<TinyVector<int,maxmatch>,1> phase;  //!< To set-up staggered sequence of symmetric passes for each group (-1 means skip)
		int& msg_phase(int grp, int match) {return(phase(grp)(match));} //!< virtual accessor
		int buffsize; //!< Size in bytes of buffer
		void *sndbuf; //!< Raw memory for outgoing message buffer
		Array<FLT,1> fsndbufarray; //!< Access to outgoing message buffer for floats
		Array<int,1> isndbufarray; //!< Access to outgoing message buffer for ints
		int msgsize; //!< Outgoing size
		boundary::msg_type msgtype; //!< Outgoing type

		/** Different types of matching boundaries,
		* local is same processor same thread
		* mpi is different processor
		* someday have threads but not yet
		*/
#ifdef MPISRC
		enum matchtype {local, mpi};
#else
		enum matchtype {local};
#endif
		int nmatch; //!< Number of matching boundaries
		int& nmatches() {return(nmatch);} //!< virtual accessor
		TinyVector<matchtype,maxmatch> mtype; //!< Local or mpi or ?
		TinyVector<boundary *,maxmatch> local_match; //!< Pointers to local matches
		TinyVector<int,maxmatch> snd_tags; //!< Identifies each connection uniquely
		TinyVector<int,maxmatch> rcv_tags; //!< Identifies each connection uniquely
		TinyVector<void *,maxmatch> rcvbuf; //!< Raw memory to store incoming messages
		TinyVector<Array<FLT,1>,maxmatch> frcvbufarray; //!< Access to incoming message buffer for floats
		TinyVector<Array<int,1>,maxmatch> ircvbufarray; //!< Access to incoming message buffer for ints

#ifdef MPISRC
		TinyVector<int,maxmatch> mpi_match; //!< Processor numbers for mpi
		TinyVector<MPI_Request,maxmatch> mpi_rcvrqst; //!< Identifier returned from mpi to monitor success of recv
		TinyVector<MPI_Request,maxmatch> mpi_sndrqst; //!< Identifier returned from mpi to monitor success of send
#endif

	public:
		comm_bdry(int inid, MESH &xin) : BASE(inid,xin), first(1), maxgroup(1), groupmask(0x3), buffsize(0), nmatch(0) {
			maxphase.resize(maxgroup+1);
			phase.resize(maxgroup+1);
			for(int m=0;m<maxmatch;++m) phase(0)(m) = 0;
		}
		comm_bdry(const comm_bdry<BASE,MESH> &inbdry, MESH& xin) : BASE(inbdry,xin), first(inbdry.first), maxgroup(inbdry.maxgroup), groupmask(inbdry.groupmask), buffsize(0), nmatch(0) {
			maxphase.resize(maxgroup+1);
			phase.resize(maxgroup+1);
			maxphase = inbdry.maxphase;
			for(int k=0;k<maxgroup+1;++k)
				phase(k) = inbdry.phase(k);

			/* COPY THESE, BUT WILL HAVE TO BE RESET TO NEW MATCHING SIDE */
			first = true; // Findmatch sets this
			mtype = inbdry.mtype;
			local_match = local_match;
			snd_tags = inbdry.snd_tags;
			rcv_tags = inbdry.rcv_tags;

#ifdef MPISRC
			mpi_match = inbdry.mpi_match;
#endif
			return;
		}

		comm_bdry<BASE,MESH>* create(MESH &xin) const {return(new comm_bdry<BASE,MESH>(*this,xin));}
		bool is_comm() {return(true);}
		bool& is_frst() {return(first);}
		int& group() {return(groupmask);}
		bool in_group(int grp) {return(((1<<grp)&groupmask));}
		int& sndsize() {return(msgsize);}
		boundary::msg_type& sndtype() {return(msgtype);}
		int& matchphase(boundary::groups group, int matchnum) {return(phase(group)(matchnum));}
		int& isndbuf(int indx) {return(isndbufarray(indx));}
		FLT& fsndbuf(int indx) {return(fsndbufarray(indx));}
		int& ircvbuf(int m,int indx) {return(ircvbufarray(m)(indx));}
		FLT& frcvbuf(int m,int indx) {return(frcvbufarray(m)(indx));}

		void resize_buffers(int nfloats) {
			if (buffsize) free(sndbuf);
			buffsize = nfloats*sizeof(FLT);
			sndbuf = malloc(buffsize);
			Array<FLT,1> temp(static_cast<FLT *>(sndbuf), buffsize/sizeof(FLT), neverDeleteData);
			fsndbufarray.reference(temp);
			Array<int,1> temp1(static_cast<int *>(sndbuf), buffsize/sizeof(int), neverDeleteData);
			isndbufarray.reference(temp1);
			
			for (int m=0;m<nmatch;++m) {
				if (rcvbuf(m)) free(rcvbuf(m));
				rcvbuf(m) = new FLT[buffsize/sizeof(FLT)]; // xmalloc(buffsize);
				Array<FLT,1> temp(static_cast<FLT *>(rcvbuf(m)), buffsize/sizeof(FLT), neverDeleteData);
				frcvbufarray(m).reference(temp);
				Array<int,1> temp1(static_cast<int *>(rcvbuf(m)), buffsize/sizeof(int), neverDeleteData);
				ircvbufarray(m).reference(temp1);
			}
		}

		void alloc(int nels) {
			BASE::alloc(nels);
			resize_buffers(nels*3);
		}

		void init(input_map& inmap) {
			int j,k,m;
			std::string keyword,val;
			std::map<std::string,std::string>::const_iterator mi;
			std::istringstream data;
			std::ostringstream nstr;

			BASE::init(inmap);

			keyword = BASE::idprefix +"_first";
			inmap.getwdefault(keyword,first,true);

			/* SET GROUP MEMBERSHIP FLAGS */
			maxgroup = 0;
			groupmask = 0;
			inmap.getlinewdefault(BASE::idprefix + "_group",val,"0 1"); // DEFAULT IS FIRST 2 GROUPS
			data.str(val);
			while(data >> m) {
				groupmask = groupmask|(1<<m);
				maxgroup = MAX(maxgroup,m);
			}
			data.clear();

			/* LOAD PHASES */
			maxphase.resize(maxgroup+1);
			maxphase = 0;
			phase.resize(maxgroup+1);
			phase = 0;
			/* SKIP GROUP 0 BECAUSE THAT GROUP IS NOT PHASED */
			for(k=1;k<maxgroup+1;++k) {
				if (!(groupmask&(1<<k))) continue;

				nstr.str("");
				nstr << BASE::idprefix << "_phase" << k << std::flush;
				if (inmap.getline(nstr.str(),val)) {
					data.str(val);
					m = 0;
					while(data >> j) {
						phase(k)(m) = j;
						maxphase(k) = MAX(maxphase(k),phase(k)(m));
						++m;
					}
					data.clear();
				}
			}
		}

		void output(std::ostream& fout) {
			BASE::output(fout);

			fout << BASE::idprefix << "_group" << ": ";
			for(int k=0;k<maxgroup+1;++k)
				if (groupmask&(1<<k)) fout << k << ' ';
			fout << std::endl;

			for(int k=0;k<maxgroup+1;++k) {
				if (groupmask&(1<<k)) {
					fout << BASE::idprefix << "_phase (not set yet so this is dumb)" << k << ": ";
					for (int m=0;m<nmatch;++m)
						fout << phase(k)(m) << " ";
					fout << std::endl;
				}
			}
		}

		int local_cnnct(boundary *bin, int snd_tag, int rcv_tag) {
			if (bin->idnum == BASE::idnum) {
				mtype(nmatch) = local;
				local_match(nmatch) = bin;
				snd_tags(nmatch) = snd_tag;
				rcv_tags(nmatch) = rcv_tag;
				rcvbuf(nmatch) = new FLT[buffsize/sizeof(FLT)]; // xmalloc(buffsize);
				Array<FLT,1> temp(static_cast<FLT *>(rcvbuf(nmatch)), buffsize/sizeof(FLT), neverDeleteData);
				frcvbufarray(nmatch).reference(temp);
				Array<int,1> temp1(static_cast<int *>(rcvbuf(nmatch)), buffsize/sizeof(int), neverDeleteData);
				ircvbufarray(nmatch).reference(temp1);
				++nmatch;
				return(0);
			}
			*BASE::x.gbl->log << "error: not local match" << BASE::idnum << bin->idnum << std::endl;
			return(1);
		}

#ifdef MPISRC
		int mpi_cnnct(int nproc, int snd_tag, int rcv_tag) {
			mtype(nmatch) = mpi;
			mpi_match(nmatch) = nproc;
			snd_tags(nmatch) = snd_tag;
			rcv_tags(nmatch) = rcv_tag;
			rcvbuf(nmatch) = new FLT[buffsize/sizeof(FLT)]; // xmalloc(buffsize);
			Array<FLT,1> temp(static_cast<FLT *>(rcvbuf(nmatch)), buffsize/sizeof(FLT), neverDeleteData);
			frcvbufarray(nmatch).reference(temp);
			Array<int,1> temp1(static_cast<int *>(rcvbuf(nmatch)), buffsize/sizeof(int), neverDeleteData);
			ircvbufarray(nmatch).reference(temp1);
			++nmatch;
			return(0);
		}
#endif

		/* MECHANISM FOR SYMMETRIC SENDING */
		void comm_prepare(boundary::groups grp, int phi, boundary::comm_type type) {
			int m;
			int nrecvs_to_post;
			int nsends_to_post;

#ifdef MPISRC
			int err;
#endif

			if (!in_group(grp)) return;

			/* SWITCHES FOR MASTER_SLAVE */
			switch(type) {
				case(boundary::master_slave): {
					if (first) {
						nsends_to_post = nmatch;
						nrecvs_to_post = 0;
					}
					else {
						nrecvs_to_post = 1;
						nsends_to_post = 0;
					}
					break;
				}
				case(boundary::slave_master): {
					if (!first) {
						nsends_to_post = 1;
						nrecvs_to_post = 0;
					}
					else {
						nrecvs_to_post = nmatch;
						nsends_to_post = 0;
					}
					break;
				}
				default: { // SYMMETRIC
					nrecvs_to_post = nmatch;
					nsends_to_post = nmatch;
					break;
				}
			}

#ifdef MPDEBUG
			if (nsends_to_post) {
				*BASE::x.gbl->log << "preparing to send these messages from "  << BASE::idprefix << "with type " << type << std::endl;
				switch(sndtype()) {
					case(boundary::flt_msg): {
						*BASE::x.gbl->log << fsndbufarray(Range(0,sndsize()-1)) << std::endl;
						break;
					}
					case(boundary::int_msg): {
						*BASE::x.gbl->log << isndbufarray(Range(0,sndsize()-1)) << std::endl;

						break;
					}
				}
			}
#endif

			/* MPI POST RECEIVES FIRST */
			for(m=0;m<nrecvs_to_post;++m) {
				if (phi != phase(grp)(m)) continue;

				switch(mtype(m)) {
					case(local):
						/* NOTHING TO DO FOR LOCAL RECEIVES */
						break;
#ifdef MPISRC
					case(mpi):
						switch(sndtype()) {
							case(boundary::flt_msg): {
#ifdef SINGLE
								err = MPI_Irecv(&frcvbuf(m,0), buffsize/sizeof(FLT), MPI_FLOAT,
									mpi_match(m), rcv_tags(m), MPI_COMM_WORLD, &mpi_rcvrqst(m));
#else
								err = MPI_Irecv(&frcvbuf(m,0), buffsize/sizeof(FLT), MPI_DOUBLE,
									mpi_match(m), rcv_tags(m), MPI_COMM_WORLD, &mpi_rcvrqst(m));
#endif
								break;
							}
							case(boundary::int_msg): {
								err = MPI_Irecv(&ircvbuf(m,0), buffsize/sizeof(int), MPI_INT,
									mpi_match(m), rcv_tags(m), MPI_COMM_WORLD, &mpi_rcvrqst(m));
								break;
							}
						}
#endif
				}
			}

			/* LOCAL POST SENDS FIRST */
			for(m=0;m<nsends_to_post;++m) {
				if (phi != phase(grp)(m)) continue;

				switch(mtype(m)) {
					case(local):
						sim::blks.notify_change(snd_tags(m),true);
						break;
#ifdef MPISRC
					case(mpi):
						/* NOTHING TO DO FOR MPI SENDS */
						break;
#endif
				}
			}
		}

		void comm_exchange(boundary::groups grp, int phi, boundary::comm_type type) {
			int i,m;
#ifdef MPISRC
			int err;
#endif
			int nrecvs_to_post = nmatch, nsends_to_post = nmatch;

			if (!in_group(grp)) return;

			switch(type) {
				case(boundary::master_slave): {
					if (!first) {
						nrecvs_to_post = 1;
						nsends_to_post = 0;
					}
					else {
						nrecvs_to_post = 0;
						nsends_to_post = nmatch;
					}
					break;
				}
				case(boundary::slave_master): {
					if (first) {
						nrecvs_to_post = nmatch;
						nsends_to_post = 0;
					}
					else {
						nrecvs_to_post = 0;
						nsends_to_post = 1;
					}
					break;
				}
				case(boundary::symmetric): {
					nrecvs_to_post = nmatch;
					nsends_to_post = nmatch;
					break;
				}
			}


			/* LOCAL PASSES */
			for(m=0;m<nrecvs_to_post;++m) {
				if (phi != phase(grp)(m) || mtype(m) != local) continue;

				sim::blks.waitforslot(rcv_tags(m),true);

				switch(sndtype()) {
					case(boundary::flt_msg): {
						for(i=0;i<local_match(m)->sndsize();++i)
							frcvbuf(m,i) = local_match(m)->fsndbuf(i);
						break;
					}
					case(boundary::int_msg): {
						for(i=0;i<local_match(m)->sndsize();++i)
							ircvbuf(m,i) = local_match(m)->isndbuf(i);
						break;
					}
				}
				sim::blks.notify_change(rcv_tags(m),false);
			}

#ifdef MPISRC
			/* MPI PASSES */
			for(m=0;m<nsends_to_post;++m) {
				if (phi != phase(grp)(m) || mtype(m) != mpi) continue;

				switch(sndtype()) {
					case(boundary::flt_msg): {
#ifdef SINGLE
						err = MPI_Isend(&fsndbufarray(0), msgsize, MPI_FLOAT,
							mpi_match(m), snd_tags(m), MPI_COMM_WORLD, &mpi_sndrqst(m));
#else
						err = MPI_Isend(&fsndbufarray(0), msgsize, MPI_DOUBLE,
							mpi_match(m), snd_tags(m), MPI_COMM_WORLD, &mpi_sndrqst(m));
#endif
						break;
					}
					case(boundary::int_msg): {
						err = MPI_Isend(&isndbufarray(0), msgsize, MPI_INT,
							mpi_match(m), snd_tags(m), MPI_COMM_WORLD, &mpi_sndrqst(m));
						break;
					}
				}
			}
#endif
			return;
		}

		int comm_wait(boundary::groups grp, int phi, boundary::comm_type type) {
			int nrecvs_to_post;
			int nsends_to_post;

			if (!in_group(grp)) return(1);

			switch(type) {
				case(boundary::master_slave): {
					if (first) {
						nrecvs_to_post = 0;
						nsends_to_post = nmatch;
					}
					else {
						nrecvs_to_post = 1;
						nsends_to_post = 0;
					}
					break;
				}

				case(boundary::slave_master): {
					if (!first) {
						nrecvs_to_post = 0;
						nsends_to_post = 1;
					}
					else {
						nrecvs_to_post = nmatch;
						nsends_to_post = 0;
					}
					break;
				}

				default: {  // SYMMETRIC
					nrecvs_to_post = nmatch;
					nsends_to_post = nmatch;
					break;
				}
			}

			for(int m=0;m<nsends_to_post;++m) {
				if (phi != phase(grp)(m)) continue;

				switch(mtype(m)) {
					case(local): {
						sim::blks.waitforslot(snd_tags(m),false);
						break;
					}
#ifdef MPISRC
					case(mpi): {
						MPI_Status status;
						MPI_Wait(&mpi_sndrqst(m), &status);
						break;
					}
#endif
				}
			}


			for(int m=0;m<nrecvs_to_post;++m) {
				if (phi != phase(grp)(m)) continue;

				switch(mtype(m)) {
					case(local): {
						break;
					}
#ifdef MPISRC
					case(mpi): {
						MPI_Status status;
						MPI_Wait(&mpi_rcvrqst(m), &status);
						break;
					}
#endif
				}
			}

			/* ONE MEANS FINISHED 0 MEANS MORE TO DO */
			return((phi-maxphase(grp) >= 0 ? 1 : 0));
		}

		bool comm_finish(boundary::groups grp, int phi, boundary::comm_type type, boundary::operation op) {

			if (!in_group(grp)) return(false);

			switch(type) {
				case(boundary::slave_master): {
					if (!is_frst()) return(false);
				}

				case(boundary::master_slave): {
					if (is_frst() || (phase(grp)(0) != phi)) return(false);

#ifdef MPDEBUG
					*BASE::x.gbl->log << "finish master_slave"  << BASE::idprefix << std::endl;
#endif
					switch(sndtype()) {
						case(boundary::int_msg):
							for(int j=0;j<sndsize();++j) {
								isndbufarray(j) = ircvbuf(0,j);
							}
#ifdef MPDEBUG
							*BASE::x.gbl->log << isndbufarray(Range(0,sndsize()-1)) << std::endl;
#endif
							break;

						case(boundary::flt_msg):
							for(int j=0;j<sndsize();++j) {
								fsndbufarray(j) = frcvbuf(0,j);
							}
#ifdef MPDEBUG
							*BASE::x.gbl->log << fsndbufarray(Range(0,sndsize()-1)) << std::endl;
#endif
							break;
					}
					return(true);

				}

				default: {
					switch(sndtype()) {
						case(boundary::flt_msg): {
							switch(op) {
								case(boundary::average): {
									int matches = 1;
									for(int m=0;m<nmatch;++m) {
										if (phase(grp)(m) != phi) continue;

										++matches;
										for(int j=0;j<sndsize();++j) {
											fsndbufarray(j) += frcvbuf(m,j);
										}
									}
									if (matches > 1 ) {
										FLT mtchinv = 1./matches;
										for(int j=0;j<sndsize();++j)
											fsndbufarray(j) *= mtchinv;
#ifdef MPDEBUG
										*BASE::x.gbl->log << "finish average"  << BASE::idprefix << std::endl;
										*BASE::x.gbl->log << fsndbufarray(Range(0,sndsize()-1)) << std::endl;
#endif
									}
									return(true);
								}

								case(boundary::sum): {
									int matches = 1;

									for(int m=0;m<nmatch;++m) {
										if (phase(grp)(m) != phi) continue;

										++matches;
										for(int j=0;j<sndsize();++j) {
											fsndbufarray(j) += frcvbuf(m,j);
										}
									}
									if (matches > 1 ) {
#ifdef MPDEBUG
										*BASE::x.gbl->log << "finish sum"  << BASE::idprefix << std::endl;
										*BASE::x.gbl->log << fsndbufarray(Range(0,sndsize()-1)) << std::endl;
#endif
										return(true);

									}
									return(false);
								}

								case(boundary::maximum): {
									int matches = 1;

									for(int m=0;m<nmatch;++m) {
										if (phase(grp)(m) != phi) continue;

										++matches;

										for(int j=0;j<sndsize();++j) {
											fsndbufarray(j) = MAX(fsndbufarray(j),frcvbuf(m,j));
										}
									}

									if (matches > 1 ) {
#ifdef MPDEBUG
										*BASE::x.gbl->log << "finish max"  << BASE::idprefix << std::endl;
										*BASE::x.gbl->log << fsndbufarray(Range(0,sndsize()-1)) << std::endl;
#endif
										return(true);
									}
									return(false);
								}

								case(boundary::minimum): {
									int matches = 1;

									for(int m=0;m<nmatch;++m) {
										if (phase(grp)(m) != phi) continue;

										++matches;

										for(int j=0;j<sndsize();++j) {
											fsndbufarray(j) = MIN(fsndbufarray(j),frcvbuf(m,j));
										}
									}

									if (matches > 1 ) {
#ifdef MPDEBUG

										*BASE::x.gbl->log << "finish min"  << BASE::idprefix << std::endl;
										*BASE::x.gbl->log << fsndbufarray(Range(0,sndsize()-1)) << std::endl;
#endif
										return(true);

									}
									return(false);
								}

								default: {
									*BASE::x.gbl->log << "replacement with symmetric sending?" << std::endl;
									exit(1);
								}
							}
						}

						case(boundary::int_msg): {
							switch(op) {
								case(boundary::average): {
									*BASE::x.gbl->log << "averaging with integer messages?" << std::endl;
									exit(1);
								}

								case(boundary::sum): {
									int matches = 1;

									for(int m=0;m<nmatch;++m) {
										if (phase(grp)(m) != phi) continue;

										++matches;
										for(int j=0;j<sndsize();++j) {
											isndbufarray(j) += ircvbuf(m,j);
										}
									}
									if (matches > 1 ) {
#ifdef MPDEBUG

										*BASE::x.gbl->log << "finish sum"  << BASE::idprefix << std::endl;
										*BASE::x.gbl->log << isndbufarray(Range(0,sndsize()-1)) << std::endl;
#endif
										return(true);
									}
									return(false);
								}

								case(boundary::maximum): {
									int matches = 1;

									for(int m=0;m<nmatch;++m) {
										if (phase(grp)(m) != phi) continue;

										++matches;

										for(int j=0;j<sndsize();++j) {
											isndbufarray(j) = MAX(isndbufarray(j),ircvbuf(m,j));
										}
									}

									if (matches > 1 ) {
#ifdef MPDEBUG

										*BASE::x.gbl->log << "finish max"  << BASE::idprefix << std::endl;
										*BASE::x.gbl->log << isndbufarray(Range(0,sndsize()-1)) << std::endl;
#endif
										return(true);
									}
									return(false);
								}

								case(boundary::minimum): {
									int matches = 1;

									for(int m=0;m<nmatch;++m) {
										if (phase(grp)(m) != phi) continue;

										++matches;

										for(int j=0;j<sndsize();++j) {
											isndbufarray(j) = MIN(isndbufarray(j),ircvbuf(m,j));
										}
									}

									if (matches > 1 ) {
#ifdef MPDEBUG

										*BASE::x.gbl->log << "finish min"  << BASE::idprefix << std::endl;
										*BASE::x.gbl->log << isndbufarray(Range(0,sndsize()-1)) << std::endl;
#endif
										return(true);
									}
									return(false);
								}

								default: {
									*BASE::x.gbl->log << "replacement with symmetric sending?" << std::endl;
									exit(1);
								}
							}
						}
					}
				}
			}

			return(false);
		}
};



/* GENERIC TEMPLATE FOR SHAPES */
template<int ND> class geometry {
	protected:
		virtual FLT hgt(TinyVector<FLT,ND> pt, FLT time = 0.0) {return(0.0);}
		virtual FLT dhgt(int dir, TinyVector<FLT,ND> pt, FLT time = 0.0) {return(1.0);}

	public:
		geometry* create() const {return(new geometry);}
		virtual void mvpttobdry(TinyVector<FLT,ND> &pt, FLT time) {
			int iter,n;
			FLT mag, delt_dist;

			/* FOR AN ANALYTIC SURFACE */
			iter = 0;
			do {
				mag = 0.0;
				for(n=0;n<ND;++n)
					mag += pow(dhgt(n,pt.data(),time),2);
				mag = sqrt(mag);
				delt_dist = -hgt(pt.data(),time)/mag;
				for(n=0;n<ND;++n)
					pt(n) += delt_dist*dhgt(n,pt.data(),time)/mag;
				if (++iter > 100) {
					std::cout << "curved iterations exceeded curved boundary " << pt(0) << ' ' << pt(1) << '\n';  // FIXME:  NEED TO FIX
					exit(1);
				}
			} while (fabs(delt_dist) > 10.*EPSILON);

			return;
		}
		virtual void init(input_map& inmap, std::string idprefix, std::ostream& log) {}
		virtual void output(std::ostream& fout, std::string idprefix) {}
		virtual ~geometry() {}
};


template<int ND> class symbolic_point : public geometry<ND> {
	protected:
		TinyVector<symbolic_function<0>,2> loc;
	public:
		symbolic_point() : geometry<ND>() {}
		symbolic_point(const symbolic_point& tgt) : geometry<ND>(tgt), loc(tgt.loc) {}

		void init(input_map& inmap, std::string idprefix, std::ostream& log) {
			geometry<ND>::init(inmap,idprefix,log);

			if (inmap.find(idprefix +"_locx0") != inmap.end()) {
				loc(0).init(inmap,idprefix+"_locx0");
			}
			else {
				log << "couldn't find shape function " << idprefix << "_locx0" << std::endl;
				exit(1);
			}

			if (inmap.find(idprefix +"_locx1") != inmap.end()) {
				loc(1).init(inmap,idprefix+"_locx1");
			}
			else {
				log << "couldn't find shape function " << idprefix << "_locx1" << std::endl;
				exit(1);
			}
		}

		void mvpttobdry(TinyVector<FLT,ND> &pt, FLT time) {
			pt(0) = loc(0).Eval(time);
			pt(1) = loc(1).Eval(time);
		}
};

template<int ND> class symbolic_shape : public geometry<ND> {
	protected:
		symbolic_function<ND> h, dhdx0, dhdx1;
		FLT hgt(TinyVector<FLT,ND> pt, FLT time = 0.0) {
			return(h.Eval(pt,time));
		}
		FLT dhgt(int dir, TinyVector<FLT,ND> pt, FLT time = 0.0) {
			if (dir) return(dhdx1.Eval(pt,time));
			return(dhdx0.Eval(pt,time));
		}
	public:
		symbolic_shape() : geometry<ND>() {}
		symbolic_shape(const symbolic_shape& tgt) : geometry<ND>(tgt), h(tgt.h), dhdx0(tgt.dhdx0), dhdx1(tgt.dhdx1) {}
		void init(input_map& inmap, std::string idprefix, std::ostream& log) {
			geometry<ND>::init(inmap,idprefix,log);
			if (inmap.find(idprefix +"_h") != inmap.end()) {
				h.init(inmap,idprefix+"_h");
			}
			else {
				log << "couldn't find shape function " << idprefix << "_h" << std::endl;
				exit(1);
			}

			if (inmap.find(idprefix +"_dhdx0") != inmap.end()) {
				dhdx0.init(inmap,idprefix+"_dhdx0");
			}
			else {
				log << "couldn't find shape function for " << idprefix << "_dhdx0" << std::endl;
				exit(1);
			}

			if (inmap.find(idprefix +"_dhdx1") != inmap.end()) {
				dhdx1.init(inmap,idprefix+"_dhdx1");
			}
			else {
				log << "couldn't find shape function " << idprefix << "_dhdx1" << std::endl;
				exit(1);
			}
		}

};


class circle : public geometry<2> {
	public:
		FLT center[2];
		FLT radius;
		FLT hgt(TinyVector<FLT,2> pt,FLT time = 0.0) {
			return(radius*radius -pow(pt[0]-center[0],2) -pow(pt[1]-center[1],2));
		}
		FLT dhgt(int dir, TinyVector<FLT,2> pt,FLT time = 0.0) {
			return(-2.*(pt[dir]-center[dir]));
		}

		circle() : geometry<2>(), radius(0.5) {center[0] = 0.0; center[1] = 0.0;}
		circle(const circle &inbdry) : geometry<2>(inbdry), radius(inbdry.radius) {center[0] = inbdry.center[0]; center[1] = inbdry.center[1];}
		circle* create() const {return(new circle(*this));}

		void output(std::ostream& fout,std::string idprefix) {
			geometry<2>::output(fout,idprefix);
			fout << idprefix << "_center: " << center[0] << '\t' << center[1] << std::endl;
			fout << idprefix << "_radius: " << radius << std::endl;
		}

		void init(input_map& inmap,std::string idprefix, std::ostream& log) {
			geometry<2>::init(inmap,idprefix,log);

			FLT dflt[2] = {0.0, 0.0};
			inmap.getwdefault(idprefix+"_center",center,2,dflt);
			inmap.getwdefault(idprefix+"_radius",radius,0.5);
		}
};

class ellipse : public geometry<2> {
	public:
		TinyVector<FLT,2> axes;
		FLT hgt(TinyVector<FLT,2> pt, FLT time = 0.0) {
			return(1 -pow(pt[0]/axes(0),2) -pow(pt[1]/axes(1),2));
		}
		FLT dhgt(int dir, TinyVector<FLT,2> pt, FLT time = 0.0) {
			return(-2.*pt[dir]/pow(axes(dir),2));
		}

			public:
		ellipse() : geometry<2>() {}
		ellipse(const ellipse &inbdry) : geometry<2>(inbdry), axes(inbdry.axes) {}
		ellipse* create() const {return(new ellipse(*this));}

		void output(std::ostream& fout,std::string idprefix) {
			geometry<2>::output(fout,idprefix);
			fout << idprefix << "_a" << axes(0) << std::endl;
			fout << idprefix << "_b" << axes(1) << std::endl;
		}
		void init(input_map& inmap, std::string idprefix,std::ostream& log) {
			geometry<2>::init(inmap,idprefix,log);
			inmap.getwdefault(idprefix+"_a",axes(0),1.0);
			inmap.getwdefault(idprefix+"_b",axes(1),1.0);
		}
};


class naca : public geometry<2> {
	public:
		FLT sign;
		TinyVector<FLT,5> coeff;
		FLT scale;
		FLT theta;
		TinyVector<FLT,2> pos;

		FLT hgt(TinyVector<FLT,2> x, FLT time = 0.0) {
			TinyVector<FLT,2> pt;
			for(int n=0;n<2;++n)
				pt[n] = x[n] -pos(n);

			FLT temp = pt[0]*cos(theta) -pt[1]*sin(theta);
			pt[1] = pt[0]*sin(theta) +pt[1]*cos(theta);
			pt[0] = temp;
			pt *= scale;

			FLT poly = coeff[1]*pt[0] +coeff[2]*pow(pt[0],2) +coeff[3]*pow(pt[0],3) +coeff[4]*pow(pt[0],4) - sign*pt[1];
			return(coeff[0]*pt[0] -poly*poly/coeff[0]);
		}
		FLT dhgt(int dir, TinyVector<FLT,2> x, FLT time = 0.0) {
			TinyVector<FLT,2> pt;
			for(int n=0;n<2;++n)
				pt[n] = x[n] -pos(n);

			FLT temp = pt[0]*cos(theta) -pt[1]*sin(theta);
			pt[1] = pt[0]*sin(theta) +pt[1]*cos(theta);
			pt[0] = temp;
			pt *= scale;

			TinyVector<FLT,2> ddx;
			FLT poly = coeff[1]*pt[0] +coeff[2]*pow(pt[0],2) +coeff[3]*pow(pt[0],3) +coeff[4]*pow(pt[0],4) - sign*pt[1];
			FLT dpolydx = coeff[1] +2*coeff[2]*pt[0] +3*coeff[3]*pow(pt[0],2) +4*coeff[4]*pow(pt[0],3);
			FLT dpolydy = -sign;
			ddx(0) = coeff[0] -2*poly*dpolydx/coeff[0];
			ddx(1) = -2*poly*dpolydy/coeff[0];
			ddx *= scale;

			if (dir == 0) return(ddx(0)*cos(theta) +ddx(1)*sin(theta));
			return(ddx(0)*(-sin(theta)) +ddx(1)*cos(theta));
		}

		naca() : geometry<2>(), sign(1.0), scale(1.0), theta(0.0) {
			/* NACA 0012 is the default */
			sign = 1;
			coeff[0] = 1.4845; coeff[1] = -0.63; coeff[2] = -1.758; coeff[3] = 1.4215; coeff[4] = -0.5180;
			coeff *= 0.12;
			pos = 0.0;
		}
		naca(const naca &inbdry) : geometry<2>(inbdry), sign(inbdry.sign), scale(inbdry.scale), theta(inbdry.theta) {
			for(int i=0;i<5;++i)
				coeff[i] = inbdry.coeff[i];

			pos = inbdry.pos;
		}
		naca* create() const {return(new naca(*this));}

		void output(std::ostream& fout,std::string idprefix) {
			geometry<2>::output(fout,idprefix);
			fout << idprefix << "_sign: " << sign << std::endl;
			fout << idprefix << "_coeff: ";
			for(int i=0;i<5;++i)
				fout << coeff[i] << " ";
			fout << std::endl;
			fout << idprefix << "_scale: " << scale << std::endl;
			fout << idprefix << "_theta: " << theta << std::endl;
			fout << idprefix << "_center: " << pos << std::endl;
		}

		void init(input_map& inmap,std::string idprefix,std::ostream& log) {
			geometry<2>::init(inmap,idprefix,log);

			FLT thickness;
			inmap.getwdefault(idprefix+"_sign",sign,1.0);
			inmap.getwdefault(idprefix+"_thickness",thickness,0.12);
			inmap.getwdefault(idprefix+"_scale",scale,1.0);
			inmap.getwdefault(idprefix+"_theta",theta,0.0);
			scale = 1./scale;
			theta = theta*M_PI/180.0;

			std::string linebuf;
			istringstream datastream;
			FLT naca_0012_dflt[5] = {1.4845, -0.63, -1.758, 1.4215, -0.5180};
			inmap.getwdefault(idprefix+"_coeff",coeff.data(),5,naca_0012_dflt);
			coeff *= thickness;

			FLT ctr_dflt[2] = {0.0, 0.0};
			inmap.getwdefault(idprefix+"_center",pos.data(),2,ctr_dflt);
		}
};
#endif



