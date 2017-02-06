//
//  hp_vrtx_bdry.cpp
//  tri_hp
//
//  Created by Brian Helenbrook on 4/25/16.
//
//

#include "hp_boundary.h"

void hp_vrtx_bdry::init(input_map& inmap,void* gbl_in) {
	std::string keyword,ibcname;
	
	/* FIND INITIAL CONDITION TYPE */
	keyword = base.idprefix + "_ibc";
	if (inmap.get(keyword,ibcname)) {
		ibc = x.getnewibc(ibcname);
		ibc->init(inmap,keyword);
	}
	
	keyword = base.idprefix + "_coupled";
	coupled = false;
	inmap.get(keyword,coupled);
	
	keyword = base.idprefix + "_frozen";
	frozen = false;
	inmap.get(keyword,frozen);
	
	keyword = base.idprefix +"_report";
	report_flag = false;
	inmap.get(keyword,report_flag);
	
	type.resize(x.NV,natural);
	Array<int,1> atemp(x.NV);
	if (inmap.get(base.idprefix+"_hp_typelist", atemp.data(), x.NV)) {
		for (int n=0;n<x.NV;++n) {
			type[n] = static_cast<bctypes>(atemp(n));
			if (type[n] == essential) {
				essential_indices.push_back(n);
			}
		}
	}
	
	std::string val;
	if (!inmap.getline(base.idprefix +"_c0_indices",val)) {
		/* Default is that all variables are continuous across boundary */
		for(int n=0;n<x.NV;++n) {
			c0_indices.push_back(n);
		}
	}
	else {
		istringstream data;
		data.str(val);
		int ind;
		while (data >> ind)
			c0_indices.push_back(ind);
		data.clear();
	}
	c0_indices_xy = c0_indices;
	if (x.mmovement == x.coupled_deformable) {
		c0_indices_xy.push_back(x.NV);
		c0_indices_xy.push_back(x.NV+1);
	}
	
#ifdef petsc
	base.resize_buffers((x.NV+x.ND)*60*(3 +3*x.sm0+x.im0));  // Allows for 10 elements of jacobian entries to be sent
#endif
}

/** This is to read solution data **/
void hp_vrtx_bdry::input(ifstream& fin,tri_hp::filetype typ,int tlvl) {
	std::string idin,mytypein;
	
	switch(typ) {
		case(tri_hp::text):
			fin >> idin >> mytypein;
			break;
		default:
			break;
	}
	return;
}

void hp_vrtx_bdry::output(const std::string& filename, tri_hp::filetype typ,int tlvl) {
	switch(typ) {
		case(tri_hp::text): {
			std::string fname;
			fname = filename +"_" +x.gbl->idprefix +".txt";
			ofstream fout;
			fout.open(fname.c_str(),std::ofstream::out | std::ofstream::app);
			fout << base.idprefix << " " << mytype << std::endl;
			fout.close();
			break;
		}
		case(tri_hp::tecplot): {
			if (report_flag) {
				streamsize oldprecision = (*x.gbl->log).precision(10);
				*x.gbl->log << base.idprefix << " position: " << x.pnts(base.pnt) << std::endl;
				*x.gbl->log << base.idprefix << " value: " << x.ug.v(base.pnt,Range::all()) << std::endl;
				(*x.gbl->log).precision(oldprecision);
			}
			break;
		}
		default:
			break;
	}
	return;
}


void hp_vrtx_bdry::rsdl(int stage) {
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	
	Array<FLT,1> lf(vdofs);
	element_rsdl(lf);
	
	for(int n=0;n<x.NV;++n)
		x.gbl->res.v(base.pnt,n) += lf(n);
	
}

void hp_vrtx_bdry::element_jacobian(Array<FLT,2>& K) {
	
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	Array<FLT,1> Rbar(vdofs),lf(vdofs);;
	
	/* Calculate and store initial residual */
	lf = 0.0;
	element_rsdl(lf);
	Rbar = lf;
	
	Array<FLT,1> dw(x.NV);
	dw = 0.0;
	for(int i=0;i<2;++i)
		for(int n=0;n<x.NV;++n)
			dw(n) = dw(n) + fabs(x.ug.v(base.pnt,n));
	
	dw *= eps_r;
	dw += eps_a;
	
	int kcol = 0;
	for(int var = 0; var < x.NV; ++var){
		x.ug.v(base.pnt,var) += dw(var);
		
		lf = 0.0;
		element_rsdl(lf);
		
		int krow = 0;
		for(int n=0;n<vdofs;++n)
			K(krow++,kcol) = (lf(n)-Rbar(n))/dw(var);
		
		++kcol;
		
		x.ug.v(base.pnt,var) -= dw(var);
	}
	
	int sind = x.ebdry(base.ebdry(1))->seg(0);
	FLT dx = eps_r*x.distance(x.seg(sind).pnt(0),x.seg(sind).pnt(1)) +eps_a;
	
	for(int var = 0; var < vdofs -x.NV; ++var){
		x.pnts(base.pnt)(var) += dx;
		
		lf = 0.0;
		element_rsdl(lf);
		
		int krow = 0;
		for(int n=0;n<vdofs;++n)
			K(krow++,kcol) = (lf(n)-Rbar(n))/dx;
		
		++kcol;
		
		x.pnts(base.pnt)(var) -= dx;
	}
}

#ifdef petsc
void hp_vrtx_bdry::non_sparse_snd(Array<int,1> &nnzero, Array<int,1> &nnzero_mpi) {
	
	if (!base.in_group(boundary::all_phased)) return;
	
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	
	/* Going to send all jacobian entries,  Diagonal entries for matching DOF's will be merged together not individual */
	/* Send number of non-zeros to matches */
	base.sndsize() = 0;
	base.sndtype() = boundary::int_msg;
	
	int pind = base.pnt*vdofs;
	for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
		base.isndbuf(base.sndsize()++) = nnzero(pind +*n);
	}
}


int hp_vrtx_bdry::non_sparse_rcv(int phase, Array<int,1> &nnzero, Array<int,1> &nnzero_mpi) {
	
	if (!base.in_group(boundary::all_phased) || base.matchphase(boundary::all_phased,0) != phase) return(0);
	
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	
	/* Reload to avoid overlap with sides */
	int pind = base.pnt*vdofs;
	int count = 0;
	for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
		nnzero(pind +*n) = base.isndbuf(count++);
		nnzero_mpi(pind +*n) = 0;
	}
	
	for(int m=0;m<base.nmatches();++m) {
		Array<int,1> target;
		if (base.is_local(m))  // only 1 matching boundary
			target.reference(nnzero);
		else
			target.reference(nnzero_mpi);
		
		int count = 0;
		int pind = base.pnt*vdofs;
		for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
			target(pind +*n) += base.ircvbuf(m,count++);
		}
	}
	return(count);
}

void hp_vrtx_bdry::petsc_jacobian() {
	
	/* Generic method for adding Jacobian terms */
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	Array<FLT,2> K(vdofs,vdofs);
	Array<int,1> rows(vdofs);
	Array<int,1> cols(vdofs);
	
	const int v0 = base.pnt;
	
	int rind = 0;
	int cind = 0;
	int gindx = vdofs*v0;
	for(int n=0;n<vdofs;++n) {
		rows(rind++) = gindx;
		cols(cind++) = gindx++;
	}
	
	element_jacobian(K);
	
#ifdef MY_SPARSE
	x.J.add_values(vdofs,rows,vdofs,cols,K);
#else
	MatSetValuesLocal(x.petsc_J,vdofs,rows.data(),vdofs,cols.data(),K.data(),ADD_VALUES);
#endif
}

void hp_vrtx_bdry::petsc_matchjacobian_snd() {
	
	if (!base.in_group(boundary::all_phased)) return;
	
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	
#ifdef MY_SPARSE
	/* I am cheating here and sending floats and int's together */
	/* Send Jacobian entries for continous variables  */
	base.sndsize() = 0;
	base.sndtype() = boundary::flt_msg;
	/* Send index of start of jacobian */
	base.fsndbuf(base.sndsize()++) = x.jacobian_start +0.1;
	/* Send index of start of surface unknowns (if they exist) */
	base.fsndbuf(base.sndsize()++) = jacobian_start+0.1;
	
	int rowbase = base.pnt*vdofs;
	/* Send continuous variables */
	for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
		int row = rowbase + *n;
		/* attach diagonal column # to allow continuity enforcement */
		base.fsndbuf(base.sndsize()++) = row +0.1;
		base.fsndbuf(base.sndsize()++) = x.J._cpt(row+1) -x.J._cpt(row) +0.1;
#ifdef MPDEBUG
		*x.gbl->log << "vertex sending " << x.J._cpt(row+1) -x.J._cpt(row) << " jacobian entries for vertex " << row/vdofs << " and variable " << *n << std::endl;
#endif
		for (int col=x.J._cpt(row);col<x.J._cpt(row+1);++col) {
#ifdef MPDEBUG
			*x.gbl->log << x.J._col(col) << ' ' << x.J._val(col) << ' ';
#endif
			base.fsndbuf(base.sndsize()++) = x.J._col(col) +0.1;
			base.fsndbuf(base.sndsize()++) = x.J._val(col);
		}
#ifdef MPDEBUG
		*x.gbl->log << std::endl;
#endif
	}
}

int hp_vrtx_bdry::petsc_matchjacobian_rcv(int phase)	{
	
	if (!base.in_group(boundary::all_phased) || base.matchphase(boundary::all_phased,0) != phase) return(0);
	
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	const int rowbase = base.pnt*vdofs;
	
	int count = 0;
	assert(x.jacobian_start == static_cast<int>(base.fsndbuf(count++)));
	assert(jacobian_start == static_cast<int>(base.fsndbuf(count++)));
	for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
		int row = static_cast<int>(base.fsndbuf(count++));
		x.J.zero_row(row);
		int ncol = static_cast<int>(base.fsndbuf(count++));
		for (int k = 0;k<ncol;++k) {
			int col = static_cast<int>(base.fsndbuf(count++));
			FLT val = base.fsndbuf(count++);
			/* This is a check that we aren't receiving unusued local degrees of freedom */
			if (col < INT_MAX-10 && col > -1) {
				x.J.set_values(row,col,val);
			}
		}
		x.J_mpi.multiply_row(row, 0.0);
	}
	
	for (int m=0;m<base.nmatches();++m) {
		int count = 0;
		int Jstart_mpi = static_cast<int>(base.frcvbuf(m, count++)); // Start of jacobian on matching block
		count++; // int Jstart_mpi_vrtx_unknowns = static_cast<int>(base.frcvbuf(m, count++)); // Start of vertex unknowns on mathcing block (not used typically)
		
		sparse_row_major *pJ_mpi;
		if (base.is_local(m)) {
			pJ_mpi = &x.J;
			Jstart_mpi = 0;
		}
		else {
			pJ_mpi = &x.J_mpi;
		}
		
		vector<int> row_mpi_storage;
		for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
			int row = rowbase + *n;
			int row_mpi = static_cast<int>(base.frcvbuf(m,count++)) +Jstart_mpi;
			row_mpi_storage.push_back(row_mpi);
			int ncol = static_cast<int>(base.frcvbuf(m,count++));
#ifdef MPDEBUG
			*x.gbl->log << "vertex receiving " << ncol << " jacobian entries for vertex " << row/vdofs << " and variable " << *n << std::endl;
#endif
			for (int k = 0;k<ncol;++k) {
				int col = static_cast<int>(base.frcvbuf(m,count++));
				FLT val = base.frcvbuf(m,count++);
				if (col < INT_MAX-10 && col > -1) {
					col += Jstart_mpi;
#ifdef MPDEBUG
					*x.gbl->log  << col << ' ' << val << ' ';
#endif
					(*pJ_mpi).add_values(row,col,val);
				}
			}
#ifdef MPDEBUG
			*x.gbl->log << std::endl;
#endif
		}
		
		/* Shift all diagonal block entries for this vertex */
		for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
			int row = rowbase + *n;
			std::vector<int>::iterator row_mpi=row_mpi_storage.begin();
			for(std::vector<int>::iterator n_mpi=c0_indices_xy.begin();n_mpi != c0_indices_xy.end();++n_mpi) {
#ifdef MPDEBUG
				*x.gbl->log << "vertex swapping " << row << ',' << *row_mpi << " for " << row << ',' << rowbase +*n_mpi << std::endl;
#endif
#ifdef DEBUG_JAC
				if (!x.gbl->jac_debug)
#endif
				{
					FLT dval = (*pJ_mpi)(row,*row_mpi);
					(*pJ_mpi)(row,*row_mpi) = 0.0;
					x.J(row,rowbase+*n_mpi) += dval;
				}
				++row_mpi;
			}
		}
		row_mpi_storage.clear();
	}
	
	for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
		x.J.multiply_row(rowbase+*n,1.0/(1.0+base.nmatches()));
		x.J_mpi.multiply_row(rowbase+*n,1.0/(1.0+base.nmatches()));
	}
	
#else
	This Part not working
#endif
	return(count);
}

void hp_vrtx_bdry::petsc_jacobian_dirichlet() {
	const int nessentials = essential_indices.size();
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	const int v0 = base.pnt;
	const int gind = v0*vdofs;
	Array<int,1> indices(nessentials);
	
	int counter = 0;
	for(std::vector<int>::const_iterator n=essential_indices.begin();n != essential_indices.end();++n) {
		indices(counter++)=gind +*n;
	}
	
#ifdef MY_SPARSE
	x.J.zero_rows(counter,indices);
	x.J_mpi.zero_rows(counter,indices);
	x.J.set_diag(counter,indices,1.0);
#else
	MatZeroRows(x.petsc_J,counter,indices.data(),1.0);
#endif
}
#endif


void multi_physics_pnt::init(input_map& inmap,void* gbl_in) {
	hp_vrtx_bdry::init(inmap,gbl_in);
	match_pairs.resize(base.nmatches());
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	denom.resize(vdofs);
	denom = 1.0;
	
	std::string match_name;
	for(int m=0;m<base.nmatches();++m) {
		base.match_name(m,match_name);
		istringstream nstr;
		std::string vals;
		if (base.idprefix <= match_name) {
			if (!inmap.getline(base.idprefix +'_' +match_name +"_matching",vals)) {
				for(int n=0;n<vdofs;++n) {
					match_pairs[m].push_back(std::pair<int,int>(n,n));
				}
				denom += 1.0;  // All matched
			}
			else {
				nstr.str(vals);
				int idx1, idx2;
				while(nstr >> idx1 >> idx2) {
					match_pairs[m].push_back(std::pair<int,int>(idx1,idx2));
					denom(idx1) += 1.0;
				}
			}
		}
		else {
			if (!inmap.getline(match_name +'_' +base.idprefix +"_matching",vals)) {
				for(int n=0;n<vdofs;++n) {
					match_pairs[m].push_back(std::pair<int,int>(n,n));
				}
				denom += 1.0;  // All matched
			}
			else {
				nstr.str(vals);
				int idx1, idx2;
				while(nstr >> idx2 >> idx1) {
					match_pairs[m].push_back(std::pair<int,int>(idx1,idx2));
					denom(idx1) += 1.0;
				}
			}
		}
	}
}

void multi_physics_pnt::pmatchsolution_rcv(int phase, FLT *pdata, int vrtstride) {
	
	if (!base.in_group(boundary::all_phased)) return;

	int offset = base.pnt*vrtstride*x.NV;
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	Array<FLT,1> counts(vdofs);
	
	counts = 1.;
	for(int m=0;m<base.nmatches();++m) {
		if(base.matchphase(boundary::all_phased,m) != phase) continue;
		
		for(std::vector<std::pair<int,int> >::const_iterator it = match_pairs[m].begin();it != match_pairs[m].end()-(vdofs-x.NV); ++it) {
#ifdef MPDEBUG
			*x.gbl->log << it->first << ' ' << base.fsndbuf(it->first) << ' ' << it->second << ' ' << base.frcvbuf(m,it->second) << std::endl;
#endif
			base.fsndbuf(it->first) += base.frcvbuf(m,it->second);
			counts(it->first) += 1.;
		}
	}
	for(int n=0;n<x.NV;++n) {
		pdata[offset+n] = base.fsndbuf(n)/counts(n);

	}
}

// FIXME: Phasing of vertices with sides is possibly messed up.

#ifdef petsc
int multi_physics_pnt::non_sparse_rcv(int phase, Array<int,1> &nnzero, Array<int,1> &nnzero_mpi) {
	
	if (!base.in_group(boundary::all_phased) || base.matchphase(boundary::all_phased,0) != phase) return(0);
	
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	
	/* Reload to avoid overlap with sides */
	int pind = base.pnt*vdofs;
	int count = 0;
	for(int n=0;n<vdofs;++n) {
		nnzero(pind +n) = base.isndbuf(count++);
		nnzero_mpi(pind +n) = 0;
	}
	
	for(int m=0;m<base.nmatches();++m) {
		Array<int,1> target;
		if (base.is_local(m))  // only 1 matching boundary
			target.reference(nnzero);
		else
			target.reference(nnzero_mpi);
		
		int pind = base.pnt*vdofs;
		for(std::vector<std::pair<int,int> >::const_iterator it = match_pairs[m].begin();it != match_pairs[m].end(); ++it) {
			target(pind +it->first) += base.ircvbuf(m,it->second);
		}
	}
	return(0);
}

int multi_physics_pnt::petsc_matchjacobian_rcv(int phase)	{
	
	if (!base.in_group(boundary::all_phased) || base.matchphase(boundary::all_phased,0) != phase) return(0);
	
	const int vdofs = x.NV +(x.mmovement == tri_hp::coupled_deformable)*x.ND;
	const int rowbase = base.pnt*vdofs;

#ifdef MY_SPARSE
	int count = 0;
	assert(x.jacobian_start == static_cast<int>(base.fsndbuf(count++)));
	assert(jacobian_start == static_cast<int>(base.fsndbuf(count++)));
	for(std::vector<int>::iterator n=c0_indices_xy.begin();n != c0_indices_xy.end();++n) {
		int row = static_cast<int>(base.fsndbuf(count++));
		x.J.zero_row(row);
		int ncol = static_cast<int>(base.fsndbuf(count++));
		for (int k = 0;k<ncol;++k) {
			int col = static_cast<int>(base.fsndbuf(count++));
			FLT val = base.fsndbuf(count++);
			/* This is a check that we aren't receiving unusued local degrees of freedom */
			if (col < INT_MAX-10 && col > -1) {
				x.J.set_values(row,col,val);
			}
		}
		x.J_mpi.multiply_row(row, 0.0);
	}

	for (int m=0;m<base.nmatches();++m) {
		int count = 0;
		int Jstart_mpi = static_cast<int>(base.frcvbuf(m, count++)); // Start of jacobian on matching block
		count++; // int Jstart_mpi_vrtx_unknowns = static_cast<int>(base.frcvbuf(m, count++)); // Start of vertex unknowns on mathcing block (not used typically)
		
		sparse_row_major *pJ_mpi;
		if (base.is_local(m)) {
			pJ_mpi = &x.J;
			Jstart_mpi = 0;
		}
		else {
			pJ_mpi = &x.J_mpi;
		}
		
		vector<int> row_mpi_storage;
		for(std::vector<std::pair<int,int> >::const_iterator it = match_pairs[m].begin();it != match_pairs[m].end(); ++it) {
			int row = rowbase + it->first;
			/* Skip communication that is unrelated */
			count = 2; // Go back to beginning and scan through list
			for(int j=0;j<it->second;++j) {
				count++;  // row_mpi
				int ncol = static_cast<int>(base.frcvbuf(m,count++));
				for (int k = 0;k<ncol;++k) {
					count++; // column
					count++; // value
				}
			}
			
			int row_mpi = static_cast<int>(base.frcvbuf(m,count++)) +Jstart_mpi;
			row_mpi_storage.push_back(row_mpi);
			int ncol = static_cast<int>(base.frcvbuf(m,count++));
#ifdef MPDEBUG
			*x.gbl->log << "vertex receiving " << ncol << " jacobian entries for vertex " << row/vdofs << " and variable " << it->first << " from mpi_row " << row_mpi << std::endl;
#endif
			for (int k = 0;k<ncol;++k) {
				int col = static_cast<int>(base.frcvbuf(m,count++));
				FLT val = base.frcvbuf(m,count++);
				if (col < INT_MAX-10 && col > -1) {
					col += Jstart_mpi;
#ifdef MPDEBUG
					*x.gbl->log  << col << ' ' << val << ' ';
#endif
					(*pJ_mpi).add_values(row,col,val);
				}
			}
#ifdef MPDEBUG
				*x.gbl->log << std::endl;
#endif
		}
		
		/* Shift all diagonal block entries for this vertex */
		for(std::vector<std::pair<int,int> >::const_iterator it = match_pairs[m].begin();it != match_pairs[m].end(); ++it) {
			int row = rowbase + it->first;
			std::vector<int>::iterator row_mpi=row_mpi_storage.begin();
			for(std::vector<std::pair<int,int> >::const_iterator it_mpi = match_pairs[m].begin();it_mpi != match_pairs[m].end(); ++it_mpi) {
#ifdef MPDEBUG
				*x.gbl->log << "vertex swapping " << row << ',' << *row_mpi << " for " << row << ',' << rowbase +it_mpi->first << std::endl;
#endif
#ifdef DEBUG_JAC
				if (!x.gbl->jac_debug)
#endif
				{
					FLT dval = (*pJ_mpi)(row,*row_mpi);
					(*pJ_mpi)(row,*row_mpi) = 0.0;
					x.J(row,rowbase+it_mpi->first) += dval;
				}
				++row_mpi;
			}
		}
		row_mpi_storage.clear();
	}

	for(int n=0;n<vdofs;++n) {
		x.J.multiply_row(rowbase+n,1./denom(n));
		x.J_mpi.multiply_row(rowbase+n,1./denom(n));
	}

#else
	This Part not working
#endif
	return(count);
}

#endif

