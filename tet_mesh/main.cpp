#include <blocks.h>
#include <string>
#include <parseargs.h>
#include "tet_mesh.h"
#include <input_map.h>

#ifdef MPISRC
#include <mpi.h>
#endif

using namespace std;

static int informat = 4;
static int outformat = 4;
static GBool Generate = gFalse;
static GBool Shift = gFalse;
static GBool Scale = gFalse;
static GBool Smooth = gFalse;
static GBool Coarsen_hp = gFalse;
static GBool Refineby2 = gFalse;
static GBool Partition = gFalse;
static GBool GMSHPartition = gFalse;
static GBool GMSHLabel = gFalse;
static GBool Format = gFalse;
static GBool Coarsen_Marks = gFalse;
static GBool Symmetrize = gFalse;
static GBool Cut = gFalse;
static GBool Vlngth = gFalse;
GBool printHelp = gFalse;


static ArgDesc argDesc[] = {
  {"-g",        argFlag,        &Generate,      0,
		"generate mesh from .d file"},
  {"-m",        argFlag,        &Shift,        0,
		"shift mesh position"},
  {"-s",        argFlag,      &Scale,        0,
		"scale mesh"},
  {"-f",        argFlag,      &Smooth,        0,
		"smooth mesh"},
  {"-h",        argFlag,      &printHelp,      0,
		"print usage information"},
  {"-help",    argFlag,      &printHelp,      0,
		"print usage information"},
  {"-c",        argFlag,      &Coarsen_hp,      0,
		"coarsen substructured mesh"},
  {"-r",        argFlag,      &Refineby2,             0,
		"refine mesh by 2"},
  {"-i",        argInt,      &informat,          0,
		"input mesh format"},
  {"-o", argInt,    &outformat,        0,
		"output format"},
  {"-p"  ,argFlag,     &Partition,            0,
		"partition mesh"},
	{"-gp"  ,argFlag,     &GMSHPartition,            0,
		"partition gmsh physical volumes"},
	{"-gl"  ,argFlag,     &GMSHLabel,            0,
		"create label of gmsh physical volumes"},
  {"-x"  ,argFlag,     &Format,            0,
		"change format"},
  {"-l"  ,argFlag,     &Coarsen_Marks,            0,
		"Coarsen vertices based on list of marks"},
  {"-y"  ,argFlag,     &Symmetrize,            0,
		"Make mesh symmetric about y = 0"},
  {"-z"  ,argFlag,     &Cut,            0,
		"Cut mesh using indicator function"},
  {"-v"  ,argFlag,     &Vlngth,            0,
		"Create a mesh resolution file"},
  {NULL}
};


int main(int argc, char *argv[]) {
	GBool ok;
	
	// tet_mesh zx;
	// zx.test();
	
#ifdef MPISRC
	int myid;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif
	//#ifdef PTH
	//	// For debugging put interrupt here
	//	// On interrupt type this into gdb console: "handle SIGUSR1 nostop print pass"
	//	// Then continue
	//	int rc = pth_init();
	//	if (!rc) {
	//		std::cerr << "couldn't start pth environment\n";
	//	}
	//#endif
	
  // parse args
	ok = parseArgs(argDesc, &argc, argv);
	if (!ok || printHelp) {
		fprintf(stderr, "mesh utility ");
		printUsage("mesh", "<inputfile> <outputfile>]", argDesc);
		exit(1);
	}
	
	tet_mesh::filetype in = static_cast<tet_mesh::filetype>(informat);
	tet_mesh::filetype out = static_cast<tet_mesh::filetype>(outformat);
	std::string bdry_nm(std::string(argv[1]) +"_bdry.inpt");
	ifstream intest;
	input_map inmap;
	intest.open(bdry_nm.c_str());
	if (intest) {
		intest.close();
		inmap.input(bdry_nm);
		inmap.echo = true;
		std::cout << "Using " << bdry_nm << std::endl;
	}
	//
	//    if (Cut) {
	//        class tet_mesh zx;
	//        zx.input(argv[1],in,1.0,inmap);
	//        for(int i=0;i<zx.npnt;++i)
	//            zx.gbl->fltwk(i) = zx.pnts(i)(0)*zx.pnts(i)(0) +zx.pnts(i)(1)*zx.pnts(i)(1) - 0.25;
	//
	//        zx.cut();
	//
	//        return(0);
	//    }
	//
	//
	//    /* TO SYMMETRIZE A MESH */
	//    if (Symmetrize) {
	//        class tet_mesh zx;
	//        zx.input(argv[1],in,8.0,inmap);
	//        zx.symmetrize();
	//        return 0;
	//    }
	//
	
	if (Vlngth) {
		class tet_mesh zx;
		
		zx.input(argv[1],in,8.0,inmap);
		std::string name;
		name = std::string(argv[1]) +".lngth";
		FILE *fp = fopen(name.c_str(),"w");
		for(int i=0;i<zx.npnt;++i) fprintf(fp,"%e\n",0.3); // 5.*zx.lngth(i));
		fclose(fp);
		return 0;
	}
	//
	//    if (Smooth) {
	//        class tet_mesh zx;
	//
	//        zx.input(argv[1],in,8.0,inmap);
	//        zx.smooth_cofa(2);
	//        zx.output(argv[2],out);
	//
	//        return 0;
	//    }
	//
	if (Refineby2) {
		class tet_mesh zx,zy;
		zx.input(argv[1],in,8.0,inmap);
		zy.refineby2(zx);
		zy.checkintegrity();
		zy.output(argv[2],out);
		return 0;
	}
	//
	//    if (Coarsen_hp) {
	//        class tet_mesh zx,zy;
	//
	//        int p;
	//        zx.input(argv[1],in,1.0,inmap);
	//        printf("input p\n");
	//        scanf("%d",&p);
	//        zy.coarsen_substructured(zx,p);
	//        zy.output(argv[2],out);
	//        return 0;
	//    }
	//
	if (Scale) {
		class tet_mesh zx;
		TinyVector<FLT,tet_mesh::ND> s;
		printf("Enter x y and z scaling\n");
		scanf("%le%le%le",&s(0),&s(1),&s(2));
		zx.input(argv[1],in,1.0,inmap);
		zx.scale(s);
		zx.output(argv[2],out);
		return 0;
	}
	
	if (Shift) {
		class tet_mesh zx;
		
		TinyVector<FLT,tet_mesh::ND> s;
		printf("Enter x y and z shift\n");
		scanf("%le%le%le",&s(0),&s(1),&s(2));
		zx.input(argv[1],in,1.0,inmap);
		zx.shift(s);
		zx.output(argv[2],out);
		return 0;
	}
	
	if (Format) {
		class tet_mesh zx;
		zx.input(argv[1],in,1.0,inmap);
		zx.output(argv[2],out);
		return(0);
	}

	
	if (GMSHPartition) {
		/* This separates different volumes in a GMSH mesh */
		/* The volumes must be numbered sequentially starting at 1 */
		class tet_mesh zx;
		int p;
		sscanf(argv[2],"%d",&p);
		std::string fname;
		ostringstream nstr;
		
		zx.input(argv[1],in,1.0,inmap);
		
		/* input calls setinfo but to make sure call it again because partition needs it to work */
		zx.tet_mesh::setinfo();
		
		for(int i=0;i<zx.ntet;++i)
			zx.tet(i).info = zx.gbl->fltwk(i)-1;
				
		tet_mesh *zpart;
		
		Array<int,2> blist;
		Array<int,1> bnum;
		
		zx.setup_partition(p,blist,bnum);
		
		for(int i=0;i<p;++i) {
			nstr << "b" << i << std::flush;
			fname = "partition_" +nstr.str();
			std::cout << nstr.str() << "_mesh: " << fname << std::endl;
			nstr.str("");
			zpart = new tet_mesh;
			zpart->partition2(zx,i,p,blist,bnum);
			zpart->output(fname,tet_mesh::gmsh);
			zpart->output(fname);
			delete zpart;
		}
		return(0);
	}
	
	if (GMSHLabel) {
		/* This creates a file containing labels of the elements in different volumes in a GMSH mesh */
		class tet_mesh zx;
		zx.input(argv[1],in,1.0,inmap);

		string gridname(argv[1]), filename;
		size_t dotloc;
		dotloc = gridname.find_last_of('.');
		
		if (dotloc != string::npos) {
			/* Found and ending */
			filename = gridname.substr(0,dotloc) +".marks";
		}
		else {
			filename = gridname +".marks";
		}
		
		
		ofstream fout;
		fout.open(filename.c_str());
		
		for(int i=0;i<zx.ntet;++i)
			fout << i << ": " << static_cast<int>(zx.gbl->fltwk(i)-1) << '\n';
		
		fout.close();
		
		return(0);
	}
	
	if (Partition) {
#ifdef METIS
		class tet_mesh zx;
		int p;
		sscanf(argv[2],"%d",&p);
		std::string fname,mname;
		ostringstream nstr;
		zx.input(argv[1],in,1.0,inmap);
		
		/* input calls setinfo but to make sure call it again because partition needs it to work */
		zx.tet_mesh::setinfo();
		
		/* Uses metis to set tet(i).info to a partition number */
		zx.setpartition(p);
		
		//Array<int,2> blist;
		//Array<int,1> bnum;
		// zx.setup_partition(p,blist,bnum);
		zx.setup_partition2(p);
		
		/* The following is to load a marker file if it exists */
		ifstream fin;
		bool marks_flag = false;
		Array<int,1> marks;
		string gridname(argv[1]), filename;
		size_t dotloc;
		
		dotloc = gridname.find_last_of('.');
		if (dotloc != string::npos) {
			/* Found and ending */
			filename = gridname.substr(0,dotloc) +".marks";
		}
		else {
			filename = gridname +".marks";
		}
		fin.open(filename.c_str());
		if (fin) {
			marks.resize(zx.ntet);
			for (int i=0;i<zx.ntet;++i) {
				fin.ignore(80,':');
				fin >> marks(i);
				marks_flag = true;
			}
		}
		
			
		/* now partition mesh and marks if found */
		for(int i=0;i<p;++i) {
			tet_mesh zpart;
			nstr << "b" << i << std::flush;
			fname = "partition_" +nstr.str();
			mname = fname +".marks";
			// std::cout << nstr.str() << "_mesh: " << fname << std::endl;
			nstr.str("");
			// zpart(i).partition2(zx,i,p,blist,bnum);
			zpart.partition3(zx,i);

			if (marks_flag) {
				ofstream fout;
				fout.open(mname.c_str());
				
				/* tet.info gets wiped out so have to redo it */
				int ntet = 0;
				for(int tind=0;tind<zx.ntet;++tind) {
					if (zx.tet(tind).info == i) {
						fout << ntet << ": " << marks(tind) << '\n';
						++ntet;
					}
				}
				fout.close();
			}
			
			// zpart(i).output(fname,tet_mesh::gmsh);
			zpart.output(fname);
		}
		
		//        for(int i=0;i<p;++i) {
		//			nstr << "b" << i << std::flush;
		//			fname = "partition_" +nstr.str();
		//			std::cout << nstr.str() << "_mesh: " << fname << std::endl;
		//			nstr.str("");
		//			zpart(i).partition(zx,i,p);
		//			zpart(i).checkintegrity();
		//			zpart(i).output(fname,out);
		//			zpart(i).output(fname,tet_mesh::gmsh);
		//
		//            //zpart(i).output(fname,tet_mesh::boundary);//temp fixme
		//        }
#else
		printf("Need metis package to partition\n");
#endif
		return(0);
	}
	
	//    if (Coarsen_Marks) {
	//        class tet_mesh zx;
	//
	//        zx.input(argv[1],in,1.0,inmap);
	//        FILE *fp = fopen(argv[3],"r");
	//
	//        for(int i=0;i<zx.npnt;++i) {
	//            fscanf(fp,"%d\n",&zx.pnt(i).info);
	//            zx.pnt(i).info = 1-zx.pnt(i).info;
	//        }
	//        zx.coarsen3();
	//        zx.output(argv[2],out);
	//        return(0);
	//    }
	
	
	if (argc == 2) {
		/* READ INPUT MAP FROM FILE */
		sim::blks.go(argv[1]);
	}
	else if (argc == 3) {
		/* READ INPUT MAP FROM FILE & OUTPUT TO FILE */
		sim::blks.go(argv[1],argv[2]);
	}
	
#ifdef PTH
	pth_kill();
#endif
#ifdef MPISRC
	MPI_Finalize();
#endif
	
	return(0);
}


multigrid_interface* block::getnewlevel(input_map& input) {
	int type;
	multigrid_interface *temp;
	
	if (!input.get(idprefix+"_type",type)) input.getwdefault("type",type,1);
	
	switch(type) {
		default: {
			temp = new tet_mesh();
			break;
		}
	} 
	
	return(temp);
}


