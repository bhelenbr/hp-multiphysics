#include "blocks.h"
#include <string>
#include <parseargs.h>
#include "tri_mesh.h"
#include "math.h"
#include <input_map.h>

#ifdef MPISRC
#include <mpi.h>
#endif

using namespace std;

static int informat = 3;
static int outformat = 3;
static GBool Generate = gFalse;
static GBool Shift = gFalse;
static GBool Scale = gFalse;
static GBool Smooth = gFalse;
static GBool Coarsen_hp = gFalse;
static GBool Refineby2 = gFalse;
static GBool Partition = gFalse;
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
	clock_t cpu_time;

#ifdef MPISRC
	int myid;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif
#ifdef PTH
	// For debugging put interrupt here
	// On interrupt type this into gdb console: "handle SIGUSR1 nostop print pass"
	// Then continue
	int rc = pth_init();
	if (!rc) {
		std::cerr << "couldn't start pth environment\n";
	}
#endif

  // parse args
	ok = parseArgs(argDesc, &argc, argv);
	if (!ok || printHelp) {
		fprintf(stderr, "mesh utility ");
		printUsage("mesh", "<inputfile> <outputfile>]", argDesc);
		exit(1);
	}
	tri_mesh::filetype in = static_cast<tri_mesh::filetype>(informat);
	tri_mesh::filetype out = static_cast<tri_mesh::filetype>(outformat);
	std::string bdry_nm(std::string(argv[1]) +"_bdry.inpt");
	ifstream intest;
	input_map bdrymap;
	intest.open(bdry_nm.c_str());
	if (intest.good()) {
		intest.close();
		bdrymap.input(bdry_nm);
		bdrymap.echo = true;
		std::cout << "Using " << bdry_nm << std::endl;
	}

	if (Cut) {
		class tri_mesh zx;
		zx.input(argv[1],in,1.0,bdrymap);
		for(int i=0;i<zx.npnt;++i)
			zx.gbl->fltwk(i) = zx.pnts(i)(0)*zx.pnts(i)(0) +zx.pnts(i)(1)*zx.pnts(i)(1) - 0.25;

		zx.cut();

		return(0);
	}


	/* TO SYMMETRIZE A MESH */
	if (Symmetrize) {
	    class tri_mesh zx;
		zx.input(argv[1],in,8.0,bdrymap);
		zx.symmetrize();
		return 0;
	}

	if (Vlngth) {
	    class tri_mesh zx;

		zx.input(argv[1],in,8.0,bdrymap);
		std::string name;
		name = std::string(argv[1]) +".lngth";
		FILE *fp = fopen(name.c_str(),"w");
		for(int i=0;i<zx.npnt;++i) fprintf(fp,"%e\n",0.3); // 5.*zx.lngth(i));
		fclose(fp);
		return 0;
	}

	if (Smooth) {
	    class tri_mesh zx;

		zx.input(argv[1],in,8.0,bdrymap);
		zx.smooth_cofa(2);
		zx.output(argv[2],out);

		return 0;
	}

	if (Refineby2) {
	    class tri_mesh zx,zy;
		zx.input(argv[1],in,8.0,bdrymap);
		zy.refineby2(zx);
		zy.checkintegrity();
		zy.output(argv[2],out);
		return 0;
	}

	if (Coarsen_hp) {
	    class tri_mesh zx,zy;

		int p;
		zx.input(argv[1],in,1.0,bdrymap);
		printf("input p\n");
		scanf("%d",&p);
		zy.coarsen_substructured(zx,p);
		zy.output(argv[2],out);
		return 0;
	}

	if (Scale) {
	    class tri_mesh zx;
		TinyVector<FLT,2> s;
		printf("Enter x and y scaling\n");
		scanf("%le%le",&s(0),&s(1));
		zx.input(argv[1],in,1.0,bdrymap);
		zx.scale(s);
		zx.output(argv[2],out);
		return 0;
	}

	if (Shift) {
		class tri_mesh zx;

		TinyVector<FLT,2> s;
		printf("Enter x and y shift\n");
		scanf("%le %le",&s(0),&s(1));
		zx.input(argv[1],in,1.0,bdrymap);
		zx.shift(s);
		zx.output(argv[2],out);
		return 0;
	}

	if (Format) {
		class tri_mesh zx;
		zx.input(argv[1],in,1.0,bdrymap);
		zx.output(argv[2],out);
		return(0);
	}

	if (Partition) {
#ifdef METIS
		class tri_mesh zx;
		int p;
		sscanf(argv[2],"%d",&p);
		std::string fname;
		ostringstream nstr;
		zx.input(argv[1],in,1.0,bdrymap);
		zx.setpartition(p);
		Array<tri_mesh,1> zpart(p);
		for(int i=0;i<p;++i) {
			nstr << "b" << i << std::flush;
			fname = "partition_" +nstr.str();
			std::cout << nstr.str() << "_mesh: " << fname << std::endl;
			nstr.str("");
			zpart(i).partition(zx,i);
			zpart(i).checkintegrity();
			zpart(i).output(fname,out);
			//zpart(i).output(fname,tri_mesh::boundary);
		}
#else
		printf("Need metis package to partition\n");
#endif
		return(0);
	}

	if (Coarsen_Marks) {
		class tri_mesh zx;

		zx.input(argv[1],in,1.0,bdrymap);
		FILE *fp = fopen(argv[3],"r");

		for(int i=0;i<zx.npnt;++i) {
			fscanf(fp,"%d\n",&zx.pnt(i).info);
			zx.pnt(i).info = 1-zx.pnt(i).info;
		}
		clock();
		zx.coarsen3();
		cpu_time = clock();
		std::cout << "that took " << cpu_time << " cpu time" << std::endl;

		zx.output(argv[2],out);
		return(0);
	}


	if (argc == 2) {
		/* READ INPUT MAP FROM FILE */
		sim::blks.go(argv[1]);
	}
	else if (argc == 3) {
		/* READ INPUT MAP FROM FILE & OUTPUT TO FILE */
		sim::blks.go(argv[1],argv[2]);
	}

#ifdef PTH
	pth_exit(NULL);
#endif
#ifdef MPISRC
	MPI_Finalize();
#endif

	return(0);
}
