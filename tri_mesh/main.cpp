#include "blocks.h"
#include <string>
#include <parseargs.h>
#include "tri_mesh.h"
#include "math.h"
#include <input_map.h>
#include <symbolic_function.h>

#ifdef MPISRC
#include <mpi.h>
#endif

using namespace std;

static int informat = 3;
static int outformat = 3;
static GBool Echo = gFalse;
static GBool Generate = gFalse;
static GBool Shift = gFalse;
static GBool Scale = gFalse;
static GBool Smooth = gFalse;
static GBool Coarsen_hp = gFalse;
static GBool Refineby2 = gFalse;
static GBool Partition = gFalse;
static GBool Format = gFalse;
static GBool Coarsenby2 = gFalse;
static GBool Symmetrize = gFalse;
static GBool Cut = gFalse;
static GBool Vlngth = gFalse;
static GBool Append = gFalse;
GBool printHelp = gFalse;

static ArgDesc argDesc[] = {
  {"-A",        argFlag,        &Append,      0,
	"append two meshes"},
	{"-e",        argFlag,        &Echo,      0,
		"display basic mesh information"},
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
	{"-c"  ,argFlag,     &Coarsenby2,            0,
	"Coarsen mesh by factor of 2"},
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
  {"-l",        argFlag,      &Coarsen_hp,      0,
	"coarsen substructured mesh"},
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
	// On interrupt type this into gdb console:
	// for gdb: handle SIGUSR1 nostop print pass
	// for lldb: pro hand -p true -s false SIGUSR1
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
		sim::abort(__LINE__,__FILE__,&std::cerr);
	}
	tri_mesh::filetype in = static_cast<tri_mesh::filetype>(informat);
	tri_mesh::filetype out = static_cast<tri_mesh::filetype>(outformat);
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

	if (Echo) {
		class tri_mesh zx;
		zx.input(argv[1],in,1.0,inmap);
		std::cout << "npnt: " << zx.npnt << " nseg: " << zx.nseg << " ntri: " << zx.ntri << std::endl;
		return(0);
	}
	

	if (Cut) {
		class tri_mesh zx;
		zx.input(argv[1],in,1.0,inmap);
		for(int i=0;i<zx.npnt;++i)
			zx.gbl->fltwk(i) = zx.pnts(i)(0)*zx.pnts(i)(0) +zx.pnts(i)(1)*zx.pnts(i)(1) - 0.25;

		zx.cut();

		return(0);
	}
	
	if (Append) {
		class tri_mesh zx,zy;
		zx.input(argv[1],in,100.0,inmap);
		zy.input(argv[2],in,1.0,inmap);
		zx.append(zy);
		zx.cleanup_after_adapt();
		zx.output(argv[3],out);
		/* Output a marks file so you know which element came from where */
		ofstream file;
		std::string filename;
		filename = std::string(argv[3]) +".marks";
		file.open(filename.c_str());
		for(int i=0;i<zx.ntri-zy.ntri;++i) {
			file << 0 << std::endl;
		}
		for (int i=0;i<zy.ntri;++i) {
			file << 1 << std::endl;
		}
		file.close();
		
		return(0);
	}

	/* TO SYMMETRIZE A MESH */
	if (Symmetrize) {
		class tri_mesh zx;
		zx.input(argv[1],in,8.0,inmap);
		zx.symmetrize();
		return 0;
	}

	if (Vlngth) {
		class tri_mesh zx;
		zx.input(argv[1],in,8.0,inmap);
		
		symbolic_function<2> length_modifier_function;
		input_map length_modifier_input;
		std::string length_modifier_string;
		std::cout << "Input function of (x0,x1) and t where t is the current mesh length" << std::endl;
		std::cin >> length_modifier_string;
		length_modifier_input["length_modifier"] = length_modifier_string;
		length_modifier_function.init(length_modifier_input,"length_modifier");
		for(int i=0;i<zx.npnt;++i) 
			zx.lngth(i) = length_modifier_function.Eval(zx.pnts(i),zx.lngth(i));
		zx.output(argv[2],out);
		return 0; 
	}

	if (Smooth) {
		class tri_mesh zx;
		zx.input(argv[1],in,8.0,inmap);
		zx.smooth_cofa(2);
		zx.output(argv[2],out);

		return 0;
	}

	if (Refineby2) {
		class tri_mesh zx,zy;
		zx.input(argv[1],in,8.0,inmap);
		zy.refineby2(zx);
		zy.checkintegrity();
		zy.output(argv[2],out);
		return 0;
	}

	if (Coarsen_hp) {
		class tri_mesh zx,zy;

		int p;
		zx.input(argv[1],in,1.0,inmap);
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
		zx.input(argv[1],in,1.0,inmap);
		zx.scale(s);
		zx.output(argv[2],out);
		return 0;
	}

	if (Shift) {
		class tri_mesh zx;

		TinyVector<FLT,2> s;
		printf("Enter x and y shift\n");
		scanf("%le %le",&s(0),&s(1));
		zx.input(argv[1],in,1.0,inmap);
		zx.shift(s);
		zx.output(argv[2],out);
		return 0;
	}

	if (Format) {
		class tri_mesh zx;
		zx.input(argv[1],in,1.0,inmap);
		zx.output(argv[2],out);
		return(0);
	}

	if (Partition) {
#ifdef METIS
		class tri_mesh zx;
		int p;
		sscanf(argv[2],"%d",&p);
		if (p > 1) {
			
			/* Load mesh */
			std::string fname;
			ostringstream nstr;
			zx.input(argv[1],in,1.0,inmap);
			
			/* If there is a marks file then try to subdivide appended multiphysics mesh */
			ifstream marks_file;
			std::string marksfilename = std::string(argv[1]) +".marks";
			marks_file.open(marksfilename.c_str());
			if (marks_file) {
				for(int i=0;i<zx.ntri;++i)
					marks_file >> zx.tri(i).info;
				zx.subpartition(p);
			}
			else {
				zx.setpartition(p);
			}

			/* Extract partitions */
			for(int i=0;i<p;++i) {
				tri_mesh zpart;
				nstr << "b" << i << std::flush;
				fname = "partition_" +nstr.str();
				std::cout << nstr.str() << "_mesh: " << fname << std::endl;
				nstr.str("");
				zpart.partition(zx,i);
				zpart.checkintegrity();
				zpart.output("partition",out);
				//zpart(i).output(fname,tri_mesh::boundary);
			}
		}
#else
		printf("Need metis package to partition\n");
#endif
		return(0);
	}
	
	if (Coarsenby2) {
		class tri_mesh zx,zy;

		zx.input(argv[1],in,1.0,inmap);
		
		for(int i=0;i<zx.npnt;++i) {
			zx.lngth(i) *= 2.0;
		}
		clock();
		zy.coarsen(1.6,zx);
		cpu_time = clock();
		std::cout << "that took " << cpu_time << " cpu time" << std::endl;

		zy.output(argv[2],out);
		return(0);
	}
	
//	if (Coarsen_Marks) {
//		class tri_mesh zx;
//
//		zx.input(argv[1],in,1.0,inmap);
//		FILE *fp = fopen(argv[3],"r");
//
//		for(int i=0;i<zx.npnt;++i) {
//			fscanf(fp,"%d\n",&zx.pnt(i).info);
//			zx.pnt(i).info = 1-zx.pnt(i).info;
//		}
//		clock();
//		zx.coarsen3();
//		cpu_time = clock();
//		std::cout << "that took " << cpu_time << " cpu time" << std::endl;
//
//		zx.output(argv[2],out);
//		return(0);
//	}


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
