#include "blocks.h"
#include <string>
#include <parseargs.h>
#include "mesh.h"
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
static GBool Coarsen_hp = gFalse;
static GBool Refineby2 = gFalse;
static GBool Partition = gFalse;
static GBool Format = gFalse;
static GBool Coarsen_Marks = gFalse;
static GBool Symmetrize = gFalse;
static GBool Vlngth = gFalse;
GBool printHelp = gFalse;


static ArgDesc argDesc[] = {
  {"-g",        argFlag,        &Generate,      0,
    "generate mesh from .d file"},
  {"-m",        argFlag,        &Shift,        0,
    "shift mesh position"},
  {"-s",        argFlag,      &Scale,        0,
    "scale mesh"},
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
  {"-v"  ,argFlag,     &Vlngth,            0,
    "Create a mesh resolution file"},
  {NULL}
};


int main(int argc, char *argv[]) {
    GBool ok;

  // parse args
    ok = parseArgs(argDesc, &argc, argv);
    if (!ok || printHelp) {
        fprintf(stderr, "mesh utility ");
        printUsage("mesh", "<inputfile> <outputfile>]", argDesc);
        exit(1);
    }
    
    class mesh zx, zy;
    mesh::filetype in = static_cast<mesh::filetype>(informat);
    mesh::filetype out = static_cast<mesh::filetype>(outformat);
    
    std::string bdry_nm(std::string(argv[1]) +"_bdry.inpt");
    ifstream intest;
    input_map bdrymap;
    intest.open(bdry_nm.c_str());
    if (intest.good()) {
        intest.close();
        bdrymap.input(bdry_nm);
    }

    /* TO SYMMETRIZE A MESH */
    if (Symmetrize) {
        zx.input(argv[1],in,8.0,bdrymap);
        zx.symmetrize();
        return 0;
    }
    
    if (Vlngth) {
        zx.input(argv[1],in,8.0,bdrymap);
        std::string name;
        name = std::string(argv[1]) +".vlngth";
        FILE *fp = fopen(name.c_str(),"w");
        for(int i=0;i<zx.nvrtx;++i) fprintf(fp,"%e\n",0.3); // 5.*zx.vlngth(i));
        fclose(fp);
        return 0;
    }
    
    if (Refineby2) {
        zx.input(argv[1],in,8.0,bdrymap);
        zy.refineby2(zx);
        zy.checkintegrity();
        zy.output(argv[2],out);
        return 0;
    }

    if (Coarsen_hp) {
        int p;
        zx.input(argv[1],in,1.0,bdrymap);
        printf("input p\n");
        scanf("%d",&p);
        zy.coarsen_substructured(zx,p);
        zy.output(argv[2],out);
        return 0;
    }

    if (Scale) {     
        TinyVector<FLT,2> s;
        printf("Enter x and y scaling\n");
        scanf("%le%le",&s(0),&s(1));
        zx.input(argv[1],in,1.0,bdrymap);
        zx.scale(s);
        zx.output(argv[2],out);
        return 0;
    }
    
    if (Shift) {     
        TinyVector<FLT,2> s;
        printf("Enter x and y shift\n");
        scanf("%le%le",&s(0),&s(1));
        zx.input(argv[1],in,1.0,bdrymap);
        zx.shift(s);
        zx.output(argv[2],out);
        return 0;
    }
    
    if (Format) {
        class mesh zx;
        zx.input(argv[1],in,1.0,bdrymap);
        zx.output(argv[2],out);
        return(0);
    }
    
    if (Partition) {
#ifdef METIS
        int p;
        std::string fname;
        ostringstream nstr;
        std::cout << "input # of partitions" << std::endl;
        std::cin >> p;
        zx.input(argv[1],in,1.0,bdrymap);
        zx.setpartition(p);
        Array<mesh,1> zpart(p);
        
        for(int i=0;i<p;++i) {
            nstr << i << std::flush;
            fname = argv[1] +nstr.str();
            nstr.str("");
            zpart(i).partition(zx,i);
            zpart(i).output(fname,out);
            zpart(i).bdry_output(fname);
        }
#else
        printf("Need metis package to partition\n");
#endif
        return(0);
    }
    
    if (Coarsen_Marks) {
        zx.input(argv[1],in,1.0,bdrymap);
        FILE *fp = fopen(argv[3],"r");
        for(int i=0;i<zx.nvrtx;++i) {
            fscanf(fp,"%d\n",&zx.vd(i).info);
            zx.vd(i).info = 1-zx.vd(i).info;
        }
        zx.coarsen3();
        zx.output(argv[2],out); 
        return(0);      
    }
        

    if (Generate) {
        sim::blks.init(argv[1]);
        for (int i=0;i<1;++i)
            sim::blks.restructure();
        sim::blks.output(argv[1]);
        return(0);
    }

#ifdef MPISRC
    int myid;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif
    
    if (argc == 2) {
        /* READ INPUT MAP FROM FILE */
        sim::blks.init(argv[1]);
    }
    else if (argc == 3) {
        /* READ INPUT MAP FROM FILE & OUTPUT TO FILE */
        sim::blks.init(argv[1],argv[2]);
    }
    
    sim::blks.go();
    
#ifdef MPISRC
    MPI_Finalize();
#endif

    return(0);
}
