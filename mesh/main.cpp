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
  {"-g",      argFlag,      &Generate,     0,
   "generate mesh from .d file"},
  {"-m",      argFlag,      &Shift,      0,
   "shift mesh position"},
  {"-s",      argFlag,     &Scale,      0,
   "scale mesh"},
  {"-h",      argFlag,     &printHelp,     0,
   "print usage information"},
  {"-help",   argFlag,     &printHelp,     0,
   "print usage information"},
  {"-c",      argFlag,     &Coarsen_hp,     0,
   "coarsen substructured mesh"}, 
  {"-r",      argFlag,     &Refineby2,          0,
   "refine mesh by 2"},
  {"-i",      argInt,     &informat,        0,
   "input mesh format"},
  {"-o", argInt,   &outformat,      0,
   "output format"},
  {"-p"  ,argFlag,    &Partition,         0,
   "partition mesh"},
  {"-x"  ,argFlag,    &Format,         0,
   "change format"},
  {"-l"  ,argFlag,    &Coarsen_Marks,         0,
   "Coarsen vertices based on list of marks"},
  {"-y"  ,argFlag,    &Symmetrize,         0,
   "Make mesh symmetric about y = 0"},
  {"-v"  ,argFlag,    &Vlngth,         0,
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
   
   /* TO SYMMETRIZE A MESH */
   if (Symmetrize) {
      zx.input(argv[1],in,8.0);
      zx.symmetrize();
      return 0;
   }
   
   if (Vlngth) {
      zx.input(argv[1],in,8.0);
      std::string name;
      name = std::string(argv[1]) +".vlngth";
      FILE *fp = fopen(name.c_str(),"w");
      for(int i=0;i<zx.nvrtx;++i) fprintf(fp,"%e\n",0.3); // 5.*zx.vlngth(i));
      fclose(fp);
      return 0;
   }
   
   if (Refineby2) {
      zx.input(argv[1],in,8.0);
      zy.refineby2(zx);
      zy.output(argv[2],out);
      return 0;
   }

   if (Coarsen_hp) {
      int p;
      zx.input(argv[1],in);
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
      zx.input(argv[1],in);
      zx.scale(s);
      zx.output(argv[2],out);
      return 0;
   }
   
   if (Shift) {    
      TinyVector<FLT,2> s;
      printf("Enter x and y shift\n");
      scanf("%le%le",&s(0),&s(1));
      zx.input(argv[1],in);
      zx.shift(s);
      zx.output(argv[2],out);
      return 0;
   }
   
   if (Format) {
      class mesh zx;
      zx.input(argv[1],in);
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
      zx.input(argv[1],in);
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
      zx.input(argv[1],in);
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
   else {
      /* CREATE INPUT MAPS HERE*/
      input_map input;
      
      input["mglvls"] = "3";
      input["ngrid"] = "3";
      input["ntstep"] = "2";
      input["fadd"] = "1.0";
      input["vnn"] = "0.9";
      input["itercrsn"] = "1";
      input["iterrfne"] = "0";
      input["njacobi"] = "1";
      input["ncycle"] = "100";
      input["vwcycle"] = "1";
      input["tolerance"] = "2.2";
      input["nologfile"] = "duh";
      input["adapt"] = "1";
      
#ifdef MPISRC
      /* LIST OF NBLOCKS FOR EACH PROCESSOR */
      input["nblock"] = "1 1";
#else
      input["nblock"] = "2";
#endif
      input["b0.type"] = "1";
      input["b0.mesh"] = "${HOME}/Codes/grids/WIND/PRDC/top8";
      input["b0.filetype"] = "3";
      input["b0.growth factor"] = "4.0";
      input["b0.bdryfile"] = "${HOME}/Codes/grids/WIND/PRDC/top_bdry.inpt";
      input["b0.bdryfile"] = "${HOME}/Codes/grids/WIND/PRDC/top_bdry_phased.inpt";

      input["b1.type"] = "1";
      input["b1.mesh"] = "${HOME}/Codes/grids/WIND/PRDC/bot8";
      input["b1.filetype"] = "3";
      input["b1.growth factor"] = "4.0";
      input["b1.bdryfile"] = "${HOME}/Codes/grids/WIND/PRDC/bot_bdry.inpt";
      input["b1.bdryfile"] = "${HOME}/Codes/grids/WIND/PRDC/bot_bdry_phased.inpt";

      std::cout << input;
      sim::blks.init(input);
   }

   sim::blks.go();
   
#ifdef MPISRC
   MPI_Finalize();
#endif

   return(0);
}
