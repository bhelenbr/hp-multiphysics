#include "blocks.h"
#include <map>
#include <string>
#include <parseargs.h>
#include "mesh.h"
#include "math.h"

#ifdef MPISRC
#include <mpi.h>
#endif

using namespace std;

/* GENERATE CONVERT MAXMIN HEIGHT SHIFT COARSEN_TECPLOT SCALE REFINEBY2 */
#define RUN

static int informat = 3;
static int outformat = 3;
static GBool Generate = gFalse;
static GBool Shift = gFalse;
static GBool Scale = gFalse;
static GBool Coarsen_hp = gFalse;
static GBool Refineby2 = gFalse;
static GBool Partition = gFalse;
static GBool Format = gFalse;
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
  {NULL}
};


int main(int argc, char *argv[]) {
   GBool ok;
   char *outfile;
   
  // parse args
   ok = parseArgs(argDesc, &argc, argv);
   if (argc == 2) {
      outfile = argv[1];
   }
   else if (argc == 3) {
      outfile = argv[2];
   }
   else if (!ok || argc < 2 || argc > 3 || printHelp) {
      fprintf(stderr, "mesh utility ");
      printUsage("mesh", "<inputfile> <outputfile>]", argDesc);
      exit(1);
   }

   class mesh zx, zy;
   ftype::name in = static_cast<ftype::name>(informat);
   ftype::name out = static_cast<ftype::name>(outformat);
   
   if (Refineby2) {
      zx.input(argv[1],in,8.0);
      zy.refineby2(zx);
      zy.output(outfile,out);
      return 0;
   }

   if (Coarsen_hp) {
      int p;
      zx.input(argv[1],in);
      printf("input # of times to unrefineby2\n");
      scanf("%d",&p);
      zy.coarsen_substructured(zx,p);
      zy.output(outfile,out);
      return 0;
   }

   if (Scale) {    
      FLT s[2];
      printf("Enter x and y scaling\n");
      scanf("%le%le",s,s+1);
      zx.input(argv[1],in);
      zx.scale(s);
      zx.output(outfile,out);
      return 0;
   }
   
   if (Shift) {    
      FLT s[2];
      printf("Enter x and y shift\n");
      scanf("%le%le",s,s+1);
      zx.input(argv[1],in);
      zx.shift(s);
      zx.output(outfile,out);
      return 0;
   }
   
   if (Format) {
      class mesh zx;
      zx.input(argv[1],in);
      zx.output(outfile,out);
      return(0);
   }
   
   if (Partition) {
      int p;
      char outappend[100];
      printf("input # of partitions\n");
      // scanf("%d",&p);
      p = 4;
      zx.input(argv[1],in);
      zx.setpartition(p);
      mesh zpart[p];
      
      for(int i=0;i<p;++i) {
         number_str(outappend,outfile,i,1);
         zpart[i].partition(zx,i);
         zpart[i].output(outappend,out);
      }
      return(0);
   }

   class blocks z;
   if (Generate) {
      z.init(argv[1]);
      for (int i=0;i<2;++i)
         z.restructure();
      z.output(outfile,out);
      return 0;
   }

#ifdef MPISRC
   int myid;
   MPI_Init(&argc,&argv);
   MPI_Comm_rank(MPI_COMM_WORLD,&myid);
#endif
   
   if (argc == 2) {
      /* READ INPUT MAP FROM FILE */
      z.init(argv[1]);
   }
   else if (argc == 3) {
      /* READ INPUT MAP FROM FILE & OUTPUT TO FILE */
      z.init(argv[1],outfile);
   }
   else {
      /* CREATE INPUT MAPS HERE*/
      map<string,string> input;
      
      input["mglvls"] = "2";
      input["ngrid"] = "2";
      input["ntstep"] = "1";
      input["fadd"] = "0.0";
      input["vnn"] = "0.5";
      input["itercrsn"] = "1";
      input["iterrfne"] = "0";
      input["njacobi"] = "1";
      input["ncycle"] = "10";
      input["vwcycle"] = "1";
      input["tolerance"] = "0.66";
      input["nologfile"] = "duh";
      
#ifdef MPISRC
      /* LIST OF NBLOCKS FOR EACH PROCESSOR */
      input["nblock"] = "1 1";
#else
      input["nblock"] = "2";
#endif
      input["b0.type"] = "1";
      input["b0.mesh"] = "/Users/helenbrk/Codes/grids/WIND/PRDC/top8";
      input["b0.filetype"] = "0";
      input["b0.growth factor"] = "4.0";
      input["b0.bdryfile"] = "/Users/helenbrk/Codes/grids/WIND/PRDC/top.bdry.inpt";
      input["b1.type"] = "1";
      input["b1.mesh"] = "/Users/helenbrk/Codes/grids/WIND/PRDC/bot8";
      input["b1.filetype"] = "0";
      input["b1.growth factor"] = "4.0";
      input["b1.bdryfile"] = "/Users/helenbrk/Codes/grids/WIND/PRDC/bot.bdry.inpt";
      z.init(input);
   }

   z.go();
   
#ifdef MPISRC
   MPI_Finalize();
#endif

   return(0);
}
