#include "blocks.h"
#include <map>
#include <string>
#include "mesh.h"
#include "math.h"

#ifdef MPISRC
#include <mpi.h>
#endif

#ifdef CAPRI
#include <capri.h>
#endif

using namespace std;

#ifdef CAPRI
CAPRI_MAIN(int argc, char *argv[]) {
#else
int main(int argc, char *argv[]) {
#endif

   class blocks z;
   
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
      z.init(argv[1],argv[2]);
   }
   else {
      /* CREATE INPUT MAPS HERE*/
      map<string,string> input[3];
      
      input[0]["mglvls"] = "1";
      input[0]["ngrid"] = "1";
      input[0]["ntstep"] = "1";
      input[0]["fadd"] = "1.0";
      input[0]["vnn"] = "0.5";
      input[0]["itercrsn"] = "1";
      input[0]["iterrfne"] = "0";
      input[0]["njacobi"] = "1";
      input[0]["ncycle"] = "1000";
      input[0]["vwcycle"] = "2";
      input[0]["tolerance"] = "0.66";
      input[0]["logfile"] = "duh";
      
#ifdef MPISRC
      /* LIST OF NBLOCKS FOR EACH PROCESSOR */
      input[0]["nblock"] = "1 1";
      if (myid == 0) {
         input[1]["blktype"] = "0";
         input[1]["mesh"] = "/Users/helenbrk/Codes/grids/WIND/PRDC/top8";
         input[1]["growth factor"] = "4.0";
         input[1]["filetype"] = "0";
      }
      else {
         input[1]["blktype"] = "0";
         input[1]["mesh"] = "/Users/helenbrk/Codes/grids/WIND/PRDC/bot8";
         input[1]["growth factor"] = "4.0";
         input[1]["filetype"] = "0";
      }
#else
      input[0]["nblock"] = "2";
      input[1]["blktype"] = "0";
      input[1]["mesh"] = "/Users/helenbrk/Codes/grids/WIND/PRDC/top8";
      input[1]["growth factor"] = "4.0";
      input[1]["filetype"] = "0";
      input[2]["blktype"] = "0";
      input[2]["mesh"] = "/Users/helenbrk/Codes/grids/WIND/PRDC/bot8";
      input[2]["growth factor"] = "4.0";
      input[2]["filetype"] = "0";
#endif
      z.init(input);
   }

   z.go();
   
#ifdef MPISRC
   MPI_Finalize();
#endif

   return(0);
}
