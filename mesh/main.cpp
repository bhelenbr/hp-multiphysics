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

/* GENERATE CONVERT MAXMIN HEIGHT SHIFT COARSEN_TECPLOT SCALE */
#define NOSCALE

#ifdef CAPRI
CAPRI_MAIN(int argc, char *argv[]) {
#else
int main(int argc, char *argv[]) {
#endif

#ifdef COARSEN_TECPLOT
   class mesh<2> zx, zy;
   int i,sind,sign,p,p2;
   
   if (argc != 3) {
      printf("usage: tectohp filename p");
      exit(1);
   }
   zx.in_mesh(argv[1],ftype::tecplot);
   sscanf(argv[2],"%d",&p);

   zy.copy(zx);
   
   p2 = p*p;
   zy.ntri = zx.ntri/p2;
   zy.nside = (zy.ntri*3 +zx.sbdry[0]->nsd()/p)/2;
   zy.nvrtx = zx.nvrtx -zy.nside*(p-1) -zy.ntri*((p-1)*(p-2))/2;
   
   for(i=0;i<zx.ntri;i+=p2) {
      zy.tvrtx[i/p2][0] = zx.tvrtx[i+2*p-2][2];
      zy.tvrtx[i/p2][1] = zx.tvrtx[i][0];
      zy.tvrtx[i/p2][2] = zx.tvrtx[i+p2-1][1];
      sind = (zx.tvrtx[i][1] -zy.nvrtx)/(p-1);
      sign = ((zx.tvrtx[i][1] -zy.nvrtx)%(p-1) == 0 ? 1 : -1);
      zy.tside[i/p2].side[0] = sind; 
      zy.tside[i/p2].sign[0] = sign;
      zy.svrtx[sind][(1-sign)/2] = zy.tvrtx[i/p2][1];
      zy.svrtx[sind][(1+sign)/2] = zy.tvrtx[i/p2][2];
      
      sind = (zx.tvrtx[i+p2-1][2] -zy.nvrtx)/(p-1);
      sign = ((zx.tvrtx[i+p2-1][2] -zy.nvrtx)%(p-1) == 0 ? 1 : -1);
      zy.tside[i/p2].side[1] = sind;
      zy.tside[i/p2].sign[1] = sign;
      zy.svrtx[sind][(1-sign)/2] = zy.tvrtx[i/p2][2];
      zy.svrtx[sind][(1+sign)/2] = zy.tvrtx[i/p2][0];
      
      sind = (zx.tvrtx[i+2*p-2][0] -zy.nvrtx)/(p-1);
      sign = ((zx.tvrtx[i+2*p-2][0] -zy.nvrtx)%(p-1) == 0 ? 1 : -1);
      zy.tside[i/p2].side[2] = sind;
      zy.tside[i/p2].sign[2] = sign;
      zy.svrtx[sind][(1-sign)/2] = zy.tvrtx[i/p2][0];
      zy.svrtx[sind][(1+sign)/2] = zy.tvrtx[i/p2][1];
   }
   zy.sbdry[0]->nsd() = zx.sbdry[0]->nsd()/p;
   i = 0;
   if (zx.svrtx[zx.sbdry[0]->sd(0)][0] < zy.nvrtx) ++i;
   for(;i<zx.sbdry[0]->nsd();i+=p) {
      sind = (zx.svrtx[zx.sbdry[0]->sd(i)][0] -zy.nvrtx)/(p-1);
      zy.sbdry[0]->sd(i/p) = sind;
   }
   zy.out_mesh(argv[1],ftype::grid);
   
   return 0;
#endif
      
   

#ifdef HEIGHT
    int count,sind;
    class mesh<2> zx;
    
    zx.in_mesh(argv[1],ftype::tecplot);
    sscanf(argv[1],"data%d",&count);
    
    double xmax = -1.0e16;
    double xmin = 1.0e16;
    
    for(int i=0;i<zx.sbdry[0]->nsd();++i) {
        sind = zx.sbdry[0]->sd(i);
        if (zx.vrtx[zx.svrtx[sind][1]][0] == 0.0) {
            xmax = max(zx.vrtx[zx.svrtx[sind][1]][1],xmax);
            xmin = min(zx.vrtx[zx.svrtx[sind][1]][1],xmin);
        }
    }
    printf("%d %e\n",count,xmax-xmin);
    return 0;
#endif

#ifdef SCALE
    class mesh<2> zx;
    double scale;
    
    sscanf(argv[1],"%le",&scale);
    zx.in_mesh(argv[2],ftype::easymesh);
    for(int i=0;i<zx.nvrtx;++i)
       zx.vrtx[i][1] *= scale;
    zx.out_mesh(argv[2],ftype::grid);
    return 0;
#endif

#ifdef SHIFT
    class mesh<2> zx;
    double offset;
    
    sscanf(argv[1],"%le",&offset);
    zx.in_mesh(argv[2],ftype::grid);
    for(int i=0;i<zx.nvrtx;++i)
       zx.vrtx[i][1] += offset;
    zx.out_mesh(argv[2],ftype::text);
    return 0;
#endif
   
#ifdef MAXMIN
   int count;
   double xpos;
   class mesh<2> zx;

   sscanf(argv[1],"data%d",&count);
   double x[3],y[3],xloc,ymax;
   double num = 0.0, den = 0.0;
   
   int sind,sind1;
   sind = zx.sbdry[0]->sd(0);
   for(int i=1;i<zx.sbdry[0]->nsd();++i) {
      sind = zx.sbdry[0]->sd(i);
#ifdef SKIP
      if(zx.vrtx[zx.svrtx[sind][1]][1] > zx.vrtx[zx.svrtx[sind][0]][1] &&
         zx.vrtx[zx.svrtx[sind1][1]][1] < zx.vrtx[zx.svrtx[sind1][0]][1]) {
            /* FIT QUADRATIC */
            x[0] = zx.vrtx[zx.svrtx[sind][0]][0];
            x[1] = zx.vrtx[zx.svrtx[sind][1]][0];
            x[2] = zx.vrtx[zx.svrtx[sind1][1]][0];
            y[0] = zx.vrtx[zx.svrtx[sind][0]][1];
            y[1] = zx.vrtx[zx.svrtx[sind][1]][1];
            y[2] = zx.vrtx[zx.svrtx[sind1][1]][1];
            
            if (x[1] <= 0.0001 || x[1] >= 5.9999) continue;
            
            y[0] /= (x[0]-x[1])*(x[0]-x[2]);
            y[1] /= (x[1]-x[0])*(x[1]-x[2]);
            y[2] /= (x[2]-x[0])*(x[2]-x[1]);
            
            xloc =  y[0]*(x[2]+x[1]) +y[1]*(x[2]+x[0]) +y[2]*(x[1]+x[0]);
            xloc /= 2.*(y[0] +y[1] +y[2]);
            
            ymax = y[0]*(xloc-x[1])*(xloc-x[2]) +y[1]*(xloc-x[0])*(xloc-x[2]) +y[2]*(xloc-x[0])*(xloc-x[1]);
                   
            printf("%e ",ymax);
      }
      sind = sind1;
#endif
      if ( zx.vrtx[zx.svrtx[sind][1]][1] > 0.1 && zx.vrtx[zx.svrtx[sind][1]][0] != zx.vrtx[zx.svrtx[sind][0]][0]) {
         
         
       //  xpos = zx.vrtx[zx.svrtx[sind][1]][0] -0.0581*2/5.0*count;
      //   xpos = zx.vrtx[zx.svrtx[sind][1]][0] -1.0/160.0*count;

         xpos = zx.vrtx[zx.svrtx[sind][1]][0] -1.053/40.0*count;

         //xpos = xpos -(int)((xpos-6.)/6.) *6.0;
         xpos = xpos -(int)((xpos-8.)/8.) *8.0;
         
      //   num += zx.vrtx[zx.svrtx[sind][1]][1]*sin(2.*M_PI/1.5*xpos);
      //   den += pow(sin(2.*M_PI/1.5*xpos),2.0);
         printf("%d %f %e\n",count,xpos,zx.vrtx[zx.svrtx[sind][1]][1]);
      }
      
   }
   
//   printf("%d %e\n",count,num/den);
  // printf("\n");

   return 0;
#endif
   
   
#ifdef CONVERT
   class mesh<2> zx;

   zx.in_mesh(argv[1],ftype::easymesh);
   zx.out_mesh(argv[1],ftype::tecplot);
   
   return(0);
#endif


class blocks z;

#ifdef GENERATE
   z.init(argv[1]);
   z.restructure();
   z.output("testing",ftype::grid);
   for (int i=0;i<10;++i)
      z.restructure();
   z.output("testing2",ftype::grid);
   
   return 0;
#endif

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
      
      input[0]["mglvls"] = "3";
      input[0]["ngrid"] = "3";
      input[0]["ntstep"] = "1";
      input[0]["fadd"] = "1.0";
      input[0]["vnn"] = "0.5";
      input[0]["itercrsn"] = "1";
      input[0]["iterrfne"] = "0";
      input[0]["njacobi"] = "1";
      input[0]["ncycle"] = "1500";
      input[0]["vwcycle"] = "2";
      input[0]["tolerance"] = "0.66";
      input[0]["logfile"] = "duh";
      
#ifdef MPISRC
      /* LIST OF NBLOCKS FOR EACH PROCESSOR */
      input[0]["nblock"] = "1 1";
      if (myid == 0) {
         input[1]["blktype"] = "3";
         input[1]["mesh"] = "/Users/helenbrk/Codes/grids/WIND/PRDC/top8";
         input[1]["growth factor"] = "4.0";
         input[1]["filetype"] = "0";
      }
      else {
         input[1]["blktype"] = "3";
         input[1]["mesh"] = "/Users/helenbrk/Codes/grids/WIND/PRDC/bot8";
         input[1]["growth factor"] = "4.0";
         input[1]["filetype"] = "0";
      }
#else
      input[0]["nblock"] = "2";
      input[1]["blktype"] = "3";
      input[1]["mesh"] = "/Users/helenbrk/Codes/grids/WIND/PRDC/top8";
      input[1]["growth factor"] = "4.0";
      input[1]["filetype"] = "0";
      input[2]["blktype"] = "3";
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
