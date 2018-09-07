#include<stdio.h>
#include<math.h>

#define FSRF_MASK (1<<0)
#define IFCE_MASK (1<<1)
#define INFL_MASK (1<<2)
#define OUTF_MASK (1<<3)
#define SYMM_MASK (1<<4)
#define EULR_MASK (1<<5)
#define PRDX_MASK (1<<6)
#define PRDY_MASK (1<<7)
#define COMX_MASK (1<<8)
#define COMY_MASK (1<<8)
#define CURV_MASK (1<<9)

main() {
	int i,j,k,l,cseg;
	double dens1,dens2,dens3,x,y;
	double front,back,width;
	
	cseg = 4;
	dens1 = 2.0/10.0;

	front = -1.5;
	back = 8.5;
	width = 1.0;

	k = 0;
	printf("%d #NUMBER OF POINTS#\n\n",4 +cseg);

	printf("#BOUNDARY POINTS#\n");	
	printf("%d:\t%lf\t%lf\t%lf\t%d\n",k,front,-width,dens1,INFL_MASK);
	++k;
	printf("%d:\t%lf\t%lf\t%lf\t%d\n",k,back,-width,dens1,OUTF_MASK);
	++k;
	printf("%d:\t%lf\t%lf\t%lf\t%d\n",k,back,width,dens1,OUTF_MASK);	
	++k;
	printf("%d:\t%lf\t%lf\t%lf\t%d\n",k,front,width,dens1,INFL_MASK);
	++k;

	printf("\n#CIRCLE POINTS#\n");		
	for(i=0;i<cseg;++i) {
		x = 0.50*cos(-(i+.5)*2*M_PI/cseg);
		y = 0.50*sin(-(i+.5)*2*M_PI/cseg);
		printf("%d:\t%lf\t%lf\t%lf\t%d\n",k,x,y,100.0,EULR_MASK);
		++k;
	}	
	k = 0;
	
	printf("\n%d #NUMBER OF SEGMENTS#\n\n",4+cseg);
	printf("\n#BOUNDARY SEGMENTS#\n");	
	printf("%d:\t%d\t%d\t%d\n",k,k,k+1,OUTF_MASK);
	++k;
	printf("%d:\t%d\t%d\t%d\n",k,k,k+1,OUTF_MASK+(1<<16));
	++k;
	printf("%d:\t%d\t%d\t%d\n",k,k,k+1,OUTF_MASK+(2<<16));
	++k;
	printf("%d:\t%d\t%d\t%d\n",k,k,0,INFL_MASK);
	++k;
	j = k;			
	printf("\n#CIRCLE SEGMENTS#\n");				
	for(i=0;i<cseg-1;++i)  {
		printf("%d:\t%d\t%d\t%d\n",k,k,k+1,EULR_MASK);
		++k;
	}
	printf("%d:\t%d\t%d\t%d\n",k,k,j,EULR_MASK);	
	++k;
}
