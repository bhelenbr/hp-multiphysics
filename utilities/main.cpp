#include "utilities.h"
#include "myblas.h"
#include <stdio.h>
#include <new>
#include <fstream>

int main (int argc, const char * argv[]) {
	char out[100];
	
	number_str(out,"name",5,3);
	
	printf("%s\n",out);
	
	sharedmem mymem1(10*sizeof(int));
	sharedmem mymem2(10*sizeof(int));
	
	mymem1.checkout();
	
	int *myints = new (mymem1.data()) int;
	
	for(int i=0; i<10; ++i)
		myints[i] = i;
	
	for(int i=0;i<10;++i)
		printf("%d %d\n",i,myints[i]);
	
	mymem1.checkin();
	
	mymem2.reference(mymem1);
	mymem2.checkout();
	myints = new (mymem2.data()) int;
	
	for(int i=0;i<10;++i)
		printf("%d %d\n",i,myints[i]);
	
	mymem2.checkin();
	
	
	
	Array<int,1> ncols(5);
	ncols = 3, 2, 4, 2, 3;
	sparse_row_major M(5,ncols);
	M.set_values(0,0,0.0);
	M.set_values(0,1,2.0);
	M.add_values(0,1,-1.);
	M(0,4) = 3.0;
	
	M.set_values(1,0,0.0);
	M.set_values(1,1,2.0);
	
	Array<FLT,2> K(2,2);
	K = 1.0,2.0,
			2.0,3.0;
	Array<int,1> rows(2);
	rows = 2,3;
	Array<int,1> cols(2);
	cols = 3,4;
	M.set_values(2,rows,2,cols,K);
	
	
	M(2,0) = 0.5;
	M(2,1) = -0.5;
	
	Array<int,1> cols2(2);
	cols2 = 0,1;
	Array<double,1> vals(2);
	vals = -1., -1.;
	M.set_values(4,3,cols2,vals);
	M.match_patterns(0,4);
	
	std::cout << M << std::endl;
	
	M.swap_rows(1, 3);
	
	K = 0.5,0.5,
			1.0,0.0;
	rows = 0,4;
	cols = 0,4;
	M.combine_rows(2, rows, 2, cols, K);
	std::cout << M << std::endl;

	
	
	return(0);
	
}
