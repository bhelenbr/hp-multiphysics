#include <stdio.h>

int main(int argc, char **argv) {
	char buffer[1000],buffer1[1000];
        FILE *fp1, *fp2;

	fp1 = fopen(argv[1],"r");
	fp2 = fopen(argv[2],"r");
		
	while (fgets(buffer,1000,fp1)) {
		sscanf(buffer,"%[^\n]\n",buffer1);
		printf("%s ",buffer1);
		fgets(buffer,1000,fp2);
		printf("%s",buffer);
	}
	return 0;
}
