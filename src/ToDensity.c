#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "Distribution.h"

#define HELPMESSAGE "\nusage: inter <file0> <file1> ...\n"


int main(int argc, char **argv) {		
	char option[256];
	FILE *fi, *fo;
	int i, j;
	for(i=0; i<256; i++)
		option[i] = 0; 
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['h']) {
			option['h'] = 0;
			printf("%s\n", HELPMESSAGE);
			exit(0);
		}
	}
	if(i<argc && (fi = fopen(argv[i], "r"))) {
		i++;
		TypeDistribution d = readDistribution(fi);
		if(i<argc && (fo = fopen(argv[i], "w"))) {
			deriveDistribution(&d);
			fprintDistribution(fo, d);
			fclose(fo);
		} else {
			fprintf(stderr, "Error while opening %s\n", argv[i]);
			exit(1);
		}
		fclose(fi);
	} else {
		fprintf(stderr, "Error while opening %s\n", argv[i]);
		exit(1);
	}
	return 0;
}
