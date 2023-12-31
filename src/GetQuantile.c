#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "Distribution.h"

#define HELPMESSAGE "\nusage: inter <file0> <file1> ...\n"


int main(int argc, char **argv) {		
	char option[256];
	FILE *fi;
	int i, j;
	double conf = 0.95;
	for(i=0; i<256; i++)
		option[i] = 0; 
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['c']) {
			option['c'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &conf) == 1)
				i++;
			else {
				fprintf(stderr, "a number is required after option -c");
				exit(1);
			}
		}
		if(option['h']) {
			option['h'] = 0;
			printf("%s\n", HELPMESSAGE);
			exit(0);
		}
	}
	if(i<argc && (fi = fopen(argv[i], "r"))) {
		i++;
		double inf, sup, mode, max;
		TypeDistribution d = readDistribution(fi);
printf("sum %le\n", sumDistribution(d));
		if(sumDistribution(d) <= 1.000001) {
			max = d.item[0].dens;
			mode = d.item[0].val;
			for(j=1; j<d.size; j++)
				if(d.item[j].dens>max) {
					max = d.item[j].dens;
					mode = d.item[j].val;
				}
			integrateDistribution(&d);
//fprintDistribution(stdout, d);
			inf = getQuantileInf(d, (1.-conf)/2.);
			sup = getQuantileSup(d, (1.-conf)/2.);
		} else {
			inf = getQuantileInf(d, (1.-conf)/2.);
			sup = getQuantileSup(d, (1.-conf)/2.);
			deriveDistribution(&d);
			max = d.item[0].dens;
			mode = d.item[0].val;
			for(j=1; j<d.size; j++)
				if(d.item[j].dens>max) {
					max = d.item[j].dens;
					mode = d.item[j].val;
				}
		}
		fprintf(stdout, "%.2lf %% confidence interval: [%.2lf, %.2lf] ;  mode: %.2lf\n", conf, inf, sup, mode);
	}
	return 0;
}
