#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include "Utils.h"
#include "MyR.h"
#include "MyRandom.h"
#include "Tree.h"
#include "Fossil.h"
#include "FossilInt.h"
#include "PiecewiseModel.h"
#include "FBDDensity.h"
#include "SimulPiecewise.h"
#include "MCMCImportanceSamplingFossilInt.h"


//./xfbd -u 40 -p -20. 0.5 0.1 0.5 1. -p -10 0.05 0.2 0.1 1. -s 5000 -b 10000 -g 50 -f 1. -a 0.33 0.33 -w 0.5 0.5 0.5 -i 5. 5. 5.  -r 0000 
//valgrind ./dfbd -s 200 ~/Dropbox/FossilPiecewise/dataMS/TreeEupely2.phy ~/Dropbox/FossilPiecewise/dataMS/CotylosauriaAgesNew.csv ~/Dropbox/FossilPiecewise/modelMS/schemeM0.txt

#define NAME_CODA "mcmc_sample"

#define MAX_PARAM 10

#define STRING_SIZE 300
#define HELP_MESSAGE "\nusage: sample [options] [<output file>]\n\nsample simulates random trees and fossils finds and saves them in Newick format\n\nOptions:\n\t-h : display help\n\t-b <birth>\t: set birth rate\n\t-d <death>\t: set death rate\n\t-f <fossil>\t: set fossil find rate\n\t-m <min>\t: set minimum number of contemporary species of a simulation to be considered\n\t-M <size>\t: set maximum size of a simulation to be considered\n\t-i <niter>\t: set the number of simulations\n\t-t <time> : the end time of the diversification (start is always 0)\n"

int main(int argc, char **argv) {	
	char *inputFileNameTree, outputName[STRING_SIZE], outputFileName[STRING_SIZE], *outputPrefix, option[256];
	double endTime;
	FILE *fit, *fo;
	int i;

	for(i=0; i<256; i++)
		option[i] = 0;
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		int j;
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['h']) {
			printf("%s\n", HELP_MESSAGE);
			return 0;
		}
	}
	if(i<argc) {
		inputFileNameTree = argv[i++];
	} else {
		fprintf(stderr, "Please provide the name of a file containing a phylogenetic tree in Newick format\n");
		exit(1);
	}
	if(i<argc) {
		endTime = atof(argv[i++]);
	} else {
		fprintf(stderr, "Please provide the ending time\n");
		exit(1);
	}
	strcpy(outputName, inputFileNameTree);
	if((outputPrefix = strrchr(outputName, '.')) != NULL)
		outputPrefix[0] = '\0';
	if((outputPrefix=strrchr(outputName, '/')) == NULL)
		outputPrefix = outputName;
	else
		outputPrefix++;
	
	if((fit = fopen(inputFileNameTree, "r"))) {
		TypeTree **tree;
		tree = readTrees(fit);
        fclose(fit);
		for(i=0; tree[i]!=NULL; i++) {
			int n;
			toBinary(tree[i]);
			if(tree[i]->name!=NULL)
				for(n=0; n<tree[i]->size; n++)
					if(tree[i]->name[n]!=NULL)
						fixSpace(tree[i]->name[n]);
			for(n=0; n<tree[i]->size; n++)
				if(tree[i]->time[n] < endTime)
					tree[i]->time[n] = NO_TIME;
				else
					tree[i]->time[n] = endTime;
		}
		sprintf(outputFileName, "%s_no_time.newick", outputName);
		if((fo = fopen(outputFileName, "w"))) {
			for(i=0; tree[i]!=NULL; i++) {
				fprintTreeNewick(fo, tree[i]);
				fprintf(fo, "\n");
			}
		}
	}
	return 0;
}
