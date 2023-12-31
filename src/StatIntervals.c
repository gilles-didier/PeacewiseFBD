#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>



#include "Utils.h"
#include "FossilInt.h"
#ifdef DO_PS
#endif


#define PREFIX "table.csv"


#define HELPMESSAGE "--------------------------\n\nNAME\n	dist - Computation of the divergence time distibution associated to a given clade from a set of possible trees, the fossil ages, and the diversification rates\n	\nSYNOPSIS\n	dist [OPTIONS] <input Tree File> <input Fossil File> <input List Clade> [output File]\n\nDESCRIPTION\n	Compute the distribution of the divergence time associated to the clade corresponding to the list of tips given in <input List Clade> by sampling into the trees contianed in <input Tree File> (in Newick format) into the fossil ranges provided in <input Fossil File> (in csv format) and output the distribution as a .csv table <output File>.csv\n\n	Options are\n	-z <input Tree File>\n		output the tree in text debug format in the console and exit \n	-o <origin bound inf> [origin bound sup]\n		set the origin range; if origin bound sup is not given it is set to the most ancient fossil age\n	-e <end bound inf> <end bound sup>\n		set the end time range\n	-p <speciation rate> <extinction rate> <fossilization rate>\n		set the speciation, extinction and fossilization rates\n	-s <number>\n		set the number of samples\n	-d\n		return the distribution (otherwise the density is returned by default)\n	-u <value>\n		set the step discretizing the time distribution to <value>\n	-s <number>\n		set the number of thread running in parallell\n	-h\n		display help\n\n--------------------------\n\n"


int main(int argc, char **argv) {	
	char *inputFileNameFossil = NULL, *outputFileName = PREFIX, option[256];
	FILE *fo, *ff;
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
	if(i<argc) {
		inputFileNameFossil = argv[i++];
	} else {
		fprintf(stderr, "Please provide the name of a file containing a fossil list\n");
		exit(1);
	}
	if(i<argc)
		outputFileName = argv[i++];
	if((ff = fopen(inputFileNameFossil, "r"))) {
		TypeNameFossilIntTab *tab;
		double mean = 0.;
		int tot = 0;
		tab = readFossilIntTab(ff);
		fclose(ff);
		if((fo = fopen(outputFileName, "w"))) {
			for(i=0; i<tab->size; i++) {
				tot += tab->fossilIntTab[i].size;
				for(j=0; j<tab->fossilIntTab[i].size; j++) {
					fprintf(fo, "%lf\n", fabs(tab->fossilIntTab[i].fossilInt[j].sup-tab->fossilIntTab[i].fossilInt[j].inf));
					fprintf(stdout, "i %d, j %d %lf\n", i, j, fabs(tab->fossilIntTab[i].fossilInt[j].sup-tab->fossilIntTab[i].fossilInt[j].inf));
					mean += fabs(tab->fossilIntTab[i].fossilInt[j].sup-tab->fossilIntTab[i].fossilInt[j].inf);
				}
			}
			fprintf(fo, "%d fossils\nmean %lf\n", tot, mean/((double) tot));
			fprintf(stdout, "%d fossils\nmean %lf\n", tot, mean/((double) tot));
			fclose(fo);
		}
			
		freeNameFossilIntTab(tab);
	} else {
		fprintf(stderr, "Cannot read %s\n", inputFileNameFossil);
		exit(1);
	}
	return 0;
}
