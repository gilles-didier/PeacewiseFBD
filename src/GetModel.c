#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "Utils.h"
#include "MyR.h"
#include "MyRandom.h"
#include "Tree.h"
#include "SimulPiecewise.h"
#include "Fossil.h"
#include "FossilInt.h"
#include "PiecewiseModel.h"
#include "MCMCImportanceSamplingFossilInt.h"

#ifdef DO_PS
#endif

#define STRING_SIZE 300
#define INC_CUT 5
#define SEQ_SIZE 30
#define EXT_OUTPUT "_added.phy"
#define MAX_PRED 7.5E+07

#define SIZE_BUFFER_CHAR 300
#define INFTY 1E99
#define RINFTY 1E99
#define DEF 10
#define MIN_VAL 0.000001
#define DELTA 0.000001

#define MINVAL 0.01
#define TRIAL 10
#define FIX_VAL(x) (((x)<=0)?MINVAL:(x))

#define MAX_ITER 1000

#define M_MAX 6
#define M_MAX_F 4
#define MIN_TREE 20
#define PREFIX "table"
#define MAX_TRIALS 1000
#define MAX_PARAM 20

#define HELPMESSAGE "\n--------------------------\n\nNAME\n\n	getPost - provides coda files for MCMC convergence diagnostics and posterior distribution of the model specification parameters\n	\nSYNOPSIS\n\n	getPost [OPTIONS] <tree(s)> <fossil ages> <model specification> [output File]\n\nDESCRIPTION\n\n	Provides coda files for MCMC convergence diagnostics and posterior distribution of the model specification parameters of  <model specification> from the tree(s) in file <tree(s)> (in Newick format) with the fossil ranges provided in <fossil ages> (in csv format). The two files <output File>.out and <output File>.ind are to be read by the R package coda.\n\n	Options are\n	-w <speciation rate width> <extinction rate width> <fossilization rate width> <sampling probability width>\n		set the widths of the sliding windows used for sampling the speciation, extinction and fossilization rates and of the sampling probability during the MCMC\n	-i <speciation rate upper bound> <extinction rate  upper bound> <fossilization rate upper bound> <sampling probability upper bound>\n		set the upper bounds of the interval in which the speciation, extinction and fossilization rates and the sampling probability are uniformly drawn at the beginning of the MCMC (lower bounds are all 0)\n	-b <number>\n		set the number of iterations for the burning stage\n	-g <number>\n		set the thinning parameter (only one iteration in <number> is considered)  \n	-f <proportion>\n		set the proportion of moves in the parameters in the MCMC proposal (other moves are drawn uniformly amont the fossil ages)\n	-a <speciation proportion> <extinction proportion> <fossilization proportion>\n		set the relation proportion of moves in the speciation, the extinction and the fossilzation rates (thus in the sampling probability)\n	-s <number>\n		set the number of samples required\n	-r <number>\n		set the random seed\n	-h\n		display help\n\nEXAMPLE\n\n./getPost -s 5000 -b 50000 -g 50 -f 0.2 -a 0.25 0.25 0.25 -w 0.05 0.05 0.05 0.25 -i 0.5 0.5 0.5 1.  ../data/Simulated_dataset_tree.newick ../data/Simulated_dataset_fossils.csv ../data/Simulated_Dataset_model_spec/Simul_model_spec_0.txt"


int main(int argc, char **argv) {	
	char *inputFileNameTree, *inputFileNameFossil = NULL, *inputFileNameScheme, outputFileName[STRING_SIZE], option[256], *outCoda = NULL;
	FILE *ft, *ff, *fs;
	int i, nSamp = 5000, nBurn = 10000, nGap = 10;
	unsigned long int seed = 0;
	double al = 0.75, probSpe = 0.25, probExt = 0.25, probFos = 0.25, propParam = 0.1;
	TypeModelParam windSize = {.birth=0.05, .death = 0.05, .fossil = 0.025, .sampling = 1.}, init = {.birth=0.5, .death = 0.5, .fossil = 0.1, .sampling = 1.};
	
	for(i=0; i<256; i++)
		option[i] = 0;
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		int j;
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['a']) {
			option['a'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &probSpe) == 1)
				i++;
			else
				error("3 values are expected after -a");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &probExt) == 1)
				i++;
			else
				error("3 values are expected after -a");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &probFos) == 1)
				i++;
			else
				error("3 values are expected after -a");
		}
		if(option['w']) {
			option['w'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &windSize.birth) == 1)
				i++;
			else
				error("4 values are expected after -w");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &windSize.death) == 1)
				i++;
			else
				error("4 values are expected after -w");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &windSize.fossil) == 1)
				i++;
			else
				error("4 values are expected after -w");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &windSize.sampling) == 1)
				i++;
			else
				error("4 values are expected after -w");
		}
		if(option['i']) {
			option['i'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &init.birth) == 1)
				i++;
			else
				error("4 values are expected after -w");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &init.death) == 1)
				i++;
			else
				error("4 values are expected after -w");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &init.fossil) == 1)
				i++;
			else
				error("4 values are expected after -w");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &init.sampling) == 1)
				i++;
			else
				error("4 values are expected after -w");
		}
		if(option['s']) {
			option['s'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &nSamp) == 1)
				i++;
			else
				error("a number is expected after -s");
		}
		if(option['r']) {
			option['r'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lu", &seed) == 1)
				i++;
			else
				error("a number is expected after -r");
		}
		if(option['b']) {
			option['b'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &nBurn) == 1)
				i++;
			else
				error("a number is expected after -b");
		}
		if(option['g']) {
			option['g'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &nGap) == 1)
				i++;
			else
				error("a number is expected after -b");
		}
		if(option['f']) {
			option['f'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &propParam) == 1)
				i++;
			else
				error("a number is required after option -t");
		}
		if(option['s']) {
			option['s'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &nSamp) == 1)
				i++;
			else
				error("a number is expected after -s");
		}
		if(option['h']) {
			option['h'] = 0;
			printf("%s\n", HELPMESSAGE);
			exit(0);
		}
	}
	if(i<argc) {
		inputFileNameTree = argv[i++];
	} else {
		fprintf(stderr, "Please provide the name of a file containing a phylogenetic tree in Newick format\n");
		exit(1);
	}
	if(i<argc) {
		inputFileNameFossil = argv[i++];
	} else {
		fprintf(stderr, "Please provide the name of a file containing a fossil list\n");
		exit(1);
	}
	if(i<argc) {
		inputFileNameScheme = argv[i++];
	} else
		error("Please provide the name of a file containing a  model scheme\n");
	if(i<argc)
		outCoda = argv[i++];
	else
		outCoda = "CodaTmp";
	if((ft = fopen(inputFileNameTree, "r")) && (ff = fopen(inputFileNameFossil, "r")) && (fs = fopen(inputFileNameScheme, "r"))) {
		TypeTree **tree;
		TypeFossilIntFeature *fint;
		TypeSamplingScheme scheme;
		FILE *find, *fout, *fo;
		char nameInd[STRING_SIZE], nameOut[STRING_SIZE];
		void *rand_data;

		int i, n, sizeTree;		
 		scheme = readSamplingScheme(fs);
		fclose(fs);

       tree = readTrees(ft);
        fclose(ft);
        if(tree[0] == NULL) {
			fprintf(stderr, "Error: no tree\n");
			return 1;
		}
		sizeTree = 0;
		for(i=0; tree[i]!=NULL; i++) {
			int n;
			toBinary(tree[i]);
			if(tree[i]->name!=NULL)
				for(n=0; n<tree[i]->size; n++)
					if(tree[i]->name[n]!=NULL)
						fixSpace(tree[i]->name[n]);
			tree[i]->minTime = scheme.param.startTime[0];
			tree[i]->maxTime = scheme.param.startTime[scheme.param.size];
			sizeTree++;
		}
		if(sizeTree == 0)
			error("No tree!\n");
		reindexTreesFromName(tree, sizeTree);
		fint = getFossilIntFeature(ff, tree[0]->name, tree[0]->size);
		fclose(ff);
		if(getMaxFossilIntTime(fint) > 0.)
			negateFossilInt(fint);
		fixStatus(tree[0], fint);
		for(n=0; n<tree[0]->size; n++)
			if(fint->status[n] == unknownNodeStatus)
				for(i=0; tree[i]!=NULL; i++)
					tree[i]->time[n] = fint->endTimeTable[n].inf;
printf("\noutput files: %s.ind, %s.out, %s.R\n\n", outCoda, outCoda, outCoda);
		my_rand_init();
		rand_data = my_rand_get_data();
		my_rand_set_seed(rand_data, seed);
		sprintf(nameInd, "%s.ind", outCoda);
		sprintf(nameOut, "%s.out", outCoda);
		if((fout = fopen(nameOut, "w")) && (find = fopen(nameInd, "w"))) {
			MCMCSamplingSamplePosteriorParametersOnly(fout, find, tree, sizeTree, fint, scheme, al, nBurn, nGap, nSamp, propParam, &windSize, &init, probSpe,probExt, probFos, rand_data);
		}
		my_rand_free_data(rand_data);
		if(outCoda != NULL) {
			sprintf(outputFileName, "%s.R", outCoda);
			if((fo = fopen(outputFileName, "w"))) {
				int n;
				fprintf(fo, "library(coda)\nlibrary(tikzDevice)\nsetwd(\"directory of the results\")\n");
				fprintf(fo, "A = read.coda(\"%s.out\", \"%s.ind\")\n", outCoda, outCoda);
				for(n=0; n<scheme.sizeParamBirth; n++) {
					fprintf(fo, "\nhist(c(as.matrix(A)[,\"birth_%d\"]", n);
					fprintf(fo, "), 45,axes=T, xlab=\"speciation rate %d\", main=NULL, yaxt='n', ylab = \"\")\n\n", n+1);
				}
				for(n=0; n<scheme.sizeParamDeath; n++) {
					fprintf(fo, "\nhist(c(as.matrix(A)[,\"death_%d\"]", n);
					fprintf(fo, "), 45,axes=T, xlab=\"extinction rate %d\", main=NULL, yaxt='n', ylab = \"\")\n\n", n+1);
				}
				for(n=0; n<scheme.sizeParamFossil; n++) {
					fprintf(fo, "hist(c(as.matrix(A)[,\"fossil_%d\"]", n);
					fprintf(fo, "), 45,axes=T, xlab=\"fossilization rate %d\", main=NULL, yaxt='n', ylab = \"\")\n\n", n+1);
				}
				for(n=0; n<scheme.sizeParamSampling; n++) {
					fprintf(fo, "\nhist(c(as.matrix(A)[,\"sampling_%d\"]", n);
					fprintf(fo, "), 45,axes=T, xlab=\"sampling probability %d\", main=NULL, yaxt='n', ylab = \"\")\n\n", n+1);
				}
				fclose(fo);
			}
		}
	} else
		error("Cannot read %s, %s or %s\n", inputFileNameTree, inputFileNameFossil, inputFileNameScheme);
	return 0;
}
