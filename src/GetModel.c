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
//./bfbd -S -320 -S -272.3 x ./TreeEupely1.phy ~/Dropbox/FF/git/data/Extinction/FossilAges.csv

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

#define HELPMESSAGE "--------------------------\n\nNAME\n	dist - Computation of the divergence time distibution associated to a given clade from a set of possible trees, the fossil ages, and the diversification rates\n	\nSYNOPSIS\n	dist [OPTIONS] <input Tree File> <input Fossil File> <input List Clade> [output File]\n\nDESCRIPTION\n	Compute the distribution of the divergence time associated to the clade corresponding to the list of tips given in <input List Clade> by sampling into the trees contianed in <input Tree File> (in Newick format) into the fossil ranges provided in <input Fossil File> (in csv format) and output the distribution as a .csv table <output File>.csv\n\n	Options are\n	-z <input Tree File>\n		output the tree in text debug format in the console and exit \n	-o <origin bound inf> [origin bound sup]\n		set the origin range; if origin bound sup is not given it is set to the most ancient fossil age\n	-e <end bound inf> <end bound sup>\n		set the end time range\n	-p <speciation rate> <extinction rate> <fossilization rate>\n		set the speciation, extinction and fossilization rates\n	-s <number>\n		set the number of samples\n	-d\n		return the distribution (otherwise the density is returned by default)\n	-u <value>\n		set the step discretizing the time distribution to <value>\n	-s <number>\n		set the number of thread running in parallell\n	-h\n		display help\n\n--------------------------\n\n"


int main(int argc, char **argv) {	
	char *inputFileNameTree, *inputFileNameFossil = NULL, *inputFileNameScheme, outputFileName[STRING_SIZE], option[256], *outCoda = NULL;
	FILE *ft, *ff, *fs;
	int i, niter = 5, nSamp = 5000, nBurn = 10000, nGap = 10;
	double al = 0.75, probSpe = 0.25, probExt = 0.25, probFos = 0.25, propParam = 0.1, step = 0.1;
	TypeModelParam windSize = {.birth=0.05, .death = 0.05, .fossil = 0.025, .sampling = 1.}, init = {.birth=0.5, .death = 0.5, .fossil = 0.1, .sampling = 1.};
	
	for(i=0; i<256; i++)
		option[i] = 0;
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		int j;
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['z']) {
			FILE *fi;
			option['z'] = 0;
			if((i+1)<argc && (fi = fopen(argv[++i], "r"))) {
				TypeTree *tree;
				tree = readTree(fi);
				toBinary(tree);
				printTreeDebug(stdout, tree->root, tree, tree->name);
				reorderTree(tree->name, tree);
				if(tree->minTime == NO_TIME || tree->minTime == 0.)
					tree->minTime = tree->time[tree->root]*0.9;
				if(tree->maxTime == NO_TIME) {
					int n;
					tree->maxTime = 0.;
					for(n=0; n<tree->size; n++)
						if(tree->time[n]>tree->maxTime)
							tree->maxTime = tree->time[n];
				}
				printTreeDebug(stdout, tree->root, tree, tree->name);
			} else {
				fprintf(stderr, "Error while reading %s.\n", argv[i]);
				exit(1);
			}
			exit(0);
		}
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
		if(option['x']) {
			option['x'] = 0;
			if((i+1)<argc)
				outCoda = argv[++i];
			else
				error("a file name is expected after -u");
		}
		if(option['u']) {
			option['u'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &niter) == 1)
				i++;
		}
		if(option['s']) {
			option['s'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &nSamp) == 1)
				i++;
			else
				error("a number is expected after -s");
		}
		if(option['u']) {
			option['u'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &(step)) == 1)
				i++;
			else
				error("a number is expected after -u");
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
		printSamplingScheme(stdout, scheme);

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
fprintTreeNewick(stdout, tree[0]);
//exit(0);
		my_rand_init();
		rand_data = my_rand_get_data();
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
				fprintf(fo, "library(coda)\nlibrary(tikzDevice)\nsetwd(\"~/Dropbox/FossilPiecewise/dev/src0.9/\")\n");
				fprintf(fo, "A = read.coda(\"%s.out\", \"%s.ind\")\n", outCoda, outCoda);
				for(n=0; n<scheme.sizeParamBirth; n++) {
					fprintf(fo, "tikz('%s_birth_%d.tex', standAlone = TRUE, width=5, height=5)\nhist(c(as.matrix(A)[,\"birth_%d\"]", outCoda, n, n);
					fprintf(fo, "), 45,axes=T, xlab=\"speciation rate %d\", main=NULL, yaxt='n', ylab = \"\")\ndev.off()\n\n", n+1);
				}
				for(n=0; n<scheme.sizeParamDeath; n++) {
					fprintf(fo, "tikz('%s_death_%d.tex', standAlone = TRUE, width=5, height=5)\nhist(c(as.matrix(A)[,\"death_%d\"]", outCoda, n, n);
					fprintf(fo, "), 45,axes=T, xlab=\"extinction rate %d\", main=NULL, yaxt='n', ylab = \"\")\ndev.off()\n\n", n+1);
				}
				for(n=0; n<scheme.sizeParamFossil; n++) {
					fprintf(fo, "tikz('%s_fossil_%d.tex', standAlone = TRUE, width=5, height=5)\nhist(c(as.matrix(A)[,\"fossil_%d\"]", outCoda, n, n);
					fprintf(fo, "), 45,axes=T, xlab=\"fossilization rate %d\", main=NULL, yaxt='n', ylab = \"\")\ndev.off()\n\n", n+1);
				}
				for(n=0; n<scheme.sizeParamSampling; n++) {
					fprintf(fo, "tikz('%s_sampling_%d.tex', standAlone = TRUE, width=5, height=5)\nhist(c(as.matrix(A)[,\"sampling_%d\"]", outCoda, n, n);
					fprintf(fo, "), 45,axes=T, xlab=\"sampling probability %d\", main=NULL, yaxt='n', ylab = \"\")\ndev.off()\n\n", n+1);
				}
				
				//fprintf(fo, "library(coda)\nlibrary(tikzDevice)\nsetwd(\"~/Dropbox/FossilPiecewise/dev/src0.9/\")\n");
				//fprintf(fo, "A = read.coda(\"%s.out\", \"%s.ind\")\n", outCoda, outCoda);
				//for(n=0; n<scheme.sizeParamBirth; n++) {
					//fprintf(fo, "\nhist(c(as.matrix(A)[,\"birth_%d\"]", n);
					//fprintf(fo, "), 45,axes=T, xlab=\"speciation rate %d\", main=NULL, yaxt='n', ylab = \"\")\n\n", n+1);
				//}
				//for(n=0; n<scheme.sizeParamDeath; n++) {
					//fprintf(fo, "\nhist(c(as.matrix(A)[,\"death_%d\"]", n);
					//fprintf(fo, "), 45,axes=T, xlab=\"extinction rate %d\", main=NULL, yaxt='n', ylab = \"\")\n\n", n+1);
				//}
				//for(n=0; n<scheme.sizeParamFossil; n++) {
					//fprintf(fo, "hist(c(as.matrix(A)[,\"fossil_%d\"]", n);
					//fprintf(fo, "), 45,axes=T, xlab=\"fossilization rate %d\", main=NULL, yaxt='n', ylab = \"\")\n\n", n+1);
				//}
				//for(n=0; n<scheme.sizeParamSampling; n++) {
					//fprintf(fo, "\nhist(c(as.matrix(A)[,\"sampling_%d\"]", n);
					//fprintf(fo, "), 45,axes=T, xlab=\"sampling probability %d\", main=NULL, yaxt='n', ylab = \"\")\n\n", n+1);
				//}
				fclose(fo);
			}
		}
	} else
		error("Cannot read %s, %s or %s\n", inputFileNameTree, inputFileNameFossil, inputFileNameScheme);
	return 0;
}
