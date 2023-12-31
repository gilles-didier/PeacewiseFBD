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

#define NAME_CODA "mcmc_sample"

#define MAX_PARAM 10

#define STRING_SIZE 300
#define HELP_MESSAGE "\nusage: sample [options] [<output file>]\n\nsample simulates random trees and fossils finds and saves them in Newick format\n\nOptions:\n\t-h : display help\n\t-b <birth>\t: set birth rate\n\t-d <death>\t: set death rate\n\t-f <fossil>\t: set fossil find rate\n\t-m <min>\t: set minimum number of contemporary species of a simulation to be considered\n\t-M <size>\t: set maximum size of a simulation to be considered\n\t-i <niter>\t: set the number of simulations\n\t-t <time> : the end time of the diversification (start is always 0)\n"

//./afbd -p -350. 0.03 0.01 0.02 0.2 -p -250 0.1 0.01 0.02 1. -p -200 0.01 0.01 0.02 1. -E 0. -f 5 -b 10000 -s 1000  -u 40 -C schemeDatasetPresent.txt 

//./afbd -p -350. 0.03 0.01 0.02 0.2 -p -250 0.8 0.01 0.02 1. -p -200 0.01 0.01 0.02 1. -E 0. -f 5 -b 10000 -s 1000  -u 40 -C schemeDatasetPresent.txt 
//./samp -b 2 -d 1 -f 1 -t 5 -i 4 simul.newick

void nameLeavesRec(int n, int *ind, int length, TypeTree *tree);
void renameLeaves(TypeTree *tree);
TypeTree *getRandomTreeFossil(TypePiecewiseModelParam *param, int minSizeTree, int maxSizeTree, int minSlice, double stopTime, void *rand_data);
TypeTree *getRandomTreeFossilInt(TypePiecewiseModelParam *param, int minSizeTree, int maxSizeTree, int minSlice, double stopTime, double fos_width, void *rand_data);


int main(int argc, char **argv) {	
	char *inputFileNameModel, *inputFileNameScheme, outputName[STRING_SIZE], outputFileName[STRING_SIZE], *outputPrefix, option[256], *outCoda = NULL;
	FILE *fim, *fis;
	int i, niter = 5, maxSizeTree = 1000, minSizeTree = 100, minSlice= 20, nSamp = 5000, nBurn = 10000, nGap = 10;
	double al = 0.75, probSpe = 0.25, probExt = 0.25, probFos = 0.25, propParam = 0.2, fos_width = 1., stopTime = DBL_MAX;
	TypeModelParam windSize = {.birth=0.02, .death = 0.02, .fossil = 0.02, .sampling = 0.2}, init = {.birth=0.5, .death = 0.5, .fossil = 0.1, .sampling = 1.};

	outCoda = NAME_CODA;
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
		if(option['m']) {
			option['m'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &minSizeTree) == 1)
				i++;
		}
		if(option['M']) {
			option['M'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &maxSizeTree) == 1)
				i++;
		}
		if(option['l']) {
			option['l'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &minSlice) == 1)
				i++;
		}
		if(option['f']) {
			option['f'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &fos_width) == 1)
				i++;
		}
		if(option['S']) {
			option['S'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &stopTime) == 1)
				i++;
		}
		if(option['h']) {
			printf("%s\n", HELP_MESSAGE);
			return 0;
		}
	}
	if(i<argc) {
		inputFileNameModel = argv[i++];
	} else
		error("Please provide the name of a file containing a piecewise model\n");
	if(i<argc) {
		inputFileNameScheme = argv[i++];
	} else
		error("Please provide the name of a file containing a  model scheme\n");
	if(!(i<argc && sscanf(argv[i++], "%s", outputName) == 1)) {
		strcpy(outputName, "Tree_simul");
	}
	if((outputPrefix = strrchr(outputName, '.')) != NULL)
		outputPrefix[0] = '\0';
	if((outputPrefix=strrchr(outputName, '/')) == NULL)
		outputPrefix = outputName;
	else
		outputPrefix++;
	outCoda = outputPrefix;

	if((fim = fopen(inputFileNameModel, "r")) && (fis = fopen(inputFileNameScheme, "r"))) {
		FILE *fo;
		TypePiecewiseModelParam param;
		TypeSamplingScheme scheme;
		param = readPiecewiseModelParam(fim);
		fclose(fim);
		printPiecewiseModel(stdout, &param);
		scheme = readSamplingScheme(fis);
		fclose(fis);
		printSamplingScheme(stdout, scheme);
		my_rand_init();
		#pragma omp parallel
		{
			int i;
			#pragma omp for
			for (i=0; i<niter; i++)	{
				void *rand_data;
				TypeTree *tree;
				char nameInd[STRING_SIZE], nameOut[STRING_SIZE];
				FILE *find, *fout;

				rand_data = my_rand_get_data();
				tree = getRandomTreeFossilInt(&param, minSizeTree, maxSizeTree, minSlice, stopTime, fos_width, rand_data);
				tree->name = nameAll("Tip_", "Int_", tree);
				sprintf(nameInd, "%s_%d.ind", outCoda, i);
				sprintf(nameOut, "%s_%d.out", outCoda, i);
				if((fout = fopen(nameOut, "w")) && (find = fopen(nameInd, "w"))) {
					MCMCSamplingSamplePosteriorParametersOnly(fout, find, &tree, 1, (TypeFossilIntFeature*) tree->info, scheme, al, nBurn, nGap, nSamp, propParam, &windSize, &init, probSpe,probExt, probFos, rand_data);
				}
				my_rand_free_data(rand_data);
			}
		}
		if(outCoda != NULL) {
			sprintf(outputFileName, "%s.R", outCoda);
			if((fo = fopen(outputFileName, "w"))) {
				int n;
				fprintf(fo, "library(coda)\nlibrary(tikzDevice)\nsetwd(\"~/Dropbox/FossilPiecewise/dev/src0.6/\")\n");
				for(i=0; i<niter; i++)
					fprintf(fo, "A%d = read.coda(\"%s_%d.out\", \"%s_%d.ind\")\n", i, outCoda, i, outCoda, i);
				for(n=0; n<scheme.sizeParamBirth; n++) {
					fprintf(fo, "tikz('%s_birth_%d.tex', standAlone = TRUE, width=5, height=5)\nhist(c(as.matrix(A0)[,\"birth_%d\"]", outCoda, n, n);
					for(i=1; i<niter; i++)
						fprintf(fo, ", as.matrix(A%d)[,\"birth_%d\"]", i, n);
					fprintf(fo, "), 45,axes=T, xlab=\"speciation rate %d (%.3lf)\", main=NULL, yaxt='n', ylab = \"\")\ndev.off()\n\n", n+1, param.param[n].birth);
				}
				for(n=0; n<scheme.sizeParamDeath; n++) {
					fprintf(fo, "tikz('%s_death_%d.tex', standAlone = TRUE, width=5, height=5)\nhist(c(as.matrix(A0)[,\"death_%d\"]", outCoda, n, n);
					for(i=1; i<niter; i++)
						fprintf(fo, ", as.matrix(A%d)[,\"death_%d\"]", i, n);
					fprintf(fo, "), 45,axes=T, xlab=\"extinction rate %d (%.3lf)\", main=NULL, yaxt='n', ylab = \"\")\ndev.off()\n\n", n+1, param.param[n].death);
				}
				for(n=0; n<scheme.sizeParamFossil; n++) {
					fprintf(fo, "tikz('%s_fossil_%d.tex', standAlone = TRUE, width=5, height=5)\nhist(c(as.matrix(A0)[,\"fossil_%d\"]", outCoda, n, n);
					for(i=1; i<niter; i++)
						fprintf(fo, ", as.matrix(A%d)[,\"fossil_%d\"]", i, n);
					fprintf(fo, "), 45,axes=T, xlab=\"fossilization rate %d (%.3lf)\", main=NULL, yaxt='n', ylab = \"\")\ndev.off()\n\n", n+1, param.param[n].fossil);
				}
				for(n=0; n<scheme.sizeParamSampling; n++) {
					fprintf(fo, "tikz('%s_sampling_%d.tex', standAlone = TRUE, width=5, height=5)\nhist(c(as.matrix(A0)[,\"sampling_%d\"]", outCoda, n, n);
					for(i=1; i<niter; i++)
						fprintf(fo, ", as.matrix(A%d)[,\"sampling_%d\"]", i, n);
					fprintf(fo, "), 45,axes=T, xlab=\"sampling probability %d (%.3lf)\", main=NULL, yaxt='n', ylab = \"\")\ndev.off()\n\n", n+1, param.param[n].sampling);
				}
				fclose(fo);
			}
		}
	}
	return 0;
}

void nameLeavesRec(int n, int *ind, int length, TypeTree *tree) {
	if(tree->node[n].child < 0) {
		int i, tmp;
		tree->name[n] = (char*) malloc((length+1)*sizeof(char));
		tmp = *ind;
		tree->name[n][length] = '\0';
		for(i=0; i<length; i++) {
			tree->name[n][length-1-i] = '0'+tmp%10;
			tmp /= 10;
		}
		(*ind)++;
	}
}

void renameLeaves(TypeTree *tree) {
	int length, ind;
	if(tree->name == NULL) {
		int i;
		tree->name = (char**) malloc(tree->sizeBuf*sizeof(char*));
		for(i=0; i<tree->sizeBuf; i++)
			tree->name[i] = NULL;
	}
	length = (int) ceil(log(tree->size)/log(10.));
	ind = 0;
	nameLeavesRec(tree->root, &ind, length, tree);
}

TypeTree *getRandomTreeFossil(TypePiecewiseModelParam *param, int minSizeTree, int maxSizeTree, int minSliceB, double stopTime, void *rand_data) {
	TypeTree *tree, *treeTmp;
	int n, *count, minSlice;
	tree = NULL;
	count = (int*) malloc(param->size*sizeof(int));
	do {
		if(tree != NULL) {
			if(tree->info != NULL) {
				freeFossilFeature((TypeFossilFeature*) tree->info);
				tree->info = NULL;	
			}
			freeTree(tree);
			tree = NULL;
		}
		treeTmp = simulTreePiecewise(param, rand_data);
		if(treeTmp != NULL && treeTmp->size > minSizeTree) {
			tree = pruneFossilBis(treeTmp, (TypeFossilFeature*)treeTmp->info);
			freeFossilFeature((TypeFossilFeature*)treeTmp->info);
			treeTmp->info = NULL;
			freeTree(treeTmp);
			treeTmp = tree;
			tree = fixBinaryFossil(treeTmp, (TypeFossilFeature*) treeTmp->info);
			freeFossilFeature((TypeFossilFeature*)treeTmp->info);
			treeTmp->info = NULL;
			freeTree(treeTmp);
		} else {
			if(treeTmp != NULL) {
				freeFossilFeature((TypeFossilFeature*)treeTmp->info);
				treeTmp->info = NULL;
				freeTree(treeTmp);
			}
			tree = NULL;
		}
		if(tree != NULL) {
			int i;
			countNodeTime(count, tree, param->size, param->startTime);
			minSlice = count[0];
			for(i=1; i<param->size; i++)
				if(count[i] < minSlice)
					minSlice = count[i];
//			printf("tree size %d minSlice %d\n", tree->size, minSlice);
		} else {
//			printf("NULL\n");
			minSlice = 0;
		}
	} while(!(tree != NULL && tree->size >= minSizeTree  && minSlice >= minSliceB && tree->size <= maxSizeTree));
//	} while(!(tree != NULL && minSlice>minSizeTree && tree->size<maxSizeTree));
	free((void*)count);
	treeTmp = tree;
	if(stopTime != DBL_MAX) {
		tree = stopTreeFossil(treeTmp, (TypeFossilFeature*) treeTmp->info, stopTime);
		freeFossilFeature((TypeFossilFeature*)treeTmp->info);
		treeTmp->info = NULL;
		freeTree(treeTmp);
	}
	tree->minTimeInt.inf = tree->minTime;
	tree->minTimeInt.sup = tree->minTime;
	for(n=0; n <tree->size; n++)
		if(tree->node[n].child != NOSUCH)
			tree->time[n] = NO_TIME;
	return tree;
}

TypeTree *getRandomTreeFossilInt(TypePiecewiseModelParam *param, int minSizeTree, int maxSizeTree, int minSlice, double stopTime, double fos_width, void *rand_data) {
	TypeTree *tree;
	TypeFossilIntFeature *fint;
	TypeSampleIntData siData;
	tree = getRandomTreeFossil(param, minSizeTree, maxSizeTree, minSlice, stopTime, rand_data);
	siData.min = param->startTime[0];
	siData.max = (stopTime!=DBL_MAX)?stopTime:param->startTime[param->size];
	siData.length = fos_width;
	fint = sampleFossilIntFromFossil((TypeFossilFeature*) tree->info, tree->size, &siData, sampleFossilIntFixed, rand_data);
	freeFossilFeature((TypeFossilFeature*) tree->info);
	tree->info = fint;
	return tree;
}
