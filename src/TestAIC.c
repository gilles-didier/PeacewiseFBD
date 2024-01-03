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

#define MAX_PARAM 20

#define STRING_SIZE 300
#define HELP_MESSAGE "\n--------------------------\n\nNAME\n\n	testAIC - simulates trees under a given model and writes files containing the Akaike weights of several model specifications over these simulated trees\n	\nSYNOPSIS\n\n	testAIC [OPTIONS] <model> <model specification 1>  <model specification 2> ... <model specification n>\n\nDESCRIPTION\n\n	 Simulate trees under the model <model> and write the files containing the Akaike weights of model specifications <model specification 1>  <model specification 2> ... <model specification n> over these simulated trees\n\n	Options are\n	-w <speciation rate width> <extinction rate width> <fossilization rate width> <sampling probability width>\n		set the widths of the sliding windows used for sampling the speciation, extinction and fossilization rates and of the sampling probability during the MCMC\n	-i <speciation rate upper bound> <extinction rate  upper bound> <fossilization rate upper bound> <sampling probability upper bound>\n		set the upper bounds of the interval in which the speciation, extinction and fossilization rates and the sampling probability are uniformly drawn at the beginning of the MCMC (lower bounds are all 0)\n	-f <proportion>\n		set the proportion of moves in the parameters in the MCMC proposal (other moves are drawn uniformly amont the fossil ages)\n	-a <speciation proportion> <extinction proportion> <fossilization proportion>\n		set the relation proportion of moves in the speciation, the extinction and the fossilzation rates (thus in the sampling probability)\n	-s <number>\n		set the number of samples required to estimate the maximum likelihood\n	-r <number>\n		set the random seed\n	-l <length>\n		set the length of the stratigraphic range of the simulated fossils to <length>\n	-m <number>\n		set the minimum size of the simulated trees to <number>\n	-M <number>\n		set the maximum size of the simulated trees to <number>\n	-x <number>\n		set the minimum number of events in a slice to <number>\n	-u <number>\n		set the number of trees to simulate to <number>\n	-h\n		display help\n\nEXAMPLE\n\n./testAIC -l 5 -s 25000 -f 0.2 -a 0.25 0.25 0.25 -w 0.05 0.05 0.05 0.25 -i 0.5 0.5 0.5 1. -u 30 -m 100 -M 500 ../data/Simulated_Dataset_model_spec/Simul_model.txt ../data/Simulated_Dataset_model_spec/Simul_model_spec_0.txt ../data/Simulated_Dataset_model_spec/Simul_model_spec_1.txt ../data/Simulated_Dataset_model_spec/Simul_model_spec_1E.txt\n"


void nameLeavesRec(int n, int *ind, int length, TypeTree *tree);
void renameLeaves(TypeTree *tree);
TypeTree *getRandomTreeFossil(TypePiecewiseModelParam *param, int minSizeTree, int maxSizeTree, int minSlice, double stopTime, void *rand_data);
TypeTree *getRandomTreeFossilInt(TypePiecewiseModelParam *param, int minSizeTree, int maxSizeTree, int minSlice, double stopTime, double fos_width, void *rand_data);


int main(int argc, char **argv) {	
	char *inputFileNameModel0, *(inputFileNameSchemeAlt[MAX_PARAM]), outputFileName[STRING_SIZE], tmpName[STRING_SIZE], *modelName, option[256];
	FILE *fim0, *fis, *fo;
	int sizeAlt = 0, i, niter = 5, maxSizeTree = 150, minSizeTree = 50, minSlice= 5, nSamp = 5000;
	unsigned long int seed = 1;
	double al = 0.75, probSpe = 0.25, probExt = 0.25, probFos = 0.25, propParam = 0.2, fos_width = 1., endTime = 0., stopTime = DBL_MAX;
	TypeModelParam windSize = {.birth=0.05, .death = 0.05, .fossil = 0.05, .sampling = 0.5}, init = {.birth=0.5, .death = 0.5, .fossil = 0.1, .sampling = 1.};

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
		if(option['f']) {
			option['f'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &propParam) == 1)
				i++;
			else
				error("a number is required after option -t");
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
		if(option['x']) {
			option['x'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &minSlice) == 1)
				i++;
		}
		if(option['l']) {
			option['l'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &fos_width) == 1)
				i++;
		}
		if(option['E']) {
			option['E'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &endTime) == 1)
				i++;
		}
		if(option['h']) {
			printf("%s\n", HELP_MESSAGE);
			return 0;
		}
	}
	if(i<argc) {
		inputFileNameModel0 = argv[i++];
	} else
		error("Please provide the name of a file containing a piecewise model\n");
	while(i<argc && sizeAlt<MAX_PARAM) {
		inputFileNameSchemeAlt[sizeAlt++] = argv[i++];
	}
	if(sizeAlt == 0)
		error("Please provide the name of a file containing a model scheme\n");
	strcpy(tmpName, inputFileNameModel0);
	if((modelName = strrchr(tmpName, '.')) != NULL)
		modelName[0] = '\0';
	if((modelName=strrchr(tmpName, '/')) == NULL)
		modelName = tmpName;
	else
		modelName++;
	if((fim0 = fopen(inputFileNameModel0, "r"))) {
		int a, nParamAlt[MAX_PARAM];
		TypePiecewiseModelParam param0;
		TypeSamplingScheme schemeAlt[MAX_PARAM];
		param0 = readPiecewiseModelParam(fim0);
		fclose(fim0);
		printPiecewiseModel(stdout, &param0);
		for(a=0; a<sizeAlt; a++) {
			char *tmp;
			if((fis = fopen(inputFileNameSchemeAlt[a], "r")))
				schemeAlt[a] = readSamplingScheme(fis);
			else
				error("Error reading %s\n", inputFileNameSchemeAlt[a]);
			fclose(fis);
			printSamplingScheme(stdout, schemeAlt[a]);
			nParamAlt[a] = schemeAlt[a].sizeParamBirth+schemeAlt[a].sizeParamDeath+schemeAlt[a].sizeParamFossil+schemeAlt[a].sizeParamSampling;
			if((tmp = strrchr(inputFileNameSchemeAlt[a], '.')) != NULL)
				tmp[0] = '\0';
			sprintf(outputFileName, "%s_weight_AIC_%s.csv", inputFileNameSchemeAlt[a], modelName);
			if((fo = fopen(outputFileName, "w")))
				fclose(fo);
		}
			void *rand_data_tree;

			my_rand_init();
			rand_data_tree = my_rand_get_data();
			my_rand_set_seed(rand_data_tree, seed);
		#pragma omp parallel
		{
			int i;
			#pragma omp for
			for (i=0; i<niter; i++)	{
				int a;
				double logLikeAlt[MAX_PARAM], den, min;
				void *rand_data;
				TypeTree *tree;
				rand_data = my_rand_get_data();
				my_rand_set_seed(rand_data, seed);
//				tree = getRandomTreeFossilInt(&param0, minSizeTree, maxSizeTree, minSlice, stopTime, fos_width, rand_data);
				#pragma omp critical
				{
					tree = getRandomTreeFossilInt(&param0, minSizeTree, maxSizeTree, minSlice, stopTime, fos_width, rand_data_tree);
				}
				min = DBL_MAX;
				for(a=0; a<sizeAlt; a++) {
					double aic;
					logLikeAlt[a] = MCMCGetMaxLogLikelihoodSC(&tree, 1, (TypeFossilIntFeature*) tree->info, schemeAlt[a], al, nSamp, propParam, &windSize, &init, probSpe,probExt, probFos, rand_data);
					aic = 2.*(((double)nParamAlt[a])-logLikeAlt[a]);
					if(aic<min)
						min = aic;
				}
				my_rand_free_data(rand_data);
				den = 0.;
				for(a=0; a<sizeAlt; a++)
					den += exp(-0.5*(2.*(((double)nParamAlt[a])-logLikeAlt[a])-min));
				#pragma omp critical
				{
					for(a=0; a<sizeAlt; a++) {
						sprintf(outputFileName, "%s_weight_AIC_%s.csv", inputFileNameSchemeAlt[a], modelName);
						if((fo = fopen(outputFileName, "a"))) {
							fprintf(fo, "%le\n", exp(-0.5*(2.*(((double)nParamAlt[a])-logLikeAlt[a])-min))/den);
							fclose(fo);
						} else
							error("Cannot open %s\n", outputFileName);
					}
				}
			}
		}
	} else
		error("Cannot open %s or %s\n", inputFileNameModel0);
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

#define MAX_TRY 10000

TypeTree *getRandomTreeFossil(TypePiecewiseModelParam *param, int minSizeTree, int maxSizeTree, int minSliceB, double stopTime, void *rand_data) {
	TypeTree *tree, *treeTmp;
	int n, *count, minSlice, try = 0;
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
		} else
			minSlice = 0;
		try++;
//	} while(!(tree != NULL && tree->size >= minSizeTree  && minSlice >= minSliceB && tree->size <= maxSizeTree));
	} while(!(tree != NULL && tree->size>=minSizeTree && tree->size<=maxSizeTree) && try < MAX_TRY);
	free((void*)count);
	if(try >= MAX_TRY)
		error("Problem in simulation: too much reject (check model and constraints)\n");
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
