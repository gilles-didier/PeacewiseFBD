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
#define HELP_MESSAGE "\n--------------------------\n\nNAME\n\n	getAIC - returns the maximum likelihood and the AIC of a dataset under a model specification\n	\nSYNOPSIS\n\n	getAIC [OPTIONS] <tree(s)> <fossil ages> <model specification> [output File]\n\nDESCRIPTION\n\n	Return the maximum likelihood and the AIC of the dataset made of <tree(s)> and <fossil ages> under the model specification <model specification>.\n\n	Options are\n	-w <speciation rate width> <extinction rate width> <fossilization rate width> <sampling probability width>\n		set the widths of the sliding windows used for sampling the speciation, extinction and fossilization rates and of the sampling probability during the MCMC\n	-i <speciation rate upper bound> <extinction rate  upper bound> <fossilization rate upper bound> <sampling probability upper bound>\n		set the upper bounds of the interval in which the speciation, extinction and fossilization rates and the sampling probability are uniformly drawn at the beginning of the MCMC (lower bounds are all 0)\n	-f <proportion>\n		set the proportion of moves in the parameters in the MCMC proposal (other moves are drawn uniformly amont the fossil ages)\n	-a <speciation proportion> <extinction proportion> <fossilization proportion>\n		set the relation proportion of moves in the speciation, the extinction and the fossilzation rates (thus in the sampling probability)\n	-s <number>\n		set the number of samples required to estimate the maximum likelihood\n	-r <number>\n		set the random seed\n	-h\n		display help\n\nEXAMPLE\n\n./getAIC -s 10000 -f 0.2  -a 0.25 0.25 0.25 -w 0.05 0.05 0.05 0.25 -i 0.5 0.5 0.5 1.  ../data/Simulated_dataset_tree.newick ../data/Simulated_dataset_fossils.csv ../data/Simulated_Dataset_model_spec/Simul_model_spec_0.txt \n"


int main(int argc, char **argv) {	
	char *inputFileNameTree, *inputFileNameFossil, *inputFileNameScheme, outputName[STRING_SIZE], outputFileName[STRING_SIZE+20], *outputPrefix, option[256];
	FILE *fit, *fif, *fis, *fo;
	int i, nSamp = 5000;
	unsigned long int seed = 0;
	double al = 0.75, probSpe = 0.25, probExt = 0.25, probFos = 0.25, propParam = 0.2;
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
		if(option['r']) {
			option['r'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lu", &seed) == 1)
				i++;
			else
				error("a number is expected after -r");
		}
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
		inputFileNameFossil = argv[i++];
	} else {
		fprintf(stderr, "Please provide the name of a file containing a fossil list\n");
		exit(1);
	}
	if(i<argc) {
		inputFileNameScheme = argv[i++];
	} else {
		fprintf(stderr, "Please provide the name of a file containing a  model scheme\n");
		exit(1);
	}
	strcpy(outputName, inputFileNameScheme);
	if((outputPrefix = strrchr(outputName, '.')) != NULL)
		outputPrefix[0] = '\0';
	//if((outputPrefix=strrchr(outputName, '/')) == NULL)
		//outputPrefix = outputName;
	//else
		//outputPrefix++;
	
	if((fit = fopen(inputFileNameTree, "r")) && (fif = fopen(inputFileNameFossil, "r")) && (fis = fopen(inputFileNameScheme, "r"))) {
		TypeTree **tree;
		TypeFossilIntFeature *fint;
		TypeSamplingScheme scheme;
		TypePiecewiseModelParam maxParam;
		int sizeTree, nParam, n;
		double logLike;
		void *rand_data;
		my_rand_init();
		rand_data = my_rand_get_data();
		my_rand_set_seed(rand_data, seed);
		scheme = readSamplingScheme(fis);
        fclose(fis);
		tree = readTrees(fit);
        fclose(fit);
		sizeTree = 0;
		for(i=0; tree[i]!=NULL; i++) {
			int n;
			toBinary(tree[i]);
			if(tree[i]->name!=NULL)
				for(n=0; n<tree[i]->size; n++)
					if(tree[i]->name[n]!=NULL)
						fixSpace(tree[i]->name[n]);
			for(n=0; n<tree[i]->size; n++)
				if(tree[i]->time[n] < scheme.param.startTime[scheme.param.size])
					tree[i]->time[n] = NO_TIME;
				else
					tree[i]->time[n] = scheme.param.startTime[scheme.param.size];
			tree[i]->minTime = scheme.param.startTime[0];
			tree[i]->maxTime = scheme.param.startTime[scheme.param.size];
			sizeTree++;
		}
		if(sizeTree == 0)
			error("No tree!\n");
		reindexTreesFromName(tree, sizeTree);
		fint = getFossilIntFeature(fif, tree[0]->name, tree[0]->size);
        fclose(fif);
		if(getMaxFossilIntTime(fint) > 0.)
			negateFossilInt(fint);
		fixStatus(tree[0], fint);
		for(n=0; n<tree[0]->size; n++)
			if(fint->status[n] == unknownNodeStatus)
				for(i=0; tree[i]!=NULL; i++)
					tree[i]->time[n] = fint->endTimeTable[n].inf;
					
		nParam = scheme.sizeParamBirth+scheme.sizeParamDeath+scheme.sizeParamFossil+scheme.sizeParamSampling;
		maxParam.param = (TypeModelParam*) malloc(scheme.param.size*sizeof(TypeModelParam));
		maxParam.startTime = scheme.param.startTime;
		logLike = MCMCFillMaxParamSampleSC(tree, sizeTree, fint, scheme, al, nSamp, propParam, &windSize, &init, probSpe,probExt, probFos, &maxParam, rand_data);
		printf("\n\nlog-Likelihood\t%le\nAIC\t%le\n", logLike, 2.*(((double)nParam+fint->sizeFossil)-logLike));
		sprintf(outputFileName, "%s_ML_model.txt", outputName);
		printf("output file: %s\n", outputFileName);
		if((fo = fopen(outputFileName, "w"))) {
			fprintf(fo, "#%d iterations\n#log-likelihood max\t%le\n#%d parameters\n#AIC\t%le (%d param and %d fossils)\n", nSamp, logLike, nParam, 2.*(((double)nParam+fint->sizeFossil)-logLike), nParam, fint->sizeFossil);
			printPiecewiseModel(fo, &maxParam);
			fclose(fo);
		} else
			error("Can't open %s\n", outputFileName);
		freeSamplingScheme(&scheme);
		free((void*)maxParam.param);
	}
	return 0;
}
