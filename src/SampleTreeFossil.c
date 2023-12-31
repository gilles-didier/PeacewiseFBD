#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>


#include "Utils.h"
#include "MyR.h"
#include "Tree.h"
#include "SimulPiecewise.h"
#include "Fossil.h"
#include "FossilInt.h"
#include "DrawTreeCairo.h"
#include "DrawTreePSTricks.h"
#include "DrawTreeTikz.h"
#include "DrawTreeGeneric.h"
#include "DrawFossil.h"
#include "DrawFossilInt.h"
#include "MyRandom.h"


//./samp -p 0. 2 1.5 3 0.5 -p 5 1 0.9 5 1. 
//./samp -p -200. 0.2 0.15 0.5 0.5 -p -100 0.1 0.09 0.5 1. 
//./samp -p -20. 0.5 0.12 0.5 1. -p -10 0.1 0.2 0.1 1. -t 0.

//./samp  -p -10 0.4 0.2 0.5 1. -t 0.

//./samp -p -320 0.2 0.18 0.5 0.2 -p -290 0.4 0.2 0.5 1. -p -270 0.2 0.18 0.5 1. -t 0.

//./samp -p -350. 0.03 0.01 0.02 0.2 -p -250 0.1 0.01 0.02 1. -p -225 0.01 0.01 0.02 1. -t 0.  -i 20
//./samp -p -350. 0.03 0.01 0.02 0.2 -p -250 0.1 0.01 0.02 1. -p -225 0.01 0.01 0.02 1. -t 0. -x 6 -y 15 -i 20

#define MAX_PARAM 10

#define STRING_SIZE 300
#define HELP_MESSAGE "\nusage: samp [options] [<output file>]\n\nsample simulates random trees and fossils finds and saves them in Newick format\n\nOptions:\n\t-h : display help\n\t-b <birth>\t: set birth rate\n\t-d <death>\t: set death rate\n\t-f <fossil>\t: set fossil find rate\n\t-m <min>\t: set minimum number of contemporary species of a simulation to be considered\n\t-M <size>\t: set maximum size of a simulation to be considered\n\t-i <niter>\t: set the number of simulations\n\t-t <time> : the end time of the diversification (start is always 0)\n"

//./samp -p 0 2 1 1 1 -p 5 3 2 5 1 -t 6 -i 4 simul.newick

void nameLeavesRec(int n, int *ind, int length, TypeTree *tree);
void renameLeaves(TypeTree *tree);
TypeTree *getRandomTreeFossil(TypePiecewiseModelParam *param, int minSizeTree, int maxSizeTree, int minSlice, double stopTime, void *rand_data);

int main(int argc, char **argv) {	
	char outputName[STRING_SIZE], outputFileName[STRING_SIZE], outputFileNameG[STRING_SIZE],*outputPrefix, option[256], format = '1';
	FILE *fi, *fo, *fs;
	int i, nb, niter = 5, number, minContemp = 10, minFossil = 5, maxSizeTree = 156, minSizeTree = 153, minSlice = 25;
	double birth = 2., death = 1., fossil = 1., maxTime = 5., figwidth = 500., fos_width = 5.;
	TypeAdditionalDrawTreeGeneric add;
	TypeDataDrawFossil data;
	TypeDataDrawFossilInt dataInt;
	TypePiecewiseModelParam param;
	TypeModelParam par[MAX_PARAM];
	double start[MAX_PARAM+1];
	void *rand_data;
	param.size = 0;
	param.param = par;
	param.startTime = start;

	for(i=0; i<256; i++)
		option[i] = 0;
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		int j;
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['p']) {
			option['p'] = 0;
			if(param.size < MAX_PARAM) {
				if((i+1)<argc && sscanf(argv[i+1], "%lf", &(param.startTime[param.size])) == 1)
					i++;
				if((i+1)<argc && sscanf(argv[i+1], "%lf", &(param.param[param.size].birth)) == 1)
					i++;
				if((i+1)<argc && sscanf(argv[i+1], "%lf", &(param.param[param.size].death)) == 1)
					i++;
				if((i+1)<argc && sscanf(argv[i+1], "%lf", &(param.param[param.size].fossil)) == 1)
					i++;
				if((i+1)<argc && sscanf(argv[i+1], "%lf", &(param.param[param.size].sampling)) == 1)
					i++;
				param.size++;
			}
		}
		if(option['d']) {
			option['d'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &death) == 1)
				i++;
		}
		if(option['f']) {
			option['f'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &fossil) == 1)
				i++;
		}
		if(option['i']) {
			option['i'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &niter) == 1)
				i++;
		}
		if(option['m']) {
			option['m'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &minContemp) == 1)
				i++;
		}
		if(option['M']) {
			option['M'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &maxSizeTree) == 1)
				i++;
		}
		if(option['t']) {
			option['t'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &maxTime) == 1)
				i++;
		}
		if(option['x']) {
			option['x'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%c", &format) == 1)
				i++;
			else
				error("a character is required after option -f");
		}
		if(option['y']) {
			option['y'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &figwidth) == 1)
				i++;
			else
				error(ErrorArgument, "a number is required after option -t");
		}
		if(option['h']) {
			printf("%s\n", HELP_MESSAGE);
			exit(EXIT_SUCCESS);
		}
	}
	param.startTime[param.size] = maxTime;
	if(!(i<argc && sscanf(argv[i++], "%s", outputName) == 1)) {
		strcpy(outputName, "Tree_simul");
	}
	if((outputPrefix = strrchr(outputName, '.')) != NULL)
		outputPrefix[0] = '\0';
	if((outputPrefix=strrchr(outputName, '/')) == NULL)
		outputPrefix = outputName;
	else
		outputPrefix++;
printPiecewiseModel(stdout, &param);
//exit(0);
		rand_data = my_rand_get_data();
	for(nb=1; nb <= niter; nb++) {
		TypeTree *tree, *tree1, *tree2;
		FILE *ft;
		char outputName[STRING_SIZE];
		TypeInfoDrawTreeGeneric info;
		tree = NULL;

		//do {
			//if(tree != NULL)
				//freeTree(tree);
			//tree = simulTreePiecewise(&param, rand_data);
////if(tree != NULL)
	////printf("tree size %d contemp %d\n", tree->size, countContemp(tree));
		//} while(!(tree != NULL && (countContemp(tree))>=minContemp && tree->size>=minSizeTree && tree->size<=maxSizeTree));
////		} while(!(tree != NULL && tree->size>=minSizeTree && tree->size<=maxSizeTree));
////		} while(!(tree != NULL && (countContemp(tree))==0 && tree->size>=minSizeTree && tree->size<=maxSizeTree));
tree = getRandomTreeFossil(&param, minSizeTree, maxSizeTree, minSlice, DBL_MAX, rand_data);
printf("tree %d size %d contemp %d\n", nb, tree->size, countContemp(tree));
		reorderTreeSize(tree);
		renameLeaves(tree);
//fprintSubtreeNewick(stdout, tree->root, tree);
//printTreeDebug(stdout, tree->root, tree, tree->name);
//exit(0);
		add.draw = drawFossil;
		data.color = (TypeRGB) {.red = 0.3, .green = 0., .blue = 0.};
		data.radius = 3.;
		data.alpha = 0.75;
		info.param.tmin = tree->minTime;
//		info.param.tmax = tree->maxTime;
		info.param.tmax = -268.8;
//		info.param.scaleStep = 25.;
		info.param.scaleStep = 5.;
		info.param.width = figwidth;
		data.fos = (TypeFossilFeature*) tree->info;
		add.data = (void*) &data;
		sprintf(outputFileName, "%s_sim_all_%d", outputPrefix, nb);
		switch(format) {
			case '1':
				sprintf(outputFileNameG, "%s_tree.pdf", outputFileName);
				setFunctPDF(&(info.funct));
				data.drawDot = drawDotCairo;
				drawTreeFileGeneric(outputFileNameG, tree, &info, &add, NULL);
				break;
			case '2':
				sprintf(outputFileNameG, "%s_tree.ps", outputFileName);
				setFunctPS(&(info.funct));
				data.drawDot = drawDotCairo;
				drawTreeFileGeneric(outputFileNameG, tree, &info, &add, NULL);
				break;
			case '3':
				sprintf(outputFileNameG, "%s_tree.png", outputFileName);
				setFunctPNG(&(info.funct));
				data.drawDot = drawDotCairo;
				drawTreeFileGeneric(outputFileNameG, tree, &info, &add, NULL);
				break;
			case '4':
				sprintf(outputFileNameG, "%s_tree.svg", outputFileName);
				setFunctSVG(&(info.funct));
				data.drawDot = drawDotCairo;
				drawTreeFileGeneric(outputFileNameG, tree, &info, &add, NULL);
				break;
			case '5':
				sprintf(outputFileNameG, "%s_tree_pst.tex", outputFileName);
				setFunctPSTricks(&(info.funct));
				data.drawDot = drawDotPSTricks;
				drawTreeFileGeneric(outputFileNameG, tree, &info, &add, NULL);
				break;
			case '6':
				sprintf(outputFileNameG, "%s_tree_tikz.tex", outputFileName);
				setFunctTikz(&(info.funct));
				data.drawDot = drawDotTikz;
				data.radius = 0.08;
				drawTreeFileGeneric(outputFileNameG, tree, &info, &add, NULL);
				break;
			default:
				;
		}
		tree->minTimeInt.inf = NO_TIME;
		tree->minTimeInt.sup = NO_TIME;
		tree->maxTimeInt.inf = NO_TIME;
		tree->maxTimeInt.sup = NO_TIME;
		tree->name = NULL;
		tree1 = pruneFossilBis(tree, (TypeFossilFeature*) tree->info);
//		freeFossilFeature((TypeFossilFeature*) tree->info);
		tree2 = fixBinaryFossil(tree1, (TypeFossilFeature*) tree1->info);
		freeFossilFeature((TypeFossilFeature*) tree1->info);
		freeTree(tree1);
		sprintf(outputFileName, "%s_sim_obs_%d", outputPrefix, nb);
		data.fos = (TypeFossilFeature*) tree2->info;
		switch(format) {
			case '1':
				sprintf(outputFileNameG, "%s_tree.pdf", outputFileName);
				setFunctPDF(&(info.funct));
				data.drawDot = drawDotCairo;
				drawTreeFileGeneric(outputFileNameG, tree2, &info, &add, NULL);
				break;
			case '2':
				sprintf(outputFileNameG, "%s_tree.ps", outputFileName);
				setFunctPS(&(info.funct));
				data.drawDot = drawDotCairo;
				drawTreeFileGeneric(outputFileNameG, tree2, &info, &add, NULL);
				break;
			case '3':
				sprintf(outputFileNameG, "%s_tree.png", outputFileName);
				setFunctPNG(&(info.funct));
				data.drawDot = drawDotCairo;
				drawTreeFileGeneric(outputFileNameG, tree2, &info, &add, NULL);
				break;
			case '4':
				sprintf(outputFileNameG, "%s_tree.svg", outputFileName);
				setFunctSVG(&(info.funct));
				data.drawDot = drawDotCairo;
				drawTreeFileGeneric(outputFileNameG, tree2, &info, &add, NULL);
				break;
			case '5':
				sprintf(outputFileNameG, "%s_tree_pst.tex", outputFileName);
				setFunctPSTricks(&(info.funct));
				data.drawDot = drawDotPSTricks;
				drawTreeFileGeneric(outputFileNameG, tree2, &info, &add, NULL);
				break;
			case '6':
				sprintf(outputFileNameG, "%s_tree_tikz.tex", outputFileName);
				setFunctTikz(&(info.funct));
				data.drawDot = drawDotTikz;
				data.radius = 0.08;
				drawTreeFileGeneric(outputFileNameG, tree2, &info, &add, NULL);
				break;
			default:
				;
		}
		TypeFossilIntFeature *fint;
		TypeSampleIntData siData;
		siData.min = param.startTime[0];
		siData.max = param.startTime[param.size];
		siData.length = fos_width;
		fint = sampleFossilIntFromFossil((TypeFossilFeature*) tree2->info, tree2->size, &siData, sampleFossilIntFixed, rand_data);
		freeFossilFeature((TypeFossilFeature*) tree2->info);
		tree2->info = fint;
		
		add.draw = drawFossilInt;
		dataInt.fos = (TypeFossilFeature*) tree2->info;
		dataInt.color = (TypeRGB) {.red = 0.3, .green = 0., .blue = 0.};
		dataInt.radius = 3.;
		dataInt.alpha = 0.75;
		add.data = (void*) &dataInt;
		switch(format) {
			case '1':
				sprintf(outputFileNameG, "%s_tree_int.pdf", outputFileName);
				setFunctPDF(&(info.funct));
				dataInt.drawLineDot = drawLineDotCairo;
				drawTreeFileGeneric(outputFileNameG, tree2, &info, &add, NULL);
				break;
			case '2':
				sprintf(outputFileNameG, "%s_tree_int.ps", outputFileName);
				setFunctPS(&(info.funct));
				dataInt.drawLineDot = drawLineDotCairo;
				drawTreeFileGeneric(outputFileNameG, tree2, &info, &add, NULL);
				break;
			case '3':
				sprintf(outputFileNameG, "%s_tree_int.png", outputFileName);
				setFunctPNG(&(info.funct));
				dataInt.drawLineDot = drawLineDotCairo;
				drawTreeFileGeneric(outputFileNameG, tree2, &info, &add, NULL);
				break;
			case '4':
				sprintf(outputFileNameG, "%s_tree_int.svg", outputFileName);
				setFunctSVG(&(info.funct));
				dataInt.drawLineDot = drawLineDotCairo;
				drawTreeFileGeneric(outputFileNameG, tree2, &info, &add, NULL);
				break;
			case '5':
				sprintf(outputFileNameG, "%s_tree_int_pst.tex", outputFileName);
				setFunctPSTricks(&(info.funct));
				dataInt.drawLineDot = drawLineDotPSTricks;
				drawTreeFileGeneric(outputFileNameG, tree2, &info, &add, NULL);
				break;
			case '6':
				sprintf(outputFileNameG, "%s_tree_int_tikz.tex", outputFileName);
				setFunctTikz(&(info.funct));
				dataInt.drawLineDot = drawLineDotTikz;
				drawTreeFileGeneric(outputFileNameG, tree2, &info, &add, NULL);
				break;
			default:
				;
		}
		sprintf(outputName, "%s_sim_%d.newick", outputPrefix, nb);
		if(ft = fopen(outputName, "w")) {
			int n;
			tree->name = nameLeaves("leaf_", tree);
			tree->comment = (char**) malloc(tree->sizeBuf*sizeof(char*));
			for(n=0; n<tree->sizeBuf; n++)
				tree->comment[n] = NULL;
//			fillCommentFossil(tree, (TypeFossilFeature*) tree->info);
			fprintSubtreeNewick(ft, tree->root, tree);
			fprintf(ft, "\n\n");
			tree2->name = nameLeaves("leaf_", tree2);
			tree2->comment = (char**) malloc(tree2->sizeBuf*sizeof(char*));
			for(n=0; n<tree2->sizeBuf; n++)
				tree2->comment[n] = NULL;
//			fillCommentFossil(tree2, (TypeFossilFeature*) tree2->info);
			fprintSubtreeNewick(ft, tree2->root, tree2);
			fprintf(ft, "\n");
			fclose(ft);
		}
		sprintf(outputName, "%s_sim_tree_%d.newick", outputPrefix, nb);
		if(ft = fopen(outputName, "w")) {
			free((void*)tree2->name);
			tree2->name = nameAll("Tip_", "Int_", tree2);
			free((void*)tree2->comment);
			tree2->comment = NULL;
			fprintSubtreeNewick(ft, tree2->root, tree2);
			fprintf(ft, "\n");
			fclose(ft);
		}
		sprintf(outputName, "%s_sim_fos_%d.csv", outputPrefix, nb);
		if(ft = fopen(outputName, "w")) {
			fprintFossilIntFeature(ft, (TypeFossilIntFeature*) tree2->info, tree2->name, tree2->size);
			fclose(ft);
		}			
		freeFossilIntFeature((TypeFossilIntFeature*) tree2->info);
		freeTree(tree2);
		freeTree(tree);
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


TypeTree *getRandomTreeFossil(TypePiecewiseModelParam *param, int minSizeTree, int maxSizeTree, int minSliceX, double stopTime, void *rand_data) {
	TypeTree *tree, *treeTmp;
	int n, *count, minSlice, imin=0;
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
		tree = simulTreePiecewise(param, rand_data);
//		treeTmp = simulTreePiecewise(param, rand_data);
		//if(treeTmp != NULL && treeTmp->size > minSizeTree) {
			//tree = pruneFossilBis(treeTmp, (TypeFossilFeature*)treeTmp->info);
			//freeFossilFeature((TypeFossilFeature*)treeTmp->info);
			//treeTmp->info = NULL;
			//freeTree(treeTmp);
			//treeTmp = tree;
			//tree = fixBinaryFossil(treeTmp, (TypeFossilFeature*) treeTmp->info);
			//freeFossilFeature((TypeFossilFeature*)treeTmp->info);
			//treeTmp->info = NULL;
			//freeTree(treeTmp);
		//} else {
			//if(treeTmp != NULL) {
				//freeFossilFeature((TypeFossilFeature*)treeTmp->info);
				//treeTmp->info = NULL;
				//freeTree(treeTmp);
			//}
			//tree = NULL;
		//}
		if(tree != NULL) {
			int i;
			countNodeTime(count, tree, param->size, param->startTime);
			minSlice = count[0];
			imin = 0;
			for(i=1; i<param->size; i++)
				if(count[i] < minSlice) {
					minSlice = count[i];
					imin = i;
				}
		} else
			minSlice = 0;
//if(tree != NULL)
	//printf("tree  size %d min %d contemp %d\n", tree->size, minSlice, countContemp(tree));
//	} while(!(tree != NULL && tree->size>minSizeTree && tree->size<maxSizeTree));
	} while(!(tree != NULL && tree->size>=minSizeTree && minSlice>=minSliceX && tree->size<maxSizeTree));
	free((void*)count);
	if(stopTime != DBL_MAX) {
		treeTmp = tree;
		tree = stopTreeFossil(treeTmp, (TypeFossilFeature*) treeTmp->info, stopTime);
		freeFossilFeature((TypeFossilFeature*)treeTmp->info);
		treeTmp->info = NULL;
		freeTree(treeTmp);
	}
	tree->minTimeInt.inf = tree->minTime;
	tree->minTimeInt.sup = tree->minTime;
	//for(n=0; n <tree->size; n++)
		//if(tree->node[n].child != NOSUCH)
			//tree->time[n] = NO_TIME;
	return tree;
}
