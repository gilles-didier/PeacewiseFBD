#ifndef MCMCImportanceSamplingFossilIntF
#define MCMCImportanceSamplingFossilIntF

#include <stdlib.h>
#include <stdio.h>

#include "Tree.h"
#include "Fossil.h"
#include "FossilInt.h"
#include "PiecewiseModel.h"

typedef struct MCMC_PARAM {
	double al;
	int burn, gap, iter;
} TypeMCMCParam;



#ifdef __cplusplus
extern "C" {
#endif

double MCMCFillMaxParamSampleSC(TypeTree **tree, int nTree, TypeFossilIntFeature *fi, TypeSamplingScheme sc, double al, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt, double probFos, TypePiecewiseModelParam *maxParam, void *rand_data);
double MCMCGetMaxLogLikelihoodSC(TypeTree **tree, int nTree, TypeFossilIntFeature *fi, TypeSamplingScheme sampScheme, double al, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt, double probFos, void *rand_data);
void MCMCSamplingSamplePosteriorParametersOnly(FILE *fout, FILE *find, TypeTree **tree, int nTree, TypeFossilIntFeature *fi, TypeSamplingScheme sampScheme, double al, int burn, int gap, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt, double probFos, void *rand_data);
void MCMCSamplingSamplePosteriorParametersFossils(FILE *fout, FILE *find, TypeTree **tree, int nTree, TypeFossilIntFeature *fi, TypeSamplingScheme sampScheme, double al, int burn, int gap, int iter, double prop, TypeModelParam *windSize, TypeModelParam *init, double probSpe, double probExt, double probFos, void *rand_data);

#ifdef __cplusplus
}
#endif



#endif
