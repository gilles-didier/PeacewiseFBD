#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <nlopt.h>

#include "Utils.h"

#include "MinimizeNLOpt.h"


//#define NLOPT_ALGO NLOPT_GN_ISRES
//#define NLOPT_ALGO NLOPT_GN_ESCH
//#define NLOPT_ALGO NLOPT_LN_BOBYQA
//#define NLOPT_ALGO NLOPT_LN_COBYLA
//#define NLOPT_ALGO NLOPT_AUGLAG
//#define NLOPT_ALGO NLOPT_GN_DIRECT_L
#define NLOPT_ALGO NLOPT_LN_SBPLX

#define MINVAL 0.01
#define INFTY 1E99
#define RINFTY 1E99
#define DEF 10
#define MIN_VAL 0.000001
#define TOLERANCE_CONSTRAINT 0.000000001
#define TOLERANCE_OPTIM 0.001


typedef struct MINIMIZATION_SET_TREE_FOS_DATA {
	TypeTree **tree;
	TypeFossilFeature *fos;
	int nTree;
	TypePiecewiseModelParam param;
	TypeSamplingCurrent sampCurrent;
	TypeSamplingScheme sampScheme;
	TypeLikelihoodSetTreeFosFunction *likelihood;
} TypeMinimizationSetTreeFossilData;



#define TAG_SPE "SPE"
#define TAG_EXT "EXT"
#define TAG_FOS "FOS"
#define TAG_TRI "TRI"
#define TAG_TOL "TOL"
#define TAG_ITE "ITE"
#define SIZE_TAG 20
#define SIZE_VAL 100

void fprintNLoptOptionTag(FILE *f, TypeNLOptOption *option) {
    fprintf(f, ":%s [%lE;%lE]\n", TAG_SPE, option->infSpe, option->supSpe);
    fprintf(f, ":%s [%lE;%lE]\n", TAG_EXT, option->infExt, option->supExt);
    fprintf(f, ":%s [%lE;%lE]\n", TAG_FOS, option->infFos, option->supFos);
    fprintf(f, ":%s %d\n", TAG_TRI, option->trials);
    fprintf(f, ":%s %lE\n", TAG_TOL, option->tolOptim);
    fprintf(f, ":%s %d\n", TAG_ITE, option->maxIter);
}

void fscanNLoptOptionTag(FILE *f, TypeNLOptOption *option) {
    char c, tag[SIZE_TAG+1], val[SIZE_VAL+1];
    for(c=fgetc(f); c!=EOF && isspace(c); c=fgetc(f));
    while(c == ':') {
        int i;
        c=fgetc(f);
        for(i=0; c!=EOF && !isspace(c) && i<SIZE_TAG; c=fgetc(f))
            tag[i++] = c;
        tag[i] = '\0';
        if(i>=SIZE_TAG) {
            fprintf(stderr, "Error when reading an optimizer options file - Tag too long:\n%s...\n", tag);
            exit(1);
        }
        for(; c!=EOF && isspace(c); c=fgetc(f));
        for(i=0; c!=EOF && !isspace(c) && i<SIZE_VAL; c=fgetc(f))
            val[i++] = c;
        val[i] = '\0';
        if(i>=SIZE_VAL) {
            fprintf(stderr, "Error when reading an optimizer options file - value too long:\n%s...\n", val);
            exit(1);
        }
        if(strcmp(tag, TAG_SPE) == 0)
            toInterval(val, &(option->infSpe), &(option->supSpe));
        if(strcmp(tag, TAG_EXT) == 0)
            toInterval(val, &(option->infExt), &(option->supExt));
        if(strcmp(tag, TAG_FOS) == 0)
            toInterval(val, &(option->infFos), &(option->supFos));
        if(strcmp(tag, TAG_TRI) == 0)
            option->trials = atoi(val);
        if(strcmp(tag, TAG_TOL) == 0)
            option->tolOptim = atof(val);
        if(strcmp(tag, TAG_ITE) == 0)
            option->maxIter = atoi(val);
        for(; c!=EOF && isspace(c); c=fgetc(f));
    }
}

void fprintNLoptOption(FILE *f, TypeNLOptOption *option) {
    fprintf(f, "Speciation rates are sampled in [%.2lE:%.2lE]\n", option->infSpe, option->supSpe);
    fprintf(f, "Extinction rates are sampled in [%.2lE:%.2lE]\n", option->infExt, option->supExt);
    fprintf(f, "Fossil rates are sampled in [%.2lE:%.2lE]\n", option->infFos, option->supFos);
    fprintf(f, "Optimizer runs %d trials and stops with tolerance %.lE or after more than %d iterations.\n", option->trials, option->tolOptim, option->maxIter);
}

void sprintNLoptOption(char *buffer, TypeNLOptOption *option) {
    buffer += sprintf(buffer, "Speciation rates are sampled in [%.2lE:%.2lE]\n", option->infSpe, option->supSpe);
    buffer += sprintf(buffer, "Extinction rates are sampled in [%.2lE:%.2lE]\n", option->infExt, option->supExt);
    buffer += sprintf(buffer, "Fossil rates are sampled in [%.2lE:%.2lE]\n", option->infFos, option->supFos);
    buffer += sprintf(buffer, "Optimizer runs %d trials and stops with tolerance %.lE or after more than %d iterations.\n", option->trials, option->tolOptim, option->maxIter);
}



void setSamplingCurrentFromVector(TypeSamplingCurrent *sampCurrent, double *x, TypeSamplingScheme sampScheme) {
	int i, ind = 0;
	for(i=0; i<sampScheme.sizeParamBirth; i++)
		sampCurrent->currentBirth[i] = x[ind++];
	for(i=0; i<sampScheme.sizeParamDeath; i++)
		sampCurrent->currentDeath[i] = x[ind++];
	for(i=0; i<sampScheme.sizeParamFossil; i++)
		sampCurrent->currentFossil[i] = x[ind++];
	for(i=0; i<sampScheme.sizeParamSampling; i++)
		sampCurrent->currentSampling[i] = x[ind++];
}

void setVectorFromSamplingCurrent(double *x, TypeSamplingCurrent sampCurrent, TypeSamplingScheme sampScheme) {
	int i, ind = 0;
	for(i=0; i<sampScheme.sizeParamBirth; i++)
		x[ind++] = sampCurrent.currentBirth[i];
	for(i=0; i<sampScheme.sizeParamDeath; i++)
		x[ind++] = sampCurrent.currentDeath[i];
	for(i=0; i<sampScheme.sizeParamFossil; i++)
		x[ind++] = sampCurrent.currentFossil[i];
	for(i=0; i<sampScheme.sizeParamSampling; i++)
		x[ind++] = sampCurrent.currentSampling[i];
}

double toMaximizeSetTreeFossil(unsigned n, const double *x, double *grad, void *data) {
	double tmp, old;
//	old = ((TypeMinimizationSetTreeFossilData*)data)->likelihood(((TypeMinimizationSetTreeFossilData*)data)->tree, ((TypeMinimizationSetTreeFossilData*)data)->nTree, ((TypeMinimizationSetTreeFossilData*)data)->fos, &(((TypeMinimizationSetTreeFossilData*)data)->param));
	setSamplingCurrentFromVector(&(((TypeMinimizationSetTreeFossilData*)data)->sampCurrent), x, ((TypeMinimizationSetTreeFossilData*)data)->sampScheme);
	setParamFromSamplingCurrent(&(((TypeMinimizationSetTreeFossilData*)data)->param), ((TypeMinimizationSetTreeFossilData*)data)->sampCurrent, ((TypeMinimizationSetTreeFossilData*)data)->sampScheme);
	tmp = ((TypeMinimizationSetTreeFossilData*)data)->likelihood(((TypeMinimizationSetTreeFossilData*)data)->tree, ((TypeMinimizationSetTreeFossilData*)data)->nTree, ((TypeMinimizationSetTreeFossilData*)data)->fos, &(((TypeMinimizationSetTreeFossilData*)data)->param));
//	printf("like %.4le -> %.4le\n", old, tmp);
	return tmp;
}

int minimizePiecewiseParamFromSetTreeFossil(TypeLikelihoodSetTreeFosFunction *f, TypeTree **tree, int nTree, TypeFossilFeature *fos,TypeSamplingCurrent sampCurrentInit,  TypeSamplingScheme sampScheme, TypeNLOptOption *option, TypeEstimation *estim) {
	double *lb, *lu, *x, minLikelihood;
	nlopt_opt opt;
	TypeMinimizationSetTreeFossilData data;
	int result, i, ind, nParam;
	//fprintNLoptOptionTag(stdout, option);
	//printf("\n");
	nParam = sampScheme.sizeParamBirth+sampScheme.sizeParamDeath+sampScheme.sizeParamFossil+sampScheme.sizeParamSampling;
	
	lb = (double*) malloc(nParam*sizeof(double));
	lu = (double*) malloc(nParam*sizeof(double));
	x = (double*) malloc(nParam*sizeof(double));
	ind = 0;
	for(i=0; i<sampScheme.sizeParamBirth; i++) {
		lb[ind] = 0.;
		lu[ind] = option->supSpe;
		ind++;
	}
	for(i=0; i<sampScheme.sizeParamDeath; i++) {
		lb[ind] = 0.;
		lu[ind] = option->supExt;
		ind++;
	}
	for(i=0; i<sampScheme.sizeParamFossil; i++) {
		lb[ind] = 0.;
		lu[ind] = option->supFos;
		ind++;
	}
	for(i=0; i<sampScheme.sizeParamSampling; i++) {
		lb[ind] = 0.;
		lu[ind] = 1.;
		ind++;
	}
	data.tree = tree;
	data.nTree = nTree;
	data.fos = fos;
	data.param = getPiecewiseParamFromSamplingScheme(sampScheme);
	data.sampScheme = sampScheme;
	data.sampCurrent = getSamplingCurrent(sampScheme);
	data.likelihood = f;
	opt = nlopt_create(NLOPT_ALGO, nParam); /* algorithm and dimensionality */
//	nlopt_add_inequality_constraint(opt, ratesConstraint, NULL, 1e-8);
	nlopt_set_lower_bounds(opt, lb);
	nlopt_set_upper_bounds(opt, lu);
	nlopt_set_max_objective(opt, toMaximizeSetTreeFossil, &data);
	nlopt_set_xtol_abs1(opt, option->tolOptim);
	nlopt_set_maxeval(opt, option->maxIter);
	setVectorFromSamplingCurrent(x, sampCurrentInit, sampScheme);
	if(((result = nlopt_optimize(opt, x, &minLikelihood)) >= 0)) {
		estim->logLikelihood = minLikelihood;
		setSamplingCurrentFromVector(&(data.sampCurrent), x, data.sampScheme);
		setParamFromSamplingCurrent(&(estim->param), data.sampCurrent, data.sampScheme);
		printf("result (%d) %.4le\n", result, minLikelihood);
	} else {
		switch(result) {
			case NLOPT_SUCCESS:
				printf("Success\n");
				break;
			case NLOPT_STOPVAL_REACHED:
				printf("Stop val reached\n");
				break;
			case NLOPT_FTOL_REACHED:
				printf("Ftol reached\n");
				break;
			case NLOPT_XTOL_REACHED:
				printf("Xtol reached\n");
				break;
			case NLOPT_MAXEVAL_REACHED:
				printf("Max eval reached\n");
				break;
			case NLOPT_MAXTIME_REACHED:
				printf("Max time reached\n");
				break;
			case NLOPT_FAILURE:
				printf("General failure\n");
				break;
			case NLOPT_INVALID_ARGS:
				printf("Invalid args\n");
				break;
			case NLOPT_OUT_OF_MEMORY:
				printf("Out of memory\n");
				break;
			case NLOPT_ROUNDOFF_LIMITED:
				printf("Roundoff limited\n");
				break;
			case NLOPT_FORCED_STOP:
				printf("Forced stop\n");
				break;
			default:
				printf("failure %d\n", result);
		}
		printf("result (%d)\n", result);
	}
	estim->logLikelihood = estim->logLikelihood;
	free((void*)lb);
	free((void*)lu);
	free((void*)x);
	freeSamplingCurrent(&(data.sampCurrent));
	freePiecewiseParam(&(data.param));
	nlopt_destroy(opt);
	return result;
}
