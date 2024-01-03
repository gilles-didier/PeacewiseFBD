#ifndef FBDDensityF
#define FBDDensityF

#include <stdlib.h>
#include <stdio.h>
#include "Utils.h"
#include "Tree.h"
#include "Fossil.h"
#include "PiecewiseModel.h"

typedef enum {
	FBD_negative_constraint_type = 0, /*div time is anterior to time*/
	FBD_positive_constraint_type /*div time is posterior to time*/
} TypeFBDConstraintTimeType;



typedef struct FBD_TIME_CONSTRAINT {
	double time;
	TypeFBDConstraintTimeType type;
	int node;
} TypeFBDTimeConstraint;

#ifdef __cplusplus
extern "C" {
#endif
void testFBD(int k, int l, double tstartA, double tstartB, double tend, TypePiecewiseModelParam *model);
void testFBDStd(int k, int l, double tstartA, double tstartB, double tend, TypePiecewiseModelParam *model);
double logLikelihoodTreeFossil(TypeTree *tree, TypeFossilFeature *fos, TypeFBDTimeConstraint *constraint, TypePiecewiseModelParam *model);

#ifdef __cplusplus
}
#endif


#endif
