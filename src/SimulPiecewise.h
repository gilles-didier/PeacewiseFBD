#ifndef SimulPiecewiseF
#define SimulPiecewiseF


#include "Tree.h"
#include "PiecewiseModel.h"

/*simulate a random tree with specified birth and death rates and fossil finds on this tree with rate "fossil"*/
TypeTree *simulTreePiecewise(TypePiecewiseModelParam *param, void *rand_data);

#endif
