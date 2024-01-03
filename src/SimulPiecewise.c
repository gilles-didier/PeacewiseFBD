#include "SimulPiecewise.h"
#include <math.h>
#include "MyR.h"
#include "MyRandom.h"
#include "Utils.h"
#include "Fossil.h"


#define INC_FOSSIL_ITEM 50
#define INFTY 9.9E99;

static char getType(double birth, double death, double fossil, void *rand_data);

TypeTree *simulTreePiecewise(TypePiecewiseModelParam *param, void *rand_data) {
    int *cur, ncur, ind, i;
	double time = 0.;
	TypeTree *tree;
	TypeFossilFeature *fos;
	if((cur = (int*) malloc((MAX_CURRENT+1)*sizeof(int))) == NULL)
		return NULL;
	tree = (TypeTree*) malloc(sizeof(TypeTree));
	tree->parent = NULL;
	tree->name = NULL;
	tree->comment = NULL;
	tree->sizeBuf = INC_SIZE;
	tree->node = (TypeNode*) malloc(tree->sizeBuf*sizeof(TypeNode));
	tree->time = (double*) malloc(tree->sizeBuf*sizeof(double));
	tree->maxTime = param->startTime[param->size];
	tree->minTime = param->startTime[0];
	tree->size = 1;
	tree->root = 0;
	tree->time[0] = INFTY;
	tree->node[0].child = NOSUCH;
	tree->node[0].sibling = NOSUCH;
	fos = (TypeFossilFeature*) malloc(sizeof(TypeFossilFeature));
	fos->sizeBuf = INC_FOSSIL_ITEM;
	fos->fossil = (int*) malloc(tree->sizeBuf*sizeof(int));
	fos->status = (TypeNodeStatus*) malloc(tree->sizeBuf*sizeof(TypeNodeStatus));
	fos->fossil[0] = NOSUCH;
	fos->status[0] = contempNodeStatus;
	fos->fossilList = (TypeFossilList*) malloc(fos->sizeBuf*sizeof(TypeFossilList));
	fos->size = 0;
//	fos->status = NULL;
	cur[0] = 0;
	ncur = 1;
	time = param->startTime[0];
	ind = 0;
	while(ind<param->size && ncur>0) {
		int j;
		while(time<param->startTime[ind+1] && ncur>0) {
			int type, which;
			double wait;
			type = getType(param->param[ind].birth, param->param[ind].death, param->param[ind].fossil, rand_data);
			which = my_rand_unif_range(rand_data, ncur-1);
			wait = -log(my_rand_unif_unit(rand_data))/(ncur*(param->param[ind].birth+param->param[ind].death+param->param[ind].fossil));
			time += wait;
			if(time < param->startTime[ind+1]) {
				switch(type) {
					case 'b':
						if(ncur > MAX_CURRENT) {
							warning("too much lineages generated during simulations (%d - max %d)", ncur, MAX_CURRENT);
							freeTree(tree);
							return NULL;
						}
						tree->time[cur[which]] = time;
						if((tree->size+1)>=tree->sizeBuf) {
							tree->sizeBuf += INC_SIZE;
							tree->node = (TypeNode*) realloc((void*)tree->node, tree->sizeBuf*sizeof(TypeNode));
							tree->time = (double*) realloc((void*)tree->time, tree->sizeBuf*sizeof(double));
							fos->fossil = (int*) realloc((void*)fos->fossil, tree->sizeBuf*sizeof(int));
							fos->status = (TypeNodeStatus*) realloc((void*)fos->status, tree->sizeBuf*sizeof(TypeNodeStatus));
						}
						fos->status[cur[which]] = noneNodeStatus;
						tree->node[cur[which]].child = tree->size;
						tree->time[tree->size] = INFTY;
						tree->node[tree->size].child = NOSUCH;
						tree->node[tree->size].sibling = tree->size+1;
						fos->fossil[tree->size] = NOSUCH;
						fos->status[tree->size] = contempNodeStatus;
						tree->time[tree->size+1] = INFTY;
						tree->node[tree->size+1].child = NOSUCH;
						tree->node[tree->size+1].sibling = NOSUCH;
						fos->fossil[tree->size+1] = NOSUCH;
						fos->status[tree->size+1] = contempNodeStatus;
						cur[which] = tree->size;
						cur[ncur] = tree->size+1;
						ncur++;
						tree->size += 2;
						break;
					case 'd':
						tree->time[cur[which]] = time;
						fos->status[cur[which]] = extinctNodeStatus;
						for(i=which+1; i<ncur; i++)
							cur[i-1] = cur[i];
						ncur--;
						break;
					case 'f':
						if(fos->size>=fos->sizeBuf) {
							fos->sizeBuf += INC_SIZE;
							fos->fossilList = (TypeFossilList*) realloc((void*)fos->fossilList, fos->sizeBuf*sizeof(TypeFossilList));
						}
						fos->fossilList[fos->size].time = time;
						fos->fossilList[fos->size].prec = fos->fossil[cur[which]];
						fos->fossil[cur[which]] = fos->size;
						fos->size++;
						break;
					default:
						break;
				}
			}
		}
		time = param->startTime[ind+1];
		j = 0;
		if(time < param->startTime[param->size])
			for(i=0; i<ncur; i++)
				if(my_rand_bern(rand_data, param->param[ind].sampling))
					cur[j++] = cur[i];
				else {
					tree->time[cur[i]] = time;
					fos->status[cur[i]] = extinctNodeStatus;
				}
		else
			for(i=0; i<ncur; i++)
				if(my_rand_bern(rand_data, param->param[ind].sampling))
					cur[j++] = cur[i];
				else {
					if(fos->fossil[cur[i]] != NOSUCH)
						tree->time[cur[i]] = param->startTime[param->size]-(param->startTime[param->size]-fos->fossilList[fos->fossil[cur[i]]].time)/10000.;
					else
						tree->time[cur[i]] = param->startTime[param->size]-(param->startTime[param->size]-param->startTime[param->size-1])/10000.;
					fos->status[cur[i]] = extinctNodeStatus;
				}
		ncur = j;
		ind++;
	}
	for(i=0; i<ncur; i++)
		tree->time[cur[i]] = param->startTime[param->size];
	free((void*)cur);
	tree->sizeBuf = tree->size;
	if(tree->size) {
		tree->node = (TypeNode*) realloc((void*)tree->node, tree->sizeBuf*sizeof(TypeNode));
		tree->time = (double*) realloc((void*)tree->time, tree->sizeBuf*sizeof(double));
		fos->fossil = (int*) realloc((void*)fos->fossil, tree->sizeBuf*sizeof(int));
		fos->status = (TypeNodeStatus*) realloc((void*)fos->status, tree->sizeBuf*sizeof(TypeNodeStatus));
	} else {
		free((void*)tree->node);
		free((void*)tree->time);
		tree->node = NULL;
		tree->time = NULL;
		fos->fossil = NULL;
		fos->status = NULL;
	}
	if(fos->size>0) {
		fos->sizeBuf = fos->size;
		fos->fossilList = (TypeFossilList*) realloc((void*)fos->fossilList, fos->sizeBuf*sizeof(TypeFossilList));
	} else {
		fos->sizeBuf = 0;
		free((void*)fos->fossilList);
		fos->fossilList = NULL;
	}
    tree->info = (void*) fos;
	return tree;
}


/* return  a random type of event wrt the rates*/
char getType(double birth, double death, double fossil, void *rand_data) {
	double uni = my_rand_unif_unit(rand_data);
	if(uni<birth/(birth+death+fossil)) {
		return 'b';
	} else {
			if(uni<(birth+death)/(birth+death+fossil))
			return 'd';
		else
			return 'f';
	}
}
