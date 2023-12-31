#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "Utils.h"
#include "Tree.h"
#include "Fossil.h"
#include "FossilInt.h"
#include "Distribution.h"
#include "DrawTreeCairo.h"
#include "DrawTreePSTricks.h"
#include "DrawTreeTikz.h"
#include "DrawTreeGeneric.h"
#include "DrawFossilInt.h"
#include "DrawDensity.h"
#include "DrawTimePosition.h"
#include "DrawDensity.h"
//./drsp -a ../../data/map_nc.txt -y 20 -x 6 ~/Dropbox/Extinction/data/Pipid_analysenoncontrainte_consensusmajoritaire.tre ~/Dropbox/Extinction/data/Pipoid_occurences.csv
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
#define PREFIX "outTree"
#define MAX_TRIALS 1000
#define MAX_NODE 100

#define HELPMESSAGE "--------------------------\n\nNAME\n	draw - Drawing a single tree with the fossil ages\n	\nSYNOPSIS\n	draw [OPTIONS] <input Tree File> <input Fossil File>\n\nDESCRIPTION\n	Output a figure in various graphic formats of the tree in <input Tree File> with the fossil ages of  <input Fossil File>\n\n	Options are\n	-z <input Tree File>\n		output the tree in text debug format in the console and exit \n	-o <origin> \n		set the origin time\n	-e <end> \n		set the end time\n	-x <number>\n		set the graphic format of the output (option is required if one wants a graphic output)\n			-x 1 -> pdf\n			-x 2 -> postscript\n			-x 3 -> png\n			-x 4 -> svg\n			-x 5 -> LaTeX (psTricks)\n			-x 6 -> LaTeX (TikZ)\n	-h\n		display help\n\n--------------------------\n"

typedef struct DATA_DRAW_FOSSIL_DENSITY {
	TypeDataDrawFossilInt *dataFossil;
	TypeDataDrawDensity *dataDensity;
} TypeDataDrawFossilDensity;


void drawFossilDensity(int n , double x, double y, TypeInfoDrawTreeGeneric *info, void *data) {
	drawFossilInt(n, x, y, info, (void*) ((TypeDataDrawFossilDensity*)data)->dataFossil);
	drawDensity(n, x, y, info, (void*) ((TypeDataDrawFossilDensity*)data)->dataDensity);
}

typedef struct MAP {
	char *key, *val;
} TypeMap;

static char **readList(FILE *f);
static TypeMap *readMap(FILE *f);
static int compareMap(const void* a, const void* b);



int main(int argc, char **argv) {	
	char *inputFileNameTree, *inputFileNameFossil = NULL, *tmp, outputFileNameG[SIZE_BUFFER_CHAR], *outputFileName, option[256], format = '1', *nameOutput = "/home/gilles/Dropbox/Extinction/data/fig.gnu", *prefixClade = "/home/gilles/Dropbox/Extinction/data/Pipid_trees_analysenoncontrainte.tre_", *suffixClade = ".txt", *prefixDist="/home/gilles/Dropbox/Extinction/data/Pipid_nc", *suffixDist=".csv";
	FILE *fi;
	int i, j, debug = 0;
	double minTimeIntInf = NO_TIME, maxTimeIntSup = 0., maxDisplayTime = 0., treeScaleStep=10., figwidth = 100.;
	TypeMap *map = NULL;
	for(i=0; i<256; i++)
		option[i] = 0;
	   
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['x']) {
			option['x'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%c", &format) == 1)
				i++;
			else
				error("a character is required after option -f");
		}
		if(option['a']) {
			option['a'] = 0;
			if((i+1)<argc && (fi = fopen(argv[++i], "r"))) {
				map = readMap(fi);
				int j;
				printf("Reading map file:\n");
				for(j=0; map[j].key!=NULL; j++)
					printf("%s -> %s\n", map[j].key, map[j].val);
				fclose(fi);
			} else {
				error("Error while reading %s.\n", argv[i]);
			}
		}
		if(option['z']) {
			option['z'] = 0;
			if((i+1)<argc && (fi = fopen(argv[++i], "r"))) {
				TypeTree *tree;
				tree = readTree(fi);
//				toBinary(tree);
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
				error("Error while reading %s.\n", argv[i]);
			}
			exit(0);
		}
		if(option['q']) {
			option['q'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &treeScaleStep) == 1)
				i++;
			else
				error("a number is required after option -t");
		}
		if(option['c']) {
			option['c'] = 0;
			if((i+1)<argc)
				prefixClade = argv[++i];
			else
				error("2 strings are expected after -c");
			if((i+1)<argc)
				suffixClade = argv[++i];
			else
				error("2 strings are expected after -c");
		}
		if(option['d']) {
			option['d'] = 0;
			if((i+1)<argc)
				prefixDist = argv[++i];
			else
				error("2 strings are expected after -d");
			if((i+1)<argc)
				suffixDist = argv[++i];
			else
				error("2 strings are expected after -d");
		}
		if(option['o']) {
			option['o'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &minTimeIntInf) == 1)
				i++;
			else
				error("1 values are expected after -o");
		}
		if(option['e']) {
			option['e'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &maxTimeIntSup) == 1)
				i++;
			else
				error("1 values are expected after -o");
		}
		if(option['m']) {
			option['m'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &(maxDisplayTime)) == 1)
				i++;
			else
				error("a number is expected after -m");
		}
		if(option['y']) {
			option['y'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &figwidth) == 1)
				i++;
			else
				error("a number is required after option -t");
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
		error("Please provide the name of a file containing a phylogenetic tree in Newick format\n");
	}
	if(i<argc) {
		inputFileNameFossil = argv[i++];
	} else {
		error("Warning: No name of file containing a fossil list was provided\n");
	}
	if((fi = fopen(inputFileNameTree, "r"))) {
		TypeTree *tree;
		TypeFossilIntFeature *fos;
		int n;
		TypeInfoDrawTreeGeneric info;
		TypeAdditionalDrawTreeGeneric add;
		TypeAdditionalDrawTreeGenericGeneral addG;
		TypeDataDrawFossilInt data;
		TypeDataDrawTimePosition timePos;
		TypeDataDrawDensity dataD;
		TypeDataDrawFossilDensity dataFD;
		TypeDistribution *d;
//printf("OK %s\n", inputFileNameTree);
		
        tree = readTree(fi);
        fclose(fi);
//		toBinary(tree);
printTreeDebug(stdout, tree->root, tree, tree->name);
		reorderTreeSize(tree);
		if(inputFileNameFossil != NULL) {
			FILE *ff;
			if((ff = fopen(inputFileNameFossil, "r"))) {
				fos = getFossilIntFeature(ff, tree->name, tree->size);
			} else {
				error("Cannot read %s\n", inputFileNameFossil);
			}
		} else
			fos = fosToFossilInt(tree);
		if(getMaxFossilIntTime(fos) > 0.)
			negateFossilInt(fos);
		fixStatus(tree, fos);
		double minFossilTime = getMinFossilIntTime(fos);
		if(minTimeIntInf == NO_TIME) {
			double tmp = getMinFossilIntTime(fos);
			if(tmp<0) {
				minTimeIntInf = 1.2*minFossilTime;
			} else {
				minTimeIntInf = 0.8*minFossilTime;
			}
		}
		if(minTimeIntInf > minFossilTime)
			minTimeIntInf = minFossilTime;
		tree->maxTime = maxTimeIntSup;
		tree->minTime = minTimeIntInf;
		tree->minTimeInt.inf = minTimeIntInf;
		tree->minTimeInt.sup = minTimeIntInf;
		if(tree->parent==NULL)
			tree->parent = getParent(tree);
		outputFileName = inputFileNameTree;
		d = (TypeDistribution*) malloc(tree->size*sizeof(TypeDistribution));
		for(n=0; n<tree->size; n++) {
			d[n].size = 0;
			d[n].item = NULL;
		}
		if(map != NULL) {
			int i;
			FILE *fo;
			if((fo = fopen(nameOutput, "w"))) {
				fprintf(fo, "set style fill transparent solid 0.25 border\nfilter(x,min,max) = (x > min && x < max) ? x : 1/0\n\nset terminal tikz color size 10,6 font \",6\" createstyle;\nset xtics 10\nunset ytics\nunset xlabel\nset key top left\nset output \"FigTikzSpe_dist.tex\"\nset xrange [175:25]\nplot");
				for(i=0; map[i].key!=NULL; i++)
					;
				qsort(map, i, sizeof(TypeMap), compareMap);
				for(i=0; map[i].key!=NULL; i++) {
					FILE *fi;
					char nameTmp[STRING_SIZE];
					int node;
					sprintf(nameTmp, "%s%s%s", prefixClade, map[i].key, suffixClade);
					if((fi = fopen(nameTmp, "r"))) {
						char **list = readList(fi);
						fclose(fi);
						node = getClade(list, tree);
	printf("Clade %s -> node %d\n", nameTmp, node);
						if(node != NOSUCH) {
							if(tree->name[node]!=NULL)
								free((void*)tree->name[node]);
							tree->name[node] = (char*) malloc((strlen(map[i].val)+1)*sizeof(char));
							strcpy(tree->name[node], map[i].val);
						}
						sprintf(nameTmp, "%s%s%s", prefixDist, map[i].key, suffixDist);
						if((fi = fopen(nameTmp, "r"))) {
							d[node] = readDistribution(fi);
							fclose(fi);
							if(i<8)
								fprintf(fo, " '%s' using (abs($1)):2 with lines lw 2 title \"%s\",", nameTmp, map[i].val);
							else
								fprintf(fo, " '%s' using (abs($1)):2 with lines lw 2 dashtype \".\" title \"%s\",", nameTmp, map[i].val);
						}
					}
				}
				fclose(fo);
			}
		}
		for(n=0; n<tree->size; n++) {
			if(tree->node[n].child == NOSUCH) {
						int f;
						double max;
				switch(fos->status[n]) {
					case contempNodeStatus:
						tree->time[n] = tree->maxTime;
					break;
					case unknownNodeStatus:
						//tree->time[n] = NO_TIME;
					break;
					case extinctNodeStatus:
					//case unknownNodeStatus:
						max = NEG_INFTY;
						for(f=fos->fossilInt[n]; f!=NOSUCH; f=fos->fossilIntList[f].prec)
							if(fos->fossilIntList[f].fossilInt.sup>max)
								max = fos->fossilIntList[f].fossilInt.sup; 
						tree->time[n] = max;
					break;
					default:
						error("Node %d has no status\n", n);
						return 1;
				}
			}
		}
		for(n=0; n<tree->size; n++) {
			if(d[n].size>0) {
				tree->time[n] = getMedianDens(d[n]);
				printf("n %d -> %lf\n", n, tree->time[n]);
			}
		}
		if(debug) {
			int n;
			if(tree->name == NULL) {
				tree->name = (char**) malloc(tree->size*sizeof(char*));
				for(n=0; n<tree->size; n++)
					tree->name[n] = NULL;
			}
			for(n=0; n<tree->size; n++) {
				char tmp[1000];
				if(tree->name[n] == NULL) {
					sprintf(tmp, "%d", n);
					tree->name[n] = (char*) malloc((strlen(tmp)+1)*sizeof(char));
					strcpy(tree->name[n], tmp);
				} else {
					sprintf(tmp, "%s (%d)", tree->name[n], n);
					free((void*)tree->name[n]);
					tree->name[n] = (char*) malloc((strlen(tmp)+1)*sizeof(char));
					strcpy(tree->name[n], tmp);
				}
			}
		}
		if((tmp = strrchr(outputFileName, '.')) != NULL)
			tmp[0] = '\0';
printf("OutputFile : %s\n", outputFileName);
		info.param.tmin = tree->minTime;
		info.param.tmax = maxDisplayTime;
		info.param.ratio = 0.55;
		info.param.scaleStep = treeScaleStep;
		info.param.width = figwidth;
		data.color = (TypeRGB) {.red = 0.63, .green = 0.32, .blue = 0.18};
		data.radius = 5.;
		data.alpha = 0.5;
		data.fos = fos;
		dataD.color = (TypeRGB) {.red = 0., .green = 0., .blue = 1.};
		dataD.alpha = 0.5;
		dataD.dens = d;
		dataFD.dataFossil = &data;
		dataFD.dataDensity = &dataD;
		add.data = &dataFD;
		add.draw = drawFossilDensity;
		timePos.color = (TypeRGB) {.red = 0., .green = 0., .blue = 1.};
		timePos.alpha = 0.5;
		timePos.pos.start = -278.;
		timePos.pos.end = -276.;
//		addG.data = &timePos;
		addG.data = NULL;
		addG.draw= drawTimePosition;
		tree->maxTime = info.param.tmax;
		switch(format) {
			case '1':
				sprintf(outputFileNameG, "%s_tree.pdf", outputFileName);
				setFunctPDF(&(info.funct));
				data.drawLineDot = drawLineDotCairo;
				dataD.fillPolygon = fillPolygonCairo;
				drawTreeFileGeneric(outputFileNameG, tree, &info, &add, &addG);
				break;
			case '2':
				sprintf(outputFileNameG, "%s_tree.ps", outputFileName);
				setFunctPS(&(info.funct));
				data.drawLineDot = drawLineDotCairo;
				dataD.fillPolygon = fillPolygonCairo;
				drawTreeFileGeneric(outputFileNameG, tree, &info, &add, &addG);
				break;
			case '3':
				sprintf(outputFileNameG, "%s_tree.png", outputFileName);
				setFunctPNG(&(info.funct));
				data.drawLineDot = drawLineDotCairo;
				dataD.fillPolygon = fillPolygonCairo;
				drawTreeFileGeneric(outputFileNameG, tree, &info, &add, &addG);
				break;
			case '4':
				sprintf(outputFileNameG, "%s_tree.svg", outputFileName);
				setFunctSVG(&(info.funct));
				data.drawLineDot = drawLineDotCairo;
				dataD.fillPolygon = fillPolygonCairo;
				drawTreeFileGeneric(outputFileNameG, tree, &info, &add, &addG);
				break;
			case '5':
				sprintf(outputFileNameG, "%s_tree_pst.tex", outputFileName);
				setFunctPSTricks(&(info.funct));
				data.drawLineDot = drawLineDotPSTricks;
				dataD.fillPolygon = fillPolygonPSTricks;
				drawTreeFileGeneric(outputFileNameG, tree, &info, &add, &addG);
				break;
			case '6':
				sprintf(outputFileNameG, "%s_tree_tikz.tex", outputFileName);
				setFunctTikz(&(info.funct));
				data.drawLineDot = drawLineDotTikz;
				dataD.fillPolygon = fillPolygonTikz;
				drawTreeFileGeneric(outputFileNameG, tree, &info, &add, &addG);
				break;
			default:
				;
		}
		freeTree(tree);
		freeFossilIntFeature(fos);
	} else {
		error("Cannot read %s or %s\n", inputFileNameTree, inputFileNameFossil);
	}
	return 0;
}

#define MAX_SIZE_TMP 50
#define INC_BUFFER 50
#define IS_SEP(c) (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == ';')

char **readList(FILE *f) {
	char c, tmp[MAX_SIZE_TMP+1], **list;
	int size, sizeBuffer;

	sizeBuffer = INC_BUFFER;
	list= (char**) malloc(sizeBuffer*sizeof(char*));
	size = 0;
	do {
		c = getc(f);
	} while(c!=EOF && IS_SEP(c)); 
	while(c != EOF) {
		int i;
		i = 0;
		while(i<MAX_SIZE_TMP && c !=EOF && !IS_SEP(c)) {
			tmp[i] = c;
			c = getc(f);
			i++;
		}
		tmp[i++] = '\0';
		if(i == MAX_SIZE_TMP) {
			fprintf(stderr, "Ident too long (%s) ...", tmp);
			exit(1);
		}
		if(i>1) {
			if(size >= sizeBuffer) {
				sizeBuffer += INC_BUFFER;
				list = (char**) realloc((void *) list, sizeBuffer*sizeof(char*));
			}
			list[size] = (char*) malloc((strlen(tmp)+1)*sizeof(char));
			strcpy(list[size], tmp);
			size++;
		}
		while(c!=EOF && IS_SEP(c))
			c=getc(f);
	}
	if(size >= sizeBuffer) {
		sizeBuffer += INC_BUFFER;
		list = (char**) realloc((void *) list, sizeBuffer*sizeof(char*));
	}
	list[size++] = NULL;
	return list;
}

int compareMap(const void* a, const void* b) {
    return strcmp(((TypeMap*)a)->val, ((TypeMap*)b)->val);
}

TypeMap *readMap(FILE *f) {
	char c, tmp[MAX_SIZE_TMP+1];
	int size, sizeBuffer;
	TypeMap *list;
	
	sizeBuffer = INC_BUFFER;
	list= (TypeMap*) malloc(sizeBuffer*sizeof(TypeMap));
	size = 0;
	do {
		c = getc(f);
	} while(c!=EOF && IS_SEP(c)); 
	while(c != EOF) {
		int i;
		i = 0;
		while(i<MAX_SIZE_TMP && c !=EOF && !IS_SEP(c)) {
			tmp[i] = c;
			c = getc(f);
			i++;
		}
		tmp[i++] = '\0';
		if(i == MAX_SIZE_TMP) {
			fprintf(stderr, "Ident too long (%s) ...", tmp);
			exit(1);
		}
		if(i>1) {
			if(size >= sizeBuffer) {
				sizeBuffer += INC_BUFFER;
				list = (TypeMap*) realloc((void *) list, sizeBuffer*sizeof(TypeMap));
			}
			list[size].key = (char*) malloc((strlen(tmp)+1)*sizeof(char));
			strcpy(list[size].key, tmp);
			while(c!=EOF && IS_SEP(c))
				c=getc(f);
			if(c==EOF) {
				fprintf(stderr, "Warning: key without val in map file\n");
			} else {
				i = 0;
				while(i<MAX_SIZE_TMP && c !=EOF && !IS_SEP(c)) {
					tmp[i] = c;
					c = getc(f);
					i++;
				}
				tmp[i++] = '\0';
				if(i == MAX_SIZE_TMP) {
					fprintf(stderr, "Ident too long (%s) ...", tmp);
					exit(1);
				}
				if(i>1) {
					if(size >= sizeBuffer) {
						sizeBuffer += INC_BUFFER;
						list = (TypeMap*) realloc((void *) list, sizeBuffer*sizeof(TypeMap));
					}
					list[size].val = (char*) malloc((strlen(tmp)+1)*sizeof(char));
					strcpy(list[size].val, tmp);
					size++;
				}
				while(c!=EOF && IS_SEP(c))
					c=getc(f);
			}
		}
		while(c!=EOF && IS_SEP(c))
			c=getc(f);
	}
	if(size >= sizeBuffer) {
		sizeBuffer += INC_BUFFER;
		list = (TypeMap*) realloc((void *) list, sizeBuffer*sizeof(TypeMap));
	}
	list[size++].key = NULL;
	list[size++].val = NULL;
	return list;
}
