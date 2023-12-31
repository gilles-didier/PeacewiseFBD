#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "Utils.h"
#include "Tree.h"
#include "Fossil.h"
#include "DrawTreeCairo.h"
#include "DrawTreePSTricks.h"
#include "DrawTreeTikz.h"
#include "DrawTreeGeneric.h"
#include "DrawFossil.h"


#ifdef DO_PS
#endif

#define STRING_SIZE 300


void drawTreeFossil(char *outputFileName, char format, double figwidth, TypeTree *tree, TypeFossilFeature *fos) {
	char outputFileNameG[STRING_SIZE];
	TypeInfoDrawTreeGeneric info;
	TypeAdditionalDrawTreeGeneric add;
	TypeDataDrawFossil data;
	add.draw = drawFossil;
	data.color = (TypeRGB) {.red = 0.3, .green = 0., .blue = 0.};
	data.radius = 3.;
	data.alpha = 0.75;
	info.param.tmin = tree->minTime;
	info.param.tmax = tree->maxTime;
	info.param.scaleStep = 1.;
	info.param.width = figwidth;
	data.fos = fos;
	add.data = (void*) &data;
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
}
