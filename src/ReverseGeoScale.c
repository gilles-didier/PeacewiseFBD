#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>


#define STRING_SIZE 300
#define PREFIX "table"
#define HELPMESSAGE "\n\nusage: geos [options] <input file> [<output file>]\n\nEstimate the diversification rates of the tree contained in the input file.\nThe input file has to be in Newick format with special tags for fossils ages and origin and end of the diversification, \nit returns a text report with the estimates.\n\nOptions are:\n\t-o <options file name>\tload the settings of the optimizer. <options file name> has to be in the format:\n\t\t:SPE [0;1] :EXT [0;1] :FOS [0:1] :TRI 10 :TOL 1.E-7 :ITE 500\n\t-h\tdisplay help\n\n"

typedef struct RGB_COLOR {
	int red, green, blue;
} TypeRGB;

/*name[i] between bound[i] bound[i+1]*/
typedef struct GEO_SCALE {
	int size;
	double *bound;
	TypeRGB *color;
	char **name;
} TypeGeoScale;

static int blackTextColor(TypeRGB back);
static char *getFontSize(char fs);
static void error(const char *message, ...);
static void mallocGeoScale(TypeGeoScale *g, int size);
static void reallocGeoScale(TypeGeoScale *g, int size);
static TypeGeoScale readGeoScale(FILE *f);
static void fprintGeoScale(FILE *f, TypeGeoScale g);
static void fprintReverseGeoScale(FILE *f, TypeGeoScale g);

int main(int argc, char **argv) {	
	char *inputFileNameGeoScale, outputFileNameGeoScale[STRING_SIZE], option[256], fontSize;
	FILE *fi, *fo;
	int i, j;
	
	for(i=0; i<256; i++)
		option[i] = 0;
	   
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['h']) {
			option['h'] = 0;
			error("%s\n", HELPMESSAGE);
		}
	}
	if(i<argc) {
		inputFileNameGeoScale = argv[i++];
	} else
		error("Please provide the name of a file containing a phylogenetic tree in Newick format\n");
	if((fi = fopen(inputFileNameGeoScale, "r"))) {
		TypeGeoScale g;
		g = readGeoScale(fi);
//		fprintGeoScale(stdout, g);
		sprintf(outputFileNameGeoScale, "%s_rev.txt", inputFileNameGeoScale);
		if((fo = fopen(outputFileNameGeoScale, "w"))) {
			fprintReverseGeoScale(fo, g);
			fclose(fo);
		} else
			error("Cannot write %s\n", outputFileNameGeoScale);
	} else
		error("Cannot read %s\n", inputFileNameGeoScale);
	return 0;
}

#define MAX_SIZE_TMP 50
#define INC_BUFFER 50
#define IS_SEP(c) (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == ';')

char *getFontSize(char fs) {
	switch(fs) {
		case 'H':
			return "Huge";
		case 'h':
			return "huge";
		case 'L':
			return "Large";
		case 'l':
			return "large";
		case 'n':
			return "normalsize"; 
		case 's':
			return "small";
		case 'f':
			return "footnotesize";
		case 'c':
			return "scriptsize";
		case 't':
			return "tiny";
		default:
			return "normalsize"; 
	}
	return "normalsize"; 
}

int blackTextColor(TypeRGB back) {
	double r, g, b, R, G, B, L;
	
	r = ((double)back.red)/255.;
	g = ((double)back.green)/255.;
	b = ((double)back.blue)/255.;
	if(r <= 0.03928)
		R = r/12.92;
	else
		R = pow((r+0.055)/1.055, 2.4);
	if(g <= 0.03928)
		G = g/12.92;
	else
		G = pow((g+0.055)/1.055, 2.4);
	if(b <= 0.03928)
		B = b/12.92;
	else
		B = pow((b+0.055)/1.055, 2.4);
	L = 0.2126*R+0.7152*G+0.0722*B;
	return(L > sqrt(1.05*0.05)-0.05);
}		

void mallocGeoScale(TypeGeoScale *g, int size) {
	g->size = size;
	g->bound = (double*) malloc((size+1)*sizeof(double));
	g->color = (TypeRGB*) malloc(size*sizeof(TypeRGB));
	g->name = (char**) malloc(size*sizeof(char*));
}

void reallocGeoScale(TypeGeoScale *g, int size) {
	g->size = size;
	g->bound = (double*) realloc((void*)g->bound, (size+1)*sizeof(double));
	g->color = (TypeRGB*) realloc((void*)g->color, size*sizeof(TypeRGB));
	g->name = (char**) realloc((void*)g->name, size*sizeof(char*));
}

TypeGeoScale readGeoScale(FILE *f) {
	char c, tmp[MAX_SIZE_TMP+1];
	int size, sizeBuffer;
	TypeGeoScale g; 
	sizeBuffer = INC_BUFFER;
	mallocGeoScale(&g, sizeBuffer);
	size = 0;
	do {
		c = getc(f);
	} while(c!=EOF && IS_SEP(c)); 
	while(c != EOF) {
		int i;
		for(i=0; i<MAX_SIZE_TMP && c !=EOF && !IS_SEP(c); i++) {
			tmp[i] = c;
			c = fgetc(f);
		}
		tmp[i++] = '\0';
		if(i == MAX_SIZE_TMP)
			error("Ident too long (%s) ...\n", tmp);
		if(i>1) {
			if(size >= sizeBuffer) {
				sizeBuffer += INC_BUFFER;
				reallocGeoScale(&g, sizeBuffer);
			}
			g.bound[size] = atof(tmp);
		} else
			error("Empty ident...\n");
		while(c!=EOF && IS_SEP(c))
			c=getc(f);
		if(c != EOF) {
			i = 0;
			if(c == '\'') {
				c = fgetc(f);
				for(i=0; i<MAX_SIZE_TMP && c !=EOF && c != '\''; i++) {
					tmp[i] = c;
					c = fgetc(f);
				}
				if(c == '\'')
					c = fgetc(f);
			} else {
				for(i=0; i<MAX_SIZE_TMP && c !=EOF && !IS_SEP(c); i++) {
					tmp[i] = c;
					c = fgetc(f);
				}
			}
			tmp[i++] = '\0';
			if(i == MAX_SIZE_TMP)
				error("Ident too long (%s) ...\n", tmp);
			if(i>1) {
				if(size >= sizeBuffer) {
					sizeBuffer += INC_BUFFER;
					reallocGeoScale(&g, sizeBuffer);
				}
				g.name[size] = (char*) malloc(i*sizeof(char));
				strcpy(g.name[size], tmp);
			} else
				error("Empty ident...\n");
		}
		while(c!=EOF && IS_SEP(c))
			c=getc(f);
		if(c != EOF) {
			i = 0;
			while(i<MAX_SIZE_TMP && c !=EOF && !IS_SEP(c)) {
				tmp[i] = c;
				c = getc(f);
				i++;
			}
			tmp[i++] = '\0';
			if(i == MAX_SIZE_TMP)
				error("Ident too long (%s) ...\n", tmp);
			if(i>1) {
				if(size >= sizeBuffer) {
					sizeBuffer += INC_BUFFER;
					reallocGeoScale(&g, sizeBuffer);
				}
				g.color[size].red = atoi(tmp);
			} else
				error("Empty color 1...\n");
		}
		while(c!=EOF && IS_SEP(c))
			c=getc(f);
		if(c != EOF) {
			i = 0;
			while(i<MAX_SIZE_TMP && c !=EOF && !IS_SEP(c)) {
				tmp[i] = c;
				c = getc(f);
				i++;
			}
			tmp[i++] = '\0';
			if(i == MAX_SIZE_TMP)
				error("Ident too long (%s) ...\n", tmp);
			if(i>1) {
				if(size >= sizeBuffer) {
					sizeBuffer += INC_BUFFER;
					reallocGeoScale(&g, sizeBuffer);
				}
				g.color[size].green = atoi(tmp);
			} else
				error("Empty color 2...\n");
		}
		while(c!=EOF && IS_SEP(c))
			c=getc(f);
		if(c != EOF) {
			i = 0;
			while(i<MAX_SIZE_TMP && c !=EOF && !IS_SEP(c)) {
				tmp[i] = c;
				c = getc(f);
				i++;
			}
			tmp[i++] = '\0';
			if(i == MAX_SIZE_TMP)
				error("Ident too long (%s) ...\n", tmp);
			if(i>1) {
				if(size >= sizeBuffer) {
					sizeBuffer += INC_BUFFER;
					reallocGeoScale(&g, sizeBuffer);
				}
				g.color[size].blue = atoi(tmp);
				size++;
			} else
				error("Empty color 3...%d (%d %d)\n", size, g.color[size].red, g.color[size].green);
		}
		while(c!=EOF && IS_SEP(c))
			c=getc(f);
	}
	reallocGeoScale(&g, size);
	return g;
}

void fprintGeoScale(FILE *f, TypeGeoScale g) {
	int i;
	for(i=0; i<g.size; i++)
		fprintf(f, "%lf\n%s\n%d %d %d\n", g.bound[i], g.name[i], g.color[i].red, g.color[i].green, g.color[i].blue);
	fprintf(f, "%lf\n", g.bound[g.size]);
}

void fprintReverseGeoScale(FILE *f, TypeGeoScale g) {
	int i;
	fprintf(f, "%lf\n", g.bound[g.size]);
	for(i=g.size-1; i>=0; i--)
		fprintf(f, "'%s'\n%d %d %d\n%lf\n", g.name[i], g.color[i].red, g.color[i].green, g.color[i].blue, g.bound[i]);
}
	
void error(const char *message, ...) {
	va_list args;

	va_start(args, message);
	vprintf(message, args);
	va_end(args);
	exit(1);
}