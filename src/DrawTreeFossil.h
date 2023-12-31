#ifndef DrawTreeFossilF
#define DrawTreeFossilF
#include <stdio.h>
#include "Tree.h"
#include "Fossil.h"



#ifdef __cplusplus
extern "C" {
#endif

void drawTreeFossil(char *outputFileName, char format, double figwidth, TypeTree *tree, TypeFossilFeature *fos);

#ifdef __cplusplus
}
#endif

#endif
