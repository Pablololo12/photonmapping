
#ifndef _INDIRECTA_
#define _INDIRECTA_

#include "kdtree.h"
#include "direct.h"

int brdf(vector, vector, vector, lista *, color *);

color luz_indirecta(vector pixel, punto lugar, vector normal, lista * objeto);

#endif