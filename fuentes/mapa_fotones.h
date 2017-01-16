#ifndef _MAPA_FOTONES_
#define _MAPA_FOTONES_

#include "tipos.h"
#include "kdtree.h"
#include "common_functions.h"
#include "direct.h"

vector global_desde_local(vector, vector, vector, vector);

double acumulativa_inversa_inclinacion(double);

double acumulativa_inversa_acimut(double);

int obtener_rayo(vector, vector *);

int enviar_photon(vector, punto, color, int);

int mapa_fotones();

#endif