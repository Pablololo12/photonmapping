
#ifndef _DIRECT
#define _DIRECT

#include "tipos.h"
#include "indirecta.h"
#include "common_functions.h"

/*
 * Método para calcular la luz directa
 */
color luz_directa(punto, lista *, luces *, vector, vector);

/*
 * Método para calcular la luz de un pixel
 */
color calcular_luz(vector, punto, int);

/*
 * Método para calcular el color de reflexión, los vectores deben estar normalizados
 */
color reflection (punto *, vector *, vector *, int);

/*
 * Método para calcular la refraccion
 */
color refraction (lista *, punto *, vector *, vector *, int, int);

#endif