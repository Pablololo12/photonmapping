/*
 * Autores: Pablo Hernandez Almudi y Mario Arcega Sanjuan
 */

#ifndef COMMON_
#define COMMON_

#include "tipos.h"

/*
 * Método para parsear el fichero de la escena
 */
int parserOBJ(FILE *);

/*
 * Método para parsear el fichero de la escena
 */
int parser(FILE *);

/*
 * Método para normalizar un vector
 */
int normalizar(vector *);

/*
 * Método para calcular el producto escalar
 */
double dotproduct (vector *, vector *);

/*
 * Método para calcular el producto vectorial
 */
int crossproduct (vector *, vector *, vector *);

/*
 * Metodo para intersectar un triangulo
 * Extraido de: https://en.wikipedia.org/wiki/Möller–Trumbore_intersection_algorithm
 */
double toca_triangulo(punto, vector, punto, punto, punto);

/*
 * Método para descubrir si un rayo toca una esfera
 * 	el vector vector debe estar normalizado
 */
double toca_esfera(punto, vector, punto, double, char *);

#endif