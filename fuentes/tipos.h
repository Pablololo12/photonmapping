/*
 * Autores: Pablo Hernandez Almudi y Mario Arcega Sanjuan
 */

#ifndef TIPOS_H
#define TIPOS_H

typedef struct punto{
	double x;
	double y;
	double z;
}punto;

typedef struct vector{
	double x;
	double y;
	double z;
}vector;

typedef struct color{
	double r;
	double g;
	double b;
}color;

typedef struct propiedades{
	color * color;
	color * Krfl;
	color * Krfr;
	double indice_ref;
	double ks;
	double alpha;
}propiedades;

typedef struct lista{
	punto * punto;
	vector * normales;
	propiedades * propiedades;
	double radio;
	struct lista * l;
}lista;

typedef struct luces{
	punto * punto;
	color * color;
	struct luces * l;
}luces;

struct kd_node_t{
	double x[MAX_DIM];
	double dist;
	color * color;
	vector * vector;
	struct kd_node_t *left, *right;
};

#endif