/*
 * Autores: Pablo Hernandez Almudi y Mario Arcega Sanjuan
 */

#ifndef TIPOS_H
#define TIPOS_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <signal.h>
#include <pthread.h>

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

typedef struct photon{
	color * color;
	vector * direccion;
}photon;

//Variables globales
#define NUM_RAYOS 1000.0
#define POT_INDIRECTA 10.0
#define RECURSIONES 5
#define NUM_THREADS 1
#define EPSILON 0.000001

extern volatile sig_atomic_t progreso;
extern volatile sig_atomic_t porcentaje;

extern unsigned char * img_buff;


// Valores de resolución y posición de la cámara
extern int ancho;
extern int alto;
extern punto camara;
extern int incrementador;
extern int total_insertado;

// Nombres por defecto de los ficheros
extern char * img;
extern char * scn;

// Punteros a las listas donde se almacenan las esferas y los focos
extern lista * l;
extern luces * lights;

extern struct kdtree * mapaFotones;

#endif