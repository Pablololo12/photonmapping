/*
 * Autores: Pablo Hernandez Almudi y Mario Arcega Sanjuan
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <signal.h>
#include <pthread.h>
#include "tipos.h"
#include "kdtree.h"

#ifdef OPENGL
	#ifdef __APPLE__
	#include <GLUT/glut.h>
	#else
	#include <GL/glut.h>
	#endif
#endif

#define NUM_RAYOS 1000.0
#define RECURSIONES 5
#define NUM_THREADS 1

volatile sig_atomic_t progreso=0;
volatile sig_atomic_t porcentaje=0;

unsigned char * img_buff;


color calcular_luz(vector pixel, punto cam, int recursivo);

// Valores de resolución y posición de la cámara
int ancho=500;
int alto=500;
punto camara={0.5,0.5,-3.0};
int incrementador;
int total_insertado=0;

double EPSILON = 0.000001;

// Nombres por defecto de los ficheros
char * img={"imagen.ppm"};
char * scn={"escenas/imagen.scn"};

// Punteros a las listas donde se almacenan las esferas y los focos
lista * l = NULL;
luces * lights = NULL;

struct kdtree * mapaFotones;


/*
 * Método para calcular el producto escalar
 */
double dotproduct (vector *a, vector *b){
	return (a->x * b->x) + (a->y * b->y) + (a->z * b->z);
}

/*
 * Método para calcular el producto vectorial
 */
int crossproduct (vector * a, vector * b, vector * c){
	c->x = a->y*b->z - a->z*b->y;
	c->y = a->z*b->x - a->x*b->z;
	c->z = a->x*b->y - a->y*b->x;
	return 1;
}

/*
 * Método para parsear el fichero de la escena
 */
int parserOBJ(FILE * f)
{
	punto vertices[50000];
	vector normales[50000];
	int i=0;
	int d=0;
	// Iteradores para ir guardando los datos
	lista * aux = NULL;
	// variables auxiliares
	double x,y,z;
	int a,b,c,dummy,a2,b2,c2;
	char buffdummy[100];
	// Se lee hasta alcanzar el final del fichero

	while(feof(f)==0)
	{		
		char op=fgetc(f);
		if(op=='\n') op=fgetc(f);
		if(op=='v'){
			op=fgetc(f);
			if(op==' '){
				fscanf(f,"%lf %lf %lf",&x,&y,&z);
				vertices[i].x=x+0.5;
				vertices[i].y=y+0.5;
				vertices[i].z=-z;//
				i++;
			} else if(op=='n'){
				fscanf(f,"%lf %lf %lf",&x,&y,&z);
				normales[d].x=x;
				normales[d].y=y;
				normales[d].z=-z;
				d++;
			}
		} else if(op=='f' && fgetc(f)==' '){
			fscanf(f,"%d/%d/%d %d/%d/%d %d/%d/%d",&a,&dummy,&a2,&b,&dummy,&b2,&c,&dummy,&c2);
			a--;b--;c--;
			a2--;b2--;c2--;
			lista * la = calloc(1,sizeof(lista));
			la->radio=-1.0;
			la->punto=calloc(3,sizeof(punto));
			la->normales=calloc(3,sizeof(vector));
			la->propiedades=calloc(1,sizeof(propiedades));
			la->propiedades->color=calloc(1,sizeof(color));
			la->propiedades->Krfl=calloc(1,sizeof(color));
			la->propiedades->Krfr=calloc(1,sizeof(color));

			la->punto[0].x=vertices[a].x; la->punto[0].y=vertices[a].y; la->punto[0].z=vertices[a].z;
			la->punto[1].x=vertices[b].x; la->punto[1].y=vertices[b].y; la->punto[1].z=vertices[b].z;
			la->punto[2].x=vertices[c].x; la->punto[2].y=vertices[c].y; la->punto[2].z=vertices[c].z;

			la->normales[0].x=normales[a2].x; la->normales[0].y=normales[a2].y; la->normales[0].z=normales[a2].z;
			la->normales[1].x=normales[b2].x; la->normales[1].y=normales[b2].y; la->normales[1].z=normales[b2].z;
			la->normales[2].x=normales[c2].x; la->normales[2].y=normales[c2].y; la->normales[2].z=normales[c2].z;

			la->propiedades->color->r=0.5; la->propiedades->color->g=0.5; la->propiedades->color->b=0.5;
			la->propiedades->Krfl->r=0.0; la->propiedades->Krfl->g=0.0; la->propiedades->Krfl->b=0.0;
			la->propiedades->Krfr->r=0.0; la->propiedades->Krfr->g=0.0; la->propiedades->Krfr->b=0.0;
			la->propiedades->indice_ref=0.0;
			la->propiedades->ks=0.5;
			la->propiedades->alpha=1;
			if(l==NULL){
				l=la;
			} else{
				aux->l=la;
			}
			aux=la;
		} else{
			fgets(buffdummy,99,f);
		}
	}
	return 0;
}

/*
 * Método para parsear el fichero de la escena
 */
int parser(FILE * f)
{
	// Iteradores para ir guardando los datos
	lista * aux = NULL;
	luces * aux2 = NULL;
	// variables auxiliares
	double x,y,z,radio,R,G,B,klr,klg,klb,krr,krg,krb,ind,ks,alpha;
	// Se lee la resolución
	fscanf(f, "%d %d",&ancho,&alto);
	// Se lee la posición de la cámara
	fscanf(f, "%lf %lf %lf",&x,&y,&z);
	camara.x=x; camara.y=y; camara.z=z;
	// Se lee hasta alcanzar el final del fichero
	while(feof(f)==0)
	{
		char c = fgetc(f);
		// Si la linea empieza por l es un foco de luz
		if(c=='l')
		{
			fscanf(f, "%lf %lf %lf %lf %lf %lf",&x,&y,&z,&R,&G,&B);

			luces * a = calloc(1, sizeof(luces));
			a->punto=calloc(1,sizeof(punto));
			a->color=calloc(1,sizeof(color));

			a->punto->x=x; a->punto->y=y; a->punto->z=z;
			a->color->r=R; a->color->g=G; a->color->b=B;

			if(lights==NULL){
				lights=a;
			} else{
				aux2->l=a;
			}
			aux2=a;
		}
		// Si la linea empieza por e se lee una esfera
		else if(c=='e')
		{
			fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&x,&y,&z,&radio,&R,&G,&B,&klr,&klg,&klb,&krr,&krg,&krb,&ind, &ks, &alpha);
			lista * a = calloc(1,sizeof(lista));
			a->radio=radio;
			a->punto=calloc(1,sizeof(punto));
			a->propiedades=calloc(1,sizeof(propiedades));
			a->propiedades->color=calloc(1,sizeof(color));
			a->propiedades->Krfl=calloc(1,sizeof(color));
			a->propiedades->Krfr=calloc(1,sizeof(color));

			a->punto->x=x; a->punto->y=y; a->punto->z=z;
			a->propiedades->color->r=R; a->propiedades->color->g=G; a->propiedades->color->b=B;
			a->propiedades->Krfl->r=klr; a->propiedades->Krfl->g=klg; a->propiedades->Krfl->b=klb;
			a->propiedades->Krfr->r=krr; a->propiedades->Krfr->g=krg; a->propiedades->Krfr->b=krb;
			a->propiedades->indice_ref=ind;
			a->propiedades->ks=ks;
			a->propiedades->alpha=alpha;
			if(l==NULL){
				l=a;
			} else{
				aux->l=a;
			}
			aux=a;
		} else if(c=='t')
		{
			double x2,y2,z2,x3,y3,z3;
			fscanf(f, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&x,&y,&z,&x2,&y2,&z2,&x3,&y3,&z3,&ks,&alpha,&R,&G,&B);
			lista * a = calloc(1,sizeof(lista));
			a->radio=-1.0;
			a->punto=calloc(3,sizeof(punto));
			a->propiedades=calloc(1,sizeof(propiedades));
			a->propiedades->color=calloc(1,sizeof(color));
			a->propiedades->Krfl=calloc(1,sizeof(color));
			a->propiedades->Krfr=calloc(1,sizeof(color));

			a->punto[0].x=x; a->punto[0].y=y; a->punto[0].z=z;
			a->punto[1].x=x2; a->punto[1].y=y2; a->punto[1].z=z2;
			a->punto[2].x=x3; a->punto[2].y=y3; a->punto[2].z=z3;
			a->normales=NULL;
			a->propiedades->color->r=R; a->propiedades->color->g=G; a->propiedades->color->b=B;
			a->propiedades->Krfl->r=0.0; a->propiedades->Krfl->g=0.0; a->propiedades->Krfl->b=0.0;
			a->propiedades->Krfr->r=0.0; a->propiedades->Krfr->g=0.0; a->propiedades->Krfr->b=0.0;
			a->propiedades->indice_ref=0.0;
			a->propiedades->ks=ks;
			a->propiedades->alpha=alpha;
			if(l==NULL){
				l=a;
			} else{
				aux->l=a;
			}
			aux=a;
		}
	}
	return 0;
}

/*
 * Método para normalizar un vector
 */
int normalizar(vector * vec)
{
	double tamanyo = sqrt(vec->x*vec->x + vec->y*vec->y + vec->z*vec->z);
	vec->x = (vec->x)/tamanyo;
	vec->y = (vec->y)/tamanyo;
	vec->z = (vec->z)/tamanyo;
	return 1;
}

/*
 * Metodo para intersectar un triangulo
 * Extraido de: https://en.wikipedia.org/wiki/Möller–Trumbore_intersection_algorithm
 */
double toca_triangulo(punto origen, vector vec, punto V1, punto V2, punto V3)
{
	vector e1 = {V2.x-V1.x, V2.y-V1.y, V2.z-V1.z};
	vector e2 = {V3.x-V1.x, V3.y-V1.y, V3.z-V1.z};

	vector P = {0.0,0.0,0.0};
	crossproduct(&vec, &e2, &P);

	double det = dotproduct(&e1, &P);

	if(det > -EPSILON && det < EPSILON) return -1.0;
	double inv_det = 1.0 / det;

	vector T = {origen.x-V1.x, origen.y-V1.y, origen.z-V1.z};

	double u = dotproduct(&T, &P) * inv_det;

	if(u < 0.0 || u > 1.0) return -1.0;

	vector Q = {0.0, 0.0, 0.0};
	crossproduct(&T, &e1, &Q);

	double v = dotproduct(&vec, &Q) * inv_det;

	if(v < 0.0 || u + v > 1.0 ) return -1.0;

	double t = dotproduct(&e2, &Q) * inv_det;

	if(t > EPSILON) {
		return t;
	}

	return -1.0;
}

/*
 * Método para descubrir si un rayo toca una esfera
 * 	el vector vector debe estar normalizado
 */
double toca_esfera(punto origen, vector vector, punto centro, double radio, char * dir)
{
	// Se calculan los operandos de la ecuacion cuadratica

	// A no hace falta calcularse si los vectores dirección son normales

	// Aqui se calcula la b (2*(O*D-D*C))
	double b = (origen.x*vector.x + origen.y*vector.y + origen.z*vector.z);
	b -= (centro.x*vector.x + centro.y*vector.y + centro.z*vector.z);
	b = b*2;

	// Aqui se calcula la c (O*O+C*C-2*O*C-r*r)
	double c = origen.x*origen.x + origen.y*origen.y + origen.z*origen.z;
	c += centro.x*centro.x + centro.y*centro.y + centro.z*centro.z;
	c -= 2*(centro.x*origen.x + centro.y*origen.y + centro.z*origen.z);
	c -= radio*radio;

	// Aqui calculo el interior de la raiz para ver si es negativo
	double aux = b*b - 4*c;

	if(aux<0.0)
	{
		// No se toca
		return -1.0;
	}

	// Se termina de realizar la ecuacion
	double t1 = 0.0 - b + sqrt(aux);
	t1 = t1 / 2.0;
	double t2 = 0.0 - b - sqrt(aux);
	t2 = t2 / 2.0;


	if(t1==t2 && t1>0.0)
	{
		// Se toca un único punto
		*dir=1;
		if(t1<EPSILON) return 0.0;
		return t1;
	} else if(t1>=0.0 && t2<0.0)
	{
		// Dentro del circulo
		*dir=0;
		if(t1<=EPSILON){*dir=1; return 0.0;} 
		return t1;
	} else if(t1>=0.0 && t2>=0.0)
	{
		//Fuera y se da el más cercano
		*dir=1;
		if(t2<=EPSILON){ *dir=0; return t1;}
		return t2;
	} else{
		//No se apunta a la esfera
		return -1.0;
	}
}

/*
 * Método para calcular la luz directa
 */
color luz_directa(punto esfera, lista * minimo, luces * luz, vector pixel, vector normal)
{
	// Con el punto calculamos Li, de momento un unico punto de luz directa
	vector inter;
	inter.x = luz->punto->x - esfera.x;
	inter.y = luz->punto->y - esfera.y;
	inter.z = luz->punto->z - esfera.z;

	double light = inter.x*inter.x + inter.y*inter.y + inter.z*inter.z;
	double dist_luz = sqrt(light);

	color power={0.0,0.0,0.0};
	power.r = luz->color->r / light;
	power.g = luz->color->g / light;
	power.b = luz->color->b / light;
	//if(light<1.0) printf("%f %f %f %f\n", power.r,power.g,power.b,light);
	
	normalizar(&inter);

	lista * aux = l;
	double dist;
	char dir=1;
	// Se comprueba si hay algun objeto entre el punto y la luz
	while(1)
	{
		if(aux->radio==-1.0)
		{
			dist = toca_triangulo(esfera, inter, aux->punto[0],aux->punto[1],aux->punto[2]);
		} else{
			dist = toca_esfera(esfera, inter,*aux->punto,aux->radio,&dir);
		}
		
		if(dist>0.0 && dist <= dist_luz ){
			color col={0.0,0.0,0.0};
			return col;
		}
		
		if(aux->l==NULL) break;
		aux = aux->l;
	}

	// Ahora calculamos la BRDF de phong
	// Primero obtenemos la normal, omega de o y de r
	
	vector omegao = {-pixel.x, -pixel.y, -pixel.z};
	normalizar(&omegao);
	
	//Obtenemos omega r
	double dotproduc = inter.x*normal.x + inter.y*normal.y + inter.z*normal.z;
	vector omegar;
	omegar.x = inter.x - 2*(inter.x - normal.x*dotproduc);
	omegar.y = inter.y - 2*(inter.y - normal.y*dotproduc);
	omegar.z = inter.z - 2*(inter.z - normal.z*dotproduc);
	normalizar(&omegar);


	double dotproductIntegral = dotproduct(&omegar, &normal);
	if (dotproductIntegral < 0) dotproductIntegral = -dotproductIntegral;

	// Por ultimo aplicamos phong
	double dotproductPhong = dotproduct(&omegao, &omegar);
	if (dotproductPhong < 0) dotproductPhong = -dotproductPhong;
	dotproductPhong = pow(dotproductPhong, minimo->propiedades->alpha);
	
	double especular = minimo->propiedades->ks * (minimo->propiedades->alpha + 2) / 2 * dotproductPhong;

	color rgb={0.0,0.0,0.0};
	rgb.r = power.r * (minimo->propiedades->color->r + especular) / M_PI * dotproductIntegral;
	rgb.g = power.g * (minimo->propiedades->color->g + especular) / M_PI * dotproductIntegral;
	rgb.b = power.b * (minimo->propiedades->color->b + especular) / M_PI * dotproductIntegral;

	if(rgb.r > 1.0) rgb.r=1.0;
	if(rgb.g > 1.0) rgb.g=1.0;
	if(rgb.b > 1.0) rgb.b=1.0;

	return rgb;
}



/*
 * Método para calcular el color de reflexión, los vectores deben estar normalizados
 */
color reflection (punto *point, vector *normal, vector *ray, int recursivo){

	// se calcula el factor para calcular el rayo reflectado y se calcula
	double factor = normal->x * -ray->x + normal->y * -ray->y + normal->z * -ray->z;
	vector reflection;
	reflection.x = -ray->x - 2 * (-ray->x - factor * normal->x);
	reflection.y = -ray->y - 2 * (-ray->y - factor * normal->y);
	reflection.z = -ray->z - 2 * (-ray->z - factor * normal->z);
	normalizar(&reflection);

	color col;
	col=calcular_luz(reflection,*point,recursivo);
	return col;
}

/*
 * Método para calcular la refraccion
 */
color refraction (lista *esfera, punto *point, vector *normal, vector *ray, int recursivo, int dentro){
	double n_refraction = esfera->propiedades->indice_ref; // Aire-Cristal 1.00/1.52=0.65
	if(dentro == 0) n_refraction = 1.0/n_refraction;
	double dot = -normal->x * ray->x + -normal->y * ray->y + -normal->z * ray->z;
	double k = 1.0 - n_refraction * n_refraction * (1.0 - dot * dot);
	vector refractado;
	if(k<0.0){
		refractado.x=0.0; refractado.y=0.0; refractado.z=0.0;
	}
	refractado.x = n_refraction * ray->x + (n_refraction * dot - sqrt(k)) * normal->x;
	refractado.y = n_refraction * ray->y + (n_refraction * dot - sqrt(k)) * normal->y;
	refractado.z = n_refraction * ray->z + (n_refraction * dot - sqrt(k)) * normal->z;
	normalizar(&refractado);

	color col; 
	col = calcular_luz(refractado,*point,recursivo);
	return col;
}

int brdf(vector omegao, vector inter, vector normal, lista * minimo, color * potencia)
{
	
	//Obtenemos omega r
	double dotproduc = inter.x*normal.x + inter.y*normal.y + inter.z*normal.z;
	vector omegar;
	omegar.x = inter.x - 2*(inter.x - normal.x*dotproduc);
	omegar.y = inter.y - 2*(inter.y - normal.y*dotproduc);
	omegar.z = inter.z - 2*(inter.z - normal.z*dotproduc);
	normalizar(&omegar);


	// Por ultimo aplicamos phong
	double dotproductPhong = dotproduct(&omegao, &omegar);
	if (dotproductPhong < 0) dotproductPhong = -dotproductPhong;
	dotproductPhong = pow(dotproductPhong, minimo->propiedades->alpha);
	
	double especular = dotproductPhong * minimo->propiedades->ks * (minimo->propiedades->alpha + 2) / 2;

	potencia->r = potencia->r * ((minimo->propiedades->color->r + especular) / M_PI);
	potencia->g = potencia->g * ((minimo->propiedades->color->g + especular) / M_PI);
	potencia->b = potencia->b * ((minimo->propiedades->color->b + especular) / M_PI);

	return 1;
}

color luz_indirecta(vector pixel, punto lugar, vector normal, lista * objeto)
{
	double pos[3] = {lugar.x,lugar.y,lugar.z};
	double i=0.05;
	struct kdres * busqueda=kd_nearest_range(mapaFotones, pos, i);

	color acumulador={0.0,0.0,0.0};
	double max_dist=0.0;

	while(kd_res_size(busqueda)<50)
	{
		i+=0.05;
		kd_res_free(busqueda);
		busqueda=kd_nearest_range(mapaFotones, pos, i);
	}

	//int num=kd_res_size(busqueda);

	while( !kd_res_end( busqueda ) ) {
		/* get the data and position of the current result item */
		photon * foton;
		foton = (photon *)kd_res_item( busqueda, pos );

		double acum = lugar.x-pos[0];
		double acum2 = acum*acum;
		acum = lugar.y-pos[1];
		acum2 += acum*acum;
		acum = lugar.z-pos[2];
		acum2 += acum*acum;

		acum = sqrt(acum2);
		
		color aux;
		aux.r = foton->color->r; aux.g = foton->color->g; aux.b = foton->color->b;
		brdf(pixel, *foton->direccion, normal, objeto, &aux);

		// Filtro gausiano
		/*acum2 = -1.953*(acum*acum)/(2*M_PI*M_PI);
		acum2 = pow(M_E,acum2);
		acum2 = (1 - acum2) / (1 - pow(M_E,-1.953));
		acum2 = 0.918*(1 - acum2);*/

		acum2 = 1 - (acum/(2.0*i));

		acumulador.r+=aux.r*acum2; acumulador.g+=aux.g*acum2; acumulador.b+=aux.b*acum2;

		/* go to the next entry */
		kd_res_next( busqueda );

		if(max_dist<acum)
			max_dist=acum;
	}

	kd_res_free(busqueda);

	max_dist=max_dist*max_dist;
	//max_dist=max_dist*max_dist*max_dist;
	double area = M_PI*max_dist;
	//double area = (4.0*M_PI*max_dist)/3.0;
	area = area * (1-(2/(3*2.0)));

	//acumulador.r=acumulador.r/area;
	//acumulador.g=acumulador.g/area;
	//acumulador.b=acumulador.b/area;

	acumulador.r=(acumulador.r/area)*10;
	acumulador.g=(acumulador.g/area)*10;
	acumulador.b=(acumulador.b/area)*10;

	if(acumulador.r>1.0) acumulador.r=1.0;
	if(acumulador.g>1.0) acumulador.g=1.0;
	if(acumulador.b>1.0) acumulador.b=1.0;

	return acumulador;
}

/*
 * Método para calcular la luz de un pixel
 */
color calcular_luz(vector pixel, punto cam, int recursivo)
{
	color rgb={0.0,0.0,0.0};
	// Primero buscamos el punto de intersección más cercano
	double min = 65535.0;
	double dist = 0.0;
	lista * aux = l;
	lista * minimo = NULL;
	char dir = 1;
	char dir_aux = 1;
	// Buscamos donde se choca el rayo
	while(1){
		if(aux->radio==-1.0)
		{
			dist = toca_triangulo(cam, pixel, aux->punto[0],aux->punto[1],aux->punto[2]);
		} else{
			dist = toca_esfera(cam,pixel,*aux->punto,aux->radio,&dir);
		}

		if(dist>0.0 && dist<min){
			min = dist;
			minimo = aux;
			dir_aux = dir;
		}
		if(aux->l==NULL) break;
		aux = aux->l;
	}
	if(min == 65535.0){
		color col={0.0,0.0,0.0};
		return col;
	}

	// Obtenemos las coordenadas del punto en el espacio
	punto esfera;
	esfera.x = cam.x + pixel.x*min;
	esfera.y = cam.y + pixel.y*min;
	esfera.z = cam.z + pixel.z*min;

	// Calculamos la normal dependiendo si se esta dentro o fuera del circulo
	vector * normal = calloc(1,sizeof(vector));
	if(minimo->radio>=0.0){
		if(dir_aux==1){
			normal->x = esfera.x - minimo->punto->x;
			normal->y = esfera.y - minimo->punto->y;
			normal->z = esfera.z - minimo->punto->z;
		} else{
			normal->x = minimo->punto->x - esfera.x;
			normal->y = minimo->punto->y - esfera.y;
			normal->z = minimo->punto->z - esfera.z;
		}
	} else{
		if(minimo->normales == NULL){
			vector e1 = {minimo->punto[1].x-minimo->punto[0].x, minimo->punto[1].y-minimo->punto[0].y, minimo->punto[1].z-minimo->punto[0].z};
			vector e2 = {minimo->punto[2].x-minimo->punto[0].x, minimo->punto[2].y-minimo->punto[0].y, minimo->punto[2].z-minimo->punto[0].z};
			crossproduct(&e1,&e2,normal);
			double dotp=dotproduct(normal, &pixel);
			if(dotp>0.0){normal->x=-normal->x;normal->y=-normal->y;normal->z=-normal->z;}
		
		}else{
			vector d1 = {minimo->punto[0].x-esfera.x, minimo->punto[0].y-esfera.y, minimo->punto[0].z-esfera.z};
			vector d2 = {minimo->punto[1].x-esfera.x, minimo->punto[1].y-esfera.y, minimo->punto[1].z-esfera.z};
			vector d3 = {minimo->punto[2].x-esfera.x, minimo->punto[2].y-esfera.y, minimo->punto[2].z-esfera.z};

			vector aux;
			crossproduct(&d1,&d2,&aux);
			double area3 = sqrt(aux.x*aux.x+aux.y*aux.y+aux.z*aux.z)/2.0;

			crossproduct(&d2,&d3,&aux);
			double area1 = sqrt(aux.x*aux.x+aux.y*aux.y+aux.z*aux.z)/2.0;

			crossproduct(&d3,&d1,&aux);
			double area2 = sqrt(aux.x*aux.x+aux.y*aux.y+aux.z*aux.z)/2.0;

			double sumatorio = area1+area2+area3;
			area1=area1/sumatorio; area2=area2/sumatorio; area3=area3/sumatorio;

			normal->x = minimo->normales[0].x*area1+minimo->normales[1].x*area2+minimo->normales[2].x*area3;
			normal->y = minimo->normales[0].y*area1+minimo->normales[1].y*area2+minimo->normales[2].y*area3;
			normal->z = minimo->normales[0].z*area1+minimo->normales[1].z*area2+minimo->normales[2].z*area3;

		}
	}
	normalizar(normal);

	// Obtenemos la luz directa
	luces * aux2 = lights;
	while(1){
		color col;
		col=luz_directa(esfera, minimo, aux2, pixel, *normal);
		rgb.r = rgb.r + col.r; rgb.g = rgb.g + col.g; rgb.b = rgb.b + col.b;
		if(aux2->l==NULL) break;
		aux2 = aux2->l;
	}
	
	// Obtenemos el porcentaje que corresponde por luz directa
	rgb.r = rgb.r * (1.0 - minimo->propiedades->Krfl->r - minimo->propiedades->Krfr->r);
	rgb.g = rgb.g * (1.0 - minimo->propiedades->Krfl->g - minimo->propiedades->Krfr->g);
	rgb.b = rgb.b * (1.0 - minimo->propiedades->Krfl->b - minimo->propiedades->Krfr->b);
	

	int flag = 0;
	if(recursivo>0){
		recursivo--;
		color color_reflexion;
		color color_refraccion;
				
		// Se calcula la reflexion solo si es necesario
		if(minimo->propiedades->Krfl->r!=0.0 && 
			minimo->propiedades->Krfl->g!=0.0 && 
			minimo->propiedades->Krfl->b!=0.0)
		{
			color_reflexion = reflection(&esfera, normal, &pixel, recursivo);
			rgb.r = rgb.r + color_reflexion.r * minimo->propiedades->Krfl->r;
			rgb.g = rgb.g + color_reflexion.g * minimo->propiedades->Krfl->g;
			rgb.b = rgb.b + color_reflexion.b * minimo->propiedades->Krfl->b;
			flag=1;
		}

		// Se calcula la refraccion solo si es necesario
		if(minimo->propiedades->Krfr->r!=0.0 && 
			minimo->propiedades->Krfr->g!=0.0 && 
			minimo->propiedades->Krfr->b!=0.0)
		{

			color_refraccion = refraction(minimo, &esfera, normal, &pixel, recursivo, dir_aux);
			rgb.r = rgb.r + color_refraccion.r * minimo->propiedades->Krfr->r;
			rgb.g = rgb.g + color_refraccion.g * minimo->propiedades->Krfr->g;
			rgb.b = rgb.b + color_refraccion.b * minimo->propiedades->Krfr->b;
			flag=1;
		}

	}

	if(flag==0){
		vector camara ={-pixel.x,-pixel.y,-pixel.z};
		normalizar(&camara);
		color indirecta = luz_indirecta(camara, esfera, *normal, minimo);

		rgb.r+=indirecta.r * (1.0 - minimo->propiedades->Krfl->r - minimo->propiedades->Krfr->r);
		rgb.g+=indirecta.g * (1.0 - minimo->propiedades->Krfl->g - minimo->propiedades->Krfr->g);
		rgb.b+=indirecta.b * (1.0 - minimo->propiedades->Krfl->b - minimo->propiedades->Krfr->b);
	}

	free(normal);

	//if(rgb.r > 1.0) rgb.r=1.0;
	//if(rgb.g > 1.0) rgb.g=1.0;
	//if(rgb.b > 1.0) rgb.b=1.0;
	return rgb;
}

vector global_desde_local(vector local, vector u, vector v, vector n){
	vector global = {0.0, 0.0, 0.0};
	global.x = u.x*local.x + v.x*local.y + n.x*local.z;
	global.y = u.y*local.x + v.y*local.y + n.y*local.z;
	global.z = u.z*local.x + v.z*local.y + n.z*local.z;
	return global;
}

double acumulativa_inversa_inclinacion(double x){
	return acos(sqrt(1.0-x));
}

double acumulativa_inversa_acimut(double x){
	return 2.0 * M_PI * x;
}

int obtener_rayo(vector n, vector * saliente)
{
	//primero se obtienen vectores perpendiculares para la geometria local
	double aleatorio=0.0;
	do{
		aleatorio = (double) (rand()%1000) / 1000.0;
	} while(aleatorio<=0.0 || aleatorio>=1.0);

	vector aleatoriov = {aleatorio, aleatorio, aleatorio};
	vector u;
	crossproduct(&n, &aleatoriov, &u);
	normalizar(&u);
	vector v;
	crossproduct(&n, &u, &v);
	normalizar(&v);

	//se eligen la inclinación y el acimut por montecarlo
	srand(clock());
	do{
		aleatorio = (double) (rand()%1000) / 1000.0;
	} while(aleatorio<0.0 || aleatorio>=1.0);
	double inclinacion = acumulativa_inversa_inclinacion(aleatorio);
	do{
		aleatorio = (double) (rand()%1000) / 1000.0;
	} while(aleatorio<0.0 || aleatorio>=1.0);
	double acimut = acumulativa_inversa_acimut(aleatorio);

	//vector reflejado en geometría local
	vector reflejado = {sin(inclinacion) * cos(acimut), sin(inclinacion) * sin(acimut), cos(inclinacion)};

	//vector reflejado en geometría global
	reflejado = global_desde_local(reflejado, u, v, n);
	normalizar(&reflejado);

	saliente->x=reflejado.x;
	saliente->y=reflejado.y;
	saliente->z=reflejado.z;

	return 1;
}

/*
 * Método para calcular la luz de un pixel
 */
int enviar_photon(vector pixel, punto cam, color potencia, int primera_vez)
{
	// Primero buscamos el punto de intersección más cercano
	double min = 65535.0;
	double dist = 0.0;
	lista * aux = l;
	lista * minimo = NULL;
	char dir = 1;
	char dir_aux = 1;
	// Buscamos donde se choca el rayo
	while(1){
		if(aux->radio==-1.0)
		{
			dist = toca_triangulo(cam, pixel, aux->punto[0],aux->punto[1],aux->punto[2]);
		} else{
			dist = toca_esfera(cam,pixel,*aux->punto,aux->radio,&dir);
		}
		if(dist>0.0 && dist<min){
			min = dist;
			minimo = aux;
			dir_aux = dir;
		}
		if(aux->l==NULL) break;
		aux = aux->l;
	}
	if(min == 65535.0){
		return 0;
	}

	// Obtenemos las coordenadas del punto en el espacio
	punto esfera;
	esfera.x = cam.x + pixel.x*min;
	esfera.y = cam.y + pixel.y*min;
	esfera.z = cam.z + pixel.z*min;

	vector inter;
	inter.x = cam.x - esfera.x;
	inter.y = cam.y - esfera.y;
	inter.z = cam.z - esfera.z;

	/*double light = inter.x*inter.x + inter.y*inter.y + inter.z*inter.z;
	if(light<1.0) light=1.0;
	potencia.r = potencia.r / (4*light);
	potencia.g = potencia.g / (4*light);
	potencia.b = potencia.b / (4*light);*/
	//printf("%f %f %f %f\n", potencia.r,potencia.g,potencia.b,light);

	normalizar(&inter);

	// Calculamos la normal dependiendo si se esta dentro o fuera del circulo
	vector * normal = calloc(1,sizeof(vector));
	if(minimo->radio>=0.0){
		if(dir_aux==1){
			normal->x = esfera.x - minimo->punto->x;
			normal->y = esfera.y - minimo->punto->y;
			normal->z = esfera.z - minimo->punto->z;
		} else{
			normal->x = minimo->punto->x - esfera.x;
			normal->y = minimo->punto->y - esfera.y;
			normal->z = minimo->punto->z - esfera.z;
		}
	} else{
		if(minimo->normales == NULL){
			vector e1 = {minimo->punto[1].x-minimo->punto[0].x, minimo->punto[1].y-minimo->punto[0].y, minimo->punto[1].z-minimo->punto[0].z};
			vector e2 = {minimo->punto[2].x-minimo->punto[0].x, minimo->punto[2].y-minimo->punto[0].y, minimo->punto[2].z-minimo->punto[0].z};
			crossproduct(&e1,&e2,normal);
			double dotp=dotproduct(normal, &pixel);
			if(dotp>0.0){normal->x=-normal->x;normal->y=-normal->y;normal->z=-normal->z;}
		
		}else{
			vector d1 = {minimo->punto[0].x-esfera.x, minimo->punto[0].y-esfera.y, minimo->punto[0].z-esfera.z};
			vector d2 = {minimo->punto[1].x-esfera.x, minimo->punto[1].y-esfera.y, minimo->punto[1].z-esfera.z};
			vector d3 = {minimo->punto[2].x-esfera.x, minimo->punto[2].y-esfera.y, minimo->punto[2].z-esfera.z};

			vector aux;
			crossproduct(&d1,&d2,&aux);
			double area3 = sqrt(aux.x*aux.x+aux.y*aux.y+aux.z*aux.z)/2.0;

			crossproduct(&d2,&d3,&aux);
			double area1 = sqrt(aux.x*aux.x+aux.y*aux.y+aux.z*aux.z)/2.0;

			crossproduct(&d3,&d1,&aux);
			double area2 = sqrt(aux.x*aux.x+aux.y*aux.y+aux.z*aux.z)/2.0;

			double sumatorio = area1+area2+area3;
			area1=area1/sumatorio; area2=area2/sumatorio; area3=area3/sumatorio;

			normal->x = minimo->normales[0].x*area1+minimo->normales[1].x*area2+minimo->normales[2].x*area3;
			normal->y = minimo->normales[0].y*area1+minimo->normales[1].y*area2+minimo->normales[2].y*area3;
			normal->z = minimo->normales[0].z*area1+minimo->normales[1].z*area2+minimo->normales[2].z*area3;

		}
	}
	normalizar(normal);

	int flag=0;

	vector salida;
				
	// Se calcula la reflexion solo si es necesario
	if(minimo->propiedades->Krfl->r!=0.0 && 
		minimo->propiedades->Krfl->g!=0.0 && 
		minimo->propiedades->Krfl->b!=0.0)
	{
			double factor = normal->x * -pixel.x + normal->y * -pixel.y + normal->z * -pixel.z;
			salida.x = -pixel.x - 2 * (-pixel.x - factor * normal->x);
			salida.y = -pixel.y - 2 * (-pixel.y - factor * normal->y);
			salida.z = -pixel.z - 2 * (-pixel.z - factor * normal->z);
			normalizar(&salida);
			flag=1;
	}

	// Se calcula la refraccion solo si es necesario
	if(minimo->propiedades->Krfr->r!=0.0 && 
		minimo->propiedades->Krfr->g!=0.0 && 
		minimo->propiedades->Krfr->b!=0.0)
	{

			double n_refraction = minimo->propiedades->indice_ref; // Aire-Cristal 1.00/1.52=0.65
			if(dir_aux == 0) n_refraction = 1.0/n_refraction;
			double dot = -normal->x * pixel.x + -normal->y * pixel.y + -normal->z * pixel.z;
			double k = 1.0 - n_refraction * n_refraction * (1.0 - dot * dot);
			
			salida.x = n_refraction * pixel.x + (n_refraction * dot - sqrt(k)) * normal->x;
			salida.y = n_refraction * pixel.y + (n_refraction * dot - sqrt(k)) * normal->y;
			salida.z = n_refraction * pixel.z + (n_refraction * dot - sqrt(k)) * normal->z;
			normalizar(&salida);
			flag=1;
	}

	if(flag==0)
	{
		if(primera_vez==0){
			photon * photon = calloc(1, sizeof(vector));
			photon->color = calloc(1,sizeof(color));
			photon->direccion = calloc(1, sizeof(vector));

			photon->color->r=potencia.r; photon->color->g=potencia.g; photon->color->b=potencia.b;
			photon->direccion->x=-pixel.x; photon->direccion->y=-pixel.y; photon->direccion->z=-pixel.z;

			double pos[3]={esfera.x,esfera.y,esfera.z};

			if(kd_insert(mapaFotones, pos, (void *) photon)==-1) printf("Error al añadir\n");;
			total_insertado++;
		}
		srand48(clock());
		double aleatorio=0.0;
		aleatorio = drand48();
		double probdifusa = (minimo->propiedades->color->r+minimo->propiedades->color->g+minimo->propiedades->color->b)/3;

		if(aleatorio<=probdifusa){
			obtener_rayo(*normal, &salida);
			brdf(salida, inter, *normal, minimo, &potencia);
			enviar_photon(salida, esfera, potencia, 0);
		} else if(aleatorio<=probdifusa+minimo->propiedades->ks){
			double factor = normal->x * -pixel.x + normal->y * -pixel.y + normal->z * -pixel.z;
			salida.x = -pixel.x - 2 * (-pixel.x - factor * normal->x);
			salida.y = -pixel.y - 2 * (-pixel.y - factor * normal->y);
			salida.z = -pixel.z - 2 * (-pixel.z - factor * normal->z);
			normalizar(&salida);
			brdf(salida, inter, *normal, minimo, &potencia);
			enviar_photon(salida, esfera, potencia, 0);
		}
		free(normal);
		return 1;
	}

	enviar_photon(salida, esfera, potencia, 0);

	free(normal);
	return 1;
}


int mapa_fotones()
{
	mapaFotones = kd_create(3);
	luces * aux = lights;

	while(1){

		printf("%g %g %g\n",aux->punto->x,aux->punto->y,aux->punto->z);

		int i;
		double aleatorio1, aleatorio2, aleatorio3;

		color potencia;
		potencia.r=aux->color->r / NUM_RAYOS;
		potencia.g=aux->color->g / NUM_RAYOS;
		potencia.b=aux->color->b / NUM_RAYOS;
		for (i = 0; i < NUM_RAYOS; i++){
			//se eligen la inclinación y el acimut por montecarlo
			srand(clock());
			do{
				aleatorio1 = (double) ((rand()%201)-100) / 100.0;
				aleatorio2 = (double) ((rand()%201)-100) / 100.0;
				aleatorio3 = (double) ((rand()%201)-100) / 100.0;
			} while(aleatorio1*aleatorio1+aleatorio2*aleatorio2+aleatorio3*aleatorio3 > 1);

			vector photon = {aleatorio1, aleatorio2, aleatorio3};
			normalizar(&photon);

			punto point;
			point.x=aux->punto->x;point.y=aux->punto->y;point.z=aux->punto->z;

			//enviar_photon(photon, point, *aux->color, 1);
			enviar_photon(photon, point, potencia, 1);
		}
		printf("%d\n", i);
		if(aux->l==NULL) break;
		aux = aux->l;
	}

	printf("Insertados: %d\n", total_insertado);
	return 1;
}

int saturacion_color(color * col)
{
	if(col->r>255.0) col->r = 255.0;
	if(col->g>255.0) col->g = 255.0;
	if(col->b>255.0) col->b = 255.0;

	return 1;
}

void* trabajador(void * argumentos)
{
	int indice = *((int *)argumentos);
	// i y d se usan para crear los rayos
	double i,d;
	// i_i y d_d se usan debido a comportamientos extraños con altas resoluciones
	int i_i,d_d;

	double i_ancho=1.0/ancho;
	double i_alto=1.0/alto;

	int index_buffer = ((ancho*alto*3)/NUM_THREADS)*indice;

	for(i = (1.0/NUM_THREADS)*(NUM_THREADS-indice), i_i=alto/NUM_THREADS-1; i_i>=0.0; i=i-i_alto,i_i--)
	{
		for (d = 0.0, d_d=0; d_d<ancho; d=d+i_ancho,d_d++)
		{

			vector pixel = {d-camara.x, i-camara.y, 0.0-camara.z};
			color col;
			normalizar(&pixel);
			col = calcular_luz(pixel,camara,RECURSIONES);
			// Desnormalizamos la luz
			col.r = col.r * 255.0; col.g = col.g * 255; col.b = col.b * 255.0;
			saturacion_color(&col);
			img_buff[index_buffer]=(unsigned char)col.r;
			img_buff[index_buffer+1]=(unsigned char)col.g;
			img_buff[index_buffer+2]=(unsigned char)col.b;
			index_buffer += 3;
		}
		progreso++;
		if(progreso>=incrementador){
			porcentaje++;progreso=0;
			int p;
			printf("\r[");
			for(p=0;p<=porcentaje/5;p++) printf("=");
			for(;p<20;p++) printf(" ");
			printf("][%02d%%]", porcentaje);
			fflush(stdout);
		}
	}
	return NULL;
}

#ifdef OPENGL
void mostrar()
{
	glDrawPixels( ancho, alto, GL_RGB, GL_UNSIGNED_BYTE, img_buff);
	glutSwapBuffers();
}

void* escritura(void * args)
{
	pthread_t * threads= (pthread_t *) args;
	// wait for each thread to complete
	for(int index = 0; index < NUM_THREADS; ++index )
	{
		pthread_join( threads[ index ], NULL );
	}

	printf("\nEscribiendo imagen...");
	FILE * imagen;
	imagen = fopen(img, "w");
	fprintf(imagen, "P3 %d %d 255\n", ancho, alto);
	int i;
	for(i=0; i<alto*ancho*3;i+=3)
		fprintf(imagen, "%d %d %d  ", img_buff[i],img_buff[i+1],img_buff[i+2]);
	fclose(imagen);
	printf("\n");

	exit(0);
	return NULL;
}
#endif

int main(int argc, char ** argv)
{
	int opt;
	int obj=0;
	// Con este bucle obtenemos los parametros dados
	while ((opt = getopt (argc, argv, "o:i:j:h")) != -1){
		switch(opt)
			{
				case 'o':
					img = optarg;
					break;
				case 'i':
					scn = optarg;
					break;
				case 'j':
					obj=1;
					scn = optarg;
					break;
				case 'h':
					printf("Use -o para nombre de imagen, -i para nombre de escena y -j para objs\n");
					return 0;
			}
	}
	FILE * escena;

	escena = fopen(scn, "r");

	if(obj==0){
		parser(escena);
	}else{
		parserOBJ(escena);
		luces * a = calloc(1, sizeof(luces));
		a->punto=calloc(1,sizeof(punto));
		a->color=calloc(1,sizeof(color));

		a->punto->x=0.5; a->punto->y=0.5; a->punto->z=0.0;
		a->color->r=5.0; a->color->g=5.0; a->color->b=5.0;
		
		lights=a;
	}
	fclose(escena);

	img_buff=malloc(ancho*alto*sizeof(char)*3);
	incrementador = alto/100;

	printf("LLenando mapa de fotones\n");
	mapa_fotones();
	printf("Mapa llenado\n");

	#ifdef OPENGL
	glutInit( &argc, argv );
	glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE );
	glutInitWindowSize( ancho, alto );
	glutSetWindow(glutCreateWindow( "Preview" ));
	glutDisplayFunc( mostrar );
	glutIdleFunc(mostrar);
	glRasterPos2f(-1,1);
	glPixelZoom( 1, -1 );
	#endif

	pthread_t threads[ NUM_THREADS ];
	int thread_args[ NUM_THREADS ];
	int index;

	for( index = 0; index < NUM_THREADS; ++index )
	{
		thread_args[ index ] = index;
		pthread_create( threads + index , NULL, trabajador, thread_args + index );
	}

	#ifdef OPENGL
		pthread_t espera;
		pthread_create(&espera, NULL, escritura, threads);
		glutMainLoop();
	
	#else

		// wait for each thread to complete
		for( index = 0; index < NUM_THREADS; ++index )
		{
			pthread_join( threads[ index ], NULL );
		}

		printf("\nEscribiendo imagen...");
		FILE * imagen;
		imagen = fopen(img, "w");
		fprintf(imagen, "P3 %d %d 255\n", ancho, alto);
		int i;
		for(i=0; i<alto*ancho*3;i+=3)
			fprintf(imagen, "%d %d %d  ", img_buff[i],img_buff[i+1],img_buff[i+2]);
		fclose(imagen);
		printf("\n");
	#endif
	return 0;
}