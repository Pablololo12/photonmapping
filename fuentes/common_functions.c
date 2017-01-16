/*
 * Autores: Pablo Hernandez Almudi y Mario Arcega Sanjuan
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tipos.h"
#include "common_functions.h"

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