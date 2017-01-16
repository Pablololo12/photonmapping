
#include "mapa_fotones.h"

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