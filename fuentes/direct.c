
#include "direct.h"

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