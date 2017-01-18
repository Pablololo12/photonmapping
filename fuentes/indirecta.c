
#include "indirecta.h"

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

	while(kd_res_size(busqueda)<NUM_BUSCAR)
	{
		i+=0.05;
		kd_res_free(busqueda);
		busqueda=kd_nearest_range(mapaFotones, pos, i);
		if(i>1.0) return acumulador;
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

		// Filtro de cono
		acum2 = 1 - (acum/(1.0*i));

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
	area = area * (1-(2/(3*1.0)));

	acumulador.r=(acumulador.r/area)*POT_INDIRECTA;
	acumulador.g=(acumulador.g/area)*POT_INDIRECTA;
	acumulador.b=(acumulador.b/area)*POT_INDIRECTA;

	if(acumulador.r>1.0) acumulador.r=1.0;
	if(acumulador.g>1.0) acumulador.g=1.0;
	if(acumulador.b>1.0) acumulador.b=1.0;

	return acumulador;
}