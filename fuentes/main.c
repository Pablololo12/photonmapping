/*
 * Autores: Pablo Hernandez Almudi y Mario Arcega Sanjuan
 */

#include "tipos.h"
#include "kdtree.h"
#include "common_functions.h"
#include "direct.h"
#include "mapa_fotones.h"

#ifdef OPENGL
	#ifdef __APPLE__
	#include <GLUT/glut.h>
	#else
	#include <GL/glut.h>
	#endif
#endif

volatile sig_atomic_t progreso=0;
volatile sig_atomic_t porcentaje=0;

unsigned char * img_buff;

// Valores de resolución y posición de la cámara
int ancho=500;
int alto=500;
punto camara={0.5,0.5,-3.0};
int incrementador;
int total_insertado=0;

// Nombres por defecto de los ficheros
char * img={"imagen.ppm"};
char * scn={"escenas/imagen.scn"};

// Punteros a las listas donde se almacenan las esferas y los focos
lista * l = NULL;
luces * lights = NULL;

struct kdtree * mapaFotones;

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