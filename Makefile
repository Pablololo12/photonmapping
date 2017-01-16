#
# Autores: Pablo Hernandez Almudi y Mario Arcega Sanjuan
#

CC = gcc
CFLAGS =-Wall -Ofast -lm -lpthread -std=gnu11
OGLFLAGS = -Wno-deprecated-declarations

SRC=fuentes/tipos.h fuentes/kdtree.c fuentes/kdtree.h fuentes/common_functions.c fuentes/common_functions.h fuentes/direct.c fuentes/direct.h fuentes/indirecta.c fuentes/indirecta.h fuentes/mapa_fotones.c fuentes/mapa_fotones.h
OBJ=fuentes/common_functions.o fuentes/direct.o fuentes/indirecta.o fuentes/mapa_fotones.o fuentes/kdtree.o

all: $(OBJ)
	$(CC) $(CFLAGS) -o ray-trace fuentes/main.c $(OBJ)

#compilar Linux
ogl: $(OBJ)
	$(CC) $(CFLAGS) $(OGLFLAGS) -lOpenGL -lglut -DOPENGL -o ray-trace fuentes/main.c $(OBJ)

#compilar Mac
oglm: $(OBJ)
	$(CC) $(CFLAGS) $(OGLFLAGS) -framework OpenGL -framework GLUT -DOPENGL -o ray-trace fuentes/main.c $(OBJ)

clean:
	$(RM) $(OBJ) ray-trace