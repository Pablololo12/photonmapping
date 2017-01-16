CC = gcc
CFLAGS =-Wall -Ofast -lm -lpthread -std=gnu11
OGLFLAGS = -Wno-deprecated-declarations

SRC=fuentes/tipos.h fuentes/kdtree.c fuentes/kdtree.h
OBJ=fuentes/kdtree.o

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