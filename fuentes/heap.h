
#ifndef HEAP_H
#define HEAP_H

#include "tipos.h"

typedef int type;

struct heap
{
	unsigned int size; // Size of the allocated memory (in number of items)
	unsigned int count; // Count of the elements in the heap
	kd_node_t *data; // Array with the elements
};

void heap_init(struct heap *restrict h);
void heap_push(struct heap *restrict h, type value);
void heap_pop(struct heap *restrict h);

// Returns the biggest element in the heap
#define heap_front(h) (*(h)->data)

// Frees the allocated memory
#define heap_term(h) (free((h)->data))

#endif