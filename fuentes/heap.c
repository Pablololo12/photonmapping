//https://gist.github.com/martinkunev/1365481

#include <unistd.h>
#include <stdlib.h>

#include "heap.h"

#define CMP(a, b) ((a) <= (b))

static const unsigned int base_size = 50;

// Prepares the heap for use
void heap_init(struct heap *restrict h)
{
    *h = (struct heap){
        .size = base_size,
        .count = 0,
        .data = malloc(sizeof(type) * base_size)
    };
    if (!h->data) _exit(1); // Exit if the memory allocation fails
}

// Inserts element to the heap
void heap_push(struct heap *restrict h, type value)
{
    unsigned int index, parent;

    // Resize the heap if it is too small to hold all the data
    if (h->count == h->size)
    {
        h->size <<= 1;
        h->data = realloc(h->data, sizeof(type) * h->size);
        if (!h->data) _exit(1); // Exit if the memory allocation fails
    }

    // Find out where to put the element and put it
    for(index = h->count++; index; index = parent)
    {
        parent = (index - 1) >> 1;
        if CMP(h->data[parent], value) break;
        h->data[index] = h->data[parent];
    }
    h->data[index] = value;
}

// Removes the biggest element from the heap
void heap_pop(struct heap *restrict h)
{
    unsigned int index, swap, other;

    // Remove the biggest element
    type temp = h->data[--h->count];

    // Resize the heap if it's consuming too much memory
    if ((h->count <= (h->size >> 2)) && (h->size > base_size))
    {
        h->size >>= 1;
        h->data = realloc(h->data, sizeof(type) * h->size);
        if (!h->data) _exit(1); // Exit if the memory allocation fails
    }

    // Reorder the elements
    for(index = 0; 1; index = swap)
    {
        // Find the child to swap with
        swap = (index << 1) + 1;
        if (swap >= h->count) break; // If there are no children, the heap is reordered
        other = swap + 1;
        if ((other < h->count) && CMP(h->data[other], h->data[swap])) swap = other;
        if CMP(temp, h->data[swap]) break; // If the bigger child is less than or equal to its parent, the heap is reordered

        h->data[index] = h->data[swap];
    }
    h->data[index] = temp;
}