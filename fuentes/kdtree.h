/*
 * Autores: Pablo Hernandez Almudi y Mario Arcega Sanjuan
 */

#ifndef KDTREE_H
#define KDTREE_H

#define MAX_DIM 3
struct kd_node_t{
    double x[MAX_DIM];
    struct kd_node_t *left, *right;
};

struct kd_node_t* make_tree(struct kd_node_t *, int, int, int);

void nearest(struct kd_node_t *, struct kd_node_t *, int, int,
        struct kd_node_t **, double *);

#endif