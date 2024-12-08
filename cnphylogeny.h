#ifndef CNPHYLOGENY_H
#define CNPHYLOGENY_H

#include <stdlib.h>

typedef unsigned char copy_num;

struct prob_matrix {
    copy_num order;
    double probs[];
};

struct prob_matrix *prob_matrix_new(copy_num order, double *probs);

struct cnp_node {
    struct cnp_node *left;
    struct cnp_node *right;
    size_t len;
    copy_num bins[];
};

struct cnp_node *cnp_node_new(
    size_t len,
    struct cnp_node *left,
    struct cnp_node *right
);

void cnp_node_free(struct cnp_node *node);

void phylogeny_optimize(
    struct cnp_node *root,
    struct prob_matrix *neighbor_probs,
    struct prob_matrix *mutation_probs
);

#endif
