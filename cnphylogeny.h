#ifndef CNPHYLOGENY_H
#define CNPHYLOGENY_H

#include <stdbool.h>
#include <stddef.h>

typedef unsigned char copy_num;

struct prob_matrix {
    copy_num order;
    double probs[];
};

struct prob_matrix *prob_matrix_new(copy_num order, double *probs);
bool prob_matrix_valid(struct prob_matrix *probs);

struct cnp_node {
    struct cnp_node *left;
    struct cnp_node *right;
    size_t len;
    copy_num bins[];
};

void phylogeny_optimize(
    struct cnp_node *root,
    struct prob_matrix *neighbor_probs,
    struct prob_matrix *mutation_probs
);

#endif
