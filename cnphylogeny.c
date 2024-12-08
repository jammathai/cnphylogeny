#include "cnphylogeny.h"

#include <math.h>
#include <stdlib.h>

struct prob_matrix *prob_matrix_new(copy_num order, double *probs)
{
    int probs_len = order * order;
    struct prob_matrix *matrix = malloc(
        sizeof(struct prob_matrix) + probs_len * sizeof(double)
    );
    matrix->order = order;
    for (int i = 0; i < probs_len; i++) matrix->probs[i] = log(probs[i]);
    return matrix;
}

struct cnp_node *cnp_node_new(
    size_t len,
    struct cnp_node *left,
    struct cnp_node *right
)
{
    struct cnp_node *node = calloc(1, sizeof(struct cnp_node) + len);
    node->len = len;
    node->left = left;
    node->right = right;
    return node;
}

void cnp_node_free(struct cnp_node *node)
{
    if (!node) return;
    cnp_node_free(node->left);
    cnp_node_free(node->right);
    free(node);
}
