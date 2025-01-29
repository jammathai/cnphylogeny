#include "../include/cnphylogeny.h"

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>


struct gibbs_node {
    copy_num *bins;
    copy_num *prev;
    int **counts;
    struct gibbs_node *parent;
    struct gibbs_node *left;
    struct gibbs_node *right;
};


static struct gibbs_node *gibbs_node_new(
    struct cnp_node *src,
    struct gibbs_node *parent
);
static void gibbs_iteration(struct gibbs_node *node, bool count);
static void gibbs_node_free(struct gibbs_node *node);
static void cnp_get_mode(struct cnp_node *node, struct gibbs_node *src);


double **prob_matrix_new(double *probs)
{
    copy_num order = max_copy_num + 1;
    double **matrix = malloc(
        order * sizeof(double *) + order * order * sizeof(double)
    );

    double *row = (double *) (matrix + order);
    for (int i = 0; i < order; i++) {
        matrix[i] = row;
        row += order;
        for (int j = 0; j < order; j++)
            matrix[i][j] = log(probs[i * order + j]);
    }

    return matrix;
}


struct cnp_node *cnp_node_new(
    int id,
    copy_num *cnp,
    struct cnp_node *left,
    struct cnp_node *right
)
{
    struct cnp_node *node = calloc(1, sizeof(struct cnp_node) + cnp_len);
    node->id = id;
    if (cnp) memcpy(node->bins, cnp, cnp_len);
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


void phylogeny_optimize(struct cnp_node *root, int burn_in, int sample_count)
{
    struct gibbs_node *gibbs_root = gibbs_node_new(root, NULL);

    for (int i = 0; i < burn_in; i++)
        gibbs_iteration(gibbs_root, false);
    for (int i = 0; i < sample_count; i++)
        gibbs_iteration(gibbs_root, true);

    cnp_get_mode(root, gibbs_root);

    gibbs_node_free(gibbs_root);
}


static struct gibbs_node *gibbs_node_new(
    struct cnp_node *src,
    struct gibbs_node *parent
)
{
    if (!src) return NULL;

    struct gibbs_node *node = malloc(sizeof(struct gibbs_node));
    node->bins = malloc(cnp_len);
    memcpy(node->bins, src->bins, cnp_len);
    node->prev = malloc(cnp_len);
    memcpy(node->prev, src->bins, cnp_len);
    node->counts = calloc(
        cnp_len * sizeof(int *) + cnp_len * (max_copy_num + 1), sizeof(int)
    );
    int *row = (int *) (node->counts + cnp_len);
    for (int i = 0; i < cnp_len; i++) {
        node->counts[i] = row;
        row += max_copy_num + 1;
    }
    node->parent = parent;
    node->left = gibbs_node_new(src->left, node);
    node->right = gibbs_node_new(src->right, node);

    return node;
}


static void gibbs_iteration(struct gibbs_node *node, bool count)
{
    if (!node || !node->left) return;

    gibbs_iteration(node->left, count);
    gibbs_iteration(node->right, count);

    if (node->parent) {
        for (int i = 0; i < cnp_len; i++) {
            double max_prob = -INFINITY;
            copy_num best = 0;
            for (copy_num s = 0; s < max_copy_num + 1; s++) {
                double prob = (
                    mutation_probs[node->parent->prev[i]][s] +
                    mutation_probs[s][node->left->prev[i]]
                );
                if (node->right)
                    prob += mutation_probs[s][node->right->prev[i]];
                if (i > 0)
                    prob += neighbor_probs[node->prev[i - 1]][s];
                if (i < cnp_len - 1)
                    prob += neighbor_probs[s][node->prev[i + 1]];

                if (prob > max_prob) {
                    best = s;
                    max_prob = prob;
                }
            }

            node->bins[i] = best;
            if (count) node->counts[i][best]++;
        }
    }

    memcpy(node->left->prev, node->left->bins, cnp_len);
    if (node->right) memcpy(node->right->prev, node->right->bins, cnp_len);
}


static void gibbs_node_free(struct gibbs_node *node)
{
    if (!node) return;

    free(node->bins);
    free(node->prev);
    free(node->counts);
    gibbs_node_free(node->left);
    gibbs_node_free(node->right);
    free(node);
}


static void cnp_get_mode(struct cnp_node *node, struct gibbs_node *src)
{
    if (!node || !node->left) return;

    cnp_get_mode(node->left, src->left);
    cnp_get_mode(node->right, src->right);

    if (!src->parent) return;

    for (int i = 0; i < cnp_len; i++) {
        int max_count = 0;
        copy_num mode;
        for (int s = 0; s < max_copy_num + 1; s++) {
            int count = src->counts[i][s];
            if (count > max_count) {
                mode = s;
                max_count = count;
            }
        }
        node->bins[i] = mode;
    }
}
