#include "cnphylogeny.h"

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>


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

static double probs_lookup(struct prob_matrix *matrix, copy_num i, copy_num j)
{
    return matrix->probs[i * matrix->order + j];
}


struct gibbs_node {
    size_t len;
    copy_num state_count;
    copy_num *bins;
    copy_num *prev;
    int *counts;
    struct gibbs_node *parent;
    struct gibbs_node *left;
    struct gibbs_node *right;
};

static struct gibbs_node *gibbs_node_new(
    struct cnp_node *src,
    copy_num state_count,
    struct gibbs_node *parent
)
{
    if (!src) return NULL;

    struct gibbs_node *node = malloc(sizeof(struct gibbs_node));
    node->len = src->len;
    node->state_count = state_count;
    node->bins = malloc(src->len);
    memcpy(node->bins, src->bins, src->len);
    node->prev = malloc(src->len);
    memcpy(node->prev, src->bins, src->len);
    node->counts = calloc(src->len * state_count, sizeof(int));
    node->parent = parent;
    node->left = gibbs_node_new(src->left, state_count, node);
    node->right = gibbs_node_new(src->right, state_count, node);

    return node;
}

static int *counts_lookup(struct gibbs_node *node, size_t bin, copy_num state)
{
    return node->counts + bin * node->state_count + state;
}

static void gibbs_iteration(
    struct gibbs_node *node,
    struct prob_matrix *neighbor_probs,
    struct prob_matrix *mutation_probs,
    bool count
)
{
    if (!node || !node->left) return;

    gibbs_iteration(node->left, neighbor_probs, mutation_probs, count);
    gibbs_iteration(node->right, neighbor_probs, mutation_probs, count);

    if (node->parent) {
        for (int i = 0; i < node->len; i++) {
            double max_prob = -INFINITY;
            copy_num best = 0;
            for (copy_num s = 0; s < node->state_count; s++) {
                double prob = (
                    probs_lookup(mutation_probs, node->parent->prev[i], s) +
                    probs_lookup(mutation_probs, s, node->left->prev[i])
                );
                if (node->right)
                    prob += probs_lookup(
                        mutation_probs,
                        s, node->right->prev[i]
                    );
                if (i > 0)
                    prob += probs_lookup(neighbor_probs, node->prev[i - 1], s);
                if (i < node->len - 1)
                    prob += probs_lookup(neighbor_probs, s, node->prev[i + 1]);

                if (prob > max_prob) {
                    best = s;
                    max_prob = prob;
                }
            }

            node->bins[i] = best;
            if (count) (*counts_lookup(node, i, best))++;
        }
    }

    memcpy(node->left->prev, node->left->bins, node->len);
    if (node->right) memcpy(node->right->prev, node->right->bins, node->len);
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

static void cnp_get_mode(struct cnp_node *node, struct gibbs_node *src)
{
    if (!node || !node->left) return;

    cnp_get_mode(node->left, src->left);
    cnp_get_mode(node->right, src->right);

    if (!src->parent) return;

    for (int i = 0; i < node->len; i++) {
        int max_count = 0;
        copy_num mode;
        for (int s = 0; s < src->state_count; s++) {
            int count = *counts_lookup(src, i, s);
            if (count > max_count) {
                mode = s;
                max_count = count;
            }
        }
        node->bins[i] = mode;
    }
}


void phylogeny_optimize(
    struct cnp_node *root,
    struct prob_matrix *neighbor_probs,
    struct prob_matrix *mutation_probs,
    int burn_in,
    int sample_rate,
    int sample_count
)
{
    struct gibbs_node *gibbs_root = gibbs_node_new(
        root,
        neighbor_probs->order,
        NULL
    );

    for (int i = 0; i < burn_in; i++)
        gibbs_iteration(gibbs_root, neighbor_probs, mutation_probs, false);
    for (int i = 0; i < sample_count; i++) {
        gibbs_iteration(gibbs_root, neighbor_probs, mutation_probs, true);
        for (int j = 0; j < sample_rate - 1; j++)
            gibbs_iteration(gibbs_root, neighbor_probs, mutation_probs, false);
    }

    cnp_get_mode(root, gibbs_root);

    gibbs_node_free(gibbs_root);
}
