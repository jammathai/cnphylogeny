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


static double phylogeny_analyze_internal(struct cnp_node *root);
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


double phylogeny_analyze(struct cnp_node *root)
{
    double result = 0;

    for (int i = 0; i < cnp_len; i++) {
        if (root->left)
            result += mutation_probs[root->bins[i]][root->left->bins[i]];
        if (root->right)
            result += mutation_probs[root->bins[i]][root->right->bins[i]];
    }

    if (root->left)
        result += phylogeny_analyze_internal(root->left);
    if (root->right)
        result += phylogeny_analyze_internal(root->right);

    return result;
}


void phylogeny_optimize(struct cnp_node *root, int burn_in, int sample_count)
{
    struct gibbs_node *gibbs_root = gibbs_node_new(root, NULL);
    double max_score = -INFINITY;

    for (int i = 0; i < burn_in; i++) gibbs_iteration(gibbs_root, false);
    for (int i = 0; i < sample_count; i++) gibbs_iteration(gibbs_root, true);

    cnp_get_mode(root, gibbs_root);

    gibbs_node_free(gibbs_root);
}


static double phylogeny_analyze_internal(struct cnp_node *root)
{
    if (!root->left) return 0;

    double result = 0;

    for (int i = 0; i < cnp_len - 1; i++)
        result += neighbor_probs[root->bins[i]][root->bins[i + 1]];

    for (int i = 0; i < cnp_len; i++) {
        result += mutation_probs[root->bins[i]][root->left->bins[i]];
        if (root->right)
            result += mutation_probs[root->bins[i]][root->right->bins[i]];
    }

    result += phylogeny_analyze(root->left);
    if (root->right) result += phylogeny_analyze(root->right);

    return result;
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
        double probs[max_copy_num + 1];

        for (int i = 0; i < cnp_len; i++) {
            for (copy_num s = 0; s <= max_copy_num; s++) {
                probs[s] = (
                    mutation_probs[node->parent->prev[i]][s] +
                    mutation_probs[s][node->left->prev[i]]
                );
                if (node->right)
                    probs[s] += mutation_probs[s][node->right->prev[i]];
                if (i > 0)
                    probs[s] += neighbor_probs[node->prev[i - 1]][s];
                if (i < cnp_len - 1)
                    probs[s] += neighbor_probs[s][node->prev[i + 1]];
            }

            double random = (double) rand() / RAND_MAX;
            double total = 0;
            for (copy_num s = 0; s <= max_copy_num; s++) {
                total += exp(probs[s]);
                if (random < total) {
                    node->bins[i] = s;
                    if (count) node->counts[i][s]++;
                    break;
                }
            }
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
