#include "../include/cnphylogeny.h"

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>


struct gibbs_node {
    copy_num *bins;
    copy_num *prev;
    copy_num *best;
    struct gibbs_node *parent;
    struct gibbs_node *left;
    struct gibbs_node *right;
};


static double phylogeny_analyze_internal(struct cnp_node *root);
static double gibbs_node_analyze(struct gibbs_node *root);
static struct gibbs_node *gibbs_node_new(
    struct cnp_node *src,
    struct gibbs_node *parent
);
static void gibbs_iteration(struct gibbs_node *node);
static copy_num sample(double *probs, double total);
static void gibbs_node_free(struct gibbs_node *node);
static void gibbs_node_set_best(struct gibbs_node *node);
static void gibbs_node_get_best(struct gibbs_node *src, struct cnp_node *dst);


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


void phylogeny_optimize(struct cnp_node *root, int sample_count)
{
    struct gibbs_node *gibbs_root = gibbs_node_new(root, NULL);
    double max_score = -INFINITY;

    for (int i = 0; i < sample_count; i++) {
        gibbs_iteration(gibbs_root);
        double score = gibbs_node_analyze(gibbs_root);
        if (score > max_score) {
            max_score = score;
            gibbs_node_set_best(gibbs_root);
            printf("Max Score: %lf\n", max_score);
        }
    }

    gibbs_node_get_best(gibbs_root, root);

    gibbs_node_free(gibbs_root);
}


static double gibbs_node_analyze(struct gibbs_node *root)
{
    if (!root->parent) {
        double result = 0;

        for (int i = 0; i < cnp_len; i++) {
            if (root->left)
                result += mutation_probs[root->bins[i]][root->left->bins[i]];
            if (root->right)
                result += mutation_probs[root->bins[i]][root->right->bins[i]];
        }

        if (root->left)
            result += gibbs_node_analyze(root->left);
        if (root->right)
            result += gibbs_node_analyze(root->right);

        return result;
    }

    if (!root->left) return 0;

    double result = 0;

    for (int i = 0; i < cnp_len - 1; i++)
        result += neighbor_probs[root->bins[i]][root->bins[i + 1]];

    for (int i = 0; i < cnp_len; i++) {
        result += mutation_probs[root->bins[i]][root->left->bins[i]];
        if (root->right)
            result += mutation_probs[root->bins[i]][root->right->bins[i]];
    }

    result += gibbs_node_analyze(root->left);
    if (root->right) result += gibbs_node_analyze(root->right);

    return result;
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

    result += phylogeny_analyze_internal(root->left);
    if (root->right) result += phylogeny_analyze_internal(root->right);

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
    node->best = malloc(cnp_len);
    node->parent = parent;
    node->left = gibbs_node_new(src->left, node);
    node->right = gibbs_node_new(src->right, node);

    return node;
}


static void gibbs_iteration(struct gibbs_node *node)
{
    if (!node || !node->left) return;

    gibbs_iteration(node->left);
    gibbs_iteration(node->right);

    if (node->parent) {
        double probs[max_copy_num + 1];
        for (int i = 0; i < cnp_len; i++) {
            double total = 0;

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
                probs[s] = exp(probs[s]);
                // logsumexp
                total += probs[s];
            }

            copy_num s = sample(probs, total);
            node->bins[i] = s;
        }
    }

    memcpy(node->left->prev, node->left->bins, cnp_len);
    if (node->right) memcpy(node->right->prev, node->right->bins, cnp_len);
}


static copy_num sample(double *probs, double total) {
    double r = (double) rand() / RAND_MAX * total;
    double sum = 0;
    for (copy_num s = 0;; s++) {
        sum += probs[s];
        if (r <= sum) return s;
    }
}


static void gibbs_node_free(struct gibbs_node *node)
{
    if (!node) return;

    free(node->bins);
    free(node->prev);
    free(node->best);
    gibbs_node_free(node->left);
    gibbs_node_free(node->right);
    free(node);
}


static void gibbs_node_set_best(struct gibbs_node *node) {
    if (!node) return;

    memcpy(node->best, node->bins, cnp_len);
    gibbs_node_set_best(node->left);
    gibbs_node_set_best(node->right);
}


static void gibbs_node_get_best(struct gibbs_node *src, struct cnp_node *dst) {
    if (!src) return;

    memcpy(dst->bins, src->best, cnp_len);
    gibbs_node_get_best(src->left, dst->left);
    gibbs_node_get_best(src->right, dst->right);
}
