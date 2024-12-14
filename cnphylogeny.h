#ifndef CNPHYLOGENY_H
#define CNPHYLOGENY_H

#include <stdlib.h>

/**
 * @brief A copy number
 */
typedef unsigned char copy_num;

/**
 * @brief A square, right stochastic matrix
 */
struct prob_matrix {
    copy_num order; /** The size of the matrix */
    double probs[]; /** The elements of the matrix */
};

/**
 * @brief Create a probability matrix
 *
 * Note that the new matrix stores a copy of `probs`, not a reference to it.
 *
 * @param order The size of the matrix
 * @param probs The elements of the matrix, flattened into an array (often, a
 *              compound literal)
 * @return A pointer to the new matrix
 */
struct prob_matrix *prob_matrix_new(copy_num order, double *probs);


/**
 * @brief A binary tree node that stores a copy number profile (CNP)
 * 
 * The fields `left` and `right` should be `NULL` for leaf nodes. For nodes with
 * only one child, `right` should be `NULL`. In other words, if `left` is
 * `NULL`, the node should have no children.
 */
struct cnp_node {
    struct cnp_node *left; /** The left child */
    struct cnp_node *right; /** The right child */
    size_t len; /** The length of the CNP */
    copy_num bins[]; /** The bins of the CNP */
};

/**
 * @brief Create a new CNP node
 *
 * This function assigns each bin a copy number of 0.
 *
 * @param len The length of the CNP
 * @param left The left child
 * @param right The right child
 * @return A pointer to the new node
 */
struct cnp_node *cnp_node_new(
    size_t len,
    struct cnp_node *left,
    struct cnp_node *right
);

/**
 * @brief Free a CNP node and its children
 *
 * @param node The node to free
 */
void cnp_node_free(struct cnp_node *node);


/**
 * @brief Use Gibbs sampling to optimize a phylogeny
 * 
 * @param root A pointer to the root of the phylogeny
 * @param neighbor_probs A probability matrix that defines the likelihood of
 *                       each possible configuration of two neighboring bins
 * @param mutation_probs A probability matrix that defines the likelihood of
 *                       each possible configuration of a parent and child bins
 * @param burn_in The number of samples to discard
 * @param sample_rate The rate at which to record samples
 * @param sample_count The number of samples to record
 */
void phylogeny_optimize(
    struct cnp_node *root,
    struct prob_matrix *neighbor_probs,
    struct prob_matrix *mutation_probs,
    int burn_in,
    int sample_rate,
    int sample_count
);


#endif
