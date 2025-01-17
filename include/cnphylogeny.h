#ifndef CNPHYLOGENY_H
#define CNPHYLOGENY_H


#include <stdlib.h>


/**
 * @brief A copy number
 */
typedef unsigned char copy_num;

/**
 * @brief A binary tree node that stores a copy number profile (CNP)
 * 
 * The fields `left` and `right` should be `NULL` for leaf nodes. For nodes with
 * only one child, `right` should be `NULL`. In other words, if a node has no
 * left child, it shouldn't have a right child either.
 */
struct cnp_node {
    struct cnp_node *left; /** The left child */
    struct cnp_node *right; /** The right child */
    copy_num bins[]; /** The bins of the CNP */
};


/**
 * @brief The length of a CNP
 *
 * This variable must be defined.
 */
extern size_t cnp_len;

/**
 * @brief The maximum possible copy number
 *
 * This variable must be defined.
 */
extern copy_num max_copy_num;

/**
 * @brief A matrix that stores the natural log probabilities of every possible
 *        configuration of two neighboring bins
 *
 * `neighbor_probs[i][j]` gives the log probability of observing a bin with copy
 * number `i` followed by a bin with copy number `j` (in the same CNP).
 *
 * This variable must be defined.
 */
extern double **neighbor_probs;

/**
 * @brief A matrix that stores the natural log probabilities of every possible
 *        configuration of a parent and child bin
 *
 * `neighbor_probs[i][j]` gives the log probability of observing a bin with copy
 * number `i` in a parent CNP and a bin with copy number `j` in the
 * corresponding position of its child.
 *
 * This variable must be defined.
 */
extern double **mutation_probs;


/**
 * @brief Allocate and assign a probability matrix
 *
 * The new matrix will have order `max_copy_num + 1`. Note that the
 * probabilities will be stored as natural log probabilities; be careful when
 * modifying an existing probability matrix.
 *
 * @param probs A pointer to the elements of the matrix, flattened into an array
 *              (often, a compound literal)
 * @return A pointer to the new matrix
 */
double **prob_matrix_new(double *probs);


/**
 * @brief Create a new CNP node
 *
 * This function assigns each bin a copy number of 0.
 *
 * @param left The left child
 * @param right The right child
 * @return A pointer to the new node
 */
struct cnp_node *cnp_node_new(
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
 * @param burn_in The number of samples to discard
 * @param sample_rate The rate at which to record samples
 * @param sample_count The number of samples to record
 */
void phylogeny_optimize(
    struct cnp_node *root,
    int burn_in,
    int sample_rate,
    int sample_count
);


#endif
