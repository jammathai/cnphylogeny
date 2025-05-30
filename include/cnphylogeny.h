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
    int id; /** The ID of the CNP */
    struct cnp_node *left; /** The left child */
    struct cnp_node *right; /** The right child */
    copy_num bins[]; /** The bins of the CNP */
};


/**
 * @brief The length of a CNP
 *
 * This variable must be defined.
 */
extern int cnp_len;

/**
 * @brief The maximum possible copy number
 *
 * This variable must be defined.
 */
extern copy_num max_copy_num;

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
 * @param id The ID of the CNP
 * @param cnp A pointer a CNP, or `NULL` to initialize all bins to zero
 * @param left The left child
 * @param right The right child
 * @return A pointer to the new node
 */
struct cnp_node *cnp_node_new(
    int id,
    copy_num *cnp,
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
 * @brief Compute a value proportional to the likelihood of a phylogeny
 * 
 * @param root A pointer to the root of the phylogeny
 * @return A value proportional to the likelihood, which can be used to compare
 *         possible phylogenies
 */
double phylogeny_analyze(struct cnp_node *root);


/**
 * @brief Use Gibbs sampling to optimize a phylogeny
 *
 * @param root A pointer to the root of the phylogeny
 * @param sample_count The number of samples to take
 */
void phylogeny_optimize(struct cnp_node *root, int sample_count);


#endif
