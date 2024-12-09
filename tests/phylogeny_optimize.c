#include "../cnphylogeny.h"

#include <assert.h>
#include <string.h>

int main()
{
    struct cnp_node *root = cnp_node_new(
        5,
        cnp_node_new(
            5,
            cnp_node_new(5, NULL, NULL),
            cnp_node_new(
                5,
                cnp_node_new(5, NULL, NULL),
                cnp_node_new(5, NULL, NULL)
            )
        ),
        NULL
    );
    struct cnp_node *leaf1 = root->left->left;
    struct cnp_node *leaf2 = root->left->right->left;
    struct cnp_node *leaf3 = root->left->right->right;

    memcpy(root->bins, (copy_num []) { 2, 2, 2, 2, 2 }, 5);
    memcpy(leaf1->bins, (copy_num []) { 2, 2, 1, 1, 1 }, 5);
    memcpy(leaf2->bins, (copy_num []) { 2, 0, 2, 2, 2 }, 5);
    memcpy(leaf3->bins, (copy_num []) { 1, 1, 1, 2, 2 }, 5);

    struct prob_matrix *neighbor_probs = prob_matrix_new(3, (double []) {
        0.9, 0.05, 0.05,
        0.05, 0.9, 0.05,
        0.05, 0.05, 0.9,
    });
    struct prob_matrix *mutation_probs = prob_matrix_new(3, (double []) {
        1, 0, 0,
        0.005, 0.99, 0.005,
        0.005, 0.005, 0.99,
    });

    phylogeny_optimize(root, neighbor_probs, mutation_probs, 0, 1, 1);

    assert(!memcmp(root->bins, (copy_num []) { 2, 2, 2, 2, 2 }, 5));
    assert(!memcmp(root->left->bins, (copy_num []) { 2, 2, 2, 2, 2 }, 5));
    assert(!memcmp(leaf1->bins, (copy_num []) { 2, 2, 1, 1, 1 }, 5));
    assert(!memcmp(root->left->right->bins, (copy_num []) { 2, 1, 2, 2, 2 }, 5));
    assert(!memcmp(leaf2->bins, (copy_num []) { 2, 0, 2, 2, 2 }, 5));
    assert(!memcmp(leaf3->bins, (copy_num []) { 1, 1, 1, 2, 2 }, 5));

    free(neighbor_probs);
    free(mutation_probs);
    cnp_node_free(root);

    return 0;
}
